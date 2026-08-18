"""Detect data leakage between the train and test halves of a dataset split.

Stages both halves to FASTA, runs an MMseqs2 all-vs-all search of test against
train, and reports how much of the test set has a near-identical counterpart in
training - sequences a model can score correctly by memorisation alone.

`run` is the entry point; the CLI is a thin wrapper around it. The report layout
is defined in `genomic_benchmarks_qc.utils.naming`.
"""

import logging
import os
import stat
import tempfile
from pathlib import Path
from typing import Optional
import shutil
import pandas as pd

from genomic_benchmarks_qc.report.report_generator import generate_splits_html_report, generate_simple_report
from genomic_benchmarks_qc.utils.mmseqs_summary import summarize_mmseqs_output, build_mmseqs_export_frame
from genomic_benchmarks_qc.utils import mmseqs_runtime
from genomic_benchmarks_qc.utils.input_utils import (
    ensure_directory,
    log_failures,
    setup_logger,
    stream_files_to_sequences,
    append_fasta_record,
    SequenceStatsAccumulator,
    read_selected_fasta_sequences,
    filter_fasta_by_ids,
)
from genomic_benchmarks_qc.utils.naming import (
    HTML_REPORT_FILE,
    MMSEQS_DIR,
    PLOTS_DIR,
    SIMPLE_REPORT_FILE,
    SPLIT_SUBDIR,
    TMP_PREFIX,
    column_dirname,
    comparison_dirname,
    strip_extensions,
    unique_slugs,
)
from genomic_benchmarks_qc.utils.split_stats import (
    get_basic_stats_from_aggregates,
    get_threshold_stats,
    flag_split_data_leakage,
)


def add_alignment_sequences(results_filt, test_fasta_path, train_fasta_path):
    """Attach the aligned query and target sequences to each MMseqs2 hit.

    The search output identifies sequences by the internal 'seq_<n>_<half>' ids
    assigned while staging, so the sequences themselves have to be read back
    from the staged FASTA files before the HTML report can render alignments.

    Raises RuntimeError if any hit id is missing from its FASTA file, which
    would silently render a blank alignment.
    """
    if results_filt.empty:
        return results_filt

    query_ids = set(results_filt['query'])
    target_ids = set(results_filt['target'])
    test_seq_by_id = read_selected_fasta_sequences(test_fasta_path, query_ids)
    train_seq_by_id = read_selected_fasta_sequences(train_fasta_path, target_ids)

    results_filt = results_filt.copy()
    results_filt['qseq'] = results_filt['query'].map(test_seq_by_id)
    results_filt['tseq'] = results_filt['target'].map(train_seq_by_id)

    missing_qseq = results_filt['qseq'].isna().sum()
    missing_tseq = results_filt['tseq'].isna().sum()
    if missing_qseq or missing_tseq:
        raise RuntimeError(
            "Failed to map MMSeqs2 hit identifiers to input sequences for alignment rendering. "
            f"Missing qseq: {missing_qseq}, missing tseq: {missing_tseq}."
        )

    return results_filt


def _stage_sequences_to_fasta(fasta_path, files, input_format, sequence_column, seq_suffix):
    """Stream one half of the split into a single FASTA, returning its statistics.

    MMseqs2 reads FASTA only, so CSV/TSV inputs are converted here; several
    sequence columns are concatenated into one record, matching the 'merged'
    analysis of the classes command. Records are numbered rather than named, and
    `seq_suffix` marks which half they came from. Sequences are streamed and
    counted incrementally so that inputs larger than memory can be handled.
    """
    acc = SequenceStatsAccumulator()

    with fasta_path.open("w", encoding="utf-8") as fasta_handle:
        for i, sequence in enumerate(stream_files_to_sequences(files, input_format, sequence_column)):
            append_fasta_record(fasta_handle, sequence, f"seq_{i}_{seq_suffix}")
            acc.add(sequence)

    return acc.finalize()


def _build_comparison_dirname(train_files, test_files):
    """Return the directory name for one train-vs-test comparison.

    Train and test files are commonly named identically and told apart only by
    their directory ('train/data.csv', 'test/data.csv'), which would give every
    comparison of a dataset the same name. The parent directories disambiguate.
    """
    train_slug, test_slug = unique_slugs(
        [strip_extensions(train_files[0]), strip_extensions(test_files[0])],
        contexts=[Path(train_files[0]).parent.name, Path(test_files[0]).parent.name],
    )
    return comparison_dirname(train_slug, test_slug)


def _build_simple_report_frame(threshold_stats):
    """Build the one-row leakage verdict written as the simple CSV report."""
    result = {
        "Data Leakage": {
            "Flag": flag_split_data_leakage(threshold_stats['perc_queries_above_thr']),
            "Percentage of leaked queries": f"{threshold_stats['perc_queries_above_thr']:.2f}%",
            "Percentage of leaked targets": f"{threshold_stats['perc_targets_above_thr']:.2f}%",
        }
    }
    df = pd.DataFrame.from_dict(result, orient='index')
    return df


def _write_mmseqs_report_bundle(
    comparison_dir,
    outfile,
    results_filt,
    train_files,
    test_files,
    train_stats,
    test_stats,
    threshold_stats,
    summary,
    train_fasta_path,
    test_fasta_path,
):
    """Generate a comprehensive report bundle for split evaluation results.

    Produces, inside `comparison_dir`:
    - mmseqs/<outfile>: TSV of the MMseqs2 hits above the similarity threshold
    - mmseqs/seq_index_mapping/: FASTA files of only the sequences involved in
      those hits, mapping the internal seq_* ids back to the input sequences
    - plots/: similarity distribution plots
    - report.html: the plots, the leakage summary and the top alignments
    """
    train_filenames = ",".join([Path(f).name for f in train_files])
    test_filenames = ",".join([Path(f).name for f in test_files])

    # Define output paths for all report components
    html_report_path = comparison_dir / HTML_REPORT_FILE
    plots_dir = comparison_dir / PLOTS_DIR
    mmseqs_dir = comparison_dir / MMSEQS_DIR
    seq_index_mapping = mmseqs_dir / 'seq_index_mapping'

    # Export MMseqs2 search results to TSV format
    # Contains hit pairs (test sequences that match training sequences) with similarity scores
    mmseqs_dir.mkdir(parents=True, exist_ok=True)
    export_results_path = mmseqs_dir / outfile
    build_mmseqs_export_frame(results_filt).to_csv(export_results_path, sep="\t", index=False)

    # Create filtered FASTA files containing only sequences involved in hits
    # Maps the internal seq_* identifiers back to original sequences for reference
    seq_index_mapping.mkdir(parents=True, exist_ok=True)
    new_test_fasta_path = seq_index_mapping / 'test_sequences.fasta'
    new_train_fasta_path = seq_index_mapping / 'train_sequences.fasta'
    query_ids = set(results_filt["query"]) if "query" in results_filt.columns else set()
    target_ids = set(results_filt["target"]) if "target" in results_filt.columns else set()
    filter_fasta_by_ids(test_fasta_path, new_test_fasta_path, query_ids)
    filter_fasta_by_ids(train_fasta_path, new_train_fasta_path, target_ids)

    # Aggregate sequence statistics (count, length, GC content, etc.) from train and test sets
    basic_stats = get_basic_stats_from_aggregates(
        train_filenames,
        train_stats,
        test_filenames,
        test_stats,
    )

    # Prepare top 100 hits with full alignment sequences for HTML visualization
    results_filt_for_html = add_alignment_sequences(
        results_filt.head(100),
        test_fasta_path,
        train_fasta_path,
    )

    # Generate interactive HTML report with:
    # - Summary statistics and leakage assessment
    # - Top sequence alignments
    # - Distribution plots comparing train vs test sets
    generate_splits_html_report(
        basic_stats,
        threshold_stats,
        results_filt_for_html,
        html_report_path,
        plots_dir,
        summary["query_similarity_max"],
        summary["target_similarity_max"],
    )

def run(
    train_files,
    test_files,
    format,
    out_folder: Optional[str] = '.',
    sequence_column: Optional[list[str]] = None,
    report_types: Optional[list[str]] = None,
    similarity_threshold: Optional[float] = 90.0,
    threads: Optional[int] = None,
    split_memory_limit: Optional[str] = None,
    keep_tmp_files: Optional[bool] = False,
    log_level: Optional[str] = 'INFO',
    log_file: Optional[str] = None,
):
    """Run the train-test split evaluation.

    This function reads sequences from the provided training and testing files, performs easy-search using MMseqs2, 
    and generates reports about potential data leakage between the training and testing datasets.

    Reports are written into a 'split/<column>/<train>_vs_<test>/'
    sub-directory of `out_folder`, matching the layout of the classes command.
    FASTA inputs have no sequence column and use 'sequence'; several columns are
    searched concatenated and land in 'merged'. Any grouping above that -
    collection, dataset - is left to the caller, who expresses it through
    `out_folder`.

    @param train_files: List of paths to training files.
    @param test_files: List of paths to testing files.
    @param format: Format of the input files (fasta, csv, csv.gz, tsv, tsv.gz).
    @param out_folder: Path to the output folder; reports go into '<out_folder>/split/'. Default: '.'.
    @param sequence_column: Name of the columns with sequences to analyze for datasets in CSV/TSV format.
                            Several columns are concatenated per row and searched together.
                            Default: ['sequence'].
    @param report_types: Types of reports to generate. Default: ['html', 'simple'].
    @param similarity_threshold: Similarity threshold for flagging potential data leakage (between 0 and 100). Default: 90.0.
    @param threads: Maximum number of threads MMseqs2 will use. Default: None.
    @param split_memory_limit: Upper RAM limit for MMseqs2 prefilter structures (e.g., 10G, 1T). Default: None.
    @param keep_tmp_files: Keep temporary files generated for MMseqs2 debugging. Default: False.
    @param log_level: Logging level, default to INFO.
    @param log_file: Path to the log file. If provided, logs will be written to this file as well as to the console.
    @return: None
    """

    if sequence_column is None:
        sequence_column = ['sequence']
    if report_types is None:
        report_types = ['html', 'simple']

    setup_logger(log_level, log_file)
    logging.info("Starting train-test split evaluation.")

    # One directory per train-vs-test comparison, so different comparisons of the
    # same dataset can run concurrently into one output folder without
    # overwriting each other. The scratch directory below lives inside it and is
    # unique per run, so even two runs of the *same* comparison keep their
    # working files apart - though they would still write the same reports.
    comparison_dir = ensure_directory(
        Path(out_folder)
        / SPLIT_SUBDIR
        / column_dirname(format, sequence_column)
        / _build_comparison_dirname(train_files, test_files)
    )

    # `mkdtemp` rather than a fixed name: it creates the directory atomically, so
    # this run is provably its only owner and the cleanup below cannot delete
    # files that belong to anything else. Re-using a fixed name would mean
    # adopting whatever was already there and then removing it - which would
    # take out the scratch files of a concurrent run of the same comparison, or
    # the files a previous `--keep-tmp-files` run was asked to preserve.
    tmp_dir = Path(tempfile.mkdtemp(dir=comparison_dir, prefix=TMP_PREFIX))

    # `mkdtemp` hardcodes mode 0700, which would make kept scratch files
    # unreadable to everyone but the owner even where the reports beside them
    # are group-readable. Match the comparison directory instead, so a shared
    # output tree stays shared.
    try:
        os.chmod(tmp_dir, stat.S_IMODE(os.stat(comparison_dir).st_mode))
    except OSError as mode_error:
        logging.debug(f"Could not match permissions of the temporary directory: {mode_error}")

    train_fasta_path = tmp_dir / 'train_sequences.fasta'
    test_fasta_path = tmp_dir / 'test_sequences.fasta'

    try:
        with log_failures("Train-test split evaluation"):
            train_stats = _stage_sequences_to_fasta(train_fasta_path, train_files, format, sequence_column, 'train')
            test_stats = _stage_sequences_to_fasta(test_fasta_path, test_files, format, sequence_column, 'test')
            num_train_seqs = train_stats["count"]
            num_test_seqs = test_stats["count"]

            if num_train_seqs == 0 or num_test_seqs == 0:
                raise ValueError(
                    "Both train and test inputs must contain at least one sequence before running split evaluation."
                )

            logging.info(f"Read {num_train_seqs} sequences from training files.")
            logging.info(f"Read {num_test_seqs} sequences from testing files.")
        
            # Run MMseqs2 search and keep raw TSV path for chunked post-processing.
            outfile = "mmseqs2_search_result.tsv"
            results_path = mmseqs_runtime.run_search(
                test_fasta_path,
                train_fasta_path,
                tmp_dir / outfile,
                tmp_dir,
                threads=threads,
                split_memory_limit=split_memory_limit,
            )

            summary = summarize_mmseqs_output(results_path, similarity_threshold)
            results_filt = summary["results_filt"]

            # Get threshold stats
            threshold_stats = get_threshold_stats(
                summary = summary,
                similarity_threshold = similarity_threshold, 
                num_train_seqs = num_train_seqs,
                num_test_seqs = num_test_seqs
            )

            if 'simple' in report_types:
                df = _build_simple_report_frame(threshold_stats)
                generate_simple_report(df, comparison_dir / SIMPLE_REPORT_FILE)

            if 'html' in report_types:
                _write_mmseqs_report_bundle(
                    comparison_dir,
                    outfile,
                    results_filt,
                    train_files,
                    test_files,
                    train_stats,
                    test_stats,
                    threshold_stats,
                    summary,
                    train_fasta_path,
                    test_fasta_path,
                )

        logging.info("Train-test split evaluation successfully completed.")
    finally:
        if keep_tmp_files:
            logging.info(f"Keeping temporary files for debugging at: {tmp_dir}")
        else:
            logging.debug("Removing temporary files.")
            try:
                shutil.rmtree(tmp_dir)
            except Exception as cleanup_error:
                logging.warning(f"Failed to remove temporary directory: {cleanup_error}")
