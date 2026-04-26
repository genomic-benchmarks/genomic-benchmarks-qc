import logging
from pathlib import Path
from typing import Optional
import shutil
import pandas as pd

from genbenchQC.report.report_generator import generate_splits_html_report, generate_simple_report
from genbenchQC.utils.mmseqs_summary import summarize_mmseqs_output
from genbenchQC.utils import mmseqs_runtime
from genbenchQC.utils.input_utils import (
    setup_logger,
    stream_files_to_sequences,
    append_fasta_record,
    init_sequence_stats,
    update_sequence_stats,
    finalize_sequence_stats,
    read_selected_fasta_sequences,
)
from genbenchQC.utils.split_stats import (
    get_basic_stats_from_aggregates,
    get_threshold_stats,
    filter_fasta_by_ids,
)

def add_alignment_sequences(results_filt, test_fasta_path, train_fasta_path):
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
    sequence_stats_acc = init_sequence_stats()

    with fasta_path.open("w", encoding="utf-8") as fasta_handle:
        for i, sequence in enumerate(stream_files_to_sequences(files, input_format, sequence_column)):
            append_fasta_record(fasta_handle, sequence, f"seq_{i}_{seq_suffix}")
            update_sequence_stats(sequence_stats_acc, sequence)

    return finalize_sequence_stats(sequence_stats_acc)


def _build_split_report_stem(train_files, test_files):
    train_stem = Path(train_files[0]).name.replace("".join(Path(train_files[0]).suffixes), "")
    test_stem = Path(test_files[0]).name.replace("".join(Path(test_files[0]).suffixes), "")
    return f"split_check_{train_stem}_vs_{test_stem}"


def _build_simple_report_frame(threshold_stats, has_leakage):
    result = {
        "Data Leakage": {
            "Flag": "Fail" if has_leakage else "Pass",
            "Percentage of leaked queries": f"{threshold_stats['perc_queries_above_thr']:.2f}%",
            "Percentage of leaked targets": f"{threshold_stats['perc_targets_above_thr']:.2f}%",
        }
    }
    df = pd.DataFrame.from_dict(result, orient='index')
    df.index.name = "Statistic"
    return df


def _write_mmseqs_report_bundle(
    out_folder,
    report_stem,
    results_path,
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
    train_filenames = ",".join([Path(f).name for f in train_files])
    test_filenames = ",".join([Path(f).name for f in test_files])

    html_report_path = out_folder / f"{report_stem}_report.html"
    plots_dir = out_folder / f"{report_stem}_plots"
    mmseqs_dir = out_folder / f"{report_stem}_mmseqs"
    seq_index_mapping = mmseqs_dir / 'seq_index_mapping'

    mmseqs_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(results_path, mmseqs_dir / outfile)

    seq_index_mapping.mkdir(parents=True, exist_ok=True)
    new_test_fasta_path = seq_index_mapping / 'test_sequences.fasta'
    new_train_fasta_path = seq_index_mapping / 'train_sequences.fasta'
    query_ids = set(results_filt["query"]) if "query" in results_filt.columns else set()
    target_ids = set(results_filt["target"]) if "target" in results_filt.columns else set()
    filter_fasta_by_ids(test_fasta_path, new_test_fasta_path, query_ids)
    filter_fasta_by_ids(train_fasta_path, new_train_fasta_path, target_ids)

    basic_stats = get_basic_stats_from_aggregates(
        train_filenames,
        train_stats,
        test_filenames,
        test_stats,
    )

    results_filt_for_html = add_alignment_sequences(
        results_filt.head(100),
        test_fasta_path,
        train_fasta_path,
    )

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
    similarity_threshold: Optional[float] = 80.0,
    threads: Optional[int] = None,
    split_memory_limit: Optional[str] = None,
    keep_tmp_files: Optional[bool] = False,
    log_level: Optional[str] = 'INFO',
    log_file: Optional[str] = None,
):
    """Run the train-test split evaluation.

    This function reads sequences from the provided training and testing files, performs easy-search using MMseqs2, 
    and generates reports about potential data leakage between the training and testing datasets.

    @param train_files: List of paths to training files.
    @param test_files: List of paths to testing files.
    @param format: Format of the input files (fasta, csv, csv.gz, tsv, tsv.gz).
    @param out_folder: Path to the output folder. Default: '.'.
    @param sequence_column: Name of the columns with sequences to analyze for datasets in CSV/TSV format. 
                            Default: ['sequence'].
    @param report_types: Types of reports to generate. Default: ['html', 'simple'].
    @param similarity_threshold: Similarity threshold for flagging potential data leakage (between 0 and 100). Default: 80.0.
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

    out_folder = Path(out_folder)
    if not out_folder.exists():
        logging.info(f"Output folder {out_folder} does not exist. Creating it.")
        out_folder.mkdir(parents=True, exist_ok=True)

    tmp_dir =  out_folder / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    train_fasta_path = tmp_dir / 'train_sequences.fasta'
    test_fasta_path = tmp_dir / 'test_sequences.fasta'

    try:
        train_stats = _stage_sequences_to_fasta(train_fasta_path, train_files, format, sequence_column, 'train')
        test_stats = _stage_sequences_to_fasta(test_fasta_path, test_files, format, sequence_column, 'test')
        num_train_seqs = train_stats["count"]
        num_test_seqs = test_stats["count"]

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

        report_stem = _build_split_report_stem(train_files, test_files)

        # Get threshold stats
        threshold_stats = get_threshold_stats(
            summary = summary,
            similarity_threshold = similarity_threshold, 
            num_train_seqs = num_train_seqs,
            num_test_seqs = num_test_seqs
        )

        if 'simple' in report_types:
            simple_report_path = out_folder / f"{report_stem}.csv"
            df = _build_simple_report_frame(threshold_stats, not results_filt.empty)
            generate_simple_report(df, simple_report_path)
        
        if 'html' in report_types:
            _write_mmseqs_report_bundle(
                out_folder,
                report_stem,
                results_path,
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
    except Exception as exc:
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.exception("Train-test split evaluation failed.")
        else:
            logging.error(f"Train-test split evaluation failed: {exc}")
        raise
    finally:
        if keep_tmp_files:
            logging.info(f"Keeping temporary files for debugging at: {tmp_dir}")
        else:
            logging.debug("Removing temporary files.")
            try:
                shutil.rmtree(tmp_dir)
            except Exception as cleanup_error:
                logging.warning(f"Failed to remove temporary directory: {cleanup_error}")
