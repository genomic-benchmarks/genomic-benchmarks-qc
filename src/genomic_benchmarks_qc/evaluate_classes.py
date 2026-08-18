"""Compare sequence characteristics between the classes of one dataset.

Reads the input as a set of classes - one FASTA file per class, or one CSV/TSV
label column - computes per-class statistics, and reports on every pair of
classes: a feature that separates two classes tells a model which class a
sequence belongs to without it having to learn anything about the biology.

`run` is the entry point; the CLI is a thin wrapper around it. The report layout
is defined in `genomic_benchmarks_qc.utils.naming`.
"""

import logging
from pathlib import Path
from itertools import combinations
from typing import Optional
import pandas as pd

from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics
from genomic_benchmarks_qc.utils.testing import flag_significant_differences
from genomic_benchmarks_qc.report.report_generator import generate_json_report, generate_simple_report, generate_dataset_html_report
from genomic_benchmarks_qc.utils.input_utils import (
    ensure_directory,
    log_failures,
    read_fasta,
    read_sequences_from_df,
    read_csv_file,
    setup_logger,
)
from genomic_benchmarks_qc.utils.naming import (
    CLASS_SUBDIR,
    DEFAULT_COLUMN_DIR,
    HTML_REPORT_FILE,
    PER_CLASS_DIR,
    PLOTS_DIR,
    SIMPLE_REPORT_FILE,
    comparison_dirname,
    per_column_dirnames,
    strip_extensions,
    unique_slugs,
)


def run_analysis(input_statistics, report_dir, report_types, plot_type):
    """Analyse each class, then every pair of classes, writing reports under `report_dir`.

    Layout produced, one directory per comparison so that every report type has
    a fixed, predictable name inside it:

        <report_dir>/
            per-class/<class>.json
            <classA>_vs_<classB>/
                report.csv
                report.html
                duplicates.txt
                plots/

    Classes are compared in the order they arrive, which `run` has already
    sorted by path name.
    """
    if report_types is None:
        report_types = ['html', 'simple']

    report_dir = Path(report_dir)

    # run individual analysis
    for s in input_statistics:
        stats, _ = s.compute()

        if "json" in report_types:
            per_class_dir = report_dir / PER_CLASS_DIR
            per_class_dir.mkdir(parents=True, exist_ok=True)
            generate_json_report(stats, per_class_dir / f'{s.slug}.json')

    if len(input_statistics) < 2:
        return

    # run pair comparison analysis with all combinations
    for stat1, stat2 in combinations(input_statistics, 2):
        comparison_dir = report_dir / comparison_dirname(stat1.slug, stat2.slug)

        if stat1.seq_column is not None:
            logging.info(f"Comparing classes for sequence column: {stat1.seq_column}")
        logging.info(f"Comparing classes: {stat1.label} vs {stat2.label}")

        logging.debug(f"Running significant differences analysis for {comparison_dir}.")
        results, failed_by_feature = flag_significant_differences(
            stat1, stat2
        )

        # Only create the directory for report types that were actually asked for,
        # so a json-only run does not leave empty comparison directories behind.
        if 'simple' in report_types or 'html' in report_types:
            comparison_dir.mkdir(parents=True, exist_ok=True)

        if 'simple' in report_types:
            generate_simple_report(results, comparison_dir / SIMPLE_REPORT_FILE)

        if 'html' in report_types:
            # Convert results dict to DataFrame using all available result fields
            results_df = pd.DataFrame.from_dict(results, orient='index')
            generate_dataset_html_report(
                stat1, stat2,
                comparison_dir / HTML_REPORT_FILE,
                plots_path=comparison_dir / PLOTS_DIR,
                end_position=min(stat1.end_position, stat2.end_position),
                plot_type=plot_type,
                results=results_df,
                failed_by_feature=failed_by_feature
            )


def _regression_labels(df, label_column):
    """Split a numeric target at its median into 'high' and 'low' classes.

    Returns the dataframe with `label_column` replaced by the two class names,
    and the class names themselves. Rows that are not numeric are dropped.

    Raises ValueError if the target does not yield two classes, which is not a
    recoverable condition: there is nothing to compare, and continuing would
    fail later with an unrelated "no sequences found for label 'low'".
    """
    df[label_column] = pd.to_numeric(df[label_column], errors='coerce')

    nan_count = df[label_column].isna().sum()
    if nan_count > 0:
        logging.warning(f"Dropped {nan_count} rows with non-numeric values in '{label_column}'.")
        # Copied, not just filtered: the classes are written back into the column
        # below, and assigning into a slice of the original frame only warns.
        df = df.dropna(subset=[label_column]).copy()

    if len(df) == 0:
        raise ValueError(
            f"Regression target '{label_column}' contains no numeric values, so the dataset "
            f"cannot be split into high and low classes."
        )

    threshold = df[label_column].median()
    logging.debug(f"Inferred threshold for regression: {threshold}")
    df[label_column] = df[label_column].apply(lambda x: 'high' if x >= threshold else 'low')
    labels = ['high', 'low']

    # Values are split at the median, so a target whose median equals its
    # minimum - a constant or heavily zero-inflated column, or a single
    # remaining row - puts every value in 'high' and leaves 'low' empty.
    missing_classes = [label for label in labels if label not in set(df[label_column])]
    if missing_classes:
        raise ValueError(
            f"Regression target '{label_column}' cannot be split into two classes: all "
            f"{len(df)} value(s) fall on one side of the median ({threshold}), leaving "
            f"'{missing_classes[0]}' empty. Check the target for constant or heavily "
            f"zero-inflated values."
        )

    return df, labels


def _resolve_labels(df, label_column, label_list, regression):
    """Work out which classes to compare, returning (dataframe, labels).

    The dataframe is returned because a regression target both rewrites and
    filters it; for the other paths it comes back unchanged.
    """
    if regression:
        # Regression derives its own classes, so anything given in --label-list
        # is silently dropped; say so rather than reporting on classes the user
        # did not ask for.
        if not (len(label_list) == 1 and label_list[0] == 'infer'):
            logging.warning(
                f"Ignoring label_list {label_list}: regression splits '{label_column}' "
                f"into 'high' and 'low' classes."
            )
        return _regression_labels(df, label_column)

    if len(label_list) == 1 and label_list[0] == 'infer':
        labels = sorted(df[label_column].unique().tolist())
        if not labels:
            raise ValueError(
                f"No usable labels found in column '{label_column}'."
            )
        logging.debug(f"Inferred labels: {labels}")
        return df, labels

    # dict.fromkeys drops repeated labels; a repeated label would otherwise be
    # compared against itself. The order given here does not survive - classes
    # are sorted by path name by the caller.
    return df, list(dict.fromkeys(str(label) for label in label_list))


def _evaluate_fasta_classes(input, out_folder, report_types, plot_type, end_position):
    """Evaluate classes given as one FASTA file per class."""
    # Each file is one class, named after its stem. Stems can repeat across
    # directories, so the parent directory disambiguates the report paths.
    labels = [strip_extensions(f) for f in input]
    slugs = unique_slugs(labels, contexts=[Path(f).parent.name for f in input])

    # Classes are ordered by their path name, so the comparison directories
    # do not depend on the order the files happened to be given in, and match
    # what a CSV input with the same class names would produce.
    ordered_inputs = sorted(zip(input, labels, slugs), key=lambda item: item[2])

    seq_stats = []
    for input_file, label, slug in ordered_inputs:
        sequences = read_fasta(input_file)
        # A class with no sequences has nothing to compare, and is almost always
        # a wrong path or a truncated file rather than something to report on.
        # Rejected here so the run fails on the input, not partway through the
        # statistics - as the CSV path does for a label with no rows.
        if not sequences:
            logging.error(f"No sequences found in FASTA file '{input_file}'.")
            raise ValueError(f"No sequences found in FASTA file '{input_file}'.")
        logging.debug(f"Read {len(sequences)} sequences from FASTA file {input_file}.")
        seq_stats += [SequenceStatistics(sequences, filename=Path(input_file).name, filepath=input_file,
                                         label=label, slug=slug, end_position=end_position)]

    run_analysis(
        input_statistics=seq_stats,
        report_dir=Path(out_folder) / CLASS_SUBDIR / DEFAULT_COLUMN_DIR,
        report_types=report_types,
        plot_type=plot_type
    )


def _evaluate_table_classes(input, format, out_folder, sequence_column, label_column, label_list,
                            regression, report_types, plot_type, end_position):
    """Evaluate classes given as a label column of one or more CSV/TSV files."""
    # read all files into one dataframe (single file wrapped in list)
    if len(input) > 1:
        logging.info(f"Merging {len(input)} input files: {', '.join(input)}")
    dfs = [read_csv_file(f, format, sequence_column, label_column) for f in input]
    df = pd.concat(dfs, ignore_index=True)
    logging.debug(f"Read {len(df)} rows from {len(dfs)} file(s)")

    df, labels = _resolve_labels(df, label_column, label_list, regression)

    # single source filename for reports
    # Shown in the report, not used as a path - hence a literal rather than the
    # directory name that happens to read the same.
    filename = Path(input[0]).name if len(input) == 1 else 'merged'
    filepath = input[0] if len(input) == 1 else ", ".join(input)

    # Labels come straight from the data, so they may contain characters
    # that are unusable in a path ('/', spaces) or that differ only in case.
    label_slugs = dict(zip(labels, unique_slugs(labels)))

    # Same ordering rule as for FASTA inputs, applied whether the labels were
    # inferred, derived from a regression target or listed explicitly.
    labels = sorted(labels, key=label_slugs.__getitem__)

    # One directory per sequence column, plus a trailing one for the merged
    # analysis when there is more than one column.
    column_dirs, merged_dir = per_column_dirnames(sequence_column)
    class_dir = Path(out_folder) / CLASS_SUBDIR

    def statistics_for(seq_col, reported_column):
        seq_stats = []
        for label in labels:
            sequences = read_sequences_from_df(df, seq_col, label_column, label)
            logging.debug(f"Read {len(sequences)} sequences for label '{label}' from column '{reported_column}'.")
            seq_stats += [SequenceStatistics(sequences, filename=filename, filepath=filepath, label=label,
                                             slug=label_slugs[label], seq_column=reported_column,
                                             end_position=end_position)]
        return seq_stats

    # loop over individual sequence columns
    for seq_col, column_dir in zip(sequence_column, column_dirs):
        run_analysis(
            input_statistics=statistics_for(seq_col, seq_col),
            report_dir=class_dir / column_dir,
            report_types=report_types,
            plot_type=plot_type,
        )

    # if multiple sequence columns, also evaluate merged sequences
    if merged_dir is not None:
        run_analysis(
            input_statistics=statistics_for(sequence_column, '_'.join(sequence_column)),
            report_dir=class_dir / merged_dir,
            report_types=report_types,
            plot_type=plot_type,
        )


def run(input, 
        format, 
        out_folder='.', 
        sequence_column: Optional[list[str]] = None, 
        label_column='label', 
        label_list: Optional[list[str]] = None,
        regression: Optional[bool] = False,
        report_types: Optional[list[str]] = None,
        end_position: Optional[int] = None,
        plot_type: Optional[str] = 'boxen',
        log_level: Optional[str] = 'INFO',
        log_file: Optional[str] = None
    ):
    """Run the dataset evaluation.

    This function reads sequences from the provided input files, performs analysis, and generates reports about the sequences.

    Reports are written to '<out_folder>/class/<column>/<classA>_vs_<classB>/',
    one directory per compared pair of classes. FASTA inputs have no sequence
    column and use 'sequence', so the layout is the same for every input format.
    Any grouping above that - collection, dataset, split - is left to the caller,
    who expresses it through `out_folder`.

    @param input: List of paths to input files. Can be a list of files, each containing sequences from one class.
    @param format: Format of the input files (fasta, csv, csv.gz, tsv, tsv.gz).
    @param out_folder: Path to the output folder; reports go into '<out_folder>/class/'. Default: '.'.
    @param sequence_column: Name of the columns with sequences to analyze for datasets in CSV/TSV format. 
                            Either one column or list of columns. Each column is analyzed separately, and
                            all of them together in an extra 'merged' report. Default: ['sequence']
    @param label_column: Name of the label column for datasets in CSV/TSV format. Default: 'label'.
    @param label_list: List of label classes to consider or "infer" to parse different labels automatically from label column.
                      For datasets in CSV/TSV format. Default: ['infer'].
    @param regression: If True, label column is considered as a regression target and values are split into 2 classes
                       at the median. Raises ValueError if that does not produce two non-empty classes.
    @param report_types: Types of reports to generate. Default: ['html', 'simple'].
    @param end_position: End position of the sequences to consider in per position statistics. 
                         If not provided, 75th percentile of sequence lengths will be used. Default: None.
    @param plot_type: Type of plot to use for visualizations. For bigger datasets, "boxen" is recommended. Default: 'boxen'.
    @param log_level: Logging level, default to INFO.
    @param log_file: Path to the log file. If provided, logs will be written to this file as well as to the console.
    @return: None
    """

    if sequence_column is None:
        sequence_column = ['sequence']
    if report_types is None:
        report_types = ['html', 'simple']
    if label_list is None:
        label_list = ['infer']

    setup_logger(log_level, log_file)
    logging.info("Starting classes evaluation.")

    ensure_directory(out_folder)

    with log_failures("Classes evaluation"):
        # we have multiple fasta files with one label each
        if format.startswith('fa'):
            _evaluate_fasta_classes(input, out_folder, report_types, plot_type, end_position)
        # we have CSV/TSV
        else:
            _evaluate_table_classes(input, format, out_folder, sequence_column, label_column,
                                    label_list, regression, report_types, plot_type, end_position)

    logging.info("Classes evaluation successfully completed.")
