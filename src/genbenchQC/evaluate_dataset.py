import argparse
import logging
from pathlib import Path
from itertools import combinations
from typing import Optional
import pandas as pd

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.utils.testing import flag_significant_differences
from genbenchQC.report.report_generator import generate_json_report, generate_sequence_html_report, generate_simple_report, generate_dataset_html_report
from genbenchQC.utils.input_utils import read_fasta, read_sequences_from_df, read_multisequence_df, read_csv_file, setup_logger

def run_analysis(input_statistics, out_folder, report_types, seq_report_types, plot_type):
   
    """
    Generate sequence-level and dataset-level reports and plots for the given SequenceStatistics objects.
    
    Creates per-sequence reports (JSON and/or HTML with plots) for each statistics object when requested, and for every pair of statistics objects generates comparison reports (CSV "simple" and/or HTML with plots). All output files and plot directories are written under out_folder.
    
    Parameters:
        input_statistics (Iterable[SequenceStatistics]): SequenceStatistics objects to analyze. Each object must provide computed statistics and an end_position via its compute() method, and expose filename, seq_column, label, and end_position attributes.
        out_folder (str | pathlib.Path): Destination directory for reports and plot directories; will be created if it does not exist.
        report_types (Iterable[str]): Dataset-level report types to produce (e.g., 'simple', 'html').
        seq_report_types (Iterable[str]): Sequence-level report types to produce (e.g., 'json', 'html').
        plot_type (str | None): Plot style/type to use when generating HTML reports and plots; passed through to report generators.
    """
    out_folder = Path(out_folder)

    # run individual analysis
    for s in input_statistics:
        stats, end_position = s.compute()

        if seq_report_types:

            filename = Path(s.filename).stem
            if s.seq_column is not None:
                filename += f'_{s.seq_column}'
            if s.label is not None:
                filename += f'_{s.label}'

            if 'json' in seq_report_types:
                json_report_path = out_folder / Path(filename + '_report.json')
                generate_json_report(stats, json_report_path)

            if 'html' in seq_report_types:
                html_report_path = out_folder / Path(filename + '_report.html')
                plots_path = out_folder / Path(filename + '_plots')
                generate_sequence_html_report(stats, html_report_path, plots_path, end_position, plot_type)

    if len(input_statistics) < 2:
        return

    # run pair comparison analysis with all combinations
    for stat1, stat2 in combinations(input_statistics, 2):
        filename = "dataset_report"
        if stat1.seq_column is not None:
            filename += f'_{stat1.seq_column}'
        if stat1.label is not None and stat2.label is not None:
            filename += f'_label_{stat1.label}_vs_{stat2.label}'
            logging.debug(f"Comparing datasets label: {stat1.label} vs {stat2.label}")
        else:
            filename += f'_{Path(stat1.filename).stem}_{Path(stat2.filename).stem}'
            logging.debug(
                f"Comparing datasets: {stat1.filename} vs {stat2.filename}")

        results = flag_significant_differences(
            stat1, stat2, 
        )
        
        if 'simple' in report_types:
            simple_report_path = out_folder / Path(f'{filename}.csv')
            generate_simple_report(results, simple_report_path)

        if 'html' in report_types:
            html_report_path = out_folder / Path(f'{filename}.html')
            plots_path = out_folder / Path(f'{filename}_plots')
            generate_dataset_html_report(
                stat1, stat2, 
                html_report_path, 
                plots_path=plots_path, 
                end_position=min(stat1.end_position, stat2.end_position),
                plot_type=plot_type,
                results=results, threshold=None
            )

def run(input, 
        format, 
        out_folder='.', 
        sequence_column: Optional[list[str]] = ['sequences'], 
        label_column='label', 
        label_list: Optional[list[str]] = ['infer'],
        regression: Optional[bool] = False,
        report_types: Optional[list[str]] = ['html', 'simple'],
        seq_report_types: Optional[list[str]] = None,
        end_position: Optional[int] = None,
        plot_type: Optional[str] = 'boxen',
        flag_threshold: Optional[float] = 0.015,
        log_level: Optional[str] = 'INFO',
        log_file: Optional[str] = None
    ):
    """
        Orchestrates reading sequence input(s), building SequenceStatistics for groups or labels, running analyses, and writing the requested reports to the output folder.
        
        Parameters:
            input (list[str]): Paths to one or more input files. For FASTA mode, each file is treated as a separate label; for tabular formats a single file may contain multiple labels.
            format (str): Input file format: 'fasta', 'csv', 'csv.gz', 'tsv', or 'tsv.gz'.
            out_folder (str): Destination folder for generated reports and plots.
            sequence_column (list[str] | str): Column name or list of column names in tabular inputs that contain sequences to analyze.
            label_column (str): Column name in tabular inputs that contains label values.
            label_list (list[str] | ['infer']): Explicit list of labels to analyze, or ['infer'] to derive labels from the label column.
            regression (bool): If True, treat the label column as a numeric regression target and split it into 'high'/'low' classes using the median.
            report_types (list[str]): Dataset-level report formats to generate (e.g., 'html', 'simple').
            seq_report_types (list[str] | None): Sequence-group-level report formats to generate for each label/group.
            end_position (int | None): Maximum sequence position to include in per-position statistics; if None a length-based heuristic is used.
            plot_type (str): Plot style used for visualizations (e.g., 'boxen').
            log_level (str): Logging level passed to the logger.
            log_file (str | None): Optional path to a log file to also write logs.
        
        """

    setup_logger(log_level, log_file)
    logging.info("Starting dataset evaluation.")

    if not Path(out_folder).exists():
        logging.info(f"Output folder {out_folder} does not exist. Creating it.")
        Path(out_folder).mkdir(parents=True, exist_ok=True)

    # we have multiple fasta files with one label each
    if format == 'fasta':
        seq_stats = []
        for input_file in input:
            sequences = read_fasta(input_file)
            logging.debug(f"Read {len(sequences)} sequences from FASTA file {input_file}.")
            seq_stats += [SequenceStatistics(sequences, filename=Path(input_file).name, 
                                             label=Path(input_file).stem, end_position=end_position)]
        run_analysis(
            input_statistics = seq_stats, 
            out_folder = out_folder, 
            report_types = report_types, 
            seq_report_types = seq_report_types, 
            plot_type = plot_type
        )

    # we have CSV/TSV
    else:
        # we have one file with multiple labels or regression target
        if len(input) == 1:
            df = read_csv_file(input[0], format, sequence_column, label_column)

            # if regression is True, we split the label column into two classes
            if regression:
                # convert the label column to numeric if it is not already
                if not pd.api.types.is_numeric_dtype(df[label_column]):
                    logging.debug(f"Converting label column '{label_column}' to numeric type for regression.")
                    df[label_column] = pd.to_numeric(df[label_column], errors='coerce')
                # infer the threshold as the median of the label column
                threshold = df[label_column].median()
                logging.debug(f"Inferred threshold for regression: {threshold}")
                df[label_column] = df[label_column].apply(lambda x: 'high' if x >= threshold else 'low')
                labels = ['high', 'low']

            # get the list of labels to consider
            elif len(label_list) == 1 and label_list[0] == 'infer':
                labels = df[label_column].unique().tolist()
                logging.debug(f"Inferred labels: {labels}")
            else:
                labels = [str(label) for label in label_list]

            # loop over sequences with specific label and run statistics
            for seq_col in sequence_column:
                seq_stats = []
                for label in labels:
                    sequences = read_sequences_from_df(df, seq_col, label_column, label)
                    logging.debug(f"Read {len(sequences)} sequences for label '{label}' from column '{seq_col}'.")
                    seq_stats += [SequenceStatistics(sequences, filename=Path(input[0]).name, label=label, 
                                                     seq_column=seq_col, end_position=end_position)]
                run_analysis(
                    input_statistics = seq_stats, 
                    out_folder = out_folder, 
                    report_types = report_types, 
                    seq_report_types = seq_report_types, 
                    plot_type = plot_type, 
                )

            # handle multiple sequence columns by concatenating sequences and running statistics on them
            if len(sequence_column) > 1:
                seq_stats = []
                for label in labels:
                    sequences = read_multisequence_df(df, sequence_column, label_column, label)
                    seq_stats += [SequenceStatistics(sequences, filename=Path(input[0]).name, label=label,
                                                     seq_column='_'.join(sequence_column))]
                run_analysis(
                    input_statistics = seq_stats, 
                    out_folder = out_folder, 
                    report_types = report_types, 
                    seq_report_types = seq_report_types, 
                    plot_type = plot_type, 
                )

        # we have multiple files with one label each
        else:
            # run statistics across input files
            for seq_col in sequence_column:
                seq_stats = []
                for input_file in input:
                    sequences = read_sequences_from_df(read_csv_file(input_file, format, seq_col), seq_col)
                    logging.debug(f"Read {len(sequences)} sequences from file {input_file} in column '{seq_col}'.")
                    seq_stats += [SequenceStatistics(sequences, filename=Path(input_file).name, 
                                                     label=Path(input_file).stem, seq_column=seq_col,
                                                     end_position=end_position)]
                run_analysis(
                    input_statistics = seq_stats, 
                    out_folder = out_folder, 
                    report_types = report_types, 
                    seq_report_types = seq_report_types, 
                    plot_type = plot_type, 
                )

            # handle multiple sequence columns
            if len(sequence_column) > 1:
                seq_stats = []
                for input_file in input:
                    sequences = read_multisequence_df(read_csv_file(input_file, format, sequence_column), sequence_column)
                    seq_stats += [SequenceStatistics(sequences, filename=Path(input_file).name, label=Path(input_file).stem,
                                                     seq_column='_'.join(sequence_column), end_position=end_position)]
                run_analysis(
                    input_statistics = seq_stats, 
                    out_folder = out_folder, 
                    report_types = report_types, 
                    seq_report_types = seq_report_types, 
                    plot_type = plot_type, 
                )

    logging.info("Dataset evaluation successfully completed.")


def parse_args():
    """
    Parse command-line arguments for the dataset evaluation tool.
    
    Performs validation (for example, when format is 'fasta' at least two input files are required) and returns the populated argument namespace.
    
    Returns:
        argparse.Namespace: Parsed CLI options with fields:
            - input: list[str] paths to input files
            - format: input file format ('fasta', 'csv', 'csv.gz', 'tsv', 'tsv.gz')
            - sequence_column: list[str] sequence column name(s)
            - label_column: str label column name
            - label_list: list[str] labels to consider or ['infer']
            - regression: bool whether to treat label column as regression target
            - out_folder: str output directory
            - report_types: list[str] dataset report types to generate
            - seq_report_types: list[str] per-sequence-group report types to generate
            - end_position: int|None end position to consider for per-position stats
            - plot_type: str plot style ('boxen' or 'violin')
            - log_level: str logging level
            - log_file: str|None path to optional log file
    """
    parser = argparse.ArgumentParser(description='A tool for evaluating sequence datasets.')
    parser.add_argument('--input', type=str, help='Path to the dataset file. '
                                                  'Can be a list of files, each containing sequences from one class.', nargs='+', required=True)
    parser.add_argument('--format', help="Format of the input files.", choices=['fasta', 'csv', 'csv.gz', 'tsv', 'tsv.gz'], required=True) # potentially add HF support
    parser.add_argument('--sequence_column', type=str, help='Name of the columns with sequences to analyze for datasets in CSV/TSV format. '
                                                            'Either one column or list of columns.', nargs='+', default=['sequence'])
    parser.add_argument('--label_column', type=str, help='Name with the label column for datasets in CSV/TSV format.', default='label')
    parser.add_argument('--label_list', type=str, nargs='+', help='List of label classes to consider or "infer" to parse different labels automatically from label column.'
                                                       ' For datasets in CSV/TSV format.', default=['infer'])
    parser.add_argument('--regression', action='store_true', help='If True, label column is considered as a regression target and values are split into 2 classes')
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    parser.add_argument('--report_types', type=str, nargs='+', choices=['json', 'html', 'simple'], default=['html', 'simple'],
                        help='Types of reports to generate. Default: [html, simple].')
    parser.add_argument('--seq_report_types', type=str, nargs='+', choices=['json', 'html'], default=[],
                        help='Types of reports to generate for individual groups of sequences. Default: [].')
    parser.add_argument('--end_position', type=int, default=None,
                        help='End position of the sequences to consider in per position statistics. If not provided, 75th percentile of sequence lengths will be used.')
    parser.add_argument('--plot_type', type=str, help='Type of plot to use for visualizations. For bigger datasets, "boxen" in recommended. Default: boxen.',
                        choices=['boxen', 'violin'], default='boxen')
    parser.add_argument('--log_level', type=str, help='Logging level, default to INFO.', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO')
    parser.add_argument('--log_file', type=str, help='Path to the log file. If provided, logs will be written to this file as well as to the console.', default=None)
    args = parser.parse_args()

    if args.format == 'fasta' and len(args.input) < 2:
        parser.error("When format is 'fasta', the input must contain individual files for each class.")

    return args

def main():
    """
    Entry point for the command-line tool that parses CLI arguments and invokes `run` with those arguments.
    
    This function reads command-line options, validates them via `parse_args()`, and dispatches execution to `run()` using the parsed values.
    """
    args = parse_args()
    run(input = args.input, 
        format = args.format, 
        out_folder = args.out_folder, 
        sequence_column = args.sequence_column, 
        label_column = args.label_column, 
        label_list = args.label_list,
        regression = args.regression,
        report_types = args.report_types,
        seq_report_types = args.seq_report_types,
        end_position = args.end_position,
        plot_type = args.plot_type,
        log_level = args.log_level,
        log_file = args.log_file
    )

if __name__ == '__main__':
    main()