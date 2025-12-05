import argparse
from pathlib import Path
from typing import Optional
import logging

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.report.report_generator import generate_json_report, generate_sequence_html_report
from genbenchQC.utils.input_utils import read_fasta, read_sequences_from_df, read_multisequence_df, read_csv_file, setup_logger

def run_analysis(seq_stats, out_folder, report_types, plot_type):

    if not Path(out_folder).exists():
        logging.info(f"Output folder {out_folder} does not exist. Creating it.")
        Path(out_folder).mkdir(parents=True, exist_ok=True)

    stats, end_position = seq_stats.compute()

    filename = Path(seq_stats.filename).stem
    if seq_stats.seq_column is not None:
        filename += f'_{seq_stats.seq_column}'
    if seq_stats.label is not None:
        filename += f'_{seq_stats.label}'

    if 'json' in report_types:
        json_report_path = Path(out_folder, filename + '_report.json')
        generate_json_report(stats, json_report_path)
    if 'html' in report_types:
        html_report_path = Path(out_folder, filename + '_report.html')
        plots_path = out_folder / Path(filename + '_plots')
        generate_sequence_html_report(
            stats, 
            html_report_path, 
            plots_path, 
            end_position=end_position, 
            plot_type=plot_type
        )

def run(input, format, 
        out_folder: Optional[str] = '.', 
        sequence_column: Optional[list[str]] = ['sequences'], 
        label_column: Optional[str] = None, 
        label: Optional[str] = None,
        report_types: Optional[list[str]] = ['html'],
        end_position: Optional[int] = None,
        plot_type: Optional[str] = 'boxen',
        log_level: Optional[str] = 'INFO',
        log_file: Optional[str] = None
    ):
    """
    Run sequence evaluation on the provided input file.

    This function reads sequences from the input file, performs analysis, and generates reports.

    @param input: Path to the input file containing sequences.
    @param format: Format of the input file (fasta, csv, csv.gz, tsv, tsv.gz).
    @param out_folder: Path to the output folder. Default: '.'.
    @param sequence_column: Name of the columns with sequences to analyze for datasets in CSV/TSV format. 
                            Default: ['sequence'].
    @param label_column: Name of the label column for datasets in CSV/TSV format. Needed only if you want to select a specific class from the dataset.
    @param label: Label of the class to select from the whole dataset. If not specified, the whole dataset is taken and analyzed as one piece.
    @param report_types: Types of reports to generate. Default: ['html'].
    @param end_position: End position of the sequences to plot in the per position plots. 
                         If not provided, 75th percentile of sequence lengths will be used. Default: None.
    @param plot_type: Type of the plot to generate for per sequence nucleotide content. For bigger datasets, "boxen" is recommended. Default: 'boxen'.
    @param log_level: Logging level, default to INFO.
    @param log_file: Path to the log file. If provided, logs will be written to this file as well as to the console.
    @return: None
    """

    setup_logger(log_level, log_file)
    logging.info("Starting sequence evaluation.")

    if format == 'fasta':
        seqs = read_fasta(input)
        logging.debug(f"Read {len(seqs)} sequences from FASTA file.")
        run_analysis(
            SequenceStatistics(seqs, filename=Path(input).name, filepath=input, 
                               label=label, end_position=end_position),
            out_folder, report_types=report_types, plot_type=plot_type
        )
    else:
        df = read_csv_file(input, format, sequence_column, label_column)

        for seq_col in sequence_column:
            sequences = read_sequences_from_df(df, seq_col, label_column, label)
            logging.debug(f"Read {len(sequences)} sequences from CSV/TSV file.")
            run_analysis(
                SequenceStatistics(sequences, filename=Path(input).name, filepath=input,
                                   seq_column=seq_col, label=label, end_position=end_position), 
                out_folder, report_types=report_types, plot_type=plot_type
            )

        if len(sequence_column) > 1:
            sequences = read_multisequence_df(df, sequence_column, label_column, label)
            run_analysis(
                SequenceStatistics(sequences, filename=Path(input).name, 
                                   filepath=input, seq_column='_'.join(sequence_column), 
                                   label=label, end_position=end_position), 
                out_folder, report_types=report_types, plot_type=plot_type
            )

    logging.info("Sequence evaluation successfully completed.")

def parse_args():
    parser = argparse.ArgumentParser(description='A tools for evaluating sequence data.')
    parser.add_argument('--input', type=str, help='Path to the input file.', required=True)
    parser.add_argument('--format', help="Format of the input file.", choices=['fasta', 'csv', 'csv.gz', 'tsv', 'tsv.gz'], required=True)
    parser.add_argument('--sequence_column', type=str,
                        help='Name of the columns with sequences to analyze for datasets in CSV/TSV format. '
                             'Either one column or list of columns.', nargs='+', default=['sequence'])
    parser.add_argument('--label_column', type=str, help='Name with the label column for datasets in CSV/TSV format. Needed only if you want to select a specific class from the dataset.',
                        default=None)
    parser.add_argument('--label', type=str,
                        help='Label of the class to select from the whole dataset. If not specified, the whole dataset is taken and analyzed as one piece.', default=None)
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    parser.add_argument('--report_types', type=str, nargs='+', choices=['json', 'html'],
                        help='Types of reports to generate. Default: [html]', default=['html'])
    parser.add_argument('--end_position', type=int, default=None,
                        help='End position of the sequences to plot in the per position plots. If not provided, 75th percentile of sequence lengths will be used.')
    parser.add_argument('--plot_type', type=str, help='Type of the plot to generate for per sequence nucleotide content. For bigger datasets, "boxen" is recommended. Default: boxen.',
                        choices=['boxen', 'violin'], default='boxen')
    parser.add_argument('--log_level', type=str, help='Logging level, default to INFO.',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO')
    parser.add_argument('--log_file', type=str, help='Path to the log file. If provided, logs will be written to this file as well as to the console.', default=None)

    args = parser.parse_args()

    if (args.label_column is not None and args.label is None) or (args.label_column is None and args.label is not None):
        parser.error("--label_column and --label must be provided together.")

    return args

def main():
    args = parse_args()
    run(input = args.input, 
        format = args.format, 
        out_folder = args.out_folder, 
        sequence_column = args.sequence_column, 
        label_column = args.label_column, 
        label = args.label, 
        report_types = args.report_types,
        end_position = args.end_position,
        plot_type = args.plot_type,
        log_level = args.log_level,
        log_file = args.log_file
    )

if __name__ == '__main__':
    main()