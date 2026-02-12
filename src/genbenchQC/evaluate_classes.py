import logging
from pathlib import Path
from itertools import combinations
from typing import Optional
import pandas as pd

from genbenchQC.utils.statistics import SequenceStatistics
from genbenchQC.utils.testing import flag_significant_differences
from genbenchQC.report.report_generator import generate_json_report, generate_simple_report, generate_dataset_html_report
from genbenchQC.utils.input_utils import read_fasta, read_sequences_from_df, read_multisequence_df, read_csv_file, setup_logger

def run_analysis(input_statistics, out_folder, report_types, plot_type):
   
    if report_types is None:
        report_types = ['html', 'simple']

    out_folder = Path(out_folder)

    # run individual analysis
    for s in input_statistics:
        stats, _ = s.compute()

        if "json" in report_types:

            filename = Path(s.filename).stem
            if s.seq_column is not None:
                filename += f'_{s.seq_column}'
            if s.label is not None:
                filename += f'_{s.label}'

            json_report_path = out_folder / Path(filename + '_report.json')
            generate_json_report(stats, json_report_path)

    if len(input_statistics) < 2:
        return

    # run pair comparison analysis with all combinations
    for stat1, stat2 in combinations(input_statistics, 2):
        filename = "dataset_report"
        if stat1.seq_column is not None:
            filename += f'_{stat1.seq_column}'
            logging.info(f"Comparing classes for sequence column: {stat1.seq_column}")
        if stat1.label is not None and stat2.label is not None:
            filename += f'_label_{stat1.label}_vs_{stat2.label}'
            logging.info(f"Comparing classes: {stat1.label} vs {stat2.label}")
        else:
            filename += f'_{Path(stat1.filename).stem}_{Path(stat2.filename).stem}'
            logging.info(
                f"Comparing classes: {stat1.filename} vs {stat2.filename}")

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
                results=results
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

    @param input: List of paths to input files. Can be a list of files, each containing sequences from one class.
    @param format: Format of the input files (fasta, csv, csv.gz, tsv, tsv.gz).
    @param out_folder: Path to the output folder. Default: '.'.
    @param sequence_column: Name of the columns with sequences to analyze for datasets in CSV/TSV format. 
                            Either one column or list of columns. Default: ['sequence']
    @param label_column: Name of the label column for datasets in CSV/TSV format. Default: 'label'.
    @param label_list: List of label classes to consider or "infer" to parse different labels automatically from label column.
                      For datasets in CSV/TSV format.
    @param regression: If True, label column is considered as a regression target and values are split into 2 classes.
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

    if not Path(out_folder).exists():
        logging.info(f"Output folder {out_folder} does not exist. Creating it.")
        Path(out_folder).mkdir(parents=True, exist_ok=True)

    # we have multiple fasta files with one label each
    if format == 'fasta':
        seq_stats = []
        for input_file in input:
            sequences = read_fasta(input_file)
            logging.debug(f"Read {len(sequences)} sequences from FASTA file {input_file}.")
            seq_stats += [SequenceStatistics(sequences, filename=Path(input_file).name, filepath=input_file,
                                             label=Path(input_file).stem, end_position=end_position)]
        run_analysis(
            input_statistics=seq_stats,
            out_folder=out_folder,
            report_types=report_types,
            plot_type=plot_type
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
                # drop rows with NaN values in label column
                nan_count = df[label_column].isna().sum()
                if nan_count > 0:
                    logging.warning(f"Dropped {nan_count} rows with non-numeric values in '{label_column}'.")
                    df = df.dropna(subset=[label_column])

                if len(df) == 0:
                    logging.error(f"No valid numeric values found in '{label_column}' for regression. Skipping analysis.")
                    return
                elif len(df) < 2:
                    logging.warning(f"Only {len(df)} sample(s) remaining after dropping non-numeric values. Results may not be meaningful.")

                # infer the threshold as the median of the label column
                threshold = df[label_column].median()
                logging.debug(f"Inferred threshold for regression: {threshold}")
                df[label_column] = df[label_column].apply(lambda x: 'high' if x >= threshold else 'low')
                labels = ['high', 'low']

            # get the list of labels to consider
            elif len(label_list) == 1 and label_list[0] == 'infer':
                labels = sorted(df[label_column].unique().tolist())
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
                                                     filepath=input[0], seq_column=seq_col, end_position=end_position)]
                run_analysis(
                    input_statistics=seq_stats,
                    out_folder=out_folder,
                    report_types=report_types,
                    plot_type=plot_type,
                )

            # handle multiple sequence columns by concatenating sequences and running statistics on them
            if len(sequence_column) > 1:
                seq_stats = []
                for label in labels:
                    sequences = read_multisequence_df(df, sequence_column, label_column, label)
                    seq_stats += [SequenceStatistics(sequences, filename=Path(input[0]).name, filepath=input[0],
                                                     label=label, seq_column='_'.join(sequence_column))]
                
                run_analysis(
                    input_statistics=seq_stats,
                    out_folder=out_folder,
                    report_types=report_types,
                    plot_type=plot_type,
                )

        # we have multiple files with one label each
        else:
            # run statistics across input files
            for seq_col in sequence_column:
                seq_stats = []
                for input_file in input:
                    sequences = read_sequences_from_df(read_csv_file(input_file, format, seq_col), seq_col)
                    logging.debug(f"Read {len(sequences)} sequences from file {input_file} in column '{seq_col}'.")
                    seq_stats += [SequenceStatistics(sequences, filename=Path(input_file).name, filepath=input_file,
                                                     label=Path(input_file).stem, seq_column=seq_col,
                                                     end_position=end_position)]
                run_analysis(
                    input_statistics=seq_stats,
                    out_folder=out_folder,
                    report_types=report_types,
                    plot_type=plot_type,
                )

            # handle multiple sequence columns
            if len(sequence_column) > 1:
                seq_stats = []
                for input_file in input:
                    sequences = read_multisequence_df(read_csv_file(input_file, format, sequence_column), sequence_column)
                    seq_stats += [SequenceStatistics(sequences, filename=Path(input_file).name, filepath=input_file,
                                                     label=Path(input_file).stem, seq_column='_'.join(sequence_column), 
                                                     end_position=end_position)]
                run_analysis(
                    input_statistics=seq_stats,
                    out_folder=out_folder,
                    report_types=report_types,
                    plot_type=plot_type,
                )

    logging.info("Classes evaluation successfully completed.")