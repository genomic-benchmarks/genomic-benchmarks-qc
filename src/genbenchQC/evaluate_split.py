import argparse
import logging
from pathlib import Path
from typing import Optional
import os
import shutil
from cdhit_reader import read_cdhit
import pandas as pd

from genbenchQC.report.report_generator import generate_json_report, generate_train_test_html_report, generate_simple_report
from genbenchQC.utils.input_utils import setup_logger, read_files_to_sequence_list, write_fasta

def run_clustering(train_fasta_file, test_fasta_file, clustered_file, identity_threshold, alignment_coverage):
    logging.info("Running CD-HIT clustering.")

    # Choose the word size for cd-hit based on the identity threshold
    if identity_threshold > 0.90:
        n = 10
    elif identity_threshold > 0.88:
        n = 7
    elif identity_threshold > 0.85:
        n = 6
    elif identity_threshold > 0.80:
        n = 5
    elif identity_threshold > 0.75:
        n = 4
    elif identity_threshold > 0.5:
        n = 3
    else:
        n = 2

    logging.debug("Running CD-HIT with the following parameters:")
    logging.debug(f"Input train file: {train_fasta_file}")
    logging.debug(f"Input test file: {test_fasta_file}")
    logging.debug(f"Output clustered file: {clustered_file}")
    logging.debug(f"Identity threshold: {identity_threshold}")
    logging.debug(f"Word size (n): {n}")
    logging.debug(f"Alignment coverage: {alignment_coverage}")

    errcode = os.system(f"cd-hit-est-2d -i {train_fasta_file} -i2 {test_fasta_file} -o {clustered_file} -c {identity_threshold} -n {n} -aS {alignment_coverage} -aL {alignment_coverage} -r 0 >/dev/null 2>&1")
    if errcode != 0:
        logging.error(f"CD-HIT clustering failed with error code {errcode}.")
        raise RuntimeError(f"CD-HIT clustering failed with error code {errcode}.")
    clusters = read_cdhit(f"{clustered_file}.clstr").read_items()
    logging.debug(f"CD-HIT clustering completed. {len(clusters)} clusters found.")

    # Collect clusters with >1 sequence
    mixed_clusters = []

    for cluster in clusters:
        if len(cluster.sequences) == 1:
            continue
        seq_ids = [seq.name for seq in cluster.sequences]
        mixed_clusters.append(seq_ids)
    return mixed_clusters

def process_mixed_clusters(clusters, train_sequences, test_sequences):
    sequence_clusters = []
    for i in range(len(clusters)):
        sequences = {"cluster": i, "train": [], "test": []}
        for seq_id in clusters[i]:
            seq_id = seq_id.split("_")
            if seq_id[2] == "train":
                sequences["train"].append(train_sequences[int(seq_id[1])])
            elif seq_id[2] == "test":
                sequences["test"].append(test_sequences[int(seq_id[1])])
            else:
                logging.warning(f"Unexpected sequence ID format: {seq_id}")
        sequence_clusters.append(sequences)

    return sequence_clusters

def run(train_files, test_files, format, 
        out_folder: Optional[str] = '.', 
        sequence_column: Optional[list[str]] = ['sequence'], 
        report_types: Optional[list[str]] = ['html', 'simple'], 
        identity_threshold: Optional[float] = 0.8, 
        alignment_coverage: Optional[float] = 0.8,
        log_level: Optional[str] = 'INFO',
        log_file: Optional[str] = None
    ):
    """Run the train-test split evaluation.

    This function reads sequences from the provided training and testing files, performs clustering using CD-HIT, 
    and generates reports about potential data leakage between the training and testing datasets.

    @param train_files: List of paths to training files.
    @param test_files: List of paths to testing files.
    @param format: Format of the input files (fasta, csv, csv.gz, tsv, tsv.gz).
    @param out_folder: Path to the output folder. Default: '.'.
    @param sequence_column: Name of the columns with sequences to analyze for datasets in CSV/TSV format. 
                            Default: ['sequence'].
    @param report_types: Types of reports to generate. Default: ['html', 'simple'].
    @param identity_threshold: Identity threshold for clustering. Default: 0.8.
    @param alignment_coverage: Alignment coverage for clustering. Default: 0.8.
    @param log_level: Logging level, default to INFO.
    @param log_file: Path to the log file. If provided, logs will be written to this file as well as to the console.
    @return: None
    """

    setup_logger(log_level, log_file)
    logging.info("Starting train-test split evaluation.")

    if not Path(out_folder).exists():
        logging.info(f"Output folder {out_folder} does not exist. Creating it.")
        Path(out_folder).mkdir(parents=True, exist_ok=True)

    Path(out_folder, "tmp").mkdir(parents=True, exist_ok=True)

    train_sequences = read_files_to_sequence_list(train_files, format, sequence_column)
    train_index = [f"{i}_train" for i in range(len(train_sequences))]
    logging.info(f"Read {len(train_sequences)} sequences from training files.")

    test_sequences = read_files_to_sequence_list(test_files, format, sequence_column)
    test_index = [f"{i}_test" for i in range(len(test_sequences))]
    logging.info(f"Read {len(test_sequences)} sequences from testing files.")

    train_fasta_path = Path(out_folder, "tmp") / 'train_sequences.fasta'
    write_fasta(train_sequences, train_fasta_path, train_index)
    test_fasta_path = Path(out_folder, "tmp") / 'test_sequences.fasta'
    write_fasta(test_sequences, test_fasta_path, test_index)

    clusters = run_clustering(train_fasta_path, test_fasta_path, Path(out_folder, "tmp/clustered_sequences"), identity_threshold, alignment_coverage)
    logging.debug(f"Having {len(clusters)} mixed clusters: {clusters}")

    filename = "split_check_" + Path(train_files[0]).stem + "_vs_" + Path(test_files[0]).stem

    if 'simple' in report_types:
        simple_report_path = Path(out_folder, filename + '.csv')
        result = {"Data leakage": "Pass" if not clusters else "Fail"}
        df = pd.DataFrame.from_dict(result, orient='index', columns=['Flag'])
        df.index.name = "Statistic"
        generate_simple_report(df, simple_report_path)

    if 'json' in report_types or 'html' in report_types:
        sequence_clusters = process_mixed_clusters(clusters, train_sequences, test_sequences)
        logging.debug(f"Transformed cluster sequence IDs to sequences: {sequence_clusters}")

    if 'json' in report_types:
        json_report_path = Path(out_folder, filename + '_report.json')
        generate_json_report({"mixed train-test clusters": sequence_clusters}, json_report_path)
    if 'html' in report_types:
        train_filenames = ",".join([Path(f).name for f in train_files])
        test_filenames = ",".join([Path(f).name for f in test_files])
        html_report_path = Path(out_folder, filename + '_report.html')
        generate_train_test_html_report(sequence_clusters, train_filenames, train_sequences, test_filenames, test_sequences, html_report_path, identity_threshold, alignment_coverage)

    # Clean up temporary files
    logging.debug("Removing temporary files.")
    shutil.rmtree(Path(out_folder, "tmp"))

    logging.info("Train-test split evaluation successfully completed.")

def parse_args():
    parser = argparse.ArgumentParser(description='Check data leakage in dataset train-test split.')
    parser.add_argument('--train_input', type=str, help='Path to the dataset file with training data. Can be multiple files that will be evaluated as one dataset part.', nargs='+', required=True)
    parser.add_argument('--test_input', type=str, help='Path to the dataset file with testing data. Can be multiple files that will be evaluated as one dataset part.', nargs='+',
                        required=True)
    parser.add_argument('--format', help="Format of the input files.", choices=['fasta', 'csv', 'csv.gz', 'tsv', 'tsv.gz'], required=True)
    parser.add_argument('--sequence_column', type=str, help='Name of the columns with sequences to analyze for datasets in CSV/TSV format. '
                                                            'Either one column or list of columns.', nargs='+', default=['sequence'])
    parser.add_argument('--out_folder', type=str, help='Path to the output folder.', default='.')
    parser.add_argument('--report_types', type=str, nargs='+', choices=['json', 'html', 'simple'],
                        help='Types of reports to generate. Default: [html]', default=['html', 'simple'])
    parser.add_argument('--identity_threshold', type=float, help='Identity threshold for clustering. Default: 0.8', default=0.8)
    parser.add_argument('--alignment_coverage', type=float, help='Alignment coverage for clustering. Default: 0.8', default=0.8)
    parser.add_argument('--log_level', type=str, help='Logging level, default to INFO.', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], default='INFO')
    parser.add_argument('--log_file', type=str, help='Path to the log file. If provided, logs will be written to this file as well as to the console.', default=None)
    args = parser.parse_args()

    return args

def main():
    args = parse_args()
    run(train_files = args.train_input, 
        test_files = args.test_input, 
        format = args.format, 
        out_folder = args.out_folder, 
        sequence_column = args.sequence_column, 
        report_types = args.report_types, 
        identity_threshold = args.identity_threshold, 
        alignment_coverage = args.alignment_coverage,
        log_level = args.log_level,
        log_file = args.log_file
    )

if __name__ == '__main__':
    main()