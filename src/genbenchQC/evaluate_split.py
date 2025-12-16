import argparse
import logging
from pathlib import Path
from typing import Optional
import os
import shutil
import pandas as pd

from genbenchQC.report.report_generator import generate_json_report, generate_split_html_report, generate_simple_report
from genbenchQC.utils.input_utils import setup_logger, read_files_to_sequence_list, write_fasta
from genbenchQC.utils.similarity_threshold import dinucleotide_shuffle_list, compute_threshold, get_basic_stats, get_threshold_stats

def run_hashFrag_stratify_test_split(train_fasta_file, test_fasta_file, out_folder): #, shuffled=False):
    logging.info("Running hashFrag stratify test split.")

    logging.debug("Running hashFrag stratify test split with the following parameters:")
    logging.debug(f"Input train file: {train_fasta_file}")
    logging.debug(f"Input test file: {test_fasta_file}")
    logging.debug(f"Output directory: {out_folder}")

    out_folder.mkdir(parents=True, exist_ok=True)
    
    errcode = os.system(f"hashFrag stratify_test_split --train-fasta-path {train_fasta_file} --test-fasta-path {test_fasta_file} -o {out_folder} >/dev/null 2>&1")
    if errcode != 0:
        logging.error(f"hashFrag stratify test split failed with error code {errcode}.")
        raise RuntimeError(f"hashFrag stratify test split failed with error code {errcode}.")

    stratified_test_split = pd.read_csv(out_folder / "hashFrag.stratified_test_split.tsv", sep="\t")
    logging.info("hashFrag stratify test split completed successfully.")

    # if shuffled:
    #     return stratified_test_split
    
    # test_seqs_processed = pd.read_csv(out_folder / "test_sequencessta.blastn.processed.tsv", sep="\t", header=None, names=["test_id", "train_id", "score"])
    # return stratified_test_split, test_seqs_processed

    return stratified_test_split

def run(train_files, test_files, format, 
        out_folder: Optional[str] = '.', 
        sequence_column: Optional[list[str]] = ['sequence'], 
        report_types: Optional[list[str]] = ['html', 'simple'], 
        log_level: Optional[str] = 'INFO',
        log_file: Optional[str] = None
    ):
    """Run the train-test split evaluation.

    This function reads sequences from the provided training and testing files into fasta format, runs hashFrag, 
    and generates reports about potential data leakage between the training and testing datasets.

    @param train_files: List of paths to training files.
    @param test_files: List of paths to testing files.
    @param format: Format of the input files (fasta, csv, csv.gz, tsv, tsv.gz).
    @param out_folder: Path to the output folder. Default: '.'.
    @param sequence_column: Name of the columns with sequences to analyze for datasets in CSV/TSV format. 
                            Default: ['sequence'].
    @param report_types: Types of reports to generate. Default: ['html', 'simple'].
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
    train_fasta_path = Path(out_folder, "tmp") / 'train_sequences.fasta'
    write_fasta(train_sequences, train_fasta_path, train_index)

    test_sequences = read_files_to_sequence_list(test_files, format, sequence_column)
    test_index = [f"{i}_test" for i in range(len(test_sequences))]
    logging.info(f"Read {len(test_sequences)} sequences from testing files.")
    test_fasta_path = Path(out_folder, "tmp") / 'test_sequences.fasta'
    write_fasta(test_sequences, test_fasta_path, test_index)
    
    out_folder_hashFrag_genomic = Path(out_folder) / "hashFrag_genomic"
    stratified_test_split = run_hashFrag_stratify_test_split(train_fasta_path, test_fasta_path, Path(out_folder_hashFrag_genomic))
    
    train_seqs_shuffled = dinucleotide_shuffle_list(train_sequences, seed=42)
    train_fasta_path_shuffled = Path(out_folder, "tmp") / 'train_sequences_shuffled.fasta'
    write_fasta(train_seqs_shuffled, train_fasta_path_shuffled, train_index)
    
    test_seqs_shuffled = dinucleotide_shuffle_list(test_sequences, seed=42)
    test_fasta_path_shuffled = Path(out_folder, "tmp") / 'test_sequences_shuffled.fasta'
    write_fasta(test_seqs_shuffled, test_fasta_path_shuffled, test_index)

    out_folder_hashFrag_shuffled = Path(out_folder) / "hashFrag_shuffled"
    stratified_test_split_shuffled = run_hashFrag_stratify_test_split(train_fasta_path_shuffled, test_fasta_path_shuffled, Path(out_folder_hashFrag_shuffled)) #, shuffled=True)

    threshold = compute_threshold(stratified_test_split_shuffled)

    filename = "split_check_" + Path(train_files[0]).stem + "_vs_" + Path(test_files[0]).stem
    # if 'simple' in report_types:
    #     simple_report_path = Path(out_folder, filename + '.csv')
    #     result = {"Data leakage": "Pass" if not clusters else "Fail"}
    #     df = pd.DataFrame.from_dict(result, orient='index', columns=['Flag'])
    #     df.index.name = "Statistic"
    #     generate_simple_report(df, simple_report_path)
    # if 'json' in report_types or 'html' in report_types:
    #     sequence_clusters = process_mixed_clusters(clusters, train_sequences, test_sequences)
    #     logging.debug(f"Transformed cluster sequence IDs to sequences: {sequence_clusters}")
    # if 'json' in report_types:
    #     json_report_path = Path(out_folder, filename + '_report.json')
    #     generate_json_report({"mixed train-test clusters": sequence_clusters}, json_report_path)
    if 'html' in report_types:
        train_filenames = ",".join([Path(f).name for f in train_files])
        test_filenames = ",".join([Path(f).name for f in test_files])
        html_report_path = Path(out_folder, filename + '_report.html')
        plots_path = Path(out_folder, filename + '_plots')

        basic_stats = get_basic_stats(train_filenames, train_sequences, test_filenames, test_sequences)
        threshold_stats = get_threshold_stats(stratified_test_split, threshold)

        generate_split_html_report(stratified_test_split, stratified_test_split_shuffled, basic_stats, threshold_stats, html_report_path, plots_path)

    # Clean up temporary files
    logging.debug("Removing temporary files.")
    # shutil.rmtree(Path(out_folder, "tmp"))

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
        log_level = args.log_level,
        log_file = args.log_file
    )

if __name__ == '__main__':
    main()