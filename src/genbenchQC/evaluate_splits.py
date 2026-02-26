import logging
from pathlib import Path
from typing import Optional
import subprocess
import shutil
import pandas as pd

from genbenchQC.report.report_generator import generate_splits_html_report, generate_simple_report
from genbenchQC.utils.input_utils import setup_logger, read_files_to_sequence_list, write_fasta
from genbenchQC.utils.data_leakage_utils import get_basic_stats, get_threshold_stats, filter_fasta_by_ids

def run_search(test_fasta_file, train_fasta_file, out_file, tmp_dir):
    logging.info(
        "Running MMSeqs2, an ultrafast and sensitive search, for test sequences (query) against train sequences (db)."
    )

    if shutil.which("mmseqs") is None:
        raise RuntimeError(
            "MMSeqs2 executable not found in PATH. "
            "Please install MMSeqs2 and ensure it is available in your environment."
        )

    cmd = [
        "mmseqs", "easy-search",
        str(test_fasta_file),
        str(train_fasta_file),
        str(out_file),
        str(tmp_dir),
        "--format-output",
        "query,target,qcov,tcov,pident,evalue,qstart,qend,tstart,tend,alnlen,qseq,tseq,qaln,taln",
        "--format-mode", "4",
        "--search-type", "3",
        "--strand", "1"
    ]

    logging.debug(f"Running command: {' '.join(cmd)}")

    try:
        subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

    except subprocess.CalledProcessError as e:
        logging.error("MMSeqs2 search failed.")
        logging.error(f"Return code: {e.returncode}")
        if e.stderr:
            logging.error(f"STDERR:\n{e.stderr.strip()}")
        if e.stdout:
            logging.debug(f"STDOUT:\n{e.stdout.strip()}")
        raise RuntimeError("MMSeqs2 search failed.") from e

    logging.debug("MMSeqs2 easy-search completed.")

    return pd.read_csv(out_file, sep="\t", header=0)

def run(train_files, test_files, format, 
        out_folder: Optional[str] = '.', 
        sequence_column: Optional[list[str]] = None, 
        report_types: Optional[list[str]] = None, 
        similarity_threshold: Optional[float] = 80.0, 
        log_level: Optional[str] = 'INFO',
        log_file: Optional[str] = None
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

    try: 
        # Read sequences from training and testing files
        train_sequences = read_files_to_sequence_list(train_files, format, sequence_column)
        train_index = [f"{i}_train" for i in range(len(train_sequences))]
        logging.info(f"Read {len(train_sequences)} sequences from training files.")
        test_sequences = read_files_to_sequence_list(test_files, format, sequence_column)
        test_index = [f"{i}_test" for i in range(len(test_sequences))]
        logging.info(f"Read {len(test_sequences)} sequences from testing files.")

        # Write sequences to temporary fasta files for MMseqs2
        train_fasta_path = tmp_dir / 'train_sequences.fasta'
        write_fasta(train_sequences, train_fasta_path, train_index)
        test_fasta_path = tmp_dir / 'test_sequences.fasta'
        write_fasta(test_sequences, test_fasta_path, test_index)
        
        # Run MMseqs2 search and get results as a DataFrame
        outfile ="mmseqs2_search_result.tsv"
        results = run_search(test_fasta_path, train_fasta_path, tmp_dir / outfile, tmp_dir)

        # Add columns to results for similarity and whether the hit is above the threshold for potential leakage
        results['min_cov'] = results[['qcov', 'tcov']].min(axis=1)
        results['min_cov*pident'] = results['min_cov'] * results['pident']
        results['Leaked'] = results.apply(lambda row: 'True' if row['min_cov*pident'] >= similarity_threshold else 'False', axis=1)

        # Build filename for reports based on input file names
        filename = ("split_check_" 
                    + Path(train_files[0]).name.replace("".join(Path(train_files[0]).suffixes), "") 
                    + "_vs_" 
                    + Path(test_files[0]).name.replace("".join(Path(test_files[0]).suffixes), ""))


        # Filter results for hits above threshold and sort by similarity
        results_filt = results[results['Leaked'] == 'True'].sort_values(by=['min_cov*pident'], ascending=False)

        # Get threshold stats
        num_train_seqs = len(train_sequences)
        num_test_seqs = len(test_sequences)
        threshold_stats = get_threshold_stats(results, results_filt, similarity_threshold, num_train_seqs, num_test_seqs)

        if 'simple' in report_types:
            simple_report_path = out_folder / (filename + '.csv')
            has_leakage = (results['Leaked'] == 'True').any()
            result = {}
            result["Data Leakage"] = {"Flag": "Fail" if has_leakage else "Pass"}
            result["Data Leakage"]["Percentage of leaked queries"] = f"{threshold_stats['perc_queries_above_thr']:.2f}%"
            result["Data Leakage"]["Percentage of leaked targets"] = f"{threshold_stats['perc_targets_above_thr']:.2f}%"
            df = pd.DataFrame.from_dict(result, orient='index')
            df.index.name = "Statistic"
            generate_simple_report(df, simple_report_path)
        
        if 'html' in report_types:
            # Build paths for output files
            train_filenames = ",".join([Path(f).name for f in train_files])
            test_filenames = ",".join([Path(f).name for f in test_files])
            html_report_path = out_folder / (filename + '_report.html')
            plots_dir = out_folder / (filename + '_plots')
            mmseqs_dir = out_folder / (filename + '_mmseqs')
            seq_index_mapping = mmseqs_dir / 'seq_index_mapping'

            # Write mmseqs2 results, including flag for hits above threshold, to file 
            mmseqs_dir.mkdir(parents=True, exist_ok=True)
            results.to_csv(mmseqs_dir / outfile, sep='\t', index=False)

            # Write filtered fasta files for mapping hits seq IDs back to seqs
            seq_index_mapping.mkdir(parents=True, exist_ok=True)
            new_test_fasta_path = seq_index_mapping / 'test_sequences.fasta'
            new_train_fasta_path = seq_index_mapping / 'train_sequences.fasta'
            filter_fasta_by_ids(test_fasta_path, new_test_fasta_path, set(results_filt["query"]))
            filter_fasta_by_ids(train_fasta_path, new_train_fasta_path, set(results_filt["target"]))

            # # Prepare results (hits over threshold) for HTML report
            # results_filt['qcov'] = results_filt['qcov'].round(2)
            # results_filt['tcov'] = results_filt['tcov'].round(2)
            # results_filt['evalue'] = results_filt['evalue'].apply(lambda x: f"{x:.2e}")
            
            # Get basic stats for HTML report
            basic_stats = get_basic_stats(train_filenames, train_sequences, test_filenames, test_sequences)

            # Generate HTML report
            generate_splits_html_report(basic_stats, threshold_stats, results, results_filt, html_report_path, plots_dir)

            logging.info("Train-test split evaluation successfully completed.")
    except Exception:
        logging.exception(f"Train-test split evaluation failed. ")
        raise
    finally: 
        logging.debug("Removing temporary files.")
        try:
            shutil.rmtree(tmp_dir)
        except Exception as cleanup_error:
            logging.warning(f"Failed to remove temporary directory: {cleanup_error}")