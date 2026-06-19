from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import logging
import json


class SequenceStatsAccumulator:
    """Streaming accumulator for sequence length statistics."""

    def __init__(self):
        """Initialize empty counters for sequence length aggregation."""
        self.count = 0
        self.total_length = 0
        self.min_length = None
        self.max_length = None

    def add(self, sequence):
        """Update aggregate length metrics with one sequence."""
        seq_len = len(sequence)
        self.count += 1
        self.total_length += seq_len
        if self.min_length is None or seq_len < self.min_length:
            self.min_length = seq_len
        if self.max_length is None or seq_len > self.max_length:
            self.max_length = seq_len

    def finalize(self):
        """Return finalized count/min/mean/max statistics as a dictionary."""
        count = self.count
        return {
            "count": count,
            "min_length": self.min_length if self.min_length is not None else 0,
            "mean_length": (self.total_length / count) if count > 0 else 0.0,
            "max_length": self.max_length if self.max_length is not None else 0,
        }

def stream_fasta_sequences(fasta_file):
    """Yield uppercase sequences one-by-one from a FASTA file."""
    logging.debug(f"Streaming FASTA file: {fasta_file}")
    for record in SeqIO.parse(fasta_file, 'fasta'):
        yield str(record.seq).upper()


def read_fasta(fasta_file):
    """Read all sequences from a FASTA file into a list."""
    logging.debug(f"Reading FASTA file: {fasta_file}")
    return [str(record.seq).upper() for record in SeqIO.parse(fasta_file, 'fasta')]


def read_selected_fasta_sequences(fasta_file, ids_to_keep):
    """Read a FASTA file and return uppercase sequences for selected IDs."""
    ids_to_keep = set(ids_to_keep)
    if not ids_to_keep:
        return {}

    out = {}
    for record in SeqIO.parse(str(fasta_file), 'fasta'):
        seq_id = record.id
        if seq_id in ids_to_keep:
            out[seq_id] = str(record.seq).upper()
            if len(out) == len(ids_to_keep):
                break

    return out


def stream_fasta_records_by_ids(fasta_path, ids_to_keep):
    """Stream FASTA records whose IDs are in ids_to_keep."""
    ids_to_keep = set(ids_to_keep)
    if not ids_to_keep:
        return

    for record in SeqIO.parse(str(fasta_path), 'fasta'):
        if record.id in ids_to_keep:
            yield record


def filter_fasta_by_ids(fasta_path, new_fasta_path, ids_to_keep):
    """Write only selected FASTA records to a new output file."""
    logging.debug(f"Filtering FASTA file: {fasta_path} -> {new_fasta_path}, keeping {len(ids_to_keep)} IDs")
    records = stream_fasta_records_by_ids(fasta_path, ids_to_keep)
    SeqIO.write(records, str(new_fasta_path), 'fasta')


def append_fasta_record(file_handle, sequence, seq_id):
    """Append one sequence to an open FASTA handle as a SeqRecord."""
    SeqIO.write(
        SeqRecord(Seq(sequence), id=seq_id, description=""),
        file_handle,
        "fasta"
    )


def read_csv_file(file_path, input_format, seq_columns, label_column=None):
    """Read CSV/TSV data and normalize sequence columns to uppercase strings."""
    delim = '\t' if input_format == 'tsv' or input_format == 'tsv.gz' else ','
    compression = 'gzip' if file_path.endswith('.gz') else None

    columns = seq_columns.copy()
    if label_column is not None:
        columns += [label_column]

    df = pd.read_csv(file_path, delimiter=delim, usecols=columns, dtype=str, compression=compression)
    
    # Drop rows with missing labels
    if label_column is not None:
        # check if label column contains any missing values
        if df[label_column].isnull().any():
            logging.warning(f"Label column '{label_column}' contains missing values. Dropping rows with missing labels.")
        df = df.dropna(subset=[label_column])
        logging.debug(f"Dropped rows with missing labels, new shape: {df.shape}")

    # Replace NaN values in sequence columns with empty strings
    df[seq_columns] = df[seq_columns].fillna('')

    # Convert sequences to uppercase
    df[seq_columns] = df[seq_columns].apply(lambda col: col.str.upper())

    logging.debug(f"Read CSV/TSV file: {file_path}, shape: {df.shape}, columns: {columns}")

    return df


def read_sequences_from_df(df, seq_columns, label_column=None, label=None):
    """Extract sequences from a dataframe with optional label filtering.

    `seq_columns` may be a single column name (str) or a list of column names.
    If multiple columns are provided the columns are concatenated per-row.
    If `label_column` is provided, rows are filtered by `label`.
    """
    # Normalize seq_columns to a list
    if isinstance(seq_columns, str):
        seq_columns = [seq_columns]

    if label_column is not None:
        logging.debug(f"Filtering sequences by label: {label} in column: {label_column}")
        df_parsed = df[df[label_column] == label]
        if df_parsed.empty:
            logging.error(f"No sequences found for label '{label}' in column '{label_column}'.")
            raise ValueError(f"No sequences found for label '{label}' in column '{label_column}'.")
    else:
        df_parsed = df

    for col in seq_columns:
        if col not in df_parsed.columns:
            logging.error(f"Sequence column '{col}' not found in dataframe.")
            raise KeyError(f"Sequence column '{col}' not found in dataframe.")

    if len(seq_columns) == 1:
        return df_parsed[seq_columns[0]].tolist()
    else:
        return df_parsed[seq_columns].agg(''.join, axis=1).tolist()

def stream_table_sequences(file_path, input_format, seq_columns, chunksize=10000):
    """Yield sequences from CSV/TSV files in chunks to limit memory usage."""
    delim = '\t' if input_format == 'tsv' or input_format == 'tsv.gz' else ','
    compression = 'gzip' if file_path.endswith('.gz') else None

    reader = pd.read_csv(
        file_path,
        delimiter=delim,
        usecols=seq_columns,
        dtype=str,
        compression=compression,
        chunksize=chunksize,
    )

    for chunk in reader:
        for col in seq_columns:
            chunk[col] = chunk[col].fillna('').str.upper()

        if len(seq_columns) == 1:
            for seq in chunk[seq_columns[0]].tolist():
                yield seq
        else:
            for seq in chunk[seq_columns].agg(''.join, axis=1).tolist():
                yield seq


def stream_files_to_sequences(files, input_format, sequence_column, chunksize=10000):
    """Yield sequences from multiple input files based on declared format."""
    for file in files:
        if input_format == 'fasta':
            yield from stream_fasta_sequences(file)
        elif input_format.startswith('csv') or input_format.startswith('tsv'):
            yield from stream_table_sequences(file, input_format, sequence_column, chunksize=chunksize)
        else:
            logging.error(f"Unsupported input format: {input_format}")
            raise ValueError(f"Unsupported input format: {input_format}")


def setup_logger(level=logging.INFO, file=None):
    """Configure console/file logging and suppress noisy matplotlib debug logs."""
    if file:
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(file, mode='w'),
                logging.StreamHandler()
            ]
        )
    else:
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.StreamHandler()
            ]
        )

    # Suppress matplotlib debug messages
    logging.getLogger("matplotlib").setLevel(logging.WARNING)


def write_stats_json(stats, stats_json_file):
    """Serialize computed statistics to JSON, converting DataFrames to plain dicts."""
    stats_dict = {}
    for key, value in stats.items():
        if isinstance(value, pd.DataFrame):
            stats_dict[key] = value.to_dict(orient='list')
        else:
            stats_dict[key] = value

    with open(stats_json_file, 'w') as file:
        json.dump(stats_dict, file, indent=4)
