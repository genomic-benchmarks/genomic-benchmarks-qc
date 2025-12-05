from Bio import SeqIO, SeqRecord, Seq
import pandas as pd
import logging
import json

def read_fasta(fasta_file):
    logging.debug(f"Reading FASTA file: {fasta_file}")
    return [str(record.seq).upper() for record in SeqIO.parse(fasta_file, 'fasta')]

def write_fasta(sequences, output_file, indices=None):
    if indices is None:
        records = [SeqRecord.SeqRecord(Seq.Seq(seq), id=f'seq_{i}', description="") for i, seq in enumerate(sequences)]
    else:
        records = [SeqRecord.SeqRecord(Seq.Seq(sequences[i]), id=f'seq_{indices[i]}', description="") for i in range(len(sequences))]
    logging.debug(f"Writing FASTA file: {output_file} with {len(sequences)} sequences")
    SeqIO.write(records, output_file, 'fasta')

def read_csv_file(file_path, input_format, seq_columns, label_column=None):
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

    # Convert sequences to uppercase
    df[seq_columns] = df[seq_columns].apply(lambda col: col.str.upper())

    logging.debug(f"Read CSV/TSV file: {file_path}, shape: {df.shape}, columns: {columns}")

    return df

def read_sequences_from_df(df, seq_column, label_column=None, label=None):
    if label_column is None:
        return df[seq_column].tolist()

    logging.debug(f"Filtering sequences by label: {label} in column: {label_column}")
    df_parsed = df[df[label_column] == label]

    if df_parsed.empty:
        logging.error(f"No sequences found for label '{label}' in column '{label_column}'.")
        raise ValueError(f"No sequences found for label '{label}' in column '{label_column}'.")

    return df_parsed[seq_column].tolist()

def read_multisequence_df(df, seq_columns, label_column=None, label=None):
    if len(seq_columns) > 1:
        logging.debug(f"Concatenating sequences from multiple columns: {seq_columns}")
    all_sequences = [read_sequences_from_df(df, seq_column, label_column, label) for seq_column in seq_columns]
    concatenated_sequences = [''.join(seqs) for seqs in zip(*all_sequences)]
    return concatenated_sequences

def setup_logger(level=logging.INFO, file=None):
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

    stats_dict = {}
    for key, value in stats.items():
        if isinstance(value, pd.DataFrame):
            stats_dict[key] = value.to_dict(orient='list')
        else:
            stats_dict[key] = value

    with open(stats_json_file, 'w') as file:
        json.dump(stats_dict, file, indent=4)

def read_files_to_sequence_list(files, input_format, sequence_column):
    sequences = []
    for file in files:
        if input_format == 'fasta':
            sequences += read_fasta(file)
        elif input_format.startswith('csv') or input_format.startswith('tsv'):
            df = read_csv_file(file, input_format, sequence_column)
            sequences += read_multisequence_df(df, sequence_column)
        else:
            logging.error(f"Unsupported input format: {input_format}")
            raise ValueError(f"Unsupported input format: {input_format}")
    return sequences
