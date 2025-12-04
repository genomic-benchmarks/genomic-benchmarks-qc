import logging
from collections import Counter
from typing import Optional
import numpy as np
import pandas as pd

class SequenceStatistics:
    def __init__(self, sequences: list[str], filename: str, filepath: str, label: str, seq_column: Optional[str] = None, end_position: Optional[int] = None):
        self.filename = filename
        self.filepath = filepath
        self.label = label
        self.seq_column = seq_column
        self.sequences = sequences
        self.end_position = end_position
        self.stats = {}

    def compute(self):
        """
        Compute various statistics from the given list of sequences.
        @return: A dictionary containing the following statistics:
            - Filename: str
            - Number of sequences: int
            - Number of bases: int
            - Unique bases: list of str
            - %GC content: float
            - Number of sequences left after deduplication: int
            - Per sequence nucleotide content: pd.DataFrame
              (index: sequence_id, columns: nucleotides, values: frequency)
            - Per sequence dinucleotide content: pd.DataFrame
              (index: sequence_id, columns: dinucleotides, values: frequency)
            - Per position nucleotide content: pd.DataFrame
              (index: position, columns: nucleotides, values: frequency)
            - Per position reversed nucleotide content: pd.DataFrame
              (index: position, columns: nucleotides, values: frequency)
            - Per sequence GC content: dict pd.DataFrame
              (index: sequence_id, columns: GC content (%), values: GC content)
            - Sequence lengths: pd.DataFrame
              (index: sequence_id, columns: Length, values: length of the sequence)
            - Sequence duplication levels: dict {sequence: count}
        """
        message = f"Computing statistics for {self.filename}"
        if self.label is not None:
            message += f", label {self.label}"
        if self.seq_column is not None:
            message += f", sequence column: {self.seq_column}"
        logging.info(message)

        self._compute_basic_statistics()
        self._compute_per_sequence_statistics()
        self._compute_sequence_duplication_levels()

        self._adjust_end_position()

        return self.stats, self.end_position

    def _adjust_end_position(self):

        col_info = f" for {self.seq_column} comparison" if self.seq_column is not None else ""

        if self.end_position is None:

            # get second end position - where one of the stats contains less then 75% values
            lengths = self.stats['Sequence lengths'].values.flatten()
            lengths_75th = np.percentile(lengths, 75)
            # round to nearest integer
            self.end_position = int(np.round(lengths_75th))

            logging.info(
                f"End position not provided. Using end position: {self.end_position}{col_info}. "
                 "This is the 75th percentile of sequence lengths."
            )
        else:
            # Ensure end_position is not greater than the maximum sequence length
            lengths = self.stats['Sequence lengths'].values.flatten()
            max_length = int(max(lengths))
            if self.end_position > max_length:
                logging.warning(f"end_position {self.end_position} is greater than the maximum sequence length {max_length}. Setting end_position to {max_length}.")
                self.end_position = max_length

            logging.info(f"Using end position: {self.end_position}{col_info}.")

    def _compute_basic_statistics(self):
        self.stats['Filename'] = self.filename
        self.stats['Filepath'] = self.filepath
        self.stats['Label'] = self.label if self.label is not None else 'N/A'
        self.stats['Sequence column'] = self.seq_column if self.seq_column is not None else 'N/A'
        self.stats['Number of sequences'] = len(self.sequences)
        self.stats['Number of bases'] = sum(len(sequence) for sequence in self.sequences)
        self.stats['Unique bases'] = list(set(''.join(self.sequences)))
        total_bases = sum(len(sequence) for sequence in self.sequences)
        self.stats['%GC content'] = sum(sequence.count('G') + sequence.count('C') for sequence in self.sequences) / total_bases if total_bases > 0 else 0.0
        self.stats['Number of sequences left after deduplication'] = len(set(self.sequences))
        self.stats['Empty sequences'] = sum(1 for sequence in self.sequences if len(sequence) == 0)

    def _compute_per_sequence_statistics(self):

        nucleotides = self.stats['Unique bases'] if 'Unique bases' in self.stats else list(set(''.join(self.sequences)))
        dinucleotides = [n1 + n2 for n1 in nucleotides for n2 in nucleotides]

        nucleotides_per_sequence = {}
        dinucleotides_per_sequence = {}
        nucleotides_per_position = {}
        nucleotides_per_position_reversed = {}
        gc_content_per_sequence = np.zeros(len(self.sequences))
        lengths_per_sequence = np.zeros(len(self.sequences))

        for id, sequence in enumerate(self.sequences):
            nucleotides_per_sequence[id] = self._compute_nucleotide_content(sequence, nucleotides)
            dinucleotides_per_sequence[id] = self._compute_dinucleotide_content(sequence, dinucleotides)
            self._compute_per_position_nucleotide_content(nucleotides_per_position, sequence)
            self._compute_per_position_nucleotide_content(nucleotides_per_position_reversed, sequence[::-1])
            seq_len = len(sequence)
            gc_content_per_sequence[id] = (sequence.count('G') + sequence.count('C')) / seq_len * 100 if seq_len > 0 else 0.0
            lengths_per_sequence[id] = len(sequence)

        self.stats['Per sequence nucleotide content'] = pd.DataFrame(nucleotides_per_sequence).T
        self.stats['Per sequence dinucleotide content'] = pd.DataFrame(dinucleotides_per_sequence).T
        self.stats['Per position nucleotide content']= pd.DataFrame(
            self._normalize_per_position(nucleotides_per_position, nucleotides)).T
        self.stats['Per position reversed nucleotide content'] = pd.DataFrame(
            self._normalize_per_position(nucleotides_per_position_reversed, nucleotides)).T
        self.stats['Per sequence GC content'] = pd.DataFrame(gc_content_per_sequence, columns=['Per sequence GC content'])
        self.stats['Sequence lengths'] = pd.DataFrame(lengths_per_sequence, columns=['Sequence lengths'])

    def _compute_nucleotide_content(self, sequence, nucleotides):
        seq_len = len(sequence)
        if seq_len == 0:
            return {nucleotide: 0 for nucleotide in nucleotides}
        return {nucleotide: sequence.count(nucleotide) / seq_len for nucleotide in nucleotides}   
     
    def _compute_dinucleotide_content(self, sequence, dinucleotides):
        dinucleotides_per_sequence = {dinucleotide: 0 for dinucleotide in dinucleotides}
        seq_len = len(sequence)

        # No dinucleotides possible for sequences shorter than 2 -> return zeros
        if seq_len < 2:
            return dinucleotides_per_sequence

        for i in range(seq_len - 1):
            dinucleotide = sequence[i:i + 2]
            # increment only known dinucleotides, but track unexpected ones too
            if dinucleotide in dinucleotides_per_sequence:
                dinucleotides_per_sequence[dinucleotide] += 1
            else:
                dinucleotides_per_sequence[dinucleotide] = dinucleotides_per_sequence.get(dinucleotide, 0) + 1

        total = sum(dinucleotides_per_sequence.values())
        if total == 0:
            return dinucleotides_per_sequence
        
        dinucleotides_per_sequence = {dinucleotide: count / total for dinucleotide, count in dinucleotides_per_sequence.items()}

        return dinucleotides_per_sequence
    
    def _compute_per_position_nucleotide_content(self, nucleotides_per_position, sequence):
        for i, nucleotide in enumerate(sequence):
            if i in nucleotides_per_position:
                nucleotides_per_position[i][nucleotide] = nucleotides_per_position[i].get(nucleotide, 0) + 1
            else:
                nucleotides_per_position[i] = {nucleotide: 1}

    def _normalize_per_position(self, nucleotides_per_position, nucleotides):
        for position in nucleotides_per_position:
            total = sum(nucleotides_per_position[position].values())
            nucleotides_per_position[position] = {nucleotide: count / total for nucleotide, count in nucleotides_per_position[position].items()}
            # add zeros for missing nucleotides
            for nucleotide in nucleotides:
                if nucleotide not in nucleotides_per_position[position]:
                    nucleotides_per_position[position][nucleotide] = 0
        return nucleotides_per_position
    
    def _compute_sequence_duplication_levels(self):
        """
        Compute the duplication levels for each sequence in the given list of sequences.
        @param sequences: A list of sequences.
        @return: A dictionary containing the duplication levels for duplicated sequences. Unique sequences are not included.
        """

        sequence_counts = Counter(self.sequences)
        # remove sequences that are not duplicated
        sequence_counts = {sequence: count for sequence, count in sequence_counts.items() if count > 1}
        # sort the sequences by their counts
        sequence_counts = dict(sorted(sequence_counts.items(), key=lambda item: item[1], reverse=True))
        
        self.stats['Sequence duplication levels'] = sequence_counts