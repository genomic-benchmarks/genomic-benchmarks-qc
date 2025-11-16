import logging
from collections import Counter
from typing import Optional
import numpy as np
import pandas as pd

class SequenceStatistics:
    def __init__(self, sequences: list[str], filename: str, label: str, seq_column: Optional[str] = None, end_position: Optional[int] = None):
        """
        Initialize a SequenceStatistics instance with input sequences and associated metadata.
        
        Parameters:
            sequences (list[str]): List of sequence strings to analyze.
            filename (str): Source filename associated with the sequences.
            label (str): Human-readable label for the sequence set.
            seq_column (Optional[str]): Optional name of the sequence column (if sequences were read from tabular data).
            end_position (Optional[int]): Optional end position to use for per-position comparisons; if None, an appropriate end position will be determined later.
        
        Initializes:
            stats (dict): Empty dictionary to be populated with computed statistics.
        """
        self.filename = filename
        self.label = label
        self.seq_column = seq_column
        self.sequences = sequences
        self.end_position = end_position
        self.stats = {}

    def compute(self):
        """
        Aggregate sequence-level and position-level statistics for the stored sequences.
        
        Populates self.stats with computed values (examples and typical types shown):
        - 'Filename', 'Label', 'Sequence column' (str)
        - 'Number of sequences' (int)
        - 'Number of bases' (int)
        - 'Unique bases' (list[str])
        - '%GC content' (float)
        - 'Number of sequences left after deduplication' (int)
        - 'Per sequence nucleotide content' (pd.DataFrame): frequencies per nucleotide for each sequence
        - 'Per sequence dinucleotide content' (pd.DataFrame): frequencies per dinucleotide for each sequence
        - 'Per position nucleotide content' (pd.DataFrame) and 'Per position reversed nucleotide content' (pd.DataFrame): nucleotide frequencies by position
        - 'Per sequence GC content' (pd.DataFrame)
        - 'Sequence lengths' (pd.DataFrame)
        - 'Sequence duplication levels' (dict): mapping duplicated sequence -> count
        
        Returns:
            tuple: (stats, end_position)
            - stats (dict): dictionary of computed statistics described above.
            - end_position (int | None): finalized end position used for position-based comparisons (may be None if not set).
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
        if self.end_position is None:

            # get second end position - where one of the stats contains less then 75% values
            lengths = self.stats['Sequence lengths'].values.flatten()
            lengths_75th = np.percentile(lengths, 75)
            # round to nearest integer
            self.end_position = int(np.round(lengths_75th))

            logging.info(
                f"End position not provided. Using end position: {self.end_position} for {self.seq_column} comparison. "
                 "This is the 75th percentile of sequence lengths."
            )
        else:
            # Ensure end_position is not greater than the maximum sequence length
            lengths = self.stats['Sequence lengths'].values.flatten()
            max_length = int(max(lengths))
            if self.end_position > max_length:
                logging.warning(f"end_position {self.end_position} is greater than the maximum sequence length {max_length}. Setting end_position to {max_length}.")
                self.end_position = max_length

            logging.info(f"Using end position: {self.end_position} for {self.seq_column} comparison.")

    def _compute_basic_statistics(self):
        self.stats['Filename'] = self.filename
        self.stats['Label'] = self.label if self.label is not None else 'N/A'
        self.stats['Sequence column'] = self.seq_column if self.seq_column is not None else 'N/A'
        self.stats['Number of sequences'] = len(self.sequences)
        self.stats['Number of bases'] = sum(len(sequence) for sequence in self.sequences)
        self.stats['Unique bases'] = list(set(''.join(self.sequences)))
        self.stats['%GC content'] = sum(sequence.count('G') + sequence.count('C') for sequence in self.sequences) / sum(len(sequence) for sequence in self.sequences)
        self.stats['Number of sequences left after deduplication'] = len(set(self.sequences))

    def _compute_per_sequence_statistics(self):

        """
        Compute per-sequence and per-position nucleotide and dinucleotide statistics and store them in self.stats.
        
        Populates self.stats with:
        - 'Per sequence nucleotide content': DataFrame of nucleotide frequency per sequence.
        - 'Per sequence dinucleotide content': DataFrame of dinucleotide frequency per sequence.
        - 'Per position nucleotide content': DataFrame of nucleotide frequency at each position (forward orientation).
        - 'Per position reversed nucleotide content': DataFrame of nucleotide frequency at each position for reversed sequences.
        - 'Per sequence GC content': DataFrame with per-sequence GC percentage.
        - 'Sequence lengths': DataFrame with per-sequence lengths.
        
        The method derives the nucleotide alphabet from existing stats ('Unique bases') if available, otherwise from the input sequences. It updates internal arrays and maps and writes the resulting pandas DataFrames into self.stats; it does not return a value.
        """
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
            gc_content_per_sequence[id] = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
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
        """
        Compute per-nucleotide frequencies for a single sequence.
        
        Parameters:
            sequence (str): The nucleotide sequence to analyze.
            nucleotides (Iterable[str]): The nucleotides to report frequencies for.
        
        Returns:
            dict: Mapping of each nucleotide from `nucleotides` to its frequency (count divided by sequence length). If `sequence` is empty, all frequencies are 0. Values are in the range [0.0, 1.0].
        """
        seq_len = len(sequence)
        if seq_len == 0:
            return {nucleotide: 0 for nucleotide in nucleotides}
        return {nucleotide: sequence.count(nucleotide) / seq_len for nucleotide in nucleotides}   
     
    def _compute_dinucleotide_content(self, sequence, dinucleotides):
        """
        Compute the frequency distribution of dinucleotides observed in a sequence.
        
        This returns normalized frequencies for dinucleotides found in `sequence`. The result always includes entries for the provided `dinucleotides` (initialized to 0). If the sequence is shorter than 2, the initialized mapping is returned unchanged. Any dinucleotide not present in the initial `dinucleotides` iterable but observed in `sequence` will also be included in the returned mapping.
        
        Parameters:
            sequence (str): Nucleotide sequence to analyze.
            dinucleotides (iterable[str]): Known dinucleotides to initialize in the result.
        
        Returns:
            dict[str, float]: Mapping from dinucleotide string to its relative frequency (count / total). Values sum to 1 when at least one dinucleotide is counted; otherwise the mapping contains zeros for the initialized dinucleotides.
        """
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