"""Per-class sequence statistics: the features the class comparison compares.

One `SequenceStatistics` holds the sequences of a single class and the features
computed from them - lengths, GC content, nucleotide and dinucleotide
composition, per-position composition, duplication levels. The comparison in
`genomic_benchmarks_qc.utils.testing` reads these; the plots and HTML report
read them too.

Two windows govern the per-position features, and they answer different
questions.

The scored window is where the flags may be set. Per-position statistics
compare only the sequences long enough to reach a position, so a position needs
a cohort that is both large enough to measure and representative enough to
speak for the class: at least
[MIN_SEQUENCES_PER_CLASS][genomic_benchmarks_qc.utils.testing.MIN_SEQUENCES_PER_CLASS]
sequences, and at least
[DEFAULT_MIN_COVERAGE][genomic_benchmarks_qc.utils.seq_stats.DEFAULT_MIN_COVERAGE]
of the class. The count is what binds on small and mid-sized classes; the
fraction binds on large ones, where a tail cohort can clear the count many times
over and still be nothing but the class's longest sequences.

The reported window is how far the per-position checks are named at all, and it
reaches much further. Past the scored window they can only say Unknown, and
saying it is the point: a report that stopped at the scored window would not
distinguish a dataset whose sequences end there from one whose tail was too thin
to compare. It stops where cohorts fall below
[MIN_SEQUENCES_PER_REPORTED_POSITION][genomic_benchmarks_qc.utils.seq_stats.MIN_SEQUENCES_PER_REPORTED_POSITION],
past which a position is reached by so few sequences that there is nothing left
to report about it either.

The figures draw the scored window and stop there. The tail is described in the
report's prose instead of drawn, because a curve is read as a measurement and
nothing out there was measured. The exception is a comparison with no scored
window at all: there the figures draw the reported one, every position Unknown,
which is what the rest of the report does with a comparison too small to score -
plot it, flag nothing.
"""

import logging
from collections import Counter

import numpy as np
import pandas as pd

from genomic_benchmarks_qc.utils.naming import slugify
from genomic_benchmarks_qc.utils.testing import MIN_SEQUENCES_PER_CLASS

# The two windows these bound are explained in the module docstring.
DEFAULT_MIN_COVERAGE = 0.25
MIN_SEQUENCES_PER_REPORTED_POSITION = 50

def cohort_floor(stats1, stats2) -> float:
    """The cohort floor a comparison of two classes runs under.

    Each class requires its own number of sequences behind a position - the
    larger of `MIN_SEQUENCES_PER_CLASS` and `min_coverage` of the class - so as a
    share of a class the floor differs between the two, and the binding one is
    the larger share: a position has to clear the floor in both classes.

    Not drawn. It sets the compared window, and the window is what the figures
    show: they stop where the floor stops them. Saying it a second time as a line
    across the coverage panel only added a rule to a panel that already carries
    two curves, so what the floor is gets said once, in the section's `?`
    explanation.

    Returns:
        Fraction of a class, or 0.0 when neither class has a floor worth drawing.
    """
    binding = 0.0
    for stats in (stats1, stats2):
        count = stats.stats['Number of sequences']
        if not count:
            continue
        needed = stats._required_cohort(count)
        binding = max(binding, needed / count)
    return binding


class SequenceStatistics:
    """The sequences of one class, and the statistics computed from them.

    Statistics are computed on demand by `compute`, not on construction, and
    cached in `stats`.
    """

    def __init__(self, sequences: list[str], filename: str, filepath: str, label: str,
                 seq_column: str | None = None, end_position: int | None = None,
                 slug: str | None = None, min_coverage: float = DEFAULT_MIN_COVERAGE):
        """Hold one class's sequences together with how to identify it.

        Args:
            sequences: The sequences of this class, uppercased by the reader.
            filename: Name of the file they came from, shown in the report.
            filepath: Full path they came from, shown in the report.
            label: The class name, shown verbatim in reports and plots.
            seq_column: Sequence column they came from, or None for FASTA.
            end_position: Last position the per-position checks reach, 1-based and
                inclusive. Defaults to the last position at least
                [MIN_SEQUENCES_PER_REPORTED_POSITION][genomic_benchmarks_qc.utils.seq_stats.MIN_SEQUENCES_PER_REPORTED_POSITION]
                of these sequences reach. It does not decide which positions may be
                flagged - the scored window does - and it is not what the figures draw,
                which is the scored window.
            slug: Path form of `label`; derived from it when not given, but normally
                passed in by the caller, which is the only place that can tell whether it
                collides with another class.
            min_coverage: Fraction of these sequences that must reach a position before it
                may set a flag, on top of the
                [MIN_SEQUENCES_PER_CLASS][genomic_benchmarks_qc.utils.testing.MIN_SEQUENCES_PER_CLASS]
                sequences every scored position needs. 0 leaves only that count.
                Default: `0.25`
                ([DEFAULT_MIN_COVERAGE][genomic_benchmarks_qc.utils.seq_stats.DEFAULT_MIN_COVERAGE]).
        """
        self.filename = filename
        self.filepath = filepath
        self.label = label
        # `label` is shown verbatim in reports and plots; `slug` is the
        # filesystem-safe, collision-free form used to build report paths.
        self.slug = slug if slug is not None else slugify(label)
        self.seq_column = seq_column
        self.sequences = sequences
        self.end_position = end_position
        self.min_coverage = min_coverage
        # Resolved by `compute`, alongside `end_position`.
        self.scored_end_position = None
        self.stats = {}

    def compute(self) -> tuple[dict, int]:
        """
        Compute various statistics from the given list of sequences.

        Results are also cached on `self.stats`, and the per-position windows
        (`self.end_position`, `self.scored_end_position`) are resolved here.

        The statistics dictionary holds:

        - Filename: str
        - Filepath: str
        - Label: str, or 'N/A'
        - Sequence column: str, or 'N/A'
        - Number of sequences: int
        - Number of bases: int
        - Unique bases: list of str
        - %GC content: float
        - Number of sequences left after deduplication: int
        - Empty sequences: int
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
        - Sequence duplication levels: dict {sequence: extra copies},
          holding only the sequences that occur more than once

        Returns:
            A tuple of that statistics dictionary and `end_position`.
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

        self._resolve_position_windows()

        return self.stats, self.end_position

    def _resolve_position_windows(self):
        """Resolve the reported window and the scored window from the lengths.

        `end_position` bounds the positions the per-position checks are named for;
        `scored_end_position` bounds the positions allowed to set a flag, and is
        the window the figures draw. The module docstring says why those are two
        different numbers. Both are 1-based and inclusive, and both are 0 when
        there is nothing to report or nothing to score.

        An explicit `end_position` is honoured as given, capped at the longest
        sequence. It cannot widen what gets flagged - the required cohort decides
        that - so an explicit window only ever drops positions off the tail of the
        results table, and drops them from the figures only where it cuts into the
        scored window.
        """

        col_info = f" for {self.seq_column} comparison" if self.seq_column is not None else ""

        lengths = self.stats['Sequence lengths'].values.flatten()

        if len(lengths) == 0:
            # No sequences means no positions at all. Guard here rather than
            # letting max() raise on the empty array, so the rest of the report
            # still gets written for a degenerate class.
            self.end_position = 0
            self.scored_end_position = 0
            logging.warning(
                f"No sequences{col_info}, so no positions are analysed in per-position statistics."
            )
            return

        scored_window = self._scored_window(lengths)

        if self.end_position is None:
            # The scored window asks for a much larger cohort than the reported
            # one, so the checks always reach at least as far as the flags. The
            # max is belt and braces for a class small enough to fall back on its
            # longest sequence: a position that sets a flag has to be reported.
            self.end_position = max(self._reported_window(lengths), scored_window)
            logging.info(
                f"No end position given, so per-position checks cover positions "
                f"1-{self.end_position}{col_info}, as far as at least "
                f"{MIN_SEQUENCES_PER_REPORTED_POSITION} sequences reach."
            )
        else:
            max_length = int(max(lengths))
            if self.end_position > max_length:
                logging.warning(
                    f"end_position {self.end_position} is greater than the maximum "
                    f"sequence length {max_length}. Setting end_position to {max_length}."
                )
                self.end_position = max_length

            logging.info(f"Using end position: {self.end_position}{col_info}.")

        # Never past the reported window: a position no check is named for has
        # nothing to flag, so trimming the checks with an explicit `end_position`
        # trims the scoring too. Left to itself it is never the binding
        # constraint.
        self.scored_end_position = min(scored_window, self.end_position)

        required = self._required_cohort(len(lengths))
        if self.scored_end_position < 1:
            logging.warning(
                f"Not enough sequences{col_info} for per-position statistics: no position is "
                f"reached by {required:,} of them, so every position is reported as Unknown. "
                "The figures are still drawn, with no flag on any position."
            )
            return

        floor = (f"{required:,} sequences ({self.min_coverage:.0%} of this class)"
                 if required > MIN_SEQUENCES_PER_CLASS
                 else f"{required:,} sequences")
        logging.info(
            f"Positions 1-{self.scored_end_position}{col_info} may be flagged, as far as "
            f"{floor} reach."
        )
        if self.end_position > self.scored_end_position:
            logging.info(
                f"Positions {self.scored_end_position + 1}-{self.end_position}{col_info} are "
                "reached by too few sequences to compare, so they are reported as Unknown and "
                "are not drawn."
            )

    def _reported_window(self, lengths) -> int:
        """Last position `MIN_SEQUENCES_PER_REPORTED_POSITION` sequences reach.

        Counted from the long end: a position is reached by at least k sequences
        exactly when it is no further out than the kth longest sequence. Past
        that a position belongs to a handful of sequences, and naming a check
        after it says more about those sequences than about the class.

        A class holding fewer sequences than the floor has no position that
        clears it, and stopping at 0 would leave it with no per-position checks at
        all. It falls back to its longest sequence: nothing there can be scored at
        that size anyway, so Unknown for every position is all it has to say.

        This floor is far below the one the scored window uses, so the reported
        window reaches past the flags rather than stopping short of them.
        """
        if len(lengths) < MIN_SEQUENCES_PER_REPORTED_POSITION:
            return int(max(lengths))
        return int(np.sort(lengths)[-MIN_SEQUENCES_PER_REPORTED_POSITION])

    def _required_cohort(self, count: int) -> int:
        """Sequences a position needs behind it before it may be flagged.

        The two floors are one number: a cohort has to be large enough to measure
        a difference and large enough to stand for the class, so it has to clear
        both `MIN_SEQUENCES_PER_CLASS` and `min_coverage` of the class. Rounded up,
        because a fraction of a sequence is not a sequence.
        """
        return max(MIN_SEQUENCES_PER_CLASS, int(np.ceil(self.min_coverage * count)))

    def _scored_window(self, lengths) -> int:
        """Last position the required cohort still reaches.

        Counted from the long end: a position is reached by at least k sequences
        exactly when it is no further out than the kth longest sequence. Because
        cohorts only shrink with position, this one number bounds the whole scored
        window and not just its final position.

        A class holding fewer sequences than it requires has no scorable position
        at all, and neither does one whose kth longest sequence is empty.
        """
        required = self._required_cohort(len(lengths))
        if len(lengths) < required:
            return 0
        return int(np.sort(lengths)[-required])

    def coverage_at(self, position: int) -> float:
        """Fraction of this class's sequences that reach `position` (1-based).

        This is the denominator behind every per-position statistic at that
        position, and it is what the report shows so a reader can tell how much
        data stands behind the far end of the per-position plots.
        """
        lengths = self.stats['Sequence lengths'].values.flatten()
        if len(lengths) == 0:
            return 0.0
        return float(np.mean(lengths >= position))

    def _compute_basic_statistics(self):
        """Compute the whole-class counts shown in the report header."""
        self.stats['Filename'] = self.filename
        self.stats['Filepath'] = self.filepath
        self.stats['Label'] = self.label if self.label is not None else 'N/A'
        self.stats['Sequence column'] = self.seq_column if self.seq_column is not None else 'N/A'
        self.stats['Number of sequences'] = len(self.sequences)
        total_bases = sum(len(sequence) for sequence in self.sequences)
        self.stats['Number of bases'] = total_bases
        unique_bases = set()
        for sequence in self.sequences:
            unique_bases.update(sequence)
        self.stats['Unique bases'] = sorted(unique_bases)
        gc_bases = sum(sequence.count('G') + sequence.count('C') for sequence in self.sequences)
        self.stats['%GC content'] = gc_bases / total_bases if total_bases > 0 else 0.0
        self.stats['Number of sequences left after deduplication'] = len(set(self.sequences))
        self.stats['Empty sequences'] = sum(1 for sequence in self.sequences if len(sequence) == 0)

    def _compute_per_sequence_statistics(self):
        """Compute every per-sequence and per-position feature in one pass."""

        # `compute` always runs `_compute_basic_statistics` first, which is what
        # puts 'Unique bases' in `stats`.
        nucleotides = self.stats['Unique bases']
        dinucleotides = [n1 + n2 for n1 in nucleotides for n2 in nucleotides]

        nucleotides_per_sequence = {}
        dinucleotides_per_sequence = {}
        nucleotides_per_position = {}
        nucleotides_per_position_reversed = {}
        gc_content_per_sequence = np.zeros(len(self.sequences))
        lengths_per_sequence = np.zeros(len(self.sequences))

        for id, sequence in enumerate(self.sequences):
            nucleotides_per_sequence[id] = self._compute_nucleotide_content(sequence, nucleotides)
            dinucleotides_per_sequence[id] = self._compute_dinucleotide_content(
                sequence, dinucleotides)
            self._compute_per_position_nucleotide_content(nucleotides_per_position, sequence)
            self._compute_per_position_nucleotide_content(
                nucleotides_per_position_reversed, sequence[::-1])
            seq_len = len(sequence)
            gc_bases = sequence.count('G') + sequence.count('C')
            gc_content_per_sequence[id] = gc_bases / seq_len * 100 if seq_len > 0 else 0.0
            lengths_per_sequence[id] = len(sequence)

        self.stats['Per sequence nucleotide content'] = pd.DataFrame(nucleotides_per_sequence).T
        self.stats['Per sequence dinucleotide content'] = pd.DataFrame(dinucleotides_per_sequence).T
        self.stats['Per position nucleotide content']= pd.DataFrame(
            self._normalize_per_position(nucleotides_per_position, nucleotides)).T
        self.stats['Per position reversed nucleotide content'] = pd.DataFrame(
            self._normalize_per_position(nucleotides_per_position_reversed, nucleotides)).T
        self.stats['Per sequence GC content'] = pd.DataFrame(
            gc_content_per_sequence, columns=['Per sequence GC content'])
        self.stats['Sequence lengths'] = pd.DataFrame(
            lengths_per_sequence, columns=['Sequence lengths'])

    def _compute_nucleotide_content(self, sequence, nucleotides):
        """Return the frequency of each nucleotide within one sequence."""
        seq_len = len(sequence)
        if seq_len == 0:
            return dict.fromkeys(nucleotides, 0)
        return {nucleotide: sequence.count(nucleotide) / seq_len for nucleotide in nucleotides}

    def _compute_dinucleotide_content(self, sequence, dinucleotides):
        """Return the frequency of each overlapping dinucleotide in one sequence.

        Dinucleotides not in `dinucleotides` are counted too rather than
        dropped, so an unexpected base cannot silently vanish from the totals.
        """
        dinucleotides_per_sequence = dict.fromkeys(dinucleotides, 0)
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
                dinucleotides_per_sequence[dinucleotide] = (
                    dinucleotides_per_sequence.get(dinucleotide, 0) + 1)

        total = sum(dinucleotides_per_sequence.values())
        if total == 0:
            return dinucleotides_per_sequence

        return {
            dinucleotide: count / total
            for dinucleotide, count in dinucleotides_per_sequence.items()
        }


    def _compute_per_position_nucleotide_content(self, nucleotides_per_position, sequence):
        """Add one sequence's bases to the per-position counts, in place."""
        for i, nucleotide in enumerate(sequence):
            if i in nucleotides_per_position:
                nucleotides_per_position[i][nucleotide] = (
                    nucleotides_per_position[i].get(nucleotide, 0) + 1)
            else:
                nucleotides_per_position[i] = {nucleotide: 1}

    def _normalize_per_position(self, nucleotides_per_position, nucleotides):
        """Turn per-position counts into frequencies, filling absent bases with 0.

        Each position is normalized by its own total, because positions beyond
        the shortest sequences are covered by fewer sequences than position 1.
        """
        for position in nucleotides_per_position:
            total = sum(nucleotides_per_position[position].values())
            nucleotides_per_position[position] = {
                nucleotide: count / total
                for nucleotide, count in nucleotides_per_position[position].items()
            }
            # add zeros for missing nucleotides
            for nucleotide in nucleotides:
                if nucleotide not in nucleotides_per_position[position]:
                    nucleotides_per_position[position][nucleotide] = 0
        return nucleotides_per_position

    def _compute_sequence_duplication_levels(self) -> None:
        """
        Compute the duplication levels for each sequence in the given list of sequences.

        Stores 'Sequence duplication levels' in `self.stats` as a dictionary of
        `{sequence: number of extra copies}`, ordered most duplicated first. Sequences
        occurring once are not included.
        """

        sequence_counts = Counter(self.sequences)
        # remove sequences that are not duplicated and decrement counts by 1 to
        # reflect number of duplications
        sequence_counts = {
            sequence: (count - 1)
            for sequence, count in sequence_counts.items() if count > 1
        }
        # sort the sequences by their counts
        sequence_counts = dict(
            sorted(sequence_counts.items(), key=lambda item: item[1], reverse=True))

        self.stats['Sequence duplication levels'] = sequence_counts
