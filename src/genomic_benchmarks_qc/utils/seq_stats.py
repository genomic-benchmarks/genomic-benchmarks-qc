"""Per-class sequence statistics: the features the class comparison compares.

One `SequenceStatistics` holds the sequences of a single class and the features
computed from them - lengths, GC content, nucleotide and dinucleotide
composition, per-position composition, duplication levels. The comparison in
`genomic_benchmarks_qc.utils.testing` reads these, and so do the plots and the
HTML report.

Two windows bound the per-position features, because a position is compared only
on the sequences long enough to reach it. `scored_end_position` is as far as a
position may be flagged: it needs a cohort of at least
[MIN_SEQUENCES_PER_CLASS][genomic_benchmarks_qc.utils.testing.MIN_SEQUENCES_PER_CLASS]
sequences and at least
[DEFAULT_MIN_COVERAGE][genomic_benchmarks_qc.utils.seq_stats.DEFAULT_MIN_COVERAGE]
of the class. It is also the window the figures draw. `end_position` reaches
further, as far as
[MIN_SEQUENCES_PER_REPORTED_POSITION][genomic_benchmarks_qc.utils.seq_stats.MIN_SEQUENCES_PER_REPORTED_POSITION]
sequences: the positions in between are named in the report as Unknown, so a
dataset whose sequences simply end is not read as one whose tail was too thin to
compare.

What the two windows mean for reading a report is on the
[per-position checks](../../guide/per-position.md) page.
"""

import logging
from collections import Counter

import numpy as np
import pandas as pd

from genomic_benchmarks_qc.utils.naming import slugify
from genomic_benchmarks_qc.utils.testing import MIN_SEQUENCES_PER_CLASS

logger = logging.getLogger(__name__)

# The two windows these bound are explained in the module docstring.
#
# Unlike MIN_SEQUENCES_PER_CLASS, this floor was not chosen by simulation: it is
# not a question about statistical power. A cohort far out along the sequence can
# clear the count many times over and still be nothing but the class's longest
# sequences, and no sample size fixes that - only stopping does.
DEFAULT_MIN_COVERAGE = 0.25
MIN_SEQUENCES_PER_REPORTED_POSITION = 50

# One vectorised pass over a class works on a block of its sequences rather than
# all of them, so that its working set does not grow with the class. Two things
# bound a block. The first is the per-character index arrays, which are eight
# bytes each per character.
STATS_BLOCK_CHARACTERS = 4_000_000
# The second is the per-sequence dinucleotide counts, one cell per sequence per
# dinucleotide: on short sequences, or on an alphabet larger than the four
# bases, they are what would otherwise be the largest array here.
STATS_BLOCK_CELLS = 4_000_000


def _sequence_blocks(lengths, cells_per_sequence):
    """Split a class into index ranges small enough for one vectorised pass.

    Args:
        lengths: Length of every sequence in the class, in order.
        cells_per_sequence: How many array cells one sequence costs, which is
            what `STATS_BLOCK_CELLS` is a budget of.

    Yields:
        `(start, stop)` pairs covering every sequence exactly once. A single
        sequence longer than the character budget still forms a block of its
        own, because the alternative is not making progress.
    """
    ends = np.cumsum(lengths)
    sequence_cap = max(1, STATS_BLOCK_CELLS // max(1, cells_per_sequence))
    start = 0
    while start < len(lengths):
        consumed = ends[start - 1] if start else 0
        stop = int(np.searchsorted(ends, consumed + STATS_BLOCK_CHARACTERS, side='right'))
        stop = min(max(stop, start + 1), start + sequence_cap, len(lengths))
        yield start, stop
        start = stop


def _as_frequencies(counts, totals):
    """Divide counts by their own totals, reading a zero total as zero.

    A sequence with no characters, and a position no sequence reaches, have
    nothing to be a fraction of. Both come back as zero everywhere rather than
    as a division by zero, which is what the dictionaries this replaced did by
    returning their zeros unchanged.

    Args:
        counts: Counts to normalize, one row per sequence or per position.
        totals: The total behind each row.

    Returns:
        Float array the shape of `counts`.
    """
    totals = np.asarray(totals, dtype=float)
    if counts.ndim == 2:
        totals = totals[:, None]
    return np.divide(counts, totals, out=np.zeros(counts.shape, dtype=float),
                     where=totals > 0)


def cohort_floor(stats1, stats2) -> float:
    """The share of a class a position's cohort has to reach to be flagged.

    Each class requires its own number of sequences behind a position - the
    larger of `MIN_SEQUENCES_PER_CLASS` and `min_coverage` of the class - so as a
    share of a class the floor differs between the two, and the binding one is
    the larger share: a position has to clear the floor in both classes.

    Returns:
        Fraction of a class, or 0.0 when neither class has a floor.
    """
    binding = 0.0
    for stats in (stats1, stats2):
        count = stats.stats['Number of sequences']
        if not count:
            continue
        needed = stats._required_cohort(count)
        binding = max(binding, needed / count)
    return binding


def _coverage_from_lengths(lengths, positions):
    """Fraction of `lengths` reaching each of `positions`, which are 1-based.

    A sequence reaches position p exactly when its length is at least p, so the
    count is the size of the sorted tail from p onwards - one binary search per
    position, rather than a positions-by-sequences comparison for every one of
    them. A class with no sequences reaches nothing.
    """
    lengths = np.asarray(lengths)
    positions = np.asarray(positions)
    if lengths.size == 0:
        return np.zeros(positions.shape, dtype=float)
    ordered = np.sort(lengths)
    return (ordered.size - np.searchsorted(ordered, positions, side='left')) / ordered.size


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
                of these sequences reach. It cannot widen what gets flagged - the
                scored window decides that - so an explicit value only ever trims.
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
        self.min_coverage = min_coverage

        self.end_position = end_position
        """Last position the per-position checks reach, resolved by `compute`."""

        self.scored_end_position = None
        """Last position that may be flagged, and the last one the figures draw.

        Resolved by `compute`, alongside `end_position`.
        """

        self.stats = {}
        """The statistics `compute` produces, empty until it has run."""

    def compute(self) -> tuple[dict, int]:
        """Compute this class's statistics and resolve its per-position windows.

        Results are cached on `self.stats`. The statistics dictionary holds:

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
            That statistics dictionary, and `end_position` - the last position
            the per-position checks reach.
        """
        message = f"Computing statistics for {self.filename}"
        if self.label is not None:
            message += f", label {self.label}"
        if self.seq_column is not None:
            message += f", sequence column: {self.seq_column}"
        logger.info(message)

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
            logger.warning(
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
            logger.info(
                f"No end position given, so per-position checks cover positions "
                f"1-{self.end_position}{col_info}, as far as at least "
                f"{MIN_SEQUENCES_PER_REPORTED_POSITION} sequences reach."
            )
        else:
            max_length = int(max(lengths))
            if self.end_position > max_length:
                logger.warning(
                    f"end_position {self.end_position} is greater than the maximum "
                    f"sequence length {max_length}. Setting end_position to {max_length}."
                )
                self.end_position = max_length

            logger.info(f"Using end position: {self.end_position}{col_info}.")

        # Never past the reported window: a position no check is named for has
        # nothing to flag, so trimming the checks with an explicit `end_position`
        # trims the scoring too. Left to itself it is never the binding
        # constraint.
        self.scored_end_position = min(scored_window, self.end_position)

        required = self._required_cohort(len(lengths))
        if self.scored_end_position < 1:
            logger.warning(
                f"Not enough sequences{col_info} for per-position statistics: no position is "
                f"reached by {required:,} of them, so every position is reported as Unknown. "
                "The figures are still drawn, with no flag on any position."
            )
            return

        floor = (f"{required:,} sequences ({self.min_coverage:.0%} of this class)"
                 if required > MIN_SEQUENCES_PER_CLASS
                 else f"{required:,} sequences")
        logger.info(
            f"Positions 1-{self.scored_end_position}{col_info} may be flagged, as far as "
            f"{floor} reach."
        )
        if self.end_position > self.scored_end_position:
            logger.info(
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

    def coverage_curve(self, end_position: int) -> np.ndarray:
        """Fraction reaching each 1-based position up to and including `end_position`.

        The denominator behind every per-position statistic, as a curve. The
        report needs it three times over - the figure draws it, the interactive
        viewer carries it, and the prose quotes single points of it - and it used
        to be worked out three ways, one of them a Python loop over the class per
        position. One answer, so the figure and the viewer cannot come to draw
        different curves.
        """
        return _coverage_from_lengths(
            self.stats['Sequence lengths'].values.flatten(),
            np.arange(1, end_position + 1),
        )

    def coverage_at(self, position: int) -> float:
        """Fraction of this class's sequences that reach `position` (1-based).

        This is the denominator behind every per-position statistic at that
        position, and it is what the report shows so a reader can tell how much
        data stands behind the far end of the per-position plots.
        """
        lengths = self.stats['Sequence lengths'].values.flatten()
        return float(_coverage_from_lengths(lengths, np.array([position]))[0])

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
        """Compute every per-sequence and per-position feature in one pass.

        The pass is vectorised. A block of sequences is joined into one string,
        read as an array of code points and turned into one column index per
        character, after which every feature here is a `bincount` over that
        array divided by its own total. Walking the characters in Python was
        71% of the statistics for a class of 10,000 500 bp sequences, nearly
        all of it the two per-position dictionaries.
        """

        # `compute` always runs `_compute_basic_statistics` first, which is what
        # puts 'Unique bases' in `stats`.
        nucleotides = self.stats['Unique bases']
        dinucleotides = [n1 + n2 for n1 in nucleotides for n2 in nucleotides]
        # A dinucleotide's column is `first * base_count + second`, which is the
        # order `dinucleotides` is built in, so the counts need no lookup table.
        base_count = len(nucleotides)

        count = len(self.sequences)
        lengths = np.fromiter((len(sequence) for sequence in self.sequences),
                              dtype=np.int64, count=count)
        max_length = int(lengths.max()) if count else 0

        # Counts first and frequencies at the end: every feature below is one of
        # these four tables divided by a total the table itself carries.
        per_sequence_bases = np.zeros((count, base_count), dtype=np.int64)
        per_sequence_pairs = np.zeros((count, base_count ** 2), dtype=np.int64)
        per_position = np.zeros((max_length, base_count), dtype=np.int64)
        per_position_reversed = np.zeros((max_length, base_count), dtype=np.int64)

        # 'Unique bases' is sorted, so it is sorted by code point too, and
        # `searchsorted` turns a character into its column in one step. Every
        # character of every sequence is in it by construction - that is what
        # makes it 'unique bases' - so there is no miss to handle, and no
        # dinucleotide can fall outside the ones named above either.
        codes = np.array([ord(base) for base in nucleotides], dtype=np.uint32)

        for start, stop in _sequence_blocks(lengths, base_count ** 2):
            block = lengths[start:stop]
            text = ''.join(self.sequences[start:stop])
            if not text:
                continue
            # UTF-32 is four bytes per code point whatever the code point, so
            # the array lines up with the string however exotic a base is.
            column = np.searchsorted(
                codes, np.frombuffer(text.encode('utf-32-le'), dtype=np.uint32))
            # Which sequence each character belongs to, and how far into it it
            # sits - the two indices every count below is grouped by.
            row = np.repeat(np.arange(stop - start), block)
            position = np.arange(len(column)) - np.repeat(np.cumsum(block) - block, block)

            per_sequence_bases[start:stop] = np.bincount(
                row * base_count + column,
                minlength=(stop - start) * base_count,
            ).reshape(stop - start, base_count)
            per_position += np.bincount(
                position * base_count + column,
                minlength=max_length * base_count,
            ).reshape(max_length, base_count)
            # The reversed feature counts from the far end of each sequence,
            # which is the same characters read under a different index.
            per_position_reversed += np.bincount(
                (np.repeat(block, block) - 1 - position) * base_count + column,
                minlength=max_length * base_count,
            ).reshape(max_length, base_count)

            # Overlapping pairs, less the ones that would straddle two
            # sequences: a pair starting at the last character of a sequence is
            # exactly the one whose second character sits at position 0.
            inside = position[1:] > 0
            per_sequence_pairs[start:stop] = np.bincount(
                row[:-1][inside] * base_count ** 2
                + column[:-1][inside] * base_count + column[1:][inside],
                minlength=(stop - start) * base_count ** 2,
            ).reshape(stop - start, base_count ** 2)

        sequence_ids = np.arange(count)
        self.stats['Per sequence nucleotide content'] = pd.DataFrame(
            _as_frequencies(per_sequence_bases, lengths),
            index=sequence_ids, columns=nucleotides)
        # A sequence's dinucleotides are a fraction of its own pairs, of which
        # there is one fewer than it has characters.
        self.stats['Per sequence dinucleotide content'] = pd.DataFrame(
            _as_frequencies(per_sequence_pairs, per_sequence_pairs.sum(axis=1)),
            index=sequence_ids, columns=dinucleotides)
        # Each position is normalized by its own total, because positions beyond
        # the shortest sequences are covered by fewer sequences than position 1.
        positions = np.arange(max_length)
        self.stats['Per position nucleotide content'] = pd.DataFrame(
            _as_frequencies(per_position, per_position.sum(axis=1)),
            index=positions, columns=nucleotides)
        self.stats['Per position reversed nucleotide content'] = pd.DataFrame(
            _as_frequencies(per_position_reversed, per_position_reversed.sum(axis=1)),
            index=positions, columns=nucleotides)

        gc_bases = np.zeros(count, dtype=np.int64)
        for base in ('G', 'C'):
            if base in nucleotides:
                gc_bases += per_sequence_bases[:, nucleotides.index(base)]
        self.stats['Per sequence GC content'] = pd.DataFrame(
            _as_frequencies(gc_bases, lengths) * 100,
            columns=['Per sequence GC content'])
        self.stats['Sequence lengths'] = pd.DataFrame(
            lengths.astype(float), columns=['Sequence lengths'])

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
