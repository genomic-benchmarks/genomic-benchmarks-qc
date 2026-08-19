"""Tests for the per-position window.

The window over which per-position statistics are computed ends where most
sequences still reach: every position is compared on the sequences long enough
to have it, so the further out a position lies the fewer sequences stand behind
it.
"""

import numpy as np
import pytest

from genomic_benchmarks_qc.utils.seq_stats import (
    DEFAULT_MIN_COVERAGE,
    SequenceStatistics,
)


def make_stats(sequences, end_position=None, label='cls'):
    stats = SequenceStatistics(
        sequences=sequences, filename='f.fa', filepath='/f.fa',
        label=label, end_position=end_position,
    )
    stats.compute()
    return stats


def random_sequences(count, length, seed):
    """Sequences of one fixed length with uniform, seed-reproducible composition."""
    rng = np.random.default_rng(seed)
    return [''.join(rng.choice(list('ACGT'), size=length)) for _ in range(count)]


class TestEndPositionWindow:
    def test_default_keeps_the_required_coverage(self):
        # A quarter of the sequences stop at 20, the rest run to 100.
        sequences = random_sequences(400, 100, seed=10)
        sequences = [seq[:20] if i % 4 == 0 else seq for i, seq in enumerate(sequences)]

        stats = make_stats(sequences)

        assert stats.coverage_at(stats.end_position) >= DEFAULT_MIN_COVERAGE

    def test_default_is_the_complement_percentile_of_the_lengths(self):
        lengths = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        sequences = ['A' * length for length in lengths]

        stats = make_stats(sequences)

        expected = int(np.floor(np.percentile(lengths, 100 * (1 - DEFAULT_MIN_COVERAGE))))
        assert stats.end_position == expected

    def test_uniform_lengths_use_the_whole_sequence(self):
        stats = make_stats(random_sequences(50, 60, seed=11))

        assert stats.end_position == 60
        assert stats.coverage_at(60) == 1.0

    def test_explicit_end_position_is_capped_at_the_longest_sequence(self):
        stats = make_stats(random_sequences(20, 30, seed=12), end_position=500)

        assert stats.end_position == 30

    def test_explicit_end_position_below_the_cap_is_kept(self):
        stats = make_stats(random_sequences(20, 30, seed=13), end_position=10)

        assert stats.end_position == 10

    def test_thin_explicit_end_position_warns(self, caplog):
        # Only one sequence in ten reaches position 100.
        sequences = ['A' * 100] + ['A' * 10] * 9

        with caplog.at_level('WARNING'):
            make_stats(sequences, end_position=100)

        assert any('reach position 100' in record.message for record in caplog.records)

    def test_default_window_is_logged_at_info(self, caplog):
        """The auto-chosen window used to be DEBUG-only, invisible at the default level."""
        with caplog.at_level('INFO'):
            stats = make_stats(random_sequences(20, 45, seed=14))

        assert any(f'end position: {stats.end_position}' in record.message
                   for record in caplog.records)


class TestCoverageAt:
    @pytest.mark.parametrize("position, expected", [(1, 1.0), (10, 1.0), (11, 0.5), (20, 0.5), (21, 0.0)])
    def test_counts_sequences_reaching_the_position(self, position, expected):
        stats = make_stats(['A' * 10] * 5 + ['A' * 20] * 5)

        assert stats.coverage_at(position) == pytest.approx(expected)

    def test_empty_class_has_no_coverage(self):
        stats = make_stats([])

        assert stats.coverage_at(1) == 0.0
