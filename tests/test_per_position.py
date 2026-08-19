"""Tests for the per-position window and the scoring inside it.

Three decisions are pinned here. A position is compared only on the sequences
that reach it, so that per-position flags answer a question about composition
rather than re-answering the one `Sequence lengths` already asks. The window
ends where most sequences still reach. And within the window, the difference a
position must show before it may set the check's verdict widens as its cohort
shrinks -- without that, the worst case over hundreds of positions reports
sampling noise as the check's verdict.
"""

import numpy as np
import pytest

from genomic_benchmarks_qc.utils.seq_stats import (
    DEFAULT_MIN_COVERAGE,
    SequenceStatistics,
)
from genomic_benchmarks_qc.utils.testing import (
    MIN_SEQUENCES_PER_POSITION,
    _adaptive_threshold,
    _compute_position_binary_scores,
    _extract_failed_features,
    _position_cohorts,
    _score_position_features,
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


class TestPositionBinaryScores:
    def test_only_sequences_reaching_the_position_contribute(self):
        # Position 3 (0-based) exists in the first two sequences only.
        scores = _compute_position_binary_scores(['AAAA', 'AAAC', 'AA'], 'A', 3, reverse=False)

        assert scores.tolist() == [1.0, 0.0]

    def test_short_sequences_are_dropped_not_scored_as_absent(self):
        """The distinction the old zero-fill lost: absent base vs absent position."""
        reaching = _compute_position_binary_scores(['CCCC'], 'A', 3, reverse=False)
        not_reaching = _compute_position_binary_scores(['CC'], 'A', 3, reverse=False)

        assert reaching.tolist() == [0.0]
        assert not_reaching.size == 0

    def test_reverse_counts_from_the_end(self):
        # One position in from the 3' end: 'T' in the first, 'A' in the second.
        scores = _compute_position_binary_scores(['GGTA', 'GGAC'], 'T', 1, reverse=True)

        assert scores.tolist() == [1.0, 0.0]

    def test_reverse_also_drops_sequences_that_are_too_short(self):
        scores = _compute_position_binary_scores(['ACGT', 'AC'], 'A', 3, reverse=True)

        assert scores.tolist() == [1.0]


class TestPositionCohorts:
    def test_counts_sequences_long_enough_for_each_position(self):
        cohorts = _position_cohorts(['AAA', 'AA', 'A'], end_position=3)

        assert cohorts.tolist() == [3, 2, 1]

    def test_no_window_gives_no_cohorts(self):
        assert _position_cohorts(['AAA'], end_position=0).tolist() == []


class TestScoringIsNotDrivenByLength:
    def test_identical_composition_with_different_lengths_passes(self):
        """The bias the old zero-fill introduced: length alone must not flag.

        Both classes share one composition; only their lengths differ. Scoring
        short sequences as "base absent" would make the second class look
        depleted at every position past its own lengths.
        """
        long_sequences = random_sequences(600, 80, seed=1)
        # Same generator, then half the sequences truncated to a third the length.
        mixed = random_sequences(600, 80, seed=2)
        mixed = [seq[:26] if i % 2 else seq for i, seq in enumerate(mixed)]

        stats1 = make_stats(long_sequences)
        stats2 = make_stats(mixed)
        end_position = min(stats1.end_position, stats2.end_position)

        results, _ = _score_position_features(
            stats1.sequences, stats2.sequences, ['A', 'C', 'G', 'T'],
            'Per position nucleotide content', end_position=end_position,
        )

        flags = {name: metrics['Flag'] for name, metrics in results.items()}
        assert set(flags.values()) <= {'Pass'}, [n for n, f in flags.items() if f != 'Pass']

    def test_a_real_composition_difference_is_still_flagged(self):
        """The counterpart: conditioning must not blunt a genuine difference."""
        sequences_1 = ['A' + seq[1:] for seq in random_sequences(600, 40, seed=3)]
        sequences_2 = ['C' + seq[1:] for seq in random_sequences(600, 40, seed=4)]

        results, _ = _score_position_features(
            sequences_1, sequences_2, ['A'],
            'Per position nucleotide content', end_position=40,
        )

        assert results['Per position nucleotide content - A position 1']['Flag'] == 'Fail'


class TestThinPositionsAreNotScored:
    def test_position_below_the_floor_is_unknown(self):
        sequences_1 = random_sequences(MIN_SEQUENCES_PER_POSITION - 1, 10, seed=5)
        sequences_2 = random_sequences(MIN_SEQUENCES_PER_POSITION + 50, 10, seed=6)

        results, _ = _score_position_features(
            sequences_1, sequences_2, ['A'],
            'Per position nucleotide content', end_position=10,
        )

        metrics = results['Per position nucleotide content - A position 1']
        assert metrics['Flag'] == 'Unknown'
        assert np.isnan(metrics['AU-ROC'])

    def test_positions_past_the_cohort_floor_stop_being_scored(self):
        """Within one window, the far positions can fall below the floor.

        Every sequence reaches position 1, but only a handful reach position 20,
        so the same comparison scores the near positions and reports the far
        ones as Unknown.
        """
        sequences = ['A' * 20] * 10 + ['A' * 5] * 200

        results, _ = _score_position_features(
            sequences, sequences, ['A'],
            'Per position nucleotide content', end_position=20,
        )

        assert results['Per position nucleotide content - A position 1']['Flag'] != 'Unknown'
        assert results['Per position nucleotide content - A position 20']['Flag'] == 'Unknown'

    def test_an_all_unknown_feature_aggregates_to_unknown(self):
        sequences = random_sequences(10, 10, seed=7)

        _, per_base = _score_position_features(
            sequences, sequences, ['A'],
            'Per position nucleotide content', end_position=10,
        )

        assert per_base['A']['Flag'] == 'Unknown'
        assert np.isnan(per_base['A']['AU-ROC'])

    def test_a_thin_cohort_cannot_manufacture_a_flag(self):
        """Small classes with identical composition used to reach Fail by chance.

        30 sequences per class over enough positions produces a worst-case
        AU-ROC around 0.67 under the null. Below the cohort floor those
        positions become Unknown instead of a Warning on the aggregate.
        """
        sequences_1 = random_sequences(30, 100, seed=8)
        sequences_2 = random_sequences(30, 100, seed=9)

        _, per_base = _score_position_features(
            sequences_1, sequences_2, ['A', 'C', 'G', 'T'],
            'Per position nucleotide content', end_position=100,
        )

        assert {aggregate['Flag'] for aggregate in per_base.values()} == {'Unknown'}


class TestAdaptiveThreshold:
    def test_tightens_as_the_number_of_tests_grows(self):
        values_1 = np.array([1.0] * 30 + [0.0] * 70)
        values_2 = np.array([1.0] * 40 + [0.0] * 60)

        one_test = _adaptive_threshold(values_1, values_2, n_tests=1)
        many_tests = _adaptive_threshold(values_1, values_2, n_tests=1600)

        assert 0.5 < one_test < many_tests

    def test_relaxes_below_the_fixed_boundary_on_large_cohorts(self):
        """On a normally sized dataset the rule must reduce to the tool's own
        boundaries, so that a decision only ever changes where data are thin."""
        small_1 = np.array([1.0] * 30 + [0.0] * 70)
        small_2 = np.array([1.0] * 40 + [0.0] * 60)
        large_1 = np.array([1.0] * 1500 + [0.0] * 3500)
        large_2 = np.array([1.0] * 2000 + [0.0] * 3000)

        assert _adaptive_threshold(small_1, small_2, n_tests=1600) > 0.6
        assert _adaptive_threshold(large_1, large_2, n_tests=1600) < 0.6

    def test_a_difference_within_sampling_noise_is_reported_at_chance(self):
        """A hundred sequences per class over a hundred positions crosses 0.6 on
        noise alone; the threshold turns those positions back into chance."""
        sequences_1 = random_sequences(100, 100, seed=20)
        sequences_2 = random_sequences(100, 100, seed=21)

        results, per_base = _score_position_features(
            sequences_1, sequences_2, ['A', 'C', 'G', 'T'],
            'Per position nucleotide content', end_position=100,
        )

        assert {aggregate['Flag'] for aggregate in per_base.values()} == {'Pass'}
        # Chance is reported as chance, not as the value that happened to come out.
        assert max(metrics['AU-ROC'] for metrics in results.values()) == pytest.approx(0.5)

    def test_a_real_difference_survives_the_same_cohort_size(self):
        """The threshold must cost specificity, not the signal it was built around."""
        sequences_1 = ['A' + seq[1:] for seq in random_sequences(100, 100, seed=22)]
        sequences_2 = ['C' + seq[1:] for seq in random_sequences(100, 100, seed=23)]

        results, _ = _score_position_features(
            sequences_1, sequences_2, ['A'],
            'Per position nucleotide content', end_position=100,
        )

        assert results['Per position nucleotide content - A position 1']['Flag'] == 'Fail'

    def test_a_scored_position_carries_a_consistent_metric_set(self):
        """A position reported at chance reports chance for every metric, so the
        row cannot show an AU-ROC of 0.5 next to an accuracy that suggests
        otherwise."""
        sequences = random_sequences(100, 20, seed=24)

        results, _ = _score_position_features(
            sequences, sequences, ['A'],
            'Per position nucleotide content', end_position=20,
        )

        metrics = results['Per position nucleotide content - A position 1']
        assert metrics['AU-ROC'] == pytest.approx(0.5)
        assert metrics['AU-PR'] == pytest.approx(0.5)
        assert metrics['Accuracy'] == pytest.approx(0.5)


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


class TestFailedFeatureExtraction:
    def test_a_clean_comparison_lists_no_positions(self):
        """An entry per base regardless of its flag leaves the dict truthy, and
        the report then renders a second, identical copy of every per-position
        plot."""
        results = {
            'Per position nucleotide content - A position 1': {'Flag': 'Pass'},
            'Per position nucleotide content - A': {'Flag': 'Pass'},
        }

        failed = _extract_failed_features(results)

        assert failed['Per position nucleotide content'] == {}

    def test_flagged_positions_are_kept_under_their_base(self):
        results = {
            'Per position nucleotide content - A position 1': {'Flag': 'Pass'},
            'Per position nucleotide content - G position 7': {'Flag': 'Fail'},
            'Per position reversed nucleotide content - T position 2': {'Flag': 'Warning'},
        }

        failed = _extract_failed_features(results)

        assert failed['Per position nucleotide content'] == {'G': {7: 'Fail'}}
        assert failed['Per position reversed nucleotide content'] == {'T': {2: 'Warning'}}
