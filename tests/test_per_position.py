"""Tests for the per-position windows and the scoring inside them.

Three decisions are pinned here. A position is compared only on the sequences
that reach it, so that per-position flags answer a question about composition
rather than re-answering the one `Sequence lengths` already asks. A position is
compared only where enough of them do: the larger of MIN_SEQUENCES_PER_CLASS
sequences and a fraction of the class, one guard against sampling noise and one
against a cohort that is large but is only the class's longest sequences. And the
reported window runs further than either, as far as a cohort worth naming a check
after survives -- what the report accounts for and what it compares are separate
questions, and the figures draw the second.
"""

import numpy as np
import pytest

from genomic_benchmarks_qc.report.per_position_payload import _coverage
from genomic_benchmarks_qc.utils.seq_stats import (
    DEFAULT_MIN_COVERAGE,
    MIN_SEQUENCES_PER_REPORTED_POSITION,
    SequenceStatistics,
    cohort_floor,
)
from genomic_benchmarks_qc.utils.testing import (
    MIN_SEQUENCES_PER_CLASS,
    _compute_position_binary_scores,
    _extract_failed_features,
    _position_cohorts,
    _score_position_features,
)


def make_stats(sequences, end_position=None, label='cls', min_coverage=DEFAULT_MIN_COVERAGE):
    stats = SequenceStatistics(
        sequences=sequences, filename='f.fa', filepath='/f.fa',
        label=label, end_position=end_position, min_coverage=min_coverage,
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
            scored_end_position=end_position,
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
            scored_end_position=40,
        )

        assert results['Per position nucleotide content - A position 1']['Flag'] == 'Fail'


class TestThinPositionsAreNotScored:
    def test_position_below_the_floor_is_unknown(self):
        sequences_1 = random_sequences(MIN_SEQUENCES_PER_CLASS - 1, 10, seed=5)
        sequences_2 = random_sequences(MIN_SEQUENCES_PER_CLASS + 50, 10, seed=6)

        results, _ = _score_position_features(
            sequences_1, sequences_2, ['A'],
            'Per position nucleotide content', end_position=10,
            scored_end_position=10,
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
        sequences = ['A' * 20] * 100 + ['A' * 5] * 200

        results, _ = _score_position_features(
            sequences, sequences, ['A'],
            'Per position nucleotide content', end_position=20,
            scored_end_position=20,
        )

        assert results['Per position nucleotide content - A position 1']['Flag'] != 'Unknown'
        assert results['Per position nucleotide content - A position 20']['Flag'] == 'Unknown'

    def test_an_all_unknown_feature_aggregates_to_unknown(self):
        sequences = random_sequences(10, 10, seed=7)

        _, per_base = _score_position_features(
            sequences, sequences, ['A'],
            'Per position nucleotide content', end_position=10,
            scored_end_position=10,
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
            scored_end_position=100,
        )

        assert {aggregate['Flag'] for aggregate in per_base.values()} == {'Unknown'}


class TestTheFixedBoundaryHoldsAtTheCohortFloor:
    """Why the per-position checks can share one boundary with the rest.

    A per-position check reduces hundreds of (position, base) pairs to their worst
    case, so the boundary it aggregates against has to survive that maximum under
    the null. It does, at the cohort floor and above: the floor is set where a
    difference of the size the boundary asks for stops turning up by chance. That
    is what lets these checks use the same fixed 0.6/0.7 as every other check,
    instead of a threshold that moves with each cohort.
    """

    def test_a_null_comparison_at_the_floor_does_not_flag(self):
        """Two classes from one process, at exactly the floor, over a hundred
        positions and four bases: the worst of those four hundred tests still has
        to come out Pass."""
        sequences_1 = random_sequences(MIN_SEQUENCES_PER_CLASS, 100, seed=20)
        sequences_2 = random_sequences(MIN_SEQUENCES_PER_CLASS, 100, seed=21)

        results, per_base = _score_position_features(
            sequences_1, sequences_2, ['A', 'C', 'G', 'T'],
            'Per position nucleotide content', end_position=100,
            scored_end_position=100,
        )

        assert {aggregate['Flag'] for aggregate in per_base.values()} == {'Pass'}
        worst = max(metrics['AU-ROC'] for metrics in results.values()
                    if not np.isnan(metrics['AU-ROC']))
        assert worst < 0.6, worst

    def test_a_real_difference_at_the_floor_is_still_flagged(self):
        """The floor must cost specificity, not the signal it was built around."""
        sequences_1 = ['A' + seq[1:]
                       for seq in random_sequences(MIN_SEQUENCES_PER_CLASS, 100, seed=22)]
        sequences_2 = ['C' + seq[1:]
                       for seq in random_sequences(MIN_SEQUENCES_PER_CLASS, 100, seed=23)]

        results, _ = _score_position_features(
            sequences_1, sequences_2, ['A'],
            'Per position nucleotide content', end_position=100,
            scored_end_position=100,
        )

        assert results['Per position nucleotide content - A position 1']['Flag'] == 'Fail'

    def test_a_scored_position_reports_what_was_measured(self):
        """A difference too small to flag is reported as the number it was.

        The threshold this replaced overwrote such a position with chance, so the
        table showed 0.500 where the observed value was something else. Nothing
        rewrites a scored row now: below the boundary it reads Pass at its own
        AU-ROC.
        """
        sequences_1 = random_sequences(MIN_SEQUENCES_PER_CLASS, 40, seed=24)
        sequences_2 = random_sequences(MIN_SEQUENCES_PER_CLASS, 40, seed=25)

        results, _ = _score_position_features(
            sequences_1, sequences_2, ['A'],
            'Per position nucleotide content', end_position=40,
            scored_end_position=40,
        )

        scores = [metrics['AU-ROC'] for name, metrics in results.items() if 'position' in name]
        assert all(0.5 <= score < 0.6 for score in scores), scores
        assert any(score > 0.5 for score in scores), 'every position came out at exactly chance'


class TestScoredWindow:
    """The window is where the required cohort still reaches.

    The requirement is one number made of two: a count that guards against
    sampling noise, and a share of the class that guards against a cohort made
    only of its longest sequences. Which one binds depends on the class size, and
    the tests below pin both regimes.
    """

    def test_the_count_binds_on_a_class_this_size(self):
        stats = make_stats(random_sequences(600, 20, seed=10))

        assert stats._required_cohort(600) == MIN_SEQUENCES_PER_CLASS

    def test_the_fraction_binds_once_the_class_is_large_enough(self):
        # A quarter of 2000 is 500 sequences, well past the count floor.
        stats = make_stats(['A' * 20] * 2000)

        assert stats._required_cohort(2000) == 500

    def test_the_window_ends_where_the_required_cohort_stops(self):
        # 300 sequences run to 100, the other 300 stop at 20, so exactly 300
        # reach anything past 20 and the floor of 250 is still met there.
        sequences = ['A' * 100] * 300 + ['A' * 20] * 300

        stats = make_stats(sequences)

        assert stats._required_cohort(600) == MIN_SEQUENCES_PER_CLASS
        assert stats.scored_end_position == 100

    def test_a_position_the_cohort_does_not_reach_is_left_out(self):
        # Only 200 sequences run past 20, which is short of the floor.
        sequences = ['A' * 100] * 200 + ['A' * 20] * 400

        stats = make_stats(sequences)

        assert stats.scored_end_position == 20

    def test_the_window_keeps_the_required_coverage(self):
        sequences = random_sequences(2000, 100, seed=10)
        sequences = [seq[:20] if i % 4 == 0 else seq for i, seq in enumerate(sequences)]

        stats = make_stats(sequences)

        assert stats.coverage_at(stats.scored_end_position) >= DEFAULT_MIN_COVERAGE

    def test_min_coverage_moves_the_window_and_leaves_the_plots_alone(self):
        sequences = ['A' * 100] * 700 + ['A' * 20] * 300

        strict = make_stats(sequences, min_coverage=0.9)
        loose = make_stats(sequences, min_coverage=0.5)

        assert strict.scored_end_position == 20
        assert loose.scored_end_position == 100
        assert strict.end_position == loose.end_position

    def test_min_coverage_cannot_go_below_the_count_floor(self):
        """`--min-coverage 0` asks for no share of the class, not for no floor:
        the count is the guard against sampling noise and is not negotiable."""
        stats = make_stats(['A' * 100] * 200 + ['A' * 20] * 400, min_coverage=0)

        assert stats._required_cohort(600) == MIN_SEQUENCES_PER_CLASS
        assert stats.scored_end_position == 20

    def test_a_class_below_the_count_floor_scores_nothing(self):
        stats = make_stats(['A' * 100] * (MIN_SEQUENCES_PER_CLASS - 1))

        assert stats.scored_end_position == 0
        assert stats.end_position == 100

    def test_too_few_sequences_to_score_still_gets_a_reported_window(self, caplog):
        sequences = ['A' * 100] * 10

        with caplog.at_level('WARNING'):
            stats = make_stats(sequences)

        assert stats.scored_end_position == 0
        assert stats.end_position == 100
        assert any('Not enough sequences' in record.message for record in caplog.records)


class TestCohortFloorForAComparison:
    """The cohort a comparison of two classes has to have behind a position.

    It is what ends the compared window the figures draw, so it is a number the
    window is derived from rather than something drawn on top of it.
    """

    def test_the_class_with_the_larger_share_sets_the_floor(self):
        """The same count is half of a 500-sequence class and a sixteenth of a
        4000-sequence one, so the smaller class sets the higher floor."""
        small = make_stats(['A' * 50] * 500, label='small')
        large = make_stats(['A' * 50] * 4000, label='large')

        assert cohort_floor(small, large) == pytest.approx(MIN_SEQUENCES_PER_CLASS / 500)

    def test_the_share_binds_once_the_class_is_large_enough(self):
        """Past 1000 sequences a quarter of the class is the larger of the two,
        so the floor stops being the sequence count and grows with the class."""
        stats1 = make_stats(['A' * 50] * 4000, label='a')
        stats2 = make_stats(['A' * 50] * 4000, label='b')

        assert cohort_floor(stats1, stats2) == pytest.approx(DEFAULT_MIN_COVERAGE)

    def test_an_empty_class_has_no_floor(self):
        empty = make_stats([], label='empty')

        assert cohort_floor(empty, empty) == 0.0


class TestPlotWindow:
    def test_default_ends_where_the_drawing_floor_does(self):
        # 60 sequences reach position 40; only 50 of them reach position 30.
        sequences = ['A' * 40] * 10 + ['A' * 30] * 40 + ['A' * 10] * 10

        stats = make_stats(sequences)

        assert stats.end_position == 30
        assert stats.coverage_at(30) == pytest.approx(50 / 60)

    def test_a_class_below_the_floor_reports_to_its_longest_sequence(self):
        """No position clears the floor, and stopping at 0 would leave the class
        with no per-position checks at all."""
        stats = make_stats(['A' * 30] * (MIN_SEQUENCES_PER_REPORTED_POSITION - 2) + ['A' * 10])

        assert stats.end_position == 30

    def test_uniform_lengths_use_the_whole_sequence(self):
        stats = make_stats(random_sequences(MIN_SEQUENCES_PER_CLASS, 60, seed=11))

        assert stats.end_position == 60
        assert stats.scored_end_position == 60
        assert stats.coverage_at(60) == 1.0

    def test_reaches_past_the_scored_window(self):
        """60 sequences run to 100, which is enough to draw and far short of the
        250 needed to compare."""
        sequences = random_sequences(600, 100, seed=15)
        sequences = [seq if i < 60 else seq[:30] for i, seq in enumerate(sequences)]

        stats = make_stats(sequences)

        assert stats.end_position == 100
        assert stats.scored_end_position == 30

    def test_explicit_end_position_is_capped_at_the_longest_sequence(self):
        stats = make_stats(random_sequences(20, 30, seed=12), end_position=500)

        assert stats.end_position == 30

    def test_explicit_end_position_below_the_cap_is_kept(self):
        stats = make_stats(random_sequences(20, 30, seed=13), end_position=10)

        assert stats.end_position == 10

    def test_an_explicit_window_shorter_than_the_scored_one_clamps_scoring(self):
        """A position that is not drawn has nothing to flag."""
        stats = make_stats(['A' * 100] * 600, end_position=10)

        assert stats.end_position == 10
        assert stats.scored_end_position == 10

    def test_an_explicit_window_cannot_widen_what_gets_scored(self):
        """The required cohort decides that, and no window argument moves it."""
        sequences = ['A' * 100] * 200 + ['A' * 20] * 400

        default = make_stats(sequences)
        widened = make_stats(sequences, end_position=100)

        assert widened.end_position == 100
        assert widened.scored_end_position == default.scored_end_position == 20

    def test_both_windows_are_logged_at_info(self, caplog):
        """The auto-chosen windows used to be DEBUG-only, invisible at the default level."""
        with caplog.at_level('INFO'):
            stats = make_stats(random_sequences(600, 45, seed=14))

        messages = [record.message for record in caplog.records]
        assert any(f'checks cover positions 1-{stats.end_position}' in message for message in messages)
        assert any(f'Positions 1-{stats.scored_end_position} may be flagged' in message
                   for message in messages)

    def test_the_reported_window_never_stops_before_the_scored_one(self):
        """A flag on a position no check is named for would have nowhere to go.
        The scored window asks for far more sequences than the reported one, so
        this holds by construction - it is pinned because the two floors are set
        apart."""
        for count in (MIN_SEQUENCES_PER_CLASS, 300, 600, 2000):
            sequences = random_sequences(count, 100, seed=40)
            sequences = [seq[:20] if i % 4 == 0 else seq for i, seq in enumerate(sequences)]

            stats = make_stats(sequences)

            assert stats.end_position >= stats.scored_end_position, count


class TestTheTailIsReportedButNotScored:
    def test_positions_past_the_scored_window_are_unknown(self):
        sequences = random_sequences(300, 40, seed=30)

        results, _ = _score_position_features(
            sequences, sequences, ['A'],
            'Per position nucleotide content', end_position=40,
            scored_end_position=20,
        )

        assert results['Per position nucleotide content - A position 20']['Flag'] != 'Unknown'
        assert results['Per position nucleotide content - A position 21']['Flag'] == 'Unknown'
        assert np.isnan(results['Per position nucleotide content - A position 40']['AU-ROC'])

    def test_the_tail_cannot_change_a_check_verdict(self):
        """The property the split rests on: widening the reported range moves no
        flag.

        The scored window, the number of tests the threshold corrects for, and
        every verdict inside the window have to come out the same whether the
        checks stop at the scored window or run to the end of the sequences.
        """
        sequences_1 = ['A' + seq[1:] for seq in random_sequences(300, 60, seed=31)]
        sequences_2 = ['C' + seq[1:] for seq in random_sequences(300, 60, seed=32)]
        sequences_2 = [seq if i % 3 else seq[:20] for i, seq in enumerate(sequences_2)]

        narrow, narrow_aggregate = _score_position_features(
            sequences_1, sequences_2, ['A', 'C', 'G', 'T'],
            'Per position nucleotide content', end_position=20, scored_end_position=20,
        )
        wide, wide_aggregate = _score_position_features(
            sequences_1, sequences_2, ['A', 'C', 'G', 'T'],
            'Per position nucleotide content', end_position=60, scored_end_position=20,
        )

        assert len(wide) > len(narrow)
        for name, metrics in narrow.items():
            assert wide[name] == pytest.approx(metrics, nan_ok=True), name
        for base, aggregate in narrow_aggregate.items():
            assert wide_aggregate[base] == pytest.approx(aggregate, nan_ok=True), base


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

    def test_the_curve_is_the_same_answer_at_every_point(self):
        """One question, one answer.

        The figure draws the curve, the interactive viewer carries it, and the
        prose quotes single points of it. These were three implementations, one
        of them a Python loop over the class per position; what a test can hold
        is that they cannot come apart again.
        """
        stats = make_stats(['A' * 7] * 3 + ['A' * 15] * 4 + ['A' * 40] * 2)

        curve = stats.coverage_curve(40)

        assert len(curve) == 40
        for position in range(1, 41):
            assert curve[position - 1] == pytest.approx(stats.coverage_at(position)), position

    def test_the_payload_carries_the_curve_the_figure_draws(self):
        """Rounded for JSON, and otherwise the same numbers."""
        stats = make_stats(['A' * 5] * 4 + ['A' * 12] * 6)

        curve = stats.coverage_curve(12)
        carried = _coverage(stats, 12)

        assert carried == [pytest.approx(round(float(value), 4)) for value in curve]

    def test_an_empty_class_has_a_curve_of_zeros_rather_than_no_curve(self):
        """The panel still gets drawn for a degenerate class."""
        assert list(make_stats([]).coverage_curve(3)) == [0.0, 0.0, 0.0]
