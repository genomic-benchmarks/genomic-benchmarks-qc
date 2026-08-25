"""Unit tests for the metric helpers in genomic_benchmarks_qc.utils.testing."""

import numpy as np
import pytest

from genomic_benchmarks_qc.utils.testing import (
    METRICS_TO_COMPUTE,
    _binary_feature_metrics,
    _compute_best_threshold_accuracy,
    _compute_metrics_from_arrays,
    _compute_position_binary_scores,
    _position_base_counts,
    _score_position_features,
)


def brute_force_best_accuracy(labels, scores):
    """Reference implementation of the metric, written straight from its definition.

    Scans every threshold that can produce a distinct set of predictions -- the
    midpoints between consecutive unique scores, plus both infinities -- and
    returns the best accuracy any of them achieves. Deliberately naive: it exists
    to pin down the specification independently of how the real function is
    optimised.
    """
    labels = np.asarray(labels)
    scores = np.asarray(scores, dtype=float)
    unique_scores = np.unique(scores)

    thresholds = np.concatenate([
        [-np.inf],
        (unique_scores[:-1] + unique_scores[1:]) / 2,
        [np.inf],
    ])

    return max(float((labels == (scores >= t)).mean()) for t in thresholds)


class TestBestThresholdAccuracy:
    @pytest.mark.parametrize(
        "labels, scores, expected",
        [
            # Perfectly separable, either orientation of the scores.
            ([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], 1.0),
            ([1, 1, 0, 0], [0.9, 0.8, 0.2, 0.1], 1.0),
            # A single sample is always classifiable.
            ([1], [5.0], 1.0),
            ([0], [5.0], 1.0),
            # One misplaced sample out of four.
            ([0, 1, 0, 1], [0.1, 0.2, 0.8, 0.9], 0.75),
            # Interleaved labels: no threshold beats the prevalence baseline.
            ([0, 1, 0, 1], [1.0, 2.0, 3.0, 4.0], 0.75),
        ],
    )
    def test_known_cases(self, labels, scores, expected):
        assert _compute_best_threshold_accuracy(
            np.asarray(labels), np.asarray(scores, dtype=float)
        ) == pytest.approx(expected)

    @pytest.mark.parametrize(
        "labels, scores, expected",
        [
            # The regression case: a threshold cannot separate the two samples
            # scoring 1.0, so the best available accuracy is 2/3. Splitting that
            # tied group would wrongly report 1.0.
            ([0, 0, 1], [0.0, 1.0, 1.0], 2 / 3),
            # Two samples with opposite labels and identical scores: 0.5 is the
            # ceiling, and tie-splitting would claim 1.0.
            ([1, 0], [5.0, 5.0], 0.5),
            # A tied group straddling the decision boundary in the middle.
            ([0, 0, 1, 0, 1, 1], [0.0, 1.0, 1.0, 1.0, 1.0, 2.0], 4 / 6),
        ],
    )
    def test_tied_scores_are_never_split(self, labels, scores, expected):
        assert _compute_best_threshold_accuracy(
            np.asarray(labels), np.asarray(scores, dtype=float)
        ) == pytest.approx(expected)

    @pytest.mark.parametrize(
        "labels, expected",
        [
            # With no score variation, only "all positive" and "all negative"
            # are available, so the answer is the prevalence baseline.
            ([1, 1, 1, 0], 0.75),
            ([0, 0, 0, 1], 0.75),
            ([0, 1], 0.5),
            ([1, 1, 1], 1.0),
        ],
    )
    def test_constant_scores_fall_back_to_prevalence(self, labels, expected):
        scores = np.zeros(len(labels), dtype=float)
        assert _compute_best_threshold_accuracy(
            np.asarray(labels), scores
        ) == pytest.approx(expected)

    def test_constant_nonzero_scores(self):
        """The constant-score shortcut must not assume the constant is zero."""
        labels = np.array([1, 1, 0])
        scores = np.full(3, -7.5)
        assert _compute_best_threshold_accuracy(labels, scores) == pytest.approx(2 / 3)

    def test_matches_brute_force_on_random_inputs(self):
        """Randomised comparison against the reference implementation."""
        rng = np.random.default_rng(20250818)

        for _ in range(500):
            n_samples = int(rng.integers(1, 15))
            labels = rng.integers(0, 2, n_samples)
            # A small score alphabet, so ties -- the interesting case -- are common.
            scores = rng.integers(0, 4, n_samples).astype(float)

            assert _compute_best_threshold_accuracy(labels, scores) == pytest.approx(
                brute_force_best_accuracy(labels, scores)
            ), f"disagreement on labels={labels.tolist()} scores={scores.tolist()}"

    def test_matches_brute_force_on_continuous_scores(self):
        """Same comparison without ties, covering the fully distinct-score path."""
        rng = np.random.default_rng(7)

        for _ in range(200):
            n_samples = int(rng.integers(2, 30))
            labels = rng.integers(0, 2, n_samples)
            scores = rng.normal(size=n_samples)

            assert _compute_best_threshold_accuracy(labels, scores) == pytest.approx(
                brute_force_best_accuracy(labels, scores)
            )

    def test_never_below_prevalence_baseline(self):
        """Predicting the majority class everywhere is always available."""
        rng = np.random.default_rng(11)

        for _ in range(200):
            n_samples = int(rng.integers(1, 20))
            labels = rng.integers(0, 2, n_samples)
            scores = rng.integers(0, 3, n_samples).astype(float)

            prevalence = labels.mean()
            baseline = max(prevalence, 1 - prevalence)
            assert _compute_best_threshold_accuracy(labels, scores) >= baseline - 1e-12

    def test_invariant_to_input_order(self):
        """The metric depends on the label/score pairing, not on row order."""
        rng = np.random.default_rng(3)
        labels = rng.integers(0, 2, 40)
        scores = rng.integers(0, 5, 40).astype(float)

        expected = _compute_best_threshold_accuracy(labels, scores)
        for _ in range(10):
            order = rng.permutation(labels.size)
            assert _compute_best_threshold_accuracy(
                labels[order], scores[order]
            ) == pytest.approx(expected)

    def test_invariant_to_increasing_rescaling(self):
        """Only the ranking of the scores matters, not their scale or offset."""
        labels = np.array([0, 0, 1, 1, 0, 1])
        scores = np.array([0.1, 0.4, 0.4, 0.9, 0.2, 0.7])

        expected = _compute_best_threshold_accuracy(labels, scores)
        for rescaled in (scores * 1000, scores * 3 - 5, np.exp(scores)):
            assert _compute_best_threshold_accuracy(
                labels, rescaled
            ) == pytest.approx(expected)

    def test_handles_negative_and_large_scores(self):
        labels = np.array([0, 0, 1, 1])
        scores = np.array([-1e6, -3.5, 2.5, 1e6])
        assert _compute_best_threshold_accuracy(labels, scores) == pytest.approx(1.0)

    def test_accepts_python_lists(self):
        """Callers pass array-likes; the function converts them itself."""
        assert _compute_best_threshold_accuracy(
            [0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9]
        ) == pytest.approx(1.0)

    def test_returns_plain_float(self):
        result = _compute_best_threshold_accuracy(
            np.array([0, 1]), np.array([0.0, 1.0])
        )
        assert type(result) is float

    def test_result_is_in_unit_interval(self):
        rng = np.random.default_rng(99)

        for _ in range(100):
            n_samples = int(rng.integers(1, 25))
            labels = rng.integers(0, 2, n_samples)
            scores = rng.integers(0, 4, n_samples).astype(float)

            result = _compute_best_threshold_accuracy(labels, scores)
            assert 0.0 <= result <= 1.0


class TestClosedFormBinaryMetrics:
    """The closed forms against the score arrays they replaced.

    `_binary_feature_metrics` reads three metrics off a 2x2 table that
    `_compute_metrics_from_arrays` gets by building one score per sequence and
    handing it to sklearn. The second is the definition; these assert the first
    is the same answer, which is the only thing that licenses the speed.
    """

    @staticmethod
    def by_the_definition(matches_1, cohort_1, matches_2, cohort_2):
        """Score the table the long way: one array entry per sequence."""
        values_1 = np.concatenate([np.ones(matches_1), np.zeros(cohort_1 - matches_1)])
        values_2 = np.concatenate([np.ones(matches_2), np.zeros(cohort_2 - matches_2)])
        return _compute_metrics_from_arrays(values_1, values_2)

    @staticmethod
    def closed_form(matches_1, cohort_1, matches_2, cohort_2):
        """Score the same table from its four counts."""
        auroc, aupr, accuracy = _binary_feature_metrics(
            np.array([matches_1]), np.array([cohort_1]),
            np.array([matches_2]), np.array([cohort_2]))
        return {'AU-ROC': auroc[0], 'AU-PR': aupr[0], 'Accuracy': accuracy[0]}

    def test_every_small_table_agrees(self):
        """Exhaustively, so the corners are covered rather than sampled.

        Both cohorts from 1 to 8, every count of matches in each: the base absent
        from one class, present in all of the other, equally common in both, and
        every table in between.
        """
        worst = dict.fromkeys(METRICS_TO_COMPUTE, 0.0)
        tables = 0

        for cohort_1 in range(1, 9):
            for cohort_2 in range(1, 9):
                for matches_1 in range(cohort_1 + 1):
                    for matches_2 in range(cohort_2 + 1):
                        fast = self.closed_form(matches_1, cohort_1, matches_2, cohort_2)
                        slow = self.by_the_definition(matches_1, cohort_1, matches_2, cohort_2)
                        tables += 1
                        for metric in METRICS_TO_COMPUTE:
                            worst[metric] = max(worst[metric], abs(fast[metric] - slow[metric]))

        assert tables == 1936
        # AU-PR comes out bit-identical. The other two differ by at most an ulp,
        # and in both cases it is the array path that rounds: sklearn's AU-ROC
        # accumulates a sum, and `_chance_metrics` reaches the accuracy of the
        # larger class through `1 - prevalence` rather than by dividing.
        assert worst['AU-PR'] == 0.0
        assert worst['Accuracy'] < 1e-12
        assert worst['AU-ROC'] < 1e-12

    def test_larger_tables_agree(self):
        """Cohort sizes around the ones a real per-position check runs at."""
        rng = np.random.default_rng(20260825)

        for _ in range(300):
            cohort_1 = int(rng.integers(1, 3000))
            cohort_2 = int(rng.integers(1, 3000))
            matches_1 = int(rng.integers(0, cohort_1 + 1))
            matches_2 = int(rng.integers(0, cohort_2 + 1))

            fast = self.closed_form(matches_1, cohort_1, matches_2, cohort_2)
            slow = self.by_the_definition(matches_1, cohort_1, matches_2, cohort_2)
            for metric in METRICS_TO_COMPUTE:
                assert fast[metric] == pytest.approx(slow[metric], abs=1e-12), (
                    f"{metric} at {matches_1}/{cohort_1} vs {matches_2}/{cohort_2}")

    @pytest.mark.parametrize("matches_1, matches_2", [(0, 0), (40, 60), (20, 30)])
    def test_a_base_that_says_nothing_scores_chance(self, matches_1, matches_2):
        """Absent from both, or equally common in both, is no information.

        These reach `_chance_metrics` in the array path only by accident of which
        branch catches them - the first because there is one score value, the
        others through the general computation. The closed form has no branch for
        any of them and lands on chance regardless, which is the check.
        """
        cohort_1, cohort_2 = 40, 60
        expected_prevalence = cohort_1 / (cohort_1 + cohort_2)

        fast = self.closed_form(matches_1, cohort_1, matches_2, cohort_2)

        assert fast['AU-ROC'] == pytest.approx(0.5)
        assert fast['AU-PR'] == pytest.approx(expected_prevalence)
        assert fast['Accuracy'] == pytest.approx(max(expected_prevalence,
                                                     1 - expected_prevalence))


class TestPositionScoringAgreesWithTheDefinition:
    """`_score_position_features` against the per-sequence score arrays.

    The closed forms are checked above on tables made up by hand; this checks the
    counting that produces those tables from real sequences, in both directions
    and with sequences that stop at different points.
    """

    @staticmethod
    def sequences(count, rng, min_length, max_length):
        lengths = rng.integers(min_length, max_length + 1, size=count)
        return [''.join(rng.choice(list('ACGT'), size=length)) for length in lengths]

    @pytest.mark.parametrize("reverse", [False, True])
    @pytest.mark.parametrize("min_length, max_length", [(30, 30), (12, 30)])
    def test_every_scored_position_matches(self, reverse, min_length, max_length):
        rng = np.random.default_rng(7)
        bases = list('ACGT')
        sequences_1 = self.sequences(120, rng, min_length, max_length)
        sequences_2 = self.sequences(90, rng, min_length, max_length)
        end_position = 30

        results, _ = _score_position_features(
            sequences_1, sequences_2, bases, 'Prefix',
            end_position=end_position, scored_end_position=end_position,
            reverse=reverse, min_cohort=10,
        )

        scored = 0
        for base in bases:
            for position in range(end_position):
                metrics = results[f'Prefix - {base} position {position + 1}']
                if np.isnan(metrics['AU-ROC']):
                    continue
                scored += 1
                values_1 = _compute_position_binary_scores(
                    sequences_1, base, position, reverse)
                values_2 = _compute_position_binary_scores(
                    sequences_2, base, position, reverse)
                expected = _compute_metrics_from_arrays(values_1, values_2)
                for metric in METRICS_TO_COMPUTE:
                    assert metrics[metric] == pytest.approx(expected[metric], abs=1e-12), (
                        f"{metric} for {base} at position {position + 1}")

        # A guard on the guard: a run that scored nothing would pass vacuously.
        assert scored > 40

    def test_the_counts_are_of_the_sequences_that_reach_the_position(self):
        """A sequence that stops early contributes to neither count.

        Read in reverse this is the whole subtlety of the check: position 1 is
        the last base of every sequence whatever its length, so the cohorts stay
        full while the bases at them come from different offsets.
        """
        sequences = ['AAAA', 'AAC', 'AG']
        counts = _position_base_counts(sequences, ['A', 'C', 'G'], 4, reverse=False)

        assert list(counts[0]) == [3, 2, 1, 1]   # A: all three, then 'AG' turns to G
        assert list(counts[1]) == [0, 0, 1, 0]   # C, only in 'AAC'
        assert list(counts[2]) == [0, 1, 0, 0]   # G, only in 'AG'

        reversed_counts = _position_base_counts(sequences, ['A', 'C', 'G'], 4, reverse=True)

        assert list(reversed_counts[0]) == [1, 3, 2, 1]  # only 'AAAA' ends in A
        assert list(reversed_counts[1]) == [1, 0, 0, 0]  # 'AAC' ends in C
        assert list(reversed_counts[2]) == [1, 0, 0, 0]  # 'AG' ends in G
