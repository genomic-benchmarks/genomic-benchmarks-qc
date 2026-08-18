"""Unit tests for the metric helpers in genomic_benchmarks_qc.utils.testing."""

import numpy as np
import pytest

from genomic_benchmarks_qc.utils.testing import _compute_best_threshold_accuracy


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
