"""
Testing utilities for Genomic Benchmarks QC.

This module provides functions to compute quality metrics (AU-ROC, AU-PR, Accuracy)
and generate flags for comparing two datasets.

Every check here asks the same question: can a single feature tell the two
classes apart? The answer is an AU-ROC, and a fixed boundary on it decides the
flag. That boundary only means what it says when there are enough sequences
behind it, and when those sequences are the class rather than a corner of it, so
two guards sit in front of the scoring:

- Nothing is scored on fewer than `MIN_SEQUENCES_PER_CLASS` sequences. For the
  per-sequence checks that is a floor on the smaller class; for the per-position
  checks it is a floor on the cohort reaching each position, since a position is
  compared only on the sequences long enough to have it.
- The per-position checks additionally stop where the cohort reaching a position
  falls below `seq_stats.DEFAULT_MIN_COVERAGE` of its class, well before the
  plots stop drawing it. A cohort far out along the sequence can be large and
  still not stand for the class - it is all of the class's long sequences, and a
  difference there can be a difference between those subsets rather than between
  the classes. Size does not bound that; only stopping does.

Together the two put the same requirement on every comparison: at least
`MIN_SEQUENCES_PER_CLASS` sequences, and at least `DEFAULT_MIN_COVERAGE` of the
class. The count binds on small and mid-sized classes, where sampling noise is
the risk; the fraction binds on large ones, where a tail cohort can clear the
count many times over and still describe only the longest sequences.

The size floor was chosen by simulation against the class comparisons this tool
was built for; the study is in `analyses/`. The coverage floor is not a question
about power and was not simulated - no sample size makes a subset of the class
stand for the class.

An underpowered comparison reports Unknown rather than Pass: no evidence of a
difference is not evidence of no difference.
"""

import logging
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

if TYPE_CHECKING:
    # Imported for annotations only: seq_stats imports this module, so a real
    # import would be circular.
    from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics

METRICS_TO_COMPUTE = ['AU-ROC', 'AU-PR', 'Accuracy']

# Sequences a comparison needs before it is made at all: in the smaller class
# for the per-sequence checks, and in each position's cohort for the per-position
# ones. Below it a check reports Unknown.
#
# Under a null - two classes drawn from the same process - the worst per-sequence
# check, dinucleotide content, crosses the 0.6 Warning boundary in 70.2% of
# replicates at 50 sequences per class and 19.4% at 100, against 0.2% at 250 and
# 0.9% at 200. The per-position checks reach 0.0% at 250 across every class size
# simulated, which is what lets them share one fixed boundary with the
# per-sequence checks instead of correcting for the hundreds of tests they
# aggregate. `analyses/` holds both studies.
MIN_SEQUENCES_PER_CLASS = 250

def _compute_best_threshold_accuracy(labels: np.ndarray, scores: np.ndarray) -> float:
    """Find the threshold that maximizes accuracy on given labels/scores.

    Args:
        labels: Binary ground truth labels.
        scores: Continuous prediction scores.

    Returns:
        Best achievable accuracy (or prevalence-based baseline if single score).
    """
    labels = np.asarray(labels)
    scores = np.asarray(scores, dtype=float)
    n_samples = scores.size

    # No score variation: no threshold separates anything, so the best accuracy
    # is the prevalence baseline. An O(n) check, cheaper than the sort below.
    if n_samples > 0 and (scores == scores[0]).all():
        return float(max(labels.mean(), 1 - labels.mean()))

    # Candidate thresholds are the boundaries between groups of equal scores:
    # a single sorted pass evaluates all of them in O(n log n).
    order = np.argsort(scores, kind='mergesort')
    sorted_scores = scores[order]
    sorted_labels = labels[order].astype(np.int64)

    # Positives among the lowest `cut` scores, for cut = 0 .. n_samples.
    positives_below = np.concatenate([[0], np.cumsum(sorted_labels)])
    total_positives = positives_below[-1]
    cuts = np.arange(n_samples + 1)

    # Predicting positive for everything at or above the cut: true positives
    # above the cut plus true negatives below it.
    correct = (total_positives - positives_below) + (cuts - positives_below)

    # Only cuts that do not split a group of tied scores are reachable, plus the
    # two extremes (predict all positive / all negative).
    reachable = np.ones(n_samples + 1, dtype=bool)
    if n_samples > 1:
        reachable[1:-1] = sorted_scores[:-1] < sorted_scores[1:]

    return float(correct[reachable].max() / n_samples)

def _compute_metrics_from_arrays(values_1: np.ndarray, values_2: np.ndarray) -> dict:
    """Compute AU-ROC, AU-PR, and Accuracy from two arrays of scores.

    Args:
        values_1: Scores from dataset 1 (treated as positive class).
        values_2: Scores from dataset 2 (treated as negative class).

    Returns:
        Dictionary with AU-ROC, AU-PR, Accuracy, and Flag.
    """
    values_1 = np.asarray(values_1, dtype=float)
    values_2 = np.asarray(values_2, dtype=float)

    # Empty input case
    if values_1.size == 0 or values_2.size == 0:
        return dict.fromkeys(METRICS_TO_COMPUTE, np.nan)

    combined = np.concatenate([values_1, values_2])

    # Single unique value case
    if np.unique(combined).size < 2:
        return _chance_metrics(values_1.size, values_2.size)

    labels = np.concatenate([
        np.ones(values_1.size, dtype=int),
        np.zeros(values_2.size, dtype=int),
    ])
    scores = combined

    auroc = roc_auc_score(labels, scores)

    # Invert scores if AU-ROC < 0.5 (indicates reverse correlation)
    if auroc < 0.5:
        scores = -scores
        auroc = 1 - auroc

    # Average precision, not the trapezoidal area under the precision-recall
    # curve. They differ by how the gap between two curve points is filled, and
    # trapezoidal interpolation fills it with a straight line drawn from a point
    # no classifier achieves - which for the per-position checks, where the score
    # takes two values and the curve has two points, is most of the curve. A
    # position whose base is equally common in both classes scores 0.625 that
    # way, next to an AU-ROC of 0.500. Average precision gives it 0.500, the same
    # answer `_chance_metrics` gives a feature that says nothing.
    aupr = average_precision_score(labels, scores)
    accuracy = _compute_best_threshold_accuracy(labels, scores)

    return {
        'AU-ROC': float(auroc),
        'AU-PR': float(aupr),
        'Accuracy': float(accuracy),
    }


def _chance_metrics(size_1: int, size_2: int) -> dict:
    """Metrics for a feature that says nothing about which class a sequence is in.

    The values a perfectly uninformative feature would produce: chance AU-ROC,
    an AU-PR at the prevalence, and the accuracy of always guessing the larger
    class. Used for a feature that takes one value across both classes, where the
    honest summary is "no information" rather than whatever the degenerate ROC
    computation would return.

    Args:
        size_1: Number of values from dataset 1.
        size_2: Number of values from dataset 2.

    Returns:
        Dictionary with AU-ROC, AU-PR, and Accuracy.
    """
    prevalence = size_1 / (size_1 + size_2)
    return {
        'AU-ROC': 0.5,
        'AU-PR': prevalence,
        'Accuracy': max(prevalence, 1 - prevalence),
    }


def _unknown_metrics() -> dict:
    """Metrics for a comparison that was not made, which flag as Unknown."""
    return dict.fromkeys(METRICS_TO_COMPUTE, np.nan)


def _flag_metrics(metrics: dict) -> dict:
    """Add the 'Flag' key implied by a metrics dict's AU-ROC, in place."""
    metrics['Flag'] = _flag_on_score(metrics['AU-ROC'])
    return metrics


def _compute_flagged_metrics(values_1: np.ndarray, values_2: np.ndarray) -> dict:
    """Compute metrics and add flag based on AU-ROC.

    Args:
        values_1: Scores from dataset 1.
        values_2: Scores from dataset 2.

    Returns:
        Metrics dictionary with added 'Flag' key.
    """
    return _flag_metrics(_compute_metrics_from_arrays(values_1, values_2))


def _aggregate_worst_case_metrics(metric_dicts: list[dict]) -> dict:
    """Aggregate multiple metric dicts by taking worst-case (max) per metric.

    Args:
        metric_dicts: List of metric dictionaries.

    Returns:
        Aggregated dictionary with max values across inputs.
    """
    aggregate = {}

    for metric_name in METRICS_TO_COMPUTE:
        values = np.asarray([m.get(metric_name, np.nan) for m in metric_dicts], dtype=float)
        if values.size == 0 or np.all(np.isnan(values)):
            aggregate[metric_name] = np.nan
        else:
            aggregate[metric_name] = float(np.nanmax(values))

    aggregate['Flag'] = _flag_on_score(aggregate.get('AU-ROC', np.nan))
    return aggregate

def _compute_position_binary_scores(
    sequences: list[str], base: str, position: int, reverse: bool
) -> np.ndarray:
    """Compute binary scores for a specific base at a position across sequences.

    Only the sequences long enough to have this position contribute: a sequence
    that ends before it is left out rather than scored as "not this base". That
    keeps the comparison a question about composition -- P(base at position |
    the sequence reaches the position) -- instead of one that also answers
    itself from sequence length, which `Sequence lengths` already scores as a
    feature of its own. It is also the definition the per-position plots use,
    since `SequenceStatistics._normalize_per_position` normalizes each position
    by its own total.

    The returned array is therefore shorter than `sequences` wherever some
    sequences stop before `position`, and how much shorter is what decides
    whether the position can be scored at all and how large a difference it has
    to show; `_score_position_features` applies both rules.

    Args:
        sequences: List of DNA/protein sequences.
        base: Nucleotide/amino acid to check for.
        position: Position index (0-based).
        reverse: If True, count from end of sequence.

    This is the definition, not the route production takes. For a score with two
    values the metrics are functions of the 2x2 table alone, so
    `_score_position_features` counts the table and closes the form; building
    these arrays walks the class once per base and position. It stays because a
    definition worth trusting is worth being able to run: `test_metrics.py` scores
    both ways and asserts they agree.

    Returns:
        Array holding 1.0 where the base matches and 0.0 where it does not, with
        one entry per sequence that reaches `position`.
    """
    values = []

    for seq in sequences:
        # Reading from either end, a sequence has this position exactly when it
        # is longer than the 0-based index.
        if len(seq) <= position:
            continue
        index = len(seq) - 1 - position if reverse else position
        values.append(1.0 if seq[index] == base else 0.0)

    return np.asarray(values, dtype=float)


# How much of a character matrix `_position_base_counts` holds at a time. NumPy's
# fixed-width text is 4 bytes per character, so this bounds the working set no
# matter how many sequences there are.
POSITION_COUNT_CHUNK_BYTES = 64 * 1024 * 1024


def _position_base_counts(
    sequences: list[str], bases: list[str], end_position: int, reverse: bool
) -> np.ndarray:
    """Count how many sequences carry each base at each position.

    One pass over the sequences answers this for every base and every position at
    once. That is what makes the per-position checks affordable: asking per
    (base, position) instead walks the whole class `len(bases) * end_position`
    times over, which was 82% of a run.

    A sequence that stops before a position does not count towards it, the same
    rule `_compute_position_binary_scores` follows - short sequences pad with
    NUL here, and NUL is not a base.

    Args:
        sequences: The class's sequences.
        bases: Bases to count, one character each. Anything longer cannot sit at
            a single position and counts zero, which is what comparing it
            against one character already did.
        end_position: How many positions to count.
        reverse: If True, position 0 is the last character of each sequence.

    Returns:
        Integer array of shape (len(bases), end_position).
    """
    counts = np.zeros((len(bases), end_position), dtype=np.int64)
    if end_position <= 0 or not sequences:
        return counts

    wanted = [(index, ord(base)) for index, base in enumerate(bases) if len(base) == 1]
    if not wanted:
        return counts

    rows_per_chunk = max(1, POSITION_COUNT_CHUNK_BYTES // (4 * end_position))
    for start in range(0, len(sequences), rows_per_chunk):
        block = sequences[start:start + rows_per_chunk]
        if reverse:
            windows = [seq[max(0, len(seq) - end_position):][::-1] for seq in block]
        else:
            windows = [seq[:end_position] for seq in block]
        # Fixed-width text viewed as integers is a rectangle of code points, NUL
        # where a sequence had already ended.
        codes = np.array(windows, dtype=f'U{end_position}').view(np.uint32)
        codes = codes.reshape(len(windows), end_position)
        for index, code in wanted:
            counts[index] += np.count_nonzero(codes == code, axis=0)

    return counts


def _binary_feature_metrics(matches_1, cohorts_1, matches_2, cohorts_2):
    """AU-ROC, AU-PR and best-threshold accuracy for a 0/1 feature, in closed form.

    A per-position check gives every sequence one of two scores - it has this
    base here, or it does not - and for a two-valued score all three metrics are
    functions of the 2x2 table. Writing p1 and p2 for the fraction of each class
    carrying the base:

    - AU-ROC is `(1 + |p1 - p2|) / 2`. Every cross-class pair of sequences either
      ties, counting a half, or separates the one way the difference in rates
      allows.
    - Average precision is one term per point of the precision-recall curve, and
      two distinct scores make two points.
    - The best threshold is one of three cuts: everything positive, everything
      negative, or the boundary between the two scores.

    The degenerate cases need no branch of their own. A base that never occurs,
    one that always occurs, and one equally common in both classes all come out
    at `_chance_metrics` - AU-ROC 0.5, AU-PR at the prevalence - because that is
    what the formulas say, not because they are checked for.

    Args:
        matches_1: Sequences of class 1 carrying the base, per position.
        cohorts_1: Sequences of class 1 reaching the position, per position.
        matches_2: As `matches_1`, for class 2.
        cohorts_2: As `cohorts_1`, for class 2.

    Returns:
        Tuple of (AU-ROC, AU-PR, Accuracy) arrays over the broadcast inputs.
        Every cohort must be non-empty; the caller decides which positions have
        enough sequences to score at all.
    """
    matches_1 = np.asarray(matches_1, dtype=float)
    matches_2 = np.asarray(matches_2, dtype=float)
    cohorts_1 = np.asarray(cohorts_1, dtype=float)
    cohorts_2 = np.asarray(cohorts_2, dtype=float)

    total = cohorts_1 + cohorts_2
    prevalence = cohorts_1 / total
    rate_1 = matches_1 / cohorts_1
    rate_2 = matches_2 / cohorts_2

    # `_compute_metrics_from_arrays` inverts a score that runs the wrong way, so
    # that a difference either way reads the same. Inverted, "carries the base"
    # becomes "does not", which is this swap.
    inverted = rate_1 < rate_2
    high_1 = np.where(inverted, cohorts_1 - matches_1, matches_1)
    high_2 = np.where(inverted, cohorts_2 - matches_2, matches_2)

    auroc = 0.5 + 0.5 * np.abs(rate_1 - rate_2)

    # The curve's two points are the cut keeping only the higher score and the
    # cut keeping everything; average precision weights each point's precision by
    # the recall it adds. Nothing scores high only when the base is absent from
    # both classes, and there the second point is the whole curve.
    scored = high_1 + high_2
    recall = high_1 / cohorts_1
    precision = np.divide(high_1, scored, out=np.zeros_like(scored), where=scored > 0)
    aupr = np.where(scored > 0,
                    precision * recall + prevalence * (1.0 - recall),
                    prevalence)

    accuracy = np.maximum(np.maximum(cohorts_1, cohorts_2),
                          high_1 + cohorts_2 - high_2) / total

    return auroc, aupr, accuracy


def _position_cohorts(sequences: list[str], end_position: int) -> np.ndarray:
    """Number of sequences reaching each 0-based position below `end_position`.

    A sequence has a position, read from either end, exactly when it is longer
    than the 0-based index, so one array of lengths answers this for the forward
    and the reverse pass alike.
    """
    lengths = np.sort(np.fromiter((len(seq) for seq in sequences), dtype=int, count=len(sequences)))
    # Sorted lengths turn every count into a boundary lookup: everything past
    # the last length that stops at or before a position still reaches it.
    positions = np.arange(end_position)
    return len(lengths) - np.searchsorted(lengths, positions, side='right')


def position_windows(stats1, stats2) -> tuple[int, int]:
    """The per-position windows a comparison of two classes runs in.

    Each class resolves its own windows from its own sequence lengths, and a
    position belongs to a window only where it belongs to it in both: a position
    that half of one class does not reach cannot be flagged on the strength of
    the other class reaching it.

    Returns:
        Tuple of (end_position, scored_end_position), 1-based and inclusive - the
        last position reported on and the last position allowed to set a flag,
        which is also the last position drawn.
    """
    return (min(stats1.end_position, stats2.end_position),
            min(stats1.scored_end_position, stats2.scored_end_position))


def _score_position_features(
    sequences_1: list[str],
    sequences_2: list[str],
    bases: list[str],
    prefix: str,
    end_position: int,
    scored_end_position: int,
    reverse: bool = False,
    min_cohort: int = MIN_SEQUENCES_PER_CLASS,
) -> tuple[dict, dict]:
    """Compute per-position metrics for all bases.

    Every position up to `end_position` gets a row, because a position the
    report is silent about is indistinguishable from one the sequences never
    reach, but only the positions up to `scored_end_position` are compared: past
    it too small a fraction of each class reaches the position for a difference
    there to be about the position rather than about the longest sequences, and
    the rows say Unknown.

    A position whose cohort has fallen below `min_cohort` in either class is
    reported as Unknown as well. `SequenceStatistics` already ends the scored
    window where that happens, so within a window it derived this is a redundant
    guard; it matters when the window was narrowed by hand, and it keeps the
    function honest on its own terms.

    Args:
        sequences_1: Sequences from dataset 1.
        sequences_2: Sequences from dataset 2.
        bases: Bases to analyze.
        prefix: Name prefix for result keys.
        end_position: Last position to report on, 1-based and inclusive. Chosen
            by `SequenceStatistics._reported_window`.
        scored_end_position: Last position allowed to set a flag, 1-based and
            inclusive. Chosen by `SequenceStatistics._scored_window`, and the
            window the per-position figures draw.
        reverse: If True, analyze from end of sequences.
        min_cohort: Sequences a position needs in both classes to be scored.

    Returns:
        Tuple of (detailed results dict, per-base aggregates dict).
    """
    results = {}
    per_base_aggregates = {}

    # Which positions can be scored at all depends only on how many sequences
    # reach them, not on the base, so this is settled once. Both gates are
    # monotone in the position, so the scorable positions are always a prefix of
    # the plotted ones.
    cohorts_1 = _position_cohorts(sequences_1, end_position)
    cohorts_2 = _position_cohorts(sequences_2, end_position)
    within_window = np.arange(end_position) < scored_end_position
    scorable = within_window & (np.minimum(cohorts_1, cohorts_2) >= max(min_cohort, 1))

    # Every scored position, for every base, in one pass over each class and one
    # vectorised evaluation of the closed forms in `_binary_feature_metrics`.
    # Counting stops at the last position anything is scored at - the rows past
    # it read Unknown whatever the sequences hold there.
    scored_positions = np.flatnonzero(scorable)
    counted_end = int(scored_positions[-1]) + 1 if scored_positions.size else 0
    matches_1 = _position_base_counts(sequences_1, bases, counted_end, reverse)
    matches_2 = _position_base_counts(sequences_2, bases, counted_end, reverse)
    auroc, aupr, accuracy = _binary_feature_metrics(
        matches_1[:, scored_positions], cohorts_1[scored_positions],
        matches_2[:, scored_positions], cohorts_2[scored_positions],
    )
    # Which column of those arrays each scorable position landed in.
    column_of = np.full(end_position, -1, dtype=int)
    column_of[scored_positions] = np.arange(scored_positions.size)

    for base_index, base in enumerate(bases):
        pos_metrics_list = []
        for position in range(end_position):
            if not scorable[position]:
                metrics = _unknown_metrics()
            else:
                column = column_of[position]
                metrics = {
                    'AU-ROC': float(auroc[base_index, column]),
                    'AU-PR': float(aupr[base_index, column]),
                    'Accuracy': float(accuracy[base_index, column]),
                }
            metrics['Flag'] = _flag_on_score(metrics['AU-ROC'])
            result_name = f'{prefix} - {base} position {position + 1}'
            results[result_name] = metrics
            pos_metrics_list.append(metrics)

        # Aggregate across positions (worst-case)
        if pos_metrics_list:
            agg = _aggregate_worst_case_metrics(pos_metrics_list)
        else:
            agg = _unknown_metrics()
            agg['Flag'] = _flag_on_score(np.nan)

        results[f'{prefix} - {base}'] = agg
        per_base_aggregates[base] = agg

    return results, per_base_aggregates

def _score_scalar_feature(
    frame_1: pd.DataFrame,
    frame_2: pd.DataFrame,
    column_name: str,
    indices_1: np.ndarray,
    indices_2: np.ndarray,
    min_class_size: int = MIN_SEQUENCES_PER_CLASS,
) -> dict:
    """Score a single scalar feature column between two datasets.

    Args:
        frame_1: DataFrame for dataset 1.
        frame_2: DataFrame for dataset 2.
        column_name: Column to score.
        indices_1: Row indices for dataset 1.
        indices_2: Row indices for dataset 2.
        min_class_size: Sequences the smaller class needs before the comparison
            is made at all; below it the feature reports Unknown. Pass 0 to score
            regardless of size.

    Returns:
        Flagged metrics dictionary.
    """
    if min(len(indices_1), len(indices_2)) < min_class_size:
        return _flag_metrics(_unknown_metrics())
    values_1 = frame_1.iloc[indices_1][column_name].fillna(0).to_numpy(dtype=float)
    values_2 = frame_2.iloc[indices_2][column_name].fillna(0).to_numpy(dtype=float)
    return _compute_flagged_metrics(values_1, values_2)


def _score_dataframe_features(
    frame_1: pd.DataFrame,
    frame_2: pd.DataFrame,
    prefix: str,
    indices_1: np.ndarray,
    indices_2: np.ndarray,
    min_class_size: int = MIN_SEQUENCES_PER_CLASS,
) -> dict:
    """Score all numeric columns in dataframes.

    Args:
        frame_1: DataFrame for dataset 1.
        frame_2: DataFrame for dataset 2.
        prefix: Name prefix for result keys.
        indices_1: Row indices for dataset 1.
        indices_2: Row indices for dataset 2.
        min_class_size: Sequences the smaller class needs before the comparison
            is made at all; below it every column, and the aggregate over them,
            reports Unknown. Pass 0 to score regardless of size.

    Returns:
        Dictionary with per-column and aggregate metrics.
    """
    results = {}
    columns = sorted(set(frame_1.columns) | set(frame_2.columns))

    # The worst case over several columns is more exposed to sampling noise than
    # any one of them, so a class too small for one column is too small for the
    # aggregate as well: the whole check reports Unknown together.
    if min(len(indices_1), len(indices_2)) < min_class_size:
        for column in columns:
            results[f'{prefix} - {column}'] = _flag_metrics(_unknown_metrics())
        if columns:
            results[prefix] = _flag_metrics(_unknown_metrics())
        return results

    for column in columns:
        if column in frame_1.columns:
            values_1 = frame_1.iloc[indices_1][column].fillna(0).to_numpy(dtype=float)
        else:
            values_1 = np.zeros(len(indices_1), dtype=float)

        if column in frame_2.columns:
            values_2 = frame_2.iloc[indices_2][column].fillna(0).to_numpy(dtype=float)
        else:
            values_2 = np.zeros(len(indices_2), dtype=float)

        results[f'{prefix} - {column}'] = _compute_flagged_metrics(values_1, values_2)

    # Compute worst-case aggregate across columns
    if columns:
        results[prefix] = _aggregate_worst_case_metrics(
            [results[f'{prefix} - {col}'] for col in columns]
        )

    return results

def direct_feature_model(stats1: 'SequenceStatistics',
                         stats2: 'SequenceStatistics') -> dict:
    """Compute all direct feature comparison metrics between two datasets.

    Uses full datasets (no subsampling) for raw-feature metrics.

    Args:
        stats1: Statistics for the first dataset.
        stats2: Statistics for the second dataset.

    Returns:
        Nested dictionary with all computed metrics and flags.
    """
    indices_1 = np.arange(len(stats1.sequences))
    indices_2 = np.arange(len(stats2.sequences))

    results = {}

    # Scalar features
    results['Sequence lengths'] = _score_scalar_feature(
        stats1.stats['Sequence lengths'],
        stats2.stats['Sequence lengths'],
        'Sequence lengths',
        indices_1, indices_2,
    )

    results['Per sequence GC content'] = _score_scalar_feature(
        stats1.stats['Per sequence GC content'],
        stats2.stats['Per sequence GC content'],
        'Per sequence GC content',
        indices_1, indices_2,
    )

    # DataFrame features
    results.update(_score_dataframe_features(
        stats1.stats['Per sequence nucleotide content'],
        stats2.stats['Per sequence nucleotide content'],
        'Per sequence nucleotide content',
        indices_1, indices_2,
    ))

    results.update(_score_dataframe_features(
        stats1.stats['Per sequence dinucleotide content'],
        stats2.stats['Per sequence dinucleotide content'],
        'Per sequence dinucleotide content',
        indices_1, indices_2,
    ))

    # Position features (forward)
    bases = sorted(set(stats1.stats['Unique bases']) | set(stats2.stats['Unique bases']))
    end_position, scored_end_position = position_windows(stats1, stats2)

    pos_results, per_base_agg = _score_position_features(
        stats1.sequences, stats2.sequences, bases,
        'Per position nucleotide content',
        end_position=end_position, scored_end_position=scored_end_position,
        reverse=False,
    )
    results.update(pos_results)
    if per_base_agg:
        results['Per position nucleotide content'] = _aggregate_worst_case_metrics(
            per_base_agg.values())

    # Position features (reverse)
    pos_results_rev, per_base_agg_rev = _score_position_features(
        stats1.sequences, stats2.sequences, bases,
        'Per position reversed nucleotide content',
        end_position=end_position, scored_end_position=scored_end_position,
        reverse=True,
    )
    results.update(pos_results_rev)
    if per_base_agg_rev:
        results['Per position reversed nucleotide content'] = (
            _aggregate_worst_case_metrics(per_base_agg_rev.values()))

    return results


def flag_significant_differences(stats1: 'SequenceStatistics',
                                 stats2: 'SequenceStatistics') -> tuple[dict, dict]:
    """Generate comprehensive QC comparison between two datasets.

    Args:
        stats1: Statistics for the first dataset.
        stats2: Statistics for the second dataset.

    Returns:
        Tuple of `(summary_statuses, failed_by_feature)`, where
        `summary_statuses` is an ordered dictionary of every flag and metric,
        and `failed_by_feature` carries the per-feature detail the plots need:

            {
                'Per sequence nucleotide content': {'A': 'Warning', 'G': 'Fail', ...},
                'Per sequence dinucleotide content': {'AA': 'Pass', 'GG': 'Fail', ...},
                'Per position nucleotide content':
                    {'A': {52: 'Warning'}, 'G': {66: 'Fail', 70: 'Fail'}, ...},
                'Per position reversed nucleotide content': {...},
            }
    """
    results = {}

    ordered_stats = [
        'Unique bases',
        'Sequence Duplications within Labels',
        'Duplicate Sequences between Labels',
        'Sequence lengths',
        'Per sequence GC content',
        'Per sequence nucleotide content',
        'Per sequence dinucleotide content',
        'Per position nucleotide content',
        'Per position reversed nucleotide content',
    ]

    all_results = {}

    all_results['Unique bases'] = {'Flag': _flag_unique_bases(stats1, stats2)}
    all_results['Sequence Duplications within Labels'] = _flag_duplicate_sequences(stats1, stats2)
    all_results['Duplicate Sequences between Labels'] = {
        'Flag': _flag_duplication_between_datasets(stats1.sequences, stats2.sequences)
    }

    model_results = direct_feature_model(stats1, stats2)
    all_results.update(model_results)

    _warn_about_unscored_checks(stats1, stats2, all_results, ordered_stats)

    # Order: aggregates first, then details
    for stat_name in ordered_stats:
        if stat_name in all_results:
            results[stat_name] = all_results[stat_name]

    for stat_name in ordered_stats:
        for key in all_results:
            if key.startswith(f"{stat_name} - ") and key != stat_name:
                results[key] = all_results[key]

    # Build failed_by_feature dict for visualization
    failed_by_feature = _extract_failed_features(all_results)

    return results, failed_by_feature


def _warn_about_unscored_checks(stats1, stats2, all_results: dict, check_names: list[str]):
    """Say on the terminal which checks were not scored, and why.

    A check that reports Unknown looks much like one that reports Pass in a
    sidebar full of icons, and the difference matters: Unknown means the
    comparison was not made, not that it came out clean. Anyone running the tool
    on a small dataset should hear that from the terminal rather than having to
    infer it from the report.

    Args:
        stats1, stats2: SequenceStatistics objects for the two classes.
        all_results: Every computed check, keyed by name.
        check_names: The top-level check names, in report order.
    """
    unknown = [name for name in check_names
               if isinstance(all_results.get(name), dict)
               and all_results[name].get('Flag') == 'Unknown']
    if not unknown:
        return

    positional = [name for name in unknown if 'position' in name.lower()]
    per_sequence = [name for name in unknown if name not in positional]

    label_1 = stats1.label if stats1.label is not None else stats1.filename
    label_2 = stats2.label if stats2.label is not None else stats2.filename
    size_1 = stats1.stats['Number of sequences']
    size_2 = stats2.stats['Number of sequences']

    if per_sequence:
        logging.warning(
            f"Not enough sequences to score {len(per_sequence)} check(s): "
            f"{', '.join(per_sequence)}. '{label_1}' has {size_1:,} sequences and "
            f"'{label_2}' has {size_2:,}, and below {MIN_SEQUENCES_PER_CLASS} in the "
            "smaller class these checks flag a difference on sampling noise alone too "
            "often to be trusted, so they are reported as Unknown rather than Pass. "
            "The plots and the descriptive statistics are still produced, so the "
            "distributions can be compared by eye."
        )

    if positional:
        logging.warning(
            f"Not enough sequences to score {len(positional)} check(s): "
            f"{', '.join(positional)}. Each position is compared on the sequences long "
            f"enough to reach it, and no position has at least "
            f"{MIN_SEQUENCES_PER_CLASS} sequences in both classes. The figures are still "
            "drawn, over every position the checks are named for, so the frequencies can "
            "be compared by eye - but nothing in them carries a flag."
        )


def _extract_failed_features(all_results: dict) -> dict:
    """Extract failure information organized by feature type for plotting.

    Args:
        all_results: Dictionary from direct_feature_model + manual flags.

    Returns:
        Nested dict structured as:
        {
            'Per sequence nucleotide content': {'A': 'Warning', 'G': 'Fail', ...},
            'Per sequence dinucleotide content': {'AA': 'Pass', 'GG': 'Fail', ...},
            'Per position nucleotide content':
                {'A': {52: 'Warning'}, 'G': {66: 'Fail', 70: 'Fail'}, ...},
            'Per position reversed nucleotide content': {'A': {10: 'Fail'}, ...},
        }
    """
    failed_by_feature = {
        'Per sequence nucleotide content': {},
        'Per sequence dinucleotide content': {},
        'Per position nucleotide content': {},
        'Per position reversed nucleotide content': {},
    }

    for key, value in all_results.items():
        flag = value.get('Flag', 'Unknown') if isinstance(value, dict) else 'Unknown'

        if key.startswith('Per sequence nucleotide content - '):
            nucleotide = key.replace('Per sequence nucleotide content - ', '')
            if flag in ('Fail', 'Warning'):
                failed_by_feature['Per sequence nucleotide content'][nucleotide] = flag

        elif key.startswith('Per sequence dinucleotide content - '):
            dinucleotide = key.replace('Per sequence dinucleotide content - ', '')
            if flag in ('Fail', 'Warning'):
                failed_by_feature['Per sequence dinucleotide content'][dinucleotide] = flag

        elif key.startswith('Per position nucleotide content - ') and ' position ' in key:
            # Parse "Per position nucleotide content - G position 52"
            # Split from the right to handle "Per position" in the prefix
            parts = key.rsplit(' position ', 1)
            if len(parts) == 2:
                base = parts[0].replace('Per position nucleotide content - ', '')
                try:
                    position = int(parts[1])
                    # Only flagged positions go in. An entry per base regardless
                    # would leave the dict truthy for a comparison with nothing
                    # to shade, and the report would draw a second, identical
                    # copy of every per-position plot.
                    if flag in ('Fail', 'Warning'):
                        per_position = failed_by_feature['Per position nucleotide content']
                        per_position.setdefault(base, {})[position] = flag
                except ValueError:
                    # Not a position entry (e.g., aggregate "Per position nucleotide content - A")
                    pass

        elif key.startswith('Per position reversed nucleotide content - ') and ' position ' in key:
            # Parse "Per position reversed nucleotide content - G position 52"
            parts = key.rsplit(' position ', 1)
            if len(parts) == 2:
                base = parts[0].replace('Per position reversed nucleotide content - ', '')
                try:
                    position = int(parts[1])
                    # Only flagged positions go in. An entry per base regardless
                    # would leave the dict truthy for a comparison with nothing
                    # to shade, and the report would draw a second, identical
                    # copy of every per-position plot.
                    if flag in ('Fail', 'Warning'):
                        per_position = failed_by_feature[
                            'Per position reversed nucleotide content']
                        per_position.setdefault(base, {})[position] = flag
                except ValueError:
                    # Not a position entry (e.g., the aggregate
                    # "Per position reversed nucleotide content - A")
                    pass

    return failed_by_feature


def _flag_on_score(score: float) -> str:
    """Assign Pass/Warning/Fail flag based on AU-ROC score.

    Args:
        score: AU-ROC value (must be finite).

    Returns:
        'Pass' if score <= 0.6, 'Warning' if <= 0.7, 'Fail' otherwise.
        Returns 'Unknown' for non-finite scores instead of raising.
    """
    if score is None or not np.isfinite(score):
        return 'Unknown'
    if score > 0.7:
        return 'Fail'
    if score > 0.6:
        return 'Warning'
    return 'Pass'


def _flag_unique_bases(stats1, stats2) -> str:
    """Check if both datasets have identical base sets."""
    same_bases = set(stats1.stats['Unique bases']) == set(stats2.stats['Unique bases'])
    return 'Pass' if same_bases else 'Fail'


def _flag_duplicate_sequences(stats1, stats2) -> dict:
    """Check for duplicates within each dataset.

    Computes combined deduplication ratio across both datasets.
    Returns:
        Dict with 'Flag' and 'percent_remaining' keys.
        Flag is 'Fail' if < 98% sequences remain, 'Warning' if >= 98% but < 100%, 'Pass' otherwise.
    """
    total_sequences = 0
    unique_sequences = 0
    for stats in [stats1, stats2]:
        total_sequences += stats.stats['Number of sequences']
        unique_sequences += stats.stats['Number of sequences left after deduplication']

    percent_remaining = unique_sequences / total_sequences if total_sequences > 0 else 1.0

    if percent_remaining < 0.98:
        flag = 'Fail'
    elif percent_remaining < 1.0:
        flag = 'Warning'
    else:
        flag = 'Pass'

    return {'Flag': flag, 'Percent Remaining': percent_remaining}


def _flag_duplication_between_datasets(sequences1: list[str], sequences2: list[str]) -> str:
    """Check for overlapping sequences between datasets."""
    return "Fail" if bool(set(sequences1) & set(sequences2)) else "Pass"
