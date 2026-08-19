"""
Testing utilities for Genomic Benchmarks QC.

This module provides functions to compute quality metrics (AU-ROC, AU-PR, Accuracy)
and generate flags for comparing two datasets.

Every check here asks the same question: can a single feature tell the two
classes apart? The answer is an AU-ROC, and a fixed boundary on it decides the
flag. That boundary only means what it says when there are enough sequences
behind it, so two guards sit in front of the scoring:

- The per-sequence checks are not scored at all below
  `MIN_SEQUENCES_PER_CLASS` sequences in the smaller class, and report Unknown.
- The per-position checks reduce hundreds of positions to their worst case, and
  the number of sequences behind a position falls as the position grows, so
  they use a threshold that widens with the sampling noise of that position's
  own cohort instead of a size floor on the whole class.

Both rules were chosen by simulation against the class comparisons this tool
was built for; the study is in the companion manuscript's `analyses/` directory.

An underpowered comparison reports Unknown rather than Pass: no evidence of a
difference is not evidence of no difference.
"""

import logging

import numpy as np
import pandas as pd
from scipy.special import ndtri
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

METRICS_TO_COMPUTE = ['AU-ROC', 'AU-PR', 'Accuracy']

# Below this many sequences in the smaller class, a per-sequence check is not
# scored and reports Unknown. Under a null - two classes drawn from the same
# process - the worst of these checks, dinucleotide content, crosses the 0.6
# Warning boundary in 70.2% of replicates at 50 sequences per class and 19.4% at
# 100, against 1.2% at 200. Only 2.5% of the class comparisons in the benchmark
# audit this tool was built for fall below it.
MIN_SEQUENCES_PER_CLASS = 200

# Per-position checks score each position on the sequences that reach it, so
# their cohorts shrink along the sequence and a floor of MIN_SEQUENCES_PER_CLASS
# would silence the far end of every variable-length dataset. They use the
# adaptive threshold below instead, and this is only the point beneath which even
# that threshold is not trusted.
MIN_SEQUENCES_PER_POSITION = 50

# Family-wise error rate the per-position threshold targets across all
# (position, base) tests within one per-position check.
ADAPTIVE_ALPHA = 0.05

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
        return {metric: np.nan for metric in METRICS_TO_COMPUTE}

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

    precision, recall, _ = precision_recall_curve(labels, scores)
    aupr = auc(recall, precision)
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
    class. Used both for a feature that takes one value across both classes and
    for a position whose observed difference does not exceed sampling noise -
    in both cases the honest summary is "no information", not the number that
    happened to come out.

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
    return {metric: np.nan for metric in METRICS_TO_COMPUTE}


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

def _compute_position_binary_scores(sequences: list[str], base: str, position: int, reverse: bool) -> np.ndarray:
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


def _position_cohorts(sequences: list[str], end_position: int) -> np.ndarray:
    """Number of sequences reaching each 0-based position below `end_position`.

    A sequence has a position, read from either end, exactly when it is longer
    than the 0-based index, so one array of lengths answers this for the forward
    and the reverse pass alike.
    """
    lengths = np.fromiter((len(seq) for seq in sequences), dtype=int, count=len(sequences))
    return np.array([int(np.sum(lengths > position)) for position in range(end_position)])


def _adaptive_threshold(values_1: np.ndarray, values_2: np.ndarray, n_tests: int) -> float:
    """Smallest AU-ROC a position must reach before it may set a check's flag.

    A per-position check reduces every (position, base) pair to its worst case.
    Each pair is a deliberately weak test - for a binary indicator the AU-ROC
    this module computes is exactly 0.5 + |p1 - p2| / 2 - so the maximum over
    hundreds of them crosses a fixed boundary on sampling noise alone once the
    cohorts behind those positions are small. The boundary is therefore moved
    out to the largest difference sampling would be expected to produce across
    `n_tests` tests, at a family-wise error rate of `ADAPTIVE_ALPHA`:

        0.5 + z(1 - alpha / 2 n_tests) * SE / 2,

    with SE the pooled standard error of the difference in rates in this
    position's own cohorts. The correction is two-sided because AU-ROC is folded
    to at least 0.5, and the bound shrinks towards zero as cohorts grow - on a
    large dataset it falls below 0.6 and the check reduces to the fixed
    boundaries, differing only where the data are thin.

    Args:
        values_1: Binary indicator values from dataset 1 at this position.
        values_2: Binary indicator values from dataset 2 at this position.
        n_tests: Number of (position, base) tests the check aggregates.

    Returns:
        The AU-ROC a difference here must exceed to count as evidence.
    """
    size_1, size_2 = values_1.size, values_2.size
    pooled = (values_1.sum() + values_2.sum()) / (size_1 + size_2)
    standard_error = np.sqrt(pooled * (1 - pooled) * (1 / size_1 + 1 / size_2))
    critical = ndtri(1 - ADAPTIVE_ALPHA / (2 * n_tests))
    return float(0.5 + critical * standard_error / 2)


def _score_position_features(
    sequences_1: list[str],
    sequences_2: list[str],
    bases: list[str],
    prefix: str,
    end_position: int,
    reverse: bool = False,
    min_cohort: int = MIN_SEQUENCES_PER_POSITION,
) -> tuple[dict, dict]:
    """Compute per-position metrics for all bases.

    Each position is scored on the sequences that reach it, and two things then
    stand between a position and the check's verdict. A position whose cohort has
    fallen below `min_cohort` in either class is reported as Unknown rather than
    scored at all. A position that is scored must additionally show a difference
    larger than `_adaptive_threshold` allows for its own cohort size; below that
    it is reported at chance, because the alternative is to let the worst case
    over hundreds of weak tests be set by whichever of them was luckiest.

    Args:
        sequences_1: Sequences from dataset 1.
        sequences_2: Sequences from dataset 2.
        bases: Bases to analyze.
        prefix: Name prefix for result keys.
        end_position: Last position to analyze, 1-based and inclusive. Chosen by
            `SequenceStatistics._adjust_end_position` so that most sequences
            reach every position in the window.
        reverse: If True, analyze from end of sequences.
        min_cohort: Sequences a position needs in both classes to be scored.

    Returns:
        Tuple of (detailed results dict, per-base aggregates dict).
    """
    results = {}
    per_base_aggregates = {}

    # Which positions can be scored at all depends only on how many sequences
    # reach them, not on the base, so this is settled once - and it has to be,
    # because the threshold each position must clear depends on how many tests
    # the check ends up aggregating.
    cohorts_1 = _position_cohorts(sequences_1, end_position)
    cohorts_2 = _position_cohorts(sequences_2, end_position)
    scorable = np.minimum(cohorts_1, cohorts_2) >= max(min_cohort, 1)
    n_tests = int(scorable.sum()) * len(bases)

    for base in bases:
        pos_metrics_list = []
        for position in range(end_position):
            if not scorable[position]:
                metrics = _unknown_metrics()
            else:
                vals1 = _compute_position_binary_scores(sequences_1, base, position, reverse)
                vals2 = _compute_position_binary_scores(sequences_2, base, position, reverse)
                metrics = _compute_metrics_from_arrays(vals1, vals2)
                threshold = _adaptive_threshold(vals1, vals2, n_tests)
                if not metrics['AU-ROC'] > threshold:
                    metrics = _chance_metrics(vals1.size, vals2.size)
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

def direct_feature_model(stats1, stats2):
    """Compute all direct feature comparison metrics between two datasets.

    Uses full datasets (no subsampling) for raw-feature metrics.

    Args:
        stats1: SequenceStatistics object for dataset 1.
        stats2: SequenceStatistics object for dataset 2.

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
    end_position = min(stats1.end_position, stats2.end_position)

    pos_results, per_base_agg = _score_position_features(
        stats1.sequences, stats2.sequences, bases,
        'Per position nucleotide content',
        end_position=end_position, reverse=False,
    )
    results.update(pos_results)
    if per_base_agg:
        results['Per position nucleotide content'] = _aggregate_worst_case_metrics(per_base_agg.values())

    # Position features (reverse)
    pos_results_rev, per_base_agg_rev = _score_position_features(
        stats1.sequences, stats2.sequences, bases,
        'Per position reversed nucleotide content',
        end_position=end_position, reverse=True,
    )
    results.update(pos_results_rev)
    if per_base_agg_rev:
        results['Per position reversed nucleotide content'] = _aggregate_worst_case_metrics(per_base_agg_rev.values())

    return results


def flag_significant_differences(stats1, stats2):
    """Generate comprehensive QC comparison between two datasets.

    Args:
        stats1: SequenceStatistics object for dataset 1.
        stats2: SequenceStatistics object for dataset 2.

    Returns:
        Tuple of (summary_statuses, failed_by_feature) where:
        - summary_statuses: Ordered dictionary with all flags and metrics.
        - failed_by_feature: Nested dict with failure info for visualization:
          {
            'Per sequence nucleotide content': {'A': 'Warning', 'G': 'Fail', ...},
            'Per sequence dinucleotide content': {'AA': 'Pass', 'GG': 'Fail', ...},
            'Per position nucleotide content': {'A': {52: 'Warning'}, 'G': {66: 'Fail', 70: 'Fail'}, ...},
            'Per position reversed nucleotide content': {...}
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
            f"enough to reach it, and no position in the analysed window has at least "
            f"{MIN_SEQUENCES_PER_POSITION} sequences in both classes. The per-position "
            "plots are still produced."
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
            'Per position nucleotide content': {'A': {52: 'Warning'}, 'G': {66: 'Fail', 70: 'Fail'}, ...},
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
                        failed_by_feature['Per position nucleotide content'].setdefault(base, {})[position] = flag
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
                        failed_by_feature['Per position reversed nucleotide content'].setdefault(base, {})[position] = flag
                except ValueError:
                    # Not a position entry (e.g., aggregate "Per position reversed nucleotide content - A")
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
    elif score > 0.6:
        return 'Warning'
    else:
        return 'Pass'


def _flag_unique_bases(stats1, stats2) -> str:
    """Check if both datasets have identical base sets."""
    return 'Pass' if set(stats1.stats['Unique bases']) == set(stats2.stats['Unique bases']) else 'Fail'


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
