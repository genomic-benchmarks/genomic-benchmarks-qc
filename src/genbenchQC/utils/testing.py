"""
Testing utilities for GenBenchQC.

This module provides functions to compute quality metrics (AU-ROC, AU-PR, Accuracy)
and generate flags for comparing two datasets.
"""

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, roc_auc_score, precision_recall_curve, auc

METRICS_TO_COMPUTE = ['AU-ROC', 'AU-PR', 'Accuracy']

def _compute_best_threshold_accuracy(labels: np.ndarray, scores: np.ndarray) -> float:
    """Find the threshold that maximizes accuracy on given labels/scores.

    Args:
        labels: Binary ground truth labels.
        scores: Continuous prediction scores.

    Returns:
        Best achievable accuracy (or prevalence-based baseline if single score).
    """
    unique_scores = np.unique(scores)
    if unique_scores.size == 1:
        return float(max(labels.mean(), 1 - labels.mean()))

    # Midpoints between consecutive unique scores + infinities
    thresholds = np.concatenate([
        [-np.inf],
        (unique_scores[:-1] + unique_scores[1:]) / 2,
        [np.inf],
    ])

    best_accuracy = 0.0
    for threshold in thresholds:
        predicted = scores >= threshold
        best_accuracy = max(best_accuracy, accuracy_score(labels, predicted))

    return float(best_accuracy)

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
        prevalence = values_1.size / (values_1.size + values_2.size)
        return {
            'AU-ROC': 0.5,
            'AU-PR': prevalence,
            'Accuracy': max(prevalence, 1 - prevalence),
        }

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


def _compute_flagged_metrics(values_1: np.ndarray, values_2: np.ndarray) -> dict:
    """Compute metrics and add flag based on AU-ROC.

    Args:
        values_1: Scores from dataset 1.
        values_2: Scores from dataset 2.

    Returns:
        Metrics dictionary with added 'Flag' key.
    """
    metrics = _compute_metrics_from_arrays(values_1, values_2)
    metrics['Flag'] = _flag_on_score(metrics['AU-ROC'])
    return metrics


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

    Args:
        sequences: List of DNA/protein sequences.
        base: Nucleotide/amino acid to check for.
        position: Position index (0-based).
        reverse: If True, count from end of sequence.

    Returns:
        Array of 1.0 where base matches at position, 0.0 otherwise.
    """
    values = np.zeros(len(sequences), dtype=float)

    for i, seq in enumerate(sequences):
        if reverse:
            index = len(seq) - 1 - position
            if index >= 0 and index < len(seq) and seq[index] == base:
                values[i] = 1.0
        elif position < len(seq) and seq[position] == base:
            values[i] = 1.0

    return values


def _score_position_features(
    sequences_1: list[str],
    sequences_2: list[str],
    bases: list[str],
    prefix: str,
    reverse: bool = False,
    end_position: int | None = None,
) -> tuple[dict, dict]:
    """Compute per-position metrics for all bases.

    Args:
        sequences_1: Sequences from dataset 1.
        sequences_2: Sequences from dataset 2.
        bases: Bases to analyze.
        prefix: Name prefix for result keys.
        reverse: If True, analyze from end of sequences.
        end_position: Max position to analyze (None = auto-detect).

    Returns:
        Tuple of (detailed results dict, per-base aggregates dict).
    """
    results = {}
    per_base_aggregates = {}

    if end_position is None:
        if not sequences_1 or not sequences_2:
            return results, per_base_aggregates
        end_position = min(
            max(len(s) for s in sequences_1),
            max(len(s) for s in sequences_2),
        )

    for base in bases:
        pos_metrics_list = []
        for position in range(end_position):
            vals1 = _compute_position_binary_scores(sequences_1, base, position, reverse)
            vals2 = _compute_position_binary_scores(sequences_2, base, position, reverse)
            metrics = _compute_metrics_from_arrays(vals1, vals2)
            metrics['Flag'] = _flag_on_score(metrics['AU-ROC'])
            result_name = f'{prefix} - {base} position {position + 1}'
            results[result_name] = metrics
            pos_metrics_list.append(metrics)

        # Aggregate across positions (worst-case)
        if pos_metrics_list:
            agg = _aggregate_worst_case_metrics(pos_metrics_list)
        else:
            agg = {metric: np.nan for metric in METRICS_TO_COMPUTE}
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
) -> dict:
    """Score a single scalar feature column between two datasets.

    Args:
        frame_1: DataFrame for dataset 1.
        frame_2: DataFrame for dataset 2.
        column_name: Column to score.
        indices_1: Row indices for dataset 1.
        indices_2: Row indices for dataset 2.

    Returns:
        Flagged metrics dictionary.
    """
    values_1 = frame_1.iloc[indices_1][column_name].fillna(0).to_numpy(dtype=float)
    values_2 = frame_2.iloc[indices_2][column_name].fillna(0).to_numpy(dtype=float)
    return _compute_flagged_metrics(values_1, values_2)


def _score_dataframe_features(
    frame_1: pd.DataFrame,
    frame_2: pd.DataFrame,
    prefix: str,
    indices_1: np.ndarray,
    indices_2: np.ndarray,
) -> dict:
    """Score all numeric columns in dataframes.

    Args:
        frame_1: DataFrame for dataset 1.
        frame_2: DataFrame for dataset 2.
        prefix: Name prefix for result keys.
        indices_1: Row indices for dataset 1.
        indices_2: Row indices for dataset 2.

    Returns:
        Dictionary with per-column and aggregate metrics.
    """
    results = {}
    columns = sorted(set(frame_1.columns) | set(frame_2.columns))

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
        reverse=False, end_position=end_position,
    )
    results.update(pos_results)
    if per_base_agg:
        results['Per position nucleotide content'] = _aggregate_worst_case_metrics(per_base_agg.values())

    # Position features (reverse)
    pos_results_rev, per_base_agg_rev = _score_position_features(
        stats1.sequences, stats2.sequences, bases,
        'Per reverse position nucleotide content',
        reverse=True, end_position=end_position,
    )
    results.update(pos_results_rev)
    if per_base_agg_rev:
        results['Per reverse position nucleotide content'] = _aggregate_worst_case_metrics(per_base_agg_rev.values())

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
            'Per reverse position nucleotide content': {...}
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
        'Per reverse position nucleotide content',
    ]

    all_results = {}

    all_results['Unique bases'] = {'Flag': _flag_unique_bases(stats1, stats2)}
    all_results['Sequence Duplications within Labels'] = _flag_duplicate_sequences(stats1, stats2)
    all_results['Duplicate Sequences between Labels'] = {
        'Flag': _flag_duplication_between_datasets(stats1.sequences, stats2.sequences)
    }

    model_results = direct_feature_model(stats1, stats2)
    all_results.update(model_results)

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


def _extract_failed_features(all_results: dict) -> dict:
    """Extract failure information organized by feature type for plotting.

    Args:
        all_results: Dictionary from direct_feature_model + manual flags.

    Returns:
        Nested dict structured as:
        {
            'Sequence lengths': {'Pass'},
            'Per sequence GC content': {'Warning'},
            'Per sequence nucleotide content': {'A': 'Warning', 'G': 'Fail', ...},
            'Per sequence dinucleotide content': {'AA': 'Pass', 'GG': 'Fail', ...},
            'Per position nucleotide content': {'A': {52: 'Warning'}, 'G': {66: 'Fail', 70: 'Fail'}, ...},
            'Per reverse position nucleotide content': {'A': {10: 'Fail'}, ...},
        }
    """
    failed_by_feature = {
        'Sequence lengths': {},
        'Per sequence GC content': {},
        'Sequence Duplications within Labels': {},
        'Per sequence nucleotide content': {},
        'Per sequence dinucleotide content': {},
        'Per position nucleotide content': {},
        'Per reverse position nucleotide content': {},
    }

    for key, value in all_results.items():
        flag = value.get('Flag', 'Unknown') if isinstance(value, dict) else 'Unknown'

        if key == 'Sequence lengths':
            if flag in ('Fail', 'Warning'):
                failed_by_feature['Sequence lengths'] = flag

        elif key == 'Sequence Duplications within Labels':
            if flag in ('Fail', 'Warning'):
                failed_by_feature['Sequence Duplications within Labels'] = flag

        elif key == 'Per sequence GC content':
            if flag in ('Fail', 'Warning'):
                failed_by_feature['Per sequence GC content'] = flag

        elif key.startswith('Per sequence nucleotide content - '):
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
                    if base not in failed_by_feature['Per position nucleotide content']:
                        failed_by_feature['Per position nucleotide content'][base] = {}
                    if flag in ('Fail', 'Warning'):
                        failed_by_feature['Per position nucleotide content'][base][position] = flag
                except ValueError:
                    # Not a position entry (e.g., aggregate "Per position nucleotide content - A")
                    pass

        elif key.startswith('Per reverse position nucleotide content - ') and ' position ' in key:
            # Parse "Per reverse position nucleotide content - G position 52"
            parts = key.rsplit(' position ', 1)
            if len(parts) == 2:
                base = parts[0].replace('Per reverse position nucleotide content - ', '')
                try:
                    position = int(parts[1])
                    if base not in failed_by_feature['Per reverse position nucleotide content']:
                        failed_by_feature['Per reverse position nucleotide content'][base] = {}
                    if flag in ('Fail', 'Warning'):
                        failed_by_feature['Per reverse position nucleotide content'][base][position] = flag
                except ValueError:
                    # Not a position entry (e.g., aggregate "Per reverse position nucleotide content - A")
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
