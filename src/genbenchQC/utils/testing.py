import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, roc_auc_score, precision_recall_curve, auc

METRICS_TO_COMPUTE = ['AU-ROC', 'AU-PR', 'Accuracy']

def _score_position_base(sequences, base, position, reverse):
    values = np.zeros(len(sequences), dtype=float)

    for i, seq in enumerate(sequences):
        if reverse:
            index = len(seq) - 1 - position
            if index >= 0 and index < len(seq) and seq[index] == base:
                values[i] = 1.0
        elif position < len(seq) and seq[position] == base:
            values[i] = 1.0

    return values

def _best_threshold_accuracy(labels, scores):
    unique_scores = np.unique(scores)
    if unique_scores.size == 1:
        return float(max(labels.mean(), 1 - labels.mean()))

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

def _metric_bundle(values_1, values_2):
    values_1 = np.asarray(values_1, dtype=float)
    values_2 = np.asarray(values_2, dtype=float)

    if values_1.size == 0 or values_2.size == 0:
        return {metric: np.nan for metric in METRICS_TO_COMPUTE}

    combined = np.concatenate([values_1, values_2])
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
    if auroc < 0.5:
        scores = -scores
        auroc = 1 - auroc

    precision, recall, _ = precision_recall_curve(labels, scores)
    aupr = auc(recall, precision)
    accuracy = _best_threshold_accuracy(labels, scores)

    return {
        'AU-ROC': float(auroc),
        'AU-PR': float(aupr),
        'Accuracy': float(accuracy),
    }

def _flagged_metrics(values_1, values_2):
    metrics = _metric_bundle(values_1, values_2)
    metrics['Flag'] = flag_on_score(metrics['AU-ROC'])
    return metrics

def _worst_case_metrics(metric_dicts):
    aggregate = {}

    for metric_name in METRICS_TO_COMPUTE:
        values = np.asarray([metrics.get(metric_name, np.nan) for metrics in metric_dicts], dtype=float)
        if values.size == 0 or np.all(np.isnan(values)):
            aggregate[metric_name] = np.nan
        else:
            aggregate[metric_name] = float(np.nanmax(values))

    aggregate['Flag'] = flag_on_score(aggregate.get('AU-ROC', np.nan))
    return aggregate

def _score_scalar_feature(frame_1, frame_2, column_name, indices_1, indices_2):
    values_1 = frame_1.iloc[indices_1][column_name].fillna(0).to_numpy(dtype=float)
    values_2 = frame_2.iloc[indices_2][column_name].fillna(0).to_numpy(dtype=float)
    return _flagged_metrics(values_1, values_2)

def _score_dataframe_features(frame_1, frame_2, prefix, indices_1, indices_2):
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

        results[f'{prefix} - {column}'] = _flagged_metrics(values_1, values_2)

    # compute worst-case aggregate across columns (max over metrics)
    if columns:
        results[f'{prefix}'] = _worst_case_metrics(
            [results[f'{prefix} - {column}'] for column in columns]
        )

    return results

def _score_position_features(sequences_1, sequences_2, bases, prefix, reverse=False, end_position=None):
    results = {}
    per_base_aggregates = {}

    if end_position is None:
        if not sequences_1 or not sequences_2:
            return results, per_base_aggregates
        end_position = min(max(len(sequence) for sequence in sequences_1), max(len(sequence) for sequence in sequences_2))

    for base in bases:
        # collect per-position raw metric dicts to compute per-base aggregate
        pos_metrics_list = []
        for position in range(end_position):
            vals1 = _score_position_base(sequences_1, base, position, reverse)
            vals2 = _score_position_base(sequences_2, base, position, reverse)
            metrics = _metric_bundle(vals1, vals2)
            metrics['Flag'] = flag_on_score(metrics['AU-ROC'])
            result_name = f'{prefix} - {base} position {position + 1}'
            results[result_name] = metrics
            pos_metrics_list.append(metrics)

        # compute per-base aggregated metrics (worst-case across positions)
        if pos_metrics_list:
            agg = _worst_case_metrics(pos_metrics_list)
        else:
            agg = {metric: np.nan for metric in METRICS_TO_COMPUTE}
            agg['Flag'] = flag_on_score(np.nan)

        results[f'{prefix} - {base}'] = agg
        per_base_aggregates[base] = agg

    return results, per_base_aggregates

def direct_feature_model(stats1, stats2):
    # Use full datasets (no balancing/subsampling) when computing raw-feature metrics.
    indices_1 = np.arange(len(stats1.sequences))
    indices_2 = np.arange(len(stats2.sequences))

    results = {}

    results['Sequence lengths'] = _score_scalar_feature(
        stats1.stats['Sequence lengths'],
        stats2.stats['Sequence lengths'],
        'Sequence lengths',
        indices_1,
        indices_2,
    )

    results['Per sequence GC content'] = _score_scalar_feature(
        stats1.stats['Per sequence GC content'],
        stats2.stats['Per sequence GC content'],
        'Per sequence GC content',
        indices_1,
        indices_2,
    )

    results.update(_score_dataframe_features(
        stats1.stats['Per sequence nucleotide content'],
        stats2.stats['Per sequence nucleotide content'],
        'Per sequence nucleotide content',
        indices_1,
        indices_2,
    ))

    results.update(_score_dataframe_features(
        stats1.stats['Per sequence dinucleotide content'],
        stats2.stats['Per sequence dinucleotide content'],
        'Per sequence dinucleotide content',
        indices_1,
        indices_2,
    ))

    bases = sorted(list(set(stats1.stats['Unique bases']) | set(stats2.stats['Unique bases'])))
    end_position = min(stats1.end_position, stats2.end_position)

    pos_results, per_base_agg = _score_position_features(
        stats1.sequences,
        stats2.sequences,
        bases,
        'Per position nucleotide content',
        reverse=False,
        end_position=end_position,
    )
    results.update(pos_results)

    # add overall worst-case across bases for forward positions
    if per_base_agg:
        results['Per position nucleotide content'] = _worst_case_metrics(per_base_agg.values())

    pos_results_rev, per_base_agg_rev = _score_position_features(
        stats1.sequences,
        stats2.sequences,
        bases,
        'Per reverse position nucleotide content',
        reverse=True,
        end_position=end_position,
    )
    results.update(pos_results_rev)

    # add overall worst-case across bases for reverse positions
    if per_base_agg_rev:
        results['Per reverse position nucleotide content'] = _worst_case_metrics(per_base_agg_rev.values())

    return results

def flag_significant_differences(stats1, stats2):

    results = {}

    # Define the desired order of statistics
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

    # Compute all results
    all_results = {}

    all_results['Unique bases'] = {}
    all_results['Unique bases']['Flag'] = flag_unique_bases(stats1, stats2)

    all_results['Sequence Duplications within Labels'] = {}
    all_results['Sequence Duplications within Labels']['Flag'] = flag_duplicate_sequences(stats1, stats2)

    all_results['Duplicate Sequences between Labels'] = {}
    all_results['Duplicate Sequences between Labels']['Flag'] = flag_duplication_between_datasets(
        stats1.sequences, stats2.sequences
    )

    model_results = direct_feature_model(stats1, stats2)
    all_results.update(model_results)

    # Reorganize: add aggregates in desired order first, then all details
    for stat_name in ordered_stats:
        if stat_name in all_results:
            results[stat_name] = all_results[stat_name]
    
    # Add all details for this statistic
    for stat_name in ordered_stats:
        for key in all_results:
            if key.startswith(f"{stat_name} - ") and key != stat_name:
                results[key] = all_results[key]

    return results

def flag_on_score(score):
    if score is None or not np.isfinite(score):
        raise ValueError(f"Score must be a finite number, got {score}")
    if score > 0.7:
        return "Fail"
    elif score > 0.6:
        return "Warning"
    else:
        return "Pass"

def flag_unique_bases(stats1, stats2):
    if set(stats1.stats['Unique bases']) == set(stats2.stats['Unique bases']):
        return 'Pass'
    else:
        return 'Fail'
    
def flag_duplicate_sequences(stats1, stats2):
    if stats1.stats['Number of sequences'] != stats1.stats['Number of sequences left after deduplication']:
        return 'Warning'
    if stats2.stats['Number of sequences'] != stats2.stats['Number of sequences left after deduplication']:
        return 'Warning'
    return 'Pass'

def flag_duplication_between_datasets(sequences1, sequences2):
    return "Fail" if bool(set(sequences1) & set(sequences2)) else "Pass"
