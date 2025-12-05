import logging
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score, precision_recall_curve, auc
import numpy as np
import pandas as pd
from typing import List, Tuple


def extract_per_position_base(sequences, base, reverse, end_position=None):
    if not sequences:  
        return np.zeros((0, 0))
    if end_position is not None:
        sequences = [seq[:end_position] for seq in sequences]
        max_len = end_position
    else:
        max_len = max([len(seq) for seq in sequences])
    features = np.zeros((len(sequences), max_len))
    for i, seq in enumerate(sequences):
        if reverse:
            seq = seq[::-1]
        for j, nt in enumerate(seq):
            features[i, j] = 1 if nt == base else 0
    return features


STATS_TO_TRAIN_PRECOMPUTED = [
    'Sequence lengths',
    'Per sequence GC content',
    'Per sequence nucleotide content',
    'Per sequence dinucleotide content'
]

METRICS_TO_COMPUTE = ['AU-ROC', 'AU-PR', 'Accuracy']

def flag_on_score(score):
    if score > 0.7:
        return "Fail"
    elif score > 0.6:
        return "Warning"
    else:
        return "Pass"

def model(stats1, stats2, max_class_size=None, metric_to_flag='AU-ROC'):

    if metric_to_flag not in METRICS_TO_COMPUTE:
        raise ValueError(f"metric_to_flag must be one of {METRICS_TO_COMPUTE}, got {metric_to_flag}")

    results = {}

    for stat in STATS_TO_TRAIN_PRECOMPUTED:

        logging.info(f"Training bias detection model for statistic: {stat}")

        X = pd.concat([stats1.stats[stat], stats2.stats[stat]], axis=0)
        X.fillna(0, inplace=True)

        if stat == 'Sequence lengths':
            X = np.log1p(X)

        y = pd.Series([1] * len(stats1.stats[stat]) + [0] * len(stats2.stats[stat]))

        metrics = cross_validation(X, y, cv=5, max_size=max_class_size)

        logging.debug(f"{metric_to_flag} scores for {stat}: {metrics[metric_to_flag]}")
        results = add_result(results, stat, metrics, metric_to_flag)

    common_nts = list(set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases']))

    for reverse in [False, True]:
        # Compute aggregate metrics for per position nucleotides - take worst case
        worse_metrics = {metric: [] for metric in METRICS_TO_COMPUTE}
        for nt in common_nts:
            log_msg = "Training bias detection model for per position nucleotide" if not reverse else "Training bias detection model for per reverse position nucleotide"
            logging.info(f"{log_msg}: {nt}")

            end_position=min(stats1.end_position, stats2.end_position)
            features = extract_per_position_base(stats1.sequences + stats2.sequences, base=nt, reverse=reverse, end_position=end_position)

            X = pd.DataFrame(features)
            y = pd.Series([1] * len(stats1.sequences) + [0] * len(stats2.sequences))

            metrics = cross_validation(X, y, cv=5, max_size=max_class_size)

            flag_name = f"Per position nucleotide content - {nt}" if not reverse else f"Per reverse position nucleotide content - {nt}"
            logging.debug(f"{metric_to_flag} scores for {flag_name}: {metrics[metric_to_flag]}")
            results = add_result(results, flag_name, metrics, metric_to_flag)

            for metric_name in METRICS_TO_COMPUTE:
                worse_metrics[metric_name].append(metrics[metric_name].mean())

        # Add worst case metrics
        worst_case_metrics = {metric: np.array([max(worse_metrics[metric])]) for metric in METRICS_TO_COMPUTE}
        flag_name = "Per position nucleotide content" if not reverse else "Per reverse position nucleotide content"
        results = add_result(results, flag_name, worst_case_metrics, metric_to_flag)

    return results

def add_result(results, key, metrics, metric_to_flag):
    results[key] = {}
    for metric_name, scores in metrics.items():
        avg_score = scores.mean()
        results[key][metric_name] = avg_score
    results[key]['Flag'] = flag_on_score(metrics[metric_to_flag].mean())

    return results

def train_model(X, y, use_dual=False, C=1.0):

    model = LogisticRegression(random_state=42, solver='liblinear', dual=use_dual, max_iter=200, C=C)
    model.fit(X, y)

    return model

def eval_model(model, X, y):
    y_prob = model.predict_proba(X)[:, 1]
    metrics = {}

    metrics['AU-ROC'] = roc_auc_score(y, y_prob)

    precision, recall, _ = precision_recall_curve(y, y_prob)
    metrics['AU-PR'] = auc(recall, precision)

    y_pred = y_prob >= 0.5
    metrics['Accuracy'] = accuracy_score(y, y_pred)

    return metrics

def balanced_kfold_splits(y, cv=5, max_size=None) -> List[Tuple[np.ndarray, np.ndarray]]:
   
    # Create custom balanced k-fold splits
    unique_labels = np.unique(y)
    min_class_size = min([sum(y == label) for label in unique_labels])
    if max_size is not None:
        min_class_size = min(min_class_size, max_size)

    balanced_indices = []
    for label in unique_labels:
        label_indices = np.where(y == label)[0]
        np.random.seed(42)  # For reproducibility
        # Subsample to match minimum class size
        if len(label_indices) > min_class_size:
            label_indices = np.random.choice(label_indices, min_class_size, replace=False)
        balanced_indices.extend(label_indices)

    # Convert to numpy array and shuffle
    balanced_indices = np.array(balanced_indices)
    np.random.shuffle(balanced_indices)

    # Create k folds
    fold_size = len(balanced_indices) // cv
    folds = []
    for i in range(cv):
        start_idx = i * fold_size
        end_idx = start_idx + fold_size if i < cv - 1 else len(balanced_indices)
        val_idx = balanced_indices[start_idx:end_idx]
        train_idx = np.array([idx for idx in balanced_indices if idx not in val_idx])
        folds.append((train_idx, val_idx))

    return folds

def cross_validation(X, y, cv=5, max_size=None):

    # in n_sample < n_features, use dual formulation and stronger regularization
    if X.shape[0] < X.shape[1]:
        use_dual = True
        logging.debug(f"Using dual formulation for Logistic Regression as n_samples < n_features ({X.shape[0]} < {X.shape[1]})")
        C = 0.1
        logging.debug(f"Using stronger regularization (C={C})")
    else:
        use_dual = False
        logging.debug(f"Using primal formulation for Logistic Regression as n_samples >= n_features ({X.shape[0]} >= {X.shape[1]})")
        C = 1.0
        logging.debug(f"Using weaker regularization (C={C})")

    metrics = {metric: [] for metric in METRICS_TO_COMPUTE}
    for train_idx, val_idx in balanced_kfold_splits(y, cv=cv, max_size=max_size):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

        model = train_model(X_train, y_train, use_dual=use_dual, C=C)
        eval_metrics = eval_model(model, X_val, y_val)
        
        for metric_name, score in eval_metrics.items():
            metrics[metric_name].append(score)

    for metric_name in metrics:
        metrics[metric_name] = np.array(metrics[metric_name])

    return metrics
