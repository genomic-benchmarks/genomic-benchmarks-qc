import logging
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score, precision_recall_curve, auc
import numpy as np
import pandas as pd
from typing import List, Tuple


def extract_per_position_base(sequences, base, reverse):
    """
    Build a per-position binary feature matrix indicating where a given nucleotide/base occurs in each sequence.
    
    Parameters:
        sequences (Iterable[str]): Collection of sequences to encode. Sequences may have varying lengths.
        base (str): Single-character nucleotide/base to mark in the output features.
        reverse (bool): If True, treat each sequence in reverse-complement orientation by reversing its order before encoding.
    
    Returns:
        numpy.ndarray: 2D array of shape (n_sequences, max_sequence_length) where each element is 1 if the sequence has `base` at that position (after optional reversal) and 0 otherwise. Positions beyond a sequence's length are 0.
    """
    max_len = max([len(seq) for seq in sequences])
    features = np.zeros((len(sequences), max_len))
    for i, seq in enumerate(sequences):
        if reverse:
            seq = seq[::-1]
        for j, nt in enumerate(seq):
            features[i, j] = 1 if nt == base else 0
    return features


STATS_TO_TRAIN_PRECOMPUTED = [
    'Per sequence nucleotide content',
    'Per sequence dinucleotide content',
    'Per sequence GC content',
    'Sequence lengths',
]

METRICS_TO_COMPUTE = ['AU-ROC', 'AU-PR', 'Accuracy']

def flag_on_score(score):
    """
    Map a numeric metric score to a categorical flag indicating bias level.
    
    Parameters:
        score (float): Metric value (expected between 0 and 1) used to determine the flag.
    
    Returns:
        str: `"Fail"` if score > 0.7, `"Warning"` if score > 0.6, `"Pass"` otherwise.
    """
    if score > 0.7:
        return "Fail"
    elif score > 0.6:
        return "Warning"
    else:
        return "Pass"

def model(stats1, stats2, max_class_size=None, metric_to_flag='AU-ROC'):

    """
    Train bias-detection classifiers comparing two datasets and return aggregated cross-validated metrics and flags.
    
    This function trains logistic regression models to detect dataset-specific signal for each statistic named in STATS_TO_TRAIN_PRECOMPUTED and for per-position nucleotide presence for bases common to both inputs (both forward and reversed). For each model it computes cross-validated metrics via cross_validation, aggregates mean scores and a flag using add_result, and returns a dictionary of results keyed by statistic or feature name.
    
    Parameters:
        stats1: An object exposing `.stats` (a mapping of statistic name -> Series/DataFrame) and `.sequences` (list of sequences) for the first dataset.
        stats2: An object exposing `.stats` and `.sequences` for the second dataset.
        max_class_size: Optional maximum number of samples per class to use when balancing folds; if None no additional cap is applied.
        metric_to_flag: Name of the metric from METRICS_TO_COMPUTE whose mean value is used to derive the 'Flag' entry for each result.
    
    Raises:
        ValueError: If `metric_to_flag` is not one of METRICS_TO_COMPUTE.
    
    Returns:
        dict: Mapping from statistic/feature name to a dict containing mean metric values and a 'Flag' derived from the mean of `metric_to_flag`.
    """
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

    for nt in common_nts:
        logging.info(f"Training bias detection model for per position nucleotide: {nt}")
        for reverse in [False, True]:

            features = extract_per_position_base(stats1.sequences + stats2.sequences, base=nt, reverse=reverse)

            X = pd.DataFrame(features)
            y = pd.Series([1] * len(stats1.sequences) + [0] * len(stats2.sequences))

            metrics = cross_validation(X, y, cv=5, max_size=max_class_size)

            flag_name = f"Per position nucleotide content - {nt}" if not reverse else f"Per reverse position nucleotide content - {nt}"
            logging.debug(f"{metric_to_flag} scores for {flag_name}: {metrics[metric_to_flag]}")
            results = add_result(results, flag_name, metrics, metric_to_flag)

    return results

def add_result(results, key, metrics, metric_to_flag):
    """
    Aggregate per-fold metric arrays into mean scores under the given key and attach a quality flag.
    
    Parameters:
        results (dict): Mutable mapping to update with aggregated metrics for `key`.
        key (hashable): Dictionary key under which the aggregated metrics and flag will be stored.
        metrics (Mapping[str, array-like]): Mapping from metric names to numeric sequences (e.g., NumPy arrays)
            containing per-fold scores; each value must support computing the mean.
        metric_to_flag (str): Name of the metric in `metrics` whose mean will be used to derive the 'Flag'.
    
    Returns:
        dict: The same `results` mapping after inserting or replacing results[key] with the aggregated
        metric means and a 'Flag' computed from the mean of `metrics[metric_to_flag]`.
    """
    results[key] = {}
    for metric_name, scores in metrics.items():
        avg_score = scores.mean()
        results[key][metric_name] = avg_score
    results[key]['Flag'] = flag_on_score(metrics[metric_to_flag].mean())

    return results

def train_model(X, y, use_dual=False, C=1.0):

    """
    Train a logistic regression classifier on the provided feature matrix and labels.
    
    Parameters:
        X (array-like or pandas.DataFrame): Feature matrix with shape (n_samples, n_features).
        y (array-like or pandas.Series): Target labels for each sample.
        use_dual (bool): Whether to use the dual formulation of the solver.
        C (float): Inverse of regularization strength; smaller values specify stronger regularization.
    
    Returns:
        sklearn.linear_model.LogisticRegression: Fitted LogisticRegression model.
    """
    model = LogisticRegression(random_state=42, solver='liblinear', dual=use_dual, max_iter=200, C=C)
    model.fit(X, y)

    return model

def eval_model(model, X, y):
    """
    Compute discrimination and classification metrics for a fitted binary classifier on the provided dataset.
    
    Parameters:
        model: A fitted classifier implementing `predict_proba(X)` and returning probability estimates for the positive class.
        X: Feature matrix for evaluation (array-like or DataFrame).
        y: True binary labels corresponding to X (array-like, values interpreted as 0/1).
    
    Returns:
        dict: Mapping metric names to float scores:
            - 'AU-ROC': Area under the receiver operating characteristic curve.
            - 'AU-PR': Area under the precision–recall curve.
            - 'Accuracy': Fraction of correct predictions using a 0.5 probability threshold.
    """
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
    """
    Generate balanced k-fold train/validation index pairs by subsampling classes to equal sizes.
    
    Parameters:
        y (array-like): Array of class labels for each sample; used to construct balanced splits.
        cv (int): Number of folds to produce.
        max_size (int | None): If provided, cap the per-class sample count to this value before folding.
    
    Returns:
        List[Tuple[np.ndarray, np.ndarray]]: A list of length `cv` where each element is a tuple (train_idx, val_idx).
        Each `train_idx` and `val_idx` is a NumPy array of sample indices. Classes are balanced by
        subsampling to the minimum class size (or `max_size` if smaller), the indices are shuffled,
        and the final fold receives any remainder when samples do not divide evenly.
    """
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
    """
    Perform cross-validated evaluation of logistic regression models using balanced folds.
    
    Parameters:
        X (pd.DataFrame): Feature matrix where rows are samples and columns are features.
        y (pd.Series): Binary labels aligned with X's rows.
        cv (int): Number of cross-validation folds to produce.
        max_size (int, optional): Maximum number of samples per class to use when balancing folds; if None, use full class sizes.
    
    Returns:
        dict: Mapping metric names (e.g., 'AU-ROC', 'AU-PR', 'Accuracy') to NumPy arrays of scores across the CV folds.
    """
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