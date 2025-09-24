import logging
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, accuracy_score
import numpy as np
import pandas as pd
from typing import List, Tuple


def extract_per_position_base(sequences, base, reverse):
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

def flag_on_score(score):
    if score > 0.8:
        return "Fail"
    elif score > 0.7:
        return "Major Warning"
    elif score > 0.6:
        return "Warning"
    else:
        return "Pass"

def model(stats1, stats2, max_class_size=None):

    results = {}

    for stat in STATS_TO_TRAIN_PRECOMPUTED:

        logging.info(f"Training bias detection model for statistic: {stat}")

        X = pd.concat([stats1.stats[stat], stats2.stats[stat]], axis=0)
        X.fillna(0, inplace=True)
        y = pd.Series([1] * len(stats1.stats[stat]) + [0] * len(stats2.stats[stat]))

        avg_precisions, accuracies = cross_validation(X, y, cv=5, max_size=max_class_size)
        avg_score = avg_precisions.mean()
        acc_score = accuracies.mean()

        logging.debug(f"Accuracy scores for {stat}: {accuracies}")
        results = add_result(results, stat, avg_score, acc_score)

    common_nts = list(set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases']))

    for nt in common_nts:
        logging.info(f"Training bias detection model for per position nucleotide: {nt}")
        for reverse in [False, True]:

            features = extract_per_position_base(stats1.sequences + stats2.sequences, base=nt, reverse=reverse)

            X = pd.DataFrame(features)
            y = pd.Series([1] * len(stats1.sequences) + [0] * len(stats2.sequences))

            avg_precisions, accuracies = cross_validation(X, y, cv=5, max_size=max_class_size)
            avg_score = avg_precisions.mean()
            acc_score = accuracies.mean()

            flag_name = f"Per position nucleotide content - {nt}" if not reverse else f"Per reverse position nucleotide content - {nt}"
            logging.debug(f"Accuracy scores for {flag_name}: {accuracies}")
            results = add_result(results, flag_name, avg_score, acc_score)

    return results

def add_result(results, key, avg_score, acc_score):
    results[key] = {}
    results[key]['Average Precision'] = avg_score
    results[key]['Accuracy'] = acc_score
    results[key]['Flag'] = flag_on_score(acc_score)

    return results

def train_model(X, y):
    model = LogisticRegression(random_state=42, solver='sag')
    model.fit(X, y)

    return model

def eval_model(model, X, y):
    y_prob = model.predict_proba(X)[:, 1]
    avg_precision = average_precision_score(y, y_prob)
    y_pred = y_prob >= 0.5
    accuracy = accuracy_score(y, y_pred)

    return avg_precision, accuracy

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

    avg_precisions = []
    accuracies = []
    for train_idx, val_idx in balanced_kfold_splits(y, cv=cv, max_size=max_size):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

        model = train_model(X_train, y_train)
        avg_precision, accuracy = eval_model(model, X_val, y_val)
        avg_precisions.append(avg_precision)
        accuracies.append(accuracy)

    return np.array(avg_precisions), np.array(accuracies)
