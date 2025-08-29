import logging
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
import numpy as np
import pandas as pd


def extract_per_position_base(sequence, reverse):
    features = {}
    if reverse:
        sequence = sequence[::-1]
    for i, base in enumerate(sequence):
        features[f'pos_{i}_{base}'] = 1
    return features

STATS_TO_TRAIN_PRECOMPUTED = [
    'Per sequence nucleotide content',
    'Per sequence dinucleotide content',
    'Per sequence GC content',
    'Sequence lengths',
    #'Per position nucleotides',
    #'Per position reversed nucleotides'
]

def flag_on_avg_precision(avg_precision):
    if avg_precision > 0.8:
        return "Fail"
    elif avg_precision > 0.7:
        return "Major Warning"
    elif avg_precision > 0.6:
        return "Warning"
    else:
        return "Pass"

def model(stats1, stats2):

    results = {}

    for stat in STATS_TO_TRAIN_PRECOMPUTED:

        logging.info(f"Training bias detection model for statistic: {stat}")

        X = pd.concat([stats1.stats[stat], stats2.stats[stat]], axis=0)
        X.fillna(0, inplace=True)
        y = pd.Series([1] * len(stats1.stats[stat]) + [0] * len(stats2.stats[stat]))

        scores = cross_validation(X, y, cv=5)
        avg_score = scores.mean()

        logging.debug(f"Scores for {stat}: {scores}")
        results[stat] = {}
        results[stat]['Average Precision'] = avg_score
        results[stat]['Flag'] = flag_on_avg_precision(avg_score)

    common_nts = list(set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases']))

    for nt in common_nts:
        logging.info(f"Training bias detection model for per position nucleotide: {nt}")

        max_len = int(max(max(stats1.stats['Sequence lengths'].values.flatten()), max(stats2.stats['Sequence lengths'].values.flatten())))

        # Per position nucleotide content
        features = np.zeros((len(stats1.sequences) + len(stats2.sequences), max_len))
        for i, seq in enumerate(stats1.sequences + stats2.sequences):
            for j, base in enumerate(seq):
                if base == nt:
                    features[i, j] = 1

        X = pd.DataFrame(features)
        y = pd.Series([1] * len(stats1.sequences) + [0] * len(stats2.sequences))

        scores = cross_validation(X, y, cv=5)
        avg_score = scores.mean()

        logging.debug(f"Scores for per position nucleotide {nt}: {scores}")
        results[f'Per position nucleotide content - {nt}'] = {}
        results[f'Per position nucleotide content - {nt}']['Average Precision'] = avg_score
        results[f'Per position nucleotide content - {nt}']['Flag'] = flag_on_avg_precision(avg_score)

        # Per reverse position nucleotide content
        features = np.zeros((len(stats1.sequences) + len(stats2.sequences), max_len))
        for i, seq in enumerate(stats1.sequences + stats2.sequences):
            for j, base in enumerate(seq[::-1]):
                if base == nt:
                    features[i, j] = 1

        X = pd.DataFrame(features)
        y = pd.Series([1] * len(stats1.sequences) + [0] * len(stats2.sequences))

        scores = cross_validation(X, y, cv=5)
        avg_score = scores.mean()

        logging.debug(f"Scores for per reverse position nucleotide {nt}: {scores}")
        results[f'Per reverse position nucleotide content - {nt}'] = {}
        results[f'Per reverse position nucleotide content - {nt}']['Average Precision'] = avg_score
        results[f'Per reverse position nucleotide content - {nt}']['Flag'] = flag_on_avg_precision(avg_score)

    return results

def train_model(X, y):
    model = LogisticRegression(random_state=42, solver='sag')
    model.fit(X, y)

    return model

def eval_model(model, X, y):
    y_pred = model.predict(X)
    avg_precision = average_precision_score(y, y_pred)

    return avg_precision

def cross_validation(X, y, cv=5):
    # Create custom balanced k-fold splits
    unique_labels = np.unique(y)
    min_class_size = min([sum(y == label) for label in unique_labels])
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
    scores = []
    for i in range(cv):
        start_idx = i * fold_size
        end_idx = start_idx + fold_size if i < cv - 1 else len(balanced_indices)
        val_idx = balanced_indices[start_idx:end_idx]
        train_idx = np.array([idx for idx in balanced_indices if idx not in val_idx])
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]

        model = train_model(X_train, y_train)
        avg_precision = eval_model(model, X_val, y_val)
        scores.append(avg_precision)

    return np.array(scores)
