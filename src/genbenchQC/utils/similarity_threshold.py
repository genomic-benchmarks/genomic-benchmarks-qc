import random
import pandas as pd
import numpy as np

def dinucleotide_shuffle(seq):
    dinucleotides = [seq[i:i+2] for i in range(len(seq) - 1)]
    random.shuffle(dinucleotides)
    
    shuffled_seq = dinucleotides[0]
    for dinucleotide in dinucleotides[1:]:
        if shuffled_seq[-1] == dinucleotide[0]:
            shuffled_seq += dinucleotide[1]
        else:
            candidates = [d for d in dinucleotides if d[0] == shuffled_seq[-1] and d != dinucleotide]
            if candidates:
                choice = random.choice(candidates)
                shuffled_seq += choice[1]
                dinucleotides.remove(choice)
            else:
                shuffled_seq += dinucleotide
    
    return shuffled_seq

def dinucleotide_shuffle_list(seq_list, seed=None):
    
    if seed is not None:
        random.seed(seed)
    
    shuffled = []
    for seq in seq_list:
        shuffled.append(dinucleotide_shuffle(seq))
    
    return shuffled

def compute_threshold(df_dinucleotide_shuffled: pd.DataFrame,
                      add_margin: float = 0.0) -> int:

    max_score = df_dinucleotide_shuffled["score"].max()
    threshold = int(max_score + add_margin)
    return threshold

# def compute_threshold(df_dinucleotide_shuffled: pd.DataFrame) -> int:
#     threshold = df_dinucleotide_shuffled["score"].quantile(0.95)
#     return int(threshold)

# def compute_threshold(df_dinucleotide_shuffled, method="percentile", percentile=0.95, B=1000):
#     """
#     Bootstrap-based threshold using shuffled alignment scores.

#     Parameters
#     ----------
#     df_dinucleotide_shuffled : pandas.DataFrame
#         Must contain a numeric column 'score'.
#     method : {"percentile", "max"}
#         - "percentile": bootstrap a percentile (default 95th)
#         - "max": bootstrap the maximum
#     percentile : float
#         Percentile used if method="percentile" (e.g. 0.95).
#     B : int
#         Number of bootstrap resamples.

#     Returns
#     -------
#     float
#         Bootstrap-estimated threshold.
#     """
#     scores = df_dinucleotide_shuffled["score"].to_numpy()
#     n = len(scores)

#     boot_vals = []
#     for _ in range(B):
#         sample = np.random.choice(scores, size=n, replace=True)

#         if method == "percentile":
#             val = np.percentile(sample, percentile * 100)
#         elif method == "max":
#             val = sample.max()
#         else:
#             raise ValueError("method must be 'percentile' or 'max'")

#         boot_vals.append(val)

#     return np.mean(boot_vals)

def get_basic_stats(filename_train, train_sequences, filename_test, test_sequences):
    basic_stats = {
        "train_filename": str(filename_train),
        "test_filename": str(filename_test),
        "number_of_sequences_train": len(train_sequences),
        "number_of_sequences_test": len(test_sequences),
        "min_length_train": min(len(seq) for seq in train_sequences),
        "mean_length_train": sum(len(seq) for seq in train_sequences) // len(train_sequences),
        "max_length_train": max(len(seq) for seq in train_sequences),
        "min_length_test": min(len(seq) for seq in test_sequences),
        "mean_length_test": sum(len(seq) for seq in test_sequences) // len(test_sequences),
        "max_length_test": max(len(seq) for seq in test_sequences),
    }
    return basic_stats

def get_threshold_stats(stratified_test_split, threshold):
    threshold_stats = {
        "threshold": threshold,
        "num_below_threshold": sum(1 for score in stratified_test_split['score'] if score <= threshold),
        "num_above_threshold": sum(1 for score in stratified_test_split['score'] if score > threshold),
        "total_alignments": len(stratified_test_split),
    }
    return threshold_stats