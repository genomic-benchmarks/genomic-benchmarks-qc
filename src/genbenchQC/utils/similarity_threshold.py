import pandas as pd
import numpy as np

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

def get_threshold_stats(results, results_filt, coverage_threshold):
    threshold_stats = {
        "coverage_threshold": coverage_threshold,
        "perc_above_threshold": len(results_filt) / len(results) * 100 if len(results) > 0 else 0,
        "perc_below_threshold": (len(results) - len(results_filt)) / len(results) * 100 if len(results) > 0 else 0,
    }
    return threshold_stats