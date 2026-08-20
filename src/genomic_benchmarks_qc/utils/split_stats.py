"""Turning the raw search summary into the numbers the split report shows."""


def _calculate_percentage(part, whole):
    """Return part as a percentage of whole, treating an empty whole as 0%."""
    return (part / whole) * 100 if whole > 0 else 0.0


def flag_split_data_leakage(perc_queries_above_thr, fail_threshold=2.0):
    """Grade a split by the share of test sequences with a near-identical match.

    Any leakage at all is worth knowing about, so a non-zero share is a Warning
    and only `fail_threshold` percent or more is a Fail.
    """
    if perc_queries_above_thr >= fail_threshold:
        return "Fail"
    if perc_queries_above_thr > 0:
        return "Warning"
    return "Pass"


def get_basic_stats_from_aggregates(filename_train, train_stats, filename_test, test_stats):
    """Combine the two halves' sequence counts and lengths for the report header."""
    return {
        "train_filename": str(filename_train),
        "test_filename": str(filename_test),
        "number_of_sequences_train": int(train_stats["count"]),
        "number_of_sequences_test": int(test_stats["count"]),
        "min_length_train": int(train_stats["min_length"]),
        "mean_length_train": float(train_stats["mean_length"]),
        "max_length_train": int(train_stats["max_length"]),
        "min_length_test": int(test_stats["min_length"]),
        "mean_length_test": float(test_stats["mean_length"]),
        "max_length_test": int(test_stats["max_length"]),
    }


def get_threshold_stats(summary, similarity_threshold, num_train_seqs, num_test_seqs):
    """Turn the search summary into the leakage counts and percentages.

    Queries are test sequences and targets are train sequences, so the two
    percentages answer different questions: how much of the test set is
    compromised, and how much of the training set is responsible.
    """
    if summary is None:
        raise ValueError("summary is required for threshold statistics computation")

    num_queries_with_hits = len(summary["query_ids_with_hits"])
    num_queries_above_thr = len(summary["query_ids_above_threshold"])
    num_targets_with_hits = len(summary["target_ids_with_hits"])
    num_targets_above_thr = len(summary["target_ids_above_threshold"])

    return {
        "similarity_threshold": similarity_threshold,
        "perc_queries_above_thr": _calculate_percentage(num_queries_above_thr, num_test_seqs),
        "num_queries_above_thr": num_queries_above_thr,
        "num_all_queries": num_test_seqs,
        "num_queries_without_hits": max(num_test_seqs - num_queries_with_hits, 0),
        "perc_targets_above_thr": _calculate_percentage(num_targets_above_thr, num_train_seqs),
        "num_targets_above_thr": num_targets_above_thr,
        "num_all_targets": num_train_seqs,
        "num_targets_without_hits": max(num_train_seqs - num_targets_with_hits, 0),
        "hits": summary["total_hits"],
        "total_combinations": num_train_seqs * num_test_seqs,
    }
