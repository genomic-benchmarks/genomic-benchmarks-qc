import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_similarity_histograms(results, threshold_stats):

    # Get the maximum signal (min_cov*pident) for each unique query and target
    def _max_signal_per_id(df: pd.DataFrame, id_col: str, metric_col: str):
        if df.empty:
            return pd.Series(dtype=float)
        idx = df.groupby(id_col)[metric_col].idxmax()
        return df.loc[idx, f"{id_col[0]}cov*pident"]
    
    results = results.copy()
    if not results.empty:
        results['qcov*pident'] = results['qcov'] * results['pident']
        results['tcov*pident'] = results['tcov'] * results['pident']
        unique_queries = _max_signal_per_id(results, 'query', 'min_cov*pident')
        unique_targets = _max_signal_per_id(results, 'target', 'min_cov*pident')
    else:
        unique_queries = pd.Series(dtype=float)
        unique_targets = pd.Series(dtype=float)

    # Count combinations without any hits
    missing_data = max(threshold_stats['total_combinations'] - threshold_stats['hits'], 0)

    def _build_hist_arrays(values: pd.Series):
        data = values.to_numpy(dtype=float)
        weights = np.ones_like(data, dtype=float)
        if missing_data > 0:
            data = np.concatenate([data, np.array([0.0], dtype=float)])
            weights = np.concatenate([weights, np.array([missing_data], dtype=float)])
        return data, weights

    qcov, qweights = _build_hist_arrays(unique_queries)
    tcov, tweights = _build_hist_arrays(unique_targets)

    # Define bins for histogram
    bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    
    # Create histogram plot
    fig, ax = plt.subplots()
    if qcov.size:
        ax.hist(qcov, bins=bins, weights=qweights, alpha=0.5, label="Test")
    else:
        ax.hist([], bins=bins, alpha=0.5, label="Test")
    if tcov.size:
        ax.hist(tcov, bins=bins, weights=tweights, alpha=0.5, label="Train")
    else:
        ax.hist([], bins=bins, alpha=0.5, label="Train")
    ax.axvline(threshold_stats["similarity_threshold"], linestyle="--", linewidth=1.2, color="red", label="Threshold")
    ax.set_yscale("log")
    ax.set_xlabel("Sequence Similarity") # Max of min(qcov, tcov) per unique query/target
    ax.set_ylabel("Count (log scale)")
    ax.legend()
    return fig
