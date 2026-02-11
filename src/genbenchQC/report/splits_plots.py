import matplotlib.pyplot as plt
import numpy as np

def plot_coverage_histograms(results, threshold_stats):

    # Get row indices of max coverage per query/target
    query_idx = results.groupby('query')['min_cov'].idxmax()
    target_idx = results.groupby('target')['min_cov'].idxmax()

    # Select the actual rows
    unique_queries = results.loc[query_idx, 'qcov']
    unique_targets = results.loc[target_idx, 'tcov']

    # Add zeros for missing combinations (i.e. those with no hits, which have 0 coverage)
    missing_data = threshold_stats['total_combinations'] - threshold_stats['hits']

    # Create arrays for plotting, filling missing combinations with zeros
    qcov = np.concatenate([unique_queries.to_numpy(), np.zeros(missing_data)])
    tcov = np.concatenate([unique_targets.to_numpy(), np.zeros(missing_data)])

    # Define bins for histogram
    bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    
    # Create histogram plot
    fig, ax = plt.subplots()
    ax.hist(qcov, bins=bins, alpha=0.5, label="Query (Test)")
    ax.hist(tcov, bins=bins, alpha=0.5, label="Target (Train)")
    ax.axvline(threshold_stats["coverage_threshold"], linestyle="--", linewidth=1.2, color="red", label="Threshold")
    ax.set_yscale("log")
    ax.set_xlabel("Max. Coverage Per Unique Sequence") # Max of min(qcov, tcov) per unique query/target
    ax.set_ylabel("Count (log scale)")
    ax.legend()
    return fig