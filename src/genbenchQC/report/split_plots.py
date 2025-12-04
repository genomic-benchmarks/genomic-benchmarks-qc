import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_distribution_max_corr_alignment_scores(stratified_test_split):
    """
    Plot the distribution of maximum corrected Smith–Waterman alignment scores
    from a hashFrag stratified_test_split dataframe.
    """

    # Extract score column
    scores = stratified_test_split["score"]

    # Compute histogram (counts and bin edges)
    counts, bin_edges = np.histogram(scores, bins=30)

    # Convert counts to percentage
    total = counts.sum()
    percentages = (counts / total) * 100

    # Use bin centers for plotting
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(bin_centers, percentages, marker="o", linestyle="-", color="#f28e2b")

    ax.set_xlabel("Maximum Pairwise SW Alignment Score")
    ax.set_ylabel("Frequency (%)")
    ax.set_title("Test set Maximum SW Alignment Scores Distribution")

    ax.grid(alpha=0.3)
    fig.tight_layout()

    return fig
