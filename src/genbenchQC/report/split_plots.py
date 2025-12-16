import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from genbenchQC.utils.similarity_threshold import compute_threshold

# def plot_distribution_max_corr_alignment_scores(stratified_test_split, threshold):
#     """
#     Plot the distribution of maximum corrected Smith–Waterman alignment scores
#     from a hashFrag stratified_test_split dataframe.
#     """

#     # Extract score column
#     scores = stratified_test_split["score"]

#     # Define bin edges with bin size = 1
#     min_score = int(scores.min())
#     max_score = int(scores.max())
#     bin_edges = range(min_score, max_score + 2, 1)

#     # Compute histogram (counts and bin edges)
#     counts, bin_edges = np.histogram(scores, bins=bin_edges)

#     # Convert counts to percentage
#     total = counts.sum()
#     percentages = (counts / total) * 100

#     # Use bin centers for plotting
#     bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

#     # Create plot
#     fig, ax = plt.subplots(figsize=(8, 5))

#     ax.plot(bin_centers, percentages, marker="o", linestyle="-", color="#f28e2b")

#     ax.axvline(
#         threshold,
#         color="red",
#         linestyle="--",
#         linewidth=1.2,
#         label=f"Threshold = {threshold}"
#     )

#     ax.set_xlabel("Maximum Pairwise SW Alignment Score")
#     ax.set_ylabel("Frequency (%)")
#     ax.set_title("Test set Maximum SW Alignment Scores Distribution")
#     ax.legend()
    
#     ax.grid(alpha=0.3)
    
#     fig.tight_layout()

#     return fig

def plot_distribution_max_corr_alignment_scores(
    stratified_test_split, threshold
):
    """
    Plot the cumulative distribution (CDF) of maximum corrected
    Smith–Waterman alignment scores.
    """

    # Extract score column
    scores = stratified_test_split["score"]

    # Sort scores
    sorted_scores = np.sort(scores)

    # Empirical CDF (in percentage)
    cumulative_percent = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores) * 100

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(
        sorted_scores,
        cumulative_percent,
        linewidth=1.2,
        label="Cumulative distribution"
    )

    ax.axvline(
        threshold,
        color="red",
        linestyle="--",
        linewidth=1.2,
        label=f"Threshold = {threshold}"
    )

    ax.set_xlabel("Maximum Pairwise Alignment Score")
    ax.set_ylabel("Cumulative frequency (%)")
    ax.set_title("Cumulative Distribution of Test Set Alignment Scores")

    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()

    return fig

def plot_genomic_vs_shuffled_hist(df_genomic, df_shuffled, bin_width=5):
    genomic = df_genomic["score"]
    shuffled = df_shuffled["score"]
    threshold = compute_threshold(df_shuffled)

    all_scores = pd.concat([genomic, shuffled])
    bins = np.arange(all_scores.min(), all_scores.max() + bin_width, bin_width)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.hist(
        genomic,
        bins=bins,
        color="cornflowerblue",
        alpha=0.8,
        edgecolor="cornflowerblue",
        linewidth=0.5,
        label="Genomic"
    )
    ax.hist(
        shuffled,
        bins= bins,
        color="red",           
        alpha=0.25,           
        edgecolor="red",
        linewidth=0.8,
        label="Dinucleotide shuffled"
    )
    ax.axvline(
        threshold,
        color="red",
        linestyle="--",
        linewidth=1.2,
        label=f"Threshold = {threshold}"
    )

    ax.set_yscale("log")
    ax.set_xlabel("Maximum pairwise SW alignment score")
    ax.set_ylabel("Count (log scale)")
    ax.legend()

    fig.tight_layout()
    return fig


