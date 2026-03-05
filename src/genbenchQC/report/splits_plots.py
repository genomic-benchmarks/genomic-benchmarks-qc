import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import pandas as pd
import seaborn as sns

from genbenchQC.report.classes_plots import HuePalette, prepare_legend

def plot_similarity_histograms(results, threshold_stats):

    # Get unique queries and targets with their max similarity
    unique_queries = results.groupby('query')['min_cov*pident'].max()
    unique_targets = results.groupby('target')['min_cov*pident'].max()

    # Get counts of missing data (queries/targets without hits) to add to histogram
    missing_data_query = threshold_stats["num_queries_without_hits"]
    missing_data_target = threshold_stats["num_targets_without_hits"]

    def _build_hist_arrays(values: pd.Series, missing_data: int):
        data = values.to_numpy(dtype=float)
        weights = np.ones_like(data, dtype=float)
        if missing_data > 0:
            data = np.concatenate([data, np.array([0.0], dtype=float)])
            weights = np.concatenate([weights, np.array([missing_data], dtype=float)])
        return data, weights

    qcov, qweights = _build_hist_arrays(unique_queries, missing_data_query)
    tcov, tweights = _build_hist_arrays(unique_targets, missing_data_target)

    bins = np.arange(0, 110, 10)

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(12, 4), dpi=300)
    palette = HuePalette()
    hist_kwargs = dict(alpha=0.65, edgecolor="black")

    if qcov.size:
        ax.hist(qcov, bins=bins, weights=qweights, label="Test", color=palette[0], **hist_kwargs)
    if tcov.size:
        ax.hist(tcov, bins=bins, weights=tweights, label="Train", color=palette[1], **hist_kwargs)

    ax.axvline(
        threshold_stats["similarity_threshold"],
        linestyle="--",
        linewidth=1.2,
        color="red",
        label="Threshold",
    )
    ax.set_xlim(0, 100)
    ax.set_xticks(np.arange(0, 101, 10))
    ax.set_yscale("log")
    ax.set_xlabel("Sequence similarity (%)", fontsize=14)
    ax.set_ylabel("Count (log scale)", fontsize=14)
    ax.tick_params(axis='both', labelsize=12)
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10))
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation(base=10))
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.0)
        spine.set_edgecolor("black")

    prepare_legend(ax, box_to_anchor=(0.5, -0.2))
    fig.tight_layout()
    return fig
