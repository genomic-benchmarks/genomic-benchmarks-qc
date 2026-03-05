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

    bins = list(range(0, 110, 10))

    fig, ax = plt.subplots(figsize=(12, 4), dpi=300)
    palette = HuePalette()

    if qcov.size:
        qhist_df = pd.DataFrame({"similarity": qcov, "weight": qweights})
        sns.histplot(
            data=qhist_df,
            x="similarity",
            weights="weight",
            bins=bins,
            stat="count",
            element="bars",
            color=palette[0],
            label="Test",
            ax=ax,
            alpha=0.7,
        )
    if tcov.size:
        thist_df = pd.DataFrame({"similarity": tcov, "weight": tweights})
        sns.histplot(
            data=thist_df,
            x="similarity",
            weights="weight",
            bins=bins,
            stat="count",
            element="bars",
            color=palette[1],
            label="Train",
            ax=ax,
            alpha=0.7,
        )

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

    ax = prepare_legend(ax, box_to_anchor=(0.5, -0.2))
    return fig
