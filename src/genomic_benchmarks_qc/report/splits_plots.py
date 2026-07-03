import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns

from genomic_benchmarks_qc.report.classes_plots import HuePalette, prepare_legend
from genomic_benchmarks_qc.report.utils import FAIL_COLOR

def plot_similarity_histograms(query_similarity_max, target_similarity_max, threshold_stats):

    unique_queries = pd.Series(query_similarity_max, dtype=float)
    unique_targets = pd.Series(target_similarity_max, dtype=float)

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

    # Combine data into single DataFrame for hue-based plotting
    hist_data = []
    if tcov.size:
        for val, w in zip(tcov, tweights):
            hist_data.append({"similarity": val, "type": "Train", "weight": w})
    if qcov.size:
        for val, w in zip(qcov, qweights):
            hist_data.append({"similarity": val, "type": "Test", "weight": w})

    if hist_data:
        df = pd.DataFrame(hist_data)
        # bars next to each other instead of stacked
        sns.histplot(
            data=df,
            x="similarity",
            weights="weight",
            bins=bins,
            stat="count",
            hue="type",
            hue_order=["Train", "Test"],
            multiple="dodge",
            element="bars",
            shrink=0.8,
            ax=ax,
            palette=[palette[0], palette[1]],
            legend=False,
            linewidth=0,
        )

    # Build legend handles for hue categories (Train, Test)
    legend_handles = [
        Patch(facecolor=palette[0], label="Train"),
        Patch(facecolor=palette[1], label="Test"),
    ]
    legend_labels = ["Train", "Test"]

    ax.axvline(
        threshold_stats["similarity_threshold"],
        linestyle="--",
        linewidth=1.2,
        color=FAIL_COLOR,
        label="Threshold",
    )
    legend_handles.append(
        Line2D([0], [0], linestyle="--", linewidth=1.2, color=FAIL_COLOR, label="Threshold")
    )
    legend_labels.append("Threshold")

    ax.set_xlim(0, 100)
    # Set x-tick labels as half-open bin intervals; last bin is closed on both ends
    ax.set_xticks(np.arange(5, 101, 10))
    ax.set_xticklabels(
        [f"[{i}, {i+10})" for i in range(0, 90, 10)] + ["[90, 100]"]
    )
    ax.set_yscale("log")
    ax.set_xlabel("Sequence similarity (%)", fontsize=14)
    ax.set_ylabel("Count (log scale)", fontsize=14)
    ax.tick_params(axis='both', labelsize=12)

    ax = prepare_legend(
        ax,
        box_to_anchor=(0.5, -0.2),
        legend_handles=legend_handles,
        legend_labels=legend_labels
    )
    
    return fig
