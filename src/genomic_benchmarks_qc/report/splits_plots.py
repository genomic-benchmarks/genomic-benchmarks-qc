"""The figures for the split report."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from genomic_benchmarks_qc.report.classes_plots import HuePalette, prepare_legend
from genomic_benchmarks_qc.report.utils import FAIL_COLOR

# The similarity axis: ten bins of ten percent, the last closed at both ends so
# a perfect 100% match lands in it rather than off the end.
BIN_WIDTH = 10
BIN_EDGES = np.arange(0, 110, BIN_WIDTH)


def _binned_counts(similarity_max, without_hits):
    """Count one half's sequences into the similarity bins.

    Sequences with no hit at all are absent from the search output and from
    `similarity_max`, but leaving them out would hide how much of the half is
    unrelated to the other one, so they are counted into the first bin.

    Args:
        similarity_max: Best similarity per sequence, NaN where there was no hit.
        without_hits: How many sequences the search returned nothing for.

    Returns:
        Array of `len(BIN_EDGES) - 1` counts.
    """
    values = np.asarray(similarity_max, dtype=float)
    values = values[~np.isnan(values)]
    counts, _ = np.histogram(values, bins=BIN_EDGES)
    counts = counts.astype(float)
    counts[0] += without_hits
    return counts


def plot_similarity_histograms(query_similarity_max, target_similarity_max, threshold_stats):
    """Plot how similar each sequence's best match in the other half is.

    One bar pair per similarity bin, with the leakage threshold marked. The
    count axis is logarithmic because leakage is usually a small tail next to a
    large bulk of unrelated sequences, which a linear axis would flatten away.

    Ten bars are counted rather than drawn from the observations. Handing the
    per-sequence maxima to a plotting library meant one Python object per
    sequence on the way in - 257 MB and a second at 300,000 sequences a half,
    for a figure that is ten numbers wide.
    """
    counts_train = _binned_counts(target_similarity_max,
                                  threshold_stats["num_targets_without_hits"])
    counts_test = _binned_counts(query_similarity_max,
                                 threshold_stats["num_queries_without_hits"])

    fig, ax = plt.subplots(figsize=(12, 4), dpi=300)
    palette = HuePalette()

    # Train and test share each bin rather than stacking, so each takes half of
    # it, and each bar is drawn at 0.8 of the half it has.
    centres = BIN_EDGES[:-1] + BIN_WIDTH / 4
    bar_width = BIN_WIDTH / 2 * 0.8
    ax.bar(centres, counts_train, width=bar_width, color=palette[0], linewidth=0)
    ax.bar(centres + BIN_WIDTH / 2, counts_test, width=bar_width,
           color=palette[1], linewidth=0)

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
