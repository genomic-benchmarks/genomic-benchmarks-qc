import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import pandas as pd
import seaborn as sns

from genbenchQC.report.classes_plots import HuePalette, prepare_legend

def plot_similarity_histograms(results, threshold_stats):

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

    bins = np.arange(0, 110, 10)

    fig, ax = plt.subplots(figsize=(12, 4), dpi=300)
    sns.set_style("white")
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
