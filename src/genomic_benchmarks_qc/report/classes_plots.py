"""The figures comparing two classes.

Each `plot_*` builds and returns a matplotlib figure; saving them is the
report generator's job. Features that were flagged are marked in the figure
itself - underlines under failing positions, colored bands over them - so a
plot carries the same verdict as the table it sits next to in the report.

Marking is a step of its own rather than an argument to the drawing, because
the report keeps both versions of every flagged figure. The generator draws
once, saves the clean file, calls the matching `mark_*` and saves again;
passing the flags in meant building the whole figure a second time to get the
second file, on exactly the datasets that have one.
"""

import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch

from genomic_benchmarks_qc.report.colors import CLASS_COLORS
from genomic_benchmarks_qc.report.utils import FAIL_COLOR, WARN_COLOR
from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics

logger = logging.getLogger(__name__)


def plot_lengths(stats1, stats2, plot_type='boxen'):
    """
    Plot the sequence lengths of two sequences.
    """
    return plot_one_stat(
        stats1, stats2,
        'Sequence lengths',
        plot_type=plot_type
    )

def plot_gc_content(stats1, stats2, plot_type='boxen'):
    """
    Plot the GC content of two sequences.
    """
    return plot_one_stat(
        stats1, stats2,
        'Per sequence GC content',
        plot_type=plot_type
    )

def reserve_flag_margin(ax: Axes, n_positions: int) -> None:
    """Reserve the bottom margin used by flag underlines so geometry is stable.

    The colored underlines from add_failed_outline are drawn with clip_on=False,
    which enlarges the tight bounding box used by savefig(bbox_inches='tight').
    Without this reservation, flagged and unflagged plots crop to different
    bounding boxes -> boxen boxes render at different widths/spacing. Drawing
    invisible stubs at every position makes the tight bbox identical whether or
    not any position is actually flagged.

    Args:
        ax: Matplotlib axis object.
        n_positions: Number of x categories to reserve space for.
    """
    xaxis_transform = ax.get_xaxis_transform()
    xticks = list(ax.get_xticks())
    for pos in range(min(n_positions, len(xticks))):
        x_center = xticks[pos]
        ax.plot(
            [x_center - 0.35, x_center + 0.35],
            [0, 0],
            color='none',
            linewidth=6,
            alpha=0,
            transform=xaxis_transform,
            clip_on=False,
            solid_capstyle='round',
            zorder=10,
        )

def add_failed_outline(ax: Axes, failed_positions: dict[int, str],
                       x_left_custom: float | None = None,
                       x_right_custom: float | None = None) -> None:
    """Add colored underlines to failed positions.

    Args:
        ax: Matplotlib axis object.
        failed_positions: Dict mapping position index -> flag status ('Fail' or
            'Warning').
        x_left_custom: Optional custom left x-coordinate for the underline.
        x_right_custom: Optional custom right x-coordinate for the underline.
    """

    xaxis_transform = ax.get_xaxis_transform()
    xticks = list(ax.get_xticks())

    for pos, flag in failed_positions.items():
        edgecolor = FAIL_COLOR if flag == 'Fail' else WARN_COLOR
        if pos < len(xticks):
            x_center = xticks[pos]
            if x_left_custom is None or x_right_custom is None:
                x_left = x_center - 0.35
                x_right = x_center + 0.35
            else:
                x_left = x_left_custom
                x_right = x_right_custom
            ax.plot(
                [x_left, x_right],
                [0, 0],
                color=edgecolor,
                linewidth=6,
                alpha=0.8,
                transform=xaxis_transform,
                clip_on=False,
                solid_capstyle='round',
                zorder=10,
            )

def plot_nucleotides(stats1: SequenceStatistics, stats2: SequenceStatistics,
                     nucleotides: list[str], plot_type: str) -> Figure:
    """
    Plot the nucleotide content of two sets of sequences.

    Args:
        stats1: Statistics for the first set of sequences.
        stats2: Statistics for the second set of sequences.
        nucleotides: List of nucleotides to plot.
        plot_type: Type of plot to create (e.g., 'boxen', 'violin').

    Returns:
        Matplotlib figure object. `mark_failed_nucleotides` adds the flags to it.
    """

    df = melt_stats(stats1, stats2, 'Per sequence nucleotide content',
                    var_name='Nucleotide', value_name='Frequency')

    fig, ax = plt.subplots(1, 1, figsize=(12, 4), dpi=300)
    if plot_type == 'violin':
        sns.violinplot(
            x='Nucleotide',
            y='Frequency',
            hue="label",
            split=True,
            data=df[df['Nucleotide'].isin(nucleotides)],
            gap=.1,
            order=nucleotides,
            hue_order=[str(stats1.label), str(stats2.label)],
            density_norm='width',
            palette=HuePalette(),
            cut=0
        )

    elif plot_type == 'boxen':
        sns.boxenplot(
            data=df[df['Nucleotide'].isin(nucleotides)],
            x='Nucleotide',
            y='Frequency',
            hue="label",
            order=nucleotides,
            hue_order=[str(stats1.label), str(stats2.label)],
            ax=ax,
            palette=HuePalette(),
            width=0.8
        )

    else:
        raise ValueError(f"Unknown plot type: {plot_type}. Supported types: 'violin', 'boxen'")

    ax.set_ylim(-0.1, 1.1)

    # Reserve the underline margin so flagged/unflagged plots share geometry
    reserve_flag_margin(ax, len(nucleotides))

    ax.set_xlabel('Nucleotide', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax = prepare_legend(ax, box_to_anchor=(0.5, -0.15))

    return fig

def mark_failed_nucleotides(fig: Figure, nucleotides: list[str],
                            failed_nucleotides: dict[str, str]) -> None:
    """Underline the flagged nucleotides on a figure `plot_nucleotides` drew.

    Args:
        fig: The figure to mark, modified in place.
        nucleotides: The nucleotides it was drawn for, in the order it drew them.
        failed_nucleotides: Dict mapping nucleotide -> flag status
            (e.g., `{'A': 'Warning', 'G': 'Fail'}`).
    """
    add_failed_outline(fig.axes[0], {index: failed_nucleotides[nt]
                                     for index, nt in enumerate(nucleotides)
                                     if nt in failed_nucleotides})

def plot_dinucleotides(stats1: SequenceStatistics, stats2: SequenceStatistics,
                       nucleotides: list[str], plot_type: str) -> Figure:
    """
    Plot the dinucleotide content of two sets of sequences.

    Args:
        stats1: Statistics for the first set of sequences.
        stats2: Statistics for the second set of sequences.
        nucleotides: List of nucleotides to generate dinucleotides from.
        plot_type: Type of plot to create (e.g., 'boxen', 'violin').

    Returns:
        Matplotlib figure object. `mark_failed_dinucleotides` adds the flags to it.
    """

    df = melt_stats(stats1, stats2, 'Per sequence dinucleotide content',
                    var_name='Dinucleotide', value_name='Frequency')

    fig, axs = plt.subplots(
        len(nucleotides), 1, figsize=(12, len(nucleotides) * 3 + 2), sharey=True, dpi=300,
        squeeze=False
    )
    axs = axs[:, 0]
    for index, nt in enumerate(nucleotides):
        dinucleotides = [nt + nt2 for nt2 in nucleotides]
        row = df[df['Dinucleotide'].isin(dinucleotides)]

        if plot_type == 'violin':
            sns.violinplot(
                x='Dinucleotide',
                y='Frequency',
                hue="label",
                split=True,
                data=row,
                gap=.1,
                order=dinucleotides,
                hue_order=[str(stats1.label), str(stats2.label)],
                ax=axs[index],
                density_norm='width',
                palette=HuePalette(),
                cut=0
            )

        elif plot_type == 'boxen':
            sns.boxenplot(
                x='Dinucleotide',
                y='Frequency',
                hue="label",
                data=row,
                order=dinucleotides,
                hue_order=[str(stats1.label), str(stats2.label)],
                ax=axs[index],
                palette=HuePalette(),
                width=0.8
            )
        else:
            logger.error(f"Unknown plot type: {plot_type}")

        axs[index].set_ylim(-0.1, 1.1)
        # Reserve the underline margin so flagged/unflagged plots share geometry
        reserve_flag_margin(axs[index], len(dinucleotides))

        axs[index].set_xlabel('')
        axs[index].legend().set_visible(False)
        axs[index].set_ylabel('Frequency', fontsize=14)
        axs[index].tick_params(axis='x', labelsize=12)
        axs[index].tick_params(axis='y', labelsize=12)

    axs[index] = prepare_legend(axs[index])
    axs[index].set_xlabel('Dinucleotide', fontsize=14)

    return fig

def mark_failed_dinucleotides(fig: Figure, nucleotides: list[str],
                              failed_dinucleotides: dict[str, str]) -> None:
    """Underline the flagged dinucleotides on a `plot_dinucleotides` figure.

    One panel per first base, which is the order `plot_dinucleotides` lays them
    out in, so a panel's axis is `fig.axes[index]` for that base.

    Args:
        fig: The figure to mark, modified in place.
        nucleotides: The nucleotides it was drawn for, in the order it drew them.
        failed_dinucleotides: Dict mapping dinucleotide -> flag status
            (e.g., `{'AA': 'Warning', 'GG': 'Fail'}`).
    """
    for index, nt in enumerate(nucleotides):
        dinucleotides = [nt + nt2 for nt2 in nucleotides]
        add_failed_outline(fig.axes[index], {i: failed_dinucleotides[dn]
                                             for i, dn in enumerate(dinucleotides)
                                             if dn in failed_dinucleotides})

def plot_one_stat(stats1, stats2, stats_name, plot_type, x_label='', title=''):
    """
    Plot a single statistic from two stats objects.
    """

    # make dataframe with two columns: label and values
    df1 = stats1.stats[stats_name]
    df2 = stats2.stats[stats_name]
    df = pd.concat([
        df1.assign(label=str(stats1.label)),
        df2.assign(label=str(stats2.label))
    ], ignore_index=True)

    min_y = df[stats_name].min()
    max_y = df[stats_name].max()

    fig, ax = plt.subplots(1, 1, figsize=(6, 6), dpi=300)
    if plot_type == 'violin':
        sns.violinplot(
            y=stats_name,
            hue="label",
            split=True,
            data=df,
            gap=.1,
            hue_order=[str(stats1.label), str(stats2.label)],
            ax=ax,
            density_norm='width',
            palette=HuePalette(),
            cut=0
        )
        if min_y != max_y:
            ax.set_ylim(min_y - 0.1 * abs(max_y - min_y), max_y + 0.1 * abs(max_y - min_y))
    elif plot_type == 'boxen':
        sns.boxenplot(
            y=stats_name,
            hue="label",
            data=df,
            hue_order=[str(stats1.label), str(stats2.label)],
            ax=ax,
            palette=HuePalette(),
            width=0.8
        )
        if min_y != max_y:
            ax.set_ylim(min_y - 0.1 * abs(max_y - min_y), max_y + 0.1 * abs(max_y - min_y))
    else:
        logger.error(f"Unknown plot type: {plot_type}")

    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(stats_name, fontsize=14)
    # remove x ticks
    ax.set_xticks([])
    ax.tick_params(axis='y', labelsize=12)
    if min_y == max_y:
        # show only one value
        ax.set_yticks([min_y])
    ax.ticklabel_format(axis='y', style='plain')
    ax = prepare_legend(ax, box_to_anchor=(0.5, -0.05))

    return fig


def plot_per_base_sequence_comparison(stats1, stats2, stats_name, nucleotides, end_position,
                                      x_label='', title=''):
    """Plot the per-position composition of two classes, base by base.

    The curves run to `end_position`, which callers set to the compared window:
    every position drawn is a position that was scored, so an unflagged stretch
    reads as a stretch that passed. The one exception is a comparison that scored
    nothing, where the caller passes the reported window instead and no position
    is flagged at all. The coverage panel below says how much data stands behind
    each position; where the window ends is said in words in the report, not as a
    line across that panel.

    Args:
        stats1, stats2: Statistics objects for two datasets.
        stats_name: Name of the stats dataframe to plot.
        nucleotides: List of nucleotides to plot (e.g., ['A', 'C', 'G', 'T']).
        end_position: Maximum position to plot.
        x_label, title: Plot labels.

    Returns:
        Matplotlib figure object. `mark_failed_positions` adds the flags to it.
    """

    df1 = stats1.stats[stats_name]
    df2 = stats2.stats[stats_name]

    fig, axs = plt.subplots(
        len(nucleotides) + 1, 1,
        figsize=(12, len(nucleotides) * 2 + 2),
        height_ratios=[1] * len(nucleotides) + [0.5],
        sharey=True, dpi=300
    )

    for index, nt in enumerate(nucleotides):

        df1_base = df1[nt][:end_position]
        df2_base = df2[nt][:end_position]
        # +1 so the position starts from 1
        axs[index].plot(df1.index[:end_position] + 1, df1_base, label=f"{stats1.label}",
                        color=HuePalette()[0], alpha=0.7)
        axs[index].plot(df2.index[:end_position] + 1, df2_base, label=f"{stats2.label}",
                        color=HuePalette()[1], alpha=0.7)

        axs[index].set_ylim(-0.1, 1.1)
        axs[index].set_ylabel('Frequency', fontsize=14)
        axs[index].legend().set_visible(False)
        axs[index].tick_params(axis='x', labelsize=12)
        axs[index].tick_params(axis='y', labelsize=12)
        axs[index].ticklabel_format(axis='both', style='plain')

        # set x ticks
        ticks = [1]
        step = max(1, end_position // 10)
        ticks.extend(range(step, end_position + 1, step))
        axs[index].set_xticks(ticks)

        # add text to the plot with the nucleotide name
        axs[index].text(0.9, 0.8, f'Nucleotide: {nt}', ha='center', va='bottom', fontsize=14,
                        transform=axs[index].transAxes)

        if title and index == 0:
            axs[index].set_title(f'{title}', fontsize=16)

    # The proportion of sequences reaching each position is the denominator behind
    # the curves above, since a position is scored only on the sequences that have
    # it. One curve per class rather than one pooled curve: a position is compared
    # only where *both* classes clear the coverage floor, so the class that runs
    # out first is what ends the compared window, and pooling hides it. Positions
    # are 1-based here to line up with the x axis of the panels above; counting
    # from 0 would shift the whole curve one position right and start it at a
    # trivial 100%.
    positions = list(range(1, end_position + 1))
    # Normalized against every sequence in the class, not against the first
    # plotted position, so the axis reads as a proportion of the class even when
    # sequences are empty. `coverage_curve` is what the interactive viewer draws
    # from as well, so the two figures cannot disagree.
    coverage = [stats.coverage_curve(end_position) for stats in (stats1, stats2)]
    # The cohort a comparison at a position actually has is the smaller of the two.
    binding = np.minimum(coverage[0], coverage[1])

    # plot coverage in the last subplot
    last_index = len(nucleotides)
    axs[last_index].fill_between(positions, binding, color='lightgray', alpha=0.5)
    axs[last_index].plot(positions, coverage[0], color=HuePalette()[0], alpha=0.7, linewidth=2)
    axs[last_index].plot(positions, coverage[1], color=HuePalette()[1], alpha=0.7, linewidth=2)
    axs[last_index].set_xlabel(f"{x_label}", fontsize=14)
    axs[last_index].set_ylabel('Proportion of\neach class\nreaching position', fontsize=14)
    axs[last_index].yaxis.set_label_position("right")
    axs[last_index].yaxis.tick_right()
    axs[last_index].set_ylim(0, 1.1)
    axs[last_index].set_yticks([0, 0.5, 1])
    axs[last_index].set_yticklabels(['0', '0.5', '1'], fontsize=12)
    axs[last_index].tick_params(axis='x', labelsize=12)
    axs[last_index].tick_params(axis='y', labelsize=12)

    # set x ticks
    ticks = [1]
    ticks.extend(range(max(1, end_position // 10), end_position + 1, max(1, end_position // 10)))
    axs[last_index].set_xticks(ticks)
    axs[last_index] = prepare_legend(axs[index], box_to_anchor=(0.5, -1))

    return fig

def mark_failed_positions(fig, nucleotides, end_position, failed_positions):
    """Shade the flagged positions on a `plot_per_base_sequence_comparison` figure.

    One panel per nucleotide, which is the order the figure lays them out in, so
    a panel's axis is `fig.axes[index]` for that nucleotide. The coverage panel
    below them takes no bands: it says how much data reaches a position, which
    is the same whether or not the position was flagged.

    The bands go behind the curves whenever they are added, because `axvspan`
    draws at a lower z-order than a line.

    Args:
        fig: The figure to mark, modified in place.
        nucleotides: The nucleotides it was drawn for, in the order it drew them.
        end_position: The last position it drew; a flag past it has no panel to
            sit on and is left off.
        failed_positions: Dict mapping nucleotide -> dict of position -> flag status
            (e.g., {'A': {52: 'Warning', 66: 'Fail'}, 'G': {70: 'Fail'}}).
    """
    for index, nt in enumerate(nucleotides):
        for position, flag in failed_positions.get(nt, {}).items():
            if 1 <= position <= end_position:
                color = FAIL_COLOR if flag == 'Fail' else WARN_COLOR
                # axvspan uses inclusive bounds, so shade from pos-0.5 to pos+0.5
                # for a 1-wide band
                fig.axes[index].axvspan(position - 0.5, position + 0.5,
                                        color=color, alpha=0.5, linewidth=0)

def _compute_duplication_bins(stats1, stats2):
    """Compute normalized duplication bin distributions for both datasets.

    Returns:
        List of (bin_dict, stats) tuples.
    """
    bins_list = []

    for stats in [stats1, stats2]:
        # count_distribution_bins for bins 1, 2, 3 .. >10, >50, >100, >500, >1000
        count_distribution_bins = {
            1: 0,
            2: 0,
            3: 0,
            4: 0,
            5: 0,
            6: 0,
            7: 0,
            8: 0,
            9: 0,
            10: 0,
            '>10': 0,
            '>50': 0,
            '>100': 0,
            '>500': 0,
            '>1000': 0
        }
        for v in stats.stats['Sequence duplication levels'].values():
            if v == 0:
                continue
            if v <= 10:
                count_distribution_bins[v] = count_distribution_bins.get(v, 0) + v
            elif v <= 50:
                count_distribution_bins['>10'] = count_distribution_bins.get('>10', 0) + v
            elif v <= 100:
                count_distribution_bins['>50'] = count_distribution_bins.get('>50', 0) + v
            elif v <= 500:
                count_distribution_bins['>100'] = count_distribution_bins.get('>100', 0) + v
            elif v <= 1000:
                count_distribution_bins['>500'] = count_distribution_bins.get('>500', 0) + v
            else:
                count_distribution_bins['>1000'] = count_distribution_bins.get('>1000', 0) + v

        # normalize counts to frequencies
        total = stats1.stats['Number of sequences'] + stats2.stats['Number of sequences']
        if total > 0:
            for k in count_distribution_bins:
                count_distribution_bins[k] /= total
                count_distribution_bins[k] *= 100  # convert to percentage

        bins_list.append((count_distribution_bins, stats))

    return bins_list


def plot_sequence_duplications_within_classes(stats1, stats2, percent_remaining_after_dedup=None):
    """Plot how much of each class is made up of sequences duplicated N times."""

    fig, ax = plt.subplots(figsize=(12, 4), dpi=300)
    palette = HuePalette()

    bins_list = _compute_duplication_bins(stats1, stats2)

    # Extract bin keys and values for both stats
    bin_keys = list(bins_list[0][0].keys())
    values1 = [bins_list[0][0][k] for k in bin_keys]
    values2 = [bins_list[1][0][k] for k in bin_keys]

    label1 = str(bins_list[0][1].label)
    label2 = str(bins_list[1][1].label)
    plot_data = []
    for bin_key, v1, v2 in zip(bin_keys, values1, values2, strict=True):
        plot_data.append({"duplication_bin": str(bin_key), "label": label1, "value": v1})
        plot_data.append({"duplication_bin": str(bin_key), "label": label2, "value": v2})

    df = pd.DataFrame(plot_data)
    sns.barplot(
        data=df,
        x="duplication_bin",
        y="value",
        hue="label",
        hue_order=[label1, label2],
        order=[str(k) for k in bin_keys],
        dodge=True,
        ax=ax,
        palette=[palette[0], palette[1]],
        saturation=1,
        edgecolor='none',
        errorbar=None,
    )

    # Set y-limits with padding
    all_values = [v for bins, _ in bins_list for v in bins.values()]
    if all_values:
        max_y = max(all_values)
        ax.set_ylim(0, max_y * 1.1 if max_y > 0 else 1.0)

    if percent_remaining_after_dedup is not None:
        ax.set_title(
            "Percent of sequences remaining after deduplication: "
            f"{percent_remaining_after_dedup:.2%}")

    ax.set_xlabel('Sequence Duplication Level', fontsize=14)
    ax.set_ylabel('% Total Sequences', fontsize=14)
    ax.tick_params(axis='both', labelsize=12)

    legend_handles = [
        Patch(facecolor=palette[0], label=label1),
        Patch(facecolor=palette[1], label=label2),
    ]
    legend_labels = [label1, label2]

    ax = prepare_legend(
        ax,
        box_to_anchor=(0.5, -0.2),
        legend_handles=legend_handles,
        legend_labels=legend_labels,
    )

    return fig

def melt_stats(stats1, stats2, stats_name, var_name='Metric', value_name='Value'):
    """
    Melt the stats DataFrame to long format and add a label column.
    """
    df1 = stats1.stats[stats_name]
    df1 = df1.melt(value_vars=df1.columns, var_name=var_name, value_name=value_name)

    df2 = stats2.stats[stats_name]
    df2 = df2.melt(value_vars=df2.columns, var_name=var_name, value_name=value_name)

    return pd.concat([
        df1.assign(label=str(stats1.label)),
        df2.assign(label=str(stats2.label))
    ], ignore_index=True)


def prepare_legend(ax, box_to_anchor=(0.5, -0.2), legend_handles=None, legend_labels=None):
    """
    Prepare the legend for the plot.
    """
    if legend_handles is None:
        legend_handles = ax.get_legend_handles_labels()[0]
    if legend_labels is None:
        legend_labels = ax.get_legend_handles_labels()[1]
    ax.legend(
        handles = legend_handles,
        labels = legend_labels,
        # A number, not a string: matplotlib takes either a point size or one of
        # its size keywords here, so a quoted '12' is silently dropped and every
        # legend falls back to the 10pt rcParams default.
        fontsize=12,
        loc='upper center',
        bbox_to_anchor=box_to_anchor,
        ncol=3,
    )
    return ax

class HuePalette:
    """The two-color palette shared by every plot, so the classes look the same everywhere.

    Instantiating returns the seaborn palette itself rather than an instance,
    built once and reused: `HuePalette()[0]` is the first class, `[1]` the second.
    """

    _palette = None

    def __new__(palette):
        """Return the shared palette, building it on first use."""
        if palette._palette is None:
            palette._palette = sns.color_palette(list(CLASS_COLORS))

        return palette._palette

