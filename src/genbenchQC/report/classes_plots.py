import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import logging
from genbenchQC.report.utils import FAIL_COLOR, WARN_COLOR

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

def reserve_flag_margin(ax, n_positions):
    """Reserve the bottom margin used by flag underlines so geometry is stable.

    The colored underlines from add_failed_outline are drawn with clip_on=False,
    which enlarges the tight bounding box used by savefig(bbox_inches='tight').
    Without this reservation, flagged and unflagged plots crop to different
    bounding boxes -> boxen boxes render at different widths/spacing. Drawing
    invisible stubs at every position makes the tight bbox identical whether or
    not any position is actually flagged.

    @param ax: Matplotlib axis object.
    @param n_positions: Number of x categories to reserve space for.
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

def add_failed_outline(ax, failed_positions, x_left_custom=None, x_right_custom=None):
    """Add colored underlines to failed positions.

    @param ax: Matplotlib axis object.
    @param failed_positions: Dict mapping position index -> flag status ('Fail' or 'Warning').
    @param x_left_custom: Optional custom left x-coordinate for the underline.
    @param x_right_custom: Optional custom right x-coordinate for the underline.
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

def plot_nucleotides(stats1, stats2, nucleotides, plot_type, failed_nucleotides=None):
    """
    Plot the nucleotide content of two sets of sequences.

    @param stats1: Statistics for the first set of sequences.
    @param stats2: Statistics for the second set of sequences.
    @param nucleotides: List of nucleotides to plot.
    @param plot_type: Type of plot to create (e.g., 'boxen', 'violin').
    @param failed_nucleotides: Dict mapping nucleotide -> flag status
        (e.g., {'A': 'Warning', 'G': 'Fail'}). Currently not applied as
        seaborn boxenplot requires special handling for outline coloring.
    @return: Matplotlib figure object.
    """
    """
    Plot the nucleotide content of two sets of sequences.

    @param stats1: Statistics for the first set of sequences.
    @param stats2: Statistics for the second set of sequences.
    @param nucleotides: List of nucleotides to plot.
    @param plot_type: Type of plot to create (e.g., 'boxen', 'violin').
    @param failed_nucleotides: Dict mapping nucleotide -> flag status
        (e.g., {'A': 'Warning', 'G': 'Fail'}). Used to add red outlines
        to failed nucleotides in boxen plots.
    @return: Matplotlib figure object.
    """

    df = melt_stats(stats1, stats2, 'Per sequence nucleotide content', var_name='Nucleotide', value_name='Frequency')

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
    # Add colored outlines to failed nucleotides if provided
    if failed_nucleotides:
        add_failed_outline(ax, {i: failed_nucleotides[nt] for i, nt in enumerate(nucleotides) if nt in failed_nucleotides})

    ax.set_xlabel('Nucleotide', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax = prepare_legend(ax, box_to_anchor=(0.5, -0.15))

    return fig

def plot_dinucleotides(stats1, stats2, nucleotides, plot_type, failed_dinucleotides=None):
    """
    Plot the dinucleotide content of two sets of sequences.

    @param stats1: Statistics for the first set of sequences.
    @param stats2: Statistics for the second set of sequences.
    @param nucleotides: List of nucleotides to generate dinucleotides from.
    @param plot_type: Type of plot to create (e.g., 'boxen', 'violin').
    @param failed_dinucleotides: Dict mapping dinucleotide -> flag status
        (e.g., {'AA': 'Warning', 'GG': 'Fail'}). Used to add red outlines
        to failed dinucleotides in boxen plots.
    @return: Matplotlib figure object.
    """

    df = melt_stats(stats1, stats2, 'Per sequence dinucleotide content', var_name='Dinucleotide', value_name='Frequency')

    fig, axs = plt.subplots(len(nucleotides), 1, figsize=(12, len(nucleotides) * 3 + 2), sharey=True, dpi=300)
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
            logging.error(f"Unknown plot type: {plot_type}")

        axs[index].set_ylim(-0.1, 1.1)
        # Reserve the underline margin so flagged/unflagged plots share geometry
        reserve_flag_margin(axs[index], len(dinucleotides))
        # Add colored outlines to failed dinucleotides if provided
        if failed_dinucleotides:
            add_failed_outline(axs[index], {i: failed_dinucleotides[dn] for i, dn in enumerate(dinucleotides) if dn in failed_dinucleotides})

        axs[index].set_xlabel('')
        axs[index].legend().set_visible(False)
        axs[index].set_ylabel('Frequency', fontsize=14)
        axs[index].tick_params(axis='x', labelsize=12)
        axs[index].tick_params(axis='y', labelsize=12)

    axs[index] = prepare_legend(axs[index])
    axs[index].set_xlabel('Dinucleotide', fontsize=14)

    return fig

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
        logging.error(f"Unknown plot type: {plot_type}")

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

def plot_per_base_sequence_comparison(stats1, stats2, stats_name, nucleotides, end_position, x_label='', title='', failed_positions=None):
    """Plot per-base sequence comparison with optional failure shading.

    Args:
        stats1, stats2: Statistics objects for two datasets.
        stats_name: Name of the stats dataframe to plot.
        nucleotides: List of nucleotides to plot (e.g., ['A', 'C', 'G', 'T']).
        end_position: Maximum position to plot.
        x_label, title: Plot labels.
        failed_positions: Dict mapping nucleotide -> dict of position -> flag status
            (e.g., {'A': {52: 'Warning', 66: 'Fail'}, 'G': {70: 'Fail'}}).
            If None or empty, no shading applied.

    Returns:
        Matplotlib figure object.
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
        axs[index].plot(df1.index[:end_position] + 1, df1_base, label=f"{stats1.label}", color=HuePalette()[0], alpha=0.7)
        axs[index].plot(df2.index[:end_position] + 1, df2_base, label=f"{stats2.label}", color=HuePalette()[1], alpha=0.7)

        # Add shading for failed positions with color based on flag type
        if failed_positions and nt in failed_positions:
            for pos, flag in failed_positions[nt].items():
                if 1 <= pos <= end_position:
                    color = FAIL_COLOR if flag == 'Fail' else WARN_COLOR
                    # axvspan uses inclusive bounds, so shade from pos-0.5 to pos+0.5 for 1-wide band
                    axs[index].axvspan(pos - 0.5, pos + 0.5, color=color, alpha=0.5, linewidth=0)

        axs[index].set_ylim(-0.1, 1.1)
        axs[index].set_ylabel('Frequency', fontsize=14)
        axs[index].legend().set_visible(False)
        axs[index].tick_params(axis='x', labelsize=12)
        axs[index].tick_params(axis='y', labelsize=12)
        axs[index].ticklabel_format(axis='both', style='plain')

        # set x ticks
        ticks = [1]
        ticks.extend(range(max(1, end_position // 10), end_position + 1, max(1, end_position // 10)))
        axs[index].set_xticks(ticks)

        # add text to the plot with the nucleotide name
        axs[index].text(0.9, 0.8, f'Nucleotide: {nt}', ha='center', va='bottom', fontsize=14, transform=axs[index].transAxes)

        if title and index == 0:
            axs[index].set_title(f'{title}', fontsize=16)

    seq_lengths = list(stats1.stats['Sequence lengths'].values.flatten()) + list(stats2.stats['Sequence lengths'].values.flatten())
    # Plot the number of sequences with length at least that position
    length_counts = [sum(1 for length in seq_lengths if length >= pos) for pos in range(end_position)]
    # normalize length_counts to [0, 1]
    if length_counts:
        length_counts = [count / max(length_counts) for count in length_counts]

    # plot length counts in the last subplot
    last_index = len(nucleotides)
    axs[last_index].fill_between(range(1, end_position + 1), length_counts, color='lightgray', alpha=0.5)
    axs[last_index].plot(range(1, end_position + 1), length_counts, color='lightgray', linewidth=2)
    axs[last_index].set_xlabel(f"{x_label}", fontsize=14)
    axs[last_index].set_ylabel('Proportion of\nsequences', fontsize=14)
    axs[last_index].yaxis.set_label_position("right")
    axs[last_index].yaxis.tick_right()
    axs[last_index].set_ylim(0, max(length_counts) * 1.1)  
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
        for k, v in stats.stats['Sequence duplication levels'].items():
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
    for bin_key, v1, v2 in zip(bin_keys, values1, values2):
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
    all_values = [v for bins, _ in bins_list for v in bins.values() if v > 0]
    if all_values:
        min_y, max_y = min(all_values), max(all_values)
        if min_y != max_y:
            ax.set_ylim(min_y - 0.1 * abs(max_y - min_y), max_y + 0.1 * abs(max_y - min_y))
        else:
            ax.set_ylim(-0.1, 1.1)

    if percent_remaining_after_dedup is not None:
        ax.set_title(f"Percent of sequences remaining after deduplication: {percent_remaining_after_dedup:.2%}")

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

    df = pd.concat([
        df1.assign(label=str(stats1.label)),
        df2.assign(label=str(stats2.label))
    ], ignore_index=True)

    return df

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
        fontsize='12',
        loc='upper center', 
        bbox_to_anchor=box_to_anchor, 
        ncol=3,
    )
    return ax

class HuePalette:

    _palette = None

    def __new__(palette):
        if palette._palette is None:
            palette._palette = sns.color_palette(['#003D99', '#66A3FF'])

        return palette._palette
                
