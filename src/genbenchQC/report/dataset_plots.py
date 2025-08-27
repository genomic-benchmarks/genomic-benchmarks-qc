import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import logging

def plot_lengths(stats1, stats2, plot_type='boxen', result=None, dist_thresh=None):
    """
    Plot the sequence lengths of two sequences.
    """
    return plot_one_stat(
        stats1, stats2, 
        'Sequence lengths', 
        plot_type=plot_type,
        x_label='Sequence length', 
        title='Sequence Length Distribution', 
        result=result, 
        dist_thresh=dist_thresh
    )

def plot_gc_content(stats1, stats2, plot_type='boxen', result=None, dist_thresh=None):
    """
    Plot the GC content of two sequences.
    """
    return plot_one_stat(
        stats1, stats2, 
        'Per sequence GC content', 
        x_label='GC content', 
        title='GC Content Distribution',
        plot_type=plot_type,
        result=result, 
        dist_thresh=dist_thresh, 
    )

def plot_nucleotides(stats1, stats2, nucleotides, plot_type, result=None, dist_thresh=None):
    """
    Plot the nucleotide content of two sets of sequences.

    @param stats1: Statistics for the first set of sequences.
    @param stats2: Statistics for the second set of sequences.
    @param nucleotides: List of nucleotides to plot.
    @param plot_type: Type of plot to create (e.g., 'boxen', 'violin').
    @param result: Optional results from a previous analysis. Tuple (distances, passed), where distances is a dict with nucleotides as keys
                   and distances between their distributions as values. Passed is bool indicating if the test passed.
                   If provided, nucleotide plots with distance > dist_thresh will be flagged.
    @param dist_thresh: Optional distance threshold for flagging significant differences.
    @return: Matplotlib figure object.
    """

    df = melt_stats(stats1, stats2, 'Per sequence nucleotide content', var_name='Nucleotide', value_name='Frequency')
    min_y = df['Frequency'].min()
    max_y = df['Frequency'].max()

    fig, ax = plt.subplots(1, 1, figsize=(10, 4), dpi=300)
    if plot_type == 'violin':
        sns.violinplot(
            x='Nucleotide', 
            y='Frequency', 
            hue="label", 
            split=True, 
            data=df[df['Nucleotide'].isin(nucleotides)],
            gap=.1,
            order=nucleotides,
            hue_order=[stats1.label, stats2.label],
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
            hue_order=[stats1.label, stats2.label],
            ax=ax,
            palette=HuePalette(),
            width=0.8
        )
    else:
        logging.error(f"Unknown plot type: {plot_type}")

    red_flag = False
    if result and dist_thresh:
        for index, nt in enumerate(nucleotides):
            if result[0][nt] > dist_thresh:
                red_flag = True
                # draw a red rectangle around the violins
                ax.add_patch(make_red_flag_rectangle(index, min_y, max_y))

    ax.set_title('Nucleotide content', fontsize=16)
    ax.set_xlabel('Nucleotide', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax = prepare_legend(ax, red_flag, dist_thresh)

    return fig

def plot_dinucleotides(stats1, stats2, nucleotides, plot_type, result=None, dist_thresh=None):

    df = melt_stats(stats1, stats2, 'Per sequence dinucleotide content', var_name='Dinucleotide', value_name='Frequency')
    min_y = df['Frequency'].min()
    max_y = df['Frequency'].max()
    
    fig, axs = plt.subplots(len(nucleotides), 1, figsize=(10, len(nucleotides) * 3 + 2), sharey=True, dpi=300)
    red_flag = False
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
                hue_order=[stats1.label, stats2.label],
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
                hue_order=[stats1.label, stats2.label],
                ax=axs[index],
                palette=HuePalette(),
                width=0.8
            )
        else:
            logging.error(f"Unknown plot type: {plot_type}")
        
        if index == 0:
            axs[index].set_title('Dinucleotide content', fontsize=16)

        axs[index].set_xlabel('')
        axs[index].legend().set_visible(False)
        axs[index].set_ylabel('Frequency', fontsize=14)
        axs[index].tick_params(axis='x', labelsize=12)
        axs[index].tick_params(axis='y', labelsize=12)

        if result and dist_thresh:
            for di_index, dint in enumerate(dinucleotides):
                if result[0][dint] > dist_thresh:
                    red_flag = True
                    # draw a red rectangle around the violins, put it behind the violins
                    axs[index].add_patch(make_red_flag_rectangle(di_index, min_y, max_y))

    axs[index] = prepare_legend(axs[index], red_flag, dist_thresh)
    axs[index].set_xlabel('Dinucleotide', fontsize=14)

    return fig

def plot_one_stat(stats1, stats2, stats_name, plot_type, x_label='', title='', result=None, dist_thresh=None):
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

    fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=300)
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
    else:
        logging.error(f"Unknown plot type: {plot_type}")
    
    # result is a tuple of (distances, passed)
    red_flag = False
    if result and dist_thresh:
        if result[0] > dist_thresh:
            red_flag = True
            # draw a red rectangle around the violins
            ax.add_patch(make_red_flag_rectangle(0, min_y, max_y))

    ax.set_title(title, fontsize=16)
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(stats_name, fontsize=14)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    if min_y == max_y:
        # show only one value
        ax.set_yticks([min_y])
    ax.ticklabel_format(axis='y', style='plain')
    ax = prepare_legend(ax, red_flag, dist_thresh) 

    return fig

def plot_per_base_sequence_comparison(stats1, stats2, stats_name, nucleotides, end_position, x_label='', title='', result=None, p_value_thresh=None):

    df1 = stats1.stats[stats_name]
    df2 = stats2.stats[stats_name]

    fig, axs = plt.subplots(
        len(nucleotides) + 1, 1, 
        figsize=(10, len(nucleotides) * 2 + 2), 
        height_ratios=[1] * len(nucleotides) + [0.5],
        sharey=True, dpi=300
    )

    red_flag = False
    for index, nt in enumerate(nucleotides):

        df1_base = df1[nt][:end_position]
        df2_base = df2[nt][:end_position]
        axs[index].plot(df1.index[:end_position], df1_base, label=f"{stats1.label}", color=HuePalette()[0], alpha=0.7)
        axs[index].plot(df2.index[:end_position], df2_base, label=f"{stats2.label}", color=HuePalette()[1], alpha=0.7)

        axs[index].set_ylim(-0.1, 1.1)
        axs[index].set_ylabel('Frequency', fontsize=14)
        axs[index].legend().set_visible(False)
        axs[index].tick_params(axis='x', labelsize=12)
        axs[index].tick_params(axis='y', labelsize=12)

        # add text to the plot with the nucleotide name
        axs[index].text(0.9, 0.8, f'Nucleotide: {nt}', ha='center', va='bottom', fontsize=14, transform=axs[index].transAxes)

        if title and index == 0:
            axs[index].set_title(f'{title}', fontsize=16)


        if result and p_value_thresh:
            for i in range(end_position):
                
                p_value = result[0][nt][i]
                if p_value < p_value_thresh:
                    axs[index].axvspan(i-0.45, i+0.45, facecolor='red', alpha=0.2)
                    red_flag = True

    axs[index].ticklabel_format(axis='both', style='plain')

    seq_lengths = list(stats1.stats['Sequence lengths'].values.flatten()) + list(stats2.stats['Sequence lengths'].values.flatten())
    # Plot the number of sequences with length at least that position
    length_counts = [sum(1 for length in seq_lengths if length >= pos) for pos in range(end_position)]
    # normalize length_counts to [0, 1]
    if length_counts:
        length_counts = [count / max(length_counts) for count in length_counts]

    last_index = len(nucleotides)
    axs[last_index].fill_between(range(end_position), length_counts, color='lightblue', alpha=0.5)
    axs[last_index].plot(range(end_position), length_counts, color='lightblue', linewidth=2)
    axs[last_index].set_xlabel(f"{x_label}", fontsize=14)
    axs[last_index].set_ylabel('Proportion of\nsequences', fontsize=14)
    axs[last_index].yaxis.set_label_position("right")
    axs[last_index].yaxis.tick_right()
    axs[last_index].set_ylim(0, max(length_counts) * 1.1)  
    axs[last_index].set_yticks([0, 0.5, 1])
    axs[last_index].set_yticklabels(['0', '0.5', '1'], fontsize=12)
    axs[last_index].tick_params(axis='x', labelsize=12)
    axs[last_index].tick_params(axis='y', labelsize=12)

    axs[last_index] = prepare_legend(axs[index], red_flag, p_value_thresh, box_to_anchor=(0.5, -1), metric='p-value <')

    return fig

def melt_stats(stats1, stats2, stats_name, var_name='Metric', value_name='Value', keep_positions=False):
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

def prepare_legend(ax, red_flag, dist_thresh, box_to_anchor=(0.5, -0.2), metric='Distance >'):
    """
    Prepare the legend for the plot.
    """
    legend_handles = ax.get_legend_handles_labels()[0]
    legend_labels = ax.get_legend_handles_labels()[1]
    if red_flag:
        legend_handles += [plt.Rectangle((0, 0), 1, 1, color='red', alpha=0.2)]
        legend_labels += [f'{metric} {dist_thresh}']
    ax.legend(
        handles = legend_handles,
        labels = legend_labels,
        title='Label',
        title_fontsize='14',
        fontsize='12',
        loc='upper center', 
        bbox_to_anchor=box_to_anchor, 
        ncol=3,
    )
    return ax

def make_red_flag_rectangle(index, min_y, max_y, margin=0.02):
    """
    Create a red rectangle to highlight the violins that are above the threshold.
    """
    flag_box_width = 1 - 2 * margin
    flag_box_hight = max_y - min_y + 2 * margin
    return plt.Rectangle((index - 0.5 + margin, min_y - margin), flag_box_width, flag_box_hight, color='red', alpha=0.2, zorder=-1)

class HuePalette:

    _palette = None

    def __new__(palette):
        if palette._palette is None:
            palette._palette = sns.color_palette()[:2]

        return palette._palette
                
