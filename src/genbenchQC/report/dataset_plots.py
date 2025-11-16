import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import logging

def plot_lengths(stats1, stats2, plot_type='boxen', result=None, dist_thresh=None):
    """
    Plot comparative sequence-length distributions for two stats objects.
    
    Parameters:
        stats1: First stats object; its `.stats` and `.label` are used to build the plot.
        stats2: Second stats object; its `.stats` and `.label` are used to build the plot.
        plot_type (str): Plot style to use, e.g. 'boxen' or 'violin'.
        result (optional): Metric values used to mark significant differences; when provided with `dist_thresh`, values above the threshold will be highlighted.
        dist_thresh (optional): Threshold for `result` above which parts of the plot are flagged as significant.
    
    Returns:
        figure: Matplotlib Figure containing the sequence length comparison plot.
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
    Plot GC content distributions for two stats objects.
    
    Parameters:
        stats1: Stats-like object containing a 'Per sequence GC content' dataframe and a `label` attribute.
        stats2: Stats-like object containing a 'Per sequence GC content' dataframe and a `label` attribute.
        plot_type (str): Plot style to use; supported values include 'boxen' and 'violin'.
        result (optional): Optional analysis result used to flag significant differences (expects indexable structure where relevant metric values are at result[0]).
        dist_thresh (optional): Numeric threshold; when provided with `result`, values above this threshold are highlighted on the plot.
    
    Returns:
        matplotlib.figure.Figure: Figure containing the GC content comparison plot.
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
    Plot per-sequence nucleotide frequency distributions for two stats objects.
    
    Plots either violin or boxen distributions for the specified nucleotides from stats1 and stats2, optionally highlighting nucleotides whose distance in `result` exceeds `dist_thresh` with a red flag and adding a legend entry.
    
    Parameters:
        nucleotides (Sequence[str]): Ordered list of nucleotide labels to include on the x-axis.
        plot_type (str): Plot style to use; recognized values are 'violin' and 'boxen'.
        result (tuple|None): Optional analysis result where result[0] is a mapping of nucleotide -> distance. If provided together with `dist_thresh`, distances greater than the threshold will be flagged.
        dist_thresh (float|None): Distance threshold used to decide which nucleotides to flag when `result` is provided.
    
    Returns:
        matplotlib.figure.Figure: Figure containing the nucleotide-content plot.
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

    """
    Plot dinucleotide frequency distributions for each nucleotide as stacked subplots comparing two stats objects.
    
    Each subplot shows frequencies for dinucleotides starting with a given nucleotide using either violin or boxen plots colored by the two input labels. If both `result` and `dist_thresh` are provided, dinucleotides whose corresponding value in `result[0]` exceeds `dist_thresh` are highlighted with a translucent red rectangle and an entry is added to the legend.
    
    Parameters:
        stats1: Object providing a `.stats` mapping and a `.label` used to source and label the first dataset.
        stats2: Object providing a `.stats` mapping and a `.label` used to source and label the second dataset.
        nucleotides (list[str]): Ordered list of single-character nucleotides to build dinucleotide groups (each subplot corresponds to one nucleotide as the first base).
        plot_type (str): Plot style to use; supported values are `'violin'` and `'boxen'`. Unknown values are logged as errors.
        result (optional): Mapping-like structure where `result[0][dinucleotide]` yields a numeric metric used for significance highlighting.
        dist_thresh (optional, numeric): Threshold applied to `result[0][dinucleotide]` to determine which dinucleotides receive a red highlight.
    
    Returns:
        matplotlib.figure.Figure: Figure containing one subplot per input nucleotide with dinucleotide frequency comparisons.
    """
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
    Plot the distribution of a named statistic for two stats objects.
    
    Creates a single matplotlib figure showing the distribution of stats_name from both stats1 and stats2 using either a violin or boxen plot. If result and dist_thresh are provided and the reported distance exceeds the threshold, the plot is annotated with a translucent red highlight and the legend is updated to indicate the significance.
    
    Parameters:
        stats1: object
            First stats object containing a mapping-like attribute `stats` and a `label` used for the plot legend.
        stats2: object
            Second stats object containing a mapping-like attribute `stats` and a `label` used for the plot legend.
        stats_name (str):
            Key/name of the statistic in each stats object's `stats` mapping to plot.
        plot_type (str):
            Plot style to use; supported values are `'violin'` and `'boxen'`. Unknown values result in an error being logged.
        x_label (str, optional):
            Label for the x-axis.
        title (str, optional):
            Plot title.
        result (optional):
            Optional tuple-like result where the first element is a distance value to compare against dist_thresh; when provided and greater than dist_thresh, the plot is highlighted.
        dist_thresh (optional):
            Threshold used with `result` to determine whether to add a red highlight for significant difference.
    
    Returns:
        matplotlib.figure.Figure: The created figure containing the plotted statistic.
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

    """
    Plot per-base sequence comparisons for the given nucleotides and include a bottom subplot showing the proportion of sequences at each position.
    
    Parameters:
        stats1: object
            First statistics container; must expose .stats (a mapping of DataFrames) and .label (string).
        stats2: object
            Second statistics container; must expose .stats and .label.
        stats_name: str
            Key in each .stats mapping whose DataFrame contains per-base columns named by `nucleotides`.
        nucleotides: Sequence[str]
            Sequence of nucleotide column names to plot (one subplot per nucleotide).
        end_position: int
            Number of positions (columns/indices) from the start to include in the plots.
        x_label: str, optional
            Label for the x-axis of the bottom subplot.
        title: str, optional
            Title to place above the first nucleotide subplot when provided.
        result: optional
            Optional structure providing per-position p-values accessible as result[0][nt][i]; when provided, positions with p < p_value_thresh are shaded.
        p_value_thresh: float, optional
            Threshold for p-value significance; positions with p-values below this value are highlighted and cause an entry to be added to the legend.
    
    Returns:
        matplotlib.figure.Figure
            Figure containing vertically stacked subplots: one per nucleotide showing frequency curves for each stats object, and a bottom subplot showing the normalized proportion of sequences with length at least each position.
    """
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
    Convert two stats DataFrames into a long-form DataFrame suitable for plotting.
    
    Each input's stats table is taken from statsX.stats[stats_name], melted into long format using the provided `var_name` and `value_name`, then concatenated into a single DataFrame with an added `label` column set to str(statsX.label) for each row.
    
    Parameters:
        stats1: Object with a `stats` mapping and a `label` attribute; its stats table at `stats_name` will be melted.
        stats2: Object with a `stats` mapping and a `label` attribute; its stats table at `stats_name` will be melted.
        stats_name (str): Key name of the stats table to extract from each object's `stats`.
        var_name (str): Column name to use for variable names in the melted output (default 'Metric').
        value_name (str): Column name to use for values in the melted output (default 'Value').
        keep_positions (bool): Present for API compatibility but not used by this implementation.
    
    Returns:
        pandas.DataFrame: Long-form DataFrame containing the melted rows from both inputs with columns
        [var_name, value_name, 'label'] where 'label' identifies the source (stats1 or stats2).
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
                