import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd
import os
from pathlib import Path
import logging

from genbenchQC.report.sequence_html_report import get_sequence_html_template
from genbenchQC.report.dataset_html_report import get_dataset_html_template
from genbenchQC.report.split_html_report import get_train_test_html_template
from genbenchQC.utils.input_utils import write_stats_json
from genbenchQC.report import dataset_plots
from genbenchQC.report import sequences_plots

def generate_sequence_plots(stats_dict, output_path, end_position, plot_type='boxen'):
    """
    Generate a plots from the given statistics dictionary.
    """

    logging.info(f"Generating PNG plots at: {output_path}")

    nucleotides = sorted(stats_dict['Unique bases'])

    plots_paths = {}
    fig = sequences_plots.hist_plot_one_stat(
        stats_dict,
        stats_name='Sequence lengths',
        x_label='Sequence lengths',
        title='Sequence lengths'
    )

    plots_paths['Sequence lengths'] = Path(output_path.name) / 'sequence_lengths.png'
    fig.savefig(output_path / 'sequence_lengths.png', bbox_inches='tight')
    plt.close(fig)

    fig = sequences_plots.hist_plot_one_stat(
        stats_dict,
        stats_name='Per sequence GC content',
        x_label='GC content (%)',
        title='Per sequence GC content'
    )
    plots_paths['Per sequence GC content'] = Path(output_path.name) / 'per_sequence_gc_content.png'
    fig.savefig(output_path / 'per_sequence_gc_content.png', bbox_inches='tight')
    plt.close(fig)

    fig = sequences_plots.plot_nucleotides(stats_dict, nucleotides=nucleotides, plot_type=plot_type)
    plots_paths['Per sequence nucleotide content'] = Path(output_path.name) / 'per_sequence_nucleotide_content.png'
    fig.savefig(output_path / 'per_sequence_nucleotide_content.png', bbox_inches='tight')
    plt.close(fig)  

    fig = sequences_plots.plot_dinucleotides(stats_dict, nucleotides=nucleotides, plot_type=plot_type)
    plots_paths['Per sequence dinucleotide content'] = Path(output_path.name) / 'per_sequence_dinucleotide_content.png'
    fig.savefig(output_path / 'per_sequence_dinucleotide_content.png', bbox_inches='tight')
    plt.close(fig)

    fig = sequences_plots.plot_per_position_nucleotide_content(
        stats_dict,
        stat_name='Per position nucleotide content',
        nucleotides=nucleotides,
        end_position=end_position,
    )
    plots_paths['Per position nucleotide content'] = Path(output_path.name) / 'per_position_nucleotide_content.png'
    fig.savefig(output_path / 'per_position_nucleotide_content.png', bbox_inches='tight')
    plt.close(fig)

    fig = sequences_plots.plot_per_position_nucleotide_content(
        stats_dict,
        stat_name='Per position reversed nucleotide content',
        nucleotides=nucleotides,
        end_position=end_position,
    )
    plots_paths['Per position reversed nucleotide content'] = Path(output_path.name) / 'per_position_reversed_nucleotide_content.png'
    fig.savefig(output_path / 'per_position_reversed_nucleotide_content.png', bbox_inches='tight')
    plt.close(fig)

    return plots_paths

def generate_sequence_html_report(stats_dict, output_path, plots_path, end_position, plot_type):
    """
    Generate an HTML report from the given statistics dictionary.
    Plots are generated using the Plotly library.
    """
    logging.info(f"Generating HTML report: {output_path}")

    plots_path.mkdir(parents=True, exist_ok=True)

    # generate 
    plots_paths = generate_sequence_plots(stats_dict, plots_path, end_position=end_position, plot_type=plot_type)

    # Load the HTML template
    template = get_sequence_html_template(stats_dict, plots_paths)

    with open(output_path, 'w') as file:
        file.write(template)

def generate_train_test_html_report(clusters, train_filename, train_seq, test_filename, test_seq, output_path, identity_threshold, alignment_coverage):
    """
    Create an HTML report that lists mixed clusters between training and test sets.
    
    Parameters:
        clusters (Sequence[Mapping]): Collection of cluster records that include membership from both datasets.
        train_filename (str): Filename or label for the training dataset.
        train_seq (Mapping): Mapping of training sequence identifiers to sequence data or metadata.
        test_filename (str): Filename or label for the test dataset.
        test_seq (Mapping): Mapping of test sequence identifiers to sequence data or metadata.
        output_path (str or Path): Path to write the rendered HTML report.
        identity_threshold (float): Sequence identity threshold used to define cluster membership.
        alignment_coverage (float): Alignment coverage threshold used to define cluster membership.
    """
    # Load the HTML template
    template = get_train_test_html_template(clusters, train_filename, train_seq, test_filename, test_seq, identity_threshold, alignment_coverage)

    with open(output_path, 'w') as file:
        file.write(template)

def generate_dataset_html_report(stats1, stats2, output_path, plots_path, end_position, plot_type, results=None, threshold=None):
    """
    Generate an HTML report that compares two dataset statistics and writes the report (and any duplicate sequences) to disk.
    
    Parameters:
        stats1: Statistics object for the first dataset containing at least `.sequences` and metadata used by the template.
        stats2: Statistics object for the second dataset containing at least `.sequences` and metadata used by the template.
        output_path (str | Path): File path where the rendered HTML report will be written.
        plots_path (str | Path): Directory where comparison plots will be generated and saved.
        end_position (int): End position index used when generating per-position plots.
        plot_type (str): Plot style to use when generating figures (for example, 'boxen').
        results (optional): Optional results object used to annotate plots (may be None).
        threshold (optional): Optional threshold value passed to plotting routines (may be None).
    
    Notes:
        - This function creates plots in `plots_path`, renders an HTML template that embeds those plots, and writes the HTML to `output_path`.
        - If any sequence labels appear in both `stats1` and `stats2`, those duplicate sequence identifiers are written to a separate file next to `output_path` with a "_duplicates.txt" suffix.
    """
    plots_path.mkdir(parents=True, exist_ok=True)

    # generate 
    plots_paths = generate_dataset_plots(stats1, stats2, plots_path, end_position, plot_type, results=None, threshold=None)

    # find duplicate sequences between labels
    duplicate_seqs = list(set(stats1.sequences).intersection(stats2.sequences))
    # Load the HTML template
    template = get_dataset_html_template(stats1, stats2, plots_paths, duplicate_seqs)

    with open(output_path, 'w') as file:
        file.write(template)

    if len(duplicate_seqs) > 0:
        # remove extension from output path, add '_duplicates.txt'
        duplicate_seqs_path = os.path.splitext(output_path)[0] + '_duplicates.txt'
        with open(duplicate_seqs_path, 'w') as f:
            for seq in duplicate_seqs:
                f.write(f"{seq}\n")
        logging.info(f"Duplicate sequences saved to {duplicate_seqs_path}")

def generate_json_report(stats_dict, output_path):
    write_stats_json(stats_dict, output_path)

def generate_simple_report(results, output_path):

    """
    Write a simple CSV report from the provided results.
    
    If `results` is a dict, it is converted to a pandas DataFrame using orient='index' before writing. The resulting DataFrame is written to `output_path` as a CSV file.
    
    Parameters:
        results (dict or pandas.DataFrame): Data to write to CSV. If a dict, keys become rows (index).
        output_path (str or pathlib.Path): Destination file path for the CSV.
    """
    logging.info(f"Generating simple report: {output_path}")

    if isinstance(results, dict):
        results = pd.DataFrame.from_dict(results, orient='index')
    results.to_csv(output_path)

def generate_dataset_plots(stats1, stats2, output_path, end_position, plot_type='boxen', results=None, threshold=None):

    """
    Generate a set of comparison PNG plots for two sequence datasets and return their file paths.
    
    Parameters:
        stats1: An object containing summary statistics for the first dataset (must expose .stats with 'Unique bases').
        stats2: An object containing summary statistics for the second dataset (must expose .stats with 'Unique bases').
        output_path (Path): Directory where generated PNG files will be written.
        end_position (int): Last sequence position to include in per-position plots.
        plot_type (str): Plot style to use for distribution plots (default: 'boxen').
        results (dict, optional): Optional dictionary of prior statistical results keyed by plot name; specific entries
            (e.g., 'Per sequence nucleotide content') will be passed to plotting helpers when present.
        threshold (float, optional): Optional significance or distance threshold forwarded to plotting helpers.
    
    Returns:
        dict: Mapping from descriptive plot names to Path objects pointing to the generated PNG files (paths are based on output_path.name).
    """
    logging.info(f"Generating PNG plots at: {output_path}")

    plots_paths = {}

    bases_overlap = list(set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases']))

    # Plot per sequence nucleotide content
    fig = dataset_plots.plot_nucleotides(
        stats1,
        stats2,
        nucleotides = bases_overlap,
        result=results['Per sequence nucleotide content'] if results else None,
        dist_thresh=threshold,
        plot_type=plot_type
    )
    plots_paths['Per sequence nucleotide content'] = Path(output_path.name) / 'per_sequence_nucleotide_content.png'
    fig.savefig(output_path / 'per_sequence_nucleotide_content.png', bbox_inches='tight')
    plt.close(fig)

    # Plot per sequence dinucleotide content
    fig = dataset_plots.plot_dinucleotides(
        stats1,
        stats2,
        nucleotides = bases_overlap,
        result=results['Per sequence dinucleotide content'] if results else None,
        dist_thresh=threshold,
        plot_type=plot_type
    )
    plots_paths['Per sequence dinucleotide content'] = Path(output_path.name) / 'per_sequence_dinucleotide_content.png'
    fig.savefig(output_path / 'per_sequence_dinucleotide_content.png', bbox_inches='tight')
    plt.close(fig)
    
    # Plot per position nucleotide content
    fig = dataset_plots.plot_per_base_sequence_comparison(
        stats1,
        stats2,
        stats_name='Per position nucleotide content',
        result=results['Per position nucleotide content'] if results else None,
        p_value_thresh=threshold,
        nucleotides = bases_overlap,
        end_position=end_position,
        x_label='Position in sequence',
        title='Nucleotide composition per position',
    )
    plots_paths['Per position nucleotide content'] = Path(output_path.name) / 'per_position_nucleotide_content.png'
    fig.savefig(output_path / 'per_position_nucleotide_content.png', bbox_inches='tight')   
    plt.close(fig)

    # Plot per reversed position nucleotide content
    fig = dataset_plots.plot_per_base_sequence_comparison(
        stats1,
        stats2,
        stats_name='Per position reversed nucleotide content',
        result=results['Per position reversed nucleotide content'] if results else None,
        p_value_thresh=threshold,
        nucleotides = bases_overlap,
        end_position=end_position,
        x_label='Position in reversed sequence',
        title='Reversed nucleotide composition per position',
    )
    plots_paths['Per position reversed nucleotide content'] = Path(output_path.name) / 'per_position_reversed_nucleotide_content.png'
    fig.savefig(output_path / 'per_position_reversed_nucleotide_content.png', bbox_inches='tight')
    plt.close(fig)

    # Plot length distribution
    fig = dataset_plots.plot_lengths(
        stats1,
        stats2,
        result=results['Sequence lengths'] if results else None,
        dist_thresh=threshold,
        plot_type=plot_type,
    )
    plots_paths['Sequence lengths'] = Path(output_path.name) / 'sequence_lengths.png'
    fig.savefig(output_path / 'sequence_lengths.png', bbox_inches='tight')
    plt.close(fig)

    # Plot per sequence GC content
    fig = dataset_plots.plot_gc_content(
        stats1,
        stats2,
        result=results['Per sequence GC content'] if results else None,
        dist_thresh=threshold,
        plot_type=plot_type,
    )
    plots_paths['Per sequence GC content'] = Path(output_path.name) / 'per_sequence_gc_content.png'
    fig.savefig(output_path / 'per_sequence_gc_content.png', bbox_inches='tight')
    plt.close(fig)

    return plots_paths