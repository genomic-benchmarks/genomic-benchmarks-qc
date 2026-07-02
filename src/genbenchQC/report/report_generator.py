import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path
import logging

from genbenchQC.report.sequence_html_report import get_sequence_html_template
from genbenchQC.report.classes_html_report import get_dataset_html_template
from genbenchQC.report.split_html_report import get_splits_html_template
from genbenchQC.utils.input_utils import write_stats_json
from genbenchQC.report import classes_plots
from genbenchQC.report import sequences_plots
from genbenchQC.report import splits_plots

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

def generate_splits_html_report(basic_stats, threshold_stats, results_filt, output_path, plots_dir, query_similarity_max, target_similarity_max):
    """
    Generate an HTML report visualising data leakage. 
    """
    plots_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"Generating HTML report: {output_path}")

    plots_paths_dict = generate_split_plots(query_similarity_max, target_similarity_max, threshold_stats, plots_dir)

    template = get_splits_html_template(basic_stats, threshold_stats, results_filt, plots_paths_dict)
    with open(output_path, 'w') as file:
        file.write(template)
        
def generate_split_plots(query_similarity_max, target_similarity_max, threshold_stats, plots_dir):

    plots_paths_dict = {}

    fig = splits_plots.plot_similarity_histograms(query_similarity_max, target_similarity_max, threshold_stats)
    plots_paths_dict['Similarity histograms'] = plots_dir / 'similarity_histograms.png'
    fig.savefig(plots_dir / 'similarity_histograms.png', bbox_inches='tight')
    plt.close(fig)

    return plots_paths_dict

def generate_dataset_html_report(stats1, stats2, output_path, plots_path, end_position, plot_type, results, failed_by_feature):
    """Generate HTML report comparing two dataset statistics.

    Generates plots with colored failure indicators (red #c62828 for Fail, orange #f57f17 for Warning),
    writes HTML report and duplicate sequences file.

    Args:
        stats1, stats2: Statistics objects for two datasets.
        output_path: Path to save HTML report (.html).
        plots_path: Directory to save plot images.
        end_position: Maximum position for per-position plots.
        plot_type: Plot type ('boxen' or 'violin').
        results: Optional DataFrame with pre-computed flags. Column 'Flag' used for summary statuses.
        failed_by_feature: Dict with failure info for plot shading:
            {
                'Sequence lengths': {'Pass'},
                'Per sequence GC content': {'Warning'},
                'Per sequence nucleotide content': {'A': 'Warning', 'G': 'Fail'},
                'Per sequence dinucleotide content': {'AA': 'Fail'},
                'Per position nucleotide content': {'A': {52: 'Warning'}, 'G': {66: 'Fail'}},
                'Per reverse position nucleotide content': {...},
                'Sequence Duplications within Labels': {'Pass'}
            }
    """
    plots_path.mkdir(parents=True, exist_ok=True)

    if results is not None:
        summary_statuses = results['Flag'].to_dict()
    else:
        summary_statuses = None

    # Extract percent remaining from results if available
    percent_remaining = None
    if results is not None and 'Sequence Duplications within Labels' in results.index:
        dup_info = results.loc['Sequence Duplications within Labels']
        if isinstance(dup_info, pd.Series) and 'Percent Remaining' in dup_info:
            percent_remaining = dup_info['Percent Remaining']
    # generate plots (with failure shading if failed_by_feature available)
    plots_paths = generate_dataset_plots(
        stats1, stats2, plots_path, end_position, plot_type,
        failed_by_feature=failed_by_feature,
        percent_remaining=percent_remaining
    )

    # find duplicate sequences between labels
    duplicate_seqs = list(set(stats1.sequences).intersection(stats2.sequences))
    # remove extension from output path, add '_duplicates.txt'
    duplicate_seqs_path = os.path.splitext(output_path)[0] + '_duplicates.txt'

    # Load the HTML template
    template = get_dataset_html_template(stats1, stats2, plots_paths, summary_statuses, duplicate_seqs, duplicate_seqs_file=duplicate_seqs_path)

    with open(output_path, 'w') as file:
        file.write(template)

    if len(duplicate_seqs) > 0:
        with open(duplicate_seqs_path, 'w') as f:
            for seq in sorted(duplicate_seqs):
                f.write(f"{seq}\n")
        logging.info(f"Duplicate sequences saved to {duplicate_seqs_path}")

def generate_json_report(stats_dict, output_path):
    logging.info(f"Generating JSON report: {output_path}")
    write_stats_json(stats_dict, output_path)

def generate_simple_report(results, output_path):

    logging.info(f"Generating simple report: {output_path}")

    if isinstance(results, dict):
        results = pd.DataFrame.from_dict(results, orient='index')
    results.to_csv(output_path)

def generate_dataset_plots(stats1, stats2, output_path, end_position, plot_type='boxen', failed_by_feature=None, percent_remaining=None):
    """Generate comparison plots between two datasets.

    Args:
        stats1, stats2: Statistics objects for two datasets.
        output_path: Directory to save plots.
        end_position: Maximum position for per-position plots.
        plot_type: Plot type ('boxen' or 'violin').
        failed_by_feature: Dict with failure info for shading (optional).
            When a feature has failures, the corresponding _with_flags plot
            is generated and returned; otherwise only the _no_flags plot is
            generated and returned.
        percent_remaining: Optional float with the percentage of sequences remaining after deduplication.

    Returns:
        Dictionary mapping plot names to file paths.
    """

    logging.info(f"Generating PNG plots at: {output_path}")

    plots_paths = {}

    bases_overlap = sorted(list(set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases'])))

    # Get failed nucleotides and dinucleotides for boxen plot outlines
    failed_nucleotides = failed_by_feature.get('Per sequence nucleotide content', {}) if failed_by_feature else None
    failed_dinucleotides = failed_by_feature.get('Per sequence dinucleotide content', {}) if failed_by_feature else None
    failed_pos_forward = failed_by_feature.get('Per position nucleotide content', {}) if failed_by_feature else None
    failed_pos_reverse = failed_by_feature.get('Per reverse position nucleotide content', {}) if failed_by_feature else None

    # Handle disjoint bases case - no common bases between datasets
    if not bases_overlap:
        
        # Create placeholder plots for nucleotide/dinucleotide/per-position content
        placeholder_nuc = output_path / 'per_sequence_nucleotide_content_with_flags.png'
        classes_plots.create_placeholder_plot(placeholder_nuc,
                                              'No common bases between datasets\nCannot compare nucleotide content')
        plots_paths['Per sequence nucleotide content'] = placeholder_nuc

        placeholder_dinuc = output_path / 'per_sequence_dinucleotide_content_with_flags.png'
        classes_plots.create_placeholder_plot(placeholder_dinuc,
                                              'No common bases between datasets\nCannot compare dinucleotide content')
        plots_paths['Per sequence dinucleotide content'] = placeholder_dinuc

        placeholder_pos_fwd = output_path / 'per_position_nucleotide_content_with_flags.png'
        classes_plots.create_placeholder_plot(placeholder_pos_fwd,
                                              'No common bases between datasets\nCannot compare per-position content')
        plots_paths['Per position nucleotide content'] = placeholder_pos_fwd

        placeholder_pos_rev = output_path / 'per_position_reversed_nucleotide_content_with_flags.png'
        classes_plots.create_placeholder_plot(placeholder_pos_rev,
                                              'No common bases between datasets\nCannot compare per-position content')
        plots_paths['Per position reversed nucleotide content'] = placeholder_pos_rev

    else:

        # Plot per sequence nucleotide content - create both versions
        # No-flags version
        fig = classes_plots.plot_nucleotides(
            stats1,
            stats2,
            nucleotides = bases_overlap,
            plot_type=plot_type,
            failed_nucleotides=None
        )
        no_flags_path = output_path / 'per_sequence_nucleotide_content_no_flags.png'
        fig.savefig(no_flags_path, bbox_inches='tight')
        plt.close(fig)
        if failed_nucleotides:
            fig = classes_plots.plot_nucleotides(
                stats1,
                stats2,
                nucleotides = bases_overlap,
                plot_type=plot_type,
                failed_nucleotides=failed_nucleotides
            )
            with_flags_path = output_path / 'per_sequence_nucleotide_content_with_flags.png'
            fig.savefig(with_flags_path, bbox_inches='tight')
            plt.close(fig)
            plots_paths['Per sequence nucleotide content'] = with_flags_path
        else:
            plots_paths['Per sequence nucleotide content'] = no_flags_path

        # Plot per sequence dinucleotide content - create both versions
        # No-flags version
        fig = classes_plots.plot_dinucleotides(
            stats1,
            stats2,
            nucleotides = bases_overlap,
            plot_type=plot_type,
            failed_dinucleotides=None
        )
        no_flags_path = output_path / 'per_sequence_dinucleotide_content_no_flags.png'
        fig.savefig(no_flags_path, bbox_inches='tight')
        plt.close(fig)
        if failed_dinucleotides:
            fig = classes_plots.plot_dinucleotides(
                stats1,
                stats2,
                nucleotides = bases_overlap,
                plot_type=plot_type,
                failed_dinucleotides=failed_dinucleotides
            )
            with_flags_path = output_path / 'per_sequence_dinucleotide_content_with_flags.png'
            fig.savefig(with_flags_path, bbox_inches='tight')
            plt.close(fig)
            plots_paths['Per sequence dinucleotide content'] = with_flags_path
        else:
            plots_paths['Per sequence dinucleotide content'] = no_flags_path

        # Plot per position nucleotide content (forward) - create both versions
        # No-flags version
        fig = classes_plots.plot_per_base_sequence_comparison(
            stats1,
            stats2,
            stats_name='Per position nucleotide content',
            nucleotides = bases_overlap,
            end_position=end_position,
            x_label='Position in sequence',
            failed_positions=None
        )
        no_flags_path = output_path / 'per_position_nucleotide_content_no_flags.png'
        fig.savefig(no_flags_path, bbox_inches='tight')
        plt.close(fig)
        if failed_pos_forward:
            fig = classes_plots.plot_per_base_sequence_comparison(
                stats1,
                stats2,
                stats_name='Per position nucleotide content',
                nucleotides = bases_overlap,
                end_position=end_position,
                x_label='Position in sequence',
                failed_positions=failed_pos_forward
            )
            with_flags_path = output_path / 'per_position_nucleotide_content_with_flags.png'
            fig.savefig(with_flags_path, bbox_inches='tight')
            plt.close(fig)
            plots_paths['Per position nucleotide content'] = with_flags_path
        else:
            plots_paths['Per position nucleotide content'] = no_flags_path

        # Plot per reversed position nucleotide content - create both versions
        # No-flags version
        fig = classes_plots.plot_per_base_sequence_comparison(
            stats1,
            stats2,
            stats_name='Per position reversed nucleotide content',
            nucleotides = bases_overlap,
            end_position=end_position,
            x_label='Position in reversed sequence',
            failed_positions=None
        )
        no_flags_path = output_path / 'per_position_reversed_nucleotide_content_no_flags.png'
        fig.savefig(no_flags_path, bbox_inches='tight')
        plt.close(fig)
        if failed_pos_reverse:
            fig = classes_plots.plot_per_base_sequence_comparison(
                stats1,
                stats2,
                stats_name='Per position reversed nucleotide content',
                nucleotides = bases_overlap,
                end_position=end_position,
                x_label='Position in reversed sequence',
                failed_positions=failed_pos_reverse
            )
            with_flags_path = output_path / 'per_position_reversed_nucleotide_content_with_flags.png'
            fig.savefig(with_flags_path, bbox_inches='tight')
            plt.close(fig)
            plots_paths['Per position reversed nucleotide content'] = with_flags_path
        else:
            plots_paths['Per position reversed nucleotide content'] = no_flags_path

    # Plot length distribution
    fig = classes_plots.plot_lengths(
        stats1,
        stats2,
        plot_type=plot_type
    )
    lengths_path = output_path / 'sequence_lengths.png'
    fig.savefig(lengths_path, bbox_inches='tight')
    plt.close(fig)
    plots_paths['Sequence lengths'] = lengths_path

    # Plot per sequence GC content
    fig = classes_plots.plot_gc_content(
        stats1,
        stats2,
        plot_type=plot_type
    )
    gc_path = output_path / 'per_sequence_gc_content.png'
    fig.savefig(gc_path, bbox_inches='tight')
    plt.close(fig)
    plots_paths['Per sequence GC content'] = gc_path

    # Plot sequence duplications within classes - only if some sequences were removed by deduplication
    if percent_remaining is not None and percent_remaining < 1.0:
        fig = classes_plots.plot_sequence_duplications_within_classes(
            stats1,
            stats2,
            percent_remaining_after_dedup=percent_remaining
        )
        dup_path = output_path / 'sequence_duplications_within_labels.png'
        fig.savefig(dup_path, bbox_inches='tight')
        plt.close(fig)
        plots_paths['Sequence Duplications within Labels'] = dup_path

    return plots_paths
