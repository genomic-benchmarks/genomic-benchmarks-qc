"""Writing the reports: which files each report type produces, and their content.

The entry points both commands call. Everything about *what* a report looks like
lives further down - the templates in `classes_html_report`/`split_html_report`,
the figures in `classes_plots`/`splits_plots` - and the names of the files are
defined in `genomic_benchmarks_qc.utils.naming`.

The figures are imported inside the two functions that draw them rather than at
the top of this module. seaborn brings scipy.stats and matplotlib with it, which
is a second of startup and 130 MB of memory before Typer has read an argument -
paid by `gb-qc --help`, and by every run that asks only for `simple` or `json`
reports and never draws anything. `test_startup.py` holds that line.
"""

import logging
from pathlib import Path

import pandas as pd

from genomic_benchmarks_qc.report.classes_html_report import get_dataset_html_template
from genomic_benchmarks_qc.report.per_position_payload import X_LABELS, build_payload, drawn_window
from genomic_benchmarks_qc.report.split_html_report import get_splits_html_template
from genomic_benchmarks_qc.utils.input_utils import write_stats_json
from genomic_benchmarks_qc.utils.naming import DUPLICATES_FILE


def save_plot(fig, path):
    """Write one figure to `path`, and return where it went.

    One call per file a figure becomes, which is two for a figure that flags:
    the callers below draw it, save it, mark the flags on it and save it again.

    Returns:
        `path`, so the caller can record it as the plot the page will show.
    """
    fig.savefig(path, bbox_inches='tight')
    return path


def validate_report_types(report_types, valid_types, command):
    """Reject a report type a command does not produce.

    Which types exist is not a property of the tool but of each command - the
    classes command writes per-class JSON and the splits command has nothing to
    put in one - so the set is declared next to the code that implements it and
    passed in here.

    Raises ValueError rather than skipping the unknown type. Skipping is what the
    splits command used to do with `json`: a full MMseqs2 search, a "successfully
    completed" line, exit 0, and no file anywhere.
    """
    unknown = [str(report_type) for report_type in report_types
               if report_type not in valid_types]
    if unknown:
        raise ValueError(
            f"{command} cannot produce report type(s): {', '.join(unknown)}. "
            f"Valid types are: {', '.join(valid_types)}."
        )


def generate_splits_html_report(basic_stats, threshold_stats, results_filt, output_path,
                                plots_dir, query_similarity_max, target_similarity_max,
                                leaked_hits=None):
    """
    Generate an HTML report visualising data leakage.

    `results_filt` holds the hits the page lists; `leaked_hits` is how many were
    at or above the similarity threshold, so the listing can say what it is not
    showing. It defaults to the number of rows given, for a caller that does not
    cap.
    """
    plots_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"Generating HTML report: {output_path}")

    plots_paths_dict = generate_split_plots(
        query_similarity_max, target_similarity_max, threshold_stats, plots_dir)

    template = get_splits_html_template(basic_stats, threshold_stats, results_filt,
                                        plots_paths_dict, leaked_hits=leaked_hits)
    with open(output_path, 'w') as file:
        file.write(template)

def generate_split_plots(query_similarity_max, target_similarity_max, threshold_stats, plots_dir):
    """Save the split figures and return {plot title: path} for the template."""
    # Deferred, for the reason in the module docstring.
    import matplotlib.pyplot as plt

    from genomic_benchmarks_qc.report import splits_plots

    plots_paths_dict = {}

    fig = splits_plots.plot_similarity_histograms(
        query_similarity_max, target_similarity_max, threshold_stats)
    plots_paths_dict['Similarity histograms'] = save_plot(
        fig, plots_dir / 'similarity_histograms.png')
    plt.close(fig)

    return plots_paths_dict

def generate_dataset_html_report(stats1, stats2, output_path, plots_path, plot_type, results,
                                 failed_by_feature):
    """Generate HTML report comparing two dataset statistics.

    Generates plots with colored failure indicators (red #c62828 for Fail, orange #f57f17
    for Warning), writes HTML report and duplicate sequences file.

    Args:
        stats1, stats2: Statistics objects for two datasets.
        output_path: Path to save HTML report (.html).
        plots_path: Directory to save plot images.
        plot_type: Plot type ('boxen' or 'violin').
        results: DataFrame with pre-computed flags, as produced by
            `flag_significant_differences`. Column 'Flag' is used for the summary
            statuses the template needs, so it is required, not optional.
        failed_by_feature: Dict with failure info for plot shading:
            {
                'Sequence lengths': {'Pass'},
                'Per sequence GC content': {'Warning'},
                'Per sequence nucleotide content': {'A': 'Warning', 'G': 'Fail'},
                'Per sequence dinucleotide content': {'AA': 'Fail'},
                'Per position nucleotide content': {'A': {52: 'Warning'}, 'G': {66: 'Fail'}},
                'Per position reversed nucleotide content': {...},
                'Sequence Duplications within Labels': {'Pass'}
            }
    """
    plots_path.mkdir(parents=True, exist_ok=True)

    # The per-position window comes from the two classes' sequence lengths, so it
    # is read off the statistics rather than passed in: a window that disagreed
    # with the statistics it is drawn from would misplace every flag. It and the
    # bases the two classes share are resolved once here and handed to the plots
    # and to the payload alike, so the PNGs in plots/ and the interactive
    # figures cannot end up drawn over different windows or different bases.
    end_position, compared = drawn_window(stats1, stats2)
    bases_overlap = sorted(set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases']))

    summary_statuses = results['Flag'].to_dict()

    # Extract percent remaining, which only the duplication check computes
    percent_remaining = None
    if 'Sequence Duplications within Labels' in results.index:
        dup_info = results.loc['Sequence Duplications within Labels']
        if isinstance(dup_info, pd.Series) and 'Percent Remaining' in dup_info:
            percent_remaining = dup_info['Percent Remaining']
    # generate plots (with failure shading if failed_by_feature available)
    plots_paths = generate_dataset_plots(
        stats1, stats2, plots_path, plot_type,
        failed_by_feature=failed_by_feature,
        percent_remaining=percent_remaining,
        end_position=end_position,
        bases_overlap=bases_overlap
    )

    # find duplicate sequences between labels
    duplicate_seqs = list(set(stats1.sequences).intersection(stats2.sequences))
    # the report lives in a directory of its own, so this is a plain sibling
    duplicate_seqs_path = Path(output_path).parent / DUPLICATES_FILE

    # The data behind the interactive per-position figures. Everything it needs
    # was computed for the flags and the plots; this only reshapes it. The window
    # is the compared one, or the reported one when nothing could be compared -
    # see `drawn_window`. Both directions come back None when the two classes
    # share no bases, which is the same condition that leaves those plots unmade.
    per_position_payloads = {
        direction: build_payload(stats1, stats2, bases_overlap, end_position,
                                 results, direction, compared)
        for direction in ('forward', 'reversed')
    }

    # Load the HTML template
    template = get_dataset_html_template(stats1, stats2, plots_paths, summary_statuses,
                                         duplicate_seqs,
                                         duplicate_seqs_file=duplicate_seqs_path,
                                         per_position_payloads=per_position_payloads)

    with open(output_path, 'w') as file:
        file.write(template)

    if len(duplicate_seqs) > 0:
        with open(duplicate_seqs_path, 'w') as f:
            for seq in sorted(duplicate_seqs):
                f.write(f"{seq}\n")
        logging.info(f"Duplicate sequences saved to {duplicate_seqs_path}")

def generate_json_report(stats_dict, output_path):
    """Dump one class's computed statistics to JSON."""
    logging.info(f"Generating JSON report: {output_path}")
    write_stats_json(stats_dict, output_path)

def generate_simple_report(results, output_path):
    """Write the per-check flags as a small CSV, one row per check."""

    logging.info(f"Generating simple report: {output_path}")

    if isinstance(results, dict):
        results = pd.DataFrame.from_dict(results, orient='index')
    # 'Percent Remaining' is computed only for plotting, exclude it from the CSV report
    results = results.drop(columns=['Percent Remaining'], errors='ignore')
    results.index.name = 'Check'
    results.to_csv(output_path)

def generate_dataset_plots(stats1, stats2, output_path, plot_type='boxen',
                           failed_by_feature=None, percent_remaining=None,
                           end_position=None, bases_overlap=None):
    """Generate comparison plots between two datasets.

    Args:
        stats1, stats2: Statistics objects for two datasets.
        output_path: Directory to save plots.
        plot_type: Plot type ('boxen' or 'violin').
        failed_by_feature: Dict with failure info for shading (optional).
        percent_remaining: Optional float with the percentage of sequences remaining
            after deduplication.
        end_position: Optional last per-position position to draw, 1-based and
            inclusive, as `drawn_window` resolves it. Read off the statistics
            when not given.
        bases_overlap: Optional sorted list of the bases the two classes share.
            Read off the statistics when not given.

    Returns:
        Dictionary mapping plot names to file paths, or to None for a figure
        that could not be drawn.
    """

    logging.info(f"Generating PNG plots at: {output_path}")

    # Deferred, for the reason in the module docstring.
    import matplotlib.pyplot as plt

    from genomic_benchmarks_qc.report import classes_plots

    # The per-position figures draw only what was compared, so the scored window
    # is the one they are plotted over; the wider window the checks themselves
    # run over is not drawn, except when nothing was compared and it is all
    # there is - see `drawn_window`. A caller that already resolved the window
    # and the shared bases hands them in, so the figures it saves and the flags
    # it draws over them come from one answer rather than two.
    if end_position is None:
        end_position, _ = drawn_window(stats1, stats2)
    if bases_overlap is None:
        bases_overlap = sorted(
            set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases'])
        )

    plots_paths = {}

    # Get failed nucleotides and dinucleotides for boxen plot outlines
    def failed_for(feature):
        return failed_by_feature.get(feature, {}) if failed_by_feature else None

    failed_nucleotides = failed_for('Per sequence nucleotide content')
    failed_dinucleotides = failed_for('Per sequence dinucleotide content')
    failed_pos_forward = failed_for('Per position nucleotide content')
    failed_pos_reverse = failed_for('Per position reversed nucleotide content')

    # Handle disjoint bases case - no common bases between datasets.
    # No plots are produced; the HTML report renders an explanatory message
    # in place of each affected plot (see get_dataset_html_template).
    if not bases_overlap:

        plots_paths['Per sequence nucleotide content'] = None
        plots_paths['Per sequence dinucleotide content'] = None
        plots_paths['Per position nucleotide content'] = None
        plots_paths['Per position reversed nucleotide content'] = None

    else:

        # Both files come from one figure: it is saved clean, the flags are
        # marked on it, and it is saved again. `reserve_flag_margin` keeps the
        # two crops identical, so the pair differs only by the marks.
        fig = classes_plots.plot_nucleotides(
            stats1,
            stats2,
            nucleotides = bases_overlap,
            plot_type=plot_type
        )
        plots_paths['Per sequence nucleotide content'] = save_plot(
            fig, output_path / 'per_sequence_nucleotide_content.png')
        if failed_nucleotides:
            classes_plots.mark_failed_nucleotides(fig, bases_overlap, failed_nucleotides)
            plots_paths['Per sequence nucleotide content'] = save_plot(
                fig, output_path / 'per_sequence_nucleotide_content_with_flags.png')
        plt.close(fig)

        fig = classes_plots.plot_dinucleotides(
            stats1,
            stats2,
            nucleotides = bases_overlap,
            plot_type=plot_type
        )
        plots_paths['Per sequence dinucleotide content'] = save_plot(
            fig, output_path / 'per_sequence_dinucleotide_content.png')
        if failed_dinucleotides:
            classes_plots.mark_failed_dinucleotides(fig, bases_overlap, failed_dinucleotides)
            plots_paths['Per sequence dinucleotide content'] = save_plot(
                fig, output_path / 'per_sequence_dinucleotide_content_with_flags.png')
        plt.close(fig)

        # The per-position figures, one per direction, drawn over the window
        # `drawn_window` resolves: what was compared, or - when nothing was - the
        # positions the checks are named for, every one of them Unknown. Only a
        # comparison with no per-position checks at all leaves them unmade, and
        # the HTML renders a message in their place the way it does for
        # disjoint bases.
        # The axis wording comes from X_LABELS, the same place the interactive
        # figure and its flag table take it from, so the PNG in plots/ and the
        # figure on the page cannot end up counting positions differently.
        for stats_name, x_label, file_stem, failed_positions in (
            ('Per position nucleotide content', X_LABELS['forward'],
             'per_position_nucleotide_content', failed_pos_forward),
            ('Per position reversed nucleotide content', X_LABELS['reversed'],
             'per_position_reversed_nucleotide_content', failed_pos_reverse),
        ):
            if end_position < 1:
                plots_paths[stats_name] = None
                continue

            fig = classes_plots.plot_per_base_sequence_comparison(
                stats1,
                stats2,
                stats_name=stats_name,
                nucleotides = bases_overlap,
                end_position=end_position,
                x_label=x_label
            )
            plots_paths[stats_name] = save_plot(
                fig, output_path / f'{file_stem}.png')

            if failed_positions:
                classes_plots.mark_failed_positions(
                    fig, bases_overlap, end_position, failed_positions)
                plots_paths[stats_name] = save_plot(
                    fig, output_path / f'{file_stem}_with_flags.png')
            plt.close(fig)

    # Plot length distribution
    fig = classes_plots.plot_lengths(
        stats1,
        stats2,
        plot_type=plot_type
    )
    plots_paths['Sequence lengths'] = save_plot(fig, output_path / 'sequence_lengths.png')
    plt.close(fig)

    # Plot per sequence GC content
    fig = classes_plots.plot_gc_content(
        stats1,
        stats2,
        plot_type=plot_type
    )
    plots_paths['Per sequence GC content'] = save_plot(
        fig, output_path / 'per_sequence_gc_content.png')
    plt.close(fig)

    # Plot sequence duplications within classes - only if some sequences were removed
    # by deduplication
    if percent_remaining is not None and percent_remaining < 1.0:
        fig = classes_plots.plot_sequence_duplications_within_classes(
            stats1,
            stats2,
            percent_remaining_after_dedup=percent_remaining
        )
        plots_paths['Sequence Duplications within Labels'] = save_plot(
            fig, output_path / 'sequence_duplications_within_labels.png')
        plt.close(fig)

    return plots_paths
