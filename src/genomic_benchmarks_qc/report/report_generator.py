"""Writing the reports: which files each report type produces, and their content.

The entry points both commands call. Everything about *what* a report looks like
lives further down - the pages in `classes_html_report`/`split_html_report`, and
each check's own section, figures included, in its module under
`genomic_benchmarks_qc.checks` - and the names of the files are defined in
`genomic_benchmarks_qc.utils.naming`.

The figures are imported inside the functions that draw them - each check's
`render` - rather than at the top of any module on the way there. seaborn brings
scipy.stats and matplotlib with it, which is a second of startup and 130 MB of
memory before Typer has read an argument - paid by `gb-qc --help`, and by every
run that asks only for `simple` or `json` reports and never draws anything.
`test_startup.py` holds that line.
"""

import logging
from pathlib import Path

import pandas as pd

from genomic_benchmarks_qc.report.classes_html_report import render_dataset_html
from genomic_benchmarks_qc.report.per_position_payload import drawn_window
from genomic_benchmarks_qc.report.sections import ClassReportContext
from genomic_benchmarks_qc.report.split_html_report import get_splits_html_template
from genomic_benchmarks_qc.report.utils import save_plot
from genomic_benchmarks_qc.utils.input_utils import write_stats_json

logger = logging.getLogger(__name__)


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

    logger.info(f"Generating HTML report: {output_path}")

    plots_paths_dict = generate_split_plots(
        query_similarity_max, target_similarity_max, threshold_stats, plots_dir)

    template = get_splits_html_template(basic_stats, threshold_stats, results_filt,
                                        plots_paths_dict, leaked_hits=leaked_hits)
    with open(output_path, 'w') as file:
        file.write(template)

def generate_split_plots(query_similarity_max, target_similarity_max, threshold_stats, plots_dir):
    """Save the split figures and return {plot title: SavedPlot} for the template."""
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

    Every check draws its own figures into `plots_path` as its section is
    rendered, marking what it flagged (red #c62828 for Fail, orange #f57f17 for
    Warning), and the check of duplicates between the classes writes them to a
    file beside the report.

    Args:
        stats1, stats2: Statistics objects for two datasets.
        output_path: Path to save HTML report (.html).
        plots_path: Directory to save plot images.
        plot_type: Plot type ('boxen' or 'violin').
        results: DataFrame with pre-computed flags, as produced by
            `flag_significant_differences`. Column 'Flag' is used for the summary
            statuses the template needs, so it is required, not optional.
        failed_by_feature: What each figure shades, as
            `flag_significant_differences` returns it:
            {
                'Per sequence nucleotide content': {'A': 'Warning', 'G': 'Fail'},
                'Per sequence dinucleotide content': {'AA': 'Fail'},
                'Per position nucleotide content': {'A': {52: 'Warning'}, 'G': {66: 'Fail'}},
                'Per position reversed nucleotide content': {},
            }
    """
    plots_path.mkdir(parents=True, exist_ok=True)
    logger.info(f"Generating PNG plots at: {plots_path}")

    # The per-position window comes from the two classes' sequence lengths, so it
    # is read off the statistics rather than passed in: a window that disagreed
    # with the statistics it is drawn from would misplace every flag. It and the
    # bases the two classes share are resolved once here and handed to every
    # check alike, so the PNGs in plots/ and the interactive figures cannot end
    # up drawn over different windows or different bases.
    drawn_end, compared = drawn_window(stats1, stats2)
    context = ClassReportContext(
        stats1=stats1, stats2=stats2,
        results=results,
        summary_statuses=results['Flag'].to_dict(),
        failed_by_feature=failed_by_feature,
        plots_dir=plots_path,
        plot_type=plot_type,
        report_dir=Path(output_path).parent,
        bases_overlap=sorted(set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases'])),
        drawn_end=drawn_end,
        compared=compared,
    )

    page = render_dataset_html(context)
    with open(output_path, 'w') as file:
        file.write(page)

def generate_json_report(stats_dict, output_path):
    """Dump one class's computed statistics to JSON."""
    logger.info(f"Generating JSON report: {output_path}")
    write_stats_json(stats_dict, output_path)

def generate_simple_report(results, output_path):
    """Write the per-check flags as a small CSV, one row per check."""

    logger.info(f"Generating simple report: {output_path}")

    if isinstance(results, dict):
        results = pd.DataFrame.from_dict(results, orient='index')
    # 'Percent Remaining' is computed only for plotting, exclude it from the CSV report
    results = results.drop(columns=['Percent Remaining'], errors='ignore')
    # The flag first, whatever the other columns are. A table built from a dict
    # takes its columns in the order they first appear, so left to that the
    # header would depend on which check happens to open the table.
    if 'Flag' in results.columns:
        results = results[['Flag', *(column for column in results.columns if column != 'Flag')]]
    results.index.name = 'Check'
    results.to_csv(output_path)
