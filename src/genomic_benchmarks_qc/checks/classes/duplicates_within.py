"""Sequence Duplications within Labels: how much of the data survives deduplication.

Decided by a rule, not a score, pooled over both classes: Pass only when every
sequence is unique, Warning while at least 98% are, Fail below that.
"""

import pandas as pd

from genomic_benchmarks_qc.checks import Check, CheckResult
from genomic_benchmarks_qc.report.sections import Section, SectionContent
from genomic_benchmarks_qc.report.utils import encode_image_to_base64, save_figure

NAME = 'Sequence Duplications within Labels'


def score(stats1, stats2) -> CheckResult:
    """Flag the combined share of sequences left after deduplication.

    The row also carries that share as 'Percent Remaining', which the report's
    figure is drawn from and `gb-qc-report.csv` leaves out.
    """
    total_sequences = 0
    unique_sequences = 0
    for stats in [stats1, stats2]:
        total_sequences += stats.stats['Number of sequences']
        unique_sequences += stats.stats['Number of sequences left after deduplication']

    percent_remaining = unique_sequences / total_sequences if total_sequences > 0 else 1.0

    if percent_remaining < 0.98:
        flag = 'Fail'
    elif percent_remaining < 1.0:
        flag = 'Warning'
    else:
        flag = 'Pass'

    return CheckResult({NAME: {'Flag': flag, 'Percent Remaining': percent_remaining}})


def _percent_remaining(results):
    """The share left after deduplication, as this check's row recorded it."""
    if NAME in results.index:
        row = results.loc[NAME]
        if isinstance(row, pd.Series) and 'Percent Remaining' in row:
            return row['Percent Remaining']
    return None


def render(context) -> SectionContent:
    """How many copies the duplicated sequences have, when deduplication removed any."""
    # Deferred, like every figure: see `report_generator`.
    from genomic_benchmarks_qc.report import classes_plots

    # Drawn only when deduplication removed something, which is exactly when the
    # flag is not a Pass - so the figure and the flag cannot disagree.
    plot = None
    percent_remaining = _percent_remaining(context.results)
    if percent_remaining is not None and percent_remaining < 1.0:
        fig = classes_plots.plot_sequence_duplications_within_classes(
            context.stats1, context.stats2, percent_remaining_after_dedup=percent_remaining)
        plot = save_figure(fig, context.plots_dir, 'sequence_duplications_within_labels')

    if context.summary_statuses[NAME].lower() in ('pass', 'ok', 'good', 'success'):
        body = """
        <p>No duplicate sequences were found in either class.</p>
        """
    else:
        body = (f'<img src="data:image/png;base64, {encode_image_to_base64(plot)}" '
                'alt="Sequence Duplications within Labels Plot" class="plot-wide">')
    return SectionContent({'{{duplications}}': body})


CHECK = Check(
    name=NAME,
    score=score,
    section=Section(
        title='Sequence Duplications within Labels',
        anchor='sequence-duplications-within-labels',
        explanation_id='within-dup-explanation',
        template='check_duplicates_within.html',
        render=render,
        docs=(('checks', 'sequence-duplications-within-labels', 'What to do about it'),),
    ),
)
