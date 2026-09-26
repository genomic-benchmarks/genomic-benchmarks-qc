"""Data Leakage: how much of the test set has a near-identical match in training.

The one check `evaluate-splits` runs. Every test sequence is searched against
the train set with MMseqs2 and scored on its best hit; a test sequence at or
above the similarity threshold is one the model has already seen, so the
accuracy it earns there is not evidence. Decided by a rule on the share of test
sequences above the threshold: Pass at none, Warning below 2%, Fail from 2% on.

Its section lists the leaked alignments as evidence, the way a per-position
check lists its flagged positions: in a collapsible panel inside the check's own
section, not in a section of their own with no flag and no navigation entry.
"""

import logging

from genomic_benchmarks_qc.checks import Check, CheckResult
from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.alignment_rendering import (
    build_alignment_string,
    has_reversed_coordinates,
)
from genomic_benchmarks_qc.report.sections import Section, SectionContent
from genomic_benchmarks_qc.report.utils import (
    encode_image_to_base64,
    escape_html_text,
    put_data,
    save_figure,
)

logger = logging.getLogger(__name__)

NAME = 'Data Leakage'

# Rows of the alignment listing that the page carries. The listing is evidence a
# reader spot-checks, not a data file - every hit above the threshold is exported
# to mmseqs/mmseqs2_search_result.tsv - and each row here embeds two coloured
# sequences, which is what the page weighs.
ROW_CAP = 100

RESULTS_TABLE = assets.template('split_results_table.html')


def flag_split_data_leakage(perc_queries_above_thr, fail_threshold=2.0):
    """Grade a split by the share of test sequences with a near-identical match.

    Any leakage at all is worth knowing about, so a non-zero share is a Warning
    and only `fail_threshold` percent or more is a Fail.
    """
    if perc_queries_above_thr >= fail_threshold:
        return "Fail"
    if perc_queries_above_thr > 0:
        return "Warning"
    return "Pass"


def score(threshold_stats) -> CheckResult:
    """Flag the share of test sequences with a near-identical train sequence.

    The row carries both shares, as the percentages `gb-qc-report.csv` shows:
    how much of the test set is compromised, and how much of the training set is
    responsible.
    """
    return CheckResult({NAME: {
        "Flag": flag_split_data_leakage(threshold_stats['perc_queries_above_thr']),
        "Percentage of leaked queries": f"{threshold_stats['perc_queries_above_thr']:.2f}%",
        "Percentage of leaked targets": f"{threshold_stats['perc_targets_above_thr']:.2f}%",
    }})


def alignments_count_text(total, shown):
    """The panel's summary line: how many pairs leaked, and how many are listed.

    The same job the flagged-position panel's count does in the class report -
    say what is inside before it is opened - and it has to distinguish a clean
    split from a listing that was capped.
    """
    if total == 0:
        return 'No high-similarity alignments'
    plural = '' if total == 1 else 's'
    capped = f' (first {shown} shown)' if shown < total else ''
    return f'{total} high-similarity alignment{plural}{capped}'


def alignment_error_html(error, reversed_coords):
    """The cell shown in place of an alignment that could not be drawn.

    Two cases, because only one of them has a known cause. Coordinates running
    backwards come from the MMseqs2 build rather than from the data, and the
    row's scores are still right, so the cell says so and names what avoids it.
    Anything else is the validator refusing to draw a hit that disagrees with
    the sequences it was matched back to, and the honest thing there is to say
    only that.
    """
    if reversed_coords:
        detail = (
            'MMseqs2 reported this hit on the reverse strand, which the '
            'forward-strand-only search did not ask for, so its alignment cannot '
            'be drawn. The scores in this row are unaffected. Seen with '
            'conda/bioconda MMseqs2 above one thread \u2014 re-running with '
            '--threads 1, or installing the upstream precompiled release, avoids it.'
        )
    else:
        detail = (
            'The hit and the sequences it was matched back to disagree, so drawing '
            'it would produce a confident-looking but wrong alignment.'
        )
    return (
        '<span class="aln-error">ALIGNMENT VISUALISATION ERROR</span><br>'
        f'{escape_html_text(type(error).__name__)}: {escape_html_text(str(error))}<br>'
        f'<span class="aln-error-detail">{detail}</span>'
    )


def _listing_html(results_filt):
    """The panel's listing: one row per leaked alignment, each with its alignment."""
    results_display = results_filt.head(ROW_CAP).copy()  # the page lists at most this many hits
    if len(results_display) == 0:
        # An empty panel says the same thing the flagged-position panel does when
        # a check found nothing: the listing is empty because there is nothing in
        # it, not because it failed to build.
        return ('<p class="qc-empty">No test sequence aligns to a train sequence at or above '
                'the similarity threshold.</p>')

    rows = []
    alignment_error_logged = False

    for i, row in results_display.iterrows():
        try:
            alignment_str_color = build_alignment_string(row, color=True)
        except Exception as e:
            # A row whose coordinates run backwards has a known cause, already
            # counted and explained once by `log_reversed_hit_warning`. Only the
            # failures that do not have one are logged here, so a build emitting
            # hundreds of backwards rows does not say it twice.
            reversed_coords = has_reversed_coordinates(row)
            if not reversed_coords and not alignment_error_logged:
                alignment_error_logged = True
                logger.warning(
                    "Alignment visualisation failed for query=%s target=%s: %s: %s. "
                    "Further alignment errors in this report are not logged.",
                    row['query'], row['target'], type(e).__name__, e,
                )
            alignment_str_color = alignment_error_html(e, reversed_coords)
        rows.append(f"""
    <tr>
        <td class="qc-mono">{escape_html_text(str(row['query']))}</td>
        <td class="qc-mono">{escape_html_text(str(row['target']))}</td>
        <td class="qc-num">{row['qcov']:.2f}</td>
        <td class="qc-num">{row['tcov']:.2f}</td>
        <td class="qc-num">{row['pident']:.1f}</td>
        <td class="qc-key">{row['min_cov*pident']:.2f}</td>
        <td class="qc-num">{row['evalue']:.2e}</td>
        <td>
            <button type="button" class="qc-btn" onclick="toggleAlignment('aln-{i}', this)">
                Show
            </button>
        </td>
    </tr>

    <tr id="aln-{i}" style="display:none;">
        <td class="aln-cell" colspan="8">
            <pre class="alignment-block">{alignment_str_color}</pre>
        </td>
    </tr>
    """)

    return put_data(RESULTS_TABLE, "{{results_rows}}", "\n".join(rows))


def render(context) -> SectionContent:
    """The leakage figures, the similarity histogram, and the listing of what leaked."""
    # Deferred, like every figure: see `report_generator`.
    from genomic_benchmarks_qc.report import splits_plots

    threshold_stats = context.threshold_stats
    fig = splits_plots.plot_similarity_histograms(
        context.query_similarity_max, context.target_similarity_max, threshold_stats)
    plot = save_figure(fig, context.plots_dir, 'similarity_histograms')

    results_filt = context.results_filt
    shown = min(len(results_filt), ROW_CAP)
    total = len(results_filt) if context.leaked_hits is None else context.leaked_hits

    return SectionContent({
        '{{row_cap}}': str(ROW_CAP),
        '{{test_filename}}': escape_html_text(context.basic_stats['test_filename']),
        '{{train_filename}}': escape_html_text(context.basic_stats['train_filename']),
        '{{perc_queries_above_thr}}': f"{threshold_stats['perc_queries_above_thr']:.2f}",
        '{{num_queries_above_thr}}': f"{threshold_stats['num_queries_above_thr']}",
        '{{perc_targets_above_thr}}': f"{threshold_stats['perc_targets_above_thr']:.2f}",
        '{{num_targets_above_thr}}': f"{threshold_stats['num_targets_above_thr']}",
        '{{histogram_similarity_base64}}': encode_image_to_base64(plot),
        '{{alignments_count}}': alignments_count_text(total, shown),
        '{{results_body}}': _listing_html(results_filt),
    }, (assets.script('split_report.js'),))


CHECK = Check(
    name=NAME,
    score=score,
    section=Section(
        title='Data Leakage',
        anchor='similarity-section',
        explanation_id='data-leakage-explanation',
        template='check_data_leakage.html',
        render=render,
        docs=(('leakage', None, 'What to do about it'),),
        css_class='table-section',
        stylesheets=('split_report.css',),
    ),
)
