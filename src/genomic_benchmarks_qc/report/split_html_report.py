"""The HTML template for the train-test split report.

Self-contained in the same way as the class report: embedded plots, inlined
styling, and the top alignments rendered inline as text.

The page is built from the same components as the class report - the cards, the
key/value tables, and the collapsible findings panel with a listing table inside
it (.qc-panel / .qc-listing in report_design.css). The high-similarity pairs are
evidence for the Data Leakage check in the way the flagged positions are
evidence for a per-position check, so they sit in that check's card in the same
kind of panel, and not in a card of their own with no flag and no nav entry.
"""

import logging
from datetime import datetime

from genomic_benchmarks_qc import __version__
from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.alignment_rendering import (
    build_alignment_string,
    has_reversed_coordinates,
)
from genomic_benchmarks_qc.report.utils import (
    COMMON_CSS,
    LOGO_CSS,
    REPORT_HEADER_HTML,
    REPORT_LOGO_HTML,
    SIDEBAR_LINKS_HTML,
    TOOL_DESCRIPTION,
    TOOL_TAGLINE,
    docs_link,
    encode_image_to_base64,
    escape_html_text,
    icon_html,
    put_data,
    put_text,
    verdict_html,
)
from genomic_benchmarks_qc.utils.split_stats import flag_split_data_leakage

logger = logging.getLogger(__name__)

# Rows of the alignment listing that the page carries. The listing is evidence a
# reader spot-checks, not a data file - every hit above the threshold is exported
# to mmseqs/mmseqs2_search_result.tsv - and each row here embeds two coloured
# sequences, which is what the page weighs.
ROW_CAP = 100

HTML_TEMPLATE = assets.template('split_report_page.html')

RESULTS_TABLE = assets.template('split_results_table.html')


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


def get_splits_html_template(basic_stats, threshold_stats, results_filt, plots_paths_dict,
                             tool_description=None, leaked_hits=None):
    """Build Train-Test Split Check HTML using shared helpers.

    Args:
        results_filt: The hits to list, already capped by the caller or not.
        leaked_hits: How many hits were at or above the similarity threshold, so
            the panel can say what it is not showing. Not the number of
            alignments the search found, which is far larger and is not what the
            panel counts. Defaults to the number of rows given.
    """

    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment; split_report.css is layered on top,
    # the way per_position_viewer.css is in the class report
    html_template = put_data(html_template, "{{common_css}}",
                             COMMON_CSS + assets.stylesheet('split_report.css') + LOGO_CSS)
    # The shared behaviour, plus this report's own alignment toggles. At the end
    # of the body so it runs against a complete page.
    html_template = put_data(html_template, "{{report_scripts}}",
                             assets.script('report_ui.js', 'split_report.js'))
    html_template = put_data(html_template, "{{report_header}}", REPORT_HEADER_HTML)

    # The nav's logo is in the template; this is the header's. Both read the one
    # copy of the picture that LOGO_CSS put in the stylesheet above.
    html_template = put_data(html_template, "{{logo}}", REPORT_LOGO_HTML)

    # header info
    if tool_description is None:
        tool_description = TOOL_DESCRIPTION
    generated_on = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_template = put_data(html_template, "{{tool_tagline}}", TOOL_TAGLINE)
    html_template = put_data(html_template, "{{tool_description}}", tool_description)
    html_template = put_data(html_template, "{{generated_on}}", generated_on)
    html_template = put_data(html_template, "{{version}}", __version__)
    html_template = put_data(html_template, "{{row_cap}}", str(ROW_CAP))

    # The subject line here is the split itself: which file was searched against
    # which, in that order, because the flag is about the test side. Only names
    # are available - evaluate_splits keeps Path(f).name, not the path - so there
    # is nothing further to add after them.
    test_name = escape_html_text(str(basic_stats['test_filename']))
    train_name = escape_html_text(str(basic_stats['train_filename']))
    html_template = put_data(
        html_template, "{{report_subject}}",
        f'Searching test set <strong>{test_name}</strong> against train set '
        f'<strong>{train_name}</strong>')
    leakage_flags = {'Data Leakage': flag_split_data_leakage(
        threshold_stats['perc_queries_above_thr'])}
    html_template = put_data(html_template, "{{report_verdict}}",
                             verdict_html(leakage_flags, ('Data Leakage',)))
    html_template = put_data(html_template, "{{sidebar_links}}", SIDEBAR_LINKS_HTML)
    html_template = put_data(
        html_template, "{{page_title}}",
        escape_html_text(f'{basic_stats["test_filename"]} vs '
                    f'{basic_stats["train_filename"]}: leakage - gb-qc'))
    html_template = put_data(html_template, "{{link_leakage}}",
                             docs_link('leakage', text='What to do about it'))

    html_template = put_text(html_template, "{{train_filename}}", basic_stats['train_filename'])
    html_template = put_text(html_template, "{{test_filename}}", basic_stats['test_filename'])
    html_template = put_data(html_template, "{{number_of_sequences_train}}", str(basic_stats['number_of_sequences_train']))
    html_template = put_data(html_template, "{{number_of_sequences_test}}", str(basic_stats['number_of_sequences_test']))
    html_template = put_data(html_template, "{{min_length_train}}", str(basic_stats['min_length_train']))
    html_template = put_data(html_template, "{{mean_length_train}}", f"{basic_stats['mean_length_train']:.2f}")
    html_template = put_data(html_template, "{{max_length_train}}", str(basic_stats['max_length_train']))
    html_template = put_data(html_template, "{{min_length_test}}", str(basic_stats['min_length_test']))
    html_template = put_data(html_template, "{{mean_length_test}}", f"{basic_stats['mean_length_test']:.2f}")
    html_template = put_data(html_template, "{{max_length_test}}", str(basic_stats['max_length_test']))
    html_template = put_data(html_template, "{{perc_queries_above_thr}}", f"{threshold_stats['perc_queries_above_thr']:.2f}")
    html_template = put_data(html_template, "{{num_queries_above_thr}}", f"{threshold_stats['num_queries_above_thr']}")
    html_template = put_data(html_template, "{{perc_targets_above_thr}}", f"{threshold_stats['perc_targets_above_thr']:.2f}")
    html_template = put_data(html_template, "{{num_targets_above_thr}}", f"{threshold_stats['num_targets_above_thr']}")

    # insert plots as base64-encoded images
    html_template = put_data(html_template, "{{histogram_similarity_base64}}",
                             encode_image_to_base64(plots_paths_dict['Similarity histograms']))

    html_template = put_data(
        html_template,
        "{{icon_data_leakage}}",
        icon_html(leakage_flags, "Data Leakage")
    )

    results_display = results_filt.head(ROW_CAP).copy()  # the page lists at most this many hits
    total = len(results_filt) if leaked_hits is None else leaked_hits
    html_template = put_data(html_template, "{{alignments_count}}",
                             alignments_count_text(total, len(results_display)))

    if len(results_display) == 0:
        # An empty panel says the same thing the flagged-position panel does when
        # a check found nothing: the listing is empty because there is nothing in
        # it, not because it failed to build.
        return put_data(
            html_template,
            "{{results_body}}",
            '<p class="qc-empty">No test sequence aligns to a train sequence at or above '
            'the similarity threshold.</p>'
        )

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

    results_table = put_data(RESULTS_TABLE, "{{results_rows}}", "\n".join(rows))
    return put_data(html_template, "{{results_body}}", results_table)
