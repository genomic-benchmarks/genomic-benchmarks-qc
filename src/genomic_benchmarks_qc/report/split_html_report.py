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

import html
import logging
from datetime import datetime

from genomic_benchmarks_qc import __version__
from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.alignment_rendering import build_alignment_string
from genomic_benchmarks_qc.report.utils import (
    COMMON_CSS,
    LOGO_BASE64,
    REPORT_HEADER_HTML,
    SIDEBAR_LINKS_HTML,
    TOOL_DESCRIPTION,
    TOOL_TAGLINE,
    docs_link,
    encode_image_to_base64,
    icon_html,
    put_data,
    put_text,
    verdict_html,
)
from genomic_benchmarks_qc.utils.split_stats import flag_split_data_leakage

# Rows of the alignment listing that the page carries. The listing is evidence a
# reader spot-checks, not a data file - every hit is exported to
# mmseqs/mmseqs2_search_result.tsv - and each row here embeds two coloured
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


def get_splits_html_template(basic_stats, threshold_stats, results_filt, plots_paths_dict,
                             tool_description=None, total_hits=None):
    """Build Train-Test Split Check HTML using shared helpers.

    Args:
        results_filt: The hits to list, already capped by the caller or not.
        total_hits: How many hits there were before any cap, so the panel can say
            what it is not showing. Defaults to the number of rows given.
    """

    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment; split_report.css is layered on top,
    # the way per_position_viewer.css is in the class report
    html_template = put_data(html_template, "{{common_css}}",
                             COMMON_CSS + assets.stylesheet('split_report.css'))
    # The shared behaviour, plus this report's own alignment toggles. At the end
    # of the body so it runs against a complete page.
    html_template = put_data(html_template, "{{report_scripts}}",
                             assets.script('report_ui.js', 'split_report.js'))
    html_template = put_data(html_template, "{{report_header}}", REPORT_HEADER_HTML)

    # insert logo
    html_template = put_data(html_template, "{{logo_base64}}", LOGO_BASE64)

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
    test_name = html.escape(str(basic_stats['test_filename']))
    train_name = html.escape(str(basic_stats['train_filename']))
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
        html.escape(f'{basic_stats["test_filename"]} vs '
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
    total = len(results_filt) if total_hits is None else total_hits
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
            if not alignment_error_logged:
                alignment_error_logged = True
                logging.warning(
                    "Alignment visualisation failed for query=%s target=%s: %s: %s. "
                    "This is often caused by a conda-installed mmseqs2 build. See the "
                    "mmseqs2 installation notes in README.md — installing precompiled "
                    "binaries is recommended over compiling from source. "
                    "Further alignment errors in this report are not logged.",
                    row['query'], row['target'], type(e).__name__, e,
                )
            alignment_str_color = (
                '<span class="aln-error">ALIGNMENT VISUALISATION ERROR</span><br>'
                f'{html.escape(type(e).__name__)}: {html.escape(str(e))}<br>'
                '<span class="aln-error-detail">This can happen with a '
                'conda-installed mmseqs2 — see the installation notes in README.md.</span>'
            )
        rows.append(f"""
    <tr>
        <td class="qc-mono">{html.escape(str(row['query']))}</td>
        <td class="qc-mono">{html.escape(str(row['target']))}</td>
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
