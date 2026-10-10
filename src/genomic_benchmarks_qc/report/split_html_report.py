"""The HTML page for the train-test split report.

Self-contained in the same way as the class report: embedded plots, inlined
styling, and the top alignments rendered inline as text.

The page is built from the same components as the class report - the cards, the
key/value tables, and the collapsible findings panel with a listing table inside
it (.qc-panel / .qc-listing in report_design.css) - and in the same way: chrome
around the checks, whose navigation entries and sections come from the registry
in `genomic_benchmarks_qc.checks.splits`. The one check there is Data Leakage.
"""

from datetime import datetime

from genomic_benchmarks_qc import __version__
from genomic_benchmarks_qc.checks import run_checks
from genomic_benchmarks_qc.checks import splits as split_checks
from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.sections import (
    render_nav,
    render_sections,
    section_stylesheets,
)
from genomic_benchmarks_qc.report.utils import (
    COMMON_CSS,
    LOGO_CSS,
    REPORT_HEADER_HTML,
    REPORT_LOGO_HTML,
    SIDEBAR_LINKS_HTML,
    TOOL_DESCRIPTION,
    TOOL_TAGLINE,
    escape_html_text,
    put_data,
    put_text,
    verdict_html,
)

HTML_TEMPLATE = assets.template('split_report_page.html')


def render_splits_html(context, tool_description=None):
    """Build the split report page, drawing every check's figures on the way.

    Args:
        context: The `SplitReportContext` every check's section is drawn from.
        tool_description: Short description of the tool for the header;
            `TOOL_DESCRIPTION` when not given.

    Returns:
        The page, as one string.
    """
    basic_stats, threshold_stats = context.basic_stats, context.threshold_stats
    # Looked up here, not at import, so what the page shows is the registry as it
    # stands when the page is built. The flags are the ones the simple report
    # writes, scored again: a check's score is cheap and has no side effects.
    checks = split_checks.SPLIT_CHECKS
    results, _ = run_checks(checks, threshold_stats)
    summary_statuses = {name: row['Flag'] for name, row in results.items()}

    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment; the sections' stylesheets are
    # layered on top - split_report.css for the leakage listing - the way
    # per_position_viewer.css is in the class report
    html_template = put_data(html_template, "{{common_css}}",
                             COMMON_CSS + section_stylesheets(checks) + LOGO_CSS)
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
    html_template = put_data(html_template, "{{report_verdict}}",
                             verdict_html(summary_statuses, [check.name for check in checks]))
    html_template = put_data(html_template, "{{sidebar_links}}", SIDEBAR_LINKS_HTML)
    html_template = put_data(
        html_template, "{{page_title}}",
        escape_html_text(f'{basic_stats["test_filename"]} vs '
                    f'{basic_stats["train_filename"]}: leakage - gb-qc'))

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

    # The checks go in last, once every placeholder of the page itself is filled,
    # so nothing a check wrote can be mistaken for one of them.
    sections, section_scripts = render_sections(checks, context, summary_statuses)

    # The shared behaviour, then the sections' own - the alignment toggles. At
    # the end of the body so it runs against a complete page. This page puts the
    # shared script first, where the class report puts it last; each is what
    # its page has always done.
    html_template = put_data(html_template, "{{report_scripts}}",
                             assets.script('report_ui.js')
                             + ''.join('\n' + script for script in section_scripts))
    html_template = put_data(html_template, "{{check_nav}}", render_nav(checks, summary_statuses))
    return put_data(html_template, "{{check_sections}}", sections)
