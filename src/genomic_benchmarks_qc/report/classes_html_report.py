"""The HTML page for the class comparison report.

Builds a single self-contained page: the plots are embedded as base64 data URIs
and the styling is inlined, so a report can be moved, zipped or served from
anywhere without losing anything.

The page is chrome around the checks. Its navigation entries and sections come
from the registry in `genomic_benchmarks_qc.checks.classes`, one per check, in
registry order - see `report.sections` for the shape they share.
"""

from datetime import datetime

from genomic_benchmarks_qc import __version__
from genomic_benchmarks_qc.checks import classes as class_checks
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
from genomic_benchmarks_qc.utils.testing import MIN_SEQUENCES_PER_CLASS

HTML_TEMPLATE = assets.template('classes_report.html')

def generate_not_scored_html(stats1, stats2, summary_statuses):
    """Explain, at the top of the report, any check that was not scored.

    A grey "?" in the sidebar is easy to read as a pass, and the difference is
    the whole point: Unknown means the comparison was not made, because there
    were not enough sequences behind it for its result to mean anything. Saying
    so once, in full, is what lets a reader trust the flags that are there.

    Returns an empty string when every check was scored, so the note costs
    nothing on a normally sized dataset.
    """
    if not summary_statuses:
        return ''

    # Only the headline checks; the per-base and per-position detail rows carry
    # their own Unknowns and would bury the message. Each check's floor says
    # which of the two reasons below is its own.
    unscored = [check for check in class_checks.CLASS_CHECKS
                if str(summary_statuses.get(check.name, '')).strip().lower() == 'unknown']
    if not unscored:
        return ''

    positional = [check.name for check in unscored if check.floor == 'per_position']
    per_sequence = [check.name for check in unscored if check.floor != 'per_position']

    def names(items):
        return ', '.join(f'<em>{name}</em>' for name in items)

    reasons = []
    if per_sequence:
        smaller = min(stats1.stats['Number of sequences'], stats2.stats['Number of sequences'])
        reasons.append(
            f'{names(per_sequence)} &mdash; the smaller class holds {smaller:,} sequences, and below '
            f'{MIN_SEQUENCES_PER_CLASS} these checks report a difference on sampling noise alone too '
            'often for a flag to be informative.'
        )
    if positional:
        reasons.append(
            f'{names(positional)} &mdash; every position is compared on the sequences long enough to '
            f'reach it, and no scored position has {MIN_SEQUENCES_PER_CLASS} sequences in both '
            'classes.'
        )

    items = ''.join(f'<li>{reason}</li>' for reason in reasons)
    return (
        '<div class="not-scored-note">'
        f'<strong>{len(unscored)} check(s) were not scored.</strong> They are marked '
        '<span class="status-icon status-unknown">?</span> rather than passed, because not enough '
        'data stands behind them to tell a real difference from sampling noise &mdash; which is not '
        'the same as having found no difference.'
        f'<ul>{items}</ul>'
        'The plots and the descriptive statistics below are computed from all the data and are '
        'unaffected, so the distributions can still be compared by eye.'
        '</div>'
    )


def render_dataset_html(context, tool_description=None):
    """Build the class comparison page, drawing every check's figures on the way.

    Args:
        context: The `ClassReportContext` every check's section is drawn from.
        tool_description: Short description of the tool for the header;
            `TOOL_DESCRIPTION` when not given.

    Returns:
        The page, as one string.
    """
    stats1, stats2 = context.stats1, context.stats2
    summary_statuses = context.summary_statuses
    # Looked up here, not at import, so what the page shows is the registry as it
    # stands when the page is built.
    checks = class_checks.CLASS_CHECKS

    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment
    html_template = put_data(html_template, "{{common_css}}",
                             COMMON_CSS + section_stylesheets(checks) + LOGO_CSS)
    html_template = put_data(html_template, "{{report_header}}", REPORT_HEADER_HTML)

    # The nav's logo is in the template; this is the header's. Both read the one
    # copy of the picture that LOGO_CSS put in the stylesheet above.
    html_template = put_data(html_template, "{{logo}}", REPORT_LOGO_HTML)

    # populate header placeholders: tool description, generated timestamp and input paths
    # Provide sensible defaults when values are not supplied
    if tool_description is None:
        tool_description = TOOL_DESCRIPTION
    if stats1.filepath == stats2.filepath:
        input_paths = stats1.filepath
    else:
        input_paths = f"{stats1.filepath}, {stats2.filepath}"
    generated_on = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_template = put_data(html_template, "{{tool_tagline}}", TOOL_TAGLINE)
    html_template = put_data(html_template, "{{tool_description}}", tool_description)
    html_template = put_data(html_template, "{{generated_on}}", generated_on)
    html_template = put_data(html_template, "{{version}}", __version__)

    # What this particular report looked at, and how it came out. The subject
    # line is the one thing a reader handed the file cannot work out for
    # themselves, and the verdict saves them counting circles in the navigation.
    label1, label2 = context.label1, context.label2
    subject = (f'Comparing label <strong>{escape_html_text(label1)}</strong> with label '
               f'<strong>{escape_html_text(label2)}</strong> in '
               f'{escape_html_text(str(input_paths))}')
    html_template = put_data(html_template, "{{report_subject}}", subject)
    html_template = put_data(html_template, "{{report_verdict}}",
                             verdict_html(summary_statuses, [check.name for check in checks]))
    html_template = put_data(html_template, "{{sidebar_links}}", SIDEBAR_LINKS_HTML)
    # Long enough to identify the report in a tab or a bookmark, ordered so that
    # a truncated tab keeps the part that distinguishes it from its neighbours.
    # The sequence column is last and still has to be there: a dataset with
    # several of them produces one report per column, and paired-sequences
    # produces three whose file and labels are identical.
    seq_col = str(stats1.seq_column) if stats1.seq_column is not None else None
    column = f' ({seq_col})' if seq_col else ''
    html_template = put_data(
        html_template, "{{page_title}}",
        escape_html_text(f'{stats1.filename}: {label1} vs {label2}{column} - gb-qc'))

    html_template = put_text(html_template, "{{filename1}}", stats1.filename)
    html_template = put_text(html_template, "{{filename2}}", stats2.filename)
    html_template = put_text(html_template, "{{label1}}", label1)
    html_template = put_text(html_template, "{{label2}}", label2)
    html_template = put_text(html_template, "{{seq_col1}}", str(stats1.seq_column) if stats1.seq_column is not None else "N/A")
    html_template = put_text(html_template, "{{seq_col2}}", str(stats2.seq_column) if stats2.seq_column is not None else "N/A")
    html_template = put_data(html_template, "{{number_of_sequences1}}", str(stats1.stats['Number of sequences']))
    html_template = put_data(html_template, "{{number_of_sequences2}}", str(stats2.stats['Number of sequences']))
    html_template = put_data(html_template, "{{dedup_sequences1}}", str(stats1.stats['Number of sequences left after deduplication']))
    html_template = put_data(html_template, "{{dedup_sequences2}}", str(stats2.stats['Number of sequences left after deduplication']))
    html_template = put_data(html_template, "{{min_length1}}", str(int(stats1.stats['Sequence lengths']['Sequence lengths'].min())))
    html_template = put_data(html_template, "{{min_length2}}", str(int(stats2.stats['Sequence lengths']['Sequence lengths'].min())))
    html_template = put_data(html_template, "{{mean_length1}}", f"{stats1.stats['Sequence lengths']['Sequence lengths'].mean():.2f}")
    html_template = put_data(html_template, "{{mean_length2}}", f"{stats2.stats['Sequence lengths']['Sequence lengths'].mean():.2f}")
    html_template = put_data(html_template, "{{max_length1}}", str(int(stats1.stats['Sequence lengths']['Sequence lengths'].max())))
    html_template = put_data(html_template, "{{max_length2}}", str(int(stats2.stats['Sequence lengths']['Sequence lengths'].max())))
    html_template = put_data(html_template, "{{number_of_bases1}}", str(stats1.stats['Number of bases']))
    html_template = put_data(html_template, "{{number_of_bases2}}", str(stats2.stats['Number of bases']))
    html_template = put_data(html_template, "{{gc_content1}}", f"{(stats1.stats['%GC content']*100):.2f}")
    html_template = put_data(html_template, "{{gc_content2}}", f"{(stats2.stats['%GC content']*100):.2f}")
    html_template = put_data(html_template, "{{not_scored_note}}",
                             generate_not_scored_html(stats1, stats2, summary_statuses))

    # The checks go in last, once every placeholder of the page itself is filled,
    # so nothing a check wrote can be mistaken for one of them. Each section
    # draws its own figures as it is rendered.
    sections, section_scripts = render_sections(checks, context, summary_statuses)

    # The behaviour goes in at the end of the body so it runs against a complete
    # page: the per-position viewer measures the width of the box it is drawn
    # into, which is only known once the layout exists. On this page the
    # sections' scripts come before the shared one, and a script no section asked
    # for - the viewer's, when there is nothing to draw - is left out entirely.
    html_template = put_data(html_template, "{{report_scripts}}",
                             ''.join(section_scripts) + assets.script('report_ui.js'))
    html_template = put_data(html_template, "{{check_nav}}", render_nav(checks, summary_statuses))
    return put_data(html_template, "{{check_sections}}", sections)
