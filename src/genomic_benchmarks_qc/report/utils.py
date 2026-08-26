"""Shared pieces of the HTML reports: templating, colors, embedded images."""

import base64
import html
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from genomic_benchmarks_qc.report import assets

# Re-exported: the colors are defined in report.colors, which the stylesheets
# also read, and imported from here by the plotting modules.
from genomic_benchmarks_qc.report.colors import (  # noqa: F401
    FAIL_COLOR,
    PASS_COLOR,
    UNKNOWN_COLOR,
    WARN_COLOR,
)

# A figure is written to plots/ at print resolution and embedded in the page at
# half of it. The images are most of a report's bytes - the four a class report
# shows were 65% of an 890 KB one at 300 dpi - and nothing on the page is drawn
# wider than 92% of a column, so the embedded copy was carrying about four times
# the pixels it could use. The file on disk keeps the resolution a figure in a
# paper needs, and is the one to reuse elsewhere.
FIGURE_DPI = 300
DISPLAY_DPI = 150


@dataclass(frozen=True)
class SavedPlot:
    """One figure, as the file it was written to and as the copy the page shows.

    The two are the same figure at the two resolutions above. Everything that
    embeds an image takes either this or a path, so a plot that is only ever
    read back off disk - the split report's histogram, say - can stay a path.

    Attributes:
        path: The PNG in plots/, at `FIGURE_DPI`.
        embedded: The same figure at `DISPLAY_DPI`, as PNG bytes, or None for a
            figure the page does not show - the file is read instead, so a
            figure that later grows a place on the page still appears.
    """

    path: Path
    embedded: bytes | None


def put_file_details(html_template, filename):
    """
    Populates the placeholders {{filename}} and {{date}} in the HTML template.
    """
    # Replace {{filename}} with the stripped filename
    html_template = html_template.replace("{{filename}}", str(filename))

    # Replace {{date}} with the current date and time
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return html_template.replace("{{date}}", current_time)


def put_data(html_template, placeholder, data):
    """
    Replaces all occurrences of a placeholder in the HTML template with the provided data.

    The data goes in as-is, so this is for HTML the report built itself. Anything
    that came out of the input files goes through `put_text` instead.
    """
    if placeholder not in html_template:
        raise ValueError(f"Placeholder not found: {placeholder}")
    return html_template.replace(placeholder, str(data))


def escape_html_text(value):
    """Return `value` as HTML text that is inert both as markup and as a placeholder.

    Two things have to be neutralised in text the report did not write itself.

    Markup, which `html.escape` handles: a label of `<img src=x onerror=...>`
    would otherwise be markup by the time the page was opened, and a sequence
    containing `<` would silently swallow the rest of the cell.

    And placeholders. A page is built by filling `{{name}}` one name at a time,
    so a value that lands early and happens to contain `{{label2}}` is still
    sitting in the page when `{{label2}}` is filled, and gets filled too - a
    file called `{{label2}}.csv` puts the class name in the filename cell.
    `html.escape` leaves braces alone, so `{` is written as `&#123;` here. It
    renders as `{` and matches no placeholder.
    """
    return html.escape(str(value)).replace('{', '&#123;')


def put_text(html_template, placeholder, data):
    """Fill a placeholder with text from the data, escaped.

    Class labels, column names, file names and the bases themselves are read out
    of files the report has no say over, and a report is a page someone else
    opens. `escape_html_text` says what that costs them.

    The counterpart of `put_data`: which one a call site uses says where the
    value came from.
    """
    return put_data(html_template, placeholder, escape_html_text(data))


def escape_str(s):
    """Return a JSON/JS quoted string for safe embedding in inline scripts.

    `{` goes in as its JSON escape for the reason `escape_html_text` gives:
    the element this lands in is inserted into a page that still has
    placeholders left to fill. The escape is invisible to the reader - the
    string JavaScript parses out is the one that went in.
    """
    return (
        json.dumps(str(s))
        .replace("&", "\\u0026")
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("'", "\\u0027")
        .replace("{", "\\u007b")
    )


def icon_html(summary_statuses, key):
    """
    Return small circular icon HTML for a summary status mapping.

    summary_statuses: mapping or None
    key: lookup key
    """
    if summary_statuses is None:
        return ''
    val = summary_statuses.get(key, '')
    if not val:
        return ''
    s = str(val).strip()
    lv = s.lower()
    if lv in ('pass', 'ok', 'good', 'success'):
        return '<span class="status-icon status-pass">✔</span>'
    if lv in ('warn', 'warning'):
        return '<span class="status-icon status-warn">!</span>'
    if lv in ('fail', 'failed', 'error'):
        return '<span class="status-icon status-fail">✖</span>'
    if lv == 'unknown':
        # Not scored: too few sequences behind the check for its result to mean
        # anything. Deliberately distinct from Pass, which is a verdict.
        return ('<span class="status-icon status-unknown" '
                'title="Not enough sequences to score this check">?</span>')
    # Otherwise assume the value is an HTML snippet or a custom symbol and return as-is
    return s

def encode_image_to_base64(image):
    """
    Return an image's base64-encoded string, ready for a data URI.

    Args:
        image: A `SavedPlot`, whose display-resolution copy is used, or a path
            to an image file, which is read.

    Returns:
        Base64-encoded string of the image

    Raises:
        FileNotFoundError: If a path was given and there is no file at it
    """
    if isinstance(image, SavedPlot):
        if image.embedded is not None:
            return base64.b64encode(image.embedded).decode('utf-8')
        image = image.path
    if not Path(image).is_file():
        raise FileNotFoundError(f"Image file not found: {image}")
    with open(image, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode('utf-8')

def image_or_message(image_path, alt, css_class, message):
    """Return an <img> element for a plot, or a text message when there is none.

    Args:
        image_path: The plot, as a `SavedPlot` or a path, or None when no plot
            was generated.
        alt: Alt text / label for the image.
        css_class: Class for the <img>: 'plot-wide' for the figures whose axes
            line up down the page, 'plot-half' for the narrower centred ones.
            Sizing lives in the stylesheet so those two families stay in step
            with the interactive figure, which has to match them.
        message: Text to show instead of the image when image_path is None.

    Returns:
        HTML snippet as a string.
    """
    if image_path is None:
        return f'<p class="no-plot-message">{escape_html_text(message)}</p>'
    return (f'<img src="data:image/png;base64, {encode_image_to_base64(image_path)}" '
            f'alt="{escape_html_text(alt)}" class="{css_class}">')

# The stylesheets, inlined into every report. report.css is the original shared
# stylesheet; the other two are layered on top of it, so their rules win where
# they overlap - see the header comment in each file.
#
# report_ui.js is the behaviour that goes with them, and is added to the page by
# the report modules rather than from here, because it has to sit at the end of
# <body>.
COMMON_CSS = assets.stylesheet('report.css', 'report_design.css')

# Where a reader who was handed the file can find the tool. A report is one
# self-contained page and the way it gets shared is by sending it, so links out
# of it are the only route from a report back to the project that made it -
# which until now the page did not offer at all.
DOCS_URL = "https://genomic-benchmarks.github.io/genomic-benchmarks-qc/"
REPO_URL = "https://github.com/genomic-benchmarks/genomic-benchmarks-qc"

# The documentation pages a report links into. Kept here, rather than spelled
# out at each use, so tests/test_report_links.py can check every page and every
# anchor against the headings the docs actually have - a report outlives the
# build that made it, and a link into a renamed heading fails silently.
DOCS_PAGES = {
    'checks': 'guide/checks/',
    'flags': 'guide/how-it-works/',
    'viewer': 'guide/per-position/',
    'leakage': 'guide/leakage/',
}

# The tagline and the one-line description. Both strings also appear on the README,
# the docs landing page and PyPI; this is the copy the reports use, and the one
# place to change when the wording changes.
TOOL_TAGLINE = "Find the shortcut. Learn the biology."
TOOL_DESCRIPTION = (
    "Automated quality control for genomic machine learning datasets: scores the "
    "biases, duplicates and data leakage a classifier could exploit before you "
    "train on it."
)


def docs_url(page, anchor=None):
    """Absolute URL of a documentation page, optionally one heading on it."""
    if page not in DOCS_PAGES:
        raise ValueError(f"Unknown documentation page: {page}")
    return f"{DOCS_URL}{DOCS_PAGES[page]}" + (f"#{anchor}" if anchor else "")


def docs_link(page, anchor=None, text='Read more'):
    """The link that closes a check's ? explanation.

    The explanations say what a figure shows and what decides its flag, and stop
    there. What to do about a finding takes a page, and the guide has one per
    check - so each box ends with the way through to it instead of growing to
    hold it. That is what keeps them to a few lines.
    """
    return (f'<a class="explanation-more" href="{docs_url(page, anchor)}" '
            f'target="_blank" rel="noopener">{text} &rarr;</a>')


# The header of every report: what the tool is, for a reader who has only been
# sent the file, then what this particular report looked at and what it found.
# The verdict line doubles as the key to the coloured circles used down the whole
# page - each count carries the icon it counts - so the page needs no separate
# legend to explain them.
REPORT_HEADER_HTML = """
<div class="report-header">
    {{logo}}
    <div class="report-tagline">{{tool_tagline}}</div>
    <div class="report-lede">{{tool_description}}</div>
    <div class="report-subject">{{report_subject}}</div>
    {{report_verdict}}
    <div class="report-meta">gb-qc v{{version}} &middot; run {{generated_on}}</div>
</div>
"""

# The foot of the navigation. In the nav rather than the header because the nav
# is fixed: a reader who scrolled to a figure half way down and wants to know
# what the tool is can still see it.
SIDEBAR_LINKS_HTML = f"""
            <div class="sidebar-links">
                <a href="{DOCS_URL}" target="_blank" rel="noopener">Documentation</a>
                <a href="{REPO_URL}" target="_blank" rel="noopener">GitHub</a>
            </div>
"""

# Order the verdict reads in: what needs attention first, and Unknown last
# because it is the absence of a result rather than one.
VERDICT_ORDER = (
    ('fail', 'Fail', '&#10006;'),
    ('warn', 'Warning', '!'),
    ('pass', 'Pass', '&#10004;'),
    ('unknown', 'Unknown', '?'),
)

_FLAG_BUCKETS = {
    'pass': 'pass', 'ok': 'pass', 'good': 'pass', 'success': 'pass',
    'warn': 'warn', 'warning': 'warn',
    'fail': 'fail', 'failed': 'fail', 'error': 'fail',
    'unknown': 'unknown',
}


def verdict_html(summary_statuses, check_names):
    """One line saying how the report came out, with the flag icons in it.

    Counted from the same flags the navigation icons are drawn from, over the
    headline checks only - the sub-checks behind a headline (a nucleotide, a
    dinucleotide, a position) are what its own flag already summarises, and
    counting them here would put a four-figure total on a nine-check report.

    Args:
        summary_statuses: mapping of check name to flag, or None.
        check_names: the headline checks, in the order the navigation lists them.

    Returns:
        HTML snippet, or an empty string when nothing was flagged either way.
    """
    if not summary_statuses:
        return ''
    scored = [
        (name, _FLAG_BUCKETS[key])
        for name in check_names
        for key in [str(summary_statuses.get(name, '')).strip().lower()]
        if key in _FLAG_BUCKETS
    ]
    if not scored:
        return ''
    counts = {}
    for _, bucket in scored:
        counts[bucket] = counts.get(bucket, 0) + 1
    total = len(scored)

    def part(key, label, symbol, count=None):
        counted = '' if count is None else f'{count} '
        return (f'<span class="verdict-part">'
                f'<span class="status-icon-small status-{key}">{symbol}</span>'
                f'{counted}{label}</span>')

    if total == 1:
        # A report with one check has nothing to tally, so the line names the
        # check and its flag instead of counting to one.
        name, bucket = scored[0]
        key, label, symbol = next(e for e in VERDICT_ORDER if e[0] == bucket)
        lead = (f'<span class="verdict-total">{escape_html_text(name)}:</span>'
                + part(key, label, symbol))
    else:
        lead = (f'<span class="verdict-total">{total} checks:</span>'
                + ' '.join(part(key, label, symbol, counts[key])
                           for key, label, symbol in VERDICT_ORDER
                           if counts.get(key)))
    return ('<div class="report-verdict">' + lead
            + docs_link('flags', text='How a flag is decided') + '</div>')

# The logo, once. `put_data` fills every occurrence of a placeholder, so a data
# URI written into the markup was embedded once per <img> - the class report has
# two, the nav and the header - which was 150 KB of an 850 KB report. It goes
# into the stylesheet as `LOGO_CSS` instead, where both elements read the one
# copy, and it is stored at the size the page shows it rather than at ten times
# that: 1462x514 shown at 250 px, now 500x176.
LOGO_BASE64 = "data:image/png;base64,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"

# Only the picture. Where each of the two elements sits and how wide it is
# belongs with the rest of the page's sizing, in report_design.css.
LOGO_CSS = f'.gb-qc-logo {{ background-image: url("{LOGO_BASE64}"); }}'

# The header's element. The nav's is the same picture in the page templates,
# which is where the rest of the page's markup lives. `role` and `aria-label`
# keep it an image to a screen reader, which is the job the alt text was doing.
REPORT_LOGO_HTML = ('<div class="gb-qc-logo report-logo" role="img" '
                    'aria-label="Genomic Benchmarks QC logo"></div>')
