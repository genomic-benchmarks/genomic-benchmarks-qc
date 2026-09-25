"""A check's section of the HTML report, and the shell every section shares.

Every section has the same shape: a heading carrying the check's flag in a
circle and a ? button, the explanation the button opens, and a body - a figure,
a table, a message. The shape is written once, in `check_section.html`, and so
is the navigation entry that links to it. What differs from check to check is
the check's own: its explanation and body markup in a fragment of its own,
`check_<name>.html`, and a `render` function that draws its figures and fills
the fragment in.

That is what lets a page be assembled from its registry. The page asks each
check for its navigation entry and its section, in registry order, and a check
the page template has never heard of appears in both.
"""

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.utils import docs_link, icon_html, put_data, put_text

if TYPE_CHECKING:
    import pandas as pd

    from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics

SECTION_TEMPLATE = assets.template('check_section.html')
NAV_ITEM_TEMPLATE = assets.template('check_nav_item.html')

# Said in place of every figure that needs a base both classes have, when they
# have none: the per-sequence composition checks and the per-position ones.
DISJOINT_BASES_MESSAGE = ('No bases in common between the two labels, so this comparison '
                          'cannot be plotted.')


@dataclass(frozen=True)
class SectionContent:
    """What a check's `render` hands back for its section.

    Attributes:
        fills: Placeholder in the check's fragment -> the markup that replaces
            it. The markup goes in as it is, so anything that came out of the
            input files has to have gone through `escape_html_text` first.
        scripts: `<script>` elements the section needs on the page, placed at
            the end of the body. Each is included once, however many sections
            ask for it.
    """

    fills: dict[str, str] = field(default_factory=dict)
    scripts: tuple[str, ...] = ()


@dataclass(frozen=True)
class Section:
    """How one check appears in the report.

    Attributes:
        title: The section's heading.
        anchor: The section's id, which its navigation entry links to. Keep it
            once published: `report_ui.js`, the tests and a reader's bookmark
            into a report find a section by it.
        explanation_id: The id of the explanation block the ? button toggles.
        template: The check's fragment in `report.assets`. It holds the
            explanation block, `<div id="{{explanation_id}}"
            class="explanation-text">`, closing with `{{docs_link}}`, and after
            it the section's body. Every other placeholder in it is the
            render's to fill.
        render: Called with the report's context - a `ClassReportContext` or a
            `SplitReportContext` - to draw the check's figures, write any file
            it keeps beside the report, and return its `SectionContent`.
        docs: The links that close the explanation, as `docs_link` arguments
            `(page, anchor, text)`, in the order they appear.
        nav_label: What the navigation calls the section, where that is not
            its title.
        css_class: A class for the <section> element, e.g. 'table-section'.
        stylesheets: Stylesheets in `report.assets` the section needs on the
            page, included once however many sections ask for them.
    """

    title: str
    anchor: str
    explanation_id: str
    template: str
    render: Callable[[Any], SectionContent]
    docs: tuple[tuple[str, str | None, str], ...] = ()
    nav_label: str | None = None
    css_class: str | None = None
    stylesheets: tuple[str, ...] = ()


@dataclass(frozen=True, eq=False)
class ClassReportContext:
    """Everything a class check's section can draw on.

    `bases_overlap` and the drawn window are resolved once, by the report, so
    every figure and the interactive per-position viewer are drawn over the same
    bases and the same positions.

    Attributes:
        stats1, stats2: The two classes.
        results: Every check's metrics, indexed by check name.
        summary_statuses: Every check's flag, by check name.
        failed_by_feature: What each figure shades, as
            `flag_significant_differences` returns it, or None.
        plots_dir: The report's plots/ directory, which exists.
        plot_type: 'boxen' or 'violin', for the distribution plots.
        report_dir: The directory the report is written to.
        bases_overlap: The bases both classes have, sorted.
        drawn_end: Last position the per-position figures draw, as
            `drawn_window` resolves it.
        compared: Whether that window is the compared one, from `drawn_window`.
    """

    stats1: 'SequenceStatistics'
    stats2: 'SequenceStatistics'
    results: 'pd.DataFrame'
    summary_statuses: dict
    failed_by_feature: dict | None
    plots_dir: Path
    plot_type: str
    report_dir: Path
    bases_overlap: list
    drawn_end: int
    compared: bool

    def flagged(self, name):
        """What check `name` flagged, for its figure: None when nothing was given."""
        return self.failed_by_feature.get(name, {}) if self.failed_by_feature else None

    @property
    def label1(self):
        return str(self.stats1.label) if self.stats1.label is not None else 'N/A'

    @property
    def label2(self):
        return str(self.stats2.label) if self.stats2.label is not None else 'N/A'


@dataclass(frozen=True, eq=False)
class SplitReportContext:
    """Everything a split check's section can draw on.

    Attributes:
        basic_stats: The two halves' counts and lengths, from
            `get_basic_stats_from_aggregates`.
        threshold_stats: The leakage counts and percentages, from
            `get_threshold_stats`.
        results_filt: The hits the page lists, with their aligned sequences.
        leaked_hits: How many hits reached the threshold, listed or not.
        query_similarity_max, target_similarity_max: Each test and each train
            sequence's best similarity, NaN for one without a hit.
        plots_dir: The report's plots/ directory, which exists.
    """

    basic_stats: dict
    threshold_stats: dict
    results_filt: 'pd.DataFrame'
    leaked_hits: int
    query_similarity_max: Any
    target_similarity_max: Any
    plots_dir: Path


def render_nav(checks, summary_statuses) -> str:
    """The navigation entries of a registry, each carrying its check's flag."""
    items = []
    for check in checks:
        section = check.section
        item = put_data(NAV_ITEM_TEMPLATE, '{{icon}}', icon_html(summary_statuses, check.name))
        item = put_data(item, '{{anchor}}', section.anchor)
        items.append(put_text(item, '{{nav_label}}', section.nav_label or section.title))
    return ''.join(items)


def render_section(check, context, summary_statuses) -> tuple[str, tuple[str, ...]]:
    """One check's section, and the scripts it needs on the page.

    The fragment is filled first - its generic placeholders, then the render's -
    and goes into the shell last, so nothing a render wrote can be mistaken for
    one of the shell's placeholders.
    """
    section = check.section
    content = section.render(context)

    body = put_data(assets.template(section.template), '{{explanation_id}}',
                    section.explanation_id)
    body = put_data(body, '{{docs_link}}', ''.join(docs_link(*link) for link in section.docs))
    for placeholder, markup in content.fills.items():
        body = put_data(body, placeholder, markup)

    shell = put_data(SECTION_TEMPLATE, '{{anchor}}', section.anchor)
    shell = put_data(shell, '{{section_class}}',
                     f' class="{section.css_class}"' if section.css_class else '')
    shell = put_data(shell, '{{icon}}', icon_html(summary_statuses, check.name))
    shell = put_text(shell, '{{title}}', section.title)
    shell = put_data(shell, '{{explanation_id}}', section.explanation_id)
    return put_data(shell, '{{section_body}}', body), content.scripts


def render_sections(checks, context, summary_statuses) -> tuple[str, tuple[str, ...]]:
    """Every section of a registry, in order, and the scripts they need.

    Returns:
        Tuple of (the sections' markup, their scripts), each script once and in
        the order it was first asked for. Where the scripts go relative to the
        page's own is the page's decision.
    """
    sections, scripts = [], []
    for check in checks:
        markup, needed = render_section(check, context, summary_statuses)
        sections.append(markup)
        scripts.extend(script for script in needed if script not in scripts)
    return '\n'.join(sections), tuple(scripts)


def section_stylesheets(checks) -> str:
    """The stylesheets a registry's sections need, each once, in registry order."""
    names = dict.fromkeys(name for check in checks for name in check.section.stylesheets)
    return assets.stylesheet(*names)
