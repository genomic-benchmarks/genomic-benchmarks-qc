"""Unique bases: do the two classes write their sequences in the same characters?

Decided by a rule, not a score. A character only one class uses makes every
sequence carrying it perfectly classifiable, so any difference in the two sets
fails, and there is no Warning.
"""

from genomic_benchmarks_qc.checks import Check, CheckResult
from genomic_benchmarks_qc.report.sections import Section, SectionContent
from genomic_benchmarks_qc.report.utils import escape_html_text

NAME = 'Unique bases'


def score(stats1, stats2) -> CheckResult:
    """Fail when the two classes' sets of characters differ at all."""
    same_bases = set(stats1.stats['Unique bases']) == set(stats2.stats['Unique bases'])
    return CheckResult({NAME: {'Flag': 'Pass' if same_bases else 'Fail'}})


def render(context) -> SectionContent:
    """The two sets of characters, side by side."""
    return SectionContent({
        '{{label1}}': escape_html_text(context.label1),
        '{{label2}}': escape_html_text(context.label2),
        '{{unique_bases1}}': escape_html_text(', '.join(context.stats1.stats['Unique bases'])),
        '{{unique_bases2}}': escape_html_text(', '.join(context.stats2.stats['Unique bases'])),
    })


CHECK = Check(
    name=NAME,
    score=score,
    section=Section(
        title='Unique Bases',
        anchor='unique-bases',
        explanation_id='unique-bases-explanation',
        template='check_unique_bases.html',
        render=render,
        docs=(('checks', 'unique-bases', 'What to do about it'),),
        css_class='table-section',
    ),
)
