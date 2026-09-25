"""Duplicate Sequences between Labels: the same sequence filed under both classes.

Decided by a rule, not a score. Identical input with opposite labels is a pair
no model can get both of right, so one shared sequence fails, and there is no
Warning.
"""

import logging

from genomic_benchmarks_qc.checks import Check, CheckResult
from genomic_benchmarks_qc.report.sections import Section, SectionContent
from genomic_benchmarks_qc.report.utils import escape_html_text, escape_str
from genomic_benchmarks_qc.utils.naming import DUPLICATES_FILE

logger = logging.getLogger(__name__)

NAME = 'Duplicate Sequences between Labels'

# How many of the shared sequences the page lists; all of them go to the file.
LISTED = 10


def score(stats1, stats2) -> CheckResult:
    """Fail when any sequence occurs in both classes."""
    shared = bool(set(stats1.sequences) & set(stats2.sequences))
    return CheckResult({NAME: {'Flag': 'Fail' if shared else 'Pass'}})


def render(context) -> SectionContent:
    """The first shared sequences on the page, and every one of them in a file beside it."""
    duplicate_seqs = list(set(context.stats1.sequences).intersection(context.stats2.sequences))
    # The report lives in a directory of its own, so this is a plain sibling.
    duplicate_seqs_path = context.report_dir / DUPLICATES_FILE

    if not duplicate_seqs:
        return SectionContent({'{{duplicates}}': """
        <p>No duplicate sequences were found between classes.</p>
        """})

    with open(duplicate_seqs_path, 'w') as handle:
        for seq in sorted(duplicate_seqs):
            handle.write(f"{seq}\n")
    logger.info(f"Duplicate sequences saved to {duplicate_seqs_path}")

    rest = f"And {len(duplicate_seqs) - LISTED} more." if len(duplicate_seqs) > LISTED else ""
    saved = (f"All {len(duplicate_seqs)} duplicate sequences saved to "
             f"{escape_html_text(str(duplicate_seqs_path))}.")
    # The sequences go in as JSON data rather than as a JavaScript literal, so
    # report_ui.js stays a file that can be linted and a sequence cannot break
    # out of the script element. Escaped for the same reason as the
    # per-position payload.
    escaped_seqs = [escape_str(seq) for seq in duplicate_seqs[:LISTED]]
    table = f"""
                <table>
                    <thead>
                        <tr>
                            <th class="sequence_column">Sequence</th>
                        </tr>
                    </thead>
                        <tbody>
                            <!-- Table rows will be dynamically populated -->
                        </tbody>
                </table>
                <div id="sequence-duplication-levels-info">
                    <p>{rest} {saved}</p>
                </div>
        """
    table += ('\n<script type="application/json" id="duplicate-sequences">['
              + ', '.join(escaped_seqs) + ']</script>')
    return SectionContent({'{{duplicates}}': table})


CHECK = Check(
    name=NAME,
    score=score,
    section=Section(
        title='Duplicate Sequences between Labels',
        # report_ui.js fills the table of shared sequences through this id.
        anchor='sequence-duplication-levels',
        explanation_id='between-dup-explanation',
        template='check_duplicates_between.html',
        render=render,
        docs=(('checks', 'duplicate-sequences-between-labels', 'What to do about it'),),
    ),
)
