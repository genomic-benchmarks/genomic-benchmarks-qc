"""The checks, one module each, and the registries that list them in report order.

A check is everything the tool says about one question: the name it is reported
under, the feature it compares and how that comparison becomes a flag, and the
section of the report that shows it. All of that lives in the check's own module,
so adding a check is writing one module and listing it in a registry - not
editing every place that used to name the checks one by one.

There is one registry per command, because the two compare different things:

- `checks.classes.CLASS_CHECKS` - two classes of one dataset, each a
  `SequenceStatistics`. Nine checks.
- `checks.splits.SPLIT_CHECKS` - the test half of a split searched against the
  train half. One check, Data Leakage.

The order of a registry is the order of its report - the navigation, the
sections, and the rows of `gb-qc-report.csv`.

This module imports neither registry. The class checks are built on the scorers
in `genomic_benchmarks_qc.utils.testing`, which imports the registry back to run
it, so the registries are imported where they are used rather than here.
"""

from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from genomic_benchmarks_qc.report.sections import Section

# The flags a sub-check is shaded for on its check's figure. Pass is not drawn -
# it would cover the whole figure - and Unknown is not a finding.
FLAGGED = ('Fail', 'Warning')


@dataclass(frozen=True)
class CheckResult:
    """What one check found, on one comparison.

    Attributes:
        rows: One entry per row of `gb-qc-report.csv`, in the order the check
            produced them. The headline is keyed by the check's own name and its
            sub-checks by `'<name> - <what>'` - a base, a dinucleotide, a base at
            a position. A check can come back without a headline, when there
            was nothing to compare (no bases, no positions).
        flagged: What the check's figure shades - only the sub-checks that came
            out Warning or Fail, keyed by what was flagged: `{'A': 'Warning'}`
            for a per-sequence check, `{'G': {52: 'Fail'}}` for a per-position
            one. None for a check with nothing to shade, which keeps it out of
            `failed_by_feature` altogether.
    """

    rows: dict[str, dict]
    flagged: dict | None = None


@dataclass(frozen=True)
class Check:
    """One check: its name, how it is scored, and how the report shows it.

    Attributes:
        name: What the check is reported as - the row in `gb-qc-report.csv`, the
            key in the results, the name in a terminal warning. Never contains
            `' - '`, which is what sets a sub-check apart from its headline.
        score: Scores one comparison and returns a `CheckResult`. For a class
            check it takes the two `SequenceStatistics`; for a split check, the
            threshold statistics of the search.
        section: The check's section of the HTML report.
        floor: Which minimum a check has to clear before it is scored, for a
            check that can report Unknown: `'per_sequence'` for one compared
            only when the smaller class holds `MIN_SEQUENCES_PER_CLASS`
            sequences, `'per_position'` for one compared position by position on
            the sequences reaching each position. None for a check decided by a
            rule rather than by a score, which is always decided. It chooses how
            an Unknown is explained, on the terminal and in the report.
    """

    name: str
    score: Callable[..., CheckResult]
    section: 'Section | None' = None
    floor: Literal['per_sequence', 'per_position'] | None = None


def run_checks(checks, *inputs) -> tuple[dict[str, dict], dict[str, dict]]:
    """Run every check of a registry on one comparison.

    Args:
        checks: The registry, in report order.
        *inputs: What each check's `score` takes.

    Returns:
        Tuple of `(results, failed_by_feature)`. `results` is every row of every
        check: the headlines first, in registry order, then each check's
        sub-checks, in the order the check produced them - which is the order of
        `gb-qc-report.csv`. `failed_by_feature` maps each check that has a figure
        to shade to its `CheckResult.flagged`.
    """
    scored = [(check, check.score(*inputs)) for check in checks]

    results = {}
    for check, result in scored:
        if check.name in result.rows:
            results[check.name] = result.rows[check.name]
    for check, result in scored:
        results.update((name, row) for name, row in result.rows.items() if name != check.name)

    failed_by_feature = {check.name: result.flagged for check, result in scored
                         if result.flagged is not None}
    return results, failed_by_feature
