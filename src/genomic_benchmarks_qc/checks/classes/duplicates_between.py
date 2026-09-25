"""Duplicate Sequences between Labels: the same sequence filed under both classes.

Decided by a rule, not a score. Identical input with opposite labels is a pair
no model can get both of right, so one shared sequence fails, and there is no
Warning.
"""

from genomic_benchmarks_qc.checks import Check, CheckResult

NAME = 'Duplicate Sequences between Labels'


def score(stats1, stats2) -> CheckResult:
    """Fail when any sequence occurs in both classes."""
    shared = bool(set(stats1.sequences) & set(stats2.sequences))
    return CheckResult({NAME: {'Flag': 'Fail' if shared else 'Pass'}})


CHECK = Check(name=NAME, score=score)
