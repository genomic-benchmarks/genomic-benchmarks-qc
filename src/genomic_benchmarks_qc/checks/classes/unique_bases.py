"""Unique bases: do the two classes write their sequences in the same characters?

Decided by a rule, not a score. A character only one class uses makes every
sequence carrying it perfectly classifiable, so any difference in the two sets
fails, and there is no Warning.
"""

from genomic_benchmarks_qc.checks import Check, CheckResult

NAME = 'Unique bases'


def score(stats1, stats2) -> CheckResult:
    """Fail when the two classes' sets of characters differ at all."""
    same_bases = set(stats1.stats['Unique bases']) == set(stats2.stats['Unique bases'])
    return CheckResult({NAME: {'Flag': 'Pass' if same_bases else 'Fail'}})


CHECK = Check(name=NAME, score=score)
