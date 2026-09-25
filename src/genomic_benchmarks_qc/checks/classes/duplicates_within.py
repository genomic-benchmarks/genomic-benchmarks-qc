"""Sequence Duplications within Labels: how much of the data survives deduplication.

Decided by a rule, not a score, pooled over both classes: Pass only when every
sequence is unique, Warning while at least 98% are, Fail below that.
"""

from genomic_benchmarks_qc.checks import Check, CheckResult

NAME = 'Sequence Duplications within Labels'


def score(stats1, stats2) -> CheckResult:
    """Flag the combined share of sequences left after deduplication.

    The row also carries that share as 'Percent Remaining', which the report's
    figure is drawn from and `gb-qc-report.csv` leaves out.
    """
    total_sequences = 0
    unique_sequences = 0
    for stats in [stats1, stats2]:
        total_sequences += stats.stats['Number of sequences']
        unique_sequences += stats.stats['Number of sequences left after deduplication']

    percent_remaining = unique_sequences / total_sequences if total_sequences > 0 else 1.0

    if percent_remaining < 0.98:
        flag = 'Fail'
    elif percent_remaining < 1.0:
        flag = 'Warning'
    else:
        flag = 'Pass'

    return CheckResult({NAME: {'Flag': flag, 'Percent Remaining': percent_remaining}})


CHECK = Check(name=NAME, score=score)
