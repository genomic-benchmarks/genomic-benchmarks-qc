"""The checks `evaluate-splits` runs on a train and test split, in report order.

To add a check, write its module beside this one and list it here - see "Adding
a check" in CONTRIBUTING.md. A split check scores the threshold statistics of
the search, `get_threshold_stats`, and draws its section from a
`SplitReportContext`.
"""

from genomic_benchmarks_qc.checks.splits import data_leakage

SPLIT_CHECKS = (
    data_leakage.CHECK,
)
