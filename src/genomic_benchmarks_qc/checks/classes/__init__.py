"""The checks `evaluate-classes` runs on every pair of classes, in report order.

The order here is the order of the report's navigation and sections and of the
rows of `gb-qc-report.csv`: the checks on the data as it stands first, then the
checks that score a feature. To add a check, write its module beside these and
list it here - see "Adding a check" in CONTRIBUTING.md.
"""

from genomic_benchmarks_qc.checks.classes import (
    dinucleotide_content,
    duplicates_between,
    duplicates_within,
    gc_content,
    lengths,
    nucleotide_content,
    per_position,
    unique_bases,
)

CLASS_CHECKS = (
    unique_bases.CHECK,
    duplicates_within.CHECK,
    duplicates_between.CHECK,
    lengths.CHECK,
    gc_content.CHECK,
    nucleotide_content.CHECK,
    dinucleotide_content.CHECK,
    per_position.FORWARD,
    per_position.REVERSED,
)
