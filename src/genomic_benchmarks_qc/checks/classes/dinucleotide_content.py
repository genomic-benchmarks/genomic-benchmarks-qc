"""Per sequence dinucleotide content: the nucleotide check, for two-base combinations.

Each dinucleotide is scored separately, as the AU-ROC of its frequency in a
sequence, and the headline is the worst of them.
"""

from genomic_benchmarks_qc.checks import Check
from genomic_benchmarks_qc.checks.scorers import column_features, stat

NAME = 'Per sequence dinucleotide content'

CHECK = Check(
    name=NAME,
    score=column_features(NAME, stat('Per sequence dinucleotide content')),
    floor='per_sequence',
)
