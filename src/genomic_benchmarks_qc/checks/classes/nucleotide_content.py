"""Per sequence nucleotide content: can how often one base occurs tell the classes apart?

Each base is scored separately, as the AU-ROC of its frequency in a sequence, and
the headline is the worst of them.
"""

from genomic_benchmarks_qc.checks import Check
from genomic_benchmarks_qc.checks.scorers import column_features, stat

NAME = 'Per sequence nucleotide content'

CHECK = Check(
    name=NAME,
    score=column_features(NAME, stat('Per sequence nucleotide content')),
    floor='per_sequence',
)
