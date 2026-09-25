"""Sequence lengths: can the length of a sequence alone tell the classes apart?

The flag is the AU-ROC of length on its own - how well a model that does nothing
but count characters separates the two classes.
"""

from genomic_benchmarks_qc.checks import Check
from genomic_benchmarks_qc.checks.scorers import scalar_feature, stat

NAME = 'Sequence lengths'

CHECK = Check(
    name=NAME,
    score=scalar_feature(NAME, stat('Sequence lengths')),
    floor='per_sequence',
)
