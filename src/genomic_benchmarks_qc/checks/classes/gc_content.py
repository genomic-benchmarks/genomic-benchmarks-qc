"""Per sequence GC content: can a sequence's GC% alone tell the classes apart?

The flag is the AU-ROC of GC content on its own.
"""

from genomic_benchmarks_qc.checks import Check
from genomic_benchmarks_qc.checks.scorers import scalar_feature, stat

NAME = 'Per sequence GC content'

CHECK = Check(
    name=NAME,
    score=scalar_feature(NAME, stat('Per sequence GC content')),
    floor='per_sequence',
)
