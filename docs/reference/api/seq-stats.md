# Sequence statistics

Where the per-sequence and per-position statistics are computed, and where the
constants that decide how far the per-position checks reach are defined.

::: genomic_benchmarks_qc.utils.seq_stats
    options:
      members:
        - SequenceStatistics
        - cohort_floor
        - DEFAULT_MIN_COVERAGE
        - MIN_SEQUENCES_PER_REPORTED_POSITION
      filters:
        - "!^_"
        # Every name here is a constructor argument, already documented in the
        # Args table under the class heading. Repeating it as a bare attribute
        # costs a heading and says nothing. The attributes that survive are the
        # ones `compute` resolves, which is what a caller reads back.
        - "!^(filename|filepath|label|slug|seq_column|sequences|min_coverage)$"
