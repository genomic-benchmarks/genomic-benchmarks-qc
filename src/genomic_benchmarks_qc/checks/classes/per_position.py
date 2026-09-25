"""Per position nucleotide content, read forward and from the sequence end.

Two checks built by one factory: at every position, can the base there tell the
classes apart? Each base at each position is scored separately, on the
sequences long enough to reach the position, and the headline is the worst
single one. The reversed check counts positions from the last base, which is
what catches something anchored to the 3' end of variable-length sequences.

Which positions are compared at all is decided by the cohort floors in
`SequenceStatistics`; `position_windows` combines the two classes' windows.
"""

from genomic_benchmarks_qc.checks import FLAGGED, Check, CheckResult
from genomic_benchmarks_qc.utils.testing import (
    _aggregate_worst_case_metrics,
    _score_position_features,
    position_windows,
)

FORWARD_NAME = 'Per position nucleotide content'
REVERSED_NAME = 'Per position reversed nucleotide content'


def _flagged_positions(rows, name, bases, end_position) -> dict:
    """The positions a per-position check flagged, keyed by base then position.

    Only flagged positions go in, and a base only once it has one. An entry per
    base regardless would leave the dict truthy for a comparison with nothing to
    shade, and the report would draw a second, identical copy of every
    per-position plot.

    Args:
        rows: The check's rows, as `_score_position_features` names them.
        name: The check's name.
        bases: The bases scored, in the order they were scored.
        end_position: Last position scored, 1-based and inclusive.

    Returns:
        Dict of base -> {position: flag}, positions 1-based.
    """
    flagged = {}
    for base in bases:
        for position in range(1, end_position + 1):
            flag = rows[f'{name} - {base} position {position}']['Flag']
            if flag in FLAGGED:
                flagged.setdefault(base, {})[position] = flag
    return flagged


def position_features(name, reverse):
    """Score every base at every position, reading from the start or the end."""
    def score(stats1, stats2) -> CheckResult:
        bases = sorted(set(stats1.stats['Unique bases']) | set(stats2.stats['Unique bases']))
        end_position, scored_end_position = position_windows(stats1, stats2)

        rows, per_base = _score_position_features(
            stats1.sequences, stats2.sequences, bases, name,
            end_position=end_position, scored_end_position=scored_end_position,
            reverse=reverse,
        )
        if per_base:
            rows[name] = _aggregate_worst_case_metrics(per_base.values())
        return CheckResult(rows, _flagged_positions(rows, name, bases, end_position))
    return score


FORWARD = Check(
    name=FORWARD_NAME,
    score=position_features(FORWARD_NAME, reverse=False),
    floor='per_position',
)

REVERSED = Check(
    name=REVERSED_NAME,
    score=position_features(REVERSED_NAME, reverse=True),
    floor='per_position',
)
