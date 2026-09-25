"""How a feature becomes a scored, flagged check - the part every AU-ROC check shares.

A check that asks "does this feature tell the two classes apart?" names the
feature and one of these factories scores it, so the scoring and the flag
boundaries are the same for every such check: the AU-ROC of the feature read as
a classifier's score, flagged by `_flag_on_score`, with the worst sub-check as
the headline. The arithmetic itself is in `genomic_benchmarks_qc.utils.testing`;
this module only adapts it to the shape of a `Check`.

A feature is a function of one class's `SequenceStatistics` returning a
DataFrame with one row per sequence. `stat` reads one that `SequenceStatistics`
already computes; a new check can just as well compute its own from
`stats.sequences`.
"""

from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from genomic_benchmarks_qc.checks import FLAGGED, CheckResult
from genomic_benchmarks_qc.utils.testing import _score_dataframe_features, _score_scalar_feature

if TYPE_CHECKING:
    from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics

Feature = Callable[['SequenceStatistics'], pd.DataFrame]


def stat(key: str) -> Feature:
    """The feature `SequenceStatistics.compute` stores under `key`."""
    return lambda stats: stats.stats[key]


def _indices(stats) -> np.ndarray:
    """Every sequence of a class, as the row indices the scorers select with."""
    return np.arange(len(stats.sequences))


def scalar_feature(name: str, feature: Feature) -> Callable[..., CheckResult]:
    """Score a feature with one value per sequence: one row, the headline.

    Args:
        name: The check's name, which keys its row.
        feature: One class's values, as a DataFrame with a single column.
    """
    def score(stats1, stats2) -> CheckResult:
        frame1, frame2 = feature(stats1), feature(stats2)
        return CheckResult({name: _score_scalar_feature(
            frame1, frame2, frame1.columns[0], _indices(stats1), _indices(stats2))})
    return score


def column_features(name: str, feature: Feature) -> Callable[..., CheckResult]:
    """Score a feature with several values per sequence, one sub-check per column.

    The headline is the worst of the columns. The figure shades the columns that
    came out Warning or Fail.

    Args:
        name: The check's name; each column's row is `'<name> - <column>'`.
        feature: One class's values, one column per thing compared - a base, a
            dinucleotide. A column one class lacks counts as zero there.
    """
    def score(stats1, stats2) -> CheckResult:
        frame1, frame2 = feature(stats1), feature(stats2)
        rows = _score_dataframe_features(frame1, frame2, name,
                                         _indices(stats1), _indices(stats2))
        flagged = {}
        for column in sorted(set(frame1.columns) | set(frame2.columns)):
            flag = rows[f'{name} - {column}']['Flag']
            if flag in FLAGGED:
                flagged[column] = flag
        return CheckResult(rows, flagged)
    return score
