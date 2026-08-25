"""Tests for the shared plot furniture.

The plot functions themselves are exercised end to end by the pipeline tests;
what is pinned here is the styling that every plot shares, where a mistake is
invisible in the output and silently consistent across the whole report.
"""

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest

from genomic_benchmarks_qc.report.classes_plots import prepare_legend
from genomic_benchmarks_qc.report.splits_plots import _binned_counts

# Axis labels are 14pt and tick labels 12pt throughout, so the legend sits with
# the tick labels.
LEGEND_FONT_SIZE = 12


@pytest.fixture
def axis():
    figure, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label='class A')
    ax.plot([0, 1], [1, 0], label='class B')
    yield ax
    plt.close(figure)


class TestPrepareLegend:
    def test_legend_text_is_the_intended_point_size(self, axis):
        """matplotlib accepts a number or a size keyword here and silently drops
        anything else, so a quoted '12' left every legend in both reports at the
        10pt default. Reading the size back off the artist is the only way to
        tell, which is why it went unnoticed."""
        prepare_legend(axis)

        sizes = {text.get_fontsize() for text in axis.get_legend().get_texts()}
        assert sizes == {LEGEND_FONT_SIZE}

    def test_every_series_is_labelled(self, axis):
        prepare_legend(axis)

        assert [text.get_text() for text in axis.get_legend().get_texts()] == ['class A', 'class B']

    def test_explicit_handles_and_labels_are_used(self, axis):
        handles, _ = axis.get_legend_handles_labels()
        prepare_legend(axis, legend_handles=handles[:1], legend_labels=['only this one'])

        assert [text.get_text() for text in axis.get_legend().get_texts()] == ['only this one']
        assert axis.get_legend().get_texts()[0].get_fontsize() == LEGEND_FONT_SIZE


class TestSimilarityBinCounts:
    """What the split histogram's ten bars are counted from.

    The bars are counted rather than drawn from the observations, so the two
    rules the figure depends on - where a sequence with no hit goes, and which
    bin a perfect match lands in - are here rather than visible in a PNG.
    """

    def test_sequences_with_no_hit_are_counted_at_the_bottom(self):
        """They are absent from the search output, not absent from the half.

        Leaving them out would draw a split as more related to itself than it
        is, which is the flattering direction.
        """
        counts = _binned_counts(np.array([95.0, 95.0]), without_hits=500)

        assert counts[0] == 500
        assert counts[-1] == 2
        assert counts.sum() == 502

    def test_a_perfect_match_lands_in_the_last_bin(self):
        """The last bin is closed at both ends, as its tick label says."""
        counts = _binned_counts(np.array([100.0, 90.0, 89.999]), without_hits=0)

        assert counts[-1] == 2
        assert counts[-2] == 1

    def test_sequences_the_search_said_nothing_about_are_skipped(self):
        """NaN is 'no hit', and it is `without_hits` that counts those."""
        counts = _binned_counts(np.array([np.nan, 45.0, np.nan]), without_hits=2)

        assert counts.sum() == 3
        assert counts[4] == 1
        assert counts[0] == 2
