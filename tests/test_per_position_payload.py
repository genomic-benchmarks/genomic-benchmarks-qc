"""Tests for the data behind the interactive per-position figure.

The figure is drawn in the browser from this payload, so anything the payload
gets wrong is invisible until someone opens the report: a series in the wrong
order, positions off by one against the flags, or a not-scored position that
arrives looking like a real measurement.
"""

import json

import numpy as np
import pandas as pd
import pytest

from genomic_benchmarks_qc.report.colors import CLASS_COLORS
from genomic_benchmarks_qc.report.per_position_payload import (
    FEATURE_NAMES,
    _coverage,
    build_payload,
    payload_script,
    viewer_html,
)
from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics
from genomic_benchmarks_qc.utils.testing import flag_significant_differences

BASES = ['A', 'C', 'G', 'T']


def make_stats(sequences, label='cls'):
    stats = SequenceStatistics(sequences=sequences, filename=f'{label}.fa',
                               filepath=f'/{label}.fa', label=label)
    stats.compute()
    return stats


def random_sequences(count, length, seed, composition=None):
    rng = np.random.default_rng(seed)
    return [''.join(rng.choice(list('ACGT'), size=length, p=composition)) for _ in range(count)]


@pytest.fixture(scope='module')
def comparison():
    """A comparison with a real difference at every position, so flags exist."""
    stats1 = make_stats(random_sequences(60, 30, seed=1), label='a')
    stats2 = make_stats(random_sequences(60, 30, seed=2, composition=[.7, .1, .1, .1]), label='b')
    results, _ = flag_significant_differences(stats1, stats2)
    return stats1, stats2, pd.DataFrame.from_dict(results, orient='index')


@pytest.fixture(scope='module')
def payload(comparison):
    stats1, stats2, results = comparison
    end = min(stats1.end_position, stats2.end_position)
    return build_payload(stats1, stats2, BASES, end, results, 'forward')


class TestShape:
    def test_the_series_are_as_long_as_the_window(self, payload):
        end = payload['endPosition']

        assert len(payload['coverage']) == end
        for base in BASES:
            per_label = payload['freq'][base]
            assert len(per_label) == 2
            assert all(len(series) == end for series in per_label)
            assert len(payload['auroc'][base]) == end

    def test_the_first_series_is_the_first_label(self, payload, comparison):
        stats1, stats2, _ = comparison
        end = payload['endPosition']

        assert payload['labels'] == ['a', 'b']
        assert payload['counts'] == [len(stats1.sequences), len(stats2.sequences)]
        expected = round(float(stats1.stats['Per position nucleotide content']['A'][0]), 3)
        assert payload['freq']['A'][0][0] == expected

    def test_the_class_colors_are_the_ones_the_plots_use(self, payload):
        assert payload['colors'] == list(CLASS_COLORS)

    def test_coverage_is_a_proportion_of_all_sequences(self, payload):
        # Every sequence here is the same length, so the whole window is covered.
        assert payload['coverage'] == [1.0] * payload['endPosition']

    def test_it_round_trips_through_json(self, payload):
        assert json.loads(json.dumps(payload))['endPosition'] == payload['endPosition']

    def test_coverage_of_nothing_is_zero_rather_than_an_error(self):
        """A class with no sequences has no coverage; the curve is still drawn."""
        class Empty:
            stats = {'Sequence lengths': pd.DataFrame({'Sequence lengths': []})}

        assert _coverage(Empty(), Empty(), 4) == [0.0, 0.0, 0.0, 0.0]


class TestFlagsAndScores:
    def test_flag_positions_are_one_based_and_match_the_results(self, payload, comparison):
        _, _, results = comparison
        prefix = FEATURE_NAMES['forward']

        flagged = [(base, int(position)) for base in BASES for position in payload['flags'][base]]
        assert flagged, 'expected this comparison to flag something'
        for base, position in flagged:
            assert results.loc[f'{prefix} - {base} position {position}', 'Flag'] \
                == payload['flags'][base][str(position)]

    def test_passing_positions_are_left_out(self, payload, comparison):
        _, _, results = comparison
        prefix = FEATURE_NAMES['forward']

        for base in BASES:
            for position in range(1, payload['endPosition'] + 1):
                flag = results.loc[f'{prefix} - {base} position {position}', 'Flag']
                if flag == 'Pass':
                    assert str(position) not in payload['flags'][base]

    def test_the_scores_are_the_ones_from_the_results(self, payload, comparison):
        _, _, results = comparison
        prefix = FEATURE_NAMES['forward']

        for base in BASES:
            for position, score in enumerate(payload['auroc'][base], start=1):
                expected = results.loc[f'{prefix} - {base} position {position}', 'AU-ROC']
                if score is None:
                    assert not np.isfinite(expected)
                else:
                    assert score == pytest.approx(expected, abs=5e-4)

    def test_a_position_that_was_not_scored_is_unknown_and_carries_no_score(self):
        """Below the cohort floor a position is not measured, and must not look
        as if it had been: no AU-ROC, and Unknown rather than Pass."""
        stats1 = make_stats(random_sequences(4, 12, seed=3), label='a')
        stats2 = make_stats(random_sequences(4, 12, seed=4), label='b')
        results, _ = flag_significant_differences(stats1, stats2)
        end = min(stats1.end_position, stats2.end_position)

        payload = build_payload(stats1, stats2, BASES, end,
                                pd.DataFrame.from_dict(results, orient='index'), 'forward')

        assert payload['flags']['A'] == {str(position): 'Unknown'
                                        for position in range(1, end + 1)}
        assert payload['auroc']['A'] == [None] * end


class TestDirections:
    def test_the_reversed_payload_reads_the_reversed_statistics(self, comparison):
        stats1, stats2, results = comparison
        end = min(stats1.end_position, stats2.end_position)

        payload = build_payload(stats1, stats2, BASES, end, results, 'reversed')

        assert payload['direction'] == 'reversed'
        assert payload['xLabel'] == 'Position in reversed sequence'
        expected = round(float(stats1.stats['Per position reversed nucleotide content']['A'][0]), 3)
        assert payload['freq']['A'][0][0] == expected

    def test_an_unknown_direction_is_refused(self, comparison):
        stats1, stats2, results = comparison

        with pytest.raises(ValueError):
            build_payload(stats1, stats2, BASES, 10, results, 'sideways')

    @pytest.mark.parametrize('bases, end', [([], 30), (BASES, 0)])
    def test_nothing_to_draw_gives_no_payload(self, comparison, bases, end):
        stats1, stats2, results = comparison

        assert build_payload(stats1, stats2, bases, end, results, 'forward') is None


class TestMarkup:
    def test_the_payload_cannot_close_its_own_script_element(self):
        script = payload_script({'labels': ['</script><img src=x>', 'b & c']}, 'ppv-fwd-data')

        assert '</script>' == script[-len('</script>'):]
        assert '</script>' not in script[:-len('</script>')]
        assert '\\u003c' in script and '\\u0026' in script

    def test_the_figure_carries_everything_the_viewer_looks_for(self, payload):
        markup = viewer_html(payload, 'ppv-fwd')

        assert 'id="ppv-fwd-data"' in markup
        assert 'data-payload="ppv-fwd-data"' in markup
        for selector in ('ppv-bar', 'ppv-canvas', 'ppv-tooltip', 'ppv-readout',
                         'ppv-plot', 'ppv-flags-count', 'ppv-flags-body'):
            assert selector in markup, selector
        # and a line for readers without JavaScript, which the viewer removes
        assert 'ppv-fallback' in markup
