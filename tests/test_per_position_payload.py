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
    drawn_window,
    payload_script,
    viewer_html,
)
from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics, cohort_floor
from genomic_benchmarks_qc.utils.testing import flag_significant_differences, position_windows

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
    stats1 = make_stats(random_sequences(300, 30, seed=1), label='a')
    stats2 = make_stats(random_sequences(300, 30, seed=2, composition=[.7, .1, .1, .1]), label='b')
    results, _ = flag_significant_differences(stats1, stats2)
    return stats1, stats2, pd.DataFrame.from_dict(results, orient='index')


@pytest.fixture(scope='module')
def payload(comparison):
    stats1, stats2, results = comparison
    return build_payload(stats1, stats2, BASES, drawn_window(stats1, stats2),
                         results, 'forward')


class TestShape:
    def test_the_series_are_as_long_as_the_window(self, payload):
        end = payload['endPosition']

        assert len(payload['coverage']) == 2
        assert all(len(series) == end for series in payload['coverage'])
        for base in BASES:
            per_label = payload['freq'][base]
            assert len(per_label) == 2
            assert all(len(series) == end for series in per_label)

    def test_the_first_series_is_the_first_label(self, payload, comparison):
        stats1, stats2, _ = comparison
        end = payload['endPosition']

        assert payload['labels'] == ['a', 'b']
        assert payload['counts'] == [len(stats1.sequences), len(stats2.sequences)]
        expected = round(float(stats1.stats['Per position nucleotide content']['A'][0]), 3)
        assert payload['freq']['A'][0][0] == expected

    def test_the_class_colors_are_the_ones_the_plots_use(self, payload):
        assert payload['colors'] == list(CLASS_COLORS)

    def test_coverage_is_one_curve_per_class(self, payload):
        # Every sequence here is the same length, so both classes cover the
        # whole window.
        assert payload['coverage'] == [[1.0] * payload['endPosition']] * 2

    def test_the_cohort_floor_travels_with_the_curves(self, payload, comparison):
        """The panel draws the floor as a bare dashed line, so the number travels
        with the curves and nothing else does: naming it is the explanation's job."""
        stats1, stats2, _ = comparison

        assert payload['coverageFloor'] == pytest.approx(cohort_floor(stats1, stats2))
        assert 'coverageFloorLabel' not in payload

    def test_it_round_trips_through_json(self, payload):
        assert json.loads(json.dumps(payload))['endPosition'] == payload['endPosition']

    def test_coverage_of_nothing_is_zero_rather_than_an_error(self):
        """A class with no sequences has no coverage; the curve is still drawn."""
        class Empty:
            stats = {'Sequence lengths': pd.DataFrame({'Sequence lengths': []})}

        assert _coverage(Empty(), 4) == [0.0, 0.0, 0.0, 0.0]


class TestFlags:
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

    def test_a_missing_check_is_unknown_rather_than_a_silent_pass(self, comparison):
        """The viewer reads a position with no entry as a Pass, so a check the
        results table does not name has to arrive as Unknown instead."""
        stats1, stats2, results = comparison
        _, scored_end = position_windows(stats1, stats2)
        prefix = FEATURE_NAMES['forward']
        dropped = f'{prefix} - A position 2'
        assert results.loc[dropped, 'Flag'] != 'Unknown'

        payload = build_payload(stats1, stats2, BASES, scored_end,
                                results.drop(index=dropped), 'forward')

        assert payload['flags']['A']['2'] == 'Unknown'

    def test_no_separability_score_travels_with_the_flags(self, payload):
        """The report shows the flags and the frequencies behind them. For a
        per-position base the AU-ROC is a restatement of the gap between those
        two frequencies, so it stays in report.csv and out of the page."""
        assert 'auroc' not in payload

    def test_pass_has_a_colour_even_though_it_is_never_a_band(self, payload):
        """The hover card names the flag on every row, Pass included."""
        assert set(payload['flagColors']) == {'Fail', 'Warning', 'Pass', 'Unknown'}

    def test_a_comparison_that_scored_nothing_still_gets_a_figure(self):
        """An underpowered comparison keeps its plots everywhere else in the
        report and flags Unknown rather than Pass. The per-position panels do the
        same: with no compared window they fall back to the reported one, and
        every position in it arrives Unknown, so nothing in the figure can be
        read as having passed."""
        stats1 = make_stats(random_sequences(60, 12, seed=3), label='a')
        stats2 = make_stats(random_sequences(60, 12, seed=4), label='b')
        results, _ = flag_significant_differences(stats1, stats2)
        end, scored_end = position_windows(stats1, stats2)
        assert scored_end == 0 and end > 0
        assert drawn_window(stats1, stats2) == end

        payload = build_payload(stats1, stats2, BASES, drawn_window(stats1, stats2),
                                pd.DataFrame.from_dict(results, orient='index'), 'forward')

        assert payload['endPosition'] == end
        for base in BASES:
            assert set(payload['flags'][base]) == {str(p) for p in range(1, end + 1)}
            assert set(payload['flags'][base].values()) == {'Unknown'}


@pytest.fixture(scope='module')
def tail_payload():
    """A comparison whose sequences outrun its scored window."""
    # Three quarters of one class stop at 20, leaving too few beyond it to
    # compare, while the sequences themselves run on to 60.
    long_reads = random_sequences(400, 60, seed=5)
    stats1 = make_stats([seq[:20] if i % 4 else seq for i, seq in enumerate(long_reads)], label='a')
    stats2 = make_stats(random_sequences(400, 60, seed=6), label='b')
    results, _ = flag_significant_differences(stats1, stats2)
    end, scored_end = position_windows(stats1, stats2)
    payload = build_payload(stats1, stats2, BASES, drawn_window(stats1, stats2),
                            pd.DataFrame.from_dict(results, orient='index'), 'forward')
    return payload, end, scored_end


class TestTheComparedWindow:
    """The figure stops where the comparison did, however far the sequences run."""

    def test_the_window_ends_at_the_last_compared_position(self, tail_payload):
        payload, end, scored_end = tail_payload

        assert scored_end < end
        assert payload['endPosition'] == scored_end

    def test_the_curves_stop_there_too(self, tail_payload):
        payload, _, scored_end = tail_payload

        assert all(len(series) == scored_end for series in payload['coverage'])
        assert all(len(series) == scored_end for series in payload['freq']['A'])

    def test_nothing_in_the_window_is_unknown(self, tail_payload):
        """Every position drawn was compared, so the viewer never has to explain
        one away - that is what moving the boundary bought."""
        payload, _, _ = tail_payload

        for base in BASES:
            assert 'Unknown' not in payload['flags'][base].values()


class TestDirections:
    def test_the_reversed_payload_reads_the_reversed_statistics(self, comparison):
        stats1, stats2, results = comparison
        _, scored_end = position_windows(stats1, stats2)

        payload = build_payload(stats1, stats2, BASES, scored_end, results, 'reversed')

        assert payload['direction'] == 'reversed'
        # The axis says which end position 1 is, not that the sequence was
        # reversed to get there; the viewer's flag table reuses the same wording.
        assert payload['xLabel'] == 'Position from sequence end'
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

    def test_there_is_no_compared_region_button(self, payload):
        """The whole figure is the compared region, so zooming to it is what
        Reset zoom already does."""
        assert 'data-action="compared"' not in viewer_html(payload, 'ppv-fwd')
