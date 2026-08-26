"""Tests for the minimum class size the per-sequence checks require.

Below `MIN_SEQUENCES_PER_CLASS` sequences in the smaller class, the worst case
over a check's features crosses the flag boundary on sampling noise often enough
that a flag carries no information, so these checks report Unknown instead of a
verdict. Two things must remain true alongside that: the integrity checks are
not statistical comparisons and are never suppressed by size, and everything
that is descriptive rather than inferential is still computed and plotted.
"""

import numpy as np
import pandas as pd
import pytest

from genomic_benchmarks_qc.report.classes_html_report import generate_not_scored_html
from genomic_benchmarks_qc.report.report_generator import generate_dataset_html_report
from genomic_benchmarks_qc.report.utils import icon_html
from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics
from genomic_benchmarks_qc.utils.testing import (
    MIN_SEQUENCES_PER_CLASS,
    _score_dataframe_features,
    _score_scalar_feature,
    flag_significant_differences,
)

# The checks that ask "does one feature separate the classes?", which is the
# question a small class cannot answer.
PER_SEQUENCE_CHECKS = [
    'Sequence lengths',
    'Per sequence GC content',
    'Per sequence nucleotide content',
    'Per sequence dinucleotide content',
]

# The checks that ask about the data as it stands, which any class size answers.
INTEGRITY_CHECKS = [
    'Unique bases',
    'Sequence Duplications within Labels',
    'Duplicate Sequences between Labels',
]


def random_sequences(count, length, seed, composition=None):
    rng = np.random.default_rng(seed)
    return [''.join(rng.choice(list('ACGT'), size=length, p=composition)) for _ in range(count)]


def make_stats(sequences, label='cls'):
    stats = SequenceStatistics(sequences=sequences, filename=f'{label}.fa',
                               filepath=f'/{label}.fa', label=label)
    stats.compute()
    return stats


def frame(values, column='feature'):
    return pd.DataFrame({column: np.asarray(values, dtype=float)})


class TestScalarFeatureFloor:
    def test_a_small_class_is_not_scored(self):
        size = MIN_SEQUENCES_PER_CLASS - 1
        indices = np.arange(size)

        metrics = _score_scalar_feature(frame(np.zeros(size)), frame(np.ones(size)),
                                        'feature', indices, indices)

        assert metrics['Flag'] == 'Unknown'
        assert np.isnan(metrics['AU-ROC'])

    def test_the_smaller_class_governs(self):
        """Real comparisons are often imbalanced; a large class cannot make up
        for a small one, because the noise comes from whichever side is thin."""
        small = np.arange(MIN_SEQUENCES_PER_CLASS - 1)
        large = np.arange(MIN_SEQUENCES_PER_CLASS * 10)

        metrics = _score_scalar_feature(frame(np.zeros(small.size)), frame(np.ones(large.size)),
                                        'feature', small, large)

        assert metrics['Flag'] == 'Unknown'

    def test_a_class_at_the_floor_is_scored(self):
        indices = np.arange(MIN_SEQUENCES_PER_CLASS)

        metrics = _score_scalar_feature(frame(np.zeros(indices.size)), frame(np.ones(indices.size)),
                                        'feature', indices, indices)

        assert metrics['Flag'] == 'Fail'
        assert metrics['AU-ROC'] == pytest.approx(1.0)

    def test_the_floor_can_be_lifted_explicitly(self):
        """The study that chose the floor has to be able to measure the scorer
        without it; passing 0 is how, and it keeps that path tested."""
        indices = np.arange(10)

        metrics = _score_scalar_feature(frame(np.zeros(10)), frame(np.ones(10)),
                                        'feature', indices, indices, min_class_size=0)

        assert metrics['Flag'] == 'Fail'


class TestDataFrameFeatureFloor:
    def test_every_column_and_the_aggregate_report_unknown_together(self):
        size = MIN_SEQUENCES_PER_CLASS - 1
        indices = np.arange(size)
        frame_1 = pd.DataFrame({'A': np.zeros(size), 'C': np.ones(size)})
        frame_2 = pd.DataFrame({'A': np.ones(size), 'C': np.zeros(size)})

        results = _score_dataframe_features(frame_1, frame_2, 'nt', indices, indices)

        assert results['nt']['Flag'] == 'Unknown'
        assert results['nt - A']['Flag'] == 'Unknown'
        assert results['nt - C']['Flag'] == 'Unknown'

    def test_a_class_at_the_floor_is_scored(self):
        size = MIN_SEQUENCES_PER_CLASS
        indices = np.arange(size)
        frame_1 = pd.DataFrame({'A': np.zeros(size), 'C': np.ones(size)})
        frame_2 = pd.DataFrame({'A': np.ones(size), 'C': np.zeros(size)})

        results = _score_dataframe_features(frame_1, frame_2, 'nt', indices, indices)

        assert results['nt']['Flag'] == 'Fail'


class TestSmallDatasetComparison:
    """End-to-end, on a dataset too small for the statistical checks."""

    @pytest.fixture
    def results(self):
        # Composition differs sharply, so these classes would flag loudly at any
        # size the tool is willing to score.
        stats1 = make_stats(random_sequences(40, 60, seed=1, composition=[0.7, 0.1, 0.1, 0.1]),
                            label='a')
        stats2 = make_stats(random_sequences(40, 60, seed=2, composition=[0.1, 0.1, 0.1, 0.7]),
                            label='b')
        results, _ = flag_significant_differences(stats1, stats2)
        return results

    @pytest.mark.parametrize('check', PER_SEQUENCE_CHECKS)
    def test_per_sequence_checks_report_unknown_not_pass(self, results, check):
        assert results[check]['Flag'] == 'Unknown'

    @pytest.mark.parametrize('check', INTEGRITY_CHECKS)
    def test_integrity_checks_still_report_a_verdict(self, results, check):
        assert results[check]['Flag'] in ('Pass', 'Warning', 'Fail')

    def test_a_warning_names_the_checks_that_were_not_scored(self, caplog):
        stats1 = make_stats(random_sequences(40, 60, seed=3), label='a')
        stats2 = make_stats(random_sequences(40, 60, seed=4), label='b')

        with caplog.at_level('WARNING'):
            flag_significant_differences(stats1, stats2)

        messages = ' '.join(record.message for record in caplog.records)
        assert 'Per sequence dinucleotide content' in messages
        assert str(MIN_SEQUENCES_PER_CLASS) in messages


class TestIntegrityChecksIgnoreSize:
    """A base or a duplicate that appears once is still worth reporting.

    These checks answer a question about the data as it stands, not a question
    about whether a model could exploit it, so they are deliberately exempt from
    the size floor - however small the class, and however rare the finding.
    """

    def test_one_stray_base_in_a_tiny_class_still_fails(self):
        stats1 = make_stats(['ACGT'] * 9 + ['ACGN'], label='a')
        stats2 = make_stats(['ACGT'] * 10, label='b')

        results, _ = flag_significant_differences(stats1, stats2)

        assert results['Unique bases']['Flag'] == 'Fail'

    def test_one_shared_sequence_in_a_tiny_class_still_fails(self):
        shared = 'ACGTACGTAC'
        stats1 = make_stats([shared] + random_sequences(9, 10, seed=5), label='a')
        stats2 = make_stats([shared] + random_sequences(9, 10, seed=6), label='b')

        results, _ = flag_significant_differences(stats1, stats2)

        assert results['Duplicate Sequences between Labels']['Flag'] == 'Fail'

    def test_one_duplicate_within_a_tiny_class_still_warns(self):
        repeated = random_sequences(9, 10, seed=7)
        stats1 = make_stats(repeated + [repeated[0]], label='a')
        stats2 = make_stats(random_sequences(10, 10, seed=8), label='b')

        results, _ = flag_significant_differences(stats1, stats2)

        assert results['Sequence Duplications within Labels']['Flag'] in ('Warning', 'Fail')


class TestDescriptiveStatisticsAreUnaffected:
    def test_a_small_class_still_computes_its_features(self):
        """Nothing descriptive is suppressed: the plots are drawn from these,
        and a reader with too few sequences to flag can still look."""
        stats = make_stats(random_sequences(10, 30, seed=9))

        assert stats.stats['Number of sequences'] == 10
        assert len(stats.stats['Per sequence GC content']) == 10
        assert len(stats.stats['Per sequence nucleotide content']) == 10
        assert not stats.stats['Per position nucleotide content'].empty


class TestReportSaysWhatWasNotScored:
    """A grey icon in a sidebar is easy to read as a pass; the report has to say
    otherwise in words, or the floor quietly turns small datasets into clean
    ones."""

    def test_unknown_gets_its_own_icon_rather_than_falling_through(self):
        markup = icon_html({'Sequence lengths': 'Unknown'}, 'Sequence lengths')

        assert 'status-unknown' in markup
        assert 'Unknown' not in markup  # the word itself must not leak into the page

    def test_the_note_names_the_checks_and_the_reason(self):
        stats1 = make_stats(random_sequences(40, 60, seed=10), label='a')
        stats2 = make_stats(random_sequences(40, 60, seed=11), label='b')

        note = generate_not_scored_html(stats1, stats2, {
            'Sequence lengths': 'Unknown',
            'Per sequence GC content': 'Unknown',
            'Unique bases': 'Pass',
        })

        assert 'Sequence lengths' in note
        assert 'Per sequence GC content' in note
        assert 'Unique bases' not in note
        assert str(MIN_SEQUENCES_PER_CLASS) in note

    def test_detail_rows_do_not_swamp_the_note(self):
        """Per-base and per-position rows carry their own Unknowns; listing them
        would bury the four names that matter."""
        stats1 = make_stats(random_sequences(40, 60, seed=12), label='a')
        stats2 = make_stats(random_sequences(40, 60, seed=13), label='b')

        note = generate_not_scored_html(stats1, stats2, {
            'Per sequence nucleotide content': 'Unknown',
            'Per sequence nucleotide content - A': 'Unknown',
            'Per sequence nucleotide content - C': 'Unknown',
        })

        assert '1 check(s) were not scored' in note

    def test_a_fully_scored_comparison_gets_no_note(self):
        stats1 = make_stats(random_sequences(10, 60, seed=14), label='a')
        stats2 = make_stats(random_sequences(10, 60, seed=15), label='b')

        assert generate_not_scored_html(stats1, stats2, {'Sequence lengths': 'Pass'}) == ''

    def test_the_note_reaches_the_rendered_report(self, tmp_path):
        stats1 = make_stats(random_sequences(40, 60, seed=16), label='a')
        stats2 = make_stats(random_sequences(40, 60, seed=17), label='b')
        results, failed_by_feature = flag_significant_differences(stats1, stats2)

        generate_dataset_html_report(
            stats1, stats2, tmp_path / 'gb-qc-report.html', plots_path=tmp_path / 'plots',
            plot_type='boxen',
            results=pd.DataFrame.from_dict(results, orient='index'),
            failed_by_feature=failed_by_feature,
        )

        page = (tmp_path / 'gb-qc-report.html').read_text()
        assert 'not-scored-note' in page
        assert 'were not scored' in page
        # Every plot is still drawn from all the data, floor or no floor - the
        # per-position ones included, over the reported window, since with no
        # compared window there is no narrower one to draw. What the floor
        # changes is the flags, not the figures.
        assert (tmp_path / 'plots' / 'sequence_lengths.png').exists()
        assert (tmp_path / 'plots' / 'per_position_nucleotide_content.png').exists()
        assert 'id="ppv-fwd"' in page
        assert 'No position could be compared' in page
