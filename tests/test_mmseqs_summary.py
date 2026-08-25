"""Tests for reducing an MMseqs2 hit table to what the split report needs.

The summariser is the one place that sees every hit, so it is also the one place
that must not hold them. Sequences are carried through the search by number and
come back as positions in an array, which is what these check: that the join is
right, that a sequence the search said nothing about stays distinguishable from
one whose best hit scored nothing, and that a table from somewhere else is
refused rather than quietly mis-joined.
"""

import numpy as np
import pandas as pd
import pytest

from genomic_benchmarks_qc.utils.mmseqs_summary import (
    MMSEQS_REQUIRED_COLS,
    sequence_id,
    staged_ids,
    summarize_mmseqs_output,
)


def write_hits(path, hits):
    """Write a hit table in the columns MMseqs2 is asked for.

    `hits` is (query index, target index, coverage, percent identity); the
    similarity the summariser scores by is their product.
    """
    rows = []
    for query, target, coverage, pident in hits:
        rows.append({
            'query': sequence_id(query, 'test'),
            'target': sequence_id(target, 'train'),
            'qcov': coverage, 'tcov': coverage, 'pident': pident,
            'evalue': 1e-9, 'qstart': 1, 'qend': 100, 'tstart': 1, 'tend': 100,
            'alnlen': 100, 'qaln': 'A' * 10, 'taln': 'A' * 10,
        })
    frame = pd.DataFrame(rows, columns=MMSEQS_REQUIRED_COLS)
    frame.to_csv(path, sep='\t', index=False)
    return path


class TestTheJoinBackToSequences:

    def test_a_sequence_with_no_hit_is_not_a_sequence_that_scored_zero(self, tmp_path):
        """NaN, not 0.0 - the histogram counts the two differently."""
        hits = write_hits(tmp_path / 'hits.tsv', [(2, 0, 1.0, 95.0)])

        summary = summarize_mmseqs_output(hits, 90.0, query_count=4, target_count=3)

        maxima = summary['query_similarity_max']
        assert maxima[2] == pytest.approx(95.0)
        assert np.isnan(maxima[[0, 1, 3]]).all()

    def test_the_maximum_is_carried_across_chunks(self, tmp_path):
        """One sequence's hits can fall either side of a chunk boundary."""
        hits = write_hits(tmp_path / 'hits.tsv', [
            (0, 0, 1.0, 70.0), (1, 1, 1.0, 60.0),
            (0, 2, 1.0, 95.0), (1, 2, 1.0, 55.0),
        ])

        summary = summarize_mmseqs_output(
            hits, 90.0, query_count=2, target_count=3, chunksize=2)

        assert summary['query_similarity_max'][0] == pytest.approx(95.0)
        assert summary['query_similarity_max'][1] == pytest.approx(60.0)
        assert list(summary['query_above_threshold']) == [True, False]

    def test_only_the_leaked_sequences_are_flagged(self, tmp_path):
        hits = write_hits(tmp_path / 'hits.tsv', [
            (0, 0, 1.0, 95.0),   # leaked
            (1, 1, 0.5, 95.0),   # 47.5, half the sequence covered
            (2, 2, 1.0, 89.9),   # just under
        ])

        summary = summarize_mmseqs_output(hits, 90.0, query_count=3, target_count=3)

        assert list(summary['query_above_threshold']) == [True, False, False]
        assert list(summary['target_above_threshold']) == [True, False, False]
        assert summary['leaked_hits'] == 1
        assert summary['total_hits'] == 3

    def test_an_empty_table_leaves_every_sequence_unhit(self, tmp_path):
        hits = write_hits(tmp_path / 'hits.tsv', [])

        summary = summarize_mmseqs_output(hits, 90.0, query_count=3, target_count=2)

        assert np.isnan(summary['query_similarity_max']).all()
        assert not summary['query_above_threshold'].any()
        assert summary['total_hits'] == 0

    def test_a_table_from_somewhere_else_is_refused(self, tmp_path):
        """A name this run did not stage cannot be joined to anything.

        Silently skipping it would under-report leakage, which is the direction
        that matters: a split would read clean because the hits were dropped.
        """
        path = tmp_path / 'hits.tsv'
        write_hits(path, [(0, 0, 1.0, 95.0)])
        frame = pd.read_csv(path, sep='\t')
        frame.loc[0, 'query'] = 'chr7:120000-120500'
        frame.to_csv(path, sep='\t', index=False)

        with pytest.raises(RuntimeError, match='did not.*stage'):
            summarize_mmseqs_output(path, 90.0, query_count=1, target_count=1)

    def test_a_sequence_number_past_the_staged_count_is_refused(self, tmp_path):
        hits = write_hits(tmp_path / 'hits.tsv', [(9, 0, 1.0, 95.0)])

        with pytest.raises(RuntimeError, match='of 3 staged'):
            summarize_mmseqs_output(hits, 90.0, query_count=3, target_count=3)


class TestTheOrderOfTheTopHits:
    """The alignment table is a listing, so its order has to be reproducible.

    Every hit here scores the same, which is not a corner case: `min_cov*pident`
    is a product of two rounded percentages, so a real search returns ties by
    the hundred and the page has to put them somewhere.
    """

    def make_ties(self, tmp_path, count=40):
        """`count` hits that all score exactly 95.0, in a known reading order."""
        return write_hits(tmp_path / 'hits.tsv',
                          [(i, i, 1.0, 95.0) for i in range(count)])

    def test_ties_keep_the_order_they_were_read_in(self, tmp_path):
        hits = self.make_ties(tmp_path)

        summary = summarize_mmseqs_output(hits, 90.0, query_count=40, target_count=40)

        assert list(summary['results_filt']['query']) == [
            sequence_id(i, 'test') for i in range(40)]

    def test_the_same_table_summarises_to_the_same_order_twice(self, tmp_path):
        """The listing a reader cites has to still say that tomorrow."""
        hits = self.make_ties(tmp_path)

        first = summarize_mmseqs_output(hits, 90.0, query_count=40, target_count=40)
        second = summarize_mmseqs_output(hits, 90.0, query_count=40, target_count=40)

        assert list(first['results_filt']['query']) == list(second['results_filt']['query'])

    def test_more_similar_hits_still_come_first(self, tmp_path):
        """Dropping the second sort must not drop the sorting."""
        hits = write_hits(tmp_path / 'hits.tsv', [
            (0, 0, 1.0, 91.0), (1, 1, 1.0, 99.0), (2, 2, 1.0, 95.0),
        ])

        summary = summarize_mmseqs_output(hits, 90.0, query_count=3, target_count=3)

        assert list(summary['results_filt']['query']) == [
            sequence_id(1, 'test'), sequence_id(2, 'test'), sequence_id(0, 'test')]


class TestStagedIds:

    def test_the_ids_come_back_from_the_mask(self):
        mask = np.array([False, True, False, True])

        assert staged_ids(mask, 'test') == {'seq_1_test', 'seq_3_test'}
        assert staged_ids(np.zeros(4, dtype=bool), 'train') == set()

    def test_they_are_the_names_the_sequences_were_staged_under(self):
        """The two halves of the round trip, in one assertion."""
        mask = np.zeros(5, dtype=bool)
        mask[[0, 4]] = True

        assert staged_ids(mask, 'train') == {sequence_id(0, 'train'), sequence_id(4, 'train')}
