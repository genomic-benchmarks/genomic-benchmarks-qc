"""Tests for reducing an MMseqs2 hit table to what the split report needs.

The summariser is the one place that sees every hit, so it is also the one place
that must not hold them. Sequences are carried through the search by number and
come back as positions in an array, which is what these check: that the join is
right, that a sequence the search said nothing about stays distinguishable from
one whose best hit scored nothing, and that a table from somewhere else is
refused rather than quietly mis-joined.
"""

import io
import logging
import random

import numpy as np
import pandas as pd
import pytest

from genomic_benchmarks_qc.utils.mmseqs_summary import (
    MMSEQS_REQUIRED_COLS,
    MMSEQS_RESULT_COLUMNS,
    log_reversed_hit_warning,
    sequence_id,
    staged_ids,
    summarize_mmseqs_output,
)


def write_hits(path, hits, backwards=()):
    """Write a hit table in the columns MMseqs2 is asked for.

    `hits` is (query index, target index, coverage, percent identity); the
    similarity the summariser scores by is their product. `backwards` is the
    positions in `hits` whose target coordinates should run the wrong way,
    which is what a build reporting a hit on a strand nobody asked for does.
    """
    rows = []
    for position, (query, target, coverage, pident) in enumerate(hits):
        tstart, tend = (100, 1) if position in backwards else (1, 100)
        rows.append({
            'query': sequence_id(query, 'test'),
            'target': sequence_id(target, 'train'),
            'qcov': coverage, 'tcov': coverage, 'pident': pident,
            'evalue': 1e-9, 'qstart': 1, 'qend': 100, 'tstart': tstart, 'tend': tend,
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

    Most hits here score the same, which is not a corner case: `min_cov*pident`
    is a product of two rounded percentages, so a real search returns ties by
    the hundred and the page has to put them somewhere. Nor can it put them
    where they were read: MMseqs2 on more than one thread writes the same hits
    in a different order from one run to the next.
    """

    def make_ties(self, tmp_path, order, name='hits.tsv'):
        """Hits that all score exactly 95.0, read in `order` of their query."""
        return write_hits(tmp_path / name, [(i, i, 1.0, 95.0) for i in order])

    def test_ties_come_in_staging_order_whatever_order_they_were_read_in(self, tmp_path):
        order = list(range(40))
        random.Random(0).shuffle(order)
        hits = self.make_ties(tmp_path, order)

        summary = summarize_mmseqs_output(hits, 90.0, query_count=40, target_count=40)

        assert list(summary['results_filt']['query']) == [
            sequence_id(i, 'test') for i in range(40)]

    def test_staging_order_is_numeric_not_alphabetical(self, tmp_path):
        """seq_9 is staged before seq_10, and 'seq_10' < 'seq_9' as strings."""
        hits = self.make_ties(tmp_path, [10, 9])

        summary = summarize_mmseqs_output(hits, 90.0, query_count=11, target_count=11)

        assert list(summary['results_filt']['query']) == [
            sequence_id(9, 'test'), sequence_id(10, 'test')]

    def test_a_query_with_several_hits_lists_them_by_train_sequence(self, tmp_path):
        hits = write_hits(tmp_path / 'hits.tsv',
                          [(0, 7, 1.0, 95.0), (0, 2, 1.0, 95.0), (0, 5, 1.0, 95.0)])

        summary = summarize_mmseqs_output(hits, 90.0, query_count=1, target_count=8)

        assert list(summary['results_filt']['target']) == [
            sequence_id(i, 'train') for i in (2, 5, 7)]

    def test_the_cap_keeps_the_same_hits_however_they_arrive(self, tmp_path):
        """Which ties make the cut is part of the order too, chunk by chunk."""
        shuffled = list(range(40))
        random.Random(1).shuffle(shuffled)
        in_order = summarize_mmseqs_output(
            self.make_ties(tmp_path, range(40), 'a.tsv'), 90.0,
            query_count=40, target_count=40, top_n=10)
        out_of_order = summarize_mmseqs_output(
            self.make_ties(tmp_path, shuffled, 'b.tsv'), 90.0,
            query_count=40, target_count=40, top_n=10, chunksize=7)

        assert list(in_order['results_filt']['query']) == [
            sequence_id(i, 'test') for i in range(10)]
        pd.testing.assert_frame_equal(in_order['results_filt'], out_of_order['results_filt'])

    def test_the_rows_are_numbered_in_listing_order(self, tmp_path):
        """The report names each alignment after its row's position."""
        hits = self.make_ties(tmp_path, [3, 1, 2, 0])

        summary = summarize_mmseqs_output(hits, 90.0, query_count=4, target_count=4)

        assert list(summary['results_filt'].index) == [0, 1, 2, 3]

    def test_more_similar_hits_still_come_first(self, tmp_path):
        hits = write_hits(tmp_path / 'hits.tsv', [
            (0, 0, 1.0, 91.0), (1, 1, 1.0, 99.0), (2, 2, 1.0, 95.0),
        ])

        summary = summarize_mmseqs_output(hits, 90.0, query_count=3, target_count=3)

        assert list(summary['results_filt']['query']) == [
            sequence_id(1, 'test'), sequence_id(2, 'test'), sequence_id(0, 'test')]


class TestTheExportedHits:
    """Every leaked hit goes to a file, in the listing's order.

    So the same search always writes the same file, and the pairs the page
    lists are its first rows. The file is written without holding every hit:
    past one chunk's worth they are sorted in runs and merged.
    """

    def hits(self, count=30):
        """Leaked hits of varied similarity, several per query, in a mixed order."""
        rows = [(i % 7, i, 1.0, 90.0 + i % 4) for i in range(count)]
        random.Random(2).shuffle(rows)
        return rows

    def export(self, tmp_path, rows, name, **kwargs):
        path = write_hits(tmp_path / f'{name}.tsv', rows)
        export_path = tmp_path / f'{name}-export.tsv'
        summary = summarize_mmseqs_output(path, 90.0, query_count=7, target_count=30,
                                          export_path=export_path, **kwargs)
        return summary, export_path.read_text()

    def test_the_file_does_not_depend_on_the_order_the_hits_came_in(self, tmp_path):
        rows = self.hits()
        _, forward = self.export(tmp_path, rows, 'forward')
        _, backward = self.export(tmp_path, rows[::-1], 'backward')

        assert forward == backward

    def test_the_file_does_not_depend_on_how_it_was_read(self, tmp_path):
        """Chunks of four make several sorted runs to merge; one chunk makes none."""
        rows = self.hits()
        _, whole = self.export(tmp_path, rows, 'whole')
        _, merged = self.export(tmp_path, rows, 'merged', chunksize=4)

        assert merged == whole

    def test_the_listing_is_the_top_of_the_file(self, tmp_path):
        summary, text = self.export(tmp_path, self.hits(), 'hits', top_n=5, chunksize=4)

        exported = pd.read_csv(io.StringIO(text), sep='\t')
        assert list(exported['query'][:5]) == list(summary['results_filt']['query'])
        assert list(exported['target'][:5]) == list(summary['results_filt']['target'])
        assert len(exported) == summary['leaked_hits']

    def test_the_sorted_runs_are_cleaned_up(self, tmp_path):
        self.export(tmp_path, self.hits(), 'hits', chunksize=4)

        assert not list(tmp_path.glob('gb-qc-export-*'))

    def test_a_clean_split_still_gets_the_header(self, tmp_path):
        _, text = self.export(tmp_path, [(0, 0, 1.0, 40.0)], 'clean')

        assert text.splitlines() == ['\t'.join(MMSEQS_RESULT_COLUMNS)]


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


class TestBackwardsAlignments:
    """Hits reported on a strand the search never asked for.

    The search is run forward-strand only, so a hit whose coordinates descend
    cannot be drawn against the sequences it names. Some MMseqs2 builds emit
    them anyway. They are counted rather than dropped: what they corrupt is the
    alignment, not the scores, so the pair still belongs in the leakage numbers.
    """

    def test_a_clean_table_counts_none(self, tmp_path):
        path = write_hits(tmp_path / 'hits.tsv', [(0, 0, 1.0, 100.0), (1, 1, 1.0, 95.0)])
        summary = summarize_mmseqs_output(path, 90.0, query_count=2, target_count=2)
        assert summary['reversed_hits'] == 0
        assert summary['reversed_leaked_hits'] == 0

    def test_backwards_coordinates_are_counted(self, tmp_path):
        path = write_hits(
            tmp_path / 'hits.tsv',
            [(0, 0, 1.0, 100.0), (1, 1, 1.0, 95.0), (2, 2, 1.0, 99.0)],
            backwards=(0, 2),
        )
        summary = summarize_mmseqs_output(path, 90.0, query_count=3, target_count=3)
        assert summary['reversed_hits'] == 2

    def test_only_the_ones_that_reach_the_report_are_counted_separately(self, tmp_path):
        # 1.0 * 40.0 is below the threshold, so that backwards hit never reaches
        # the listing; the other one does.
        path = write_hits(
            tmp_path / 'hits.tsv',
            [(0, 0, 1.0, 100.0), (1, 1, 1.0, 40.0)],
            backwards=(0, 1),
        )
        summary = summarize_mmseqs_output(path, 90.0, query_count=2, target_count=2)
        assert summary['reversed_hits'] == 2
        assert summary['reversed_leaked_hits'] == 1

    def test_a_backwards_hit_still_counts_as_a_leak(self, tmp_path):
        # The scores on these rows match a good build exactly, so dropping them
        # would under-report the leakage the command exists to measure.
        straight = write_hits(tmp_path / 'a.tsv', [(0, 0, 1.0, 100.0)])
        backwards = write_hits(tmp_path / 'b.tsv', [(0, 0, 1.0, 100.0)], backwards=(0,))
        clean = summarize_mmseqs_output(straight, 90.0, query_count=1, target_count=1)
        dirty = summarize_mmseqs_output(backwards, 90.0, query_count=1, target_count=1)
        assert dirty['leaked_hits'] == clean['leaked_hits'] == 1
        assert dirty['query_above_threshold'].tolist() == [True]

    def test_they_are_counted_even_when_nothing_leaked(self, tmp_path):
        # The early exit for a chunk with no leaks used to be taken before the
        # count, which hid exactly the case worth warning about: a build that is
        # returning nonsense on a split that happens to be clean.
        path = write_hits(tmp_path / 'hits.tsv', [(0, 0, 1.0, 40.0)], backwards=(0,))
        summary = summarize_mmseqs_output(path, 90.0, query_count=1, target_count=1)
        assert summary['leaked_hits'] == 0
        assert summary['reversed_hits'] == 1

    def test_they_are_counted_across_chunk_boundaries(self, tmp_path):
        hits = [(i, i, 1.0, 100.0) for i in range(6)]
        path = write_hits(tmp_path / 'hits.tsv', hits, backwards=(0, 3, 5))
        summary = summarize_mmseqs_output(
            path, 90.0, query_count=6, target_count=6, chunksize=2)
        assert summary['reversed_hits'] == 3


class TestSayingSoToTheUser:

    def test_nothing_is_said_when_there_are_none(self, caplog):
        with caplog.at_level(logging.WARNING, logger='genomic_benchmarks_qc'):
            log_reversed_hit_warning(0, 0, 100, threads=4)
        assert caplog.records == []

    def test_the_counts_and_both_remedies_are_named(self, caplog):
        with caplog.at_level(logging.WARNING, logger='genomic_benchmarks_qc'):
            log_reversed_hit_warning(12, 3, 400, threads=8)
        message = caplog.text
        assert '12 of 400' in message
        assert '--threads 1' in message
        assert 'precompiled' in message
        assert '--threads 8' in message

    def test_a_report_that_got_away_with_it_is_told_so(self, caplog):
        # Worth separating: a user whose report is fine still wants to know the
        # build is returning nonsense, but should not go looking for damage.
        with caplog.at_level(logging.WARNING, logger='genomic_benchmarks_qc'):
            log_reversed_hit_warning(12, 0, 400, threads=2)
        assert 'this report is unaffected' in caplog.text

    def test_a_report_that_did_not_says_what_is_missing(self, caplog):
        with caplog.at_level(logging.WARNING, logger='genomic_benchmarks_qc'):
            log_reversed_hit_warning(12, 3, 400, threads=2)
        assert 'without an alignment' in caplog.text
        assert 'count as leaks' in caplog.text
        assert 'this report is unaffected' not in caplog.text

    def test_single_threaded_runs_are_not_told_to_use_one_thread_they_already_use(self, caplog):
        with caplog.at_level(logging.WARNING, logger='genomic_benchmarks_qc'):
            log_reversed_hit_warning(12, 3, 400, threads=1)
        assert 'this run used' not in caplog.text
