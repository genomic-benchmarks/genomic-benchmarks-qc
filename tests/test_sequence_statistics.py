"""Tests for the composition statistics a class is described by.

Every feature here comes out of one vectorised pass over the class: the
sequences are joined, read as an array of code points and turned into one
column index per character, and the four tables are `bincount`s over that.
The pass works a block of sequences at a time so its working set does not grow
with the class, which is the part worth pinning - a block boundary in the wrong
place is a silent off-by-one in the middle of a frequency, not an error.
"""

import numpy as np
import pytest

from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics

# Three sequences that between them cover the cases the pass has to get right:
# one full-length, one that stops early, and one of a single repeated base.
EXAMPLE = ['ACGT', 'AC', 'AAAA']


def make_stats(sequences, label='cls'):
    stats = SequenceStatistics(sequences=sequences, filename=f'{label}.fa',
                               filepath=f'/{label}.fa', label=label)
    stats.compute()
    return stats


@pytest.fixture
def example():
    return make_stats(EXAMPLE)


class TestPerSequenceComposition:
    def test_each_sequence_is_a_fraction_of_its_own_length(self, example):
        frame = example.stats['Per sequence nucleotide content']

        assert list(frame.columns) == ['A', 'C', 'G', 'T']
        assert frame.loc[0].tolist() == [0.25, 0.25, 0.25, 0.25]
        assert frame.loc[1].tolist() == [0.5, 0.5, 0.0, 0.0]
        assert frame.loc[2].tolist() == [1.0, 0.0, 0.0, 0.0]

    def test_dinucleotides_are_a_fraction_of_the_pairs_a_sequence_has(self, example):
        frame = example.stats['Per sequence dinucleotide content']

        # One fewer pair than the sequence has characters, so 'AC' is one pair
        # and all of them, and 'AAAA' is three of the same one.
        assert frame.loc[0, 'AC'] == pytest.approx(1 / 3)
        assert frame.loc[0, 'CG'] == pytest.approx(1 / 3)
        assert frame.loc[0, 'GT'] == pytest.approx(1 / 3)
        assert frame.loc[1, 'AC'] == 1.0
        assert frame.loc[2, 'AA'] == 1.0

    def test_a_pair_never_straddles_two_sequences(self, example):
        """The last character of one sequence and the first of the next sit
        beside each other in the joined text, and are not a dinucleotide."""
        frame = example.stats['Per sequence dinucleotide content']

        # 'ACGT' then 'AC' would otherwise make a TA, and 'AC' then 'AAAA' a CA.
        assert frame['TA'].sum() == 0
        assert frame['CA'].sum() == 0

    def test_gc_content_is_a_percentage(self, example):
        gc = example.stats['Per sequence GC content']['Per sequence GC content']

        assert gc.tolist() == [50.0, 50.0, 0.0]

    def test_lengths_are_recorded_per_sequence(self, example):
        lengths = example.stats['Sequence lengths']['Sequence lengths']

        assert lengths.tolist() == [4.0, 2.0, 4.0]


class TestPerPositionComposition:
    def test_a_position_is_a_fraction_of_the_sequences_that_reach_it(self, example):
        frame = example.stats['Per position nucleotide content']

        assert list(frame.index) == [0, 1, 2, 3]
        assert frame.loc[0].tolist() == [1.0, 0.0, 0.0, 0.0]
        assert frame.loc[1, 'C'] == pytest.approx(2 / 3)
        # Only two sequences reach position 2, so each base there is a half.
        assert frame.loc[2, 'A'] == 0.5 and frame.loc[2, 'G'] == 0.5
        assert frame.loc[3, 'A'] == 0.5 and frame.loc[3, 'T'] == 0.5

    def test_the_reversed_feature_counts_from_the_far_end(self, example):
        frame = example.stats['Per position reversed nucleotide content']

        # Position 0 is the last character of each sequence: T, C, A.
        assert frame.loc[0].tolist() == pytest.approx([1 / 3, 1 / 3, 0.0, 1 / 3])
        assert frame.loc[3, 'A'] == 1.0

    def test_every_position_sums_to_one(self, example):
        for name in ('Per position nucleotide content',
                     'Per position reversed nucleotide content'):
            totals = example.stats[name].sum(axis=1)

            assert totals.tolist() == pytest.approx([1.0] * len(totals)), name

    def test_the_columns_are_the_class_bases_in_order(self):
        """Not the order the bases happen to turn up in, which is what the
        dictionaries this replaced left them in - a per-position frame whose
        columns were sorted differently from every other frame in the report,
        by an accident of the data."""
        stats = make_stats(['TTTT', 'GGGA', 'CCCC', 'AAAA'])

        assert (list(stats.stats['Per position nucleotide content'].columns)
                == stats.stats['Unique bases'] == ['A', 'C', 'G', 'T'])


class TestEdges:
    def test_an_empty_sequence_reads_as_zero_rather_than_a_division_by_zero(self):
        stats = make_stats(['ACGT', '', 'TTTT'])

        assert stats.stats['Per sequence nucleotide content'].loc[1].tolist() == [0.0] * 4
        assert stats.stats['Per sequence dinucleotide content'].loc[1].sum() == 0.0
        assert stats.stats['Per sequence GC content'].iloc[1, 0] == 0.0
        # and it contributes to no position, so the ones it does not reach are
        # still fractions of the sequences that do
        assert stats.stats['Per position nucleotide content'].loc[0, 'A'] == 0.5

    def test_a_class_of_nothing_computes(self):
        stats = make_stats([])

        assert stats.stats['Unique bases'] == []
        assert stats.stats['Per position nucleotide content'].empty
        assert stats.stats['Sequence lengths'].empty

    def test_a_base_outside_ascii_still_lands_in_its_own_column(self):
        """Bases are read as code points rather than bytes, so nothing about the
        pass assumes the alphabet is ACGT, or even Latin."""
        stats = make_stats(['ACÄ', 'AÄÄ'])

        assert stats.stats['Unique bases'] == ['A', 'C', 'Ä']
        assert stats.stats['Per sequence nucleotide content'].loc[1, 'Ä'] == pytest.approx(2 / 3)
        assert stats.stats['Per position nucleotide content'].loc[2, 'Ä'] == 1.0

    @pytest.mark.parametrize('characters, cells', [(1, 1), (3, 2), (7, 1000)])
    def test_the_numbers_do_not_depend_on_where_the_blocks_fall(
        self, monkeypatch, characters, cells
    ):
        """A block boundary is an implementation detail of the working set, and
        the only thing that keeps it one is that the counts are accumulated
        across blocks rather than computed within them."""
        rng = np.random.default_rng(4)
        sequences = [''.join(rng.choice(list('ACGT'), size=rng.integers(0, 9)))
                     for _ in range(25)]
        whole = make_stats(sequences)

        monkeypatch.setattr('genomic_benchmarks_qc.utils.seq_stats.STATS_BLOCK_CHARACTERS',
                            characters)
        monkeypatch.setattr('genomic_benchmarks_qc.utils.seq_stats.STATS_BLOCK_CELLS', cells)
        blocked = make_stats(sequences)

        for name, frame in whole.stats.items():
            if hasattr(frame, 'to_numpy'):
                assert np.allclose(frame.to_numpy(dtype=float),
                                   blocked.stats[name].to_numpy(dtype=float)), name
