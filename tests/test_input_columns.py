"""Tests for what a report says when a column name is wrong.

Mistyping a column is the ordinary way to get this command wrong: the defaults
are 'sequence' and 'label', and datasets in the wild call them 'seq', 'seqs',
'Sequence', 'x'. What came back before was pandas' own message, which names the
column that is missing and nothing else - not the file it looked in, not the
option the name came from, and not the columns that are actually there. Every
one of those is in the message now, and these keep them there.
"""

import gzip

import pytest

from genomic_benchmarks_qc.utils.input_utils import read_csv_file, stream_table_sequences


@pytest.fixture
def table(tmp_path):
    path = tmp_path / 'dataset.csv'
    path.write_text('seq,class\nACGT,0\nTTTT,1\n')
    return path


class TestTheMessageNamesWhatIsNeededToFixIt:

    def test_it_names_the_file_the_column_and_the_alternatives(self, table):
        with pytest.raises(ValueError) as excinfo:
            read_csv_file(table, 'csv', ['sequence'])
        message = str(excinfo.value)

        assert 'sequence' in message
        assert str(table) in message
        assert "'seq'" in message and "'class'" in message

    def test_it_names_the_option_the_name_came_from(self, table):
        """Two options read column names, and they fail identically otherwise."""
        with pytest.raises(ValueError, match='--sequence-column'):
            read_csv_file(table, 'csv', ['nope'], label_column='class')
        with pytest.raises(ValueError, match='--label-column'):
            read_csv_file(table, 'csv', ['seq'], label_column='label')

    def test_both_wrong_names_are_reported_together(self, table):
        """One run, one message: fixing the first should not reveal the second."""
        with pytest.raises(ValueError) as excinfo:
            read_csv_file(table, 'csv', ['nope'], label_column='label')
        message = str(excinfo.value)

        assert "'nope' (--sequence-column)" in message
        assert "'label' (--label-column)" in message

    def test_the_streaming_reader_says_the_same_thing(self, table):
        """evaluate-splits never loads the table, and gets the same message."""
        with pytest.raises(ValueError, match='--sequence-column'):
            list(stream_table_sequences(table, 'csv', ['nope']))

    def test_a_gzipped_table_is_read_far_enough_to_answer(self, tmp_path):
        path = tmp_path / 'dataset.csv.gz'
        with gzip.open(path, 'wt') as handle:
            handle.write('seq,class\nACGT,0\n')

        with pytest.raises(ValueError, match="'seq'"):
            read_csv_file(path, 'csv.gz', ['sequence'])

    def test_a_tsv_is_split_on_tabs_before_the_header_is_judged(self, tmp_path):
        """On the wrong delimiter every column looks missing, including the right one."""
        path = tmp_path / 'dataset.tsv'
        path.write_text('seq\tclass\nACGT\t0\n')

        read_csv_file(path, 'tsv', ['seq'], label_column='class')


class TestTheColumnsThatAreThereStillWork:

    def test_naming_them_correctly_reads_the_table(self, table):
        frame = read_csv_file(table, 'csv', ['seq'], label_column='class')

        assert list(frame.columns) == ['seq', 'class']
        assert list(frame['seq']) == ['ACGT', 'TTTT']

    def test_the_streaming_reader_yields_the_sequences(self, table):
        assert list(stream_table_sequences(table, 'csv', ['seq'])) == ['ACGT', 'TTTT']
