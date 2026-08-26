"""Tests for evaluate_splits behaviour beyond the report layout.

The MMseqs2 binary is stubbed out, so these drive the surrounding pipeline:
the HTML report bundle, alignment lookup, input validation, and the failure and
cleanup paths.
"""

import io
import logging
import shutil

import pandas as pd
import pytest
from Bio import SeqIO
from helpers import mmseqs_hit, write_csv, write_mmseqs_output

from genomic_benchmarks_qc import evaluate_splits
from genomic_benchmarks_qc.utils.input_utils import append_fasta_record
from genomic_benchmarks_qc.utils.mmseqs_summary import MMSEQS_RESULT_COLUMNS, sequence_id
from genomic_benchmarks_qc.utils.naming import TMP_PREFIX


@pytest.fixture
def stub_mmseqs(monkeypatch):
    """Replace the MMseqs2 search with one that writes the given hits.

    Returns a factory so each test can choose what the search "found".
    """

    def _install(hits=(), error=None):
        def fake_run_search(query_fasta, target_fasta, output_path, tmp_dir, **kwargs):
            if error is not None:
                raise error
            return write_mmseqs_output(output_path, hits)

        monkeypatch.setattr(evaluate_splits.mmseqs_runtime, 'run_search', fake_run_search)

    return _install


@pytest.fixture
def split_inputs(tmp_path):
    return (
        write_csv(tmp_path / 'train.csv', ['0'], rows_per_label=5),
        write_csv(tmp_path / 'test.csv', ['0'], rows_per_label=5),
    )


class TestHtmlReportBundle:
    def test_bundle_contains_every_component(self, tmp_path, split_inputs, stub_mmseqs):
        train, test = split_inputs
        # Sequences are staged as seq_<i>_train / seq_<i>_test before searching.
        stub_mmseqs([mmseqs_hit('seq_0_test', 'seq_0_train'), mmseqs_hit('seq_1_test', 'seq_2_train')])

        evaluate_splits.run(
            train_files=[train], test_files=[test], format='csv',
            out_folder=str(tmp_path / 'out'), report_types=['html'],
        )

        comparison = tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'
        assert (comparison / 'gb-qc-report.html').is_file()
        assert (comparison / 'plots' / 'similarity_histograms.png').is_file()
        assert (comparison / 'mmseqs' / 'mmseqs2_search_result.tsv').is_file()
        assert (comparison / 'mmseqs' / 'seq_index_mapping' / 'test_sequences.fasta').is_file()
        assert (comparison / 'mmseqs' / 'seq_index_mapping' / 'train_sequences.fasta').is_file()

    def test_exported_hits_carry_the_expected_columns(self, tmp_path, split_inputs, stub_mmseqs):
        train, test = split_inputs
        stub_mmseqs([mmseqs_hit('seq_0_test', 'seq_0_train', pident=97.5)])

        evaluate_splits.run(
            train_files=[train], test_files=[test], format='csv',
            out_folder=str(tmp_path / 'out'), report_types=['html'],
        )

        exported = pd.read_csv(
            tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test' / 'mmseqs' / 'mmseqs2_search_result.tsv',
            sep='\t',
        )
        assert list(exported['query']) == ['seq_0_test']
        assert list(exported['target']) == ['seq_0_train']
        # min(qcov, tcov) * pident, the score leakage is judged on.
        assert exported['min_cov*pident'].iloc[0] == pytest.approx(0.99 * 97.5)

    def test_only_sequences_involved_in_hits_are_exported(self, tmp_path, split_inputs, stub_mmseqs):
        train, test = split_inputs
        stub_mmseqs([mmseqs_hit('seq_0_test', 'seq_3_train')])

        evaluate_splits.run(
            train_files=[train], test_files=[test], format='csv',
            out_folder=str(tmp_path / 'out'), report_types=['html'],
        )

        mapping = tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test' / 'mmseqs' / 'seq_index_mapping'
        assert mapping.joinpath('test_sequences.fasta').read_text().count('>') == 1
        assert '>seq_0_test' in mapping.joinpath('test_sequences.fasta').read_text()
        assert mapping.joinpath('train_sequences.fasta').read_text().count('>') == 1
        assert '>seq_3_train' in mapping.joinpath('train_sequences.fasta').read_text()

    def test_a_search_with_no_hits_still_produces_a_report(self, tmp_path, split_inputs, stub_mmseqs):
        train, test = split_inputs
        stub_mmseqs()

        evaluate_splits.run(
            train_files=[train], test_files=[test], format='csv',
            out_folder=str(tmp_path / 'out'), report_types=['html', 'simple'],
        )

        comparison = tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'
        assert (comparison / 'gb-qc-report.html').is_file()
        assert (comparison / 'gb-qc-report.csv').is_file()


class TestReportTypes:
    def test_a_type_this_command_cannot_write_is_refused(self, tmp_path, split_inputs):
        """The CLI checks this too, against the same constant. This is the
        library path, where the failure was a run that did all the work - a full
        MMseqs2 search - logged "successfully completed" and wrote nothing."""
        train, test = split_inputs

        with pytest.raises(ValueError, match="cannot produce report type"):
            evaluate_splits.run(
                train_files=[train], test_files=[test], format='csv',
                out_folder=str(tmp_path / 'out'), report_types=['json'],
            )

        assert not (tmp_path / 'out' / 'split').exists()


class TestEveryLeakedHitIsCountedAndExported:
    """The panel's count and the exported table both cover every leaked hit.

    The page lists at most ROW_CAP of them, and used to be handed that capped
    frame as if it were the whole finding: a split with 300 leaked pairs
    reported 100, exported 100, and never said it had dropped anything - a
    leakage report understating leakage.
    """

    def run_with_hits(self, tmp_path, stub_mmseqs, count, threshold=90.0):
        train = write_csv(tmp_path / 'train.csv', ['0'], rows_per_label=count)
        test = write_csv(tmp_path / 'test.csv', ['0'], rows_per_label=count)
        stub_mmseqs([mmseqs_hit(f'seq_{i}_test', f'seq_{i}_train', pident=100.0,
                                qcov=1.0, tcov=1.0)
                     for i in range(count)])

        evaluate_splits.run(
            train_files=[train], test_files=[test], format='csv',
            out_folder=str(tmp_path / 'out'), report_types=['html'],
            similarity_threshold=threshold,
        )
        return tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'

    def test_the_count_is_the_leaked_hits_not_the_listed_ones(
            self, tmp_path, stub_mmseqs):
        comparison = self.run_with_hits(tmp_path, stub_mmseqs, 250)

        page = (comparison / 'gb-qc-report.html').read_text()
        assert '250 high-similarity alignments (first 100 shown)' in page

    def test_the_export_holds_every_leaked_hit(self, tmp_path, stub_mmseqs):
        comparison = self.run_with_hits(tmp_path, stub_mmseqs, 250)

        exported = pd.read_csv(comparison / 'mmseqs' / 'mmseqs2_search_result.tsv', sep='\t')
        assert len(exported) == 250
        assert set(exported['query']) == {f'seq_{i}_test' for i in range(250)}

    def test_the_mapped_fastas_cover_the_export_not_the_listing(
            self, tmp_path, stub_mmseqs):
        """Both halves of mmseqs/ describe the same hits, or the ids in the TSV
        cannot be looked up in the FASTA beside it."""
        comparison = self.run_with_hits(tmp_path, stub_mmseqs, 250)

        mapping = comparison / 'mmseqs' / 'seq_index_mapping'
        assert mapping.joinpath('test_sequences.fasta').read_text().count('>') == 250
        assert mapping.joinpath('train_sequences.fasta').read_text().count('>') == 250

    def test_hits_below_the_threshold_are_left_out_of_both(
            self, tmp_path, stub_mmseqs):
        """The export is the leaked hits, not every alignment the search made -
        a search reports far more of the latter."""
        train, test = (write_csv(tmp_path / 'train.csv', ['0'], rows_per_label=5),
                       write_csv(tmp_path / 'test.csv', ['0'], rows_per_label=5))
        stub_mmseqs([
            mmseqs_hit('seq_0_test', 'seq_0_train', pident=100.0, qcov=1.0, tcov=1.0),
            mmseqs_hit('seq_1_test', 'seq_1_train', pident=50.0, qcov=1.0, tcov=1.0),
        ])

        evaluate_splits.run(
            train_files=[train], test_files=[test], format='csv',
            out_folder=str(tmp_path / 'out'), report_types=['html'],
        )

        comparison = tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'
        exported = pd.read_csv(comparison / 'mmseqs' / 'mmseqs2_search_result.tsv', sep='\t')
        assert list(exported['query']) == ['seq_0_test']
        assert '1 high-similarity alignment<' in (
            comparison / 'gb-qc-report.html').read_text()

    def test_a_clean_split_still_leaves_a_table(self, tmp_path, split_inputs, stub_mmseqs):
        """A header and no rows, so a reader who opens it sees an empty result
        rather than wondering whether the export ran."""
        train, test = split_inputs
        stub_mmseqs()

        evaluate_splits.run(
            train_files=[train], test_files=[test], format='csv',
            out_folder=str(tmp_path / 'out'), report_types=['html'],
        )

        exported = pd.read_csv(
            tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'
            / 'mmseqs' / 'mmseqs2_search_result.tsv', sep='\t')
        assert len(exported) == 0
        assert list(exported.columns) == MMSEQS_RESULT_COLUMNS


class TestStagedFasta:
    """MMseqs2 reads FASTA, so a CSV half is written out as one before the
    search starts. The records are numbered here and the sequences already
    uppercased, so they are written as their two lines rather than built into a
    Biopython record first - which was a record object per sequence, on inputs
    of several million, before any search had begun."""

    def test_a_record_is_two_lines_however_long_the_sequence(self):
        handle = io.StringIO()

        append_fasta_record(handle, 'ACGT' * 40, sequence_id(7, 'train'))

        assert handle.getvalue() == '>seq_7_train\n' + 'ACGT' * 40 + '\n'

    def test_biopython_reads_back_what_was_written(self):
        """The staged file is parsed again to attach sequences to the hits, so
        whatever is written here has to survive that round trip."""
        handle = io.StringIO()
        sequences = {sequence_id(i, 'test'): seq
                     for i, seq in enumerate(['ACGT' * 30, '', 'TT'])}
        for seq_id, sequence in sequences.items():
            append_fasta_record(handle, sequence, seq_id)

        parsed = {record.id: str(record.seq)
                  for record in SeqIO.parse(io.StringIO(handle.getvalue()), 'fasta')}

        assert parsed == sequences


class TestAddAlignmentSequences:
    def test_empty_hits_are_returned_unchanged(self, tmp_path):
        empty = pd.DataFrame(columns=['query', 'target'])

        result = evaluate_splits.add_alignment_sequences(
            empty, tmp_path / 'missing_test.fasta', tmp_path / 'missing_train.fasta'
        )

        # No file is read at all, so the paths above need not exist.
        assert result.empty

    def test_hit_sequences_are_attached_from_the_staged_fasta(self, tmp_path):
        (tmp_path / 'test.fasta').write_text('>seq_0_test\nACGT\n')
        (tmp_path / 'train.fasta').write_text('>seq_0_train\nTTTT\n')
        hits = pd.DataFrame([{'query': 'seq_0_test', 'target': 'seq_0_train'}])

        result = evaluate_splits.add_alignment_sequences(
            hits, tmp_path / 'test.fasta', tmp_path / 'train.fasta'
        )

        assert list(result['qseq']) == ['ACGT']
        assert list(result['tseq']) == ['TTTT']

    def test_unmappable_identifiers_raise_rather_than_render_blank_alignments(self, tmp_path):
        (tmp_path / 'test.fasta').write_text('>seq_0_test\nACGT\n')
        (tmp_path / 'train.fasta').write_text('>seq_0_train\nTTTT\n')
        hits = pd.DataFrame([{'query': 'seq_0_test', 'target': 'seq_99_train'}])

        with pytest.raises(RuntimeError, match='Failed to map MMSeqs2 hit identifiers'):
            evaluate_splits.add_alignment_sequences(
                hits, tmp_path / 'test.fasta', tmp_path / 'train.fasta'
            )


class TestInputValidation:
    def test_an_empty_input_side_is_rejected(self, tmp_path, stub_mmseqs):
        stub_mmseqs()
        train = write_csv(tmp_path / 'train.csv', ['0'], rows_per_label=5)
        empty = tmp_path / 'empty.csv'
        empty.write_text('sequence,label\n')

        with pytest.raises(ValueError, match='at least one sequence'):
            evaluate_splits.run(
                train_files=[train], test_files=[str(empty)], format='csv',
                out_folder=str(tmp_path / 'out'), report_types=['simple'],
            )


class TestFailureHandling:
    def test_a_failing_search_is_logged_and_re_raised(self, tmp_path, split_inputs, stub_mmseqs, caplog):
        train, test = split_inputs
        stub_mmseqs(error=RuntimeError('mmseqs exploded'))

        with (
            caplog.at_level(logging.ERROR),
            pytest.raises(RuntimeError, match='mmseqs exploded'),
        ):
            evaluate_splits.run(
                train_files=[train], test_files=[test], format='csv',
                out_folder=str(tmp_path / 'out'), report_types=['simple'],
            )

        assert 'Train-test split evaluation failed' in caplog.text

    def test_debug_logging_records_the_traceback(self, tmp_path, split_inputs, stub_mmseqs, caplog):
        """At DEBUG level the failure is logged with its traceback rather than one line."""
        train, test = split_inputs
        stub_mmseqs(error=RuntimeError('mmseqs exploded'))

        with caplog.at_level(logging.DEBUG), pytest.raises(RuntimeError):
            evaluate_splits.run(
                train_files=[train], test_files=[test], format='csv',
                out_folder=str(tmp_path / 'out'), report_types=['simple'],
                log_level='DEBUG',
            )

        assert 'Traceback' in caplog.text

    def test_temporary_files_are_removed_even_when_the_run_fails(self, tmp_path, split_inputs, stub_mmseqs):
        train, test = split_inputs
        stub_mmseqs(error=RuntimeError('mmseqs exploded'))

        with pytest.raises(RuntimeError):
            evaluate_splits.run(
                train_files=[train], test_files=[test], format='csv',
                out_folder=str(tmp_path / 'out'), report_types=['simple'],
            )

        comparison = tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'
        assert list(comparison.glob(f'{TMP_PREFIX}*')) == []

    def test_a_failed_cleanup_warns_instead_of_failing_the_run(self, tmp_path, split_inputs, stub_mmseqs, monkeypatch, caplog):
        """A successful analysis must not be lost to an unremovable temp directory."""
        train, test = split_inputs
        stub_mmseqs()

        def refuse_to_remove(path):
            raise OSError('device busy')

        monkeypatch.setattr(shutil, 'rmtree', refuse_to_remove)

        with caplog.at_level(logging.WARNING):
            evaluate_splits.run(
                train_files=[train], test_files=[test], format='csv',
                out_folder=str(tmp_path / 'out'), report_types=['simple'],
            )

        assert 'Failed to remove temporary directory' in caplog.text
        assert (tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test' / 'gb-qc-report.csv').is_file()


class TestDefaults:
    def test_optional_arguments_default_when_omitted(self, tmp_path, stub_mmseqs):
        stub_mmseqs()
        train = write_csv(tmp_path / 'train.csv', ['0'], rows_per_label=5)
        test = write_csv(tmp_path / 'test.csv', ['0'], rows_per_label=5)

        evaluate_splits.run(
            train_files=[train], test_files=[test], format='csv',
            out_folder=str(tmp_path / 'out'),
            sequence_column=None,
            report_types=None,
        )

        comparison = tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'
        assert (comparison / 'gb-qc-report.csv').is_file()
        assert (comparison / 'gb-qc-report.html').is_file()

    def test_an_explicit_sequence_column_is_used(self, tmp_path, stub_mmseqs):
        stub_mmseqs()
        train = write_csv(tmp_path / 'train.csv', ['0'], columns=('seq',), rows_per_label=5)
        test = write_csv(tmp_path / 'test.csv', ['0'], columns=('seq',), rows_per_label=5)

        evaluate_splits.run(
            train_files=[train], test_files=[test], format='csv',
            out_folder=str(tmp_path / 'out'),
            sequence_column=['seq'],
            report_types=['simple'],
        )

        # The column directory is named after the column that was searched.
        assert (tmp_path / 'out' / 'split' / 'seq' / 'train_vs_test' / 'gb-qc-report.csv').is_file()

    def test_an_existing_comparison_folder_is_reused(self, tmp_path, split_inputs, stub_mmseqs):
        train, test = split_inputs
        stub_mmseqs()
        comparison = tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'
        comparison.mkdir(parents=True)

        evaluate_splits.run(
            train_files=[train], test_files=[test], format='csv',
            out_folder=str(tmp_path / 'out'), report_types=['simple'],
        )

        assert (comparison / 'gb-qc-report.csv').is_file()
