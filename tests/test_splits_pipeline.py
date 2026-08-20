"""Tests for evaluate_splits behaviour beyond the report layout.

The MMseqs2 binary is stubbed out, so these drive the surrounding pipeline:
the HTML report bundle, alignment lookup, input validation, and the failure and
cleanup paths.
"""

import logging
import shutil

import pandas as pd
import pytest
from helpers import mmseqs_hit, write_csv, write_mmseqs_output

from genomic_benchmarks_qc import evaluate_splits
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
