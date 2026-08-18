"""Tests for the on-disk layout of generated reports.

The tool owns everything from the 'class/' and 'split/' directories downwards;
any grouping above that (collection, dataset, split) is expressed by the caller
through the output folder. These tests pin the part the tool owns, and cover the
inputs whose names used to collide into a shared report path.
"""

import pathlib

import pandas as pd
import pytest

from helpers import sequences as _sequences, write_csv as _write_csv, write_fasta as _write_fasta, write_mmseqs_output

from genomic_benchmarks_qc import evaluate_classes, evaluate_splits


class TestClassesLayout:
    def test_reports_land_in_a_directory_per_comparison(self, tmp_path):
        inputs = [
            _write_fasta(tmp_path / 'pos.fa', 40, seed=1),
            _write_fasta(tmp_path / 'neg.fa', 40, seed=2),
        ]

        evaluate_classes.run(
            input=inputs,
            format='fasta',
            out_folder=str(tmp_path / 'out'),
            report_types=['simple', 'json'],
        )

        column_dir = tmp_path / 'out' / 'class' / 'sequence'
        assert (column_dir / 'neg_vs_pos' / 'report.csv').is_file()
        assert (column_dir / 'per-class' / 'pos.json').is_file()
        assert (column_dir / 'per-class' / 'neg.json').is_file()

    def test_fasta_files_sharing_a_stem_get_separate_directories(self, tmp_path):
        """Same file name in different directories must not share a report path."""
        inputs = [
            _write_fasta(tmp_path / 'train' / 'pos.fa', 40, seed=1),
            _write_fasta(tmp_path / 'test' / 'pos.fa', 40, seed=2),
            _write_fasta(tmp_path / 'neg.fa', 40, seed=3),
        ]

        evaluate_classes.run(
            input=inputs,
            format='fasta',
            out_folder=str(tmp_path / 'out'),
            report_types=['simple', 'json'],
        )

        class_dir = tmp_path / 'out' / 'class' / 'sequence'
        # Three inputs give three pairwise comparisons; all must be distinct.
        comparisons = sorted(d.name for d in class_dir.iterdir() if d.name != 'per-class')
        assert comparisons == [
            'neg_vs_test-pos',
            'neg_vs_train-pos',
            'test-pos_vs_train-pos',
        ]
        for comparison in comparisons:
            assert (class_dir / comparison / 'report.csv').is_file()

        assert sorted(p.name for p in (class_dir / 'per-class').iterdir()) == [
            'neg.json',
            'test-pos.json',
            'train-pos.json',
        ]

    def test_labels_unusable_in_a_path_are_slugged(self, tmp_path):
        """A label containing '/' used to abort the run with FileNotFoundError."""
        csv_path = _write_csv(tmp_path / 'data.csv', ["5' UTR", 'A/B'])

        evaluate_classes.run(
            input=[csv_path],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            report_types=['simple'],
        )

        assert (tmp_path / 'out' / 'class' / 'sequence' / '5-utr_vs_a-b' / 'report.csv').is_file()

    def test_labels_keep_their_original_form_in_the_report(self, tmp_path):
        """Slugging is for paths only; the report still shows the real label."""
        csv_path = _write_csv(tmp_path / 'data.csv', ["5' UTR", 'A/B'])

        evaluate_classes.run(
            input=[csv_path],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            report_types=['json'],
        )

        per_class = tmp_path / 'out' / 'class' / 'sequence' / 'per-class'
        assert pd.read_json(per_class / '5-utr.json', typ='series')['Label'] == "5' UTR"
        assert pd.read_json(per_class / 'a-b.json', typ='series')['Label'] == 'A/B'

    def test_each_sequence_column_gets_its_own_directory(self, tmp_path):
        csv_path = _write_csv(tmp_path / 'data.csv', ['0', '1'], columns=('seq_a', 'seq_b'))

        evaluate_classes.run(
            input=[csv_path],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            sequence_column=['seq_a', 'seq_b'],
            report_types=['simple'],
        )

        class_dir = tmp_path / 'out' / 'class'
        assert (class_dir / 'seq_a' / '0_vs_1' / 'report.csv').is_file()
        assert (class_dir / 'seq_b' / '0_vs_1' / 'report.csv').is_file()
        # Multiple columns are additionally analysed concatenated together.
        assert (class_dir / 'merged' / '0_vs_1' / 'report.csv').is_file()

    def test_a_column_named_merged_does_not_take_the_merged_directory(self, tmp_path):
        csv_path = _write_csv(tmp_path / 'data.csv', ['0', '1'], columns=('merged', 'seq_b'))

        evaluate_classes.run(
            input=[csv_path],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            sequence_column=['merged', 'seq_b'],
            report_types=['simple'],
        )

        class_dir = tmp_path / 'out' / 'class'
        column_dirs = sorted(d.name for d in class_dir.iterdir())
        assert len(column_dirs) == 3, column_dirs
        assert 'merged' in column_dirs

    def test_duplicate_sequences_are_written_beside_the_html_report(self, tmp_path):
        shared = _sequences(10, seed=7)
        csv_path = tmp_path / 'data.csv'
        pd.DataFrame(
            {
                'sequence': shared + shared + _sequences(20, seed=8),
                'label': ['a'] * 10 + ['b'] * 10 + ['b'] * 20,
            }
        ).to_csv(csv_path, index=False)

        evaluate_classes.run(
            input=[str(csv_path)],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            report_types=['html'],
        )

        comparison = tmp_path / 'out' / 'class' / 'sequence' / 'a_vs_b'
        assert (comparison / 'report.html').is_file()
        assert (comparison / 'plots').is_dir()
        assert (comparison / 'duplicates.txt').is_file()
        assert set((comparison / 'duplicates.txt').read_text().split()) == set(shared)


class TestSplitsLayout:
    @pytest.fixture
    def stub_mmseqs(self, monkeypatch):
        """Stand in for MMseqs2, which is an external binary the tests cannot run.

        Returns an empty hit table, which is enough to exercise every output path
        the command builds.
        """

        def fake_run_search(query_fasta, target_fasta, output_path, tmp_dir, **kwargs):
            # The real search writes its raw TSV here; downstream code only reads it.
            return write_mmseqs_output(output_path)

        monkeypatch.setattr(evaluate_splits.mmseqs_runtime, 'run_search', fake_run_search)

    def test_reports_land_in_a_directory_per_comparison(self, tmp_path, stub_mmseqs):
        train = _write_csv(tmp_path / 'enhancers_train.csv', ['0'])
        test = _write_csv(tmp_path / 'enhancers_test.csv', ['0'])

        evaluate_splits.run(
            train_files=[train],
            test_files=[test],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            report_types=['simple'],
        )

        comparison = tmp_path / 'out' / 'split' / 'sequence' / 'enhancers_train_vs_enhancers_test'
        assert (comparison / 'report.csv').is_file()

    def test_train_and_test_sharing_a_stem_get_separate_directories(self, tmp_path, stub_mmseqs):
        """'train/data.csv' vs 'val/data.csv' and vs 'test/data.csv' must differ.

        Both used to reduce to a single 'data_vs_data' stem, so the second run
        overwrote the first one's reports.
        """
        train = _write_csv(tmp_path / 'train' / 'data.csv', ['0'])
        out_folder = str(tmp_path / 'out')

        for split in ('val', 'test'):
            evaluate_splits.run(
                train_files=[train],
                test_files=[_write_csv(tmp_path / split / 'data.csv', ['0'])],
                format='csv',
                out_folder=out_folder,
                report_types=['simple'],
            )

        comparisons = sorted(d.name for d in (tmp_path / 'out' / 'split' / 'sequence').iterdir())
        assert comparisons == ['train-data_vs_test-data', 'train-data_vs_val-data']

    def test_temporary_files_are_removed_from_the_comparison_directory(self, tmp_path, stub_mmseqs):
        train = _write_csv(tmp_path / 'train.csv', ['0'])
        test = _write_csv(tmp_path / 'test.csv', ['0'])

        evaluate_splits.run(
            train_files=[train],
            test_files=[test],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            report_types=['simple'],
        )

        comparison = tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'
        assert not (comparison / 'tmp').exists()

    def test_temporary_files_are_kept_inside_the_comparison_directory(self, tmp_path, stub_mmseqs):
        """Per-comparison tmp keeps concurrent runs sharing an output folder apart."""
        train = _write_csv(tmp_path / 'train.csv', ['0'])
        test = _write_csv(tmp_path / 'test.csv', ['0'])

        evaluate_splits.run(
            train_files=[train],
            test_files=[test],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            report_types=['simple'],
            keep_tmp_files=True,
        )

        comparison = tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'
        assert (comparison / 'tmp').is_dir()
        assert not (tmp_path / 'out' / 'tmp').exists()

    def test_the_searched_column_names_the_directory(self, tmp_path, stub_mmseqs):
        train = _write_csv(tmp_path / 'train.csv', ['0'], columns=('gene',))
        test = _write_csv(tmp_path / 'test.csv', ['0'], columns=('gene',))

        evaluate_splits.run(
            train_files=[train],
            test_files=[test],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            sequence_column=['gene'],
            report_types=['simple'],
        )

        assert (tmp_path / 'out' / 'split' / 'gene' / 'train_vs_test' / 'report.csv').is_file()

    def test_several_columns_are_reported_as_merged(self, tmp_path, stub_mmseqs):
        """The columns are concatenated into one search, as 'class/merged/' is."""
        train = _write_csv(tmp_path / 'train.csv', ['0'], columns=('gene', 'rna'))
        test = _write_csv(tmp_path / 'test.csv', ['0'], columns=('gene', 'rna'))

        evaluate_splits.run(
            train_files=[train],
            test_files=[test],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            sequence_column=['gene', 'rna'],
            report_types=['simple'],
        )

        assert (tmp_path / 'out' / 'split' / 'merged' / 'train_vs_test' / 'report.csv').is_file()

    def test_fasta_inputs_use_the_default_column_directory(self, tmp_path, stub_mmseqs):
        evaluate_splits.run(
            train_files=[_write_fasta(tmp_path / 'train.fa', 40, seed=1)],
            test_files=[_write_fasta(tmp_path / 'test.fa', 40, seed=2)],
            format='fasta',
            out_folder=str(tmp_path / 'out'),
            report_types=['simple'],
        )

        assert (tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test' / 'report.csv').is_file()


class TestLayoutIsFormatIndependent:
    def test_fasta_and_csv_reports_sit_at_the_same_depth(self, tmp_path):
        """FASTA has no sequence column, but its reports still go one level deep.

        A caller can then resolve any report as class/<column>/<pair>/report.html
        without branching on the input format.
        """
        evaluate_classes.run(
            input=[
                _write_fasta(tmp_path / 'pos.fa', 40, seed=1),
                _write_fasta(tmp_path / 'neg.fa', 40, seed=2),
            ],
            format='fasta',
            out_folder=str(tmp_path / 'from_fasta'),
            report_types=['simple'],
        )
        evaluate_classes.run(
            input=[_write_csv(tmp_path / 'data.csv', ['pos', 'neg'])],
            format='csv',
            out_folder=str(tmp_path / 'from_csv'),
            report_types=['simple'],
        )

        def layout(root):
            return sorted(str(p.relative_to(root)) for p in root.rglob('report.csv'))

        # Identical paths, not merely identical depth: classes are sorted by path
        # name for both formats, so 'pos' and 'neg' give 'neg_vs_pos' either way.
        assert layout(tmp_path / 'from_fasta') == ['class/sequence/neg_vs_pos/report.csv']
        assert layout(tmp_path / 'from_csv') == ['class/sequence/neg_vs_pos/report.csv']

    def test_both_commands_report_at_the_same_depth(self, tmp_path, monkeypatch):
        """<command>/<column>/<comparison>/ holds every report either command writes."""

        def fake_run_search(query_fasta, target_fasta, output_path, tmp_dir, **kwargs):
            return write_mmseqs_output(output_path)

        monkeypatch.setattr(evaluate_splits.mmseqs_runtime, 'run_search', fake_run_search)

        out_folder = tmp_path / 'out'
        train = _write_csv(tmp_path / 'train.csv', ['pos', 'neg'])
        test = _write_csv(tmp_path / 'test.csv', ['pos', 'neg'])

        evaluate_classes.run(
            input=[train, test],
            format='csv',
            out_folder=str(out_folder),
            report_types=['simple'],
        )
        evaluate_splits.run(
            train_files=[train],
            test_files=[test],
            format='csv',
            out_folder=str(out_folder),
            report_types=['simple'],
        )

        reports = sorted(str(p.relative_to(out_folder)) for p in out_folder.rglob('report.csv'))
        assert reports == [
            'class/sequence/neg_vs_pos/report.csv',
            'split/sequence/train_vs_test/report.csv',
        ]
        assert all(len(pathlib.Path(report).parts) == 4 for report in reports)


class TestClassOrdering:
    """Classes are sorted by path name, so report paths are input-order independent."""

    def test_input_file_order_does_not_change_report_paths(self, tmp_path):
        pos = _write_fasta(tmp_path / 'pos.fa', 40, seed=1)
        neg = _write_fasta(tmp_path / 'neg.fa', 40, seed=2)

        for name, inputs in (('forward', [pos, neg]), ('reversed', [neg, pos])):
            evaluate_classes.run(
                input=inputs,
                format='fasta',
                out_folder=str(tmp_path / name),
                report_types=['simple'],
            )

        def layout(name):
            root = tmp_path / name
            return sorted(str(p.relative_to(root)) for p in root.rglob('report.csv'))

        assert layout('forward') == layout('reversed') == ['class/sequence/neg_vs_pos/report.csv']

    def test_explicit_label_list_order_does_not_change_report_paths(self, tmp_path):
        csv_path = _write_csv(tmp_path / 'data.csv', ['pos', 'neg'])

        for name, label_list in (('forward', ['pos', 'neg']), ('reversed', ['neg', 'pos'])):
            evaluate_classes.run(
                input=[csv_path],
                format='csv',
                out_folder=str(tmp_path / name),
                label_list=label_list,
                report_types=['simple'],
            )

        for name in ('forward', 'reversed'):
            assert (tmp_path / name / 'class' / 'sequence' / 'neg_vs_pos' / 'report.csv').is_file()
