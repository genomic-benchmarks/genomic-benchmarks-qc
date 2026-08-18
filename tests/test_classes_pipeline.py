"""Tests for evaluate_classes behaviour that the layout tests do not reach.

Covers the regression target handling, the single-class short circuit, input
merging, and the argument defaults.
"""

import logging

import pytest

from helpers import write_csv, write_fasta, write_regression_csv

from genomic_benchmarks_qc import evaluate_classes


def _reports(root):
    return sorted(str(p.relative_to(root)) for p in root.rglob('report.csv'))


class TestRegressionTargets:
    def test_numeric_target_is_split_into_high_and_low(self, tmp_path):
        csv_path = write_regression_csv(tmp_path / 'data.csv', [str(i) for i in range(20)])

        evaluate_classes.run(
            input=[csv_path],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            regression=True,
            report_types=['simple'],
        )

        assert _reports(tmp_path / 'out') == ['class/sequence/high_vs_low/report.csv']

    def test_non_numeric_values_are_dropped_with_a_warning(self, tmp_path, caplog):
        values = [str(i) for i in range(18)] + ['not-a-number', 'also-not']
        csv_path = write_regression_csv(tmp_path / 'data.csv', values)

        with caplog.at_level(logging.WARNING):
            evaluate_classes.run(
                input=[csv_path],
                format='csv',
                out_folder=str(tmp_path / 'out'),
                regression=True,
                report_types=['simple'],
            )

        assert 'Dropped 2 rows with non-numeric values' in caplog.text
        assert _reports(tmp_path / 'out') == ['class/sequence/high_vs_low/report.csv']

    def test_a_fully_non_numeric_target_is_an_error(self, tmp_path, caplog):
        csv_path = write_regression_csv(tmp_path / 'data.csv', ['abc'] * 10)

        with caplog.at_level(logging.ERROR):
            with pytest.raises(ValueError, match='contains no numeric values'):
                evaluate_classes.run(
                    input=[csv_path],
                    format='csv',
                    out_folder=str(tmp_path / 'out'),
                    regression=True,
                    report_types=['simple'],
                )

        # Reported to the log as well as raised, so it reaches --log-file.
        assert 'contains no numeric values' in caplog.text
        assert _reports(tmp_path / 'out') == []

    @pytest.mark.parametrize(
        "values",
        [
            pytest.param(['5.0'] * 10, id='every-value-identical'),
            pytest.param(['0'] * 8 + ['1', '2'], id='zero-inflated-median-equals-minimum'),
            pytest.param(['abc'] * 9 + ['7.0'], id='single-row-survives-numeric-filter'),
        ],
    )
    def test_a_target_that_cannot_be_split_into_two_classes_is_an_error(self, tmp_path, values, caplog):
        """The split is `value >= median -> high`, so a target whose median equals
        its minimum puts everything in 'high' and leaves 'low' empty. There is no
        meaningful comparison to run, so the tool reports it and stops.
        """
        csv_path = write_regression_csv(tmp_path / 'data.csv', values)

        with caplog.at_level(logging.ERROR):
            with pytest.raises(ValueError, match='cannot be split into two classes'):
                evaluate_classes.run(
                    input=[csv_path],
                    format='csv',
                    out_folder=str(tmp_path / 'out'),
                    regression=True,
                    report_types=['simple'],
                )

        # The message must name the cause, not surface as a missing-label error.
        assert "zero-inflated" in caplog.text
        assert _reports(tmp_path / 'out') == []


class TestSingleClass:
    def test_one_class_produces_per_class_output_but_no_comparison(self, tmp_path):
        csv_path = write_csv(tmp_path / 'data.csv', ['only'])

        evaluate_classes.run(
            input=[csv_path],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            report_types=['simple', 'json'],
        )

        column_dir = tmp_path / 'out' / 'class' / 'sequence'
        assert (column_dir / 'per-class' / 'only.json').is_file()
        assert _reports(tmp_path / 'out') == []


class TestInputMerging:
    def test_multiple_input_files_are_pooled_into_one_dataset(self, tmp_path, caplog):
        first = write_csv(tmp_path / 'first.csv', ['a', 'b'], rows_per_label=15)
        second = write_csv(tmp_path / 'second.csv', ['a', 'b'], rows_per_label=15)

        with caplog.at_level(logging.INFO):
            evaluate_classes.run(
                input=[first, second],
                format='csv',
                out_folder=str(tmp_path / 'out'),
                report_types=['simple', 'json'],
            )

        assert 'Merging 2 input files' in caplog.text
        assert _reports(tmp_path / 'out') == ['class/sequence/a_vs_b/report.csv']

        # Reports name the merged source rather than either individual file.
        import json
        stats = json.loads((tmp_path / 'out' / 'class' / 'sequence' / 'per-class' / 'a.json').read_text())
        assert stats['Filename'] == 'merged'
        assert stats['Number of sequences'] == 30


class TestDefaults:
    def test_run_defaults_to_html_and_simple_reports(self, tmp_path):
        csv_path = write_csv(tmp_path / 'data.csv', ['a', 'b'])

        evaluate_classes.run(
            input=[csv_path],
            format='csv',
            out_folder=str(tmp_path / 'out'),
            report_types=None,
        )

        comparison = tmp_path / 'out' / 'class' / 'sequence' / 'a_vs_b'
        assert (comparison / 'report.csv').is_file()
        assert (comparison / 'report.html').is_file()

    def test_run_analysis_defaults_to_html_and_simple_reports(self, tmp_path):
        """run_analysis is part of the public module surface and defaults on its own."""
        csv_path = write_csv(tmp_path / 'data.csv', ['a', 'b'])
        evaluate_classes.run(
            input=[csv_path], format='csv', out_folder=str(tmp_path / 'seed'), report_types=['json'],
        )

        from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics
        from helpers import sequences as make_sequences

        stats = [
            SequenceStatistics(make_sequences(20, seed=1), filename='a.csv', filepath='a.csv', label='a'),
            SequenceStatistics(make_sequences(20, seed=2), filename='a.csv', filepath='a.csv', label='b'),
        ]
        evaluate_classes.run_analysis(
            input_statistics=stats,
            report_dir=tmp_path / 'out',
            report_types=None,
            plot_type='boxen',
        )

        assert (tmp_path / 'out' / 'a_vs_b' / 'report.csv').is_file()
        assert (tmp_path / 'out' / 'a_vs_b' / 'report.html').is_file()

    def test_an_existing_output_folder_is_reused(self, tmp_path):
        out_folder = tmp_path / 'out'
        out_folder.mkdir()
        (out_folder / 'pre-existing.txt').write_text('kept')

        evaluate_classes.run(
            input=[
                write_fasta(tmp_path / 'pos.fa', 30, seed=1),
                write_fasta(tmp_path / 'neg.fa', 30, seed=2),
            ],
            format='fasta',
            out_folder=str(out_folder),
            report_types=['simple'],
        )

        assert (out_folder / 'pre-existing.txt').read_text() == 'kept'
        assert _reports(out_folder) == ['class/sequence/neg_vs_pos/report.csv']
