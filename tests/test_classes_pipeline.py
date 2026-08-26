"""Tests for evaluate_classes behaviour that the layout tests do not reach.

Covers the regression target handling, the single-class short circuit, input
merging, and the argument defaults.
"""

import logging

import pytest
from helpers import write_csv, write_fasta, write_regression_csv

from genomic_benchmarks_qc import evaluate_classes


def _reports(root):
    return sorted(str(p.relative_to(root)) for p in root.rglob('gb-qc-report.csv'))


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

        assert _reports(tmp_path / 'out') == ['class/sequence/high_vs_low/gb-qc-report.csv']

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
        assert _reports(tmp_path / 'out') == ['class/sequence/high_vs_low/gb-qc-report.csv']

    def test_a_fully_non_numeric_target_is_an_error(self, tmp_path, caplog):
        csv_path = write_regression_csv(tmp_path / 'data.csv', ['abc'] * 10)

        with (
            caplog.at_level(logging.ERROR),
            pytest.raises(ValueError, match='contains no numeric values'),
        ):
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

        with (
            caplog.at_level(logging.ERROR),
            pytest.raises(ValueError, match='cannot be split into two classes'),
        ):
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


class TestInputValidation:
    def test_an_empty_fasta_class_is_rejected(self, tmp_path, caplog):
        empty = tmp_path / 'empty.fa'
        empty.write_text('')

        with (
            caplog.at_level(logging.ERROR),
            pytest.raises(ValueError, match='No sequences found in FASTA file'),
        ):
            evaluate_classes.run(
                input=[write_fasta(tmp_path / 'pos.fa', 30, seed=1), str(empty)],
                format='fasta',
                out_folder=str(tmp_path / 'out'),
                report_types=['simple'],
            )

        assert 'empty.fa' in caplog.text
        # Rejected while reading, so no partial reports are left behind.
        assert _reports(tmp_path / 'out') == []

    def test_an_unknown_report_type_is_refused(self, tmp_path):
        """Silently skipping one would mean a run that does all the work and
        writes nothing - which is what the splits command used to do with json."""
        with pytest.raises(ValueError, match="cannot produce report type"):
            evaluate_classes.run(
                input=[write_csv(tmp_path / 'data.csv', ['a', 'b'], rows_per_label=15)],
                format='csv',
                out_folder=str(tmp_path / 'out'),
                report_types=['pdf'],
            )

        assert _reports(tmp_path / 'out') == []


class TestInferredLabels:
    """`--label-list infer` trusts the column it is pointed at, and every pair of
    classes gets its own comparison - so the cost of pointing it at the wrong
    column is quadratic and unbounded. A continuous target without --regression
    inferred one class per row and wrote tens of thousands of directories before
    it was stopped.
    """

    def test_a_column_with_too_many_values_is_refused(self, tmp_path):
        labels = [str(value) for value in range(evaluate_classes.MAX_INFERRED_LABELS + 1)]
        data = write_csv(tmp_path / 'data.csv', labels, rows_per_label=2)

        with pytest.raises(ValueError, match='distinct values'):
            evaluate_classes.run(input=[data], format='csv',
                                 out_folder=str(tmp_path / 'out'),
                                 report_types=['simple'])

        # Refused before any comparison ran, so nothing is left to clean up.
        assert _reports(tmp_path / 'out') == []

    def test_the_message_says_what_to_do_instead(self, tmp_path):
        """The three ways out, because the column that triggers this is usually a
        regression target or the sequence column itself."""
        labels = [f'value_{value}' for value in range(evaluate_classes.MAX_INFERRED_LABELS + 1)]
        data = write_csv(tmp_path / 'data.csv', labels, rows_per_label=2)

        with pytest.raises(ValueError) as failure:
            evaluate_classes.run(input=[data], format='csv',
                                 out_folder=str(tmp_path / 'out'),
                                 report_types=['simple'])

        message = str(failure.value)
        assert '--regression' in message
        assert '--label-list' in message
        assert '--label-column' in message
        assert 'value_0' in message                     # what it actually found

    def test_a_long_label_is_shortened_in_the_message(self, tmp_path):
        """The column pointed at is often the sequence column, whose values would
        otherwise make the message unreadable."""
        labels = [f'{"ACGT" * 20}{value}'
                  for value in range(evaluate_classes.MAX_INFERRED_LABELS + 1)]
        data = write_csv(tmp_path / 'data.csv', labels, rows_per_label=2)

        with pytest.raises(ValueError) as failure:
            evaluate_classes.run(input=[data], format='csv',
                                 out_folder=str(tmp_path / 'out'),
                                 report_types=['simple'])

        assert '\u2026' in str(failure.value)
        assert 'ACGT' * 20 not in str(failure.value)

    def test_the_cap_itself_is_allowed(self, tmp_path):
        labels = [str(value) for value in range(evaluate_classes.MAX_INFERRED_LABELS)]
        data = write_csv(tmp_path / 'data.csv', labels, rows_per_label=2)

        evaluate_classes.run(input=[data], format='csv',
                             out_folder=str(tmp_path / 'out'),
                             report_types=['simple'])

        pairs = evaluate_classes.MAX_INFERRED_LABELS * (evaluate_classes.MAX_INFERRED_LABELS - 1)
        assert len(_reports(tmp_path / 'out')) == pairs // 2

    def test_an_explicit_list_is_not_capped(self, tmp_path):
        """A list the caller typed out is a decision, not a mistake to catch."""
        labels = [str(value) for value in range(evaluate_classes.MAX_INFERRED_LABELS + 1)]
        data = write_csv(tmp_path / 'data.csv', labels, rows_per_label=2)

        evaluate_classes.run(input=[data], format='csv',
                             out_folder=str(tmp_path / 'out'),
                             label_list=labels[:3], report_types=['simple'])

        assert len(_reports(tmp_path / 'out')) == 3


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
        assert _reports(tmp_path / 'out') == ['class/sequence/a_vs_b/gb-qc-report.csv']

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
        assert (comparison / 'gb-qc-report.csv').is_file()
        assert (comparison / 'gb-qc-report.html').is_file()

    def test_run_analysis_defaults_to_html_and_simple_reports(self, tmp_path):
        """run_analysis is part of the public module surface and defaults on its own."""
        csv_path = write_csv(tmp_path / 'data.csv', ['a', 'b'])
        evaluate_classes.run(
            input=[csv_path], format='csv', out_folder=str(tmp_path / 'seed'), report_types=['json'],
        )

        from helpers import sequences as make_sequences

        from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics

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

        assert (tmp_path / 'out' / 'a_vs_b' / 'gb-qc-report.csv').is_file()
        assert (tmp_path / 'out' / 'a_vs_b' / 'gb-qc-report.html').is_file()

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
        assert _reports(out_folder) == ['class/sequence/neg_vs_pos/gb-qc-report.csv']
