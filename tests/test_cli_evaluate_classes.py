"""Tests for the `gb-qc evaluate-classes` command.

The pipeline itself is stubbed out (see the ``classes_run`` fixture), so these
tests cover what the CLI layer is responsible for: argument parsing, validation
and the arguments handed to the pipeline.
"""

import pytest
import typer

from genomic_benchmarks_qc.cli import app


def invoke(runner, *args):
    return runner.invoke(app, ["evaluate-classes", *args])


class TestArgumentPassing:
    def test_defaults_are_forwarded(self, runner, make_file, classes_run):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv)

        assert result.exit_code == 0
        assert classes_run.kwargs == {
            'input': [csv],
            'format': 'csv',
            'out_folder': '.',
            'sequence_column': ['sequence'],
            'label_column': 'label',
            'label_list': ['infer'],
            'regression': False,
            'report_types': ['html', 'simple'],
            'end_position': None,
            'plot_type': 'boxen',
            'log_level': 'INFO',
            'log_file': None,
        }

    def test_explicit_options_are_forwarded(self, runner, make_file, classes_run, tmp_path):
        csv = make_file("train.csv")
        out = str(tmp_path / "out")
        log = str(tmp_path / "run.log")

        result = invoke(
            runner,
            "--input", csv,
            "--sequence-column", "seq_a",
            "--sequence-column", "seq_b",
            "--label-column", "target",
            "--label-list", "0",
            "--label-list", "1",
            "--regression",
            "--out-folder", out,
            "--report-types", "json",
            "--end-position", "200",
            "--plot-type", "violin",
            "--log-level", "DEBUG",
            "--log-file", log,
        )

        assert result.exit_code == 0
        # Repeatable options replace their defaults rather than extending them.
        assert classes_run.kwargs == {
            'input': [csv],
            'format': 'csv',
            'out_folder': out,
            'sequence_column': ['seq_a', 'seq_b'],
            'label_column': 'target',
            'label_list': ['0', '1'],
            'regression': True,
            'report_types': ['json'],
            'end_position': 200,
            'plot_type': 'violin',
            'log_level': 'DEBUG',
            'log_file': log,
        }

    def test_multiple_inputs_are_forwarded_in_order(self, runner, make_file, classes_run):
        first = make_file("positives.fasta")
        second = make_file("negatives.fasta")

        result = invoke(runner, "--input", first, "--input", second)

        assert result.exit_code == 0
        assert classes_run.kwargs['input'] == [first, second]

    def test_input_is_required(self, runner, classes_run):
        result = invoke(runner)

        assert result.exit_code == 2
        assert not classes_run.called

    @pytest.mark.parametrize("flag, expected", [("--regression", True), ("--no-regression", False)])
    def test_regression_flag(self, runner, make_file, classes_run, flag, expected):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, flag)

        assert result.exit_code == 0
        assert classes_run.kwargs['regression'] is expected


class TestFormatHandling:
    @pytest.mark.parametrize(
        "names, expected_format",
        [
            (["train.csv"], 'csv'),
            (["train.csv.gz"], 'csv'),
            (["train.tsv"], 'tsv'),
            (["train.tsv.gz"], 'tsv'),
            (["pos.fasta", "neg.fasta"], 'fasta'),
            (["pos.fa", "neg.fa"], 'fasta'),
            (["pos.fasta.gz", "neg.fa.gz"], 'fasta'),
        ],
    )
    def test_format_is_inferred_and_normalized(self, runner, make_file, classes_run, names, expected_format):
        files = [make_file(name) for name in names]

        result = invoke(runner, *[arg for f in files for arg in ("--input", f)])

        assert result.exit_code == 0
        # Gzip is re-detected per file downstream, so the normalized family is passed on.
        assert classes_run.kwargs['format'] == expected_format

    def test_mixed_formats_are_rejected(self, runner, make_file, classes_run):
        csv = make_file("train.csv")
        fasta = make_file("train.fasta")

        result = invoke(runner, "--input", csv, "--input", fasta)

        assert result.exit_code == 1
        assert "All input files must have the same format" in result.stderr
        assert not classes_run.called

    def test_unsupported_extension_is_rejected(self, runner, make_file, classes_run):
        txt = make_file("train.txt", content="ACGT\n")

        result = invoke(runner, "--input", txt)

        assert result.exit_code == 1
        assert "Invalid format 'txt'" in result.stderr
        assert not classes_run.called

    @pytest.mark.parametrize("name", ["dataset", "results.v2/dataset"])
    def test_file_without_extension_is_rejected(self, runner, make_file, classes_run, name):
        # A dot in a parent directory does not stand in for a missing extension.
        dataset = make_file(name)

        result = invoke(runner, "--input", dataset)

        assert result.exit_code == 1
        assert "Invalid format ''" in result.stderr
        assert not classes_run.called

    def test_extension_is_read_from_the_file_not_the_directory(self, runner, make_file, classes_run):
        csv = make_file("results.v2/train.csv")

        result = invoke(runner, "--input", csv)

        assert result.exit_code == 0
        assert classes_run.kwargs['format'] == 'csv'


class TestFastaMinimumFiles:
    @pytest.mark.parametrize("name", ["only.fasta", "only.fa", "only.fasta.gz", "only.fa.gz"])
    def test_single_fasta_file_is_rejected(self, runner, make_file, classes_run, name):
        fasta = make_file(name)

        result = invoke(runner, "--input", fasta)

        assert result.exit_code == 1
        assert "at least 2 input files are required" in result.stderr
        assert not classes_run.called

    def test_two_fasta_files_are_accepted(self, runner, make_file, classes_run):
        first = make_file("pos.fasta")
        second = make_file("neg.fasta")

        result = invoke(runner, "--input", first, "--input", second)

        assert result.exit_code == 0
        assert classes_run.called

    @pytest.mark.parametrize("name", ["train.csv", "train.tsv"])
    def test_single_tabular_file_is_accepted(self, runner, make_file, classes_run, name):
        # The minimum of two files applies to fasta only: labels come from a column.
        tabular = make_file(name)

        result = invoke(runner, "--input", tabular)

        assert result.exit_code == 0
        assert classes_run.called


class TestValidation:
    def test_missing_input_file_is_reported(self, runner, tmp_path, classes_run):
        missing = str(tmp_path / "absent.csv")

        result = invoke(runner, "--input", missing)

        assert result.exit_code == 1
        assert f"Input file does not exist: {missing}" in result.stderr
        assert not classes_run.called

    def test_directory_is_not_accepted_as_input(self, runner, tmp_path, classes_run):
        directory = tmp_path / "inputs.csv"
        directory.mkdir()

        result = invoke(runner, "--input", str(directory))

        assert result.exit_code == 1
        assert "Input file does not exist" in result.stderr
        assert not classes_run.called

    def test_missing_file_is_reported_before_format_problems(self, runner, tmp_path, classes_run):
        # Existence is checked first, so an absent .txt file is reported as missing.
        missing = str(tmp_path / "absent.txt")

        result = invoke(runner, "--input", missing)

        assert result.exit_code == 1
        assert "does not exist" in result.stderr
        assert "Invalid format" not in result.stderr

    def test_invalid_report_type_is_rejected(self, runner, make_file, classes_run):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--report-types", "pdf")

        assert result.exit_code == 1
        assert "Invalid report type 'pdf'" in result.stderr
        assert not classes_run.called

    def test_invalid_report_type_alongside_valid_one_is_rejected(self, runner, make_file, classes_run):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--report-types", "html", "--report-types", "pdf")

        assert result.exit_code == 1
        assert "Invalid report type 'pdf'" in result.stderr
        assert not classes_run.called

    @pytest.mark.parametrize("report_type", ["json", "html", "simple"])
    def test_valid_report_types_are_accepted(self, runner, make_file, classes_run, report_type):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--report-types", report_type)

        assert result.exit_code == 0
        assert classes_run.kwargs['report_types'] == [report_type]

    def test_invalid_plot_type_is_rejected(self, runner, make_file, classes_run):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--plot-type", "scatter")

        assert result.exit_code == 1
        assert "Invalid plot type 'scatter'" in result.stderr
        assert not classes_run.called

    @pytest.mark.parametrize("plot_type", ["boxen", "violin"])
    def test_valid_plot_types_are_accepted(self, runner, make_file, classes_run, plot_type):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--plot-type", plot_type)

        assert result.exit_code == 0
        assert classes_run.kwargs['plot_type'] == plot_type

    def test_invalid_log_level_is_rejected(self, runner, make_file, classes_run):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--log-level", "TRACE")

        assert result.exit_code == 1
        assert "Invalid log level 'TRACE'" in result.stderr
        assert not classes_run.called

    def test_log_level_is_case_sensitive(self, runner, make_file, classes_run):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--log-level", "info")

        assert result.exit_code == 1
        assert "Invalid log level 'info'" in result.stderr

    @pytest.mark.parametrize("log_level", ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    def test_valid_log_levels_are_accepted(self, runner, make_file, classes_run, log_level):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--log-level", log_level)

        assert result.exit_code == 0
        assert classes_run.kwargs['log_level'] == log_level

    @pytest.mark.parametrize("end_position", ["0", "-1"])
    def test_non_positive_end_position_is_rejected(self, runner, make_file, classes_run, end_position):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--end-position", end_position)

        assert result.exit_code == 1
        assert "end_position must be a positive integer" in result.stderr
        assert not classes_run.called

    def test_positive_end_position_is_accepted(self, runner, make_file, classes_run):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--end-position", "1")

        assert result.exit_code == 0
        assert classes_run.kwargs['end_position'] == 1

    def test_non_integer_end_position_is_a_usage_error(self, runner, make_file, classes_run):
        csv = make_file("train.csv")

        result = invoke(runner, "--input", csv, "--end-position", "abc")

        assert result.exit_code == 2
        assert not classes_run.called


class TestPipelineFailures:
    def test_pipeline_exception_becomes_exit_code_1(self, runner, make_file, classes_run):
        csv = make_file("train.csv")
        classes_run.exception = RuntimeError("pipeline blew up")

        result = invoke(runner, "--input", csv)

        assert result.exit_code == 1
        assert "Error: pipeline blew up" in result.stderr

    def test_pipeline_exit_code_is_preserved(self, runner, make_file, classes_run):
        csv = make_file("train.csv")
        classes_run.exception = typer.Exit(code=4)

        result = invoke(runner, "--input", csv)

        assert result.exit_code == 4
