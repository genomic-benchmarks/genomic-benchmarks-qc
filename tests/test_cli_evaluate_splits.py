"""Tests for the `gb-qc evaluate-splits` command.

The pipeline itself is stubbed out (see the ``splits_run`` fixture), so these
tests cover what the CLI layer is responsible for: argument parsing, validation
and the arguments handed to the pipeline.
"""

import pytest
import typer

from genomic_benchmarks_qc.cli import app


def invoke(runner, *args):
    return runner.invoke(app, ["evaluate-splits", *args])


@pytest.fixture
def split_files(make_file):
    """A minimal valid train/test pair of CSV inputs."""
    return make_file("enhancers_train.csv"), make_file("enhancers_test.csv")


class TestArgumentPassing:
    def test_defaults_are_forwarded(self, runner, split_files, splits_run):
        train, test = split_files

        result = invoke(runner, "--train-input", train, "--test-input", test)

        assert result.exit_code == 0
        assert splits_run.kwargs == {
            'train_files': [train],
            'test_files': [test],
            'format': 'csv',
            'out_folder': '.',
            'sequence_column': ['sequence'],
            'report_types': ['html', 'simple'],
            'similarity_threshold': 90.0,
            'threads': None,
            'split_memory_limit': None,
            'keep_tmp_files': False,
            'log_level': 'INFO',
            'log_file': None,
        }

    def test_explicit_options_are_forwarded(self, runner, split_files, splits_run, tmp_path):
        train, test = split_files
        out = str(tmp_path / "out")
        log = str(tmp_path / "run.log")

        result = invoke(
            runner,
            "--train-input", train,
            "--test-input", test,
            "--sequence-column", "seq_a",
            "--sequence-column", "seq_b",
            "--out-folder", out,
            "--report-types", "simple",
            "--similarity-threshold", "80.5",
            "--threads", "8",
            "--split-memory-limit", "10G",
            "--keep-tmp-files",
            "--log-level", "WARNING",
            "--log-file", log,
        )

        assert result.exit_code == 0
        assert splits_run.kwargs == {
            'train_files': [train],
            'test_files': [test],
            'format': 'csv',
            'out_folder': out,
            'sequence_column': ['seq_a', 'seq_b'],
            'report_types': ['simple'],
            'similarity_threshold': 80.5,
            'threads': 8,
            'split_memory_limit': '10G',
            'keep_tmp_files': True,
            'log_level': 'WARNING',
            'log_file': log,
        }

    def test_multiple_inputs_per_side_are_forwarded_in_order(self, runner, make_file, splits_run):
        train_a, train_b = make_file("train_a.csv"), make_file("train_b.csv")
        test_a, test_b = make_file("test_a.csv"), make_file("test_b.csv")

        result = invoke(
            runner,
            "--train-input", train_a, "--train-input", train_b,
            "--test-input", test_a, "--test-input", test_b,
        )

        assert result.exit_code == 0
        assert splits_run.kwargs['train_files'] == [train_a, train_b]
        assert splits_run.kwargs['test_files'] == [test_a, test_b]

    @pytest.mark.parametrize(
        "flag, expected", [("--keep-tmp-files", True), ("--no-keep-tmp-files", False)]
    )
    def test_keep_tmp_files_flag(self, runner, split_files, splits_run, flag, expected):
        train, test = split_files

        result = invoke(runner, "--train-input", train, "--test-input", test, flag)

        assert result.exit_code == 0
        assert splits_run.kwargs['keep_tmp_files'] is expected

    @pytest.mark.parametrize(
        "args",
        [
            [],
            ["--train-input", "train.csv"],
            ["--test-input", "test.csv"],
        ],
        ids=["neither", "train_only", "test_only"],
    )
    def test_both_inputs_are_required(self, runner, splits_run, args):
        result = invoke(runner, *args)

        assert result.exit_code == 2
        assert not splits_run.called


class TestFormatHandling:
    @pytest.mark.parametrize(
        "train_name, test_name, expected_format",
        [
            ("train.csv", "test.csv", 'csv'),
            ("train.csv.gz", "test.csv.gz", 'csv'),
            ("train.tsv", "test.tsv", 'tsv'),
            ("train.fasta", "test.fasta", 'fasta'),
            # Inputs may differ in gzip or the fa/fasta alias.
            ("train.csv", "test.csv.gz", 'csv'),
            ("train.fa", "test.fasta", 'fasta'),
            ("train.fa.gz", "test.fasta", 'fasta'),
        ],
    )
    def test_format_is_inferred_across_both_sides(
        self, runner, make_file, splits_run, train_name, test_name, expected_format
    ):
        train, test = make_file(train_name), make_file(test_name)

        result = invoke(runner, "--train-input", train, "--test-input", test)

        assert result.exit_code == 0
        assert splits_run.kwargs['format'] == expected_format

    def test_single_fasta_file_per_side_is_accepted(self, runner, make_file, splits_run):
        # Unlike evaluate-classes, splits carry no per-file class labels,
        # so one fasta file per side is enough.
        train, test = make_file("train.fasta"), make_file("test.fasta")

        result = invoke(runner, "--train-input", train, "--test-input", test)

        assert result.exit_code == 0
        assert splits_run.called

    def test_format_mismatch_between_train_and_test_is_rejected(self, runner, make_file, splits_run):
        train, test = make_file("train.csv"), make_file("test.fasta")

        result = invoke(runner, "--train-input", train, "--test-input", test)

        assert result.exit_code == 1
        assert "All input files must have the same format" in result.stderr
        assert not splits_run.called

    def test_format_mismatch_within_one_side_is_rejected(self, runner, make_file, splits_run):
        train_a, train_b = make_file("train_a.csv"), make_file("train_b.tsv")
        test = make_file("test.csv")

        result = invoke(
            runner, "--train-input", train_a, "--train-input", train_b, "--test-input", test
        )

        assert result.exit_code == 1
        assert "All input files must have the same format" in result.stderr
        assert not splits_run.called

    def test_unsupported_extension_is_rejected(self, runner, make_file, splits_run):
        train = make_file("train.txt", content="ACGT\n")
        test = make_file("test.txt", content="ACGT\n")

        result = invoke(runner, "--train-input", train, "--test-input", test)

        assert result.exit_code == 1
        assert "Invalid format 'txt'" in result.stderr
        assert not splits_run.called


class TestInputExistence:
    def test_missing_train_file_is_reported(self, runner, make_file, tmp_path, splits_run):
        missing = str(tmp_path / "absent_train.csv")
        test = make_file("test.csv")

        result = invoke(runner, "--train-input", missing, "--test-input", test)

        assert result.exit_code == 1
        assert f"Training input file does not exist: {missing}" in result.stderr
        assert not splits_run.called

    def test_missing_test_file_is_reported(self, runner, make_file, tmp_path, splits_run):
        train = make_file("train.csv")
        missing = str(tmp_path / "absent_test.csv")

        result = invoke(runner, "--train-input", train, "--test-input", missing)

        assert result.exit_code == 1
        assert f"Test input file does not exist: {missing}" in result.stderr
        assert not splits_run.called

    def test_train_files_are_checked_before_test_files(self, runner, tmp_path, splits_run):
        missing_train = str(tmp_path / "absent_train.csv")
        missing_test = str(tmp_path / "absent_test.csv")

        result = invoke(runner, "--train-input", missing_train, "--test-input", missing_test)

        assert result.exit_code == 1
        assert "Training input file does not exist" in result.stderr
        assert "Test input file does not exist" not in result.stderr


class TestValidation:
    @pytest.mark.parametrize("threshold", ["0", "0.0", "50", "100", "100.0"])
    def test_threshold_within_range_is_accepted(self, runner, split_files, splits_run, threshold):
        train, test = split_files

        result = invoke(
            runner, "--train-input", train, "--test-input", test,
            "--similarity-threshold", threshold,
        )

        assert result.exit_code == 0
        assert splits_run.kwargs['similarity_threshold'] == float(threshold)

    @pytest.mark.parametrize("threshold", ["-1", "-0.5", "100.1", "1000"])
    def test_threshold_out_of_range_is_rejected(self, runner, split_files, splits_run, threshold):
        train, test = split_files

        result = invoke(
            runner, "--train-input", train, "--test-input", test,
            "--similarity-threshold", threshold,
        )

        assert result.exit_code == 1
        assert "similarity_threshold must be between 0 and 100" in result.stderr
        assert not splits_run.called

    def test_non_numeric_threshold_is_a_usage_error(self, runner, split_files, splits_run):
        train, test = split_files

        result = invoke(
            runner, "--train-input", train, "--test-input", test,
            "--similarity-threshold", "high",
        )

        assert result.exit_code == 2
        assert not splits_run.called

    @pytest.mark.parametrize("threads", ["1", "4", "64"])
    def test_positive_thread_count_is_accepted(self, runner, split_files, splits_run, threads):
        train, test = split_files

        result = invoke(runner, "--train-input", train, "--test-input", test, "--threads", threads)

        assert result.exit_code == 0
        assert splits_run.kwargs['threads'] == int(threads)

    @pytest.mark.parametrize("threads", ["0", "-1"])
    def test_non_positive_thread_count_is_rejected(self, runner, split_files, splits_run, threads):
        train, test = split_files

        result = invoke(runner, "--train-input", train, "--test-input", test, "--threads", threads)

        assert result.exit_code == 1
        assert "threads must be a positive integer" in result.stderr
        assert not splits_run.called

    def test_invalid_report_type_is_rejected(self, runner, split_files, splits_run):
        train, test = split_files

        result = invoke(
            runner, "--train-input", train, "--test-input", test, "--report-types", "pdf"
        )

        assert result.exit_code == 1
        assert "Invalid report type 'pdf'" in result.stderr
        assert not splits_run.called

    @pytest.mark.parametrize("report_type", ["html", "simple"])
    def test_valid_report_types_are_accepted(self, runner, split_files, splits_run, report_type):
        train, test = split_files

        result = invoke(
            runner, "--train-input", train, "--test-input", test, "--report-types", report_type
        )

        assert result.exit_code == 0
        assert splits_run.kwargs['report_types'] == [report_type]

    def test_a_report_type_this_command_cannot_write_is_rejected(
            self, runner, split_files, splits_run):
        """`json` is the classes command's per-class statistics; a search that
        reports on the split as a whole has nothing to put in one. Accepting it
        here meant a full MMseqs2 run that wrote no file and exited 0."""
        train, test = split_files

        result = invoke(
            runner, "--train-input", train, "--test-input", test, "--report-types", "json"
        )

        assert result.exit_code == 1
        assert "Invalid report type 'json'" in result.stderr
        assert "html, simple" in result.stderr
        assert not splits_run.called

    def test_invalid_log_level_is_rejected(self, runner, split_files, splits_run):
        train, test = split_files

        result = invoke(runner, "--train-input", train, "--test-input", test, "--log-level", "TRACE")

        assert result.exit_code == 1
        assert "Invalid log level 'TRACE'" in result.stderr
        assert not splits_run.called

    @pytest.mark.parametrize("log_level", ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    def test_valid_log_levels_are_accepted(self, runner, split_files, splits_run, log_level):
        train, test = split_files

        result = invoke(
            runner, "--train-input", train, "--test-input", test, "--log-level", log_level
        )

        assert result.exit_code == 0
        assert splits_run.kwargs['log_level'] == log_level


class TestSplitMemoryLimit:
    @pytest.mark.parametrize("value", ["0", "512", "1B", "7K", "100M", "10G", "1T"])
    def test_valid_byte_sizes_are_forwarded(self, runner, split_files, splits_run, value):
        train, test = split_files

        result = invoke(
            runner, "--train-input", train, "--test-input", test, "--split-memory-limit", value
        )

        assert result.exit_code == 0
        assert splits_run.kwargs['split_memory_limit'] == value

    @pytest.mark.parametrize(
        "value",
        ["10GB", "10Gi", "1.5G", "-1", "0G", "abc", "G", "", "10 G", "1e3"],
    )
    def test_invalid_byte_sizes_are_rejected(self, runner, split_files, splits_run, value):
        train, test = split_files

        result = invoke(
            runner, "--train-input", train, "--test-input", test, "--split-memory-limit", value
        )

        assert result.exit_code == 1
        assert "split_memory_limit must be 0 or a positive integer" in result.stderr
        # The value never reaches MMseqs2, where it would only fail mid-run.
        assert not splits_run.called

    @pytest.mark.parametrize(
        "value, expected", [("10g", "10G"), ("1t", "1T"), (" 10G ", "10G"), ("10G\n", "10G")]
    )
    def test_case_and_surrounding_whitespace_are_normalized(
        self, runner, split_files, splits_run, value, expected
    ):
        train, test = split_files

        result = invoke(
            runner, "--train-input", train, "--test-input", test, "--split-memory-limit", value
        )

        assert result.exit_code == 0
        assert splits_run.kwargs['split_memory_limit'] == expected

    def test_omitting_the_option_leaves_mmseqs2_defaults_in_place(self, runner, split_files, splits_run):
        train, test = split_files

        result = invoke(runner, "--train-input", train, "--test-input", test)

        assert result.exit_code == 0
        assert splits_run.kwargs['split_memory_limit'] is None


class TestPipelineFailures:
    def test_pipeline_exception_becomes_exit_code_1(self, runner, split_files, splits_run):
        train, test = split_files
        splits_run.exception = RuntimeError("MMSeqs2 executable not found in PATH.")

        result = invoke(runner, "--train-input", train, "--test-input", test)

        assert result.exit_code == 1
        assert "Error: MMSeqs2 executable not found in PATH." in result.stderr

    def test_pipeline_exit_code_is_preserved(self, runner, split_files, splits_run):
        train, test = split_files
        splits_run.exception = typer.Exit(code=4)

        result = invoke(runner, "--train-input", train, "--test-input", test)

        assert result.exit_code == 4
