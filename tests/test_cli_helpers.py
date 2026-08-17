"""Unit tests for the helper functions in genomic_benchmarks_qc.cli."""

import pytest
import typer

from genomic_benchmarks_qc.cli import (
    SPLIT_MEMORY_LIMIT_RE,
    _infer_format,
    _infer_format_from_inputs,
    _normalize_format,
    _run_command,
)


class TestInferFormat:
    @pytest.mark.parametrize(
        "file_path, expected",
        [
            ("data.csv", "csv"),
            ("data.tsv", "tsv"),
            ("data.fasta", "fasta"),
            ("data.fa", "fa"),
            ("data.csv.gz", "csv.gz"),
            ("data.tsv.gz", "tsv.gz"),
            ("data.fasta.gz", "fasta.gz"),
            ("data.fa.gz", "fa.gz"),
            # Extensions are matched case-insensitively.
            ("DATA.CSV", "csv"),
            ("Data.Fasta.GZ", "fasta.gz"),
            # Only the last extension counts.
            ("train.v2.csv", "csv"),
            ("folder/data.csv", "csv"),
            ("/abs/path/to/data.csv", "csv"),
            # Only the file name is considered, never the directories above it.
            ("results.v2/data.csv", "csv"),
            ("/abs/results.v2/data.fa.gz", "fa.gz"),
        ],
    )
    def test_infers_extension(self, file_path, expected):
        assert _infer_format(file_path) == expected

    @pytest.mark.parametrize(
        "file_path",
        ["dataset", "folder/dataset", "results.v2/dataset", "/abs/results.v2/dataset"],
    )
    def test_file_without_extension_has_empty_format(self, file_path):
        assert _infer_format(file_path) == ""

    def test_bare_gz_has_no_base_extension(self):
        assert _infer_format("dataset.gz") == ".gz"


class TestNormalizeFormat:
    @pytest.mark.parametrize(
        "fmt, expected",
        [
            # Gzip is detected per-file downstream, so it does not distinguish formats.
            ("csv", "csv"),
            ("csv.gz", "csv"),
            ("tsv", "tsv"),
            ("tsv.gz", "tsv"),
            # 'fa' is only an alias for 'fasta'.
            ("fasta", "fasta"),
            ("fasta.gz", "fasta"),
            ("fa", "fasta"),
            ("fa.gz", "fasta"),
            ("", ""),
        ],
    )
    def test_reduces_to_base_family(self, fmt, expected):
        assert _normalize_format(fmt) == expected


class TestInferFormatFromInputs:
    @pytest.mark.parametrize(
        "files, expected",
        [
            (["a.csv"], "csv"),
            (["a.csv", "b.csv"], "csv"),
            (["a.tsv", "b.tsv"], "tsv"),
            (["a.fasta", "b.fasta"], "fasta"),
            # Files that differ only in gzip or the fa/fasta alias may be mixed,
            # and the returned format is the normalized family.
            (["a.fa", "b.fasta"], "fasta"),
            (["a.fasta", "b.fasta.gz"], "fasta"),
            (["a.fa.gz", "b.fasta"], "fasta"),
            (["a.csv", "b.csv.gz"], "csv"),
        ],
    )
    def test_agreeing_inputs(self, files, expected):
        assert _infer_format_from_inputs(files) == expected

    @pytest.mark.parametrize(
        "files",
        [
            ["a.csv", "b.fasta"],
            ["a.csv", "b.tsv"],
            ["a.tsv", "b.fa.gz"],
            ["a.csv", "b.csv", "c.tsv"],
        ],
    )
    def test_conflicting_inputs_exit_with_error(self, files, capsys):
        with pytest.raises(typer.Exit) as excinfo:
            _infer_format_from_inputs(files)
        assert excinfo.value.exit_code == 1
        assert "All input files must have the same format" in capsys.readouterr().err

    def test_error_lists_every_file_with_its_own_extension(self, capsys):
        with pytest.raises(typer.Exit):
            _infer_format_from_inputs(["a.csv", "b.fasta.gz"])
        err = capsys.readouterr().err
        assert "a.csv (csv)" in err
        assert "b.fasta.gz (fasta.gz)" in err


class TestRunCommand:
    def test_forwards_positional_and_keyword_arguments(self):
        seen = {}

        def command(*args, **kwargs):
            seen['args'] = args
            seen['kwargs'] = kwargs

        _run_command(command, 1, 2, key='value')
        assert seen == {'args': (1, 2), 'kwargs': {'key': 'value'}}

    def test_converts_exception_into_exit_code_1(self, capsys):
        def command():
            raise ValueError("something broke")

        with pytest.raises(typer.Exit) as excinfo:
            _run_command(command)
        assert excinfo.value.exit_code == 1
        assert "Error: something broke" in capsys.readouterr().err

    def test_propagates_typer_exit_unchanged(self, capsys):
        def command():
            raise typer.Exit(code=3)

        with pytest.raises(typer.Exit) as excinfo:
            _run_command(command)
        assert excinfo.value.exit_code == 3
        # A deliberate exit is not reported as an unexpected error.
        assert capsys.readouterr().err == ""


class TestSplitMemoryLimitPattern:
    @pytest.mark.parametrize(
        "value",
        ["0", "1", "512", "1B", "7K", "100M", "10G", "1T", "2048G", "999999"],
    )
    def test_accepts_mmseqs2_byte_sizes(self, value):
        assert SPLIT_MEMORY_LIMIT_RE.match(value)

    @pytest.mark.parametrize(
        "value",
        [
            "",  # empty
            "G",  # unit without a number
            "10GB",  # only single-letter units
            "10Gi",  # no binary-prefix spelling
            "1.5G",  # integers only
            "-1",  # no sign
            "+1",
            "0G",  # zero takes no unit
            "01",  # no leading zeros
            "10 G",  # no inner whitespace
            "10G20",  # unit must be last
            "abc",
            "1e3",
            "10G\n",
        ],
    )
    def test_rejects_everything_else(self, value):
        assert SPLIT_MEMORY_LIMIT_RE.match(value) is None
