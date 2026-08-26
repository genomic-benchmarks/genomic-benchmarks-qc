"""Shared fixtures for the test suite."""

import pytest
from typer.testing import CliRunner

from genomic_benchmarks_qc import cli

CSV_CONTENT = "sequence,label\nACGT,0\nTGCA,1\n"
TSV_CONTENT = "sequence\tlabel\nACGT\t0\nTGCA\t1\n"
FASTA_CONTENT = ">seq1\nACGT\n>seq2\nTGCA\n"


class RunRecorder:
    """Stand-in for a pipeline entry point, recording how the CLI called it.

    Set ``exception`` to make the stubbed pipeline fail, which is how the tests
    exercise the error handling in ``cli._run_command``.
    """

    def __init__(self):
        self.calls = []
        self.exception = None

    def __call__(self, **kwargs):
        self.calls.append(kwargs)
        if self.exception is not None:
            raise self.exception

    @property
    def called(self):
        return bool(self.calls)

    @property
    def kwargs(self):
        """Keyword arguments of the single expected call."""
        assert len(self.calls) == 1, f"expected exactly one call, got {len(self.calls)}"
        return self.calls[0]


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def make_file(tmp_path):
    """Create an input file and return its path as a string.

    The CLI only ever inspects file names and existence — the pipelines that
    read the contents are stubbed out — so placeholder content matching the
    extension is enough. ``.gz`` files are therefore written uncompressed.
    """

    def _make_file(name, content=None):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if content is None:
            base = name[: -len('.gz')] if name.endswith('.gz') else name
            if base.endswith(('.csv',)):
                content = CSV_CONTENT
            elif base.endswith(('.tsv',)):
                content = TSV_CONTENT
            else:
                content = FASTA_CONTENT
        path.write_text(content)
        return str(path)

    return _make_file


@pytest.fixture
def classes_run(monkeypatch):
    """Replace the evaluate-classes pipeline with a recorder."""
    recorder = RunRecorder()
    monkeypatch.setattr(cli, 'run_evaluate_classes', recorder)
    return recorder


@pytest.fixture
def splits_run(monkeypatch):
    """Replace the evaluate-splits pipeline with a recorder."""
    recorder = RunRecorder()
    monkeypatch.setattr(cli, 'run_evaluate_splits', recorder)
    return recorder
