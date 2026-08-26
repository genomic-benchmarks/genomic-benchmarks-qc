"""Tests for the single place the version is written down.

The reports print the version, and the built package declares it. Those used to
be two sources - a literal in pyproject.toml, read back at runtime through
importlib.metadata - which cost a metadata scan per report and made report
generation fail outright wherever the package was not installed, e.g. from a
plain checkout or a fresh git worktree. Now pyproject reads the package
attribute, and these tests pin that arrangement so it cannot quietly revert.
"""

import pathlib
import re
import tomllib

import genomic_benchmarks_qc
from genomic_benchmarks_qc import cli

PYPROJECT = pathlib.Path(__file__).resolve().parents[1] / 'pyproject.toml'
CONFIG = tomllib.loads(PYPROJECT.read_text())


def test_the_package_states_its_version():
    assert re.fullmatch(r'\d+\.\d+\.\d+(\.\w+)?', genomic_benchmarks_qc.__version__)


def test_pyproject_takes_the_version_from_the_package():
    """A literal version here again would be a second source to keep in step."""
    assert 'version' in CONFIG['project']['dynamic']
    assert 'version' not in CONFIG['project']
    assert CONFIG['tool']['setuptools']['dynamic']['version'] == {
        'attr': 'genomic_benchmarks_qc.__version__'
    }


def test_the_reports_do_not_read_installed_metadata():
    """Reading it back through importlib.metadata is both slower and fallible."""
    report_dir = pathlib.Path(genomic_benchmarks_qc.__file__).parent / 'report'

    offenders = [path.name for path in report_dir.glob('*.py')
                 if 'importlib.metadata' in path.read_text()]

    assert offenders == []


def test_the_cli_reports_the_version(runner):
    """`gb-qc --version` is how a reader of a report finds out what wrote it.

    It also has to answer without a command after it, which is what `is_eager`
    on the option buys - `gb-qc --version` alone would otherwise be a usage
    error.
    """
    result = runner.invoke(cli.app, ['--version'])

    assert result.exit_code == 0, result.output
    assert genomic_benchmarks_qc.__version__ in result.output
