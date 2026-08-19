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
