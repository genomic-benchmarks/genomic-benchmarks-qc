"""Tests for the report's stylesheets and scripts.

They live as files in `genomic_benchmarks_qc.report.assets` and are inlined into
every generated report, which puts two things at risk that nothing else would
catch: a file that is not declared as package data (works from a checkout,
missing from the wheel), and a colour placeholder that never gets substituted.
"""

import pathlib
import tomllib

import pytest

from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.colors import FAIL_COLOR, PASS_COLOR, UNKNOWN_COLOR, WARN_COLOR

ASSET_DIR = pathlib.Path(assets.__file__).parent
ASSET_FILES = sorted(path.name for path in ASSET_DIR.iterdir()
                     if path.suffix in ('.css', '.js'))
PYPROJECT = pathlib.Path(__file__).resolve().parents[1] / 'pyproject.toml'


def test_there_are_assets_to_check():
    """Guards the tests below, which would all pass on an empty directory."""
    assert ASSET_FILES


@pytest.mark.parametrize('name', ASSET_FILES)
def test_every_asset_loads_through_the_loader(name):
    """Reading through importlib.resources is how the installed package sees them."""
    assert assets.read_asset(name).strip()


@pytest.mark.parametrize('name', ASSET_FILES)
def test_every_asset_is_declared_as_package_data(name):
    """A file added here and not to pyproject.toml is missing from the wheel."""
    config = tomllib.loads(PYPROJECT.read_text())
    globs = config['tool']['setuptools']['package-data']['genomic_benchmarks_qc.report.assets']
    assert any(pathlib.PurePath(name).match(glob) for glob in globs)


def test_the_stylesheets_carry_the_flag_colors():
    css = assets.stylesheet('report.css')

    for color in (FAIL_COLOR, WARN_COLOR, PASS_COLOR, UNKNOWN_COLOR):
        assert color in css


@pytest.mark.parametrize('name', [name for name in ASSET_FILES if name.endswith('.css')])
def test_no_placeholder_survives_in_a_stylesheet(name):
    """The colours are the only templating the stylesheets get; none may be left."""
    assert '{{' not in assets.read_asset(name)


def test_scripts_are_wrapped_one_element_each():
    markup = assets.script('report_ui.js', 'split_report.js')

    assert markup.count('<script>') == 2
    assert markup.count('</script>') == 2


def test_a_missing_asset_is_an_error_rather_than_an_empty_page():
    with pytest.raises(FileNotFoundError):
        assets.read_asset('not_a_real_asset.css')
