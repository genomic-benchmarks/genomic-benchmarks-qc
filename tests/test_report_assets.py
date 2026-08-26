"""Tests for the report's templates, stylesheets and scripts.

They live as files in `genomic_benchmarks_qc.report.assets` and are inlined into
every generated report, which puts three things at risk that nothing else would
catch: a file that is not declared as package data (works from a checkout,
missing from the wheel), a colour placeholder that never gets substituted, and -
now that the page templates are .html files too - markup that stopped being
well-formed while nobody was rendering it.
"""

import pathlib
import tomllib
from html.parser import HTMLParser

import pytest

from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.colors import FAIL_COLOR, PASS_COLOR, UNKNOWN_COLOR, WARN_COLOR

ASSET_DIR = pathlib.Path(assets.__file__).parent
ASSET_FILES = sorted(path.name for path in ASSET_DIR.iterdir()
                     if path.suffix in ('.html', '.css', '.js'))
TEMPLATES = [name for name in ASSET_FILES if name.endswith('.html')]
# A template is either a whole report page or a fragment dropped into one. Listed
# rather than sniffed, so that adding a template is a decision about which it is.
PAGE_TEMPLATES = ('classes_report.html', 'split_report_page.html')
FRAGMENT_TEMPLATES = ('split_results_table.html',)
# Elements that carry no closing tag, so an unbalanced-tag check must not wait
# for one. <p> and <li> may also be left open in HTML; the templates close both.
VOID = {'area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link',
        'meta', 'param', 'source', 'track', 'wbr'}
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
    css = assets.stylesheet('report.css', 'report_design.css')

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


class TagBalance(HTMLParser):
    """Every open tag has to be closed, and closed in the order it was opened.

    Not a validator - it is the one class of mistake that editing markup by hand
    actually produces, and the one a browser hides by silently repairing the
    page it was handed.
    """

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []
        self.errors = []

    def handle_starttag(self, tag, attrs):
        if tag not in VOID:
            self.stack.append((tag, self.getpos()[0]))

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if not self.stack:
            self.errors.append(f'</{tag}> at line {self.getpos()[0]} closes nothing')
        elif self.stack[-1][0] != tag:
            open_tag, line = self.stack[-1]
            self.errors.append(
                f'</{tag}> at line {self.getpos()[0]} closes <{open_tag}> '
                f'opened at line {line}')
            self.stack.pop()
        else:
            self.stack.pop()


def test_there_are_templates_to_check():
    assert TEMPLATES


@pytest.mark.parametrize('name', TEMPLATES)
def test_every_template_is_well_formed(name):
    parser = TagBalance()
    parser.feed(assets.template(name))

    assert parser.errors == []
    assert parser.stack == [], f'never closed: {parser.stack}'


def test_every_template_is_accounted_for():
    """A new .html file has to be classified, or the checks below skip it."""
    assert sorted(TEMPLATES) == sorted(PAGE_TEMPLATES + FRAGMENT_TEMPLATES)


@pytest.mark.parametrize('name', PAGE_TEMPLATES)
def test_every_page_template_is_a_whole_page(name):
    """A report is one self-contained file, so a page carries its own styling."""
    markup = assets.template(name)

    assert markup.startswith('<!DOCTYPE html>')
    assert '{{common_css}}' in markup          # the styling is inlined into it
    assert '{{report_header}}' in markup       # and the shared header dropped in
    assert '{{sidebar_links}}' in markup       # and the way back to the project


@pytest.mark.parametrize('name', FRAGMENT_TEMPLATES)
def test_a_fragment_is_not_a_page(name):
    """It is pasted into one, so it must not bring a second <html> with it."""
    markup = assets.template(name)

    assert '<!DOCTYPE' not in markup
    assert '<html' not in markup
