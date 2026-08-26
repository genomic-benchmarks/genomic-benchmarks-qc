"""Tests for the cross-references the docstrings carry.

The docstrings are written for mkdocstrings, so where a number comes from a
constant they link the constant instead of only repeating its value:
`[DEFAULT_MIN_COVERAGE][genomic_benchmarks_qc.utils.seq_stats.DEFAULT_MIN_COVERAGE]`.
Two things can rot there without anything else noticing - a target that was
renamed or moved, which the docs would render as plain bracketed text rather
than a link, and a documented default that drifted from the constant it cites,
which is what sent us looking in the first place.

Only references into this package are checked; a link to an external object is
left to the docs build.
"""

import ast
import importlib
import pathlib
import re

import pytest

import genomic_benchmarks_qc

SRC = pathlib.Path(genomic_benchmarks_qc.__file__).parent

# The mkdocstrings cross-reference form, scoped to this package so that prose
# like `matrix[i][j]` cannot look like a reference.
REFERENCE = re.compile(r'\[([A-Za-z_]\w*)\]\[(genomic_benchmarks_qc[\w.]*)\]')
# A documented default that also cites the constant it comes from, e.g.
# "Default: `0.25` ([DEFAULT_MIN_COVERAGE][...])".
DOCUMENTED_DEFAULT = re.compile(
    r'Default: `([^`]+)` \(\[(\w+)\]\[(genomic_benchmarks_qc[\w.]*)\]\)')

DOCUMENTED = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)


def _docstrings():
    """Every docstring in the package, whitespace collapsed onto one line."""
    for path in sorted(SRC.rglob('*.py')):
        for node in ast.walk(ast.parse(path.read_text())):
            if not isinstance(node, DOCUMENTED):
                continue
            text = ast.get_docstring(node)
            if text:
                where = f"{path.relative_to(SRC)}:{getattr(node, 'name', '<module>')}"
                yield where, ' '.join(text.split())


DOCSTRINGS = list(_docstrings())
REFERENCES = sorted({(where, text, target) for where, doc in DOCSTRINGS
                     for text, target in REFERENCE.findall(doc)})
DEFAULTS = sorted({(where, value, target) for where, doc in DOCSTRINGS
                   for value, _, target in DOCUMENTED_DEFAULT.findall(doc)})


def _resolve(target):
    """The object a dotted cross-reference target names.

    Raises ModuleNotFoundError or AttributeError if it names nothing, which is
    exactly the state these tests exist to catch.
    """
    try:
        return importlib.import_module(target)
    except ModuleNotFoundError:
        pass

    parts = target.split('.')
    for split in range(len(parts) - 1, 0, -1):
        try:
            obj = importlib.import_module('.'.join(parts[:split]))
        except ModuleNotFoundError:
            continue
        for attribute in parts[split:]:
            obj = getattr(obj, attribute)
        return obj
    raise ModuleNotFoundError(target)


def test_there_are_references_to_check():
    """Guards the tests below, which would all pass on a package with none."""
    assert REFERENCES
    assert DEFAULTS


@pytest.mark.parametrize('where, text, target', REFERENCES,
                         ids=[f'{where} -> {target}' for where, _, target in REFERENCES])
def test_every_reference_target_exists(where, text, target):
    try:
        _resolve(target)
    except (ModuleNotFoundError, AttributeError) as error:
        pytest.fail(f"{where} links to {target}, which no longer exists ({error})")


@pytest.mark.parametrize('where, text, target', REFERENCES,
                         ids=[f'{where} -> {target}' for where, _, target in REFERENCES])
def test_every_reference_shows_the_name_it_links_to(where, text, target):
    """A link reading one name and pointing at another misleads twice over."""
    assert text == target.rsplit('.', 1)[-1]


@pytest.mark.parametrize('where, value, target', DEFAULTS,
                         ids=[f'{where} -> {target}' for where, _, target in DEFAULTS])
def test_every_documented_default_matches_its_constant(where, value, target):
    """The drift this whole convention exists to prevent."""
    assert ast.literal_eval(value) == _resolve(target), (
        f"{where} documents {value}, but {target} is {_resolve(target)}")
