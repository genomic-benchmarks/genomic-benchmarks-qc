"""Tests that the README's option tables still describe the actual CLI.

The README documents every `gb-qc` option in a table, and the docs site
generates the same information from the Typer app. Two copies of a CLI reference
is exactly the kind of thing that drifts - an option gets added, or a default
changes, and the hand-written table quietly becomes wrong while nothing fails.

So the tables are checked against the app itself. Names are compared strictly:
an option that exists must be documented, and a documented option must exist.
Defaults are compared where the README states one as a literal, which is most of
them; where it uses a word instead (`required`, `auto`, `none`, `unlimited`) the
check is looser, because "auto" is a better thing to tell a reader than `None`
and the point is that the table is honest, not that it is machine-generated.
"""

import re

import pytest
from typer.main import get_command

from genomic_benchmarks_qc.cli import app

README = 'README.md'

# The README wraps each command's table in a <details> block whose summary names
# the command, e.g. "<summary><b><code>evaluate-classes</code> options</b>".
TABLE_BLOCK = re.compile(
    r'<summary><b><code>(?P<command>[\w-]+)</code> options</b></summary>'
    r'(?P<body>.*?)</details>',
    re.DOTALL)
# | `--option` | default | description |
TABLE_ROW = re.compile(r'^\|\s*`(?P<option>--[\w-]+)`\s*\|(?P<default>[^|]*)\|')

# Words the README uses in place of a literal default. Each says something true
# that the raw value does not, so they are accepted rather than rewritten.
PROSE_DEFAULTS = {
    'required': None,     # checked against the option's `required` instead
    'auto': None,         # the tool picks a value at runtime
    'none': None,         # unset means the feature is off
    'unlimited': None,    # unset means no limit is applied
}


def _readme_tables() -> dict[str, dict[str, str]]:
    """Every documented option and its stated default, keyed by command."""
    with open(README) as handle:
        text = handle.read()

    tables = {}
    for match in TABLE_BLOCK.finditer(text):
        rows = {}
        for line in match.group('body').splitlines():
            row = TABLE_ROW.match(line.strip())
            if row:
                rows[row.group('option')] = row.group('default').strip()
        tables[match.group('command')] = rows
    return tables


def _cli_options() -> dict[str, dict[str, object]]:
    """Every real option and its click default, keyed by command."""
    command = get_command(app)
    options = {}
    for name, sub in command.commands.items():
        options[name] = {
            param.opts[0]: param
            for param in sub.params
            if param.param_type_name == 'option'
        }
    return options


README_TABLES = _readme_tables()
CLI_OPTIONS = _cli_options()
COMMANDS = sorted(CLI_OPTIONS)


def _normalize(stated: str):
    """The README's stated default as a comparable value, or None for prose.

    A list default is written space-separated (`html simple`), which is how a
    reader would type it, not how Python spells it.
    """
    bare = stated.strip().strip('`').strip()
    if bare.lower() in PROSE_DEFAULTS:
        return None
    return bare


def _spell(default) -> set[str]:
    """How the README might reasonably write a click default."""
    if isinstance(default, (list, tuple)):
        return {' '.join(str(v) for v in default), ', '.join(str(v) for v in default)}
    if isinstance(default, bool):
        return {str(default), str(default).lower()}
    if isinstance(default, float):
        # 90.0 may reasonably be written 90.0 or 90
        return {str(default), str(int(default))} if default.is_integer() else {str(default)}
    return {str(default)}


def test_every_command_has_a_table():
    """A new command must bring its own README table."""
    assert set(README_TABLES) == set(CLI_OPTIONS), (
        f"README documents {sorted(README_TABLES)}, "
        f"gb-qc has {sorted(CLI_OPTIONS)}")


@pytest.mark.parametrize('command', COMMANDS)
def test_no_undocumented_options(command):
    """Every option the command accepts appears in its README table."""
    missing = sorted(set(CLI_OPTIONS[command]) - set(README_TABLES[command]))
    assert not missing, (
        f"`gb-qc {command}` accepts options the README does not document: "
        f"{', '.join(missing)}")


@pytest.mark.parametrize('command', COMMANDS)
def test_no_invented_options(command):
    """Every option the README documents is one the command accepts."""
    extra = sorted(set(README_TABLES[command]) - set(CLI_OPTIONS[command]))
    assert not extra, (
        f"README documents options `gb-qc {command}` does not accept: "
        f"{', '.join(extra)}")


@pytest.mark.parametrize('command', COMMANDS)
def test_required_options_say_required(command):
    """`required` in the table means required in the app, and vice versa."""
    for option, stated in README_TABLES[command].items():
        param = CLI_OPTIONS[command].get(option)
        if param is None:
            continue  # reported by test_no_invented_options
        says_required = stated.strip().strip('`').lower() == 'required'
        assert says_required == bool(param.required), (
            f"`gb-qc {command} {option}`: README says "
            f"{'required' if says_required else stated!r}, "
            f"app says required={param.required}")


@pytest.mark.parametrize('command', COMMANDS)
def test_stated_defaults_match(command):
    """Where the README states a literal default, it is the real one."""
    for option, stated in README_TABLES[command].items():
        param = CLI_OPTIONS[command].get(option)
        if param is None:
            continue
        expected = _normalize(stated)
        if expected is None:
            # Prose default. Only meaningful claim: the real default is falsy,
            # i.e. there genuinely is nothing to state.
            assert not param.default or param.required, (
                f"`gb-qc {command} {option}`: README says {stated!r}, but the "
                f"default is {param.default!r}, which is worth stating")
            continue
        assert expected in _spell(param.default), (
            f"`gb-qc {command} {option}`: README says default {expected!r}, "
            f"app default is {param.default!r}")
