#!/usr/bin/env python3
"""Assemble the generated parts of the documentation site, then build it.

Four things are generated rather than committed, so that none of them can drift
from the code they describe:

    docs/reference/cli.md              the CLI reference, from the Typer app
    docs/reports/<name>/               the example reports, from examples/build.py
    docs/assets/report-screenshot.png  the landing page's shot of one of those
    site/reference/api/                the API reference, by mkdocstrings

Usage:

    python scripts/build_docs.py            # everything, then mkdocs build
    python scripts/build_docs.py --serve     # ... then mkdocs serve instead
    python scripts/build_docs.py --skip-reports    # 2.5 minutes faster

`--skip-reports` is for working on prose: the example reports take about two and
a half minutes to regenerate and rarely change while you are writing. The
published build never skips them.

Building the reports needs MMseqs2 on PATH for the four `evaluate-splits` runs.
Both workflows in .github/workflows/ install it; locally, put its `bin` on PATH
first or the split reports fail and the build stops.

Two traps when running this locally, both of which produce a site that looks
freshly built and carries pre-edit code. `check_installed_matches_tree` refuses
to build under either.

The first: uv caches the wheel it
builds from the source tree and will reuse it while the version is unchanged, so
the `gb-qc` this calls - and the package mkdocstrings imports - can be older than
the working tree. Reports come out looking freshly built and carrying stale
strings. `--refresh-package genomic-benchmarks-qc` is the documented fix and
does work - but not alongside `--with-requirements`, which makes uv reuse the
cached build and ignore both --refresh and --refresh-package. Pass the docs
dependencies as individual `--with` flags, or run this with an interpreter that
has the package installed editable. CI is unaffected: it does a plain
`pip install .` on a clean runner.

The second: setuptools reuses `build/lib/` between builds and only copies a
source file over when it looks newer, so a file it decides is unchanged stays at
whatever version it was when that directory was first populated. Editing an
asset that is package data - the report's CSS and JavaScript - and rebuilding
can therefore ship the old one. `rm -rf build/lib` fixes it. CI is again
unaffected, having no `build/` at all.
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / 'docs'
CLI_PAGE = DOCS / 'reference' / 'cli.md'
REPORTS = DOCS / 'reports'

CLI_HEADER = """# CLI reference

Everything below is generated from the `gb-qc` application itself, so it cannot
fall behind the code. `gb-qc evaluate-classes --help` prints the same thing in
your terminal.

The README carries the same options as annotated tables, with notes on when to
reach for each one; `tests/test_readme_cli_tables.py` checks the two agree.

"""


def _run(argv: list[str], **kwargs) -> subprocess.CompletedProcess:
    """Run a command from the repository root, failing loudly."""
    print(f"$ {' '.join(argv)}", flush=True)
    result = subprocess.run(argv, cwd=ROOT, **kwargs)
    if result.returncode != 0:
        raise SystemExit(f"failed ({result.returncode}): {' '.join(argv)}")
    return result


def check_installed_matches_tree():
    """Refuse to build if the importable package differs from src/.

    Every generated part of this site comes from the installed package rather
    than from the working tree - the reports from running `gb-qc`, the API pages
    from mkdocstrings importing it - so if the two disagree, the build succeeds
    while publishing something the code no longer does. That has happened twice
    here from two unrelated causes: uv reusing a wheel it had cached from an
    earlier state of the tree, and setuptools reusing a `build/lib` whose copy of
    an asset it had decided was current. Rather than test for each cause, this
    tests the property both of them break.

    An editable install points at src/ and is always in agreement, so this is a
    no-op there. It compares the package's files, which covers the report's CSS
    and JavaScript - package data, and the likeliest thing to go stale unnoticed,
    since nothing about a report looks wrong when its script is a version old.

    The limit worth knowing: this checks the package that `sys.executable` can
    import, which is also where mkdocstrings reads from and, in every
    environment used here, where the `gb-qc` on PATH comes from. It cannot see a
    `gb-qc` resolved from some other environment.
    """
    import importlib.util

    spec = importlib.util.find_spec('genomic_benchmarks_qc')
    if spec is None or not spec.origin:
        raise SystemExit(
            "genomic_benchmarks_qc is not importable, so the reports cannot be "
            "built. Install it: `pip install -e .`"
        )
    installed = Path(spec.origin).resolve().parent
    source = (ROOT / 'src' / 'genomic_benchmarks_qc').resolve()
    if installed == source:
        return  # editable install: the same files, by definition

    differing = [
        path.relative_to(source) for path in sorted(source.rglob('*'))
        if path.is_file() and path.suffix not in {'.pyc'}
        and '__pycache__' not in path.parts
        and (installed / path.relative_to(source)).is_file()
        and (installed / path.relative_to(source)).read_bytes() != path.read_bytes()
    ]
    if differing:
        raise SystemExit(
            "the installed package does not match src/, so this build would "
            "publish code the tree no longer contains:\n  "
            + '\n  '.join(str(path) for path in differing)
            + f"\n\ninstalled: {installed}"
            + "\n\nReinstall it. Under uvx, note that `--with-requirements` "
              "defeats both --refresh and --refresh-package: pass the docs "
              "dependencies as individual `--with` flags, or run this from an "
              "environment with the package installed editable."
        )


def generate_cli_page():
    """Write docs/reference/cli.md from the Typer app.

    `typer ... utils docs` renders the whole command tree as Markdown. Going
    through the installed console script rather than importing the app keeps
    this honest about what a user actually gets.
    """
    result = _run(['typer', 'genomic_benchmarks_qc.cli', 'utils', 'docs',
                   '--name', 'gb-qc'], capture_output=True, text=True)
    body = result.stdout
    if 'evaluate-classes' not in body:
        raise SystemExit("generated CLI docs look wrong: no evaluate-classes section")

    # Typer emits `# gb-qc` as the top heading; the page supplies its own, so
    # demote the command tree by one level to keep a single h1 per page.
    lines = []
    for line in body.splitlines():
        if line.startswith('# `gb-qc`'):
            continue
        lines.append('#' + line if line.startswith('#') else line)

    CLI_PAGE.parent.mkdir(parents=True, exist_ok=True)
    CLI_PAGE.write_text(CLI_HEADER + '\n'.join(lines).strip() + '\n')
    print(f"wrote {CLI_PAGE.relative_to(ROOT)} ({len(lines)} lines)")


def generate_reports():
    """Build every example's reports and copy them under docs/reports/.

    They land at reports/<example>/<class|split>/<column>/<comparison>/, which
    is `gb-qc`'s own output layout below the example name - so the URL a reader
    is looking at is the same shape as the directory the tool leaves on their
    disk.

    Note this is not docs/examples/: a report directory there would collide with
    the URL mkdocs gives the example's own page.
    """
    scratch = ROOT / 'build' / 'example-reports'
    _run([sys.executable, 'examples/build.py', '--out-folder', str(scratch), '--check'])

    if REPORTS.exists():
        shutil.rmtree(REPORTS)
    REPORTS.mkdir(parents=True)

    copied = 0
    for report in sorted(scratch.rglob('gb-qc-report.html')):
        target = REPORTS / report.relative_to(scratch)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(report, target)
        # The CSV and plots sit beside it and are worth publishing too: the CSV
        # is what a reader would put in CI, and the report links its own plots.
        for extra in ('gb-qc-report.csv',):
            source = report.parent / extra
            if source.is_file():
                shutil.copy2(source, target.parent / extra)
        plots = report.parent / 'plots'
        if plots.is_dir():
            shutil.copytree(plots, target.parent / 'plots', dirs_exist_ok=True)
        copied += 1
    print(f"copied {copied} report(s) into {REPORTS.relative_to(ROOT)}")


GENERATED = DOCS / '_generated'

FLAG_SPANS = {
    'Pass': '<span class="flag flag-pass">Pass</span>',
    'Warning': '<span class="flag flag-warn">Warning</span>',
    'Fail': '<span class="flag flag-fail">Fail</span>',
    'Unknown': '<span class="flag flag-unknown">Unknown</span>',
}


def _flag_rows(report_csv: Path) -> list[tuple[str, str, str]]:
    """(check, flag, figure) for a report's headline checks.

    The figure column is whatever that report measures: AU-ROC for the class
    checks, leaked percentages for a split report, and nothing for the checks
    that are a yes/no rather than a score.
    """
    rows = []
    with report_csv.open(newline='') as handle:
        for row in csv.DictReader(handle):
            if ' - ' in row['Check']:
                continue  # per-base or per-position breakdown
            leaked = row.get('Percentage of leaked queries')
            if leaked:
                figure = f"{leaked} of queries, {row['Percentage of leaked targets']} of targets"
            elif row.get('AU-ROC'):
                figure = f"{float(row['AU-ROC']):.3f}"
            else:
                figure = '—'
            rows.append((row['Check'], row['Flag'], figure))
    return rows


def generate_flag_tables():
    """Write one flag table per example into docs/_generated/.

    The example pages include these rather than restating the flags, so a page
    cannot claim a flag the tool no longer produces. Built from the reports where
    they exist and from each `meta.toml`'s `[expect]` table otherwise, which
    `examples/build.py --check` has already proved agree.
    """
    GENERATED.mkdir(parents=True, exist_ok=True)
    for meta in sorted((ROOT / 'examples').glob('*/meta.toml')):
        name = meta.parent.name
        expect = tomllib.loads(meta.read_text()).get('expect', {})
        blocks = []
        for report in sorted(expect):
            csv_path = REPORTS / name / report / 'gb-qc-report.csv'
            if csv_path.is_file():
                rows = _flag_rows(csv_path)
                # A split report measures leaked percentages, not AU-ROC.
                leakage = any(check == 'Data Leakage' for check, _, _ in rows)
                measure = 'Leakage' if leakage else 'AU-ROC'
            else:
                # No report built (a --skip-reports preview). The declared flags
                # are still the truth; only the figures are unavailable.
                rows = [(check, flag, '—') for check, flag in expect[report].items()]
                measure = 'AU-ROC'
            # Relative to the page that includes this snippet - an example page
            # at docs/examples/<name>.md - not to docs/_generated/ where it is
            # written. mkdocs resolves the link from the including source file.
            href = f"../reports/{name}/{report}/gb-qc-report.html"
            blocks.append(
                f"**`{report}`** — [open the report]({href})\n\n"
                f"| Check | Flag | {measure} |\n|---|---|---|\n"
                + '\n'.join(f"| {check} | {FLAG_SPANS.get(flag, flag)} | {figure} |"
                            for check, flag, figure in rows))
        (GENERATED / f'{name}-flags.md').write_text('\n\n'.join(blocks) + '\n')
    print(f"wrote {len(list(GENERATED.glob('*-flags.md')))} flag table(s) "
          f"into {GENERATED.relative_to(ROOT)}")


POSITION_CHECKS = (
    ('Per position nucleotide content', 'counted from the start of each sequence'),
    ('Per position reversed nucleotide content', 'counted from the end of each sequence'),
)
POSITION_ROW = re.compile(r' - (?P<base>[A-Z]) position (?P<pos>\d+)$')
POSITION_LIMIT = 20


def _flagged_positions(report_csv: Path, check: str) -> list[tuple[int, str, float, str]]:
    """(position, base, auroc, flag) for the positions one check flagged.

    A position is scored once per base, so a single position can appear four
    times. The headline check reports the worst of them, and so does this: one
    row per position, carrying the base that made it the worst.
    """
    worst: dict[int, tuple[str, float, str]] = {}
    with report_csv.open(newline='') as handle:
        for row in csv.DictReader(handle):
            if not row['Check'].startswith(check + ' - '):
                continue
            match = POSITION_ROW.search(row['Check'])
            if not match or row['Flag'] not in ('Warning', 'Fail'):
                continue
            position, auroc = int(match['pos']), float(row['AU-ROC'])
            if position not in worst or auroc > worst[position][1]:
                worst[position] = (match['base'], auroc, row['Flag'])
    return [(pos, *rest) for pos, rest in sorted(worst.items())]


def _position_table(positions: list[tuple[int, str, float, str]]) -> str:
    """A Markdown table of flagged positions, worst first if it has to be cut."""
    shown = sorted(positions, key=lambda row: -row[2])[:POSITION_LIMIT]
    table = ('| Position | Worst base | AU-ROC | Flag |\n|---|---|---|---|\n'
             + '\n'.join(f"| {pos} | {base} | {auroc:.3f} | "
                         f"{FLAG_SPANS.get(flag, flag)} |"
                         for pos, base, auroc, flag in sorted(shown)))
    omitted = len(positions) - len(shown)
    if omitted:
        # Say so, rather than letting a cut table read as the whole finding.
        table += (f"\n\nThe {POSITION_LIMIT} highest-scoring of "
                  f"{len(positions)} flagged positions; {omitted} more are in "
                  f"the report.")
    return table


# Values that legitimately appear in prose without coming from a report: the
# flag boundaries, the --min-coverage default, and the ends of the AU-ROC range.
# They are constants of the tool, not measurements of an example.
TOOL_CONSTANTS = frozenset({0.0, 0.25, 0.5, 0.6, 0.7, 0.9, 1.0})
# Trailing \w excludes 'mid-0.80s' and the like: prose that is deliberately
# vague is not quoting a measurement.
PROSE_FIGURE = re.compile(r'(?<![\d.])0\.\d{2,4}(?![\w%])')


def check_prose_figures():
    """Fail if an example page quotes a figure its own reports do not contain.

    The flag tables and position tables are generated, but the prose around them
    is written by hand, and a hand-typed AU-ROC is exactly the kind of thing that
    survives a change to the data it described. This catches that: every
    report-shaped decimal on an example page has to be within rounding distance
    of some value in one of that example's own reports.

    Percentages are skipped - they are matched by the trailing `%` and read off
    the leakage rows, which are not in the same columns.
    """
    wrong = []
    for page in sorted((DOCS / 'examples').glob('*.md')):
        if page.stem == 'index':
            continue
        measured = set()
        for report in (REPORTS / page.stem).rglob('gb-qc-report.csv'):
            with report.open(newline='') as handle:
                for row in csv.DictReader(handle):
                    measured.update(float(row[col]) for col in
                                    ('AU-ROC', 'AU-PR', 'Accuracy') if row.get(col))
        if not measured:
            continue  # a --skip-reports preview has nothing to check against
        for number, line in ((m, n) for n, text in
                             enumerate(page.read_text().splitlines(), 1)
                             for m in PROSE_FIGURE.finditer(text)):
            quoted = float(number.group())
            if quoted in TOOL_CONSTANTS:
                continue
            # Rounding distance at the precision the page chose to write, so
            # "0.783" matches a report's 0.7835 but not its 0.7935.
            places = len(number.group().split('.')[1])
            if not any(abs(value - quoted) <= 0.5 * 10 ** -places + 1e-12
                       for value in measured):
                wrong.append(f"{page.relative_to(ROOT)}:{line}: {number.group()} "
                             f"is in no {page.stem} report")
    if wrong:
        raise SystemExit('figures on example pages that no report produced:\n  '
                         + '\n  '.join(wrong))
    print("every figure quoted on an example page is in that example's reports")


def generate_check_coverage():
    """Write the check-to-example cross-reference into docs/_generated/.

    A check nobody has seen fail is hard to reason about, so the examples index
    lists where each one does. Hand-maintaining that table went wrong twice -
    two examples missing from rows they belong in, and the reversed per-position
    check missing altogether - so it is read off the reports instead.
    """
    GENERATED.mkdir(parents=True, exist_ok=True)
    order: list[str] = []
    where: dict[str, dict[str, set[str]]] = {}
    for meta in sorted((ROOT / 'examples').glob('*/meta.toml')):
        name = meta.parent.name
        for report in sorted(tomllib.loads(meta.read_text()).get('expect', {})):
            csv_path = REPORTS / name / report / 'gb-qc-report.csv'
            if not csv_path.is_file():
                continue
            with csv_path.open(newline='') as handle:
                for row in csv.DictReader(handle):
                    check = row['Check']
                    if ' - ' in check:
                        continue
                    if check not in where:
                        order.append(check)
                        where[check] = {'Fail': set(), 'Warning': set()}
                    if row['Flag'] in where[check]:
                        where[check][row['Flag']].add(name)

    def cell(names: set[str]) -> str:
        # Links are resolved from the including page, docs/examples/index.md.
        return ', '.join(f"[{n}]({n}.md)" for n in sorted(names)) or '—'

    if order:
        lines = ['| Check | Fails in | Warns in |', '|---|---|---|']
        lines += [f"| {check} | {cell(where[check]['Fail'])} | "
                  f"{cell(where[check]['Warning'])} |" for check in order]
    else:
        # A --skip-reports preview has nothing to read this off.
        lines = ['This table is built from the example reports, which this '
                 'build skipped.']
    (GENERATED / 'check-coverage.md').write_text('\n'.join(lines) + '\n')
    print(f"wrote the coverage of {len(order)} check(s) "
          f"into {(GENERATED / 'check-coverage.md').relative_to(ROOT)}")


def generate_position_tables():
    """Write one flagged-position table per example into docs/_generated/.

    Hand-typing these is how a page ends up claiming an AU-ROC the report never
    produced, which is what had happened to the hidden-motif page before this
    existed. Both directions are tabulated: on fixed-length data the reversed
    table is the forward one mirrored, but on variable-length data it is a
    different finding - see the enhancers example, where the forward check
    flags the first base and the reversed check the last.
    """
    GENERATED.mkdir(parents=True, exist_ok=True)
    for meta in sorted((ROOT / 'examples').glob('*/meta.toml')):
        name = meta.parent.name
        blocks = []
        for report in sorted(tomllib.loads(meta.read_text()).get('expect', {})):
            csv_path = REPORTS / name / report / 'gb-qc-report.csv'
            if not csv_path.is_file():
                continue
            for check, direction in POSITION_CHECKS:
                positions = _flagged_positions(csv_path, check)
                if positions:
                    blocks.append(f"**`{report}` — {check}**, {direction}:\n\n"
                                  + _position_table(positions))
        if not blocks:
            # Either nothing was flagged, or this is a --skip-reports preview
            # with no report to read. Write the file either way: a page that
            # includes it must still build, and check_paths would fail first.
            blocks = ['No per-position check is flagged in this example, '
                      'or the reports were not built.']
        (GENERATED / f'{name}-positions.md').write_text('\n\n'.join(blocks) + '\n')
    print(f"wrote {len(list(GENERATED.glob('*-positions.md')))} flagged-position "
          f"table(s) into {GENERATED.relative_to(ROOT)}")


PLACEHOLDER = """<!doctype html>
<title>Report not built</title>
<body style="font: 15px/1.6 system-ui; max-width: 34em; margin: 4em auto; padding: 0 1em">
<h1>Report not built</h1>
<p>This is a placeholder from a <code>--skip-reports</code> preview build. The
real report for <strong>{name}</strong> is generated by
<code>examples/build.py</code>; run <code>scripts/build_docs.py</code> without
<code>--skip-reports</code> to see it.</p>
"""


def expected_report_paths() -> list[tuple[str, str]]:
    """Every (example, report path) pair the examples declare they produce.

    Read from the `[expect]` tables in each `meta.toml`, which is the same
    declaration `examples/build.py --check` asserts against - so the two cannot
    disagree about which reports are supposed to exist.
    """
    pairs = []
    for meta in sorted((ROOT / 'examples').glob('*/meta.toml')):
        expect = tomllib.loads(meta.read_text()).get('expect', {})
        pairs += [(meta.parent.name, path) for path in sorted(expect)]
    return pairs


def generate_screenshot(require: bool):
    """Shoot the landing page's screenshot from a report this build produced.

    scripts/screenshot_report.py explains what it frames and why. `require` is
    dropped only for prose builds, where accepting a placeholder beats
    demanding a 150 MB browser download; any build that could be published
    insists on the real image, so a broken install fails the build instead of
    quietly publishing a placeholder.
    """
    argv = [sys.executable, str(ROOT / 'scripts' / 'screenshot_report.py')]
    if require:
        argv.append('--require')
    _run(argv)


def placeholder_reports():
    """Stand in for the reports so links to them still resolve.

    Without this, `--skip-reports` fails the build: mkdocs checks that a link
    points at a file it knows about, and every report link would dangle. Writing
    stubs at the declared paths keeps that check meaningful during prose work -
    a link to a report no example produces still fails, which is the case worth
    catching.
    """
    for name, report in expected_report_paths():
        target = REPORTS / name / report / 'gb-qc-report.html'
        if target.exists():
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(PLACEHOLDER.format(name=f"{name}/{report}"))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--serve', action='store_true',
                        help="Serve with live reload instead of building once")
    parser.add_argument('--skip-reports', action='store_true',
                        help="Leave docs/reports/ alone; for working on prose")
    args = parser.parse_args()

    check_installed_matches_tree()
    generate_cli_page()
    if args.skip_reports:
        REPORTS.mkdir(parents=True, exist_ok=True)
        placeholder_reports()
        print(f"skipped example reports; {len(expected_report_paths())} "
              f"placeholder(s) in place so report links still resolve")
    else:
        generate_reports()
    generate_screenshot(require=not args.skip_reports)
    generate_flag_tables()
    generate_position_tables()
    generate_check_coverage()
    check_prose_figures()

    env = {**os.environ, 'MKDOCS_CONFIG_FILE': 'mkdocs.yml'}
    _run(['mkdocs', 'serve' if args.serve else 'build'], env=env)
    return 0


if __name__ == '__main__':
    sys.exit(main())
