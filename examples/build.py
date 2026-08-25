#!/usr/bin/env python3
"""Run `gb-qc` over the examples and check the flags came out as documented.

Reports are not committed - they are built here, from the committed data, by the
same commands the docs quote. Two uses:

    python examples/build.py --out-folder build/examples
    python examples/build.py --out-folder build/examples --check

The first regenerates every report. The second does that and then compares each
report's headline flags against the `[expect]` table in the example's
`meta.toml`, exiting non-zero on any difference. That is what stops a docs page
from claiming "six checks fail" after the tool's behaviour has moved on.

Only the headline checks are compared, not the per-base and per-position
sub-checks. A page cites the nine top-level flags; asserting on the several
thousand others would fail on noise and teach everyone to ignore it.

To adopt current output as the new expectation - after a deliberate change, or
when adding an example - use --print-expect and paste the block it emits into
`meta.toml`.
"""

import argparse
import csv
import subprocess
import sys
import tomllib
from pathlib import Path

EXAMPLES = Path(__file__).resolve().parent
REPORT = 'gb-qc-report.csv'


def _load(name: str | None = None):
    """Every example's (directory, parsed meta.toml), or just the one named."""
    metas = sorted(EXAMPLES.glob('*/meta.toml'))
    if not metas:
        raise SystemExit(f"no examples found under {EXAMPLES}")
    found = [(m.parent, tomllib.loads(m.read_text())) for m in metas
             if name is None or m.parent.name == name]
    if not found:
        known = ', '.join(sorted(m.parent.name for m in metas))
        raise SystemExit(f"no such example: {name}. Known: {known}")
    return found


def run_example(directory: Path, meta: dict, out_root: Path) -> int:
    """Run every command an example declares. Returns the failure count.

    Commands run with `cwd` set to the example directory, so the relative data
    paths in `meta.toml` are the same paths the docs show a reader typing.
    """
    out = (out_root / directory.name).resolve()
    failures = 0
    for spec in meta.get('run', []):
        argv = ['gb-qc', spec['command'], *spec['args'],
                '--out-folder', str(out), '--log-level', 'WARNING']
        print(f"  {spec['id']:6s} {' '.join(argv[1:3])} ...", flush=True)
        result = subprocess.run(argv, cwd=directory, capture_output=True, text=True)
        if result.returncode != 0:
            failures += 1
            print(f"  {spec['id']:6s} FAILED (exit {result.returncode})")
            for line in (result.stderr or result.stdout).splitlines()[-12:]:
                print(f"         {line}")
    return failures


def headline_flags(report: Path) -> dict[str, str]:
    """The report's top-level check flags, keyed by check name.

    Rows whose check name contains ' - ' are a per-base or per-position
    breakdown of the check above them, and are deliberately skipped.
    """
    with report.open(newline='') as handle:
        return {row['Check']: row['Flag'] for row in csv.DictReader(handle)
                if ' - ' not in row['Check']}


def collect(directory: Path, out_root: Path) -> dict[str, dict[str, str]]:
    """Every report an example produced, keyed by its path below the example."""
    out = out_root / directory.name
    return {str(path.parent.relative_to(out)): headline_flags(path)
            for path in sorted(out.rglob(REPORT))}


def check(directory: Path, meta: dict, out_root: Path) -> list[str]:
    """Differences between what an example produced and what it expects."""
    expected = meta.get('expect')
    if not expected:
        return [f"{directory.name}: no [expect] table in meta.toml"]

    actual = collect(directory, out_root)
    problems = []
    for report_path, checks in sorted(expected.items()):
        if report_path not in actual:
            problems.append(f"{directory.name}: no report at {report_path}")
            continue
        got = actual[report_path]
        for name, flag in sorted(checks.items()):
            if name not in got:
                problems.append(f"{directory.name}/{report_path}: no check {name!r}")
            elif got[name] != flag:
                problems.append(
                    f"{directory.name}/{report_path}: {name!r} is {got[name]}, "
                    f"expected {flag}")
    for report_path in sorted(set(actual) - set(expected)):
        problems.append(f"{directory.name}: report at {report_path} is not in [expect]")
    return problems


def print_expect(directory: Path, out_root: Path):
    """Emit an [expect] block matching what the example just produced."""
    print(f"\n# ---- {directory.name}/meta.toml ----")
    for report_path, checks in sorted(collect(directory, out_root).items()):
        print(f'\n[expect."{report_path}"]')
        for name, flag in checks.items():
            print(f'"{name}" = "{flag}"')


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out-folder', required=True, type=Path,
                        help="Where reports are written; one directory per example")
    parser.add_argument('--example', default=None,
                        help="Only this example. Default: all of them")
    parser.add_argument('--check', action='store_true',
                        help="Compare the flags against [expect] and fail on a difference")
    parser.add_argument('--print-expect', action='store_true',
                        help="Print [expect] blocks for the reports just built")
    parser.add_argument('--skip-build', action='store_true',
                        help="Only check or print, reusing reports already in --out-folder")
    args = parser.parse_args()

    examples = _load(args.example)
    failures = 0
    for directory, meta in examples:
        print(directory.name, flush=True)
        if not args.skip_build:
            failures += run_example(directory, meta, args.out_folder)

    if failures:
        print(f"\n{failures} command(s) failed", file=sys.stderr)
        return 1

    if args.print_expect:
        for directory, _ in examples:
            print_expect(directory, args.out_folder)

    if args.check:
        problems = [p for directory, meta in examples
                    for p in check(directory, meta, args.out_folder)]
        if problems:
            print(f"\n{len(problems)} flag mismatch(es):", file=sys.stderr)
            for problem in problems:
                print(f"  {problem}", file=sys.stderr)
            return 1
        print("\nall examples match their expected flags")

    return 0


if __name__ == '__main__':
    sys.exit(main())
