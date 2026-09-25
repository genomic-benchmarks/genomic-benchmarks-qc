# Contributing

Contributions are welcome: bug reports, ideas, fixes and new checks.

## Reporting a bug or suggesting a feature

Open an issue at <https://github.com/genomic-benchmarks/genomic-benchmarks-qc/issues>.

- For a bug, include the output of `gb-qc --version`, your operating system, the command you
  ran and, if you can, a small input that reproduces it.
- For a feature, describe how it would work, and keep the scope as narrow as you can.

Issues tagged "help wanted" are open to anyone who wants to pick them up.

## Setting up

Fork the repository on GitHub, then:

```bash
git clone git@github.com:your_name_here/genomic-benchmarks-qc.git
cd genomic-benchmarks-qc
conda env create -f dev-requirements.yml
conda activate gb-qc-dev
pip install -e '.[develop]'
```

The tests stub out MMseqs2, so you only need it to run `gb-qc evaluate-splits` yourself or to
build the example reports. The README's Installation section says how to install it.

## Making a change

1. For anything bigger than a small fix, open an issue first so the approach can be agreed on.
2. Create a branch from `main`.
3. Make the change, with tests in `tests/`. Anything public needs a docstring, and a change to
   what a user does needs the page under `docs/` that describes it updated.
4. Run the checks below.
5. Push the branch and open a pull request against `main`.

## Checks

These are the commands CI runs, with the same settings from `pyproject.toml`, so a clean run
here means a clean run on the pull request:

```bash
pytest                    # the test suite
pytest --cov              # ...and fail if total coverage drops below 60%
ruff check .              # lint; `ruff check --fix .` fixes many findings for you
(for f in src/genomic_benchmarks_qc/report/assets/*.js; do node --check "$f" || exit 1; done)
```

The last line checks the syntax of the report's JavaScript, which nothing else parses.

A few `pytest` options worth knowing:

```bash
pytest tests/test_cli_helpers.py          # one file
pytest -k SplitMemoryLimit                # tests matching a name
pytest --cov --cov-report=html            # line-by-line coverage in htmlcov/index.html
```

If you changed a docstring, a page under `docs/` or anything under `examples/`, also build the
documentation site. Its dependencies are not in the dev environment:

```bash
pip install -r docs/requirements.txt
python scripts/build_docs.py --skip-reports
```

`--skip-reports` skips rebuilding the eight example reports, which takes about two and a half
minutes and needs MMseqs2 on `PATH`. Leave it off if you changed anything a report's numbers
depend on: the full build is what checks that each example still produces the flags its
`meta.toml` expects.

### What runs on a pull request

| Workflow and job | What it runs | When |
| --- | --- | --- |
| CI / Lint | `ruff check .` and the `node --check` loop | every pull request |
| CI / Python 3.12, 3.13, 3.14 | `pytest --cov`, then four `gb-qc` runs on the example data | every pull request |
| CI / Package | builds the sdist and wheel, runs `twine check`, and runs the suite from the unpacked sdist | every pull request |
| Docs / Build | builds every example report, compares its flags with `meta.toml`, and builds the site with `--strict` | pull requests that touch `src/`, `docs/`, `examples/`, `mkdocs.yml`, `pyproject.toml` or `scripts/build_docs.py` |

The dev environment uses Python 3.14, so 3.12 and 3.13 are only tested on the pull request.

### Things that work in your checkout but not for users

- **A new file in `report/assets/`** ships in the package only if its extension is in the
  package data in `pyproject.toml` (`*.html`, `*.css`, `*.js`). `tests/test_report_assets.py`
  fails when one is left out.
- **A new file in `tests/` that is not a `test_*.py` module** (a data file, say) is left out of
  the source distribution unless `MANIFEST.in` includes it. Python files there are already
  covered. The same goes for a file outside `tests/` that a test reads, which is why
  `MANIFEST.in` lists `docs/guide/`. Only the Package job notices a missing one.
- **A new public function or constant** appears in the API docs only once you add it to the
  `members:` list of its page under `docs/reference/api/`. Nothing checks this.
- **Docstring cross-references** are written `[NAME][genomic_benchmarks_qc.path.NAME]`.
  `tests/test_docstring_references.py` checks that the target exists and that a documented
  `Default:` matches the constant it cites. Only the docs build catches a reference that
  fails to render.
- **Lines** can be up to 100 characters. Tests are exempt, and so are the modules that hold
  the report's HTML as Python strings, listed under `[tool.ruff.lint.per-file-ignores]`.

## Adding a check

Each check is one module, listed in a registry. The rest of the tool reads the registry to
build the results table, `gb-qc-report.csv`, the shading in the figures, and the report's
navigation, sections and verdict. So a new check is a new module plus one line in the
registry, and two new checks can be written independently.

There is one registry per command:

- `CLASS_CHECKS` in `checks/classes/__init__.py`: the checks `evaluate-classes` runs on each
  pair of classes. A class check scores two `SequenceStatistics`.
- `SPLIT_CHECKS` in `checks/splits/__init__.py`: the checks `evaluate-splits` runs. A split
  check scores the threshold statistics of the MMseqs2 search.

A check is a `Check` (from `checks/__init__.py`) with four parts:

- `name`, shown in the report and the CSV;
- `score`, a function returning a `CheckResult`;
- `floor`, which minimum class size it needs before it is scored (`'per_sequence'` or
  `'per_position'`), or `None` for a check decided by a rule;
- `section`, a `Section` (from `report/sections.py`) saying how the report shows it.

Most checks ask whether one feature tells the two classes apart. For those, the feature is
computed in the check's own module, by a function that takes one class's `SequenceStatistics`
and returns a DataFrame with one row per sequence. `stats.sequences` holds the class's
sequences, so nothing in `utils/seq_stats.py` needs to change. A factory in
`checks/scorers.py` then scores and flags the feature the same way as every other check:

```python
NAME = 'Per sequence CpG count'

def cpg_count(stats):
    """One row per sequence, one column: the feature being compared."""
    return pd.DataFrame({NAME: [float(seq.count('CG')) for seq in stats.sequences]})

CHECK = Check(
    name=NAME,
    score=scalar_feature(NAME, cpg_count),      # or column_features, one sub-check per column
    floor='per_sequence',
    section=Section(
        title='Per Sequence CpG Count',
        anchor='per-sequence-cpg-count',
        explanation_id='cpg-explanation',
        template='check_cpg_count.html',
        render=render,
        docs=(('checks', 'per-sequence-cpg-count', 'What to do about it'),),
    ),
)
```

To add one:

1. **Write its module** in `checks/classes/` or `checks/splits/`. Start from the closest
   existing one: `gc_content.py` scores one value per sequence, `nucleotide_content.py` one
   column per base, and `unique_bases.py` uses a rule instead of a score. The feature
   function runs again for every pair of classes a class is compared in, so if it is slow,
   cache its result per class.
2. **Write its `render` function.** It gets the report's context (`ClassReportContext` or
   `SplitReportContext`), draws the figures into `context.plots_dir` with `save_figure`, and
   returns a `SectionContent` that fills the fragment's placeholders.
   - Pass anything that came from the input files through `escape_html_text`.
   - Import the plotting modules inside `render`, not at the top of the module.
     `gb-qc --help` imports the check modules, and `tests/test_startup.py` fails if
     matplotlib comes with them.
3. **Add its fragment**, `report/assets/check_<name>.html`. It starts with the explanation,
   `<div id="{{explanation_id}}" class="explanation-text">…{{docs_link}}</div>`, followed by
   the section's body. The heading, flag and `?` button are shared, so leave them out.
4. **List it in the registry**, at the position it should have in the report.
5. **Give it a `##` heading in `docs/guide/checks.md`.** Its docs link points there, and a
   class check's `anchor` must match that heading's anchor. `tests/test_report_links.py` and
   `tests/test_check_registry.py` check both.
6. **Add its flag to every example's `[expect]` table** in `examples/*/meta.toml`.
   `python examples/build.py --out-folder build/examples --print-expect` prints the tables.
   `--check` ignores checks a table does not name, so it will not notice one missing.
7. **Update the number of checks** ("Nine checks") in the docs and the README.

`tests/test_check_registry.py` holds the rules every check must follow, such as unique names
and section ids. It also passes a dummy check through both reports to show that nothing else
needs editing.

## Cutting a release

This is a maintainer's job. PyPI accepts each version number only once, even after the
release is deleted, so the steps put the upload last and make it wait for an approval.

1. **Bump the version** in `src/genomic_benchmarks_qc/__init__.py`, the only place it lives:

   ```python
   __version__ = "1.0.1"
   ```

   `pyproject.toml`, `gb-qc --version` and every report's footer read it from there. Merge the
   bump into `main` through a pull request.
2. **Publish a GitHub Release** on `main`, tagged `vX.Y.Z`. Publishing the Release starts the
   upload; pushing a tag on its own does not.
3. **Approve the `pypi` deployment.** The publish job waits for a required reviewer.
4. **Check that it installs from PyPI:**

   ```bash
   python -m venv /tmp/gb-qc-check
   /tmp/gb-qc-check/bin/pip install genomic-benchmarks-qc
   /tmp/gb-qc-check/bin/gb-qc --version
   ```

Before anything is uploaded, CI makes sure of three things:

- Lint, all three test jobs and the Package job must pass on the release commit.
- The release tag must match the version of the built files. A forgotten bump fails the run
  instead of taking up a version number.
- The files uploaded are the same ones the Package job built and tested.

The upload uses PyPI Trusted Publishing, so there is no API token. PyPI's side names
`.github/workflows/ci.yml` and the `pypi` environment. Renaming either breaks the upload, and
the error calls the workflow an unauthorised publisher. The reasons behind each step are in
the comments in `ci.yml`.
