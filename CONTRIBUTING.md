# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

You can contribute in many ways:

## Types of Contributions

### Report Bugs

Report bugs at <https://github.com/genomic-benchmarks/genomic-benchmarks-qc/issues>.

If you are reporting a bug, please include:

-   Your operating system name and version.
-   Any details about your local setup that might be helpful in troubleshooting.
-   Detailed steps to reproduce the bug.

### Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with “bug” and “help wanted” is open to whoever wants to implement it.

### Implement Features

Look through the GitHub issues for features. Anything tagged with “enhancement” and “help wanted” is open to whoever wants to implement it.

### Submit Feedback

The best way to send feedback is to file an issue at <https://github.com/genomic-benchmarks/genomic-benchmarks-qc/issues>.

If you are proposing a feature:

-   Explain in detail how it would work.
-   Keep the scope as narrow as possible, to make it easier to implement.
-   Remember that this is a volunteer-driven project, and that contributions are welcome :)

## Workflow

- make an issue for the thing you want to implement
- create the corresponding branch
- develop
- write unit tests in tests/
- document it: a docstring on anything public, and a page under docs/ if it
  changes what a user does
- push the changes
- make a pull request
- once the pull request is merged, the issue will be closed

## Get Started!

Ready to contribute? Here’s how to set up Genomic Benchmarks QC for local development.

1.  Fork the Genomic Benchmarks QC repo on GitHub.
2.  Clone your fork locally:

    ```bash
    git clone git@github.com:your_name_here/genomic-benchmarks-qc.git
    ```

3.  Install your local copy into a conda environment. Assuming you have conda installed, this is how you set up your fork for local development.


    ```bash
    cd genomic-benchmarks-qc/
    conda env create -f dev-requirements.yml
    conda activate gb-qc-dev
    pip install -e '.[develop]'
    ```

4.  Create a branch for local development:

    ```bash
    git checkout -b name-of-your-bugfix-or-feature
    ```

    Now you can make your changes locally.

5.  When you’re done making changes, check that your changes pass the tests:

    ```bash
    pytest
    ```

    Test paths and the `src/` import path are configured in `pyproject.toml`, so plain `pytest` works from the repository root. To run a single file or test:

    ```bash
    pytest tests/test_cli_helpers.py
    pytest -k SplitMemoryLimit
    ```

    To check unit test coverage, add `--cov`:

    ```bash
    pytest --cov
    ```

    The coverage settings live under `[tool.coverage.*]` in `pyproject.toml`. Branch coverage is
    on, and the run fails if total coverage drops below 60%, which is also enforced in CI.
    Total coverage is well above that floor - 88% at the time of writing, and CI's test job
    prints the live figure as a per-file table in its job summary. Uncovered line numbers are
    printed per file; for a browsable line-by-line report use:

    ```bash
    pytest --cov --cov-report=html
    open htmlcov/index.html
    ```

    Coverage is opt-in so that running a single test file does not trip the threshold. New code
    should come with tests, and the threshold should be raised as coverage improves.

6.  Check that your changes pass linting:

    ```bash
    ruff check .
    ```

    The rules live under `[tool.ruff]` in `pyproject.toml` and are the same ones CI runs, so a
    clean run locally means a clean run on the pull request. Many findings fix themselves:

    ```bash
    ruff check --fix .
    ```

    Two exemptions are worth knowing about. The modules that hold the report's HTML as Python
    strings are exempt from the line-length limit, because a `# noqa` cannot be placed inside a
    triple-quoted string. Tests are exempt from it too. Everything else is expected to stay within
    100 columns.

7.  If you changed the report's styling or behaviour, the stylesheets and scripts live in
    `src/genomic_benchmarks_qc/report/assets/` and are inlined into the generated HTML. They are
    real files so they can be checked without generating a report:

    ```bash
    (for f in src/genomic_benchmarks_qc/report/assets/*.js; do node --check "$f" || exit 1; done)
    ```

    CI runs the same loop in the lint job, so a stray bracket fails the pull request rather than
    the reader's browser.

    Every asset has to be declared as package data in `pyproject.toml`, otherwise it works from a
    checkout and is missing from the installed package. `tests/test_report_assets.py` fails when a
    new one is not.

8.  If you changed a docstring, a page under `docs/`, or anything under `examples/`, build the
    documentation site. mkdocs is not in the dev environment - the docs dependencies are pinned
    separately in `docs/requirements.txt` - so install them first:

    ```bash
    pip install -r docs/requirements.txt
    python scripts/build_docs.py --skip-reports
    ```

    `--skip-reports` leaves `docs/reports/` alone instead of running `gb-qc` over all eight
    examples, which takes about two and a half minutes; anything missing there gets a
    placeholder so the report links still resolve. Drop it when you changed anything a report's
    numbers depend on - that full run is what asserts each example still produces the flags its
    `meta.toml` claims.

    The site is built with `--strict`, so a cross-reference that stops resolving fails the build
    rather than rendering as bracketed text. Two things to know before writing one:

    - A docstring cross-reference is `[NAME][genomic_benchmarks_qc.path.NAME]`, and
      `tests/test_docstring_references.py` checks that the target exists and that a documented
      `Default:` still matches the constant it cites. Plain `pytest` catches a dead target; only
      the build catches a link that fails to render.
    - The pages under `docs/reference/api/` name their members explicitly, so a new public
      function or constant is not published until it is added to that page's `members:` list.

    The Docs workflow runs the same build on every pull request touching `src/`, `docs/`,
    `examples/`, `mkdocs.yml` or `pyproject.toml`.

9.  If you added a file under `tests/` that is not a `test_*.py` module - a fixture, a helper, a
    data file - add it to `MANIFEST.in`. setuptools finds the test modules by itself and nothing
    else, so without it the source distribution ships a suite that cannot run. The same trap as
    package data above, one directory further out:

    ```bash
    pip install build
    rm -rf dist                          # `build` adds to dist/, it does not replace it
    python -m build
    tar tzf dist/*.tar.gz | grep /tests/
    ```

    The same applies to anything outside `tests/` that a test reads: the sdist carries
    `docs/guide/` only because `tests/test_report_links.py` opens those pages, and that entry is
    in `MANIFEST.in` for the same reason.

    CI's package job builds both artifacts, runs `twine check` on them, and then runs the whole
    suite from the unpacked tarball. That is the only place a file missing from the sdist shows
    up: every other check here reads your working tree, where the file is still sitting.

10. Commit your changes and push your branch to GitHub:

    ```bash
    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature
    ```

11. Submit a pull request through the GitHub website.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1.  The pull request should include tests.
2.  `ruff check .` should report no errors; CI runs the same check and fails the pull request otherwise.
3.  If the pull request adds functionality, the docs should be updated: a docstring on anything public, and its name added to the relevant page under `docs/reference/api/`.
4.  The pull request should work for Python >=3.12. CI runs the suite on 3.12, 3.13 and 3.14; the dev environment pins 3.14, so the older two are only exercised on the pull request.
