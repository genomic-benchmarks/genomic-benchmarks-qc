"""What importing the CLI is allowed to cost.

`report_generator` defers seaborn and matplotlib into the two functions that
draw, and `testing` defers sklearn into the one function that scores a
continuous feature. Both are load-bearing rather than tidy: between them they
are most of a second and half the memory of starting the tool up, paid by
`gb-qc --help` and by every run that asks for `simple` or `json` reports and
never draws a figure.

A deferred import is one `from x import y` away from coming back at the top of
some module, and nothing about a passing test suite would notice. So these run
the import in a subprocess and look at what arrived with it.
"""

import subprocess
import sys

# The stacks that must stay out of the way until something actually needs them,
# and roughly what each cost when it did not.
DEFERRED = ['matplotlib', 'seaborn', 'scipy', 'sklearn']


def modules_after(statement):
    """Names in `sys.modules` after running `statement` in a fresh interpreter."""
    script = f"import sys\n{statement}\nprint('\\n'.join(sorted(sys.modules)))\n"
    result = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return set(result.stdout.split())


class TestNothingHeavyLoadsWithTheCli:

    def test_importing_the_cli_leaves_the_heavy_stacks_alone(self):
        loaded = modules_after('import genomic_benchmarks_qc.cli')
        arrived = [name for name in DEFERRED if name in loaded]
        assert arrived == [], (
            f"{', '.join(arrived)} loaded just by importing the CLI. Something took an "
            f"import back to the top of a module - see the docstrings in "
            f"report_generator and utils.testing for why they are where they are."
        )

    def test_help_does_not_load_them_either(self):
        """The command Typer answers without running anything."""
        loaded = modules_after(
            "import genomic_benchmarks_qc.cli as c\n"
            "sys.argv = ['gb-qc', '--help']\n"
            "try:\n"
            "    c.app()\n"
            "except SystemExit:\n"
            "    pass"
        )
        assert [name for name in DEFERRED if name in loaded] == []

    def test_the_plotting_stack_arrives_when_a_figure_is_drawn(self):
        """The other half of the claim: deferred, not dropped."""
        loaded = modules_after(
            "from genomic_benchmarks_qc.report import report_generator\n"
            "import inspect\n"
            "assert 'matplotlib' not in sys.modules\n"
            "from genomic_benchmarks_qc.report import classes_plots"
        )
        assert 'matplotlib' in loaded
        assert 'seaborn' in loaded

    def test_a_simple_report_run_never_loads_the_plotting_stack(self):
        """The claim in the round: a run that draws nothing pays for nothing.

        sklearn and scipy do arrive - the per-sequence checks score through them
        - which is why this names the two that should not rather than reusing
        `DEFERRED`.
        """
        loaded = modules_after(
            "import tempfile, pathlib\n"
            "from genomic_benchmarks_qc.evaluate_classes import run\n"
            "folder = pathlib.Path(tempfile.mkdtemp())\n"
            "rows = ['sequence,label']\n"
            "rows += [f'ACGTACGTAC,{i % 2}' for i in range(20)]\n"
            "data = folder / 'in.csv'\n"
            "data.write_text('\\n'.join(rows))\n"
            "run(input=[str(data)], format='csv', out_folder=str(folder / 'out'),\n"
            "    report_types=['simple', 'json'], log_level='ERROR')\n"
            "assert list((folder / 'out').rglob('*.csv'))"
        )
        assert 'matplotlib' not in loaded
        assert 'seaborn' not in loaded
