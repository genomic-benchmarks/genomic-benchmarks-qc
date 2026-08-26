# Installation

## Installing `gb-qc`

```bash
pip install genomic-benchmarks-qc
```

That is everything `evaluate-classes` needs. Python 3.12 or newer; tested on
3.12, 3.13 and 3.14. The scientific stack it depends on — NumPy, pandas,
matplotlib, seaborn, scikit-learn, SciPy, Biopython and Typer — comes in with it.

Check it landed:

```bash
gb-qc --help
```

You should see the two commands, `evaluate-classes` and `evaluate-splits`.

## MMseqs2, for the leakage check

`evaluate-splits` measures how much of your test set already appears in your
training set, and it does that with a similarity search rather than string
matching. That search is [MMseqs2](https://github.com/soedinglab/MMseqs2), which
is a separate binary: it is not a Python package and `pip` will not bring it in.
It has to be on your `PATH` under the name `mmseqs`.

Use a precompiled binary. MMseqs2 ships static builds for Linux and macOS, and
they are self-contained — pick the one for your machine from the
[installation guide](https://github.com/soedinglab/MMseqs2/wiki#installation),
unpack it, and put its `bin/` directory on your `PATH`.

Then confirm it:

```bash
mmseqs version
```

If that prints a version, `evaluate-splits` will find it. If it prints
`command not found`, see
[the troubleshooting entry](faq.md#mmseqs-command-not-found).

**`evaluate-classes` does not need MMseqs2.** If the class checks are all you
want, skip this section — nothing else in the tool touches it.

## Installing `gb-qc` from source

Back to the Python package — MMseqs2 is done with. To get the current `main`
rather than the last release:

```bash
pip install git+https://github.com/genomic-benchmarks/genomic-benchmarks-qc.git
```

For a checkout you intend to edit, the repository ships a conda environment with
the test and lint tooling in it:

```bash
git clone https://github.com/genomic-benchmarks/genomic-benchmarks-qc.git
cd genomic-benchmarks-qc
conda env create -f dev-requirements.yml
conda activate gb-qc-dev
pip install -e '.[develop]'
```

`CONTRIBUTING.md` covers the rest — running the tests, the lint rules, and how
the example reports are built.

## Next

- [The checks](guide/checks.md) — what the tool measures, and what to do about
  each flag
- [Eight worked examples](examples/index.md) with live reports, if you would
  rather see it working first
- [Troubleshooting](faq.md) — for when something does not behave
