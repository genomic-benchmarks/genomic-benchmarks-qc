# Example datasets

Eight datasets, each chosen because it is the only one that shows a particular
thing. Together they exercise every check `gb-qc` runs, at every severity that
check can reach.

The data here is committed, as plain CSV, TSV and FASTA — no compression, so you
can open any of it and look. The **reports are not** committed: they are built
from this data by [`build.py`](build.py), using the same commands the
documentation shows, so a published report is always the one the current code
produces.

```bash
# every example, into build/examples/<name>/
python examples/build.py --out-folder build/examples

# and check the flags still match what the docs claim
python examples/build.py --out-folder build/examples --check
```

A full run takes about two and a half minutes and produces 15 reports.

## What each one is for

| Example | Shows | Class flags | Leakage |
|---|---|---|---|
| [clean-dataset](clean-dataset/) | The control: what "nothing wrong" looks like. AU-ROC 0.50–0.54 throughout | 9 Pass | 0.07% Warning |
| [composition-bias](composition-bias/) | The worst case, and the only example that fails both kinds of report | 6 Fail, 1 Warning, 2 Pass | 6.04% **Fail** |
| [hidden-motif](hidden-motif/) | A bias 8 nt wide inside 398 positions. Why the per-position plot is interactive | 2 Fail, 1 Warning, 6 Pass | 0.00% Pass |
| [variable-length](variable-length/) | The `Unknown` tail and `--min-coverage`; non-ACGT bases. 20,870 positions unscored against 5,536 scored | 2 Fail, 1 Warning, 6 Pass | 0.80% Warning |
| [length-bias](length-bias/) | Length alone separating the classes, plus `--regression` and `--end-position` | 3 Warning⁻, 6 Pass | — |
| [paired-sequences](paired-sequences/) | Two sequence columns in one row, and the extra `merged` report | duplication checks Fail, statistics Pass | — |
| [fasta-classes](fasta-classes/) | One FASTA file per class, label from the filename | 3 Fail, 6 Pass | — |
| [enhancers](enhancers/) | The README's quickstart dataset. GC content alone separates the classes at AU-ROC 0.66 | 6 Warning, 3 Pass | 0.67% Warning |

⁻ `length-bias` reaches Warning on the length check, not Fail, and that is the
strongest available: across all 234 dataset splits surveyed for the paper, **no
dataset fails the length check outright**. The strongest length bias on record
among them is AU-ROC 0.61; the sample here scores 0.63. So the README's first
claim — that a length classifier can beat your model — is demonstrated at
Warning severity, because in this corpus that is as far as it goes.

## Layout

```text
<example>/
├── meta.toml     # provenance, how the data was derived, commands, expected flags
└── data/         # committed inputs
```

`meta.toml` is the single source of truth for an example. It records which
collection and dataset the data came from and how to cite it, how the committed
files were derived from it, the exact commands to run, and the flags those
commands are expected to produce. `build.py` reads it, and so should the docs
pages.

The derivation is documented, not scripted. The committed files are the
artefact: each `meta.toml` says which upstream dataset and split it came from
and how much was taken, which is what you need to go back to the source. It is
not enough to reproduce the exact draw, and deliberately so — regenerating the
data is not a thing anyone should be doing casually, because two examples have
published numbers measured on exactly these bytes (see below).

## Sources and attribution

Every dataset here is a redistribution of work published by someone else, under
terms that permit it. Please cite the original when you use one.

| Examples | Source | Cite |
|---|---|---|
| `clean-dataset`, `composition-bias` | GUE | Zhou et al. *DNABERT-2: Efficient Foundation Model and Benchmark For Multi-Species Genome.* [arXiv:2306.15006](https://arxiv.org/abs/2306.15006) |
| `hidden-motif`, `variable-length` | OmniGenBench | Yang et al. *OmniGenBench: A Modular Platform for Reproducible Genomic Foundation Models Benchmarking.* [arXiv:2505.14402](https://doi.org/10.48550/arXiv.2505.14402) |
| `enhancers`, `fasta-classes` | Genomic Benchmarks | Gresova et al. *Genomic benchmarks: a collection of datasets for genomic sequence classification.* [BMC Genomic Data 24, 25 (2023)](https://doi.org/10.1186/s12863-023-01123-8) |
| `length-bias` | mRNABench | [mRNABench](https://github.com/morrislab/mRNAbench), translation efficiency task |
| `paired-sequences` | miRBench | Hejret et al., AGO2 CLASH data via [miRBench](https://github.com/katarinagresova/miRBench) |
