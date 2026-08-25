<img src="https://github.com/genomic-benchmarks/genomic-benchmarks-qc/blob/main/assets/genomic-benchmarks-qc-text.png?raw=True" alt="genomic-benchmarks-qc logo" width="500" />

[![CI](https://github.com/genomic-benchmarks/genomic-benchmarks-qc/actions/workflows/ci.yml/badge.svg)](https://github.com/genomic-benchmarks/genomic-benchmarks-qc/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/genomic-benchmarks-qc)](https://pypi.org/project/genomic-benchmarks-qc/)
[![Python](https://img.shields.io/pypi/pyversions/genomic-benchmarks-qc)](https://pypi.org/project/genomic-benchmarks-qc/)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-genomic--benchmarks.github.io-1f6fd0)](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/)

### Find the shortcut. Learn the biology.

Automated quality control for genomic machine learning datasets: scores the biases, duplicates and data leakage a classifier could exploit before you train on it.

Everything you need is below. The [documentation site](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/) adds what a README has no room for: eight [worked examples](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/examples/) with live reports — a clean dataset, a badly biased one, one whose only flaw is six positions wide — a guide to [what to do about each check](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/guide/checks/) and [why the thresholds sit where they do](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/guide/how-it-works/), and the generated CLI and API reference.

## What it catches

When the classes differ in something trivial, a high score no longer tells you what a model learned — the biology, or the shortcut. `gb-qc` looks for the differences a classifier could exploit without understanding anything:

- **Your negatives are shorter than your positives.** A length classifier now beats your model.
- **Your classes differ in GC content, base composition or dinucleotide frequencies.** In the bundled enhancers example, GC content alone separates the two classes at AU-ROC 0.66.
- **Position 1 is `N` in one class only**, or every sequence in one class starts with the same adapter — a per-position give-away that no summary statistic would show.
- **The same sequence appears in both classes, or your test set repeats your training set.** `evaluate-splits` catches near-duplicates too, with an MMseqs2 similarity search.

Each check gets a **Pass / Warning / Fail** flag, an [HTML report](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/reports/composition-bias/class/sequence/0_vs_1/gb-qc-report.html) you can read, and a [CSV](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/reports/composition-bias/class/sequence/0_vs_1/gb-qc-report.csv) you can put in CI.

<!-- Recorded by `python scripts/screenshot_report.py --mode scroll`, from a
     report the docs build has just built; both README animations are generated
     that way rather than committed, so neither can show an older report than
     the one the site is serving. -->
<a href="https://genomic-benchmarks.github.io/genomic-benchmarks-qc/reports/composition-bias/class/sequence/0_vs_1/gb-qc-report.html"><img src="https://genomic-benchmarks.github.io/genomic-benchmarks-qc/assets/report-scroll.webp" alt="The gb-qc report for the composition-bias example, scrolled from top to bottom: the sidebar's flag summary stays in view with two green ticks, an orange warning and six red failures, while the page moves through the descriptive statistics, the duplicate-sequence tables, the GC-content, nucleotide and dinucleotide plots and the two interactive per-position panels, each headed by the flag it was given" width="1200" /></a>

The report for the bundled `composition-bias` example, scrolled end to end — click it to open the real one. It is the worst case in the gallery, and the one report that carries all three flags at once: six checks **Fail** — the classes share sequences outright, and GC content, nucleotide and dinucleotide composition and both per-position checks all sit above the 0.7 line, the worst at AU-ROC 0.717 — one **Warning**, two **Pass**.

## Installation

```bash
pip install genomic-benchmarks-qc
```

Requires Python 3.12 or newer.

<details>
<summary><b>Extra step for <code>evaluate-splits</code>: mmseqs2</b></summary>

The leakage check uses [mmseqs2](https://github.com/soedinglab/MMseqs2), which has to be installed separately — preferably from a precompiled binary or built from source. See the [mmseqs2 installation guide](https://github.com/soedinglab/MMseqs2/wiki#installation).

`evaluate-classes` does not need it.
</details>

## Quick Start

Point the tool at your dataset and give it somewhere to write:

```bash
gb-qc evaluate-classes \
  --input examples/enhancers/data/enhancers_train.csv \
  --input examples/enhancers/data/enhancers_test.csv \
  --out-folder qc-out
```

It leaves you a report to open — [this one](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/reports/enhancers/class/sequence/0_vs_1/gb-qc-report.html):

```text
qc-out/class/sequence/0_vs_1/
├── gb-qc-report.html   ← open this
├── gb-qc-report.csv
└── plots/
```

To check whether your test set leaks into your training set:

```bash
gb-qc evaluate-splits \
  --train-input examples/enhancers/data/enhancers_train.csv \
  --test-input examples/enhancers/data/enhancers_test.csv \
  --sequence-column sequence \
  --out-folder qc-out
```

with a report of its own — [this one](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/reports/enhancers/split/sequence/enhancers_train_vs_enhancers_test/gb-qc-report.html):

```text
qc-out/split/sequence/enhancers_train_vs_enhancers_test/
├── gb-qc-report.html   ← open this
├── gb-qc-report.csv
└── plots/
```

Both commands can share one `--out-folder`; they write into `class/` and `split/` respectively and never overwrite each other.

## Reading the report

The HTML report is a single standalone file — no external assets, so you can mail it or archive it. It opens with a summary table of every check and its flag, followed by the basic descriptive statistics of each class, then one section per check with the numbers and the plot behind it.

What it checks:

- Nucleotide vocabulary
- Sequence lengths distribution
- GC content per sequence
- Nucleotide composition (per sequence & per position)
- Dinucleotide frequencies
- Sequence duplication levels
- Exact duplicate sequences between classes

Flags come from the AU-ROC of a classifier that sees only that one feature:

| Flag | AU-ROC | Meaning |
|------|--------|---------|
| **Pass** | ≤ 0.6 | Classes not distinguishable by this feature |
| **Warning** | ≤ 0.7 | A model could get some traction here |
| **Fail** | > 0.7 | Significant bias detected |
| **Unknown** | — | Not enough sequences to score the check |

**`Unknown` is not `Pass`**: it says the comparison was not made, not that it came out clean. A check needs at least **250 sequences** before it is scored at all. The plots, per-class statistics and descriptive tables are computed from all the data either way, so small datasets can still be compared by eye, and the terminal says which checks were skipped and why.

### The per-position plot is interactive

Full walkthrough: [the per-position plot](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/guide/per-position/).

For long sequences a static plot cannot show much: a single flagged position is lost in a figure spanning hundreds of them. The per-position panels are therefore interactive, so you can zoom in and explore the individual failing positions.

<!-- Recorded by `python scripts/screenshot_report.py --mode panel`. -->
<img src="https://genomic-benchmarks.github.io/genomic-benchmarks-qc/assets/per-position-demo.webp" alt="The per-position panel being driven: a drag across the flagged region takes 398 positions down to 48, the tooltip at position 200 reads G at 0.217 in one class against 0.784 in the other with a Fail flag beside it, and one press of Next flag returns to the same window without the drag" width="800" />

- **Drag** to zoom, shift-drag to pan, double-click to reset.
- **Hover** for the per-class base frequencies and flags at a position.
- **Prev / Next flag** jumps to the flagged positions.
- **Save view** writes the current window as a PNG.
- The **table underneath** lists every flagged position and jumps the plot to it.

With variable-length sequences, the tail positions that too few sequences reach are not scored and are left out of the figures.

### The leakage report

Full treatment, including how similarity is computed: [train/test leakage](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/guide/leakage/).

`evaluate-splits` searches every test sequence against the training set with MMseqs2 and
flags one when it exceeds `--similarity-threshold` (90% by default). The report gives the
percentage of leaked queries and targets, a histogram of the similarity distribution, and
a panel listing the leaked pairs — up to the first 100, each expanding into its rendered
alignment, so you can see what is actually shared. Every hit is exported to
`mmseqs/mmseqs2_search_result.tsv` beside the report. The bundled enhancers example carries a little real leakage: 0.67% of queries and
0.33% of targets — [see the report](https://genomic-benchmarks.github.io/genomic-benchmarks-qc/reports/enhancers/split/sequence/enhancers_train_vs_enhancers_test/gb-qc-report.html).

## Input Formats

| Format | Description |
|--------|-------------|
| `.fa` / `.fasta` | One or more FASTA files. Multiple files = multiple classes. |
| `.csv` / `.tsv` | One or more CSV/TSV files. Multiple will be pooled and evaluated as one. |

All of the input formats are supported in `.gz` version as well.

**CSV/TSV** needs a sequence column and a label column. By default they are called `sequence` and `label`:

```csv
sequence,label
GAGTGTATGTGTCGAGGAATGTATCCATTTCTTCTAGATTTTCTAGTTT,1
GGTCACCACCACCAAGTTCATGCCTGAACCCTTCAGTGGTCCTTTGCCC,0
```

Override the names with `--sequence-column` and `--label-column`. Any other columns are ignored.

**FASTA** files carry no labels, so each file is one class and its filename stem (the name without extensions) becomes the label — `coding_seqs.fasta` becomes the label `coding_seqs`:

```bash
gb-qc evaluate-classes \
  --input examples/fasta-classes/data/coding_seqs.fasta \
  --input examples/fasta-classes/data/intergenomic_seqs.fasta \
  --out-folder qc-out
```

**Multiple sequence columns**, for datasets where one row pairs two sequences: `evaluate-classes` analyzes each column separately, then concatenates all sequences for a combined `merged` analysis. `evaluate-splits` only analyzes concatenated sequences.

```bash
gb-qc evaluate-classes \
  --input examples/paired-sequences/data/miRNA_mRNA_pairs_dataset.tsv \
  --sequence-column gene \
  --sequence-column noncodingRNA \
  --label-column label \
  --out-folder qc-out
```

**Continuous labels:** `--regression` splits the label column at its median into `high`
and `low` classes. Rows whose value is not numeric are dropped with a warning, and if the
split leaves only one class — a constant column, or one that is mostly zeros — the run
exits with an error.

> **Repeat the option name for each value.** `--input a.csv --input b.csv`, not `--input a.csv b.csv`. The same goes for `--sequence-column` and `--label-list`.

## Output Files

Each command writes into its own sub-directory of `--out-folder` — `class/` for
`evaluate-classes` and `split/` for `evaluate-splits` — then into one directory
per sequence column, then into one directory per comparison:

```text
<out-folder>/
├── class/
│   └── sequence/
│       └── negatives_vs_positives/
│           ├── gb-qc-report.csv
│           └── gb-qc-report.html
└── split/
    └── sequence/
        └── train_vs_test/
            ├── gb-qc-report.csv
            └── gb-qc-report.html
```

If you are checking several datasets, give each one its own `--out-folder` to keep the results side by side.

<details>
<summary><b>How directory names are chosen</b></summary>

Directory names come from your class labels, sequence-column names, and input file
names, lowercased and stripped of characters that are unsafe on some filesystems.
If two of them would collide, the tool makes them unique and warns you which name
it used. Labels shown *inside* the reports and plots are always the originals.

Classes are always ordered alphabetically by directory name, so the same dataset
gives you the same report paths no matter which order you listed the input files
or `--label-list` values in.
</details>

<details>
<summary><b><code>evaluate-classes</code> — output layout</b></summary>

```text
<out-folder>/class/
└── <column>/                   # one directory per sequence column, named after
    │                           # it. FASTA inputs have no columns and use
    │                           # `sequence`; if you give several columns, an
    │                           # extra `merged` report joins them all together
    ├── <classA>_vs_<classB>/
    │   ├── gb-qc-report.csv
    │   ├── gb-qc-report.html
    │   ├── gb-qc-duplicates.txt
    │   └── plots/
    └── per-class/
        └── <class>.json
```

| File | Description |
|------|-------------|
| `gb-qc-report.csv` | Simple comparison report with Pass/Warning/Fail flags. |
| `gb-qc-report.html` | Standalone HTML report: one file, no external assets. |
| `plots/` | Individual plot images (PNG). |
| `gb-qc-duplicates.txt` | Sequences appearing in both compared classes; written only when there are any. |
| `per-class/<class>.json` | Per-class statistics (count, GC%, length, base/dinucleotide frequencies); written only when `json` is in `--report-types`. |
</details>

<details>
<summary><b><code>evaluate-splits</code> — output layout</b></summary>

```text
<out-folder>/split/
└── <column>/                   # the sequence column that was searched. FASTA
    └── <train>_vs_<test>/      # inputs use `sequence`, and several columns
        ├── gb-qc-report.csv    # searched together use `merged`
        ├── gb-qc-report.html
        ├── plots/
        └── mmseqs/
```

| File | Description |
|------|-------------|
| `gb-qc-report.csv` | Simple leakage summary. |
| `gb-qc-report.html` | Interactive HTML report with alignments. |
| `mmseqs/` | Raw MMseqs2 results and filtered FASTA files. |
| `plots/` | Similarity distribution plots. |

Temporary MMseqs2 files go into a `gb-qc-mmseqs-*/` directory inside the
comparison directory and are removed at the end unless you pass
`--keep-tmp-files`. Each run gets its own directory, so several runs can share
one `--out-folder` at the same time. With `--keep-tmp-files` nothing is cleaned
up and the path of each directory is written to the log.
</details>

## CLI Reference

Run `gb-qc evaluate-classes --help` or `gb-qc evaluate-splits --help` for the same thing in your terminal.

<details>
<summary><b><code>evaluate-classes</code> options</b></summary>

| Option | Default | Description |
|--------|---------|-------------|
| `--input` | required | Input file(s)  |
| `--sequence-column` | `sequence` | Column name(s) with sequences |
| `--label-column` | `label` | Column with class labels |
| `--label-list` | `infer` | Specific labels or `infer` |
| `--regression` | False | Treat label as continuous, split high/low |
| `--out-folder` | `.` | Output directory; reports go into `<out-folder>/class/` |
| `--report-types` | `html simple` | `json`, `html`, `simple` |
| `--plot-type` | `boxen` | `boxen` or `violin` |
| `--end-position` | auto | Last position the per-position checks reach. Defaults to the last position at least 50 of each class's sequences reach |
| `--min-coverage` | `0.25` | Fraction of each class that must reach a position before it can be flagged, on top of the 250 sequences every compared position needs |
| `--log-level` | `INFO` | `DEBUG`, `INFO`, `WARNING`, `ERROR` |
| `--log-file` | none | Path to also write logs to; logs go to the console only when unset |

**Tuning the per-position window.** The two options control different things:
`--min-coverage` how far a position can be *flagged*, `--end-position` how far
the checks run at all. A position is flagged only where at least
`max(250, --min-coverage x class size)` sequences in each class reach it;
positions past that are reported as `Unknown` and are not drawn.
`--min-coverage 0` leaves only the 250, which cannot be switched off. Pass
`--end-position` to trim a long fixed-length window (where the default trims
nothing, because every sequence reaches every position), to keep the report
smaller, or to pin the same window across runs so two reports line up position by
position. It can only narrow what gets flagged, never widen it.
</details>

<details>
<summary><b><code>evaluate-splits</code> options</b></summary>

| Option | Default | Description |
|--------|---------|-------------|
| `--train-input` | required | Training file(s) |
| `--test-input` | required | Test file(s) |
| `--sequence-column` | `sequence` | Column with sequences |
| `--out-folder` | `.` | Output directory; reports go into `<out-folder>/split/` |
| `--report-types` | `html simple` | Report formats |
| `--similarity-threshold` | `90.0` | % similarity for leakage flag |
| `--threads` | auto | MMseqs2 thread count |
| `--split-memory-limit` | unlimited | Upper RAM limit for MMseqs2 prefilter structures (e.g. `10G`, `1T`) |
| `--keep-tmp-files` | False | Keep temp files for debugging |
| `--log-level` | `INFO` | `DEBUG`, `INFO`, `WARNING`, `ERROR` |
| `--log-file` | none | Path to also write logs to; logs go to the console only when unset |
</details>

## Contributions & Support

Contributions and suggestions for new features are welcome, as are bug reports! Please create a new [issue](https://github.com/genomic-benchmarks/genomic-benchmarks-qc/issues/new) for any of these, including example reports where possible. Pull-requests for fixes and additions are very welcome. See the [contributing notes](CONTRIBUTING.md) for more information about how the process works.

## Citation

If you use `genomic-benchmarks-qc` in your research, please cite this repository: <https://github.com/genomic-benchmarks/genomic-benchmarks-qc>

<!-- TODO: replace with the paper reference / BibTeX once published. -->

## License

MIT-style. See [LICENSE](LICENSE).
