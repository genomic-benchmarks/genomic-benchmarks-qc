<img src="https://github.com/genomic-benchmarks/genomic-benchmarks-qc/blob/main/assets/genomic-benchmarks-qc-text.png?raw=True" alt="genomic-benchmarks-qc logo" width="500" />

[![CI](https://github.com/genomic-benchmarks/genomic-benchmarks-qc/actions/workflows/ci.yml/badge.svg)](https://github.com/genomic-benchmarks/genomic-benchmarks-qc/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/genomic-benchmarks-qc)](https://pypi.org/project/genomic-benchmarks-qc/)
[![Python](https://img.shields.io/pypi/pyversions/genomic-benchmarks-qc)](https://pypi.org/project/genomic-benchmarks-qc/)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

Automated quality control for genomic machine learning datasets. Detects biases, inconsistencies, and data leakage before model training.

<!-- TODO(screenshot): add the file, then uncomment the <img> below.
      assets/report-summary.png — top of gb-qc-report.html: the Summary
     table with its Pass/Warning/Fail flags, plus the first plot below it. ~1200px wide. -->
<!-- <img src="assets/report-summary.png" alt="Example gb-qc HTML report" width="800" /> -->

**[Browse the example reports →](example_outputs/)** (download a `gb-qc-report.html` and open it in a browser — GitHub shows the source, not the page)

## What it catches

A model that scores well on a biased dataset has learned the bias, not the biology. `gb-qc` looks for the differences a classifier could exploit without understanding anything:

- **Your negatives are shorter than your positives.** A length classifier now beats your model.
- **Your classes differ in GC content, base composition or dinucleotide frequencies.** In the bundled enhancers example, GC content alone separates the two classes at AU-ROC 0.66.
- **Position 1 is `N` in one class only**, or every sequence in one class starts with the same adapter — a per-position give-away that no summary statistic would show.
- **The same sequence appears in both classes** — reported exactly — **or your test set repeats your training set.** `evaluate-splits` catches near-duplicates too, with an MMseqs2 similarity search.

Each check gets a **Pass / Warning / Fail** flag, an HTML report you can read, and a CSV you can put in CI.

## Installation

```bash
pip install genomic-benchmarks-qc
```

Requires Python 3.12 or newer.

<details>
<summary><b>Extra step for <code>evaluate-splits</code>: mmseqs2</b></summary>

The leakage check shells out to [mmseqs2](https://github.com/soedinglab/MMseqs2), which is not a Python dependency and must be installed separately. While conda installation works, precompiled binaries or compiling from source is recommended for better performance and to avoid alignment issues. See the [mmseqs2 installation guide](https://github.com/soedinglab/MMseqs2/wiki#installation) for details.

`evaluate-classes` does not need it.
</details>

## Quick Start

Point the tool at your dataset and give it somewhere to write:

```bash
gb-qc evaluate-classes \
  --input example_datasets/enhancers_train.csv \
  --input example_datasets/enhancers_test.csv \
  --out-folder example_outputs/enhancers_dataset
```

It tells you what it is doing as it goes, including which positions it could and could not compare:

```text
INFO - Starting classes evaluation.
INFO - Merging 2 input files: enhancers_train.csv, enhancers_test.csv
INFO - Computing statistics for merged, label 0, sequence column: sequence
INFO - End position argument not provided. Per-position checks cover positions 1-483 for
       sequence comparison, which is as far as at least 50 sequences reach.
INFO - Positions 1-359 for sequence comparison may be flagged, which is as far as 250
       sequences reach.
INFO - Positions 360-483 for sequence comparison are reported as Unknown and are not drawn:
       too few sequences reach them to say whether a difference there is a difference
       between the classes or between their longest sequences.
INFO - Comparing classes: 0 vs 1
INFO - Generating simple report: .../class/sequence/0_vs_1/gb-qc-report.csv
INFO - Generating PNG plots at: .../class/sequence/0_vs_1/plots
INFO - Classes evaluation successfully completed.
```

and leaves you a report to open:

```text
example_outputs/enhancers_dataset/class/sequence/0_vs_1/
├── gb-qc-report.html   ← open this
├── gb-qc-report.csv
└── plots/
```

To check a train/test split for leakage instead:

```bash
gb-qc evaluate-splits \
  --train-input example_datasets/enhancers_train.csv \
  --test-input example_datasets/enhancers_test.csv \
  --sequence-column sequence \
  --out-folder example_outputs/enhancers_dataset
```

Both commands can share one `--out-folder`; they write into `class/` and `split/` respectively and never overwrite each other.

## Reading the report

The HTML report is a single standalone file — no external assets, so you can mail it or archive it. It opens with a summary table of every check and its flag, followed by the basic descriptive statistics of each class, then one section per check with the numbers and the plot behind it.

Flags come from the AU-ROC of a classifier that sees only that one feature:

| Flag | AU-ROC | Meaning |
|------|--------|---------|
| **Pass** | ≤ 0.6 | Classes not distinguishable by this feature |
| **Warning** | ≤ 0.7 | A model could get some traction here |
| **Fail** | > 0.7 | Significant bias detected |
| **Unknown** | — | Not enough sequences to score the check |

**`Unknown` is not `Pass`**: it says the comparison was not made, not that it came out clean. A check needs at least **250 sequences** — and, for the per-position checks, a cohort of at least **25% of the class** — before it is scored at all. The plots, per-class statistics and descriptive tables are computed from all the data either way, so small datasets can still be compared by eye, and the terminal says which checks were skipped and why.

<details>
<summary><b>Why 250 sequences and 25% coverage</b></summary>

Each check reports the worst case over its features, and the worst case over many weak features crosses a fixed boundary on sampling noise alone when the classes are small. Every comparison therefore has to clear the same two floors before it is made at all:

- At least **250 sequences**. For the per-sequence checks that is a floor on the
  smaller class; for the per-position checks it is a floor on the cohort reaching
  each position, since a position is compared only on the sequences long enough to
  have it. Below 250 the comparison reports `Unknown` rather than a verdict: at
  100 sequences per class the dinucleotide check flags two classes drawn from the
  same process 19.4% of the time, against 0.2% at 250. The per-position checks
  reach 0.0% there, which is what lets them use the same fixed boundaries as
  everything else rather than a threshold that moves with each cohort.
- At least **25% of the class** (`--min-coverage`), for the per-position checks.
  A cohort out in the tail can clear 250 many times over and still not stand for
  the class: it is all of the class's long sequences, and a difference there can
  be a difference between those subsets rather than between the classes. No
  sample size fixes that; only stopping does.

The count binds on small and mid-sized classes, where sampling noise is the risk;
the share binds above 1000 sequences per class, where the tail is long enough for
a large cohort to be unrepresentative. Both floors were set from null and power
simulations of the checks.

The vocabulary and duplication checks are exempt from all of this. They ask
whether something *is present* — a base in one class and not the other, a
sequence in both classes — rather than whether a model could exploit it, and a
single occurrence is worth reporting however small the dataset.
</details>

### The per-position plot is interactive

A flagged position is one pixel wide in a plot spanning hundreds of positions, so the per-position panels are drawn in the browser from the numbers behind them rather than embedded as an image.

<!-- TODO(gif): add the file, then uncomment the <img> below.
      assets/per-position-demo.gif — 5-10s screen capture of the Per Position
     Nucleotide Content panel: drag to zoom into a flagged region, hover to show the
     tooltip, press `n` to step to the next flagged position. ~900px wide, no audio. -->
<!-- <img src="assets/per-position-demo.gif" alt="Interactive per-position plot" width="800" /> -->

- **Drag** across the plot to zoom.
- **Hover** for the frequencies and AU-ROC at a position.
- **`n` / `p`** or the buttons step between flagged positions.
- **Save view** writes the current window as a PNG.
- The **table underneath** lists every flagged position and jumps the plot to it.

Static PNGs of the same panels are still written to `plots/`. The data is embedded per position rather than per sequence, so the report's size follows the length of the analysed window, not the size of the dataset.

<details>
<summary><b>Where the per-position figures stop</b></summary>

They are drawn over the scored window and no further, so every position in them was
compared and an unflagged stretch is a stretch that passed, with nothing on the figure to
read past. Positions the sequences reach but the comparison could not use are still
reported as `Unknown` — in `gb-qc-report.csv` and in the count of checks that
were not scored — and the section's explanation says how far they run and why
they were left out. A comparison where *no* position could be scored is the one exception:
there the figures fall back to the full reported window, as every other plot in
an underpowered report is still drawn, with every position `Unknown` and no flag
on any of them. The panel under each figure shows how many sequences stand behind each
position and marks the floor they have to stay above: the window ends where the
lower curve crosses it.
</details>

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
  --input example_datasets/coding_seqs.fasta \
  --input example_datasets/intergenomic_seqs.fasta \
  --out-folder example_outputs/coding_vs_intergenomic_dataset
```

**Multiple sequence columns**, for datasets where one row pairs two sequences: `evaluate-classes` analyzes each column separately, then concatenates all sequences for a combined `merged` analysis. `evaluate-splits` only analyzes concatenated sequences.

```bash
gb-qc evaluate-classes \
  --input example_datasets/miRNA_mRNA_pairs_dataset.tsv \
  --sequence-column gene \
  --sequence-column noncodingRNA \
  --label-column label \
  --out-folder example_outputs/miRNA_mRNA_dataset
```

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

Every report ends up at the same kind of path whatever the input format was, and because each comparison gets a directory of its own, the files inside it always have the same names. If you are checking several datasets, give each one its own `--out-folder` to keep the results side by side.

<details>
<summary><b>How directory names are chosen</b></summary>

Directory names come from your class labels, sequence-column names, and input file
names. They are lowercased and reduced to characters that are safe on any filesystem.
If two of them would
produce the same directory name, the tool makes them unique and warns you which
name it used. Labels shown *inside* the reports and plots are always the
originals.

Classes are always ordered alphabetically by directory name, so the same dataset
gives you the same report paths no matter which order you listed the input files
or `--label-list` values in.
</details>

<details>
<summary><b><code>evaluate-classes</code> — full output layout and checks</b></summary>

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

**Evaluated features:**
- Nucleotide vocabulary
- Sequence lengths distribution
- GC content per sequence
- Nucleotide composition (per sequence & per position)
- Dinucleotide frequencies
- Sequence duplication levels
- Exact duplicate sequences between classes

**Regression targets:** `--regression` splits the label column at its median into
`high` and `low` classes. Rows whose value is not numeric are dropped with a
warning. If no numeric values are left, or if every value falls on one side of
the median — which happens with a constant column, or one that is mostly
zeros — there are no two classes to compare, so the tool tells you why and exits
with an error instead of producing a meaningless report.
</details>

<details>
<summary><b><code>evaluate-splits</code> — full output layout and leakage detection</b></summary>

```text
<out-folder>/split/
└── <column>/                   # the sequence column that was searched. FASTA
    └── <train>_vs_<test>/      # inputs use `sequence`, and several columns
        ├── gb-qc-report.csv    # searched together use `merged`
        ├── gb-qc-report.html
        ├── plots/
        └── mmseqs/
```

If you give several `--sequence-column` values, they are joined per row and
searched as one sequence rather than one column at a time, so the run produces a
single `merged/` comparison rather than one per column.

| File | Description |
|------|-------------|
| `gb-qc-report.csv` | Simple leakage summary. |
| `gb-qc-report.html` | Interactive HTML report with alignments. |
| `mmseqs/` | Raw MMseqs2 results and filtered FASTA files. |
| `plots/` | Similarity distribution plots. |

**Leakage detection:** every test sequence is searched against the training set with
MMseqs2, and one is flagged when it exceeds `--similarity-threshold` (90% by default).
The report gives the percentage of leaked queries and targets, a histogram of the
similarity distribution, and the rendered alignments of the 100 most similar pairs so
you can see what is actually shared. The bundled enhancers example carries a little real
leakage: 0.67% of queries and 0.33% of targets.

Temporary MMseqs2 files are kept in a `gb-qc-mmseqs-*/` directory inside the
comparison directory and removed at the end unless you pass `--keep-tmp-files`.
The trailing suffix is random, so each run gets a directory of its own: several
comparisons — or even repeated runs of the same one — can go into one
`--out-folder` at the same time without disturbing each other's scratch files.
With `--keep-tmp-files`, each run therefore leaves its own directory behind and
none are cleaned up later; the path of each is written to the log.
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

**Tuning the per-position window.** `--min-coverage` and `--end-position` control two
different things: how far a position can be *flagged*, and how far the checks *run* at all.
A position is flagged only where at least `max(250, --min-coverage x class size)` sequences
in each class reach it; positions past that are reported as `Unknown` and are not drawn.
`--min-coverage 0` leaves only the 250, which is the guard against sampling noise and
cannot be switched off. `--end-position` is worth passing to trim a long
fixed-length window (where the default trims nothing, because every sequence
reaches every position), to hold the report and its embedded per-position data to
a size, or to pin the same window across runs so that two reports line up
position by position. It can only narrow what gets flagged, never widen it.
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
