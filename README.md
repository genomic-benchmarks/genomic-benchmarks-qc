<img src="https://github.com/genomic-benchmarks/genomic-benchmarks-qc/blob/main/assets/genomic-benchmarks-qc-text.png?raw=True" alt="genomic-benchmarks-qc logo" width="500" />

Automated quality control for genomic machine learning datasets. Detects biases, inconsistencies, and data leakage before model training.

## Installation

```bash
pip install genomic-benchmarks-qc
```

For `evaluate-splits`, install [mmseqs2](https://github.com/soedinglab/MMseqs2). While conda installation works, precompiled binaries or compiling from source is recommended for better performance and to avoid alignment issues. See the [mmseqs2 installation guide](https://github.com/soedinglab/MMseqs2/wiki#installation) for details.

## Quick Start

### Evaluate Classes - Check Class Bias

Compare sequence characteristics across different classes/labels to detect bias.

**FASTA files (one per class):**
```bash
gb-qc evaluate-classes \
  --input example_datasets/G4_positives.fasta \
  --input example_datasets/G4_negatives.fasta \
  --out-folder example_outputs/G4_dataset
```

Each FASTA file represents one class, and its filename stem (the name without extensions) is used as the class/label name — e.g. `G4_positives.fasta` becomes the label `G4_positives`.

**Multiple CSV/TSV files:**
```bash
gb-qc evaluate-classes \
  --input example_datasets/enhancers_train.csv \
  --input example_datasets/enhancers_test.csv \
  --out-folder example_outputs/enhancers_dataset
```

When multiple CSV/TSV files are given, their rows are combined into a single
dataset, and the classes to compare are taken from the values in the label
column (`label` by default; override with `--label-column`).

**CSV/TSV with label column and multiple sequence columns:**
```bash
gb-qc evaluate-classes \
  --input example_datasets/miRNA_mRNA_pairs_dataset.tsv \
  --sequence-column gene \
  --sequence-column noncodingRNA \
  --label-column label \
  --out-folder example_outputs/miRNA_mRNA_dataset
```

Note: when you want to provide multiple values for some option, such as `--input` or `--sequence-column`, prefix each value with option name:
```bash
  --sequence-column gene \
  --sequence-column noncodingRNA \
```

### Evaluate Splits - Check Data Leakage

Detect train-test data leakage via sequence similarity.

```bash
gb-qc evaluate-splits \
  --train-input example_datasets/enhancers_train.csv \
  --test-input example_datasets/enhancers_test.csv \
  --sequence-column sequence \
  --out-folder example_outputs/enhancers_dataset
```

## Input Formats

| Format | Description |
|--------|-------------|
| `.fa/.fasta` | One or more FASTA files. Multiple files = multiple classes. |
| `.csv` / `.tsv` | One or more CSV/TSV files. Multiple will be pooled and evaluated as one. |

All of the input formats are supported in .gz version as well.

**Multiple sequence columns:** `evaluate-classes` analyzes each column separately, then concatenates all sequences for combined analysis. `evaluate-splits` only analyzes concatenated sequences.

## Output Files

Each command writes into its own sub-directory of `--out-folder` — `class/` for
`evaluate-classes` and `split/` for `evaluate-splits` — then into one directory
per sequence column, then into one directory per comparison:

```text
<out-folder>/
├── class/
│   └── sequence/
│       └── negatives_vs_positives/
│           ├── report.csv
│           └── report.html
└── split/
    └── sequence/
        └── train_vs_test/
            ├── report.csv
            └── report.html
```

Every report ends up at the same kind of path whatever the input format was, and
because each comparison gets a directory of its own, the files inside it always
have the same names — `report.html` is the HTML report for whichever comparison
the directory is named after. You can run both commands on the same dataset with
the same `--out-folder` without them overwriting each other; if you are checking
several datasets, give each one its own `--out-folder` to keep the results side
by side.

Directory names come from your class labels, sequence-column names, and input file
names. They are lowercased and reduced to characters that are safe on any filesystem.
If two of them would
produce the same directory name, the tool makes them unique and warns you which
name it used. Labels shown *inside* the reports and plots are always the
originals.

Classes are always ordered alphabetically by directory name, so the same dataset
gives you the same report paths no matter which order you listed the input files
or `--label-list` values in.

### evaluate-classes outputs

```text
<out-folder>/class/
└── <column>/                   # one directory per sequence column, named after
    │                           # it. FASTA inputs have no columns and use
    │                           # `sequence`; if you give several columns, an
    │                           # extra `merged` report joins them all together
    ├── <classA>_vs_<classB>/
    │   ├── report.csv
    │   ├── report.html
    │   ├── duplicates.txt
    │   └── plots/
    └── per-class/
        └── <class>.json
```

| File | Description |
|------|-------------|
| `report.csv` | Simple comparison report with Pass/Warning/Fail flags. |
| `report.html` | Standalone HTML report: one file, no external assets. |
| `plots/` | Individual plot images (PNG). |
| `duplicates.txt` | Sequences appearing in both compared classes; written only when there are any. |
| `per-class/<class>.json` | Per-class statistics (count, GC%, length, base/dinucleotide frequencies); written only when `json` is in `--report-types`. |

**Flag thresholds (AU-ROC):**
- **Pass**: ≤ 0.6 (classes not distinguishable by feature)
- **Warning**: ≤ 0.7
- **Fail**: > 0.7 (significant bias detected)
- **Unknown**: not enough sequences to score the check (see below)

**How small is too small.** Each check reports the worst case over its features,
and the worst case over many weak features crosses a fixed boundary on sampling
noise alone when the classes are small. Two guards keep a flag meaning what it
says:

- The per-sequence checks need at least **200 sequences in the smaller class**.
  Below that they report `Unknown` rather than a verdict, because at 100
  sequences per class the dinucleotide check flags two classes drawn from the
  same process 19.4% of the time, against 1.2% at 200.
- The per-position checks compare each position on the sequences long enough to
  reach it, so their cohorts shrink along the sequence and one class-wide floor
  would silence the far end of every variable-length dataset. Instead, the
  difference a position must show before it counts widens as its own cohort
  shrinks. On a normally sized dataset this falls below 0.6 and the fixed
  boundaries above apply unchanged.

`Unknown` is not `Pass`: it says the comparison was not made, not that it came
out clean. The plots, the per-class statistics and the descriptive tables are
computed from all the data either way, so small datasets can still be compared
by eye, and the terminal says which checks were skipped and why.

The vocabulary and duplication checks are exempt from all of this. They ask
whether something *is present* — a base in one class and not the other, a
sequence in both classes — rather than whether a model could exploit it, and a
single occurrence is worth reporting however small the dataset.

**Evaluated features:**
- Nucleotide vocabulary
- Sequence lengths distribution
- GC content per sequence
- Nucleotide composition (per sequence & per position)
- Dinucleotide frequencies
- Sequence duplication levels
- Exact duplicate sequences between classes

### evaluate-splits outputs

```text
<out-folder>/split/
└── <column>/                   # the sequence column that was searched. FASTA
    └── <train>_vs_<test>/      # inputs use `sequence`, and several columns
        ├── report.csv          # searched together use `merged`
        ├── report.html
        ├── plots/
        └── mmseqs/
```

If you give several `--sequence-column` values, they are joined per row and
searched as one sequence rather than one column at a time, so the run produces a
single `merged/` comparison rather than one per column.

| File | Description |
|------|-------------|
| `report.csv` | Simple leakage summary. |
| `report.html` | Interactive HTML report with alignments. |
| `mmseqs/` | Raw MMseqs2 results and filtered FASTA files. |
| `plots/` | Similarity distribution plots. |

Temporary MMseqs2 files are kept in a `gb-qc-mmseqs-*/` directory inside the
comparison directory and removed at the end unless you pass `--keep-tmp-files`.
The trailing suffix is random, so each run gets a directory of its own: several
comparisons - or even repeated runs of the same one - can go into one
`--out-folder` at the same time without disturbing each other's scratch files.
With `--keep-tmp-files`, each run therefore leaves its own directory behind and
none are cleaned up later; the path of each is written to the log.

**Leakage detection:** Flags when test sequences exceed similarity threshold vs training sequences.

## CLI Options

### evaluate-classes

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
| `--end-position` | auto | Last position included in per-position stats. Defaults to the last position at least 75% of each class's sequences reach |
| `--log-level` | `INFO` | `DEBUG`, `INFO`, `WARNING`, `ERROR` |
| `--log-file` | none | Path to also write logs to; logs go to the console only when unset |

**Regression targets:** `--regression` splits the label column at its median into
`high` and `low` classes. Rows whose value is not numeric are dropped with a
warning. If no numeric values are left, or if every value falls on one side of
the median — which happens with a constant column, or one that is mostly
zeros — there are no two classes to compare, so the tool tells you why and exits
with an error instead of producing a meaningless report.

### evaluate-splits

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

## Contributions & Support

Contributions and suggestions for new features are welcome, as are bug reports! Please create a new [issue](https://github.com/genomic-benchmarks/genomic-benchmarks-qc/issues/new) for any of these, including example reports where possible. Pull-requests for fixes and additions are very welcome. See the [contributing notes](CONTRIBUTING.md) for more information about how the process works.

## License

MIT-style. See [LICENSE](LICENSE).
