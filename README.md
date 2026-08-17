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

### evaluate-classes outputs

Report files are named after the compared classes, using the stem
`evaluate-classes[_col_<column>]_label_<classA>_vs_<classB>` (the `_col_<column>`
part is only added for CSV/TSV inputs with a sequence column). For example:
`evaluate-classes_label_G4_positives_vs_G4_negatives` or
`evaluate-classes_col_sequence_label_0_vs_1`.

| File | Description |
|------|-------------|
| `evaluate-classes_*.csv` | Simple comparison report with Pass/Warning/Fail flags. |
| `evaluate-classes_*.html` | Interactive HTML report with visualizations. |
| `evaluate-classes_*_plots/` | Individual plot images (PNG). |
| `evaluate-classes_*_duplicates.txt` | Sequences appearing in multiple classes. |
| `evaluate-classes_*_report.json` | Per-class statistics (count, GC%, length, base/dinucleotide frequencies); written only when `json` is in `--report-types`. |

**Flag thresholds (AU-ROC):**
- **Pass**: ≤ 0.6 (classes not distinguishable by feature)
- **Warning**: ≤ 0.7
- **Fail**: > 0.7 (significant bias detected)

**Evaluated features:**
- Nucleotide vocabulary
- Sequence lengths distribution
- GC content per sequence
- Nucleotide composition (per sequence & per position)
- Dinucleotide frequencies
- Sequence duplication levels
- Exact duplicate sequences between classes

### evaluate-splits outputs

Files use the stem `evaluate-splits_split_<train>_vs_<test>`, e.g.
`evaluate-splits_split_enhancers_train_vs_enhancers_test`.

| File | Description |
|------|-------------|
| `evaluate-splits_*.csv` | Simple leakage summary. |
| `evaluate-splits_*.html` | Interactive HTML report with alignments. |
| `evalaute-splits_*_mmseqs/` | Raw MMseqs2 results and filtered FASTA files. |
| `evalaute-splits_*_plots/` | Similarity distribution plots. |

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
| `--out-folder` | `.` | Output directory |
| `--report-types` | `html simple` | `json`, `html`, `simple` |
| `--plot-type` | `boxen` | `boxen` or `violin` |
| `--end-position` | auto | Max position for per-position stats |
| `--log-level` | `INFO` | `DEBUG`, `INFO`, `WARNING`, `ERROR` |

### evaluate-splits

| Option | Default | Description |
|--------|---------|-------------|
| `--train-input` | required | Training file(s) |
| `--test-input` | required | Test file(s) |
| `--sequence-column` | `sequence` | Column with sequences |
| `--out-folder` | `.` | Output directory |
| `--report-types` | `html simple` | Report formats |
| `--similarity-threshold` | `90.0` | % similarity for leakage flag |
| `--threads` | auto | MMseqs2 thread count |
| `--keep-tmp-files` | False | Keep temp files for debugging |

## Contributions & Support

Contributions and suggestions for new features are welcome, as are bug reports! Please create a new [issue](https://github.com/genomic-benchmarks/genomic-benchmarks-qc/issues/new) for any of these, including example reports where possible. Pull-requests for fixes and additions are very welcome. See the [contributing notes](CONTRIBUTING.md) for more information about how the process works.

## License

MIT-style. See [LICENSE](LICENSE).
