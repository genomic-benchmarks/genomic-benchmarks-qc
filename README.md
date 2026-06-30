![](https://github.com/katarinagresova/GenBenchQC/blob/main/assets/logo_with_text_transparent.png?raw=True)

Automated quality control for genomic machine learning datasets. Detects biases, inconsistencies, and data leakage before model training.

## Installation

```bash
pip install genbenchQC
```

For `evaluate-splits`, install [mmseqs2](https://github.com/soedinglab/MMseqs2). While conda installation works, precompiled binaries or compiling from source is recommended for better performance and to avoid alignment issues. See the [mmseqs2 installation guide](https://github.com/soedinglab/MMseqs2/wiki#installation) for details.

## Quick Start

### Evaluate Classes - Check Class Bias

Compare sequence characteristics across different classes/labels to detect bias.

**FASTA files (one per class):**
```bash
genbenchQC evaluate-classes \
  --input example_datasets/G4_positives.fasta \
  --input example_datasets/G4_negatives.fasta \
  --out-folder example_outputs/G4_dataset
```

**Multiple CSV/TSV files:**
```bash
genbenchQC evaluate-classes \
  --input example_datasets/enhancers_train.csv \
  --input example_datasets/enhancers_test.csv \
  --out-folder example_outputs/enhancers_dataset
```

Files are pooled; classes inferred from label column.

**CSV/TSV with label column and multiple sequence columns:**
```bash
genbenchQC evaluate-classes \
  --input example_datasets/miRNA_mRNA_pairs_dataset.tsv \
  --sequence-column gene \
  --sequence-column noncodingRNA \
  --label-column label
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
genbenchQC evaluate-splits \
  --train-input example_datasets/enhancers_train.csv \
  --test-input example_datasets/enhancers_test.csv \
  --sequence-column sequence \
  --out-folder example_outputs/enhancers_dataset
```

## Input Formats

| Format | Description |
|--------|-------------|
| `fasta` | One or more FASTA files. Multiple files = multiple classes. |
| `csv` / `tsv` | One or more CSV/TSV files. Multiple will be pooled and avaluated as one. |
| `csv.gz` / `tsv.gz` | Compressed CSV/TSV. |

**Multiple sequence columns:** `evaluate-classes` analyzes each column separately, then concatenates all sequences for combined analysis. `evaluate-splits` only analyzes concatenated sequences.

## Output Files

### evaluate-classes outputs

| File | Description |
|------|-------------|
| `*_report.json` | Per-class statistics (count, GC%, length, base/dinucleotide frequencies). |
| `dataset_report*.csv` | Simple comparison report with Pass/Warning/Fail flags. |
| `dataset_report*.html` | Interactive HTML report with visualizations. |
| `dataset_report*_plots/` | Individual plot images (PNG). |
| `*_duplicates.txt` | Sequences appearing in multiple classes. |

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

| File | Description |
|------|-------------|
| `split_check_*.csv` | Simple leakage summary. |
| `split_check_*.html` | Interactive HTML report with alignments. |
| `split_check_*_mmseqs/` | Raw MMseqs2 results and filtered FASTA files. |
| `split_check_*_plots/` | Similarity distribution plots. |

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

Contributions and suggestions for new features are welcome, as are bug reports! Please create a new [issue](https://github.com/katarinagresova/GenBenchQC/issues/new) for any of these, including example reports where possible. Pull-requests for fixes and additions are very welcome. See the [contributing notes](CONTRIBUTING.md) for more information about how the process works.

## License

MIT-style. See [LICENSE](LICENSE).
