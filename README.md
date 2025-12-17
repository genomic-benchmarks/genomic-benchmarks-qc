![](https://github.com/katarinagresova/GenBenchQC/blob/main/assets/logo_with_text_transparent.png?raw=True)

# Automated Quality Control for Genomic Machine Learning Datasets

GenBenchQC is a Python package and CLI toolkit for automated quality control of genomic datasets used in machine learning.
It helps detect biases, inconsistencies, and potential data leakage across sequences, dataset classes, and train-test splits — ensuring your datasets are reliable before model training.

## Features

### Provided Tools
- **genbenchQC evaluate-classes** – QC tool to evaluate sequence characteristics between different classes/labels in the dataset.
- **genbenchQC evaluate-splits** – QC tool to evaluate data leakage in dataset train-test splits.

### General Features
- [**Class-level QC**](https://github.com/katarinagresova/GenBenchQC/tree/main?tab=readme-ov-file#evaluate-classes) – Compare multiple classes for feature similarity or bias.
- [**Train–test split QC**](https://github.com/katarinagresova/GenBenchQC/tree/main?tab=readme-ov-file#evaluate-splits) – Detect potential data leakage through sequence similarity and clustering.
- [**Multiple input formats**](https://github.com/katarinagresova/GenBenchQC/tree/main?tab=readme-ov-file#supported-input-file-formats) – Supports FASTA, CSV, and TSV datasets.
- **Customizable reporting** – Generate JSON, HTML, or simple text summaries.
- **Integration-ready** – Available as both CLI tools and a Python API.
- **Flexible sequence handling** – Works with single or multiple sequence columns.

## Installation

Install Genomic Benchmarks QC using pip:

```bash
pip install genbenchQC
```

If you plan to use `evaluate-splits`, install [cd-hit](https://www.bioinformatics.org/cd-hit/cd-hit-user-guide):

```bash
conda install -c bioconda cd-hit
# or follow: https://github.com/weizhongli/cdhit/wiki/2.-Installation
```

## Quick Start

Clone the repository to access example datasets:

```bash
git clone https://github.com/katarinagresova/GenBenchQC.git
cd GenBenchQC
```

### Evaluate Classes

Running from CLI with fasta file:

```bash
genbenchQC evaluate-classes \
  --input example_datasets/G4_positives.fasta \
  --input example_datasets/G4_negatives.fasta \
  --format fasta \
  --out-folder example_outputs/G4_dataset
```

Outputs with their description are in [example_outputs/G4_dataset](https://github.com/katarinagresova/GenBenchQC/tree/main/example_outputs/G4_dataset).

Running from CLI with tsv file and two sequence columns:

```bash
genbenchQC evaluate-classes \
  --input example_datasets/miRNA_mRNA_pairs_dataset.tsv \
  --format tsv \
  --out-folder example_outputs/miRNA_mRNA_dataset \
  --sequence-column gene \
  --sequence-column noncodingRNA
```

Note: when you want to provide multiple values for some option, such as `--input` or `--sequence-column`, prefix each value with option name:
```bash
genbenchQC evaluate-classes \
  --input example_datasets/G4_positives.fasta \
  --input example_datasets/G4_negatives.fasta 
```

Outputs with their description are in [example_outputs/miRNA_mRNA_dataset](https://github.com/katarinagresova/GenBenchQC/tree/main/example_outputs/miRNA_mRNA_dataset).

### Evaluate Splits

```bash
genbenchQC evaluate-splits \
  --train-input example_datasets/enhancers_train.csv \
  --test-input example_datasets/enhancers_test.csv \
  --format csv \
  --sequence-column sequence \
  --out-folder example_outputs/enhancers_dataset
```

Outputs with their description are in [example_outputs/enhancers_dataset](https://github.com/katarinagresova/GenBenchQC/tree/main/example_outputs/enhancers_dataset).

## Supported input file formats

You can choose to run the tools while having different dataset formats:
- **FASTA**: The input is a FASTA file / list of FASTA files. For *evaluate-classes* each fasta file is treated as separate class/label.
- **CSV/TSV**: The input is a CSV/TSV file, and you provide the name of the column containing sequences. You can have either:
  - **multiple files**, each one containing sequences from one class (similar as with FASTA input)
  - **one file** containing sequences from multiple classes. In this case, when running *evaluate-classes* tool, you need to provide the name of the column containing class labels so the tool can split the dataset into parts. The label classes can then be inferred, or you can specify their list by yourself. The dataset will then be split into pieces containing sequences with corresponding labels and analysis will be performed similarly as with multiple files.
- **CSV.GZ/TSV.GZ**: Functionality is the same as CSV/TSV files

When having CSV/TSV/CSV.GZ/TSV.GZ input, you can also decide to provide multiple sequence columns to analyze. In this case, the tool *evaluate-classes* will be performed for each column separately and lastly for sequences made by concatenating sequences throughout all the columns. 
*evaluate-splits* tool will run only the concatenated sequences.


## Development

If you want to help with the development of Genomic Benchmarks QC, you are more than welcome to join in!

For a guidance, have a look at [CONTRIBUTING.md](https://github.com/katarinagresova/GenBenchQC/blob/main/CONTRIBUTING.md)

## License

Genomic Benchmarks QC is MIT-style licensed, as found in the [LICENSE](https://github.com/katarinagresova/GenBenchQC/blob/main/LICENSE) file.