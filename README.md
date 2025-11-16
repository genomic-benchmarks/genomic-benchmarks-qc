![](assets/logo_with_text_transparent.png)

# Automated Quality Control for Genomic Machine Learning Datasets

GenBenchQC is a Python package and CLI toolkit for automated quality control of genomic datasets used in machine learning.
It helps detect biases, inconsistencies, and potential data leakage across sequences, dataset classes, and train-test splits — ensuring your datasets are reliable before model training.

## Features

### Provided Tools
- **evaluate_sequences** – QC of a single dataset or dataset subset.
- **evaluate_dataset** – QC across multiple dataset classes.
- **evaluate_split** – Train–test split leakage detection.

### General Features
- [**Sequence-level QC**](https://github.com/katarinagresova/GenBenchQC/tree/main?tab=readme-ov-file#evaluate-sequences) – Evaluate nucleotide composition, sequence length distribution, GC content, and more.
- [**Class-level QC**](https://github.com/katarinagresova/GenBenchQC/tree/main?tab=readme-ov-file#evaluate-dataset) – Compare multiple classes for feature similarity or bias.
- [**Train–test split validation**](https://github.com/katarinagresova/GenBenchQC/tree/main?tab=readme-ov-file#evaluate-split) – Detect potential data leakage through sequence similarity and clustering.
- [**Multiple input formats**](https://github.com/katarinagresova/GenBenchQC/tree/main?tab=readme-ov-file#supported-input-file-formats) – Supports FASTA, CSV, and TSV datasets.
- **Customizable reporting** – Generate JSON, HTML, or simple text summaries.
- **Integration-ready** – Available as both CLI tools and a Python API.
- **Flexible sequence handling** – Works with single or multiple sequence columns.

## Installation

Install Genomic Benchmarks QC using pip:

```bash
pip install genbenchQC
```

If you plan to use `evaluate_split`, install [cd-hit](https://www.bioinformatics.org/cd-hit/cd-hit-user-guide):

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

### Evaluate Sequences

```bash
evaluate_sequences \
  --input example_datasets/G4_positives.fasta \
  --format fasta \
  --out_folder example_outputs/G4_dataset_positives
```

Outputs with their description are in [example_outputs/G4_dataset_positives](https://github.com/katarinagresova/GenBenchQC/tree/main/example_outputs/G4_dataset_positives).

### Evaluate Dataset

Running from CLI with fasta file:

```bash
evaluate_dataset \
  --input example_datasets/G4_positives.fasta example_datasets/G4_negatives.fasta \
  --format fasta \
  --out_folder example_outputs/G4_dataset
```

Outputs with their description are in [example_outputs/G4_dataset](https://github.com/katarinagresova/GenBenchQC/tree/main/example_outputs/G4_dataset).

### Evaluate Split

```bash
evaluate_split \
  --train_input example_datasets/enhancers_train.csv \
  --test_input example_datasets/enhancers_test.csv \
  --format csv \
  --sequence_column sequence \
  --out_folder example_outputs/enhancers_dataset
```

Outputs with their description are in [example_outputs/enhancers_dataset](https://github.com/katarinagresova/GenBenchQC/tree/main/example_outputs/enhancers_dataset).

## Python Execution

The same commands can be executed from Python.

### Evaluate Sequences

```python
from genbenchQC import evaluate_sequences

evaluate_sequences.run(
  input='example_datasets/G4_positives.fasta', 
  format='fasta',
  out_folder='example_outputs/G4_dataset_positives'
)
```

### Evaluate Dataset

Running from Python with CSV file with multiple sequence columns:

```python
from genbenchQC import evaluate_dataset

evaluate_dataset.run(
  input=['example_datasets/miRNA_mRNA_pairs_dataset.tsv'], 
  format='tsv', 
  out_folder='example_outputs/miRNA_mRNA_dataset', 
  sequence_column=['gene', 'noncodingRNA'], 
  label_column='label', 
  label_list=['0', '1']
)
```

Outputs with their description are in [example_outputs/miRNA_mRNA_dataset](https://github.com/katarinagresova/GenBenchQC/tree/main/example_outputs/miRNA_mRNA_dataset).

### Evaluate Split

```python
from genbenchQC import evaluate_split

evaluate_split.run(
  train_files=['example_datasets/enhancers_train.csv'],
  test_files=['example_datasets/enhancers_test.csv'],
  format='csv',
  sequence_column='sequence',
  out_folder='example_outputs/enhancers_dataset'
)
```

## Supported input file formats

You can choose to run the tool while having different dataset formats:
- **FASTA**: The input is a FASTA file / list of FASTA files. One file needs to contain sequences of one class if running *evaluate_sequences* mode.
- **CSV/TSV**: The input is a CSV/TSV file, and you provide the name of the column containing sequences. You can have either:
  - **multiple files**, each one containing sequences from one class (similar as with FASTA input)
  - **one file** containing sequences from multiple classes. In this case, when running *evaluate_sequences* mode, you need to provide the name of the column containing class labels so the tool can split the dataset into parts. The label classes can then be inferred, or you can specify their list by yourself. The dataset will then be split into pieces containing sequences with corresponding labels and analysis will be performed similarly as with multiple files.
- **CSV.GZ/TSV.GZ**: Functionality is the same as CSV/TSV files

When having CSV/TSV/CSV.GZ/TSV.GZ input, you can also decide to provide multiple sequence columns to analyze. In this case, the analysis in modes *evaluate_sequences* and *evaluate_dataset* will be performed for each column separately and lastly for sequences made by concatenating sequences throughout all the columns. 
*evaluate_split* mode will run only the concatenated sequences.


## Development

If you want to help with the development of Genomic Benchmarks QC, you are more than welcome to join in!

For a guidance, have a look at [CONTRIBUTING.md](https://github.com/katarinagresova/GenBenchQC/blob/main/CONTRIBUTING.md)

## License

Genomic Benchmarks QC is MIT-style licensed, as found in the [LICENSE](https://github.com/katarinagresova/GenBenchQC/blob/main/LICENSE) file.