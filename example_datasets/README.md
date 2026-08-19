# Example datasets

Small, seeded subsets of published benchmarks, kept here so that the commands in
the top-level README and the reports in [`example_outputs/`](../example_outputs)
can be reproduced from a clone.

Each class holds at least 500 sequences. That is deliberate: below 200 sequences
in the smaller class the per-sequence checks report `Unknown` rather than a
verdict, and an example that cannot be scored is not much of an example.

| File | Source | Size | What it shows |
|------|--------|------|----------------|
| `coding_seqs.fasta`, `intergenomic_seqs.fasta` | `demo_coding_vs_intergenomic_seqs`, [Genomic Benchmarks](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks) | 600 sequences per file, 200 nt | One FASTA file per class. The two classes differ in composition, so most per-sequence checks fail. |
| `enhancers_train.csv`, `enhancers_test.csv` | `human_enhancers_ensembl`, [Genomic Benchmarks](https://github.com/ML-Bioinfo-CEITEC/genomic_benchmarks) | 1,200 and 600 rows, 4–568 nt | A label column, and sequences whose lengths vary widely — so the per-position window and its thinning cohorts are visible. Also the input for `evaluate-splits`, where a little real leakage shows up. |
| `miRNA_mRNA_pairs_dataset.tsv` | `AGO2_CLASH_Hejret2023`, [miRBench](https://github.com/katarinagresova/miRBench) | 1,200 rows | Two sequence columns analysed separately and merged. Sequences repeat within and across classes, so the duplication checks fail. |

Rows were drawn per class with `random_state=0` from each source's training
split, then shuffled. The classes are balanced, which the sources mostly are
already.
