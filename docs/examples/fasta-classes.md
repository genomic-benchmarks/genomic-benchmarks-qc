# fasta-classes

**One FASTA file per class, label taken from the filename.**

600 coding and 600 intergenomic sequences, each 200 nucleotides, from
[Genomic Benchmarks](https://doi.org/10.1186/s12863-023-01123-8)'s
`demo_coding_vs_intergenomic_seqs`.

The only example that is not a table. FASTA carries no label column, so the
file *is* the class.

## Run it

```bash
gb-qc evaluate-classes \
  --input examples/fasta-classes/data/coding_seqs.fasta \
  --input examples/fasta-classes/data/intergenomic_seqs.fasta \
  --out-folder qc-out
```

Each file becomes one class, and its filename stem becomes the label:
`coding_seqs.fasta` → `coding_seqs`. The two labels, ordered alphabetically, give
the report directory `class/sequence/coding_seqs_vs_intergenomic_seqs`.

Two details that surprise people:

- **Classes are always ordered alphabetically by directory name**, so the path is
  the same whichever order you list the files in. Reports are reproducible
  without you having to remember the order you typed.
- **FASTA inputs still land under a `sequence/` directory** even though there is
  no sequence column, so the output layout is identical for every input format.

## What it produces

--8<-- "_generated/fasta-classes-flags.md"

## What you should conclude

**Coding and non-coding DNA differ in composition, strongly.** GC content scores
AU-ROC **0.844**, base composition 0.825, dinucleotide frequencies 0.843. These
are the highest per-sequence scores in any example here.

This is real biology — coding sequence is constrained by the genetic code, which
constrains its composition — and it makes the dataset an easy benchmark. A model
that only counts GC gets 0.84. If you are reporting results on a
coding-vs-noncoding task, that is the number to beat, and plenty of published
architectures beat it by less than you would hope.

**Per-position passes** (0.586) while per-sequence composition fails hard. That
is the informative contrast: the composition difference is spread evenly through
the sequence rather than concentrated at particular positions. Compare
[hidden-motif](hidden-motif.md), which is the exact opposite — clean per-sequence,
failing per-position. The two checks are asking genuinely different questions, and
these two examples are the clearest demonstration that they can disagree.

The flagged features list is long — `A`, `C`, `G`, `T` and ten dinucleotides —
because when overall composition shifts, everything shifts with it. A long list
of flagged features usually means one global difference rather than many separate
ones.

**Nothing else is wrong.** No duplicates, no non-ACGT bases, and every sequence
is exactly 200 nt so `Sequence lengths` is exactly 0.500. This is a
well-constructed dataset that happens to be easy.
