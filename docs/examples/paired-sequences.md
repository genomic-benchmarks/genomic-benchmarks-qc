# paired-sequences

**One row, two sequences — and the extra `merged` report.**

miRNA and target-mRNA pairs: 1200 rows from
[miRBench](https://github.com/katarinagresova/miRBench)'s
`AGO2_CLASH_Hejret2023`, each row carrying a `gene` sequence and a
`noncodingRNA` sequence with one shared label.

The only example where a single input file produces **three** class reports.

## Run it

```bash
gb-qc evaluate-classes \
  --input examples/paired-sequences/data/miRNA_mRNA_pairs_dataset.tsv \
  --sequence-column gene \
  --sequence-column noncodingRNA \
  --label-column label \
  --out-folder qc-out
```

Repeat `--sequence-column` per column — `--sequence-column gene noncodingRNA`
will not work. Each column is analysed on its own, and then all sequences are
pooled into an extra `merged` analysis, giving three report directories:

```text
class/
├── gene/0_vs_1/            # the mRNA column alone
├── noncodingrna/0_vs_1/    # the miRNA column alone
└── merged/0_vs_1/          # both pooled
```

Directory names are lowercased and stripped of path-unsafe characters, which is
why `noncodingRNA` becomes `noncodingrna`. The labels shown *inside* the reports
keep their original spelling.

## What it produces

--8<-- "_generated/paired-sequences-flags.md"

## What you should conclude

**Duplication fails; statistics pass.** Both real columns come back
<span class="flag flag-fail">Fail</span> on
`Sequence Duplications within Labels` and
`Duplicate Sequences between Labels`, while every statistical check sits between
0.50 and 0.53 — indistinguishable classes.

That combination is characteristic of interaction datasets, and it is worth
recognising. These are pairs: one miRNA binds many targets, one target is bound
by many miRNAs. So the same `noncodingRNA` sequence appears in many rows by
construction, and so does the same `gene`. The duplication checks are reporting
the shape of the data, not a mistake in it.

It still matters, for two reasons:

- **Your effective sample size is smaller than 1200.** If a handful of miRNAs
  account for most rows, a random train/test split puts the same miRNA on both
  sides, and your model gets credit for recognising it. Splitting by miRNA
  rather than by row is the usual fix.
- **`Duplicate Sequences between Labels` failing** means the same sequence
  appears with both labels — the same miRNA has both bound and unbound targets,
  which is true and expected here. Neither column can be predictive on its own,
  which is the point of a pair dataset.

**Read the `merged` report with care.** It comes back
<span class="flag flag-pass">Pass</span> on all nine, including the duplication
checks the individual columns fail. That is not a contradiction: pooling the two
columns means a duplicated miRNA and a duplicated mRNA are diluted against each
other, and the duplication rate across the pool falls below the threshold. The
`merged` view answers "is there a bias in the sequences taken together", which
is a different question from "is either column duplicated". When they disagree,
the per-column reports are the ones about your data's structure.

The 8 <span class="flag flag-unknown">Unknown</span> sub-checks in the merged and
noncodingrna reports are the positions past where the shorter miRNA sequences
reach — 1200 rows is also close to the 250-per-class floor, so this is a small
dataset by the checks' standards.
