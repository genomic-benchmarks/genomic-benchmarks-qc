# composition-bias

**The worst case, and the only example that fails both kinds of report.**

Human core promoters with a TATA box: 4904 training sequences, every one 70
nucleotides, from [GUE](https://arxiv.org/abs/2306.15006)'s
`core_promoter_detection_tata`. Small enough that both splits ship whole, so
every number here is the published dataset's own, not a subsample's.

Six of nine checks fail. The classes differ in GC content, in base composition,
in dinucleotide frequencies and at individual positions; they share sequences
with each other; and 6% of the test set is already in the training set.

## Run it

```bash
gb-qc evaluate-classes \
  --input examples/composition-bias/data/train.csv \
  --out-folder qc-out
```

```bash
gb-qc evaluate-splits \
  --train-input examples/composition-bias/data/train.csv \
  --test-input examples/composition-bias/data/test.csv \
  --sequence-column sequence \
  --out-folder qc-out
```

## What it produces

--8<-- "_generated/composition-bias-flags.md"

## What you should conclude

Read the failures as one finding, not six. GC content 0.701, base composition
0.717, dinucleotide frequencies 0.716, per-position 0.710 — four checks landing
within 0.02 of each other is not four independent biases. It is **one
compositional difference between the classes, seen four ways.**

The classes are TATA-box promoters and non-promoters. A TATA box is AT-rich by
definition, so the positive class is AT-rich, so GC content separates the
classes, so base composition does, so dinucleotide frequencies do, so the
positions carrying the motif do. Each check is measuring a consequence of the
same thing.

What that means for a model: **a classifier that only counts bases gets AU-ROC
0.70 on this dataset.** Whatever your model scores, subtract that as the floor.
A reported 0.75 is barely above counting nucleotides.

Two further failures are separate from the composition story and worth their own
attention:

- **`Duplicate Sequences between Labels` <span class="flag flag-fail">Fail</span>.**
  The same sequence appears in both classes, which means part of your training
  signal is contradictory — identical input, opposite label. No model can learn
  that, and it caps achievable accuracy.
- **Leakage <span class="flag flag-fail">Fail</span> at 6.04% of queries.** One
  test sequence in sixteen has a 90%-or-better match in training. A model that
  memorises will score on those for free.

The <span class="flag flag-warn">Warning</span> on
`Sequence Duplications within Labels` says sequences repeat inside a class too —
less damaging than across classes, but it means the effective dataset is smaller
than the row count suggests.

## Why both reports matter

This is the one example where `evaluate-classes` and `evaluate-splits` both come
back red, and the two findings compound rather than repeat. The class checks say
a trivial model does well; the leakage check says a memorising model also does
well. A benchmark with both is not measuring what it claims to.
