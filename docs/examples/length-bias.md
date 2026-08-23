# length-bias

**When length alone separates the classes — and the only example with a
continuous label.**

Human mRNA translation efficiency: 1000 transcripts from
[mRNABench](https://github.com/morrislab/mRNAbench), running from 425 to 17,497
nucleotides with a median of 3463. The label is a continuous value, not a class.

This example carries the README's first claim — *your negatives are shorter than
your positives, and a length classifier now beats your model* — and two options
nothing else exercises.

## Run it

```bash
gb-qc evaluate-classes \
  --input examples/length-bias/data/train.csv \
  --label-column target \
  --regression \
  --end-position 500 \
  --out-folder qc-out
```

Two options here are doing real work:

- **`--regression`** splits the continuous `target` column at its median into
  `high` and `low` classes. That is why the report directory is
  `class/sequence/high_vs_low` rather than `0_vs_1`.
- **`--end-position 500`** stops the per-position checks at position 500. These
  are whole mRNAs; left unbounded the checks would run for thousands of positions
  and produce a report too large to open, for positions only the longest few
  transcripts reach. This is what that option is for.

## What it produces

--8<-- "_generated/length-bias-flags.md"

## What you should conclude

**Length carries the label.** `Sequence lengths` scores AU-ROC **0.631** — a
classifier told nothing but how long a transcript is predicts which side of the
median its translation efficiency falls on, 63% of the time.

For this dataset that is biology: transcript length genuinely relates to
translation efficiency. But it is also the single easiest thing for a model to
learn, and a sequence model that scores 0.65 here may have learned nothing but
length. If you are benchmarking on this, that is the number to beat, and it is
worth reporting alongside your model's.

**The composition warnings are downstream of it.** Base composition 0.617 and
dinucleotide frequencies 0.618 are not separate findings: longer transcripts have
systematically different composition (GC content varies with length in mRNAs), so
a length difference shows up as a composition difference. The flagged features
are `T`, `TT`, `CA` and `GT` — consistent with a length-linked shift rather than
a specific motif.

**Per-position is clean** (0.552) even though composition is not. That fits: the
difference is in overall composition, not in any particular position, which is
exactly the pattern to expect from a length-driven bias.

## About that Warning

`Sequence lengths` is a <span class="flag flag-warn">Warning</span>, not a
<span class="flag flag-fail">Fail</span>, and this is the strongest length bias
available rather than a weak example. Across the 234 dataset splits surveyed for
the paper, **no dataset fails the length check outright.** The strongest on
record scores 0.61.

So the README leads with length bias because it is the easiest bias to introduce
accidentally and the most embarrassing to ship — not because published benchmarks
are full of it. Published benchmarks are mostly length-matched, often
deliberately. Yours might not be, and it is the cheapest thing to check.
