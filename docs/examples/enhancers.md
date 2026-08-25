# enhancers

**The quickstart dataset, and a little real train/test leakage.**

1200 training and 600 test sequences sampled from
[Genomic Benchmarks](https://doi.org/10.1186/s12863-023-01123-8)'s
`human_enhancers_ensembl`, running from 4 to 568 nucleotides.

This is the dataset the README's examples use, and the one the threshold
simulations in `analyses/` read. Both quote numbers measured on exactly these
bytes, so it ships unchanged.

## Run it

```bash
gb-qc evaluate-classes \
  --input examples/enhancers/data/enhancers_train.csv \
  --input examples/enhancers/data/enhancers_test.csv \
  --out-folder qc-out
```

```bash
gb-qc evaluate-splits \
  --train-input examples/enhancers/data/enhancers_train.csv \
  --test-input examples/enhancers/data/enhancers_test.csv \
  --sequence-column sequence \
  --out-folder qc-out
```

The two files are **pooled** by the first command and their classes taken from
the label column — passing several `--input` files to `evaluate-classes` does not
compare the files, it combines them. The second command is the one that treats
them as train and test.

## What it produces

--8<-- "_generated/enhancers-flags.md"

## What you should conclude

**Six warnings, no failures — the awkward middle.** GC content 0.662, base
composition 0.634, dinucleotide frequencies 0.688, per-position 0.631. Every one
of them sits between 0.6 and 0.7: a model could get traction, but no single
feature hands it the answer.

This is the most common real-world shape, and the hardest to act on. There is no
one thing to fix. What it tells you is that **a trivial baseline gets roughly
0.66 on this dataset**, so a model reporting 0.70 has added very little, and one
in the mid-0.80s has genuinely learned something. Warnings are most useful read as
a floor for your results, not as a defect list.

`Per sequence dinucleotide content` at 0.688 is the closest to failing, and
dinucleotide frequencies being the strongest signal is a hint about *what* the
difference is — enhancers are CpG-enriched relative to background, and that shows
up in pairs before it shows up in single-base composition.

**Variable lengths, quietly.** These sequences run 4 to 568 nt, so the
per-position checks reach position 483 but only compare positions 1–359, and
1–359 is what the figures draw. Less dramatic than
[variable-length](variable-length.md), but the same mechanism — worth noticing
that a report can be silently windowed without anything failing.

## The only two flagged positions are the two ends

Everything positional about this dataset is at its edges. The per-position checks
flag exactly one position each, and it is the same position in both:

--8<-- "_generated/enhancers-positions.md"

Forward position 1 is the first base of each sequence; reversed position 1 is the
last. In class 1 the first base is `T` in 56% of sequences and `A` in 4%, and the
last base is `A` in 51% and `T` in 5%. Class 0 is spread across all four bases at
both. **Nothing in between is flagged at all.**

That is the shape [the checks page](../guide/checks.md#per-position-nucleotide-content)
calls a technical artefact rather than biology: whatever produced these intervals
put a systematic base at each end of one class, and `T` at the start against `A`
at the end is one bias seen from two directions. It is not an artefact of the
subsample either — run the check on the full upstream `human_enhancers_ensembl`
splits and both ends flag there too, at the same strength, in each split
independently.

**This is also the clearest case for the reversed panel.** These sequences run 4
to 568 nucleotides, so "the last base" sits at a different forward index in every
one of them. The forward check cannot see an end-of-sequence bias at all; only
the reversed check can. On [hidden-motif](hidden-motif.md), where every sequence
is the same length, the reversed panel is the forward one mirrored and adds
nothing. Here it adds a finding.

!!! note "What a subsample can hide"

    The full upstream splits also fail `Unique bases`, `Sequence Duplications
    within Labels` and `Duplicate Sequences between Labels`. The 1,800 sequences
    here happen not to contain the sequences responsible, so those three come
    back clean. It is a general point worth carrying: **the three categorical
    checks can be hidden by sampling**, because one non-ACGT character or one
    shared sequence is enough to fail them and half the data may not include it.
    The AU-ROC checks are far more stable — every one of them here lands within
    a hundredth of the full dataset's.

## The leakage is small and real

<span class="flag flag-warn">Warning</span> at **0.67% of queries and 0.33% of
targets** — 4 test sequences with a 90%-or-better match in training.

This is the best example for reading the leakage report, precisely because the
number is small: the alignment panel lists every leaked pair, and with four of
them you can expand each one and see exactly what is shared, rather than
scrolling through hundreds. Start here before looking at
[composition-bias](composition-bias.md) at 6%.
