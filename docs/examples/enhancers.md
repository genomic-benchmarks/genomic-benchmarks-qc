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
reporting 0.85 has genuinely learned something. Warnings are most useful read as
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

## The leakage is small and real

<span class="flag flag-warn">Warning</span> at **0.67% of queries and 0.33% of
targets** — 4 test sequences with a 90%-or-better match in training.

This is the best example for reading the leakage report, precisely because the
number is small: the alignment panel lists every leaked pair, and with four of
them you can expand each one and see exactly what is shared, rather than
scrolling through hundreds. Start here before looking at
[composition-bias](composition-bias.md) at 6%.
