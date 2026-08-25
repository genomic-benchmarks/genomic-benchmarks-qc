# Train/test leakage

`evaluate-classes` asks whether your two classes are distinguishable.
`evaluate-splits` asks a different question entirely:

> **Is your test set already in your training set?**

A model that has seen a test sequence during training scores on it for free. If
enough of the test set leaks, a benchmark stops measuring generalisation and
starts measuring memorisation — and it will rank a memorising model above a
generalising one.

```bash
gb-qc evaluate-splits \
  --train-input examples/enhancers/data/enhancers_train.csv \
  --test-input examples/enhancers/data/enhancers_test.csv \
  --sequence-column sequence \
  --out-folder qc-out
```

Unlike `evaluate-classes`, this one needs
[MMseqs2](https://github.com/soedinglab/MMseqs2) on your `PATH` —
[how to install it](../installation.md#mmseqs2-for-the-leakage-check).

## How similarity is measured

Every test sequence is searched against the whole training set with
`mmseqs easy-search`, restricted to nucleotide search on the forward strand. Each
hit gets a similarity score:

```text
similarity = min(query_coverage, target_coverage) × percent_identity
```

Both halves matter, and the `min` is the important part:

- **Percent identity alone is not enough.** A 20-nucleotide stretch matching
  perfectly inside a 500-nucleotide sequence is 100% identical over the aligned
  region and tells you nothing about leakage.
- **Coverage of the *shorter* sequence** is what makes the score symmetric and
  honest. A short sequence fully contained in a long one is leakage — the model
  has seen all of it — and taking the minimum coverage catches that, where
  averaging the two would dilute it away.

A test sequence counts as leaked when its **best** hit exceeds
`--similarity-threshold`, 90% by default. So the reported percentage is a count of
test sequences with at least one near-duplicate in training, not a count of hits.

### Queries and targets are not the same number

The report gives two percentages, and they answer different questions:

- **Leaked queries** — what fraction of your *test* set has a near-duplicate in
  training. **This is the number that matters for your benchmark**, because it
  bounds how much of your test score could come from memorisation.
- **Leaked targets** — what fraction of your *training* set is implicated. Almost
  always the smaller number, because one popular training sequence can be the
  match for many test sequences.

In [variable-length](../examples/variable-length.md) the split is 0.80% of
queries against 0.12% of targets: a small number of training sequences accounting
for most of the leakage. That shape is worth noticing — it usually means a
repeated family rather than diffuse overlap.

## Reading the report

Three things, in the order they are useful:

1. **The two percentages and the flag.** Start here and decide whether to keep
   reading.
2. **The similarity histogram.** The shape tells you what kind of leakage you
   have. A spike at 100% is exact duplicates — a plain mistake, easily fixed. A
   broad hump in the 90s is a repeated sequence family, which is harder, because
   removing it may remove a real biological class.
3. **The leaked-pair panel.** Up to the first 100 pairs, each expanding to its
   rendered alignment, so you can see what is actually shared. Every hit,
   including those past 100, is exported to
   `mmseqs/mmseqs2_search_result.tsv` beside the report.

**Start with a small example.** [enhancers](../examples/enhancers.md) has four
leaked pairs — you can expand every one and understand the whole finding in a
minute. [composition-bias](../examples/composition-bias.md) at 6.04% has far too
many to read individually, which is the point at which the histogram becomes more
useful than the list.

## What the numbers look like in practice

Measured on the bundled examples, which between them span the range:

| Example | Leaked queries | Flag | What it is |
|---|---|---|---|
| [hidden-motif](../examples/hidden-motif.md) | 0.00% | <span class="flag flag-pass">Pass</span> | Genuinely clean — worth seeing once, so you know the check can come back empty |
| [clean-dataset](../examples/clean-dataset.md) | 0.07% | <span class="flag flag-warn">Warning</span> | 2 sequences in 3,000. Real, and negligible |
| [enhancers](../examples/enhancers.md) | 0.67% | <span class="flag flag-warn">Warning</span> | 4 pairs. The right size to learn the report on |
| [variable-length](../examples/variable-length.md) | 0.80% | <span class="flag flag-warn">Warning</span> | Concentrated in few training sequences |
| [composition-bias](../examples/composition-bias.md) | 6.04% | <span class="flag flag-fail">Fail</span> | One test sequence in sixteen. This one changes conclusions |

Two things to take from that table. First, **a small non-zero number is normal** —
real genomic data contains repeated families, and 0.1% is not worth acting on.
Second, the interesting threshold is not the flag boundary but whether the leaked
fraction is large relative to the *differences between models you are comparing*.
If two architectures are half a point apart and 6% of the test set is leaked, the
comparison is noise.

## What to do about it

**Exact duplicates**: remove them from the test set. Uncontroversial.

**Near-duplicates**: the fix is to split on something other than sequence
identity. Cluster the sequences first — MMseqs2 will do it — and assign whole
clusters to train or test. This is standard practice in protein benchmarking and
underused in genomics.

**Positional overlap**: if train and test were drawn from overlapping genomic
windows, split by chromosome or by locus rather than by window. This is the most
common cause of the broad-hump histogram.

**Nothing**: a legitimate answer when the leaked fraction is small and you say so
in the paper. The failure mode to avoid is not knowing.

## Practical notes

- **Multiple sequence columns** are concatenated and searched together, unlike
  `evaluate-classes` which analyses each column separately. There is no
  per-column split report.
- **Memory.** MMseqs2's prefilter structures can be large on big datasets.
  `--split-memory-limit 10G` caps them; without it there is no limit and a big
  run can exhaust the machine.
- **Threads.** `--threads` is passed straight through; the default is whatever
  MMseqs2 picks.
- **Temporary files** go into a `gb-qc-mmseqs-*/` directory inside the comparison
  directory and are removed at the end. `--keep-tmp-files` keeps them and logs
  each path. Each run gets its own, so several runs can share one `--out-folder`
  concurrently.
