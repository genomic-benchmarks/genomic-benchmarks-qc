# clean-dataset

**The control. What a dataset with nothing wrong with it looks like.**

Mouse transcription-factor binding sites: 3000 sequences per class, every one
101 nucleotides, from [GUE](https://arxiv.org/abs/2306.15006)'s
`transcription_factor_prediction_mouse_1`.

This example exists so the other seven mean something. A tool that flags
everything is as useless as one that flags nothing, and the only way to know
which you have is to point it at a dataset that should come back clean.

## Run it

```bash
gb-qc evaluate-classes \
  --input examples/clean-dataset/data/train.csv \
  --out-folder qc-out
```

```bash
gb-qc evaluate-splits \
  --train-input examples/clean-dataset/data/train.csv \
  --test-input examples/clean-dataset/data/test.csv \
  --sequence-column sequence \
  --out-folder qc-out
```

## What it produces

--8<-- "_generated/clean-dataset-flags.md"

## What you should conclude

Nine checks, nine <span class="flag flag-pass">Pass</span>. Look at the AU-ROC
column rather than the flags, though, because the flags alone undersell it: the
scores run from **0.500 to 0.544**, where 0.500 is exactly chance. The threshold
for a <span class="flag flag-warn">Warning</span> is 0.6. Nothing here is close
to it.

That is what "the classes are not distinguishable by this feature" looks like as
a number. A classifier given only the GC content of these sequences does no
better than a coin. Same for length, base composition, dinucleotide frequencies,
and every individual position.

`Sequence lengths` scores exactly 0.500 because every sequence is 101 nt — there
is nothing there to distinguish, not even by accident. When you see 0.500 on that
row it usually means fixed-length sequences rather than a lucky draw.

**Calibration to take away:** on a clean dataset the checks sit near 0.5, not
near 0.6. If a real dataset of yours comes back with everything in the high
0.5s — technically all <span class="flag flag-pass">Pass</span> — that is not
this. It is a dataset with a faint systematic difference in everything, and
worth a look.

## The one thing that is not clean

The split report flags <span class="flag flag-warn">Warning</span> at **0.07%**
of test queries leaking into training. That is 2 sequences out of 3000.

Two honest notes about that number:

- **It is real, and it is nothing.** Two near-duplicate pairs in 3000 will not
  move a benchmark. This is what the bottom of the leakage scale looks like, and
  it is useful to have seen it before you meet
  [composition-bias](composition-bias.md) at 6%.
- **It would be lower on the full dataset.** These 3000 training sequences are a
  sample of 27000. A test sequence has fewer sequences to match against here
  than it would upstream, so this figure understates the published dataset's
  leakage, which is 0.19%. Still a
  <span class="flag flag-warn">Warning</span>, still negligible.
