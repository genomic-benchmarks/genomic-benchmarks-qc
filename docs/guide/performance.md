# Runtime and memory

The numbers below are measured on the bundled examples, single-threaded on one
core, so they are a floor rather than a benchmark. They are here because the
scaling is not what people expect.

## Measured

`evaluate-classes`, one core:

| Example | Sequences | Median nt | Window | Time | Peak RSS |
|---|---:|---:|---:|---:|---:|
| [clean-dataset](../examples/clean-dataset.md) | 6,000 | 101 | auto | 7 s | 362 MB |
| [length-bias](../examples/length-bias.md) | 1,000 | 3,462 | 500 | 16 s | 449 MB |
| [hidden-motif](../examples/hidden-motif.md) | 4,000 | 398 | auto | 17 s | 483 MB |
| [variable-length](../examples/variable-length.md) | 7,274 | 331 | auto | 25 s | 385 MB |
| [composition-bias](../examples/composition-bias.md) | 4,904 | 70 | auto | 32 s | 529 MB |

Building all eight examples — fifteen reports, five of them `evaluate-splits`
runs with MMseqs2 — takes about **two and a half to three minutes** in total.

## The surprise: failures cost time

Read that table again. `composition-bias` has the **shortest** sequences of the
five and takes the **longest**. `clean-dataset` has more sequences, longer, and
finishes in a fifth of the time.

The reason is that `composition-bias` fails six checks and `clean-dataset` fails
none. Every flagged check gets an extra annotated figure — `composition-bias`
produces 11 plots against `clean-dataset`'s 6, five of them `_with_flags`
variants — plus a `gb-qc-duplicates.txt` listing the shared sequences.

So **runtime tracks how much is wrong, not just how much data there is.** A first
run on a bad dataset is the slow one; it gets faster as you fix things. Worth
knowing before you conclude the tool is slow on your data.

## What actually drives cost

In rough order of impact:

1. **Positions scored** — the per-position checks are the bulk of the work, and
   they scale with the number of positions in the window times the number of
   sequences reaching each. This is the term to control.
2. **Number of failures** — see above.
3. **Sequence count** — matters, but less than positions.
4. **Number of classes** — every pair of classes is a separate comparison, so
   *k* classes means *k(k−1)/2* reports.

Memory sits in the 350–550 MB range across all five examples, and is dominated by
holding the sequences plus the per-position count matrices. It did not scale
alarmingly with anything here.

## Keeping a long-sequence run sane

If your sequences are thousands of nucleotides, the per-position checks are the
whole cost, and `--end-position` is the lever:

```bash
gb-qc evaluate-classes --input mrnas.csv --end-position 500 --out-folder qc-out
```

[length-bias](../examples/length-bias.md) does this. Its transcripts run to
17,497 nt; unbounded, the per-position checks would run for thousands of
positions — most of them reached by a handful of transcripts — and produce a
report too heavy to open. Trimmed to 500 it takes 16 seconds.

Note what you keep: `--end-position` bounds only the per-position checks. Length,
GC content and composition are computed on whole sequences regardless, so
trimming the window does not weaken the checks that were going to fire anyway.
That is why `length-bias` still reports its length
<span class="flag flag-warn">Warning</span> correctly.

## `evaluate-splits`

Dominated by MMseqs2, not by anything `gb-qc` does. Two options are worth
knowing:

- **`--threads`** is passed straight through. MMseqs2 parallelises well, so this
  is the main lever.
- **`--split-memory-limit 10G`** caps the prefilter structures. Unset there is no
  limit, and a large train set can exhaust the machine. Set it before you find
  out.

For the paper's run over 234 dataset splits, the widest datasets took roughly an
hour per 1,000 nt of window — see [running at scale](at-scale.md).

## Reproducing

```bash
/usr/bin/time -f "%e s  %M KB" gb-qc evaluate-classes \
  --input examples/hidden-motif/data/train.csv --out-folder /tmp/qc
```
