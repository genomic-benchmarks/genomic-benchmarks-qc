# Troubleshooting

## `mmseqs: command not found`

`evaluate-splits` needs [MMseqs2](https://github.com/soedinglab/MMseqs2) on your
`PATH`; it is not a Python dependency and `pip` will not install it. Take a
precompiled binary — see
[installing MMseqs2](installation.md#mmseqs2-for-the-leakage-check).

`evaluate-classes` does not need it, so if you only want the class checks you can
skip this entirely.

## Everything came back `Unknown`

Your classes are smaller than 250 sequences each, which is the floor below which
nothing is scored. This is not a failure — it is the tool declining to make a
claim it cannot support. See
[how a flag is decided](guide/how-it-works.md#at-least-250-sequences-per-class).

**It does not mean the dataset is clean.** The plots, per-class statistics and
descriptive tables are still computed from all of your data, so open the report
and compare by eye. That is often enough to see an obvious problem.

The floor is not adjustable. `--min-coverage 0` removes the *fraction* rule for
per-position checks but leaves the 250 in place.

## The per-position plot stops before my sequences end

Expected on variable-length data. Positions are only flagged where at least 25% of
each class still reaches them, and the figures draw the flagged window. Everything
past it is <span class="flag flag-unknown">Unknown</span>.

`--min-coverage 0` extends it as far as the 250-sequence floor allows. Read
[the window section](guide/per-position.md#what-the-window-means) first — the
reason for the rule is not sample size, and widening it changes what the numbers
mean.

## It is slow

Drawing the figures, and computing the per-sequence statistics. On 10,000
sequences of 500 bp the whole run is about seven seconds, and roughly half of
that is the figures.

So the lever is `--report-types simple`, which writes every flag and metric as a
CSV and draws nothing:

```bash
gb-qc evaluate-classes --input long.csv --report-types simple --out-folder qc-out
```

That is 7.0 s down to 3.6 s on the same data, and 450 MB down to 170 MB.

`--end-position` is no longer a speed lever — on that run it saves about a tenth
of a second. The per-position checks used to be most of the time; they are now
scored from per-position base counts rather than from one array per sequence, and
cost milliseconds however wide the window is. It is still how you say the far end
of your sequences is not worth reporting on.

Also worth knowing: **failures cost runtime**, because each flagged check gets an
extra annotated figure. A first run on a messy dataset is the slow one.

## MMseqs2 runs out of memory

Cap the prefilter structures:

```bash
gb-qc evaluate-splits ... --split-memory-limit 10G
```

Unset, there is no limit. On a large training set this is the thing that exhausts
the machine.

## `--input a.csv b.csv` does not work

Repeat the option instead:

```bash
gb-qc evaluate-classes --input a.csv --input b.csv --out-folder qc-out
```

Same for `--sequence-column` and `--label-list`. This catches everyone once.

## Several `--input` files did not become several classes

For CSV/TSV they are **pooled** and classes come from the label column. One file
or ten, the classes are whatever is in `--label-column`.

FASTA is the opposite: it carries no labels, so each file *is* a class and its
filename stem becomes the label. See
[fasta-classes](examples/fasta-classes.md).

## It refused to infer my classes

`--label-list infer` stops at 50 distinct values. Every pair of classes gets its
own comparison, so the work is quadratic — 600 values is 179,700 reports — and a
column with that many values is nearly always a continuous target or the sequence
column, not a label column.

Three ways on, and the error names all three: `--regression` splits a numeric
target at its median, `--label-column` points at a different column, and
`--label-list` names the classes to compare. An explicit list is not capped.

## My class labels are not the directory names I expected

Directory names are lowercased and stripped of characters unsafe on some
filesystems, then ordered alphabetically — so the same dataset gives the same
paths whatever order you passed the files in. `noncodingRNA` becomes
`noncodingrna`. Labels shown *inside* the reports keep their original spelling.

If two names would collide, `gb-qc` makes them unique and warns which name it
used. Worth grepping the log for at scale.

## `--regression` exited with an error

The label column is split at its median into `high` and `low`. If that does not
produce two non-empty classes — a constant column, or one that is mostly a single
value — there is nothing to compare and the run stops. Non-numeric rows are
dropped with a warning first, so check how many survived.

## A check fails but I think it is the biology

Quite likely, and the tool cannot tell the difference. A splice-site dataset
*should* fail its per-position check at the splice site — see
[hidden-motif](examples/hidden-motif.md).

What the flag tells you is that a trivial model does well on that feature. Whether
that is signal or leakage is your call. Either way the useful move is the same:
treat the flagged feature's AU-ROC as the floor your model has to beat, and report
it.

## The report is enormous

The per-position payload is most of a report's bytes, and it scales with the
window. `--end-position` is the lever. Reports on the bundled examples run
850 KB to 1.1 MB, which is normal for a standalone file with its plots inlined.

## Can I get the raw numbers?

Add `json` to `--report-types` for per-class statistics — counts, GC, lengths,
base and dinucleotide frequencies — at
`<out-folder>/class/<column>/per-class/<class>.json`.

The `gb-qc-report.csv` beside every report carries every flag and score,
including the per-base and per-position breakdown. For leakage,
`mmseqs/mmseqs2_search_result.tsv` has every pair at or above the similarity
threshold. See
[using it in CI](guide/ci.md) for the CSV's shape.
