# How a flag is decided

Every check `gb-qc` runs asks the same question, and it is a deliberately narrow
one:

> **Could a classifier that sees only this one feature tell the two classes
> apart?**

Not "are the classes different" — they had better be, or there is nothing to
model. The question is whether they differ in a way that a model could exploit
*without understanding anything about the biology*. A benchmark where GC content
alone reaches AU-ROC 0.85 is not measuring what its authors think it is.

## The measurement

For a feature — say GC content — the procedure is:

1. Compute the feature for every sequence in both classes.
2. Label class A positive and class B negative.
3. Feed the feature values in as if they were a classifier's scores, and take
   the **AU-ROC**.

There is no model to train and nothing to tune. The feature *is* the score, so
the AU-ROC is a closed-form statement about how separable the classes are along
that one axis. 0.5 means the feature is useless; 1.0 means it settles the
question by itself.

One detail worth knowing: **direction is discarded**. If the AU-ROC comes out
below 0.5 the scores are inverted and it is reported as `1 − auroc`. "Class A has
more GC" and "class A has less GC" are the same finding — a difference either way
is exploitable — so the reported value is always ≥ 0.5.

Alongside AU-ROC the report gives AU-PR (average precision) and the accuracy at
the best threshold.
They are there for context; **the flag comes from AU-ROC alone.**

### It measures difference, not structure

Because the comparison is always between the two classes, a feature both classes
share scores 0.5 no matter how striking it is. The clearest case is
[hidden-motif](../examples/hidden-motif.md#the-one-position-in-the-motif-that-passes),
where one position holds the same base in 100% of every sequence in both classes
and the check reports it
<span class="flag flag-pass">Pass</span>. That is not the check missing
something. A position that never varies cannot tell the classes apart, so it
cannot be exploited, so there is nothing to flag.

The practical consequence is that **the flagged region is usually narrower than
the pattern causing it.** Flags mark where a difference is large enough to
measure, not where a motif starts and stops.

## The boundaries

| AU-ROC | Flag | Read it as |
|---|---|---|
| ≤ 0.6 | <span class="flag flag-pass">Pass</span> | Not distinguishable by this feature |
| 0.6 – 0.7 | <span class="flag flag-warn">Warning</span> | A model could get some traction here |
| > 0.7 | <span class="flag flag-fail">Fail</span> | Significant bias |
| — | <span class="flag flag-unknown">Unknown</span> | Not scored; see below |

The 0.6 boundary is not a round number picked for looking sensible. It comes out
of simulation against datasets with no real difference between the classes: below
250 sequences per class, noise alone pushes the worst per-sequence check over 0.6
in a substantial fraction of replicates, and at 250 it essentially stops.

What the boundaries do **not** encode is whether a finding matters for *your*
task. A <span class="flag flag-fail">Fail</span> on a splice-site dataset's
per-position check is the biology — see
[hidden-motif](../examples/hidden-motif.md). The flag says "a trivial model does
well here"; only you know whether that is the signal or the leak.

## Worst case, not average

Several headline checks are made of many sub-checks. `Per sequence nucleotide
content` scores A, C, G and T separately. `Per sequence dinucleotide content`
scores all sixteen pairs. `Per position nucleotide content` scores every base at
every position — for a 398-nucleotide dataset that is over 1,500 comparisons.

**The headline number is the maximum across them, not the mean.**

This matters for reading a report. When
[hidden-motif](../examples/hidden-motif.md) reports `Per position nucleotide
content` at 0.783, that is *one* base at *one* position — G at position 200. The
other 1,500-odd comparisons are unremarkable. The headline is telling you the
worst thing in there, which is the right summary for a QC tool: a single
give-away position is a real problem even when the average position is fine, and
an average would hide it.

The consequence is that a failing headline check is an invitation to open the
report and find *which* sub-check failed. The per-feature and per-position tables
under each section are where the actual finding is.

## Two guards before anything is scored

An AU-ROC only means what it says when there is enough behind it, and when what
is behind it is the class rather than a corner of it. So two rules sit in front
of the scoring.

### At least 250 sequences per class

Nothing is scored on fewer. For the per-sequence checks this is a floor on the
smaller class. For the per-position checks it applies per position — a position
is only compared using sequences long enough to have it, so the cohort shrinks as
you move along the sequence.

Below that floor the boundary stops meaning anything: the simulations show the
worst per-sequence check crossing 0.6 on pure noise in 19.4% of replicates at 100
sequences per class, against 0.2% at 250.

### At least 25% of the class must reach a position

This one is not about sample size, and it catches people out.

Consider a dataset of variable-length sequences and a position far along them.
The cohort reaching that position may be large — thousands of sequences, well
past the 250 floor. But it is not a sample of the class. It is a sample of *the
class's longest sequences*. If length correlates with composition, and in real
sequence data it usually does, then a difference measured there is a difference
between long sequences, which is a different claim from a difference between
classes.

No amount of data fixes that. A bigger cohort of long sequences is still only
long sequences. So the per-position checks stop flagging where the cohort falls
below `--min-coverage` (25% by default) of its class, and everything past that is
<span class="flag flag-unknown">Unknown</span>.

[variable-length](../examples/variable-length.md) is the example where all of
this is visible at once: flagging stops at position 549, the checks report as far
as 2636, and the longest sequence runs to 5964.

## `Unknown` is not `Pass`

This is the single most important thing to get right when reading a report.

<span class="flag flag-pass">Pass</span> means: *the comparison was made, and the
classes were not separable.*

<span class="flag flag-unknown">Unknown</span> means: *the comparison was not
made.*

No evidence of a difference is not evidence of no difference. A small dataset
where everything comes back
<span class="flag flag-unknown">Unknown</span> has not been given a clean bill of
health — it has not been checked. The terminal output says which checks were
skipped and why, so read it rather than assuming.

The plots, the per-class statistics and the descriptive tables are computed from
**all** the data regardless of whether a check was scored. So a dataset too small
to flag can still be compared by eye, which is often enough to see an obvious
problem.

## What is not decided by AU-ROC

Three checks are not scores and have no AU-ROC:

- **`Unique bases`** — do the two classes use the same alphabet? Yes or no.
- **`Sequence Duplications within Labels`** — a threshold on how much of the data
  survives deduplication.
- **`Duplicate Sequences between Labels`** — is any sequence in both classes at
  all?

Their exact rules are on the [checks page](checks.md).
