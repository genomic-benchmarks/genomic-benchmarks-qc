# variable-length

**Why most positions in a report go unscored.**

Tomato long non-coding RNAs: 7274 training sequences from
[OmniGenBench](https://doi.org/10.48550/arXiv.2505.14402)'s
`lncrna_s_lycopersicum`, running from **105 to 5964 nucleotides** with a median
of 331. Shipped whole, so the flags are the published dataset's.

Every other example has fixed-length sequences, where every sequence reaches
every position and the per-position checks are simple. This one is what happens
when they do not.

## Run it

```bash
gb-qc evaluate-classes \
  --input examples/variable-length/data/train.csv \
  --out-folder qc-out
```

```bash
gb-qc evaluate-splits \
  --train-input examples/variable-length/data/train.csv \
  --test-input examples/variable-length/data/test.csv \
  --sequence-column sequence \
  --out-folder qc-out
```

## What it produces

--8<-- "_generated/variable-length-flags.md"

## What you should conclude

### Most of this report is `Unknown`, and that is correct

Of the 26,408 per-position sub-checks in the class report, **20,870 are
<span class="flag flag-unknown">Unknown</span>** — 79% of them. Three different
boundaries are in play, and seeing them in one report is the point of this
example:

| Position | What happens there |
|---|---|
| **541** | As far as 25% of each class still reaches. This is the default `--min-coverage` window: the furthest a position can be *flagged*, and the last position the figures draw |
| **2636** | As far as 250 sequences per class reach. The furthest anything is *scored* at all |
| **5964** | The longest sequence in the dataset |

Past 541, positions are reported <span class="flag flag-unknown">Unknown</span>.
Not <span class="flag flag-pass">Pass</span> — the comparison was not made.

The reason is not sample size. It is that a cohort far along the sequence is not
a sample of the class; it is a sample of *the class's longest sequences*. If
length correlates with composition — and in lncRNAs it does — then a difference
found at position 3000 is a difference between long sequences, which is a
different claim from a difference between classes. No amount of data fixes that.
Declining to score is the only honest answer.

**If you need to look further out**, `--min-coverage 0` drops the fraction rule
and leaves only the 250-sequence floor, which extends flagging to position 2636.
Read what comes back with the caveat above firmly in mind.

### The vocabulary check fails

<span class="flag flag-fail">Fail</span> on `Unique bases`: some sequences
contain characters outside ACGT. In a genomic dataset that is usually `N` from
unresolved assembly, and it matters for a mundane reason — your tokeniser has to
do something with them, and whatever it does is a modelling decision you did not
consciously make. Worth finding before training rather than after.

### Sequences repeat, within and between classes

<span class="flag flag-fail">Fail</span> on
`Duplicate Sequences between Labels` and
<span class="flag flag-warn">Warning</span> within labels. Identical sequences
carrying opposite labels are unlearnable by construction and cap achievable
accuracy.

### Composition is clean

Everything statistical passes: GC 0.508, composition 0.556, dinucleotides 0.558,
per-position 0.559. So this is a dataset whose *content* looks fine and whose
*bookkeeping* does not — duplicates and non-ACGT characters. A useful shape to
recognise, because it is fixable without touching the biology.

## Leakage

<span class="flag flag-warn">Warning</span> at 0.80% of queries, 0.12% of
targets. Small, real, worth a look at the alignment panel to see what the shared
sequence actually is.
