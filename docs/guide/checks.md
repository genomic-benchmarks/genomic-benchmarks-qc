# The checks, and what to do about them

Nine checks. For each one: what it measures, the exact rule behind the flag, what
it means when it fires, and what to actually do.

For how a score becomes a flag, see [how a flag is decided](how-it-works.md).
Every check below links to an example where you can see it fire on real data.

!!! tip "Read a failure as a floor, not a defect"

    Most of these checks are telling you what a trivial model scores on your
    dataset. If GC content alone reaches AU-ROC 0.70, then 0.70 is the number
    your model has to beat to have demonstrated anything. That is more useful
    than treating the flag as a pass/fail gate — and it is why a
    <span class="flag flag-warn">Warning</span> is often more actionable than it
    looks.

## Unique bases

**Measures.** Whether the two classes use the same set of characters.

**Rule.** <span class="flag flag-fail">Fail</span> if the sets differ at all,
<span class="flag flag-pass">Pass</span> if identical. There is no
<span class="flag flag-warn">Warning</span> — it is a yes-or-no question.

Note what this is *not*: it does not flag non-ACGT characters as such. A dataset
where both classes contain `N` passes. What fails is an **asymmetry** — one class
containing a character the other does not.

**When it fires.** That asymmetry is a give-away of the worst kind. If `N`
appears only in the positive class, then every sequence containing `N` is
perfectly classifiable, and a model will find that long before it finds any
biology. It does not take many: in
[variable-length](../examples/variable-length.md) just **6 sequences out of
3,607** carry an `N` that the other class never has, and the check fails on
that.

**What to do.** Find the odd characters and decide deliberately:

- If they are `N` from unresolved assembly, either drop those sequences from both
  classes or ensure both classes are sampled the same way. Asymmetry usually
  means the two classes came through different pipelines.
- If they are IUPAC ambiguity codes, decide what your tokeniser does with them —
  silently mapping to `N`, or to a random base, is a modelling decision you
  should make on purpose.
- If they are lowercase (soft-masked repeats), normalise case before anything
  else.

**See it fail:** [variable-length](../examples/variable-length.md).

## Sequence lengths

**Measures.** AU-ROC of sequence length as a classifier.

**When it fires.** A model that counts characters beats your model. This is the
easiest bias to introduce by accident — negatives sampled from a different source
than positives, or trimmed differently — and among the most embarrassing to
publish.

**What to do.** Match the length distributions. If negatives are drawn from
background genome, draw them at the length of the positive they pair with rather
than at a fixed length. If the difference is intrinsic to the biology (mRNAs
really do differ in length by translation efficiency), report the length-only
baseline alongside your model so a reader can see what the model added.

**See it warn:** [length-bias](../examples/length-bias.md).

## Per sequence GC content

**Measures.** AU-ROC of each sequence's GC fraction.

**When it fires.** The classic compositional confound. It fires whenever the two
classes come from different genomic contexts — promoters versus background,
coding versus intergenic, CpG islands versus not.

**What to do.** GC-match your negatives. Sampling background regions at the same
GC content as your positives is standard practice and removes most of this. If
you cannot, expect a reviewer to ask what a GC-only baseline scores, and have the
answer ready.

**See it fail:** [fasta-classes](../examples/fasta-classes.md) at 0.844,
[composition-bias](../examples/composition-bias.md) at 0.701.

## Per sequence nucleotide content

**Measures.** The fraction of A, C, G and T per sequence, scored separately. The
headline is the **worst** of the four.

**When it fires.** Almost always alongside GC content, because they measure
overlapping things. When both fire at similar values, treat it as one finding.

**What to do.** As for GC content. If this fires *without* GC content firing, the
imbalance is between A and T or between C and G rather than between GC and AT,
which points at strand asymmetry — worth checking whether one class was
reverse-complemented and the other not.

**See it fail:** [composition-bias](../examples/composition-bias.md),
[fasta-classes](../examples/fasta-classes.md).

## Per sequence dinucleotide content

**Measures.** All sixteen dinucleotide frequencies per sequence, scored
separately. Headline is the worst of the sixteen.

**When it fires.** More sensitive than single-base composition, and the flagged
dinucleotide tells you something specific. `CG` is the one to look for: CpG
depletion varies hugely across genomic contexts, so a `CG`-driven failure often
means the classes differ in methylation context or in promoter proximity.

**What to do.** Check *which* dinucleotide. A single flagged pair points at a
specific motif or context; all sixteen flagged means overall composition has
shifted and this is a shadow of the GC finding.

**See it:** [enhancers](../examples/enhancers.md) at 0.688 (CpG in enhancers),
[hidden-motif](../examples/hidden-motif.md) at 0.602 (a fixed-position GT motif
leaking into whole-sequence frequency).

## Per position nucleotide content

**Measures.** Each base at each position, scored separately — for a 400-nt
dataset, over 1,500 comparisons. Headline is the worst single one.

**When it fires.** Something at a *specific location* gives the class away.
Common causes, roughly in order of how often they turn out to be the real
explanation:

- An adapter, primer or barcode left on one class.
- A padding or trimming convention applied to one class only.
- A fixed-offset biological motif — a TATA box, a splice site — which is the
  signal, not a bug.
- One class starting or ending with `N`.

**What to do.** Open the report and look at *where*. This is the check where the
headline number tells you almost nothing and the per-position detail tells you
everything — the [interactive plot](per-position.md) exists for exactly this. Six
flagged positions clustered together is a motif; position 1 alone is a technical
artefact; a broad smear is composition showing through.

One thing this check does **not** do is reward conservation. It compares the two
classes, so a position both classes share is a
<span class="flag flag-pass">Pass</span> however invariant it is — see the
[position every sequence agrees on](../examples/hidden-motif.md#the-one-position-in-the-motif-that-passes),
where 100% of both classes carry the same base and the check has nothing to say.
Expect the flagged window to be narrower than the motif that causes it.

**See it fail:** [hidden-motif](../examples/hidden-motif.md) at 0.783, from six
positions out of 398. **See the artefact shape:**
[enhancers](../examples/enhancers.md#the-only-two-flagged-positions-are-the-two-ends),
where the only flagged position is the first base.

## Per position reversed nucleotide content

**Measures.** The same thing, counted from the end of each sequence.

**Why it exists.** For fixed-length sequences it is redundant — it finds exactly
what the forward version finds, which is why the two usually report identical
numbers. It earns its place on **variable-length** data, where a feature at a
fixed distance from the *end* (a poly-A tail, a 3' adapter) sits at a different
forward position in every sequence and the forward check cannot see it.

**What to do.** If the reversed check fires and the forward one does not, look for
something anchored to the end of your sequences.

**See it earn its place:**
[enhancers](../examples/enhancers.md#the-only-two-flagged-positions-are-the-two-ends).
Both directions flag position 1 there, but they are not the same position — the
forward check finds the first base of each sequence and the reversed check the
last, and on sequences running from 4 to 568 nucleotides only the reversed check
could have found the second one.

## Sequence Duplications within Labels

**Measures.** How much of the data survives deduplication, pooled across both
classes.

**Rule.** <span class="flag flag-fail">Fail</span> below 98% remaining,
<span class="flag flag-warn">Warning</span> between 98% and 100%,
<span class="flag flag-pass">Pass</span> at exactly 100%. Note that any
duplication at all is at least a <span class="flag flag-warn">Warning</span>.

**When it fires.** Your effective dataset is smaller than your row count. Worse,
if duplicates land on both sides of a random train/test split, your test set is
partly a copy of your training set and your reported accuracy is inflated.

**What to do.** Deduplicate, or split on a grouping key rather than on rows. For
interaction datasets this is expected rather than wrong — see
[paired-sequences](../examples/paired-sequences.md), where one miRNA legitimately
appears in many rows — but it still means splitting by row is the wrong move.

**See it fail:** [paired-sequences](../examples/paired-sequences.md).

## Duplicate Sequences between Labels

**Measures.** Whether any sequence appears in both classes.

**Rule.** <span class="flag flag-fail">Fail</span> if the intersection is
non-empty. One shared sequence is enough. There is no
<span class="flag flag-warn">Warning</span>.

**When it fires.** Part of your training signal is self-contradictory: identical
input, opposite label. No model can learn it, and it puts a hard ceiling on
achievable accuracy — a ceiling you may spend weeks trying to beat.

**What to do.** The offending sequences are written to
`gb-qc-duplicates.txt` beside the report, so you can look at them. Then decide
whether they are a labelling error (remove or relabel) or genuine ambiguity
(a sequence that really is both, which means the task is not well posed at
sequence level).

**See it fail:** [composition-bias](../examples/composition-bias.md),
[paired-sequences](../examples/paired-sequences.md),
[variable-length](../examples/variable-length.md).

## Data leakage

Not part of `evaluate-classes` — this is what `evaluate-splits` does, and it has
[a page of its own](leakage.md).

## When several checks fail together

Failures are usually correlated, and reading them as one finding rather than
several is most of the skill in using this tool.

| Pattern | Usually means |
|---|---|
| GC, composition and dinucleotides all fail at similar values | **One** compositional difference, seen three ways. See [composition-bias](../examples/composition-bias.md) |
| Per-sequence fails, per-position passes | The difference is spread through the sequence. [fasta-classes](../examples/fasta-classes.md) |
| Per-position fails, per-sequence passes | The difference is localised. [hidden-motif](../examples/hidden-motif.md) |
| Length warns and composition warns | Composition is probably following length. [length-bias](../examples/length-bias.md) |
| Duplication checks fail, everything statistical passes | A bookkeeping problem, not a biology one — and fixable without touching the data's content |
| Everything <span class="flag flag-unknown">Unknown</span> | The dataset is too small to check. Not the same as clean |
