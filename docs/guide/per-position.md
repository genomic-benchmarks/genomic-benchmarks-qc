# The per-position plot

The per-position checks produce more numbers than any other part of the report,
and a static figure is the wrong way to read them. For a 400-nucleotide dataset
there are over 1,500 comparisons behind one headline flag; for a 2,000-nucleotide
one, ten times that.

So the per-position panels in the HTML report are interactive. This page is about
using them.

## The problem it solves

Open [hidden-motif](../examples/hidden-motif.md)'s report and look at the
per-position panel before doing anything else. It looks flat. Nothing obviously
wrong.

That dataset's `Per position nucleotide content` check is a
<span class="flag flag-fail">Fail</span> at AU-ROC 0.783 — the strongest
per-position finding in any of the examples. It comes from **six positions out of
398**, clustered between 199 and 205. At the default zoom that region is about
one and a half percent of the x-axis: a hairline you would not notice and could
not read a value off.

That is the case the interactivity exists for. Not for prettiness — for the fact
that the finding is genuinely invisible until you zoom.

## Controls

| Control | Does |
|---|---|
| **Drag** | Zoom into a horizontal range |
| **Shift-drag** | Pan |
| **Scroll** | Zoom in and out |
| **Double-click** | Reset the view |
| **◀ Prev flag** / **Next flag ▶** | Jump to the previous or next flagged position |
| **Reset zoom** | Same as double-click |
| **↓ Save view** | Write the current window out as a PNG |
| **?** | Inline explanation of what the panel shows |
| **Hover** | Per-class base frequencies and the flag at that position |

There is also a two-button filter above the table below the plot: **All**
positions, or only those needing attention. On a long dataset the second is
usually what you want — [variable-length](../examples/variable-length.md) has
26,408 per-position sub-checks, and scrolling all of them is not a plan.

## A way to read it

This works on any dataset and takes about a minute.

1. **Look at the panel unzoomed.** You are asking one question: is the difference
   *spread* or *localised*? A broad elevated region is composition showing
   through; a spike is a motif or an artefact. A spike at position 1, or at
   position 1 of the reversed panel, is almost always an artefact of how the
   sequences were cut.
2. **Hit Next flag.** This jumps to the first flagged position, whatever the zoom.
   It is the fastest way to find out whether the flags cluster or scatter.
3. **Drag to zoom around the first cluster** — ten or twenty positions either
   side. Now the per-class lines separate and you can see which base is doing it.
4. **Hover the worst position.** You get both classes' base frequencies. This is
   where "class 1 is 94% G here, class 0 is 26%" becomes concrete rather than an
   AU-ROC.
5. **Keep hitting Next flag** to the end. Six clustered positions is one finding.
   Forty scattered ones is composition. Position 1 by itself is a technical
   artefact.
6. **Save view** on whatever you want to keep, and use the table underneath to
   jump back to specific positions.

## What the window means

Two boundaries decide what you can see, and they are different things.

**`--min-coverage`** (25% by default) decides how far positions can be
*flagged*. A position is only flagged where at least that fraction of each class
still reaches it — and at least 250 sequences, whichever binds harder. **This is
also the window the figures draw.** So if the plot stops before your longest
sequence ends, this is why.

**`--end-position`** decides how far the checks *run at all*. By default it is
the last position that at least 50 sequences of each class reach. It can only
narrow what gets flagged, never widen it.

Everything past the flagging window is reported
<span class="flag flag-unknown">Unknown</span> and is not drawn. That is not the
plot hiding data from you — it is the tool declining to make a comparison it
cannot make honestly. The reasoning is in
[how a flag is decided](how-it-works.md#at-least-25-of-the-class-must-reach-a-position);
the short version is that a cohort far along a variable-length dataset is a
sample of the class's longest sequences, not of the class.

For a dataset where this dominates, see
[variable-length](../examples/variable-length.md): flagging stops at position 541,
scoring at 2,636, and the longest sequence runs to 5,964. **79% of its
per-position sub-checks are <span class="flag flag-unknown">Unknown</span>.**

### Widening it

`--min-coverage 0` drops the fraction rule and leaves only the 250-sequence
floor, which is not switchable off. That extends flagging as far as the counts
allow — on `variable-length`, from 541 to 2,636. Read what comes back with the
caveat above in mind: those positions are real comparisons on real sequences, but
the sequences are the long ones.

### Narrowing it

`--end-position N` is worth reaching for in three cases:

- **Long sequences.** [length-bias](../examples/length-bias.md) holds whole mRNAs,
  median 3,463 nt and up to 17,497. Unbounded, its per-position checks would run
  for thousands of positions and produce a report too heavy to open, for
  positions only a handful of transcripts reach. It uses `--end-position 500`.
- **Comparing two reports.** Pinning the same window across runs makes two reports
  line up position by position.
- **Report size.** The per-position payload is the bulk of a report's bytes.

## Both directions

Each report has two per-position panels: forward, and counted from the end of
each sequence. On fixed-length data they find the same thing and report identical
numbers — [hidden-motif](../examples/hidden-motif.md) shows 0.783 for both, and
its flagged positions mirror (199 forward ↔ 200 reversed on a 398-nt sequence).

The reversed panel earns its place on variable-length data, where anything
anchored to the *end* of a sequence — a poly-A tail, a 3' adapter — sits at a
different forward position in every sequence and the forward panel cannot see it.
If the reversed check fires and the forward one does not, that is what to look
for.

[enhancers](../examples/enhancers.md#the-only-two-flagged-positions-are-the-two-ends)
is the case where it pays off. Sequences run 4 to 568 nucleotides; both panels
flag position 1, and they are flagging *different bases* — the first of each
sequence and the last. Only the reversed panel could have found the second.
