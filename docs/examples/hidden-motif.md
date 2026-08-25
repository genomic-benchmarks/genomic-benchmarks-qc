# hidden-motif

**A bias six positions wide, inside a 398-position sequence.**

Arabidopsis splice donor sites: 2000 sequences per class, every one 398
nucleotides long, from
[OmniGenBench](https://doi.org/10.48550/arXiv.2505.14402)'s
`splicing_arabidopsis_thaliana_donor`.

Length, GC content, duplication and vocabulary are all clean. Per-sequence
composition is clean. The only thing wrong with this dataset is what happens
between positions 199 and 205 — and that is enough for a classifier looking at a
single position to reach AU-ROC **0.783**.

## Run it

```bash
gb-qc evaluate-classes \
  --input examples/hidden-motif/data/train.csv \
  --out-folder qc-out
```

```bash
gb-qc evaluate-splits \
  --train-input examples/hidden-motif/data/train.csv \
  --test-input examples/hidden-motif/data/test.csv \
  --sequence-column sequence \
  --out-folder qc-out
```

## What it produces

--8<-- "_generated/hidden-motif-flags.md"

## What you should conclude

Six of the nine checks pass, and one of the two that fail is the mirror of the
other — `Per position reversed nucleotide content` is the same comparison read
from the far end, which for fixed-length sequences finds the same thing. So the
finding is singular: **one region of these sequences gives the class away, and
nothing else does.**

That region is the splice donor site, and the report locates it to the position:

--8<-- "_generated/hidden-motif-positions.md"

Six positions out of 398. **The other 392 are clean.** The reversed table is the
same finding read from the far end of a fixed-length sequence — 398 − 200 + 1 =
199, and the AU-ROCs match exactly.

### The one position in the motif that passes

Position 201 is missing from that table, and it is the most interesting thing in
it.

Every sequence the report was built from — all 2,000 in each class — has a `G` at
position 201, and position 202 is only ever `T` or `C`. Every window here,
positive or negative, is centred on a `GT` or a `GC`. The positives are `GT` in
99.0% of cases; the negatives are `GT` in 61.8% and `GC` in 38.3%. **The
negatives are not background genome — they are decoy donor-like sites that are
not donors**, which is what makes this dataset worth benchmarking on.

And it makes position 201 carry no information whatsoever about the label. The
check reports it <span class="flag flag-pass">Pass</span>.

That is the clearest illustration in these examples of what a per-position check
actually measures. **A position can be perfectly conserved and still pass.** The
question is never "is this position informative about the sequence" but "does it
differ between the classes" — and the single most conserved base in the dataset
is the one the report has nothing to say about.

What the flags mark, then, is the flanking consensus on either side of that
shared anchor:

```text
position    198  199  200 | 201  202  203  204  205  206
positives   A/C   A    G  |  G    T    A    A    G    T
negatives    ·    ·    ·  |  G   T/C   ·    ·    ·    ·
flagged           ✓    ✓  |       ✓    ✓    ✓    ✓
                          ^ exon / intron boundary
```

Which is `(A/C)AG|GTAAGT`, the canonical U2 donor site, exactly where a dataset of
donor sites should put it. Note that it is wider than the flagged window:
position 198 and position 206 both differ between the classes, just not by enough
to cross 0.6. **Flags mark where a difference is large enough to detect, not
where the motif ends.**

Whether that is a problem is a judgement the tool cannot make for you, and this
is the example that makes the distinction concrete:

- **If you are building a splice-site classifier**, this is the signal, not a
  bias. A model that finds the donor motif has done its job. `gb-qc` is telling
  you where the learnable signal is concentrated, which is worth knowing — a
  model scoring well here may have learned six positions and nothing else.
- **If you did not expect a positional give-away**, you have just found one, and
  it took a per-position check to see it. Nothing in the summary statistics
  hints at it.

The <span class="flag flag-warn">Warning</span> on
`Per sequence dinucleotide content` (0.602) is the same motif showing up
second-hand: enough sequences carry GT at a fixed position that the whole-sequence
GT frequency shifts slightly. It is a shadow of the real finding, not a separate
one.

## Why this is the example for the interactive plot

The flagged window, positions 199 to 205, is under two percent of a figure 398
positions wide. In the static PNG it is a hairline you would not notice
and could not read.

Open the report and use the per-position panel:

1. Land on the panel. The line looks flat.
2. Hit **Next flag**. The view jumps to position 199.
3. Drag to zoom into 195–210. Now the donor site is unmistakable — one class is
   78% G at position 200 where the other is 21%.
4. Hover any position for the per-class base frequencies behind the flag.

That is the whole argument for making the plot interactive rather than shipping
another PNG. This dataset is what it looks like when the argument is true.

## No leakage here

The split report comes back <span class="flag flag-pass">Pass</span> at
**0.00%** — not one test sequence has a 90%-similar match in the training set.
Worth seeing at least once, because it is the only example that is completely
clean on leakage, and it establishes that the leakage numbers elsewhere are real
findings rather than an artefact of how the check works.
