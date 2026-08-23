# hidden-motif

**A bias eight nucleotides wide, inside a 398-position sequence.**

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

That region is the GT donor motif. These are splice sites, so of course position
200 is a G in the positive class — the biology is the label. Six positions carry
it:

| Position | 199 | 200 | 202 | 203 | 204 | 205 |
|---|---|---|---|---|---|---|
| Worst AU-ROC | 0.683 | **0.783** | 0.672 | 0.714 | 0.641 | 0.659 |
| Flag | <span class="flag flag-warn">Warning</span> | <span class="flag flag-fail">Fail</span> | <span class="flag flag-warn">Warning</span> | <span class="flag flag-fail">Fail</span> | <span class="flag flag-warn">Warning</span> | <span class="flag flag-warn">Warning</span> |

Six positions out of 398. **The other 392 are clean.**

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

Six flagged positions in a figure 398 positions wide is roughly one and a half
percent of the x-axis. In the static PNG it is a hairline you would not notice
and could not read.

Open the report and use the per-position panel:

1. Land on the panel. The line looks flat.
2. Hit **Next flag**. The view jumps to position 199.
3. Drag to zoom into 195–210. Now the donor site is unmistakable — one class is
   almost entirely G at 200, the other is background.
4. Hover any position for the per-class base frequencies behind the flag.

That is the whole argument for making the plot interactive rather than shipping
another PNG. This dataset is what it looks like when the argument is true.

## No leakage here

The split report comes back <span class="flag flag-pass">Pass</span> at
**0.00%** — not one test sequence has a 90%-similar match in the training set.
Worth seeing at least once, because it is the only example that is completely
clean on leakage, and it establishes that the leakage numbers elsewhere are real
findings rather than an artefact of how the check works.
