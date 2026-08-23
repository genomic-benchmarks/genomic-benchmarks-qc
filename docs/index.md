# genomic-benchmarks-qc

Automated quality control for genomic machine learning datasets. Finds the
biases, inconsistencies and data leakage that make a benchmark easier than it
looks — before you train on it.

A model that scores well on a biased dataset has learned the bias, not the
biology. `gb-qc` looks for the differences a classifier could exploit without
understanding anything: classes that differ in length or base composition, a
give-away at a single position, the same sequence in both classes, a test set
that repeats the training set.

```bash
pip install genomic-benchmarks-qc
```

```bash
gb-qc evaluate-classes \
  --input examples/enhancers/data/enhancers_train.csv \
  --input examples/enhancers/data/enhancers_test.csv \
  --out-folder qc-out
```

Every check gets a <span class="flag flag-pass">Pass</span>,
<span class="flag flag-warn">Warning</span> or
<span class="flag flag-fail">Fail</span> flag, a standalone HTML report you can
read or mail, and a CSV you can put in CI.

!!! note "Where things are, for now"

    This site is being built out. The [README][readme] is complete and remains
    the place to start: installation, every input format, the full output
    layout, and the CLI options in table form. What lives here is the
    [generated reference](reference/cli.md) and — as they land — the worked
    examples and the longer-form guides that a README has no room for.

  [readme]: https://github.com/genomic-benchmarks/genomic-benchmarks-qc#readme

## The examples

Eight datasets, each the only one that shows a particular thing, with flags
measured rather than asserted. They are the fastest way to see what a report
actually tells you.

| Example | What it shows | Report |
|---|---|---|
| `clean-dataset` | The control: what "nothing wrong" looks like, AU-ROC 0.50–0.54 throughout | [open](reports/clean-dataset/class/sequence/0_vs_1/gb-qc-report.html) |
| `composition-bias` | The worst case — six checks fail, and 6% of the test set is already in training | [open](reports/composition-bias/class/sequence/0_vs_1/gb-qc-report.html) |
| `hidden-motif` | A bias eight nucleotides wide inside 398 positions. Why the per-position plot is interactive | [open](reports/hidden-motif/class/sequence/0_vs_1/gb-qc-report.html) |
| `variable-length` | Sequences that stop at different places, and why most positions go unscored | [open](reports/variable-length/class/sequence/0_vs_1/gb-qc-report.html) |
| `length-bias` | Length alone separating the classes, on a continuous label | [open](reports/length-bias/class/sequence/high_vs_low/gb-qc-report.html) |
| `paired-sequences` | Two sequence columns in one row | [open](reports/paired-sequences/class/merged/0_vs_1/gb-qc-report.html) |
| `fasta-classes` | One FASTA file per class | [open](reports/fasta-classes/class/sequence/coding_seqs_vs_intergenomic_seqs/gb-qc-report.html) |
| `enhancers` | The quickstart dataset, with a little real train/test leakage | [open](reports/enhancers/class/sequence/0_vs_1/gb-qc-report.html) |

Their data and provenance are in
[`examples/`](https://github.com/genomic-benchmarks/genomic-benchmarks-qc/tree/main/examples).
The reports are built from that data by the same commands shown above, so what
you read here is always what the current code produces.

## Flags

Flags come from the AU-ROC of a classifier that sees only the one feature under
test:

| Flag | AU-ROC | Meaning |
|---|---|---|
| <span class="flag flag-pass">Pass</span> | ≤ 0.6 | Classes not distinguishable by this feature |
| <span class="flag flag-warn">Warning</span> | ≤ 0.7 | A model could get some traction here |
| <span class="flag flag-fail">Fail</span> | > 0.7 | Significant bias detected |
| <span class="flag flag-unknown">Unknown</span> | — | Not enough sequences to score the check |

<span class="flag flag-unknown">Unknown</span> is not
<span class="flag flag-pass">Pass</span>. It says the comparison was not made,
not that it came out clean — a check needs at least 250 sequences per class
before it is scored at all. The plots and descriptive statistics are computed
from all the data either way, so a small dataset can still be compared by eye.
