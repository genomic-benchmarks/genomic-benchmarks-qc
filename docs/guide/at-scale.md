# Running across many datasets

`gb-qc` handles one comparison at a time. Auditing a whole benchmark collection —
every dataset, every split — is a matter of arranging many independent runs, and
the arrangement is where the practical problems are.

The QC survey behind the paper drove **234 dataset splits** across seven
collections through this, which is where the notes below come from.

## The layout does the work

Both commands write into `<out-folder>/class/` and `<out-folder>/split/`, and
never collide. Everything above that is yours to arrange through `--out-folder`,
which is the whole mechanism you need:

```text
results/<collection>/<dataset>/<split>/
├── class/<column>/<classA>_vs_<classB>/
└── split/<column>/<train>_vs_<test>/
```

Give every run its own `--out-folder` and the results end up side by side,
comparable, with no naming scheme to invent:

```bash
for dataset in "${datasets[@]}"; do
  gb-qc evaluate-classes \
    --input "data/${dataset}_train.csv" \
    --out-folder "results/${collection}/${dataset}/train" \
    --log-file "results/${collection}/${dataset}/train/gb-qc.log"
done
```

`--log-file` per run matters more than it looks. At 234 runs you will not be
watching the terminal, and the log is where "this check was skipped, and why"
gets recorded.

## Collect the CSVs, not the HTML

The reports are for reading one at a time. For a survey, the CSVs are the data:

```python
import csv, pathlib

rows = []
for report in pathlib.Path('results').rglob('class/**/gb-qc-report.csv'):
    collection, dataset, split = report.relative_to('results').parts[:3]
    with report.open(newline='') as handle:
        for r in csv.DictReader(handle):
            if ' - ' in r['Check']:
                continue                       # skip the per-position breakdown
            rows.append({'collection': collection, 'dataset': dataset,
                         'split': split, **r})

with open('combined.csv', 'w', newline='') as out:
    writer = csv.DictWriter(out, fieldnames=rows[0])
    writer.writeheader()
    writer.writerows(rows)
```

Keep the unfiltered version too. The headline flags answer "which datasets have a
problem"; the per-position rows answer "and where", and you will want them once
something interesting turns up. In the paper's run the headline table came to
2,100 rows and the full one to 2.6 million — both worth having, for different
questions.

## Sizing the jobs

From the paper's run, in a workflow manager with per-job resources:

- **Runtime scales with the per-position window**, roughly an hour per 1,000 nt
  of window on the widest datasets. Budget on window width, not row count.
- **Memory** for `evaluate-classes` stayed under about 1 GB for datasets up to a
  few hundred nucleotides. The wide ones were not measured carefully enough to
  quote.
- **`evaluate-splits` is the long pole.** MMseqs2 dominates, and a few datasets
  needed several times the default time budget. Identify them and give them their
  own resource class rather than raising the limit for everything.
- **Pin `--end-position`** across a collection so datasets are comparable and the
  wide ones do not run away. It also makes the runtime predictable, which is what
  a scheduler wants.

## Things that will bite you

**Directory-name collisions.** Directory names come from your class labels and
column names, lowercased and stripped of filesystem-unsafe characters. Two labels
that differ only in case or punctuation collapse to the same name; `gb-qc` makes
them unique and warns which name it used. At scale you will not see the warning,
so check the logs for it.

**`Unknown` is not `Pass`.** Across many datasets some will fall below 250
sequences per class and stop being checked. An aggregation that counts
`Fail` rows will score those as clean. Count `Unknown` explicitly and report it as
its own category.

**Concurrent runs sharing an `--out-folder`.** Safe: each `evaluate-splits` run
gets its own `gb-qc-mmseqs-*` temporary directory. But give each run its own
output folder anyway, or you cannot tell which log belongs to which result.

**Version pinning.** Flags are thresholds on a computed score. A survey spanning a
version change is not internally comparable. Pin it and record the version beside
the results.
