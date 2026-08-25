# Using it in CI

`gb-qc` writes a CSV beside every report for exactly this: so a build can read
the flags and decide something. The point is not to block merges on every
<span class="flag flag-warn">Warning</span> — it is to notice the day a dataset
regenerates differently and nobody looks.

## The CSV

`gb-qc-report.csv` from `evaluate-classes`:

```csv
Check,Flag,AU-ROC,AU-PR,Accuracy
Unique bases,Pass,,,
Sequence lengths,Pass,0.5,0.5,0.5
Per sequence GC content,Fail,0.701,0.7266,0.6429
Per sequence nucleotide content - A,Fail,0.7043,0.7141,0.6531
```

and from `evaluate-splits`:

```csv
Check,Flag,Percentage of leaked queries,Percentage of leaked targets
Data Leakage,Fail,6.04%,0.98%
```

Rows whose `Check` contains ` - ` are the per-base and per-position breakdown of
the check above them. **Filter those out for a gate**.

## A gate

```python
#!/usr/bin/env python3
"""Fail the build if any headline check regressed. Usage: gate.py <report.csv>"""
import csv, sys

ALLOWED = {'Pass', 'Warning'}          # tighten to {'Pass'} when you can
report = sys.argv[1]

with open(report, newline='') as handle:
    rows = [r for r in csv.DictReader(handle) if ' - ' not in r['Check']]

bad = [r for r in rows if r['Flag'] not in ALLOWED]
for r in bad:
    print(f"::error::{r['Check']}: {r['Flag']} (AU-ROC {r.get('AU-ROC') or 'n/a'})")

unknown = [r for r in rows if r['Flag'] == 'Unknown']
if unknown:
    # Unknown is not Pass - it means the check never ran. Surface it.
    print(f"::warning::{len(unknown)} check(s) not scored: "
          + ', '.join(r['Check'] for r in unknown))

sys.exit(1 if bad else 0)
```

```yaml
- name: Dataset QC
  run: |
    pip install genomic-benchmarks-qc
    gb-qc evaluate-classes --input data/train.csv --out-folder qc-out
    python gate.py qc-out/class/sequence/0_vs_1/gb-qc-report.csv

- name: Keep the report
  if: always()
  uses: actions/upload-artifact@v4
  with:
    name: qc-report
    path: qc-out/**/gb-qc-report.html
```

Uploading the HTML on `always()` is the part people skip and then regret: when
the gate fails, the CSV tells you *that* a check failed and the report tells you
*why*.

## Pinning flags instead of thresholding them

A gate that allows <span class="flag flag-warn">Warning</span> everywhere passes
a dataset that has quietly got worse. The alternative is to record the flags you
expect and fail on any change — in either direction, so an improvement is also
something you have to acknowledge.

This repository does that for its own examples. Each one declares its flags in
`meta.toml`:

```toml
[expect."class/sequence/0_vs_1"]
"Per sequence GC content" = "Fail"
"Per position nucleotide content" = "Fail"
```

and `examples/build.py --check` compares them against the generated reports. The
docs build runs it, so a change in the tool's behaviour fails the build rather
than silently making a published page wrong. The
[implementation](https://github.com/genomic-benchmarks/genomic-benchmarks-qc/blob/main/examples/build.py)
is about forty lines and worth copying.

## Practical notes

- **Pin the version.** `pip install genomic-benchmarks-qc==x.y.z`. Flags are
  thresholds on a computed score; a change in either moves them, and you want
  that to be a deliberate upgrade rather than a surprise on a Tuesday.
- **`evaluate-classes` needs no MMseqs2.** If you only gate on class checks, skip
  installing it and save a couple of minutes per run.
- **Pin the window** with `--end-position` if you compare reports across runs, so
  two runs line up position by position.
- **Watch for `Unknown`.** A dataset that shrinks below 250 sequences per class
  stops being checked, and a gate that only looks for `Fail` will call that a
  pass. The script above prints a warning for it.
