# Python API

The CLI is a thin wrapper. Both commands are one function call, so you can drive
them from a pipeline without shelling out.

```python
from genomic_benchmarks_qc.evaluate_classes import run as evaluate_classes
from genomic_benchmarks_qc.evaluate_splits import run as evaluate_splits
```

Full signatures and every argument are in the
[API reference](../reference/api/evaluate-classes.md).

## Comparing classes

```python
evaluate_classes(
    input=['examples/enhancers/data/enhancers_train.csv',
           'examples/enhancers/data/enhancers_test.csv'],
    format='csv',
    out_folder='qc-out',
)
```

`format` is required here — the CLI infers it from the file extension, but the
function does not. It takes `'csv'`, `'tsv'` or `'fasta'`; gzip is detected
per-file from the filename regardless.

Continuous labels, several sequence columns, and a trimmed per-position window:

```python
evaluate_classes(
    input=['pairs.tsv'],
    format='tsv',
    sequence_column=['gene', 'noncodingRNA'],
    label_column='target',
    regression=True,
    end_position=500,
    out_folder='qc-out',
)
```

## Checking a split

```python
evaluate_splits(
    train_files=['train.csv'],
    test_files=['test.csv'],
    format='csv',
    sequence_column=['sequence'],
    similarity_threshold=90.0,
    out_folder='qc-out',
)
```

Note `train_files` / `test_files` rather than the CLI's `--train-input` /
`--test-input`.

## Reading the results back

Neither function returns the results — they write files, and the CSV is the
interface. Reports land at
`<out_folder>/class/<column>/<classA>_vs_<classB>/` and
`<out_folder>/split/<column>/<train>_vs_<test>/`, with directory names derived
from your labels and column names.

Rather than reconstruct those paths, glob for the report:

```python
import csv, pathlib

for report in pathlib.Path('qc-out').rglob('gb-qc-report.csv'):
    with report.open(newline='') as handle:
        flags = {r['Check']: r['Flag'] for r in csv.DictReader(handle)
                 if ' - ' not in r['Check']}          # headline checks only
    failed = [c for c, f in flags.items() if f == 'Fail']
    print(report.parent.name, failed or 'clean')
```

Adding `'json'` to `report_types` also writes per-class statistics — counts, GC,
lengths, base and dinucleotide frequencies — to
`<out_folder>/class/<column>/per-class/<class>.json`, which is the machine-readable
route to the numbers behind the flags.

## Logging and failures

Both functions call `setup_logger` and log progress at `INFO`. Pass
`log_level='WARNING'` for quiet, or `log_file='qc.log'` to tee to a file.

They raise on failure rather than exiting, so wrap them as you would any other
call. The CLI's non-zero exit codes are its own translation layer, not the
functions' behaviour.
