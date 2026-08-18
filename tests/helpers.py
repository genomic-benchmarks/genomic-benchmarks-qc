"""Builders for tests that drive the analysis pipelines end to end.

The CLI-level tests stub the pipelines out; these helpers are for the tests that
actually run them, and need real input files and realistic MMseqs2 output.
"""

import random

import pandas as pd

from genomic_benchmarks_qc.utils.mmseqs_summary import MMSEQS_REQUIRED_COLS


def sequences(count, length=60, seed=0):
    """Return `count` reproducible random nucleotide sequences."""
    rng = random.Random(seed)
    return [''.join(rng.choice('ACGT') for _ in range(length)) for _ in range(count)]


def write_fasta(path, count, seed):
    path.parent.mkdir(parents=True, exist_ok=True)
    records = [f">seq{i}\n{seq}" for i, seq in enumerate(sequences(count, seed=seed))]
    path.write_text('\n'.join(records) + '\n')
    return str(path)


def write_csv(path, labels, columns=('sequence',), rows_per_label=30):
    """Write a CSV with `rows_per_label` rows for each label in `labels`."""
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {column: [] for column in columns}
    data['label'] = []
    for seed, label in enumerate(labels):
        for column in columns:
            data[column] += sequences(rows_per_label, seed=seed + len(column))
        data['label'] += [label] * rows_per_label
    pd.DataFrame(data).to_csv(path, index=False)
    return str(path)


def write_regression_csv(path, values):
    """Write a CSV whose label column holds regression targets."""
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {'sequence': sequences(len(values), seed=1), 'label': values}
    ).to_csv(path, index=False)
    return str(path)


def mmseqs_hit(query, target, pident=96.9, qcov=0.99, tcov=0.99, alnlen=60):
    """One row of MMseqs2 easy-search output, in the columns the tool requests.

    `qcov`/`tcov` are fractions and `pident` a percentage, matching the real
    `--format-mode 4` output; the leakage score is min(qcov, tcov) * pident.
    """
    alignment = 'A' * alnlen
    return {
        'query': query, 'target': target,
        'qcov': qcov, 'tcov': tcov, 'pident': pident,
        'evalue': 1e-200,
        'qstart': 1, 'qend': alnlen, 'tstart': 1, 'tend': alnlen, 'alnlen': alnlen,
        'qaln': alignment, 'taln': alignment,
    }


def write_mmseqs_output(path, hits=()):
    """Write an MMseqs2 result table, with the header the real run produces."""
    pd.DataFrame(list(hits), columns=MMSEQS_REQUIRED_COLS).to_csv(path, sep='\t', index=False)
    return path
