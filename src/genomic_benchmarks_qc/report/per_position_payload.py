"""The data behind the interactive per-position figure, and its markup.

The per-position panels used to go into the report as a PNG. They now go in as
the numbers they were drawn from, and the page draws them in a canvas, so a
reader can zoom into a flagged position instead of squinting at a 1px band in a
picture. The numbers are all computed already - this module only reshapes what
`flag_significant_differences` and the statistics objects produced:

    freq      per label, per base, per position - the lines themselves
    coverage  the proportion of sequences reaching each position
    auroc     per base, per position - how well that base separates the labels
    flags     per base, only the positions that were flagged, or not scored

Nothing here rounds or subsets for size beyond three decimals, because the
payload grows with the *length* of the analysed window and not with the number
of sequences: it is ~75 bytes per position per direction, against ~900 KB for
the PNG it replaces.
"""

import json

import numpy as np

from genomic_benchmarks_qc.report.colors import (
    CLASS_COLORS,
    FAIL_COLOR,
    UNKNOWN_COLOR,
    WARN_COLOR,
)

# One name per direction, shared by the statistics frame, the check names in the
# results table ('{name} - {base} position {n}') and the plot file names.
FEATURE_NAMES = {
    'forward': 'Per position nucleotide content',
    'reversed': 'Per position reversed nucleotide content',
}

X_LABELS = {
    'forward': 'Position in sequence',
    'reversed': 'Position in reversed sequence',
}

# Flag colors used for the bands drawn over the panels. Fail and Warning match
# the shading of the static plots; Unknown is not a severity, it says the
# position was not scored, so it reads as a flat grey wash instead.
FLAG_COLORS = {
    'Fail': FAIL_COLOR,
    'Warning': WARN_COLOR,
    'Unknown': UNKNOWN_COLOR,
}

DECIMALS = 3


def _series(frame, base, end_position):
    """One label's frequency curve for one base, as plain floats.

    Sliced positionally, like the static plot does, so the curve in the report
    and the curve in the PNG are the same numbers.
    """
    values = np.asarray(frame[base], dtype=float)[:end_position]
    return [round(float(value), DECIMALS) for value in values]


def _coverage(stats1, stats2, end_position):
    """Proportion of all sequences reaching each position, 1-based.

    Pooled over both labels, which is what the static plot's bottom panel shows:
    it is the denominator behind the curves above it.
    """
    lengths = np.concatenate([
        np.asarray(stats1.stats['Sequence lengths']).flatten(),
        np.asarray(stats2.stats['Sequence lengths']).flatten(),
    ])
    if lengths.size == 0:
        return [0.0] * end_position
    positions = np.arange(1, end_position + 1)
    reaching = (lengths[None, :] >= positions[:, None]).sum(axis=1) / lengths.size
    return [round(float(value), DECIMALS) for value in reaching]


def _label(stats):
    """The name to show for a class: its label, or its filename when unlabelled."""
    return str(stats.label if stats.label is not None else stats.filename)


def _metrics_by_position(results, prefix, bases, end_position):
    """Pull the per-position AU-ROC and flag out of the results table.

    `results` is indexed by check name, and the per-position checks are named
    '<prefix> - <base> position <n>' with n 1-based. Looked up in one reindex per
    base rather than row by row, because a wide window is thousands of rows.

    Returns:
        Tuple of (auroc, flags) where auroc maps base -> list of AU-ROC values or
        None per position, and flags maps base -> {position: flag} holding only
        the positions that were flagged or not scored.
    """
    auroc = {}
    flags = {}
    for base in bases:
        names = [f'{prefix} - {base} position {position}'
                 for position in range(1, end_position + 1)]
        rows = results.reindex(names)
        scores = rows['AU-ROC'].to_numpy(dtype=float) if 'AU-ROC' in rows else np.full(end_position, np.nan)
        auroc[base] = [None if not np.isfinite(value) else round(float(value), DECIMALS)
                       for value in scores]

        base_flags = {}
        if 'Flag' in rows:
            for position, flag in enumerate(rows['Flag'].tolist(), start=1):
                # Pass is the common case and says nothing worth drawing.
                if flag in FLAG_COLORS:
                    base_flags[str(position)] = flag
        flags[base] = base_flags
    return auroc, flags


def build_payload(stats1, stats2, bases, end_position, results, direction):
    """Assemble the payload for one direction of the per-position figure.

    Args:
        stats1, stats2: SequenceStatistics objects for the two classes.
        bases: Bases present in both classes, in plotting order.
        end_position: Last analysed position, 1-based and inclusive.
        results: DataFrame of every check's metrics, indexed by check name.
        direction: 'forward' or 'reversed'.

    Returns:
        A JSON-serialisable dict, or None when there is nothing to draw.
    """
    if direction not in FEATURE_NAMES:
        raise ValueError(f"Unknown direction: {direction}")
    if not bases or end_position < 1:
        return None

    frame1 = stats1.stats[FEATURE_NAMES[direction]]
    frame2 = stats2.stats[FEATURE_NAMES[direction]]
    auroc, flags = _metrics_by_position(results, FEATURE_NAMES[direction], bases, end_position)

    return {
        'direction': direction,
        'endPosition': int(end_position),
        'xLabel': X_LABELS[direction],
        'labels': [_label(stats1), _label(stats2)],
        'colors': list(CLASS_COLORS),
        'flagColors': dict(FLAG_COLORS),
        'counts': [int(stats1.stats['Number of sequences']),
                   int(stats2.stats['Number of sequences'])],
        'nucleotides': list(bases),
        'freq': {base: [_series(frame1, base, end_position),
                        _series(frame2, base, end_position)] for base in bases},
        'coverage': _coverage(stats1, stats2, end_position),
        'flags': flags,
        'auroc': auroc,
    }


def payload_script(payload, dom_id):
    """Return the <script type="application/json"> element holding a payload.

    JSON goes into the page as data rather than as a JavaScript literal, so a
    class label containing markup cannot break out of the script element. The
    three characters that could close it early are escaped.
    """
    text = (json.dumps(payload, separators=(',', ':'))
            .replace('<', '\\u003c')
            .replace('>', '\\u003e')
            .replace('&', '\\u0026'))
    return f'<script type="application/json" id="{dom_id}">{text}</script>'


def viewer_html(payload, dom_id):
    """Return the payload plus the figure's markup, ready to drop into a section.

    Args:
        payload: The dict from build_payload.
        dom_id: Element id prefix, e.g. 'ppv-fwd'; the payload element gets
            '<dom_id>-data' and the viewer reads it through data-payload.

    Returns:
        HTML snippet as a string.
    """
    labels = payload['labels']
    aria = (f"Per-position nucleotide frequency for {', '.join(payload['nucleotides'])} "
            f"across positions 1 to {payload['endPosition']}, comparing "
            f"{labels[0]} with {labels[1]}. Flagged positions are listed below the plot.")
    return f'''{payload_script(payload, dom_id + '-data')}
<div class="ppv" id="{dom_id}" data-payload="{dom_id}-data">
  <div class="ppv-bar">
    <button type="button" data-action="prev" title="Previous flagged position (p)">&#9664; Prev flag</button>
    <button type="button" data-action="next" title="Next flagged position (n)">Next flag &#9654;</button>
    <span class="ppv-readout" aria-live="polite"></span>
    <span class="ppv-bar-gap"></span>
    <button type="button" data-action="reset" title="Reset zoom (double-click or Home)">Reset zoom</button>
    <button type="button" data-action="save" title="Save the current window as PNG (s)">&#8595; Save view</button>
  </div>
  <div class="ppv-plot">
    <canvas class="ppv-canvas" tabindex="0" role="img" aria-label="{aria}"></canvas>
    <div class="ppv-tooltip" hidden></div>
    <p class="ppv-fallback">This figure is drawn in the browser and needs JavaScript.
    The flagged positions are listed below, and the same plot is in this report's
    <code>plots/</code> directory as a PNG.</p>
  </div>
  <p class="ppv-hint">Drag to zoom &middot; shift-drag to pan &middot; scroll to zoom &middot;
  double-click to reset &middot; hover for the values at a position &middot;
  keys: n / p for the next and previous flag, s to save the view</p>
  <details class="ppv-flags">
    <summary><span class="ppv-flags-count">Flagged positions</span></summary>
    <div class="ppv-flags-body"></div>
  </details>
</div>'''
