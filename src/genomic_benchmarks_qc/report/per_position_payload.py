"""The data behind the interactive per-position figure, and its markup.

The per-position panels used to go into the report as a PNG. They now go in as
the numbers they were drawn from, and the page draws them in a canvas, so a
reader can zoom into a flagged position instead of squinting at a 1px band in a
picture. The numbers are all computed already - this module only reshapes what
`flag_significant_differences` and the statistics objects produced:

    freq      per label, per base, per position - the lines themselves
    coverage  per label, the proportion of that class reaching each position
    flags     per base, only the positions that were flagged

The window is the compared one, so an unflagged stretch of the figure is a
stretch that passed. Positions past it are still reported - as Unknown, in
`gb-qc-report.csv` and in the count of checks that were not scored - but they are not
drawn: a figure that mixed the two would have a tail where a missing flag means
nothing, and saying so took a wash and a caption that were most of what the
figure had to explain about itself. What is left to say about them is a sentence
in the section's explanation, where a reader who wants it can open it.

When no position clears the cohort floor there is no compared window, and the
figure falls back to the reported one - drawn, but with every position Unknown.
That is what the rest of the report does with an underpowered comparison: the
plots are still made, because a distribution is worth looking at whether or not
it was scored, and the flags say Unknown rather than Pass. The fallback keeps
that promise for the per-position panels without ever mixing the two cases in
one figure: either every position drawn was compared, or none of them was.

The report shows the flags and the frequencies behind them, and no separate
separability score: for a per-position base the AU-ROC is a restatement of the
gap between the two frequencies, so printing it beside them says the same thing
twice. It stays in `gb-qc-report.csv`, which is where a number to compute on belongs.

Nothing here rounds or subsets for size beyond three decimals, because the
payload grows with the *length* of the compared window and not with the number
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
from genomic_benchmarks_qc.utils.testing import position_windows

# One name per direction, shared by the statistics frame, the check names in the
# results table ('{name} - {base} position {n}') and the plot file names.
FEATURE_NAMES = {
    'forward': 'Per position nucleotide content',
    'reversed': 'Per position reversed nucleotide content',
}

# What the x axis counts, in both directions. The viewer puts the same wording on
# the position column of its flag table, so a reader never has to work out which
# end position 1 is: 'in reversed sequence' named the transform rather than the
# measurement, and left that open.
X_LABELS = {
    'forward': 'Position in sequence',
    'reversed': 'Position from sequence end',
}

# Flag colors, for the bands drawn over the panels and for the chips the viewer
# puts in its hover card. Fail and Warning match the shading of the static
# plots; Unknown is not a severity, it says the position was not scored, so it
# reads as a flat grey wash instead. Inside a compared window Unknown is only
# reached when the results table is missing a check the window says should be
# there; the wash covers the whole figure in the one case where nothing could be
# compared at all. Pass has no color here: it is never drawn as a band - it would
# cover the whole window - and the hover card no longer names it either, because
# every position the figure draws was scored, so what the card marks is a
# finding.
FLAG_COLORS = {
    'Fail': FAIL_COLOR,
    'Warning': WARN_COLOR,
    'Unknown': UNKNOWN_COLOR,
}

# The flags a position is worth carrying in the payload. Pass is the common
# case and would be most of the bytes, so it is left out and the viewer reads a
# position with no entry as a Pass.
STORED_FLAGS = ('Fail', 'Warning', 'Unknown')

DECIMALS = 3


def drawn_window(stats1, stats2):
    """The window the per-position figures draw, given two classes.

    The compared window, which is what makes an unflagged position on the figure
    a position that passed. When nothing could be compared there is no such
    window, and rather than drop the figure the panels fall back to the reported
    one: an underpowered comparison still gets its plots everywhere else in the
    report, and every position in the fallback is Unknown, so nothing in it can
    be misread as having passed.

    Returns:
        Last position to draw, 1-based and inclusive; 0 when the comparison has
        no per-position checks at all and there is genuinely nothing to draw.
    """
    end_position, scored_end_position = position_windows(stats1, stats2)
    return scored_end_position if scored_end_position >= 1 else end_position


def _series(frame, base, end_position):
    """One label's frequency curve for one base, as plain floats.

    Sliced positionally, like the static plot does, so the curve in the report
    and the curve in the PNG are the same numbers.
    """
    values = np.asarray(frame[base], dtype=float)[:end_position]
    return [round(float(value), DECIMALS) for value in values]


def _coverage(stats, end_position):
    """Proportion of one class's sequences reaching each position, 1-based.

    One curve per class rather than one pooled curve: this is the denominator
    behind that class's frequencies, and it is also what ends the compared
    window - which happens where the *lower* of the two curves drops through the
    coverage floor, something a pooled curve cannot show.
    """
    lengths = np.asarray(stats.stats['Sequence lengths']).flatten()
    if lengths.size == 0:
        return [0.0] * end_position
    positions = np.arange(1, end_position + 1)
    # A sequence reaches position p iff its length is >= p, so the count is the
    # size of the sorted tail from p onwards - one binary search per position
    # rather than a positions-by-sequences comparison matrix, which for a wide
    # window and a large cohort is the biggest allocation in the payload.
    ordered = np.sort(lengths)
    reaching = (ordered.size - np.searchsorted(ordered, positions, side='left')) / ordered.size
    return [round(float(value), DECIMALS) for value in reaching]


def _label(stats):
    """The name to show for a class: its label, or its filename when unlabelled."""
    return str(stats.label if stats.label is not None else stats.filename)


def _flags_by_position(results, prefix, bases, end_position):
    """Pull the per-position flag out of the results table.

    `results` is indexed by check name, and the per-position checks are named
    '<prefix> - <base> position <n>' with n 1-based. Looked up in one reindex per
    base rather than row by row, because a wide window is thousands of rows.

    Returns:
        Dict of base -> {position: flag}, holding only the positions that were
        flagged. A position missing from it passed; that is the reading the
        viewer relies on, so anything the results table does not name is recorded
        as Unknown rather than dropped.
    """
    flags = {}
    for base in bases:
        names = [f'{prefix} - {base} position {position}'
                 for position in range(1, end_position + 1)]
        rows = results.reindex(names)

        base_flags = {}
        if 'Flag' in rows:
            for position, flag in enumerate(rows['Flag'].tolist(), start=1):
                if flag == 'Pass':
                    continue
                base_flags[str(position)] = flag if flag in STORED_FLAGS else 'Unknown'
        flags[base] = base_flags
    return flags


def build_payload(stats1, stats2, bases, end_position, results, direction):
    """Assemble the payload for one direction of the per-position figure.

    Args:
        stats1, stats2: SequenceStatistics objects for the two classes.
        bases: Bases present in both classes, in plotting order.
        end_position: Last position to draw, 1-based and inclusive, as
            `drawn_window` resolves it - normally the compared window, and the
            reported one when nothing could be compared.
        results: DataFrame of every check's metrics, indexed by check name.
        direction: 'forward' or 'reversed'.

    Returns:
        A JSON-serialisable dict, or None when there is nothing to draw: no base
        in common, or no per-position check at all.
    """
    if direction not in FEATURE_NAMES:
        raise ValueError(f"Unknown direction: {direction}")
    if not bases or end_position < 1:
        return None

    frame1 = stats1.stats[FEATURE_NAMES[direction]]
    frame2 = stats2.stats[FEATURE_NAMES[direction]]
    flags = _flags_by_position(results, FEATURE_NAMES[direction], bases, end_position)

    # Whether anything in the window was compared. False only in the fallback
    # case, where the figure is the reported window with every position Unknown;
    # the viewer needs it to say 'nothing here was compared' where it would
    # otherwise say 'every position passed'.
    _, scored_end_position = position_windows(stats1, stats2)

    return {
        'direction': direction,
        'endPosition': int(end_position),
        'compared': bool(scored_end_position >= 1),
        'xLabel': X_LABELS[direction],
        'labels': [_label(stats1), _label(stats2)],
        'colors': list(CLASS_COLORS),
        'flagColors': dict(FLAG_COLORS),
        'counts': [int(stats1.stats['Number of sequences']),
                   int(stats2.stats['Number of sequences'])],
        'nucleotides': list(bases),
        'freq': {base: [_series(frame1, base, end_position),
                        _series(frame2, base, end_position)] for base in bases},
        'coverage': [_coverage(stats1, end_position), _coverage(stats2, end_position)],
        'flags': flags,
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
    <button type="button" class="qc-btn" data-action="prev" title="Previous flagged position">&#9664; Prev flag</button>
    <button type="button" class="qc-btn" data-action="next" title="Next flagged position">Next flag &#9654;</button>
    <span class="ppv-readout" aria-live="polite"></span>
    <span class="ppv-bar-gap"></span>
    <button type="button" class="qc-btn" data-action="reset" title="Reset zoom (or double-click the figure)">Reset zoom</button>
    <button type="button" class="qc-btn" data-action="save" title="Save the current window as PNG">&#8595; Save view</button>
  </div>
  <div class="ppv-plot">
    <canvas class="ppv-canvas" role="img" aria-label="{aria}"></canvas>
    <div class="ppv-tooltip" hidden></div>
    <p class="ppv-fallback">This figure is drawn in the browser and needs JavaScript.
    The flagged positions are listed below, and the same plot is in this report's
    <code>plots/</code> directory as a PNG.</p>
  </div>
  <p class="ppv-hint">Drag to zoom &middot; shift-drag to pan &middot; scroll to zoom &middot;
  double-click to reset &middot; hover for the frequencies and flags at a position</p>
  <details class="qc-panel">
    <summary><span class="ppv-flags-count">Flagged positions</span></summary>
    <div class="qc-panel-body ppv-flags-body"></div>
  </details>
</div>'''
