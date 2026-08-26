"""Rendering one MMseqs2 hit as a text alignment for the HTML report.

The search reports an alignment of a region of the query against a region of
the target; this lays the two full sequences out against each other with that
region aligned, so a reader can see how much of each sequence the match covers.
"""

from itertools import zip_longest

from genomic_benchmarks_qc.report.utils import escape_html_text


def has_reversed_coordinates(row):
    """True if the hit's coordinates run backwards on either sequence.

    The one failure of `build_alignment_string` with a known cause, so the
    report can name it instead of guessing. See
    `mmseqs_summary._reversed_coordinate_mask` for what produces these.
    """
    return int(row["qstart"]) > int(row["qend"]) or int(row["tstart"]) > int(row["tend"])


def _validate_alignment_coords(start, end, sequence_length, label):
    """Check one sequence's 1-based alignment range lies inside the sequence."""
    if start < 1:
        raise ValueError(f"{label} start < 1 (1-based expected): {start}")
    if end < start:
        raise ValueError(f"{label} end ({end}) < start ({start})")
    if end > sequence_length:
        raise ValueError(
            f"{label} coordinates [{start}, {end}] exceed sequence length {sequence_length}"
        )


def _validate_alignment_inputs(qstart, tstart, qend, tend, qseq, tseq, qaln, taln):
    """Check the reported alignment really describes these two sequences.

    The aligned strings, with their gaps removed, must be exactly the slices the
    coordinates point at. If they are not, the hit and the sequence it was
    matched back to have come apart somewhere, and rendering it would produce a
    confident-looking but wrong alignment.
    """
    _validate_alignment_coords(qstart, qend, len(qseq), "Query")
    _validate_alignment_coords(tstart, tend, len(tseq), "Target")

    if len(qaln) != len(taln):
        raise ValueError(
            f"Alignment length mismatch: "
            f"  len(qaln)={len(qaln)}, len(taln)={len(taln)}"
        )

    qstart0 = qstart - 1
    tstart0 = tstart - 1
    expected_q = qseq[qstart0:qend]
    expected_t = tseq[tstart0:tend]

    if qaln.replace("-", "") != expected_q:
        raise ValueError(
            "Query alignment mismatch:\n"
            "  Ungapped aligned sequence does not match original sequence slice."
        )

    if taln.replace("-", "") != expected_t:
        raise ValueError(
            "Target alignment mismatch:\n"
            "  Ungapped aligned sequence does not match original sequence slice."
        )


def _get_midline(qaln, taln):
    """Return the line between the two sequences: '|' match, '.' mismatch, ' ' gap."""
    return "".join(
        " " if (q == "-" or t == "-")
        else "|" if q == t
        else "."
        for q, t in zip(qaln, taln, strict=True)
    )


def _get_alignment(qstart, tstart, qseq, tseq, qaln, taln, tend, qend):
    """Lay the two full sequences out with their aligned regions lined up.

    Returns (target line, midline, query line), padded at the front so that the
    aligned region starts at the same column in all three, and with the
    unaligned flanks of each sequence kept around it.
    """
    qstart0 = qstart - 1
    tstart0 = tstart - 1

    q_left = ""
    t_left = ""

    delta = qstart0 - tstart0
    if delta > 0:
        t_left = " " * delta
    elif delta < 0:
        q_left = " " * (-delta)

    mid_left = " " * max(qstart0, tstart0)
    mid_line = mid_left + _get_midline(qaln, taln)

    new_taln = tseq[:tstart0] + taln + tseq[tend:]
    new_qaln = qseq[:qstart0] + qaln + qseq[qend:]

    return t_left + new_taln, mid_line, q_left + new_qaln


def _wrap_text(text, width):
    """Split text into fixed-width chunks, without breaking on whitespace."""
    return [text[i:i + width] for i in range(0, len(text), width)]


def _color_base(base):
    """Wrap one base in the span that colors it, escaping anything unexpected."""
    safe_base = escape_html_text(base)
    return {
        "A": '<span class="base-A">A</span>',
        "C": '<span class="base-C">C</span>',
        "G": '<span class="base-G">G</span>',
        "T": '<span class="base-T">T</span>',
        "-": '<span class="base-gap">-</span>',
    }.get(base, f'<span class="base-other">{safe_base}</span>')


def _format_sequence(seq, color):
    """Escape a sequence for HTML, coloring each base unless color is off."""
    if not color:
        return escape_html_text(seq)
    return "".join(_color_base(base) for base in seq)

def build_alignment_string(row, width=80, color=True, validate=True):
    """Render one hit as wrapped blocks of target, midline and query.

    `row` is one row of the MMseqs2 results with the two sequences attached.
    Returns HTML unless `color` is off, in which case it is plain text. Raises
    ValueError when `validate` is on and the hit does not match its sequences.
    """
    qstart = int(row["qstart"])
    tstart = int(row["tstart"])
    qend   = int(row["qend"])
    tend   = int(row["tend"])

    qseq = row["qseq"]
    tseq = row["tseq"]
    qaln = row["qaln"]
    taln = row["taln"]

    if validate:
        _validate_alignment_inputs(qstart, tstart, qend, tend, qseq, tseq, qaln, taln)

    t_line, mid_line, q_line = _get_alignment(
        qstart, tstart, qseq, tseq, qaln, taln, tend, qend
    )

    blocks = []
    for t, m, q in zip_longest(
        _wrap_text(t_line, width),
        _wrap_text(mid_line, width),
        _wrap_text(q_line, width),
        fillvalue="",
    ):
        blocks.append(f"T    {_format_sequence(t, color)}")
        blocks.append(f"     {m}")   # keep midline uncolored
        blocks.append(f"Q    {_format_sequence(q, color)}")
        blocks.append("")

    return "\n".join(blocks)
