from itertools import zip_longest
import html

# Filter fasta file, keeping only sequences with IDs in ids_to_keep (for hits)
def filter_fasta_by_ids(fasta_path, new_fasta_path, ids_to_keep):
    ids_to_keep = set(ids_to_keep)
    with fasta_path.open() as file_in, new_fasta_path.open("w") as file_out:
        write = False

        for line in file_in:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                write = seq_id in ids_to_keep

            if write:
                file_out.write(line)


def _validate_alignment_coords(start, end, sequence_length, label):
    if start < 1:
        raise ValueError(f"{label} start < 1 (1-based expected): {start}")
    if end < start:
        raise ValueError(f"{label} end ({end}) < start ({start})")
    if end > sequence_length:
        raise ValueError(
            f"{label} coordinates [{start}, {end}] exceed sequence length {sequence_length}"
        )


def _validate_alignment_inputs(qstart, tstart, qend, tend, qseq, tseq, qaln, taln):
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
    return "".join(
        " " if (q == "-" or t == "-")
        else "|" if q == t
        else "."
        for q, t in zip(qaln, taln)
    )


def _get_alignment(qstart, tstart, qseq, tseq, qaln, taln, tend, qend):
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
    return [text[i:i + width] for i in range(0, len(text), width)]


def _color_base(base):
    safe_base = html.escape(base)
    return {
        "A": '<span class="base-A">A</span>',
        "C": '<span class="base-C">C</span>',
        "G": '<span class="base-G">G</span>',
        "T": '<span class="base-T">T</span>',
        "-": '<span class="base-gap">-</span>',
    }.get(base, f'<span class="base-other">{safe_base}</span>')


def _format_sequence(seq, color):
    if not color:
        return html.escape(seq)
    return "".join(_color_base(base) for base in seq)

def build_alignment_string(row, width=80, color=True, validate=True):
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