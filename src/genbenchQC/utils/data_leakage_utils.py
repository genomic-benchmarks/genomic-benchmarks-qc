from itertools import zip_longest

def get_basic_stats_from_aggregates(filename_train, train_stats, filename_test, test_stats):
    basic_stats = {
        "train_filename": str(filename_train),
        "test_filename": str(filename_test),
        "number_of_sequences_train": int(train_stats["count"]),
        "number_of_sequences_test": int(test_stats["count"]),
        "min_length_train": int(train_stats["min_length"]),
        "mean_length_train": float(train_stats["mean_length"]),
        "max_length_train": int(train_stats["max_length"]),
        "min_length_test": int(test_stats["min_length"]),
        "mean_length_test": float(test_stats["mean_length"]),
        "max_length_test": int(test_stats["max_length"]),
    }
    return basic_stats

def get_threshold_stats(results, results_filt, similarity_threshold, num_train_seqs, num_test_seqs):

    num_queries_with_hits = len(results['query'].unique())
    num_queries_above_thr = len(results_filt['query'].unique())
    perc_queries_above_thr = (num_queries_above_thr / num_test_seqs) * 100 if num_test_seqs > 0 else 0.0

    num_targets_with_hits = len(results['target'].unique())
    num_targets_above_thr = len(results_filt['target'].unique())
    perc_targets_above_thr = (num_targets_above_thr / num_train_seqs) * 100 if num_train_seqs > 0 else 0.0

    threshold_stats = {
        "similarity_threshold": similarity_threshold,
        "perc_queries_above_thr": perc_queries_above_thr,
        "num_queries_above_thr": num_queries_above_thr,
        "num_all_queries": num_test_seqs,
        "num_queries_without_hits": max(num_test_seqs - num_queries_with_hits, 0),
        "perc_targets_above_thr": perc_targets_above_thr,
        "num_targets_above_thr": num_targets_above_thr,
        "num_all_targets": num_train_seqs,
        "num_targets_without_hits": max(num_train_seqs - num_targets_with_hits, 0),
        "hits": len(results),
        "total_combinations": num_train_seqs * num_test_seqs,
    }
    return threshold_stats

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

def build_alignment_string(row, width=80, color=True, validate=True):
    def get_midline(qaln: str, taln: str) -> str:
        # '|' = match, '.' = mismatch, ' ' = gap
        return "".join(
            " " if (q == "-" or t == "-")
            else "|" if q == t
            else "."
            for q, t in zip(qaln, taln)
        )

    def validate_inputs(qstart, tstart, qend, tend, qseq, tseq, qaln, taln):
        # ---- Coordinate sanity ---
        def validate_coords(start, end, seqlen, label):
            if start < 1:
                raise ValueError(f"{label} start < 1 (1-based expected): {start}")
            if end < start:
                raise ValueError(f"{label} end ({end}) < start ({start})")
            if end > seqlen:
                raise ValueError(
                    f"{label} coordinates [{start}, {end}] exceed sequence length {seqlen}"
                )
            
        q_len = len(qseq)
        t_len = len(tseq)
        validate_coords(qstart, qend, q_len, "Query")
        validate_coords(tstart, tend, t_len, "Target")

        # ---- Alignment consistency ----
        if len(qaln) != len(taln):
            raise ValueError(
                f"Alignment length mismatch: "
                f"  len(qaln)={len(qaln)}, len(taln)={len(taln)}"
            )

        # ---- Ungapped alignment should match original slice ----
        qstart0 = qstart - 1
        tstart0 = tstart - 1

        expected_q = qseq[qstart0:qend]
        expected_t = tseq[tstart0:tend]

        ungapped_qaln = qaln.replace("-", "")
        ungapped_taln = taln.replace("-", "")

        if ungapped_qaln != expected_q:
            raise ValueError(
                "Query alignment mismatch:\n"
                "  Ungapped aligned sequence does not match original sequence slice."
            )

        if ungapped_taln != expected_t:
            raise ValueError(
                "Target alignment mismatch:\n"
                "  Ungapped aligned sequence does not match original sequence slice."
            )    
    
    # Get aligned sequences with appropriate padding to show the alignment in the correct position relative to the original sequences
    def get_alignment(qstart, tstart, qseq, tseq, qaln, taln, tend, qend):
        # 1-based → 0-based indices
        qstart0 = qstart - 1
        tstart0 = tstart - 1

        # initialize padding
        q_left = ""
        t_left = ""    

        delta = qstart0 - tstart0
        if delta > 0:
            t_left = " " * delta   # target shifted right
        elif delta < 0:
            q_left = " " * (-delta) # query shifted right

        mid_left = " " * (max(qstart0, tstart0))
        mid_line = mid_left + get_midline(qaln, taln)

        new_taln = tseq[:tstart0] + taln + tseq[tend:]
        new_qaln = qseq[:qstart0] + qaln + qseq[qend:]

        t_line = t_left + new_taln
        q_line = q_left + new_qaln
        
        return t_line, mid_line, q_line
    
    # Wrap string s into lines of length width
    def wrap(s, width): 
        return [s[i:i+width] for i in range(0, len(s), width)]
    
    # Define color mapping for bases
    def color_base(c):
        return {
            "A": '<span class="base-A">A</span>',
            "C": '<span class="base-C">C</span>',
            "G": '<span class="base-G">G</span>',
            "T": '<span class="base-T">T</span>',
            "-": '<span class="base-gap">-</span>',
        }.get(c, f'<span class="base-other">{c}</span>')

    # Format the alignment for display
    def format_seq(seq):
        return "".join(color_base(c) for c in seq) if color else seq

    qstart = int(row["qstart"])
    tstart = int(row["tstart"])
    qend   = int(row["qend"])
    tend   = int(row["tend"])

    qseq = row["qseq"]
    tseq = row["tseq"]
    qaln = row["qaln"]
    taln = row["taln"]

    if validate:
        validate_inputs(qstart, tstart, qend, tend, qseq, tseq, qaln, taln)

    t_line, mid_line, q_line = get_alignment(
        qstart, tstart, qseq, tseq, qaln, taln, tend, qend
    )    

    blocks = []
    for t, m, q in zip_longest(wrap(t_line, width), wrap(mid_line, width), wrap(q_line, width), fillvalue=""):
        blocks.append(f"T    {format_seq(t)}")
        blocks.append(f"     {m}")   # keep midline uncolored
        blocks.append(f"Q    {format_seq(q)}")
        blocks.append("")

    return "\n".join(blocks)