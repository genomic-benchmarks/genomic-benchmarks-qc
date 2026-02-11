from itertools import zip_longest

def get_basic_stats(filename_train, train_sequences, filename_test, test_sequences):
    basic_stats = {
        "train_filename": str(filename_train),
        "test_filename": str(filename_test),
        "number_of_sequences_train": len(train_sequences),
        "number_of_sequences_test": len(test_sequences),
        "min_length_train": min(len(seq) for seq in train_sequences),
        "mean_length_train": sum(len(seq) for seq in train_sequences) // len(train_sequences),
        "max_length_train": max(len(seq) for seq in train_sequences),
        "min_length_test": min(len(seq) for seq in test_sequences),
        "mean_length_test": sum(len(seq) for seq in test_sequences) // len(test_sequences),
        "max_length_test": max(len(seq) for seq in test_sequences),
    }
    return basic_stats

def get_threshold_stats(results, results_filt, coverage_threshold, num_train_seqs, num_test_seqs):

    num_queries_above_thr = len(results_filt['query'].unique())
    perc_queries_above_thr = (num_queries_above_thr / num_test_seqs) * 100 # Assuming num_test_seqs contains unique test sequences

    num_targets_above_thr = len(results_filt['target'].unique())
    perc_targets_above_thr = (num_targets_above_thr / num_train_seqs) * 100 # Assuming num_train_seqs contains unique train sequences

    threshold_stats = {
        "coverage_threshold": coverage_threshold,
        "perc_queries_above_thr": perc_queries_above_thr,
        "perc_targets_above_thr": perc_targets_above_thr,
        "hits": len(results),
        "total_combinations": num_train_seqs * num_test_seqs,
    }
    return threshold_stats

# Filter fasta file, keeping only sequences with IDs in ids_to_keep (for hits)
def filter_fasta_by_ids(fasta_path, new_fasta_path,ids_to_keep):
    with fasta_path.open() as file_in, new_fasta_path.open("w") as file_out:
        write = False

        for line in file_in:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                write = seq_id in ids_to_keep

            if write:
                file_out.write(line)

def build_alignment_string(row, width=80, color=True):
    def get_midline(qaln: str, taln: str) -> str:
        # '|' = match, '.' = mismatch, ' ' = gap
        return "".join(
            " " if (q == "-" or t == "-")
            else "|" if q == t
            else "."
            for q, t in zip(qaln, taln)
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

    t_line, mid_line, q_line = get_alignment(
        qstart=int(row["qstart"]),
        tstart=int(row["tstart"]),
        qseq=row["qseq"],
        tseq=row["tseq"],
        qaln=row["qaln"],
        taln=row["taln"],
        tend=int(row["tend"]),
        qend=int(row["qend"]),
    )

    blocks = []
    for t, m, q in zip_longest(wrap(t_line, width), wrap(mid_line, width), wrap(q_line, width), fillvalue=""):
        blocks.append(f"T    {format_seq(t)}")
        blocks.append(f"     {m}")   # keep midline uncolored
        blocks.append(f"Q    {format_seq(q)}")
        blocks.append("")

    return "\n".join(blocks)

# Add alignment strings to the results DataFrame
def add_alignments_to_results(results):
    results["alignment"] = results.apply(lambda row: build_alignment_string(row, color=False), axis=1)
    return results