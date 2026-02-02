import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
from itertools import zip_longest

def plot_coverage_histograms(results, coverage_threshold, bins=50):
    fig, ax = plt.subplots()
    ax.hist(results["qcov"], bins=bins, alpha=0.5, label="Query (Test)")
    ax.hist(results["tcov"], bins=bins, alpha=0.5, label="Target (Train)")
    ax.axvline(coverage_threshold, linestyle="--", linewidth=1.2, color="red", label="Threshold")
    ax.set_yscale("log")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Count")
    ax.legend()
    return fig

def build_alignment_string(row, width=80):
    def get_midline(qaln: str, taln: str) -> str:
        # '|' = match, '.' = mismatch, ' ' = gap
        return "".join(
            " " if (q == "-" or t == "-")
            else "|" if q == t
            else "."
            for q, t in zip(qaln, taln)
        )
    
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

        t_line = t_left + "".join(new_taln)
        q_line = q_left + "".join(new_qaln)
        
        return t_line, mid_line, q_line

    def wrap(s, width):
        return [s[i:i+width] for i in range(0, len(s), width)]
    
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
        blocks.append(t)
        blocks.append(m)
        blocks.append(q)
        blocks.append("")

    return "\n".join(blocks)

def add_alignments_to_results(results):
    results["alignment_str"] = results.apply(build_alignment_string, axis=1)
    return results