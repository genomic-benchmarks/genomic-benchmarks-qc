"""Summarising an MMseqs2 hit table without holding it in memory.

An all-vs-all search of a large split produces far more hits than the report
needs, so the table is read in chunks and reduced on the way past: the per-
sequence maximum similarity, the sets of sequences with hits and above the
threshold, and the top hits for the alignment view.

Similarity is `min(qcov, tcov) * pident` - the coverage of the shorter of the
two sequences scaled by how identical the aligned part is - so that a short
exact match inside a long sequence does not count as a leak.

Sequences are carried through the search by number - `seq_<n>_<half>`, written
by the staging in `evaluate_splits` and read back here - so that joining a hit
to a sequence is an array index. Keyed by the id string instead, the two
per-sequence maxima were 240 MB at 300,000 sequences a half, which is the same
scale as the search itself and for two floats per sequence.

The hits above the threshold are also written out as they go past, if the caller
asks for it. That is a side effect in a module otherwise made of returns, and it
is here for the same reason the rest of it is: the export holds every leaked hit,
there can be far more of them than the report lists, and collecting them to write
at the end would give back the bound this module exists to keep.
"""

import heapq
import logging
import re

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError

# How a staged sequence is named in the FASTA both halves go into. MMseqs2
# carries only the name through the search, so this is the whole of the join
# between a hit and the sequence it is about.
SEQUENCE_ID_PATTERN = re.compile(r"seq_(\d+)_(?:train|test)\Z")


def sequence_id(index, half):
    """The FASTA id the sequence at `index` of `half` is staged under."""
    return f"seq_{index}_{half}"


def staged_ids(mask, half):
    """The FASTA ids of the sequences a boolean mask selects.

    The counterpart of `sequence_id`, for the one thing that still wants ids
    rather than positions: pulling records back out of a FASTA by name. It is
    built from the mask at the point of use rather than carried alongside it,
    because the mask is a bit per sequence and the ids are a string each.
    """
    return {sequence_id(int(index), half) for index in np.flatnonzero(mask)}


def _staged_indices(ids, count, source_name):
    """Turn staged ids into positions, refusing anything this did not stage."""
    extracted = pd.Index(ids).str.extract(SEQUENCE_ID_PATTERN, expand=False)
    if extracted.isna().any():
        unknown = pd.Index(ids)[extracted.isna()][0]
        raise RuntimeError(
            f"MMSeqs2 output in {source_name} names a sequence this run did not "
            f"stage: {unknown!r}. The hit table has to come from the FASTA files "
            f"staged alongside it."
        )
    indices = extracted.to_numpy(dtype=np.int64)
    if indices.size and (indices.max() >= count or indices.min() < 0):
        raise RuntimeError(
            f"MMSeqs2 output in {source_name} names sequence number "
            f"{int(indices.max())} of {count} staged."
        )
    return indices


MMSEQS_REQUIRED_COLS = [
    "query",
    "target",
    "qcov",
    "tcov",
    "pident",
    "evalue",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "alnlen",
    "qaln",
    "taln",
]

MMSEQS_DERIVED_COLS = ["min_cov", "min_cov*pident"]
MMSEQS_RESULT_COLUMNS = MMSEQS_REQUIRED_COLS + MMSEQS_DERIVED_COLS


def _missing_required_columns(frame_columns):
    """Return the required MMseqs2 columns absent from the given columns."""
    return [column for column in MMSEQS_REQUIRED_COLS if column not in frame_columns]


def _validate_mmseqs_frame(frame, source_name):
    """Check a chunk has every required column and drop rows with gaps in them.

    Raises RuntimeError on missing columns, which means the search was run with
    a different `--format-output` than this code expects.
    """
    missing_cols = _missing_required_columns(frame.columns)
    if missing_cols:
        raise RuntimeError(
            f"MMSeqs2 output is missing required columns in {source_name}: "
            + ", ".join(missing_cols)
        )

    invalid_rows = frame[MMSEQS_REQUIRED_COLS].isna().any(axis=1).sum()
    if invalid_rows > 0:
        logging.debug(
            "Found %s rows with missing values in MMSeqs2 output from %s. Removing them.",
            invalid_rows,
            source_name,
        )
        frame = frame.dropna(subset=MMSEQS_REQUIRED_COLS)

    return frame


def _score_mmseqs_chunk(chunk):
    """Add the coverage and similarity columns this module scores hits by."""
    min_cov = chunk[["qcov", "tcov"]].min(axis=1)
    scored_chunk = chunk.copy()
    scored_chunk["min_cov"] = min_cov
    scored_chunk["min_cov*pident"] = min_cov * scored_chunk["pident"]
    return scored_chunk


def _update_similarity_max(current_max, grouped_max, source_name):
    """Merge a chunk's per-sequence maxima into the running maxima, in place.

    Hits for one sequence can be spread over several chunks, so the maximum has
    to be carried across them rather than computed per chunk. `fmax` rather than
    `maximum` because a sequence with no hit yet holds NaN, which is what tells
    it apart from one whose best hit scored zero.
    """
    indices = _staged_indices(grouped_max.index, current_max.size, source_name)
    np.fmax.at(current_max, indices, grouped_max.to_numpy(dtype=np.float32))


def _push_top_rows(top_rows_heap, rows, row_order, top_n):
    """Keep the `top_n` most similar hits seen so far in a min-heap.

    The heap holds (similarity, arrival order, row), so the weakest hit is
    always the one dropped and equally similar hits keep the order they were
    read in. Returns the updated arrival counter.
    """
    for _, row in rows.iterrows():
        row_dict = row.to_dict()
        heapq.heappush(
            top_rows_heap,
            (float(row_dict["min_cov*pident"]), row_order, row_dict),
        )
        row_order += 1
        if len(top_rows_heap) > top_n:
            heapq.heappop(top_rows_heap)

    return row_order


def _finalize_results_frame(top_rows_heap):
    """Turn the heap of top hits into a frame sorted most-similar first.

    The sort key is (similarity descending, arrival order), which is a total
    order over the rows: two equally similar hits are separated by the counter
    `_push_top_rows` stamped them with, so the same input always produces the
    same table. That is the whole ordering - the frame is built in this order
    and left in it. Sorting the frame again afterwards was where the
    determinism went: `sort_values` defaults to quicksort, which is not stable,
    so it was free to shuffle the ties this key had just settled.
    """
    top_rows = [
        row_dict
        for _, _, row_dict in sorted(top_rows_heap, key=lambda item: (-item[0], item[1]))
    ]
    if top_rows:
        results_filt = pd.DataFrame(top_rows)
    else:
        results_filt = pd.DataFrame(columns=MMSEQS_RESULT_COLUMNS)
    return results_filt.reindex(columns=MMSEQS_RESULT_COLUMNS)


def build_mmseqs_export_frame(results_filt):
    """Return the hits as a frame with every export column, in a fixed order.

    Missing columns are filled with NA so the exported TSV has the same shape
    whether or not anything was found.
    """
    if results_filt is None or results_filt.empty:
        return pd.DataFrame(columns=MMSEQS_RESULT_COLUMNS)

    export_frame = results_filt.copy()

    for column in MMSEQS_RESULT_COLUMNS:
        if column not in export_frame.columns:
            export_frame[column] = pd.NA

    return export_frame.loc[:, MMSEQS_RESULT_COLUMNS]


def _append_leaked_hits(export_path, leaked, first_write):
    """Append one chunk's above-threshold hits to the export TSV."""
    build_mmseqs_export_frame(leaked).to_csv(
        export_path,
        sep="\t",
        index=False,
        mode="w" if first_write else "a",
        header=first_write,
    )


def summarize_mmseqs_output(results_path, similarity_threshold, query_count, target_count,
                            top_n=100, chunksize=100000, export_path=None):
    """Reduce an MMseqs2 hit table to the values the split report needs.

    Reads the table in chunks of `chunksize` rows, so memory use is bounded by
    the chunk size and by four bytes and a bit per sequence, rather than by the
    number of hits. An empty table is not an error - it means nothing matched.

    Args:
        results_path: The hit table MMseqs2 wrote.
        similarity_threshold: Percentage at or above which a hit counts as a leak.
        query_count: Sequences staged as queries, which is the test half.
        target_count: Sequences staged as targets, which is the train half.
        top_n: How many of the leaked hits to keep for the alignment view.
        chunksize: Rows read at a time.
        export_path: Where to write every leaked hit, if anywhere.

    Returns a dict with, for each of queries (test) and targets (train), the
    per-sequence maximum similarity as an array in staging order - NaN for a
    sequence the search returned nothing for - and a boolean array of which are
    at or above `similarity_threshold`; plus `total_hits`, every alignment the
    search found; `leaked_hits`, the ones at or above the threshold; and
    `results_filt`, the `top_n` most similar of those for the alignment view.

    `total_hits` and `leaked_hits` are easy to reach for interchangeably and are
    not the same number: a search reports far more alignments than it finds
    leaks, and it is the leaks the report counts.

    When `export_path` is given, every leaked hit is written there as a TSV, in
    the columns of `MMSEQS_RESULT_COLUMNS`. The file is created even when nothing
    leaked, so a clean split leaves a header rather than a missing file.
    """
    # NaN, not zero: a sequence the search said nothing about has to stay
    # distinguishable from one whose best alignment scored nothing.
    query_similarity_max = np.full(query_count, np.nan, dtype=np.float32)
    target_similarity_max = np.full(target_count, np.nan, dtype=np.float32)
    query_above_threshold = np.zeros(query_count, dtype=bool)
    target_above_threshold = np.zeros(target_count, dtype=bool)
    top_rows_heap = []
    row_order = 0
    total_hits = 0
    leaked_hits = 0
    exported_any = False

    try:
        chunk_iter = pd.read_csv(results_path, sep="\t", chunksize=chunksize)
    except EmptyDataError:
        chunk_iter = []

    for chunk in chunk_iter:
        chunk = _validate_mmseqs_frame(chunk, results_path)
        if chunk.empty:
            continue

        scored_chunk = _score_mmseqs_chunk(chunk)
        total_hits += len(scored_chunk)

        # Which sequences had a hit at all falls out of the maxima: the entries
        # that are no longer NaN.
        _update_similarity_max(
            query_similarity_max,
            scored_chunk.groupby("query")["min_cov*pident"].max(),
            results_path,
        )
        _update_similarity_max(
            target_similarity_max,
            scored_chunk.groupby("target")["min_cov*pident"].max(),
            results_path,
        )

        leaked = scored_chunk[scored_chunk["min_cov*pident"] >= similarity_threshold]
        if leaked.empty:
            continue

        leaked_hits += len(leaked)
        if export_path is not None:
            _append_leaked_hits(export_path, leaked, first_write=not exported_any)
            exported_any = True

        query_above_threshold[
            _staged_indices(pd.Index(leaked["query"].unique()), query_count, results_path)
        ] = True
        target_above_threshold[
            _staged_indices(pd.Index(leaked["target"].unique()), target_count, results_path)
        ] = True

        leaked_top = leaked.nlargest(top_n, "min_cov*pident")
        row_order = _push_top_rows(top_rows_heap, leaked_top, row_order, top_n)

    results_filt = _finalize_results_frame(top_rows_heap)

    if export_path is not None and not exported_any:
        _append_leaked_hits(export_path, None, first_write=True)

    return {
        "query_similarity_max": query_similarity_max,
        "target_similarity_max": target_similarity_max,
        "query_above_threshold": query_above_threshold,
        "target_above_threshold": target_above_threshold,
        "results_filt": results_filt,
        "total_hits": total_hits,
        "leaked_hits": leaked_hits,
    }
