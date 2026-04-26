import heapq
import logging

import pandas as pd
from pandas.errors import EmptyDataError


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

MMSEQS_RESULT_COLUMNS = MMSEQS_REQUIRED_COLS + ["min_cov", "min_cov*pident"]


def _missing_required_columns(frame_columns):
    return [column for column in MMSEQS_REQUIRED_COLS if column not in frame_columns]


def _validate_mmseqs_frame(frame, source_name):
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
    min_cov = chunk[["qcov", "tcov"]].min(axis=1)
    scored_chunk = chunk.copy()
    scored_chunk["min_cov"] = min_cov
    scored_chunk["min_cov*pident"] = min_cov * scored_chunk["pident"]
    return scored_chunk


def _update_similarity_max(current_max, grouped_max):
    for sequence_id, similarity_value in grouped_max.items():
        similarity_value = float(similarity_value)
        if similarity_value > current_max.get(sequence_id, float("-inf")):
            current_max[sequence_id] = similarity_value


def _push_top_rows(top_rows_heap, rows, row_order, top_n):
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
    top_rows = [
        row_dict
        for _, _, row_dict in sorted(top_rows_heap, key=lambda item: (-item[0], item[1]))
    ]
    if top_rows:
        results_filt = pd.DataFrame(top_rows)
    else:
        results_filt = pd.DataFrame(columns=MMSEQS_RESULT_COLUMNS)
    if not results_filt.empty:
        results_filt = (
            results_filt.sort_values(by=["min_cov*pident"], ascending=False)
            .reset_index(drop=True)
        )
    return results_filt


def summarize_mmseqs_output(results_path, similarity_threshold, top_n=100, chunksize=100000):
    query_similarity_max = {}
    target_similarity_max = {}
    query_ids_with_hits = set()
    target_ids_with_hits = set()
    query_ids_above_threshold = set()
    target_ids_above_threshold = set()
    top_rows_heap = []
    row_order = 0
    total_hits = 0

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

        query_ids_with_hits.update(scored_chunk["query"].unique().tolist())
        target_ids_with_hits.update(scored_chunk["target"].unique().tolist())

        _update_similarity_max(
            query_similarity_max,
            scored_chunk.groupby("query")["min_cov*pident"].max(),
        )
        _update_similarity_max(
            target_similarity_max,
            scored_chunk.groupby("target")["min_cov*pident"].max(),
        )

        leaked = scored_chunk[scored_chunk["min_cov*pident"] >= similarity_threshold]
        if leaked.empty:
            continue

        query_ids_above_threshold.update(leaked["query"].unique().tolist())
        target_ids_above_threshold.update(leaked["target"].unique().tolist())

        leaked_top = leaked.nlargest(top_n, "min_cov*pident")
        row_order = _push_top_rows(top_rows_heap, leaked_top, row_order, top_n)

    results_filt = _finalize_results_frame(top_rows_heap)

    return {
        "query_similarity_max": query_similarity_max,
        "target_similarity_max": target_similarity_max,
        "query_ids_with_hits": query_ids_with_hits,
        "target_ids_with_hits": target_ids_with_hits,
        "query_ids_above_threshold": query_ids_above_threshold,
        "target_ids_above_threshold": target_ids_above_threshold,
        "results_filt": results_filt,
        "total_hits": total_hits,
    }