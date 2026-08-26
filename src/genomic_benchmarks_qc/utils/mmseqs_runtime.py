"""Running the MMseqs2 search as a subprocess.

MMseqs2 is an external binary, so this module owns everything about talking to
it: checking it can run here at all, building the `easy-search` command line,
and forwarding its output into the logger instead of the terminal.
"""

import logging
import os
import platform
import shutil
import subprocess
import threading
import time
from pathlib import Path

from genomic_benchmarks_qc.utils.mmseqs_summary import MMSEQS_REQUIRED_COLS

logger = logging.getLogger(__name__)

SUPPORTED_CPU_FLAGS = ("avx2", "sse4_1", "sse2")
MMSEQS_PROGRESS_LOG_MIN_INTERVAL_SEC = 1.0


def _read_linux_cpu_flags():
    """Return the CPU feature flags from /proc/cpuinfo, or None if unreadable."""
    try:
        with open("/proc/cpuinfo", encoding="utf-8") as handle:
            for line in handle:
                if line.lower().startswith("flags"):
                    _, flags_str = line.split(":", 1)
                    return set(flags_str.strip().split())
    except FileNotFoundError:
        return None
    return None


def check_mmseqs_preflight():
    """Verify MMseqs2 can run on this machine, raising RuntimeError if not.

    The published binaries are compiled per instruction set, and one built for a
    newer CPU dies with SIGILL rather than a readable error, so the CPU flags
    are checked up front. Only x86_64 Linux can be checked this way; elsewhere
    the check is skipped with a warning.

    Skipped, not refused. The flags this looks for are x86_64 ones, so on ARM
    there is nothing here to check - but MMSeqs2 runs on ARM, and a check that
    exists to turn a SIGILL into a sentence must not become the reason the search
    never starts. On an architecture it cannot vouch for, letting MMSeqs2 speak
    for itself is the better answer.
    """
    mmseqs_path = shutil.which("mmseqs")
    if mmseqs_path is None:
        raise RuntimeError(
            "MMSeqs2 executable not found in PATH. "
            "Please install MMSeqs2 and ensure it is available in your environment."
        )
    logger.debug("Found MMSeqs2 at: %s", mmseqs_path)

    system = platform.system()
    if system != "Linux":
        logger.warning(
            "Skipping CPU feature checks for non-Linux system (%s). "
            "Ensure your MMSeqs2 binary is compatible with this platform.",
            system,
        )
        return

    arch = platform.machine().lower()
    if arch not in ("x86_64", "amd64"):
        logger.warning(
            "Skipping CPU feature checks on %s: %s are x86_64 instruction sets and "
            "this machine is not x86_64. MMSeqs2 ships builds for other "
            "architectures, so this is not a reason to stop - if the binary here is "
            "the wrong one, MMSeqs2 says so itself.",
            arch,
            ", ".join(SUPPORTED_CPU_FLAGS),
        )
        return

    flags = _read_linux_cpu_flags()
    if flags is None:
        raise RuntimeError(
            "Unable to read CPU flags from /proc/cpuinfo to verify MMSeqs2 support."
        )
    matched_flags = [flag for flag in SUPPORTED_CPU_FLAGS if flag in flags]
    if not matched_flags:
        raise RuntimeError(
            "CPU does not support any of the MMSeqs2-supported instruction set flags: "
            + ", ".join(SUPPORTED_CPU_FLAGS)
        )
    logger.debug(
        "CPU feature checks passed for MMSeqs2 using supported flags: %s",
        ", ".join(matched_flags),
    )


def run_search(
    test_fasta_file,
    train_fasta_file,
    out_file,
    tmp_dir,
    threads: int | None = None,
    split_memory_limit: str | None = None,
):
    """Search the test sequences against the train sequences, returning the hit table.

    Runs `mmseqs easy-search` with the test half as the query and the train half
    as the database, restricted to nucleotide search on the forward strand, and
    asks for `MMSEQS_REQUIRED_COLS` as a tab-separated table with a header
    (`--format-mode 4`) so it can be read back in chunks.

    Returns the path of that table. Raises RuntimeError if MMseqs2 is unusable
    or exits non-zero.
    """
    logger.info(
        "Running MMSeqs2, an ultrafast and sensitive search, for test sequences "
        "(query) against train sequences (db)."
    )

    check_mmseqs_preflight()

    cmd = [
        "mmseqs",
        "easy-search",
        str(test_fasta_file),
        str(train_fasta_file),
        str(out_file),
        str(tmp_dir),
        "--format-output",
        ",".join(MMSEQS_REQUIRED_COLS),
        "--format-mode",
        "4",
        "--search-type",
        "3",
        "--strand",
        "1",
    ]

    if threads is not None:
        cmd.extend(["--threads", str(threads)])
    if split_memory_limit is not None:
        cmd.extend(["--split-memory-limit", split_memory_limit])

    logger.debug("Running command: %s", " ".join(cmd))
    mmseqs_env = os.environ.copy()
    mmseqs_env["TTY"] = "1"
    logger.debug("Running MMSeqs2 with TTY=%s", mmseqs_env["TTY"])

    process = None
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=False,
            bufsize=0,
            env=mmseqs_env,
        )

        def _forward_stream(stream, stream_name):
            """Log one of the subprocess streams line by line as it arrives.

            MMseqs2 draws a progress bar by rewriting the same line with '\r',
            which would otherwise flood the log, so those lines are logged at
            most once per `MMSEQS_PROGRESS_LOG_MIN_INTERVAL_SEC`.
            """
            if stream is None:
                return
            last_progress_log_ts = None

            try:
                while True:
                    chunk = stream.read(4096)
                    if not chunk:
                        break

                    text = chunk.decode("utf-8", errors="replace")
                    for line in text.split("\n"):
                        if not line:
                            continue

                        # Check if this is a progress line (ends with \r before stripping)
                        is_progress = line.endswith("\r")
                        line = line.rstrip("\r").strip()
                        if not line:
                            continue

                        # Progress lines get throttled; regular lines always logged
                        if is_progress:
                            now = time.monotonic()
                            due = (
                                last_progress_log_ts is None
                                or (now - last_progress_log_ts)
                                >= MMSEQS_PROGRESS_LOG_MIN_INTERVAL_SEC
                            )
                            if due:
                                logger.debug("MMSeqs2 progress: %s", line)
                                last_progress_log_ts = now
                        else:
                            log_fn = logger.error if stream_name == "stderr" else logger.debug
                            log_fn("MMSeqs2 %s: %s", stream_name, line)
            finally:
                stream.close()

        stdout_thread = threading.Thread(
            target=_forward_stream,
            args=(process.stdout, "stdout"),
            daemon=True,
        )
        stderr_thread = threading.Thread(
            target=_forward_stream,
            args=(process.stderr, "stderr"),
            daemon=True,
        )
        stdout_thread.start()
        stderr_thread.start()

        return_code = process.wait()
        stdout_thread.join()
        stderr_thread.join()

        if return_code != 0:
            logger.error("MMSeqs2 search failed with return code: %s", return_code)
            raise RuntimeError("MMSeqs2 search failed.")

    except BaseException:
        # BaseException rather than Exception. A Ctrl-C or a cancelled task is
        # when a search left running matters most - it is the whole machine for
        # as long as it takes - and neither of those is an Exception. In a
        # terminal the signal reaches MMseqs2 anyway, because it is in the same
        # process group; called from a program, nothing else would stop it.
        if process is not None and process.poll() is None:
            process.kill()
            process.wait()
        raise

    logger.debug("MMSeqs2 easy-search completed.")

    return Path(out_file)
