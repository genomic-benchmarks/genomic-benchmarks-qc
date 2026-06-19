import logging
import os
import platform
import shutil
import subprocess
import threading
import time
from pathlib import Path
from typing import Optional

from genbenchQC.utils.mmseqs_summary import MMSEQS_REQUIRED_COLS


SUPPORTED_CPU_FLAGS = ("avx2", "sse4_1", "sse2")
MMSEQS_PROGRESS_LOG_MIN_INTERVAL_SEC = 1.0


def _read_linux_cpu_flags():
    try:
        with open("/proc/cpuinfo", "r", encoding="utf-8") as handle:
            for line in handle:
                if line.lower().startswith("flags"):
                    _, flags_str = line.split(":", 1)
                    return set(flags_str.strip().split())
    except FileNotFoundError:
        return None
    return None


def check_mmseqs_preflight():
    mmseqs_path = shutil.which("mmseqs")
    if mmseqs_path is None:
        raise RuntimeError(
            "MMSeqs2 executable not found in PATH. "
            "Please install MMSeqs2 and ensure it is available in your environment."
        )
    logging.debug("Found MMSeqs2 at: %s", mmseqs_path)

    system = platform.system()
    if system != "Linux":
        logging.warning(
            "Skipping CPU feature checks for non-Linux system (%s). "
            "Ensure your MMSeqs2 binary is compatible with this platform.",
            system,
        )
        return

    arch = platform.machine().lower()
    if arch not in ("x86_64", "amd64"):
        raise RuntimeError(
            f"Unsupported architecture for MMSeqs2 preflight checks: {arch}. "
            "Expected x86_64."
        )

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
    logging.debug(
        "CPU feature checks passed for MMSeqs2 using supported flags: %s",
        ", ".join(matched_flags),
    )


def run_search(
    test_fasta_file,
    train_fasta_file,
    out_file,
    tmp_dir,
    threads: Optional[int] = None,
    split_memory_limit: Optional[str] = None,
):
    logging.info(
        "Running MMSeqs2, an ultrafast and sensitive search, for test sequences (query) against train sequences (db)."
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

    logging.debug("Running command: %s", " ".join(cmd))
    mmseqs_env = os.environ.copy()
    mmseqs_env["TTY"] = "1"
    logging.debug("Running MMSeqs2 with TTY=%s", mmseqs_env["TTY"])

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
                            if last_progress_log_ts is None or (now - last_progress_log_ts) >= MMSEQS_PROGRESS_LOG_MIN_INTERVAL_SEC:
                                logging.debug("MMSeqs2 progress: %s", line)
                                last_progress_log_ts = now
                        else:
                            log_fn = logging.error if stream_name == "stderr" else logging.debug
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
            logging.error("MMSeqs2 search failed with return code: %s", return_code)
            raise RuntimeError("MMSeqs2 search failed.")

    except Exception:
        if "process" in locals() and process.poll() is None:
            process.kill()
            process.wait()
        raise

    logging.debug("MMSeqs2 easy-search completed.")

    return Path(out_file)