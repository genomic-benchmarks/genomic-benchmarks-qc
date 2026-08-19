"""The `gb-qc` command line interface.

Each command validates its options, reports the first problem as a message on
stderr with a non-zero exit code, and then hands over to the matching
`evaluate_*.run`. Nothing else lives here: the commands are thin wrappers, so
that the same evaluations can be called directly from Python.
"""

import logging
import re
from pathlib import Path
from typing import List, Optional

import typer

from genomic_benchmarks_qc.evaluate_classes import run as run_evaluate_classes
from genomic_benchmarks_qc.evaluate_splits import run as run_evaluate_splits

app = typer.Typer(no_args_is_help=True)

# Valid choices for validation.
# VALID_FORMATS serves two roles: it is the human-readable list of accepted file
# extensions shown in error messages, and its normalized members ('csv', 'tsv',
# 'fasta') are what inputs are actually validated against — see
# _normalize_format, which collapses the gzip and 'fa' variants before the check.
VALID_FORMATS = ['fasta', 'fasta.gz', 'fa', 'fa.gz', 'csv', 'csv.gz', 'tsv', 'tsv.gz']
VALID_REPORT_TYPES = ['json', 'html', 'simple']
VALID_PLOT_TYPES = ['boxen', 'violin']
VALID_LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

# MMseqs2 byte sizes: 0, or a non-zero integer with an optional unit suffix.
SPLIT_MEMORY_LIMIT_RE = re.compile(r'\A(?:0|[1-9][0-9]*[BKMGT]?)\Z')


def _fail(message: str):
    """Report a validation problem on stderr and exit non-zero."""
    typer.echo(f"Error: {message}", err=True)
    raise typer.Exit(code=1)


def _validate_input_files(files: List[str], description: str):
    """Check every input path exists, naming the kind of input in the message."""
    for file_path in files:
        if not Path(file_path).is_file():
            _fail(f"{description} does not exist: {file_path}")


def _validate_choice(value, valid_values: List[str], description: str):
    """Check one value is among the accepted ones."""
    if value not in valid_values:
        _fail(f"Invalid {description} '{value}'. Must be one of: {', '.join(valid_values)}")


def _validate_choices(values, valid_values: List[str], description: str):
    """Check every value of a repeatable option is among the accepted ones."""
    for value in values:
        _validate_choice(value, valid_values, description)


def _infer_format(file_path: str) -> str:
    """Infer format from file extension."""
    # Only the file name matters: a dot in a parent directory must not be
    # mistaken for an extension.
    name = Path(file_path).name.lower()
    if name.endswith('.gz'):
        base = name[:-len('.gz')]
        ext = base.rsplit('.', 1)[-1] if '.' in base else ''
        return f"{ext}.gz"
    return name.rsplit('.', 1)[-1] if '.' in name else ''


def _normalize_format(fmt: str) -> str:
    """Reduce a format to its base family, ignoring gzip and fa/fasta variants.

    Gzip is detected per-file downstream from the actual filename, and 'fa' is
    just an alias for 'fasta', so files differing only in those respects are
    considered the same format and may be mixed.
    """
    base = fmt[:-len('.gz')] if fmt.endswith('.gz') else fmt
    return 'fasta' if base == 'fa' else base


def _infer_format_from_inputs(files: List[str]) -> str:
    """Infer format across all input files, ensuring they agree."""
    formats = {_normalize_format(_infer_format(f)) for f in files}
    if len(formats) > 1:
        detected = ', '.join(f"{f} ({_infer_format(f)})" for f in files)
        _fail(f"All input files must have the same format, but got: {detected}")
    return formats.pop()


def _resolve_format(files: List[str]) -> str:
    """Infer the shared input format of all files and check it is supported."""
    format = _infer_format_from_inputs(files)
    _validate_choice(format, VALID_FORMATS, 'format')
    return format


def _run_command(command, *args, **kwargs):
    """Call an evaluation and turn any failure into a non-zero exit code.

    Only the message is echoed here, but the traceback is not thrown away with
    the exception: it goes to the log at DEBUG, so `--log-level DEBUG` records
    the cause even for failures raised before `log_failures` takes over - a
    scratch directory that cannot be created, say.
    """
    try:
        command(*args, **kwargs)
    except typer.Exit:
        raise
    except Exception as exc:
        # Failures inside the evaluation proper are logged twice at DEBUG, once
        # here and once by `log_failures`. Worth it: the alternative is guessing
        # whether the evaluation got far enough to have logged anything at all.
        logging.debug("Command failed.", exc_info=True)
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1)

@app.command()
def evaluate_classes(
    input: List[str] = typer.Option(..., help="Input file(s)."),
    sequence_column: List[str] = typer.Option(['sequence'], help="One or more sequence column names for CSV/TSV inputs."),
    label_column: str = typer.Option('label', help="Label column name for single-file CSV inputs."),
    label_list: List[str] = typer.Option(['infer'], help='List of labels to consider or "infer" to detect labels automatically.'),
    regression: bool = typer.Option(False, help="Treat label column as regression target and split into high/low."),
    out_folder: str = typer.Option('.', help="Output folder for reports."),
    report_types: List[str] = typer.Option(['html', 'simple'], help="Types of reports to generate (json, html, simple)."),
    end_position: Optional[int] = typer.Option(None, help="Last position included in per-position stats. Defaults to the last position at least 75% of each class's sequences reach."),
    plot_type: str = typer.Option('boxen', help="Plot type to use for visualizations (boxen, violin)."),
    log_level: str = typer.Option('INFO', help="Logging level."),
    log_file: Optional[str] = typer.Option(None, help="Optional path to write logs to."),
):
    """
    Evaluate sequence characteristics across different classes/labels in the dataset.
    """
    _validate_input_files(input, 'Input file')
    format = _resolve_format(input)

    # One FASTA file holds one class, so a single file has nothing to compare
    if format == 'fasta' and len(input) < 2:
        _fail("When format is 'fasta', at least 2 input files are required (one per class).")

    _validate_choices(report_types, VALID_REPORT_TYPES, 'report type')
    _validate_choice(plot_type, VALID_PLOT_TYPES, 'plot type')
    _validate_choice(log_level, VALID_LOG_LEVELS, 'log level')

    if end_position is not None and end_position <= 0:
        _fail(f"end_position must be a positive integer, got {end_position}")

    _run_command(
        run_evaluate_classes,
        input=input,
        format=format,
        out_folder=out_folder,
        sequence_column=sequence_column,
        label_column=label_column,
        label_list=label_list,
        regression=regression,
        report_types=report_types,
        end_position=end_position,
        plot_type=plot_type,
        log_level=log_level,
        log_file=log_file,
    )

@app.command()
def evaluate_splits(
    train_input: List[str] = typer.Option(..., help="Path to the dataset file(s) with training data."),
    test_input: List[str] = typer.Option(..., help="Path to the dataset file(s) with testing data."),
    sequence_column: List[str] = typer.Option(['sequence'], help="One or more sequence column names for CSV/TSV inputs."),
    out_folder: str = typer.Option('.', help="Output folder for reports."),
    report_types: List[str] = typer.Option(['html', 'simple'], help="Types of reports to generate (json, html, simple)."),
    similarity_threshold: float = typer.Option(90.0, help="Similarity threshold for data leakage detection (%)."),
    threads: Optional[int] = typer.Option(None, help="Set maximum number of threads MMseqs2 will use."),
    split_memory_limit: Optional[str] = typer.Option(None, "--split-memory-limit", help="Upper RAM limit for MMseqs2 prefilter structures (e.g., 10G, 1T)."),
    keep_tmp_files: bool = typer.Option(False, help="Keep temporary files for debugging."),
    log_level: str = typer.Option('INFO', help="Logging level."),
    log_file: Optional[str] = typer.Option(None, help="Optional path to write logs to."),
):
    """
    Evaluate data leakage in dataset train-test split.
    """
    _validate_input_files(train_input, 'Training input file')
    _validate_input_files(test_input, 'Test input file')

    # Both halves are searched against each other, so they must share a format
    format = _resolve_format(train_input + test_input)

    _validate_choices(report_types, VALID_REPORT_TYPES, 'report type')
    _validate_choice(log_level, VALID_LOG_LEVELS, 'log level')

    if not 0 <= similarity_threshold <= 100:
        _fail(f"similarity_threshold must be between 0 and 100, got {similarity_threshold}")

    if threads is not None and threads <= 0:
        _fail(f"threads must be a positive integer, got {threads}")

    # MMseqs2's byte-size grammar: 0, or a positive integer with an optional unit
    if split_memory_limit is not None:
        split_memory_limit = split_memory_limit.strip().upper()
        if not SPLIT_MEMORY_LIMIT_RE.match(split_memory_limit):
            _fail(
                f"split_memory_limit must be 0 or a positive integer optionally "
                f"followed by B, K, M, G or T (e.g., 10G, 1T), got {split_memory_limit}"
            )

    _run_command(
        run_evaluate_splits,
        train_files=train_input,
        test_files=test_input,
        format=format,
        out_folder=out_folder,
        sequence_column=sequence_column,
        report_types=report_types,
        similarity_threshold=similarity_threshold,
        threads=threads,
        split_memory_limit=split_memory_limit,
        keep_tmp_files=keep_tmp_files,
        log_level=log_level,
        log_file=log_file,
    )

def main():
    """Console-script entry point installed as `gb-qc`."""
    app()

if __name__ == "__main__":
    main()