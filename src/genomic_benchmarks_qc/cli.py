"""The `gb-qc` command line interface.

Each command validates its options, reports the first problem as a message on
stderr with a non-zero exit code, and then hands over to the matching
`evaluate_*.run`. Nothing else lives here: the commands are thin wrappers, so
that the same evaluations can be called directly from Python.
"""

import logging
import re
import traceback
from pathlib import Path

import typer

from genomic_benchmarks_qc import __version__
from genomic_benchmarks_qc.evaluate_classes import REPORT_TYPES as CLASSES_REPORT_TYPES
from genomic_benchmarks_qc.evaluate_classes import run as run_evaluate_classes
from genomic_benchmarks_qc.evaluate_splits import REPORT_TYPES as SPLITS_REPORT_TYPES
from genomic_benchmarks_qc.evaluate_splits import run as run_evaluate_splits
from genomic_benchmarks_qc.utils.input_utils import PACKAGE_LOGGER_NAME
from genomic_benchmarks_qc.utils.seq_stats import (
    DEFAULT_MIN_COVERAGE,
    MIN_SEQUENCES_PER_REPORTED_POSITION,
)
from genomic_benchmarks_qc.utils.testing import MIN_SEQUENCES_PER_CLASS

logger = logging.getLogger(__name__)

app = typer.Typer(no_args_is_help=True)

# Valid choices for validation.
# VALID_FORMATS serves two roles: it is the human-readable list of accepted file
# extensions shown in error messages, and its normalized members ('csv', 'tsv',
# 'fasta') are what inputs are actually validated against — see
# _normalize_format, which collapses the gzip and 'fa' variants before the check.
VALID_FORMATS = ['fasta', 'fasta.gz', 'fa', 'fa.gz', 'csv', 'csv.gz', 'tsv', 'tsv.gz']
VALID_PLOT_TYPES = ['boxen', 'violin']

# Each command's --report-types help. Built from the lists the pipelines declare,
# rather than spelled out here: a command that cannot write a report type must
# not be able to offer it.
CLASSES_REPORT_TYPES_HELP = f"Types of reports to generate ({', '.join(CLASSES_REPORT_TYPES)})."
SPLITS_REPORT_TYPES_HELP = f"Types of reports to generate ({', '.join(SPLITS_REPORT_TYPES)})."
VALID_LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

# MMseqs2 byte sizes: 0, or a non-zero integer with an optional unit suffix.
SPLIT_MEMORY_LIMIT_RE = re.compile(r'\A(?:0|[1-9][0-9]*[BKMGT]?)\Z')


def _report_version(value: bool):
    """Print the version and stop, for the `--version` option below."""
    if value:
        typer.echo(f"gb-qc {__version__}")
        raise typer.Exit()


def _fail(message: str):
    """Report a validation problem on stderr and exit non-zero."""
    typer.echo(f"Error: {message}", err=True)
    raise typer.Exit(code=1)


def _validate_input_files(files: list[str], description: str):
    """Check every input path exists, naming the kind of input in the message."""
    for file_path in files:
        if not Path(file_path).is_file():
            _fail(f"{description} does not exist: {file_path}")


def _validate_choice(value, valid_values: list[str], description: str):
    """Check one value is among the accepted ones."""
    if value not in valid_values:
        _fail(f"Invalid {description} '{value}'. Must be one of: {', '.join(valid_values)}")


def _validate_choices(values, valid_values: list[str], description: str):
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


def _infer_format_from_inputs(files: list[str]) -> str:
    """Infer format across all input files, ensuring they agree."""
    formats = {_normalize_format(_infer_format(f)) for f in files}
    if len(formats) > 1:
        detected = ', '.join(f"{f} ({_infer_format(f)})" for f in files)
        _fail(f"All input files must have the same format, but got: {detected}")
    return formats.pop()


def _resolve_format(files: list[str]) -> str:
    """Infer the shared input format of all files and check it is supported."""
    format = _infer_format_from_inputs(files)
    _validate_choice(format, VALID_FORMATS, 'format')
    return format


def _run_command(command, *args, **kwargs):
    """Call an evaluation and turn any failure into a non-zero exit code.

    Only the message is echoed here, but the traceback is not thrown away with
    the exception: it goes to the log at DEBUG, so `--log-level DEBUG` records
    the cause even for failures raised before `log_failures` takes over - a
    scratch directory that cannot be created, say. Setting the log up is itself
    one of those failures, and the one case where the log cannot record it, so
    the traceback falls back to stderr there.
    """
    try:
        command(*args, **kwargs)
    except typer.Exit:
        raise
    except Exception as exc:
        # `setup_logger` runs inside the command, so it can be what failed - an
        # unwritable `--log-file`, say - and then nothing ever configured the
        # package logger and the DEBUG record below has nowhere to go.
        debug_lost = (str(kwargs.get('log_level', '')).upper() == 'DEBUG'
                      and not logging.getLogger(PACKAGE_LOGGER_NAME).isEnabledFor(logging.DEBUG))
        # Failures inside the evaluation proper are logged twice at DEBUG, once
        # here and once by `log_failures`. Worth it: the alternative is guessing
        # whether the evaluation got far enough to have logged anything at all.
        logger.debug("Command failed.", exc_info=True)
        if debug_lost:
            # Asking for DEBUG has to produce a traceback somewhere, or the one
            # failure that breaks logging is the one that explains itself least.
            typer.echo(traceback.format_exc(), err=True)
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

# Nothing to configure at this level, but a Typer app without a callback has no
# place to hang an option that belongs to no command, and `--version` is the one
# thing a user asks the tool before asking it anything else - most often to say
# which version produced a report. Eager, so it answers on its own rather than
# needing a command after it.
@app.callback()
def gb_qc(
    version: bool = typer.Option(
        False, '--version', callback=_report_version, is_eager=True,
        help="Show the version and exit."),
):
    """Quality control for genomic machine learning datasets."""


@app.command()
def evaluate_classes(
    input: list[str] = typer.Option(..., help="Input file(s)."),
    sequence_column: list[str] = typer.Option(
        ['sequence'], help="One or more sequence column names for CSV/TSV inputs."),
    label_column: str = typer.Option('label', help="Label column name for single-file CSV inputs."),
    label_list: list[str] = typer.Option(
        ['infer'], help='List of labels to consider or "infer" to detect labels automatically.'),
    regression: bool = typer.Option(
        False, help="Treat label column as regression target and split into high/low."),
    out_folder: str = typer.Option('.', help="Output folder for reports."),
    report_types: list[str] = typer.Option(
        ['html', 'simple'],
        help=CLASSES_REPORT_TYPES_HELP),
    end_position: int | None = typer.Option(
        None, help=(
            "Last position the per-position checks reach. Defaults to the last position at "
            f"least {MIN_SEQUENCES_PER_REPORTED_POSITION} of each class's sequences reach. Can "
            "only narrow the flagged window, never widen it - "
            "what is flagged is decided by --min-coverage, and that window is what the figures "
            "draw."
        )),
    min_coverage: float = typer.Option(
        DEFAULT_MIN_COVERAGE, help=(
            "Fraction of each class's sequences that must reach a position before it can be "
            f"flagged, on top of the {MIN_SEQUENCES_PER_CLASS} sequences every compared position "
            "needs. This window is what the per-position figures draw; further positions are "
            f"reported as Unknown and not drawn. 0 leaves only the {MIN_SEQUENCES_PER_CLASS}."
        )),
    plot_type: str = typer.Option(
        'boxen', help="Plot type to use for visualizations (boxen, violin)."),
    log_level: str = typer.Option('INFO', help="Logging level."),
    log_file: str | None = typer.Option(None, help="Optional path to write logs to."),
):
    """Compare every pair of classes and flag what tells them apart.

    Sequence length, base and dinucleotide composition, per-position composition
    and duplicate sequences, each flagged Pass, Warning or Fail. Reports land in
    out-folder/class/, one directory per compared pair.
    """
    _validate_input_files(input, 'Input file')
    format = _resolve_format(input)

    # One FASTA file holds one class, so a single file has nothing to compare
    if format == 'fasta' and len(input) < 2:
        _fail("When format is 'fasta', at least 2 input files are required (one per class).")

    _validate_choices(report_types, CLASSES_REPORT_TYPES, 'report type')
    _validate_choice(plot_type, VALID_PLOT_TYPES, 'plot type')
    _validate_choice(log_level, VALID_LOG_LEVELS, 'log level')

    if end_position is not None and end_position <= 0:
        _fail(f"end_position must be a positive integer, got {end_position}")

    if not 0 <= min_coverage <= 1:
        _fail(f"min_coverage must be a fraction between 0 and 1, got {min_coverage}")

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
        min_coverage=min_coverage,
        plot_type=plot_type,
        log_level=log_level,
        log_file=log_file,
    )

@app.command()
def evaluate_splits(
    train_input: list[str] = typer.Option(
        ..., help="Path to the dataset file(s) with training data."),
    test_input: list[str] = typer.Option(
        ..., help="Path to the dataset file(s) with testing data."),
    sequence_column: list[str] = typer.Option(
        ['sequence'], help="One or more sequence column names for CSV/TSV inputs."),
    out_folder: str = typer.Option('.', help="Output folder for reports."),
    report_types: list[str] = typer.Option(
        ['html', 'simple'],
        help=SPLITS_REPORT_TYPES_HELP),
    similarity_threshold: float = typer.Option(
        90.0, help="Similarity threshold for data leakage detection (%)."),
    threads: int | None = typer.Option(
        None, help="Set maximum number of threads MMseqs2 will use."),
    split_memory_limit: str | None = typer.Option(
        None, "--split-memory-limit",
        help="Upper RAM limit for MMseqs2 prefilter structures (e.g., 10G, 1T)."),
    keep_tmp_files: bool = typer.Option(False, help="Keep temporary files for debugging."),
    log_level: str = typer.Option('INFO', help="Logging level."),
    log_file: str | None = typer.Option(None, help="Optional path to write logs to."),
):
    """Search the test half of a split against the train half and flag the leakage.

    An MMseqs2 similarity search reports how much of the test set already appears
    in the training set, at or above --similarity-threshold. Reports land in
    out-folder/split/, one directory per compared pair.
    """
    _validate_input_files(train_input, 'Training input file')
    _validate_input_files(test_input, 'Test input file')

    # Both halves are searched against each other, so they must share a format
    format = _resolve_format(train_input + test_input)

    _validate_choices(report_types, SPLITS_REPORT_TYPES, 'report type')
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
