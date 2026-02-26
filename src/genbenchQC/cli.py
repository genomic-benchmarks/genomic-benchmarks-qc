import typer
from typing import List, Optional
from pathlib import Path

from genbenchQC.evaluate_classes import run as run_evaluate_classes
from genbenchQC.evaluate_splits import run as run_evaluate_splits

app = typer.Typer(no_args_is_help=True)

# Valid choices for validation
VALID_FORMATS = ['fasta', 'csv', 'csv.gz', 'tsv', 'tsv.gz']
VALID_REPORT_TYPES = ['json', 'html', 'simple']
VALID_PLOT_TYPES = ['boxen', 'violin']
VALID_LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

@app.command()
def evaluate_classes(
    input: List[str] = typer.Option(..., help="Input file(s)."),
    format: str = typer.Option(..., help="Input format: fasta, csv, csv.gz, tsv, tsv.gz"),
    sequence_column: List[str] = typer.Option(['sequence'], help="One or more sequence column names for CSV/TSV inputs."),
    label_column: str = typer.Option('label', help="Label column name for single-file CSV inputs."),
    label_list: List[str] = typer.Option(['infer'], help='List of labels to consider or "infer" to detect labels automatically.'),
    regression: bool = typer.Option(False, help="Treat label column as regression target and split into high/low."),
    out_folder: str = typer.Option('.', help="Output folder for reports."),
    report_types: List[str] = typer.Option(['html', 'simple'], help="Types of reports to generate (json, html, simple)."),
    end_position: Optional[int] = typer.Option(None, help="End position for per-position stats."),
    plot_type: str = typer.Option('boxen', help="Plot type to use for visualizations (boxen, violin)."),
    log_level: str = typer.Option('INFO', help="Logging level."),
    log_file: Optional[str] = typer.Option(None, help="Optional path to write logs to."),
):
    """
    Evaluate sequence characteristics across different classes/labels in the dataset.
    """
    # Validate format
    if format not in VALID_FORMATS:
        typer.echo(f"Error: Invalid format '{format}'. Must be one of: {', '.join(VALID_FORMATS)}", err=True)
        raise typer.Exit(code=1)
    
    # Validate input files exist
    for file_path in input:
        if not Path(file_path).is_file():
            typer.echo(f"Error: Input file does not exist: {file_path}", err=True)
            raise typer.Exit(code=1)
    
    # Validate fasta format requires at least 2 files
    if format == 'fasta' and len(input) < 2:
        typer.echo("Error: When format is 'fasta', at least 2 input files are required (one per class).", err=True)
        raise typer.Exit(code=1)
    
    # Validate report_types
    for rt in report_types:
        if rt not in VALID_REPORT_TYPES:
            typer.echo(f"Error: Invalid report type '{rt}'. Must be one of: {', '.join(VALID_REPORT_TYPES)}", err=True)
            raise typer.Exit(code=1)
    
    # Validate plot_type
    if plot_type not in VALID_PLOT_TYPES:
        typer.echo(f"Error: Invalid plot type '{plot_type}'. Must be one of: {', '.join(VALID_PLOT_TYPES)}", err=True)
        raise typer.Exit(code=1)
    
    # Validate log_level
    if log_level not in VALID_LOG_LEVELS:
        typer.echo(f"Error: Invalid log level '{log_level}'. Must be one of: {', '.join(VALID_LOG_LEVELS)}", err=True)
        raise typer.Exit(code=1)
    
    # Validate end_position is positive
    if end_position is not None and end_position <= 0:
        typer.echo(f"Error: end_position must be a positive integer, got {end_position}", err=True)
        raise typer.Exit(code=1)
    
    run_evaluate_classes(
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
    format: str = typer.Option(..., help="Format of the input files: fasta, csv, csv.gz, tsv, tsv.gz"),
    sequence_column: List[str] = typer.Option(['sequence'], help="One or more sequence column names for CSV/TSV inputs."),
    out_folder: str = typer.Option('.', help="Output folder for reports."),
    report_types: List[str] = typer.Option(['html', 'simple'], help="Types of reports to generate (json, html, simple)."),
    similarity_threshold: float = typer.Option(80.0, help="Similarity threshold for data leakage detection (%)."),
    log_level: str = typer.Option('INFO', help="Logging level."),
    log_file: Optional[str] = typer.Option(None, help="Optional path to write logs to."),
):
    """
    Evaluate data leakage in dataset train-test split.
    """
    # Validate format
    if format not in VALID_FORMATS:
        typer.echo(f"Error: Invalid format '{format}'. Must be one of: {', '.join(VALID_FORMATS)}", err=True)
        raise typer.Exit(code=1)
    
    # Validate train input files exist
    for file_path in train_input:
        if not Path(file_path).is_file():
            typer.echo(f"Error: Training input file does not exist: {file_path}", err=True)
            raise typer.Exit(code=1)
    
    # Validate test input files exist
    for file_path in test_input:
        if not Path(file_path).is_file():
            typer.echo(f"Error: Test input file does not exist: {file_path}", err=True)
            raise typer.Exit(code=1)
    
    # Validate report_types
    for rt in report_types:
        if rt not in VALID_REPORT_TYPES:
            typer.echo(f"Error: Invalid report type '{rt}'. Must be one of: {', '.join(VALID_REPORT_TYPES)}", err=True)
            raise typer.Exit(code=1)
    
    # Validate log_level
    if log_level not in VALID_LOG_LEVELS:
        typer.echo(f"Error: Invalid log level '{log_level}'. Must be one of: {', '.join(VALID_LOG_LEVELS)}", err=True)
        raise typer.Exit(code=1)
    
    # Validate thresholds are in valid range [0, 1]
    if not 0 <= similarity_threshold <= 100:
        typer.echo(f"Error: similarity_threshold must be between 0 and 100, got {similarity_threshold}", err=True)
        raise typer.Exit(code=1)

    run_evaluate_splits(
        train_files=train_input,
        test_files=test_input,
        format=format,
        out_folder=out_folder,
        sequence_column=sequence_column,
        report_types=report_types,
        similarity_threshold=similarity_threshold,
        log_level=log_level,
        log_file=log_file,
    )

def main():
    app()

if __name__ == "__main__":
    main()