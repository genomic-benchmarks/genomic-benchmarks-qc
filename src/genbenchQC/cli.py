import typer
from typing import List, Optional

from genbenchQC.evaluate_classes import run as run_evaluate_classes
from genbenchQC.evaluate_splits import run as run_evaluate_splits

app = typer.Typer(no_args_is_help=True)

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
    seq_report_types: List[str] = typer.Option([], help="Sequence-level report types (json, html)."),
    end_position: Optional[int] = typer.Option(None, help="End position for per-position stats."),
    plot_type: str = typer.Option('boxen', help="Plot type to use for visualizations (boxen, violin)."),
    log_level: str = typer.Option('INFO', help="Logging level."),
    log_file: Optional[str] = typer.Option(None, help="Optional path to write logs to."),
):
    """
    Evaluate sequence characteristics across different classes/labels in the dataset.
    """
    run_evaluate_classes(
        input=input,
        format=format,
        out_folder=out_folder,
        sequence_column=sequence_column,
        label_column=label_column,
        label_list=label_list,
        regression=regression,
        report_types=report_types,
        seq_report_types=seq_report_types,
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
    identity_threshold: float = typer.Option(0.8, help="Identity threshold for clustering."),
    alignment_coverage: float = typer.Option(0.8, help="Alignment coverage for clustering."),
    log_level: str = typer.Option('INFO', help="Logging level."),
    log_file: Optional[str] = typer.Option(None, help="Optional path to write logs to."),
):
    """
    Evaluate data leakage in dataset train-test split.
    """
    run_evaluate_splits(
        train_files=train_input,
        test_files=test_input,
        format=format,
        out_folder=out_folder,
        sequence_column=sequence_column,
        report_types=report_types,
        identity_threshold=identity_threshold,
        alignment_coverage=alignment_coverage,
        log_level=log_level,
        log_file=log_file,
    )

def main():
    app()

if __name__ == "__main__":
    main()