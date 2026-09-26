"""Per sequence dinucleotide content: the nucleotide check, for two-base combinations.

Each dinucleotide is scored separately, as the AU-ROC of its frequency in a
sequence, and the headline is the worst of them.
"""

from genomic_benchmarks_qc.checks import Check
from genomic_benchmarks_qc.checks.scorers import column_features, stat
from genomic_benchmarks_qc.report.sections import DISJOINT_BASES_MESSAGE, Section, SectionContent
from genomic_benchmarks_qc.report.utils import image_or_message, save_figure

NAME = 'Per sequence dinucleotide content'


def render(context) -> SectionContent:
    """One panel per pair of shared bases, with the flagged ones marked."""
    # Deferred, like every figure: see `report_generator`.
    from genomic_benchmarks_qc.report import classes_plots

    # With no base in common there is no panel to draw, and the section says so.
    plot = None
    if context.bases_overlap:
        bases = context.bases_overlap
        fig = classes_plots.plot_dinucleotides(context.stats1, context.stats2,
                                               nucleotides=bases, plot_type=context.plot_type)
        plot = save_figure(
            fig, context.plots_dir, 'per_sequence_dinucleotide_content',
            flagged=context.flagged(NAME),
            mark=lambda fig, flagged: classes_plots.mark_failed_dinucleotides(fig, bases, flagged))
    return SectionContent({'{{plot}}': image_or_message(
        plot, 'Per Sequence Dinucleotide Content', 'plot-wide', DISJOINT_BASES_MESSAGE)})


CHECK = Check(
    name=NAME,
    score=column_features(NAME, stat('Per sequence dinucleotide content')),
    floor='per_sequence',
    section=Section(
        title='Per Sequence Dinucleotide Content',
        anchor='per-sequence-dinucleotide-content',
        explanation_id='dinucleotide-explanation',
        template='check_dinucleotide_content.html',
        render=render,
        docs=(('checks', 'per-sequence-dinucleotide-content', 'What to do about it'),),
    ),
)
