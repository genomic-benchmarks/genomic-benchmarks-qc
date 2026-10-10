"""Sequence lengths: can the length of a sequence alone tell the classes apart?

The flag is the AU-ROC of length on its own - how well a model that does nothing
but count characters separates the two classes.
"""

from genomic_benchmarks_qc.checks import Check
from genomic_benchmarks_qc.checks.scorers import scalar_feature, stat
from genomic_benchmarks_qc.report.sections import Section, SectionContent
from genomic_benchmarks_qc.report.utils import encode_image_to_base64, save_figure

NAME = 'Sequence lengths'


def render(context) -> SectionContent:
    """The two length distributions."""
    # Deferred, like every figure: see `report_generator`.
    from genomic_benchmarks_qc.report import classes_plots

    fig = classes_plots.plot_lengths(context.stats1, context.stats2, plot_type=context.plot_type)
    plot = save_figure(fig, context.plots_dir, 'sequence_lengths')
    return SectionContent({'{{plot_base64}}': encode_image_to_base64(plot)})


CHECK = Check(
    name=NAME,
    score=scalar_feature(NAME, stat('Sequence lengths')),
    floor='per_sequence',
    section=Section(
        title='Sequence Lengths',
        anchor='sequence-lengths',
        explanation_id='lengths-explanation',
        template='check_lengths.html',
        render=render,
        docs=(('checks', 'sequence-lengths', 'What to do about it'),),
    ),
)
