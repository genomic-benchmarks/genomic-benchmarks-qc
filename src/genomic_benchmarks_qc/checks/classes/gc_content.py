"""Per sequence GC content: can a sequence's GC% alone tell the classes apart?

The flag is the AU-ROC of GC content on its own.
"""

from genomic_benchmarks_qc.checks import Check
from genomic_benchmarks_qc.checks.scorers import scalar_feature, stat
from genomic_benchmarks_qc.report.sections import Section, SectionContent
from genomic_benchmarks_qc.report.utils import encode_image_to_base64, save_figure

NAME = 'Per sequence GC content'


def render(context) -> SectionContent:
    """The two GC-content distributions."""
    # Deferred, like every figure: see `report_generator`.
    from genomic_benchmarks_qc.report import classes_plots

    fig = classes_plots.plot_gc_content(context.stats1, context.stats2,
                                        plot_type=context.plot_type)
    plot = save_figure(fig, context.plots_dir, 'per_sequence_gc_content')
    return SectionContent({'{{plot_base64}}': encode_image_to_base64(plot)})


CHECK = Check(
    name=NAME,
    score=scalar_feature(NAME, stat('Per sequence GC content')),
    floor='per_sequence',
    section=Section(
        title='Per Sequence GC Content',
        anchor='per-sequence-gc-content',
        explanation_id='gc-explanation',
        template='check_gc_content.html',
        render=render,
        docs=(('checks', 'per-sequence-gc-content', 'What to do about it'),),
    ),
)
