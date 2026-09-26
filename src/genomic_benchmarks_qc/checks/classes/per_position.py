"""Per position nucleotide content, read forward and from the sequence end.

Two checks built by one factory: at every position, can the base there tell the
classes apart? Each base at each position is scored separately, on the
sequences long enough to reach the position, and the headline is the worst
single one. The reversed check counts positions from the last base, which is
what catches something anchored to the 3' end of variable-length sequences.

Which positions are compared at all is decided by the cohort floors in
`SequenceStatistics`; `position_windows` combines the two classes' windows.
"""

from genomic_benchmarks_qc.checks import FLAGGED, Check, CheckResult
from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.per_position_payload import X_LABELS, build_payload, viewer_html
from genomic_benchmarks_qc.report.sections import DISJOINT_BASES_MESSAGE, Section, SectionContent
from genomic_benchmarks_qc.report.utils import escape_html_text, save_figure
from genomic_benchmarks_qc.utils.testing import (
    MIN_SEQUENCES_PER_CLASS,
    _aggregate_worst_case_metrics,
    _score_position_features,
    position_windows,
)

FORWARD_NAME = 'Per position nucleotide content'
REVERSED_NAME = 'Per position reversed nucleotide content'

# The interactive figure's behaviour, and the call that starts it once the page
# is laid out. Asked for by both directions and put on the page once.
VIEWER_SCRIPT = (assets.script('per_position_viewer.js')
                 + '\n<script>initPerPositionViewers();</script>\n')


def _flagged_positions(rows, name, bases, end_position) -> dict:
    """The positions a per-position check flagged, keyed by base then position.

    Only flagged positions go in, and a base only once it has one. An entry per
    base regardless would leave the dict truthy for a comparison with nothing to
    shade, and the report would draw a second, identical copy of every
    per-position plot.

    Args:
        rows: The check's rows, as `_score_position_features` names them.
        name: The check's name.
        bases: The bases scored, in the order they were scored.
        end_position: Last position scored, 1-based and inclusive.

    Returns:
        Dict of base -> {position: flag}, positions 1-based.
    """
    flagged = {}
    for base in bases:
        for position in range(1, end_position + 1):
            flag = rows[f'{name} - {base} position {position}']['Flag']
            if flag in FLAGGED:
                flagged.setdefault(base, {})[position] = flag
    return flagged


def position_features(name, reverse):
    """Score every base at every position, reading from the start or the end."""
    def score(stats1, stats2) -> CheckResult:
        bases = sorted(set(stats1.stats['Unique bases']) | set(stats2.stats['Unique bases']))
        end_position, scored_end_position = position_windows(stats1, stats2)

        rows, per_base = _score_position_features(
            stats1.sequences, stats2.sequences, bases, name,
            end_position=end_position, scored_end_position=scored_end_position,
            reverse=reverse,
        )
        if per_base:
            rows[name] = _aggregate_worst_case_metrics(per_base.values())
        return CheckResult(rows, _flagged_positions(rows, name, bases, end_position))
    return score


def generate_position_window_html(stats1, stats2):
    """Describe the per-position window and how well it is covered.

    The figure stops where the comparison did, so nothing in it has to be read
    with a caveat - and the caveat itself, that there are later positions the
    report saw and could not compare, still has to be said somewhere or the
    figure quietly stands in for the whole sequence. It is said here. The one
    comparison whose figure does carry a caveat is the one where nothing could be
    compared at all: there the panels fall back to the reported window, as the
    rest of the report's plots do when a comparison is underpowered, and the
    caveat covers all of it.

    It goes inside the section's ? explanation rather than above the figure: it is
    background for reading the plot, not a finding, and on its own above the
    figure it read as a second explanation of the same check.
    """
    end_position, scored_end_position = position_windows(stats1, stats2)

    if end_position < 1:
        return ('<p>Sequences are too short for per-position '
                'statistics, so no positions were analysed.</p>')

    def coverage_at(position):
        parts = []
        for stats in (stats1, stats2):
            name = stats.label if stats.label is not None else stats.filename
            label = escape_html_text(str(name))
            coverage = stats.coverage_at(position)
            count = int(round(coverage * stats.stats['Number of sequences']))
            parts.append(f"{label}: {coverage:.1%} ({count:,} sequences)")
        return f'{parts[0]} and {parts[1]}'

    def required_cohorts():
        """What each class needed behind a position, in its own terms.

        Collapsed to one phrase when both classes require the same number, which
        is the common case and reads as padding spelled out twice.
        """
        needed = []
        for stats in (stats1, stats2):
            count = stats.stats['Number of sequences']
            needed.append(stats._required_cohort(count) if count else 0)
        if needed[0] == needed[1]:
            return f'{needed[0]:,} sequences in each class'
        parts = []
        for stats, count in zip((stats1, stats2), needed, strict=True):
            name = stats.label if stats.label is not None else stats.filename
            label = escape_html_text(str(name))
            parts.append(f'{label}: {count:,}')
        return f'{parts[0]} and {parts[1]} sequences'

    if scored_end_position < 1:
        return (f'<p>No position could be compared. A position is compared only where '
                f'{required_cohorts()} reach it, and none does, so all {end_position} positions '
                f'the report covers are reported as Unknown rather than compared. The figure is '
                f'still drawn, over all {end_position} of them, because the frequencies are worth '
                f'looking at whether or not they were scored - but nothing in it carries a flag, '
                f'and a difference visible in it is not a difference this report is standing '
                f'behind. {coverage_at(end_position)} reach the last position drawn.</p>')

    # The common case, and the one where the cohort rule never bites: fixed-length
    # sequences, so the window is the whole sequence and every sequence is behind
    # every position of it. Spelling out the floors there costs eighty words to
    # explain a boundary the data never came near. What the rule is stays one
    # link away, in the same place the rest of the explanation sends the reader.
    whole = (end_position == scored_end_position
             and stats1.coverage_at(end_position) >= 1
             and stats2.coverage_at(end_position) >= 1)
    if whole:
        return (f'<p>Positions 1&ndash;{end_position} were compared, all of them: every '
                f'sequence in both labels reaches every position. Everything the figure '
                f'draws was scored, so a stretch with no band on it is a stretch that '
                f'passed.</p>')

    # The share each class was asked for, not the share that turned out to be
    # binding: on most datasets the latter is the sequence count restated as a
    # percentage, which explains nothing. The counts themselves come from the
    # boundary position below, where the binding class sits exactly on its floor.
    requested = min(stats1.min_coverage, stats2.min_coverage)
    parts = [f'<p>Positions 1&ndash;{scored_end_position} were compared, and those are the '
             f'positions the figure draws: everything in it was scored, so a stretch with no '
             f'flag on it is a stretch that passed. '
             f'Each position is compared on the sequences that reach it, and only where enough of '
             f'them do: the larger of {MIN_SEQUENCES_PER_CLASS} sequences &mdash; below which a '
             'difference this size turns up on sampling noise alone &mdash; and '
             f'{requested:.0%} of the class, below which a cohort can be large and still describe '
             f'only the longest sequences. Position {scored_end_position} is the last that clears '
             f'both: {coverage_at(scored_end_position)}.</p>']

    if end_position > scored_end_position:
        parts.append(
            f'<p>The sequences run further, to position {end_position}, where '
            f'{coverage_at(end_position)} remain. Those later positions are not drawn and are '
            'reported as Unknown rather than compared: they are reached by too few of each class '
            'for a difference there to be a difference between the classes rather than between '
            'their longest sequences. So the figure ends before the sequences do, and says '
            'nothing either way about what happens past its right-hand edge. The panel at the '
            'bottom shows how the number of sequences behind each class falls along the window: '
            'the window ends where the lower of the two curves falls below the cohort a position '
            'has to have behind it.</p>')

    return ''.join(parts)


def position_section(name, direction, dom_id, stem):
    """Render one direction's section: the interactive figure, and its PNG in plots/.

    The page draws the panels in the browser from the numbers behind them rather
    than embedding a picture: a flagged position is one pixel wide in a
    400-position window, which is what the zoom exists for. The same figure is
    still written to plots/ as a PNG, it is just not what the page shows.
    """
    def render(context) -> SectionContent:
        # Deferred, like every figure: see `report_generator`.
        from genomic_benchmarks_qc.report import classes_plots

        bases, end_position = context.bases_overlap, context.drawn_end
        # Drawn over the window `drawn_window` resolves: what was compared, or -
        # when nothing was - the positions the checks are named for, every one of
        # them Unknown. Only a comparison with no per-position checks at all, or
        # no base in common, leaves the figure unmade. The axis wording comes
        # from X_LABELS, which the interactive figure and its flag table read
        # too, so the PNG and the page cannot count positions differently.
        if bases and end_position >= 1:
            fig = classes_plots.plot_per_base_sequence_comparison(
                context.stats1, context.stats2, stats_name=name, nucleotides=bases,
                end_position=end_position, x_label=X_LABELS[direction])
            save_figure(
                fig, context.plots_dir, stem, flagged=context.flagged(name),
                mark=lambda fig, flagged: classes_plots.mark_failed_positions(
                    fig, bases, end_position, flagged),
                embed=False)

        payload = build_payload(context.stats1, context.stats2, bases, end_position,
                                context.results, direction, context.compared)
        fills = {'{{position_window_note}}': generate_position_window_html(context.stats1,
                                                                           context.stats2)}
        if payload is not None:
            fills['{{plot}}'] = viewer_html(payload, dom_id)
            return SectionContent(fills, (VIEWER_SCRIPT,))

        # No payload for either of two reasons, and only one of them is about the
        # bases: the labels share none, or the comparison has no position to
        # draw at all. Saying the first where the second holds would explain the
        # empty section with something that is not true of it.
        message = ('There is nothing to plot for this comparison.' if bases
                   else DISJOINT_BASES_MESSAGE)
        fills['{{plot}}'] = f'<p class="no-plot-message">{message}</p>'
        return SectionContent(fills)
    return render


FORWARD = Check(
    name=FORWARD_NAME,
    score=position_features(FORWARD_NAME, reverse=False),
    floor='per_position',
    section=Section(
        title='Per Position Nucleotide Content',
        anchor='per-position-nucleotide-content',
        explanation_id='per-position-explanation',
        template='check_per_position.html',
        render=position_section(FORWARD_NAME, 'forward', 'ppv-fwd',
                                'per_position_nucleotide_content'),
        # The one explanation with two places to go on to: what the check means,
        # and how to drive the figure it is attached to.
        docs=(('checks', 'per-position-nucleotide-content', 'What to do about it'),
              ('viewer', None, 'Reading the figure')),
        stylesheets=('per_position_viewer.css',),
    ),
)

REVERSED = Check(
    name=REVERSED_NAME,
    score=position_features(REVERSED_NAME, reverse=True),
    floor='per_position',
    section=Section(
        title='Per Position Reversed Nucleotide Content',
        anchor='per-position-reversed-nucleotide-content',
        explanation_id='per-position-rev-explanation',
        template='check_per_position_reversed.html',
        render=position_section(REVERSED_NAME, 'reversed', 'ppv-rev',
                                'per_position_reversed_nucleotide_content'),
        docs=(('checks', 'per-position-reversed-nucleotide-content', 'More on this check'),),
        stylesheets=('per_position_viewer.css',),
    ),
)
