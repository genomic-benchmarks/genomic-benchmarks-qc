"""The HTML template for the class comparison report.

Builds a single self-contained page: the plots are embedded as base64 data URIs
and the styling is inlined, so a report can be moved, zipped or served from
anywhere without losing anything.
"""

import html
from datetime import datetime

from genomic_benchmarks_qc import __version__
from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.per_position_payload import viewer_html
from genomic_benchmarks_qc.report.utils import (
    COMMON_CSS,
    LOGO_BASE64,
    REPORT_HEADER_HTML,
    SIDEBAR_LINKS_HTML,
    TOOL_DESCRIPTION,
    TOOL_TAGLINE,
    docs_link,
    encode_image_to_base64,
    escape_str,
    icon_html,
    image_or_message,
    put_data,
    put_text,
    verdict_html,
)
from genomic_benchmarks_qc.utils.testing import MIN_SEQUENCES_PER_CLASS, position_windows

# The nine headline checks, in the order the navigation lists them. The flags
# behind them are what the verdict line counts; every other key in
# summary_statuses is a sub-check one of these already summarises.
CHECK_NAMES = (
    'Unique bases',
    'Sequence Duplications within Labels',
    'Duplicate Sequences between Labels',
    'Sequence lengths',
    'Per sequence GC content',
    'Per sequence nucleotide content',
    'Per sequence dinucleotide content',
    'Per position nucleotide content',
    'Per position reversed nucleotide content',
)

# The link that closes each ? explanation, keyed by the placeholder that holds
# it. Each anchor is a heading in docs/guide/checks.md; tests/test_report_links.py
# checks that they all still exist.
EXPLANATION_LINKS = {
    '{{link_unique_bases}}': ('checks', 'unique-bases', 'What to do about it'),
    '{{link_within_dup}}': ('checks', 'sequence-duplications-within-labels',
                            'What to do about it'),
    '{{link_between_dup}}': ('checks', 'duplicate-sequences-between-labels',
                             'What to do about it'),
    '{{link_lengths}}': ('checks', 'sequence-lengths', 'What to do about it'),
    '{{link_gc}}': ('checks', 'per-sequence-gc-content', 'What to do about it'),
    '{{link_nucleotide}}': ('checks', 'per-sequence-nucleotide-content',
                            'What to do about it'),
    '{{link_dinucleotide}}': ('checks', 'per-sequence-dinucleotide-content',
                              'What to do about it'),
    '{{link_per_position_rev}}': ('checks', 'per-position-reversed-nucleotide-content',
                                  'More on this check'),
}


def generate_nucleotide_flags_html(summary_statuses, flag_prefix):
    """
    Generate HTML for nucleotide-level flags.

    Args:
        summary_statuses: dict with summary statuses for various checks.
        flag_prefix: Prefix for the flag keys (e.g., "Per position nucleotide content")

    Returns:
        HTML string with nucleotide flags
    """
    if summary_statuses is None:
        return ''

    flags_html = ''

    for flag in summary_statuses:
        if not flag.startswith(f"{flag_prefix} - "):
            continue

        # ignore per-position entries like "... - A position 1" and only show
        # per-nucleotide aggregated flags like "... - A".
        remainder = flag.split(" - ", 1)[1]
        if ' position ' in remainder:
            continue

        nt = remainder
        flag_value = summary_statuses.get(flag, '')

        # Determine status class based on flag value
        status_class = ''
        symbol = '?'
        if isinstance(flag_value, str):
            flag_lower = flag_value.lower()
            if flag_lower in ('pass', 'ok', 'good', 'success'):
                status_class = 'status-pass'
                symbol = '✔'
            elif flag_lower in ('warn', 'warning'):
                status_class = 'status-warn'
                symbol = '!'
            elif flag_lower in ('fail', 'failed', 'error'):
                status_class = 'status-fail'
                symbol = '✖'
            else:
                status_class = 'status-unknown'
                symbol = '?'

        flags_html += f'''<div class="nucleotide-flag-item">
                <span class="status-icon-small {status_class}">{symbol}</span>
                <span class="nucleotide-label">{html.escape(nt)}</span>
            </div>
            '''

    return flags_html

HTML_TEMPLATE = assets.template('classes_report.html')

def generate_not_scored_html(stats1, stats2, summary_statuses):
    """Explain, at the top of the report, any check that was not scored.

    A grey "?" in the sidebar is easy to read as a pass, and the difference is
    the whole point: Unknown means the comparison was not made, because there
    were not enough sequences behind it for its result to mean anything. Saying
    so once, in full, is what lets a reader trust the flags that are there.

    Returns an empty string when every check was scored, so the note costs
    nothing on a normally sized dataset.
    """
    if not summary_statuses:
        return ''

    # Only the top-level checks; the per-base and per-position detail rows carry
    # their own Unknowns and would bury the message.
    unscored = [name for name, flag in summary_statuses.items()
                if ' - ' not in name and str(flag).strip().lower() == 'unknown']
    if not unscored:
        return ''

    positional = [name for name in unscored if 'position' in name.lower()]
    per_sequence = [name for name in unscored if name not in positional]

    def names(items):
        return ', '.join(f'<em>{name}</em>' for name in items)

    reasons = []
    if per_sequence:
        smaller = min(stats1.stats['Number of sequences'], stats2.stats['Number of sequences'])
        reasons.append(
            f'{names(per_sequence)} &mdash; the smaller class holds {smaller:,} sequences, and below '
            f'{MIN_SEQUENCES_PER_CLASS} these checks report a difference on sampling noise alone too '
            'often for a flag to be informative.'
        )
    if positional:
        reasons.append(
            f'{names(positional)} &mdash; every position is compared on the sequences long enough to '
            f'reach it, and no scored position has {MIN_SEQUENCES_PER_CLASS} sequences in both '
            'classes.'
        )

    items = ''.join(f'<li>{reason}</li>' for reason in reasons)
    return (
        '<div class="not-scored-note">'
        f'<strong>{len(unscored)} check(s) were not scored.</strong> They are marked '
        '<span class="status-icon status-unknown">?</span> rather than passed, because not enough '
        'data stands behind them to tell a real difference from sampling noise &mdash; which is not '
        'the same as having found no difference.'
        f'<ul>{items}</ul>'
        'The plots and the descriptive statistics below are computed from all the data and are '
        'unaffected, so the distributions can still be compared by eye.'
        '</div>'
    )


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
            label = html.escape(str(stats.label if stats.label is not None else stats.filename))
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
            label = html.escape(str(stats.label if stats.label is not None else stats.filename))
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


def get_dataset_html_template(stats1, stats2, plots_path, summary_statuses, duplicate_seqs, duplicate_seqs_file=None,
                             tool_description=None, per_position_payloads=None):
    """
    Returns the HTML template for the report.

    Args:
        stats1, stats2: objects containing dataset statistics (unchanged API).
        plots_path: dict with plot image paths.
        summary_statuses: dict with summary statuses for various checks.
        duplicate_seqs: list of duplicate sequences.
        duplicate_seqs_file: path to file with duplicate sequences (optional).
        tool_description: short description of the tool to include in the header (optional).
        per_position_payloads: dict keyed 'forward'/'reversed' with the data for
            the interactive per-position figures, as built by
            per_position_payload.build_payload. A direction that is absent or
            None renders the no-plot message instead.

    `summary_statuses` must contain the 'Sequence Duplications within Labels' flag,
    which selects between the duplication plot and the "no duplicates" message.
    """
    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment
    html_template = put_data(html_template, "{{common_css}}",
                             COMMON_CSS + assets.stylesheet('per_position_viewer.css'))
    html_template = put_data(html_template, "{{report_header}}", REPORT_HEADER_HTML)

    # insert logo
    html_template = put_data(html_template, "{{logo_base64}}", LOGO_BASE64)

    # populate header placeholders: tool description, generated timestamp and input paths
    # Provide sensible defaults when values are not supplied
    if tool_description is None:
        tool_description = TOOL_DESCRIPTION
    if stats1.filepath == stats2.filepath:
        input_paths = stats1.filepath
    else:
        input_paths = f"{stats1.filepath}, {stats2.filepath}"
    generated_on = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_template = put_data(html_template, "{{tool_tagline}}", TOOL_TAGLINE)
    html_template = put_data(html_template, "{{tool_description}}", tool_description)
    html_template = put_data(html_template, "{{generated_on}}", generated_on)
    html_template = put_data(html_template, "{{version}}", __version__)

    # What this particular report looked at, and how it came out. The subject
    # line is the one thing a reader handed the file cannot work out for
    # themselves, and the verdict saves them counting circles in the navigation.
    label1 = str(stats1.label) if stats1.label is not None else 'N/A'
    label2 = str(stats2.label) if stats2.label is not None else 'N/A'
    subject = (f'Comparing label <strong>{html.escape(label1)}</strong> with label '
               f'<strong>{html.escape(label2)}</strong> in '
               f'{html.escape(str(input_paths))}')
    html_template = put_data(html_template, "{{report_subject}}", subject)
    html_template = put_data(html_template, "{{report_verdict}}",
                             verdict_html(summary_statuses, CHECK_NAMES))
    html_template = put_data(html_template, "{{sidebar_links}}", SIDEBAR_LINKS_HTML)
    # Long enough to identify the report in a tab or a bookmark, ordered so that
    # a truncated tab keeps the part that distinguishes it from its neighbours.
    # The sequence column is last and still has to be there: a dataset with
    # several of them produces one report per column, and paired-sequences
    # produces three whose file and labels are identical.
    seq_col = str(stats1.seq_column) if stats1.seq_column is not None else None
    column = f' ({seq_col})' if seq_col else ''
    html_template = put_data(
        html_template, "{{page_title}}",
        html.escape(f'{stats1.filename}: {label1} vs {label2}{column} - gb-qc'))

    for placeholder, (page, anchor, text) in EXPLANATION_LINKS.items():
        html_template = put_data(html_template, placeholder,
                                 docs_link(page, anchor, text))
    # The one explanation with two places to go on to: what the check means, and
    # how to drive the figure it is attached to.
    html_template = put_data(
        html_template, "{{link_per_position}}",
        docs_link('checks', 'per-position-nucleotide-content', 'What to do about it')
        + docs_link('viewer', text='Reading the figure'))

    html_template = put_text(html_template, "{{filename1}}", stats1.filename)
    html_template = put_text(html_template, "{{filename2}}", stats2.filename)
    html_template = put_text(html_template, "{{label1}}", label1)
    html_template = put_text(html_template, "{{label2}}", label2)
    html_template = put_text(html_template, "{{seq_col1}}", str(stats1.seq_column) if stats1.seq_column is not None else "N/A")
    html_template = put_text(html_template, "{{seq_col2}}", str(stats2.seq_column) if stats2.seq_column is not None else "N/A")
    html_template = put_data(html_template, "{{number_of_sequences1}}", str(stats1.stats['Number of sequences']))
    html_template = put_data(html_template, "{{number_of_sequences2}}", str(stats2.stats['Number of sequences']))
    html_template = put_data(html_template, "{{dedup_sequences1}}", str(stats1.stats['Number of sequences left after deduplication']))
    html_template = put_data(html_template, "{{dedup_sequences2}}", str(stats2.stats['Number of sequences left after deduplication']))
    html_template = put_data(html_template, "{{min_length1}}", str(int(stats1.stats['Sequence lengths']['Sequence lengths'].min())))
    html_template = put_data(html_template, "{{min_length2}}", str(int(stats2.stats['Sequence lengths']['Sequence lengths'].min())))
    html_template = put_data(html_template, "{{mean_length1}}", f"{stats1.stats['Sequence lengths']['Sequence lengths'].mean():.2f}")
    html_template = put_data(html_template, "{{mean_length2}}", f"{stats2.stats['Sequence lengths']['Sequence lengths'].mean():.2f}")
    html_template = put_data(html_template, "{{max_length1}}", str(int(stats1.stats['Sequence lengths']['Sequence lengths'].max())))
    html_template = put_data(html_template, "{{max_length2}}", str(int(stats2.stats['Sequence lengths']['Sequence lengths'].max())))
    html_template = put_data(html_template, "{{number_of_bases1}}", str(stats1.stats['Number of bases']))
    html_template = put_data(html_template, "{{number_of_bases2}}", str(stats2.stats['Number of bases']))
    html_template = put_text(html_template, "{{unique_bases1}}", ', '.join(x for x in stats1.stats['Unique bases']))
    html_template = put_text(html_template, "{{unique_bases2}}", ', '.join(x for x in stats2.stats['Unique bases']))
    html_template = put_data(html_template, "{{gc_content1}}", f"{(stats1.stats['%GC content']*100):.2f}")
    html_template = put_data(html_template, "{{gc_content2}}", f"{(stats2.stats['%GC content']*100):.2f}")

    if summary_statuses['Sequence Duplications within Labels'].lower() in ('pass', 'ok', 'good', 'success'):
        # no duplicate sequences found
        duplication_message = """
        <p>No duplicate sequences were found in either class.</p>
        """
        html_template = put_data(html_template, "{{sequence_duplications_within_classes}}", duplication_message)
    else:
        # insert plot showing duplicate sequences
        html_template = put_data(html_template, "{{sequence_duplications_within_classes}}",
                                 f'<img src="data:image/png;base64, {encode_image_to_base64(plots_path["Sequence Duplications within Labels"])}" alt="Sequence Duplications within Labels Plot" class="plot-wide">')
    html_template = put_data(html_template, "{{sequence_length_plot_base64}}",
                             encode_image_to_base64(plots_path['Sequence lengths']))
    html_template = put_data(html_template, "{{per-sequence-gc-content_base64}}",
                             encode_image_to_base64(plots_path['Per sequence GC content']))
    # Sections that require a shared set of bases between the two labels. When
    # the datasets have no bases in common the plot path is None, so we render
    # an explanatory message instead of an image.
    disjoint_bases_message = ('No bases in common between the two labels, so this comparison '
                              'cannot be plotted.')
    html_template = put_data(html_template, "{{per-sequence-nucleotide-content}}",
                             image_or_message(plots_path['Per sequence nucleotide content'],
                                              'Per Sequence Nucleotide Content', 'plot-wide',
                                              disjoint_bases_message))
    html_template = put_data(html_template, "{{per-sequence-dinucleotide-content}}",
                             image_or_message(plots_path['Per sequence dinucleotide content'],
                                              'Per Sequence Dinucleotide Content', 'plot-wide',
                                              disjoint_bases_message))
    html_template = put_data(html_template, "{{not_scored_note}}",
                             generate_not_scored_html(stats1, stats2, summary_statuses))
    position_window_note = generate_position_window_html(stats1, stats2)
    html_template = put_data(html_template, "{{position_window_note}}", position_window_note)
    html_template = put_data(html_template, "{{position_window_note_reversed}}", position_window_note)
    # The per-position panels are drawn in the browser from the numbers behind
    # them, not embedded as a picture: a flagged position is one pixel wide in a
    # 400-position window, which is what the zoom exists for. The same figure is
    # still written to the plots directory as a PNG, it is just not what the page
    # shows.
    per_position_payloads = per_position_payloads or {}
    # A direction arrives without a payload for either of two reasons, and only
    # one of them is about the bases: the labels share none, or the comparison
    # has no position to draw at all. Saying the first where the second holds
    # would explain the empty section with something that is not true of it.
    if set(stats1.stats['Unique bases']) & set(stats2.stats['Unique bases']):
        no_payload_message = 'There is nothing to plot for this comparison.'
    else:
        no_payload_message = disjoint_bases_message
    for placeholder, direction, dom_id, _name in (
        ("{{per-position-nucleotide-content}}", 'forward', 'ppv-fwd',
         'Per Position Nucleotide Content'),
        ("{{per-position-reversed-nucleotide-content}}", 'reversed', 'ppv-rev',
         'Per Position Reversed Nucleotide Content'),
    ):
        payload = per_position_payloads.get(direction)
        if payload is None:
            html_template = put_data(html_template, placeholder,
                                     f'<p class="no-plot-message">{no_payload_message}</p>')
        else:
            html_template = put_data(html_template, placeholder, viewer_html(payload, dom_id))

    # Populate sidebar icon placeholders (if provided). summary_statuses may contain
    # simple status keywords ('pass', 'warn', 'fail') or an HTML snippet. This helper
    # returns the HTML for the small circular icon shown before each section link.
    html_template = put_data(html_template, "{{icon_unique_bases}}", icon_html(summary_statuses, 'Unique bases'))
    html_template = put_data(html_template, "{{icon_sequence_duplications_within_classes}}", icon_html(summary_statuses, 'Sequence Duplications within Labels'))
    html_template = put_data(html_template, "{{icon_sequence_lengths}}", icon_html(summary_statuses, 'Sequence lengths'))
    html_template = put_data(html_template, "{{icon_sequence_duplication_levels}}", icon_html(summary_statuses, 'Duplicate Sequences between Labels'))
    html_template = put_data(html_template, "{{icon_per_sequence_nucleotide_content}}", icon_html(summary_statuses, 'Per sequence nucleotide content'))
    html_template = put_data(html_template, "{{icon_per_sequence_dinucleotide_content}}", icon_html(summary_statuses, 'Per sequence dinucleotide content'))
    html_template = put_data(html_template, "{{icon_per_position_nucleotide_content}}", icon_html(summary_statuses, 'Per position nucleotide content'))
    html_template = put_data(html_template, "{{icon_per_position_reversed_nucleotide_content}}", icon_html(summary_statuses, 'Per position reversed nucleotide content'))
    html_template = put_data(html_template, "{{icon_per_sequence_gc_content}}", icon_html(summary_statuses, 'Per sequence GC content'))

    if duplicate_seqs == []:
        # no duplicate sequences found
        duplication_message = """
        <p>No duplicate sequences were found between classes.</p>
        """
        html_template = put_data(html_template, "{{sequence_duplication_levels}}", duplication_message)
    else:
        # insert table showing duplicate sequences
        duplication_table = """
                <table>
                    <thead>
                        <tr>
                            <th class="sequence_column">Sequence</th>
                        </tr>
                    </thead>
                        <tbody>
                            <!-- Table rows will be dynamically populated -->
                        </tbody>
                </table>
                <div id="sequence-duplication-levels-info">
                    <p>{{sequence_duplication_levels_rest}} {{sequence_duplication_levels_file}}</p>
                </div>
        """
        # The sequences go in as JSON data rather than as a JavaScript literal, so
        # report_ui.js stays a file that can be linted and a sequence cannot
        # break out of the script element. Escaped for the same reason as the
        # per-position payload.
        escaped_seqs = [escape_str(seq) for seq in duplicate_seqs[:10]]
        duplication_table += ('\n<script type="application/json" id="duplicate-sequences">['
                              + ', '.join(escaped_seqs) + ']</script>')
        html_template = put_data(html_template, "{{sequence_duplication_levels}}", duplication_table)
        if len(duplicate_seqs) > 10:
            # If there are more than 10 sequences, we show how many more there are
            # and set the rest to a placeholder
            html_template = put_data(html_template, "{{sequence_duplication_levels_rest}}", f"And {str(len(duplicate_seqs) - 10)} more.")
        else:
            html_template = put_data(html_template, "{{sequence_duplication_levels_rest}}", "")

        if duplicate_seqs_file is not None:
            html_template = put_data(
                html_template, "{{sequence_duplication_levels_file}}",
                f"All {len(duplicate_seqs)} duplicate sequences saved to "
                f"{html.escape(str(duplicate_seqs_file))}.")
        else:
            html_template = put_data(html_template, "{{sequence_duplication_levels_file}}", "")

    # The behaviour goes in at the end of the body so it runs against a complete
    # page: the viewer measures the width of the box it is drawn into, which is
    # only known once the layout exists. The viewer and its data are left out
    # entirely when there is nothing to draw.
    scripts = ''
    if any(payload is not None for payload in per_position_payloads.values()):
        scripts += assets.script('per_position_viewer.js')
        scripts += '\n<script>initPerPositionViewers();</script>\n'
    scripts += assets.script('report_ui.js')
    html_template = put_data(html_template, "{{report_scripts}}", scripts)

    # Generate unique bases flags
    unique_bases_flags = generate_nucleotide_flags_html(summary_statuses, 'Unique bases')
    return put_data(html_template, "{{unique_bases_flags}}", unique_bases_flags)

