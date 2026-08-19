"""The HTML template for the class comparison report.

Builds a single self-contained page: the plots are embedded as base64 data URIs
and the styling is inlined, so a report can be moved, zipped or served from
anywhere without losing anything.
"""

from datetime import datetime
from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.per_position_payload import viewer_html
from genomic_benchmarks_qc.report.utils import put_data, encode_image_to_base64, image_or_message, escape_str, icon_html, COMMON_CSS, REPORT_HEADER_HTML, LOGO_BASE64
from genomic_benchmarks_qc.utils.seq_stats import DEFAULT_MIN_COVERAGE
from genomic_benchmarks_qc.utils.testing import MIN_SEQUENCES_PER_CLASS, MIN_SEQUENCES_PER_POSITION
import importlib.metadata


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
    
    for flag in summary_statuses.keys():
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
                <span class="nucleotide-label">{nt}</span>
            </div>
            '''
    
    return flags_html

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HTML Report Output</title>
    <style>{{common_css}}</style>
</head>
<body>
    <div class="container">
        <div class="sidebar">
            <div class="logo" style="text-align: center; margin-bottom: 20px;">
                <img src="{{logo_base64}}" alt="Genomic Benchmarks QC Logo" style="max-width: 150px; height: auto;">
                <span style="display: block; font-size: 14px; color: #555;">v{{version}}</span>
            </div>
            <h2>Summary</h2>
            <div class="sidebar-item">
                <span style="display: inline-block; width: 40px;"></span>
                <a href="#basic-descriptive-statistics">Basic Descriptive Statistics</a>
            </div>
            <div class="sidebar-item">
                {{icon_unique_bases}}
                <a href="#unique-bases">Unique Bases</a>
            </div>
            <div class="sidebar-item">{{icon_sequence_duplications_within_classes}}<a href="#sequence-duplications-within-classes">Sequence Duplications within Labels</a></div>
            <div class="sidebar-item">{{icon_sequence_duplication_levels}}<a href="#sequence-duplication-levels">Duplicate Sequences between Labels</a></div>
            <div class="sidebar-item">{{icon_sequence_lengths}}<a href="#sequence-lengths">Sequence lengths</a></div>
            <div class="sidebar-item">{{icon_per_sequence_gc_content}}<a href="#per-sequence-gc-content">Per Sequence GC Content</a></div>
            <div class="sidebar-item">{{icon_per_sequence_nucleotide_content}}<a href="#per-sequence-nucleotide-content">Per Sequence Nucleotide Content</a></div>
            <div class="sidebar-item">{{icon_per_sequence_dinucleotide_content}}<a href="#per-sequence-dinucleotide-content">Per Sequence Dinucleotide Content</a></div>
            <div class="sidebar-item">{{icon_per_position_nucleotide_content}}<a href="#per-position-nucleotide-content">Per Position Nucleotide Content</a></div>
            <div class="sidebar-item">{{icon_per_position_reversed_nucleotide_content}}<a href="#per-position-reversed-nucleotide-content">Per Position Reversed Nucleotide Content</a></div>
        </div>

        <div class="content">

            <!-- Report header: logo, short description and generated-on/data source info -->
            {{report_header}}

            {{not_scored_note}}

            <section id="basic-descriptive-statistics" class="table-section">
                <div class="section-header">
                    <h2>Basic Descriptive Statistics</h2>
                </div>
                <table style="width: 100%; border-collapse: collapse; margin: 20px 0;">
                    <tr>
                        <td><span>Filename</span></td>
                        <td style="text-align: center;">{{filename1}}</td>
                        <td style="text-align: center;">{{filename2}}</td>
                    </tr>
                    <tr>
                        <td><span>Label</span></td>
                        <td style="text-align: center;">{{label1}}</td>
                        <td style="text-align: center;">{{label2}}</td>
                    </tr>
                    <tr>
                        <td><span>Sequence column</span></td>
                        <td style="text-align: center;">{{seq_col1}}</td>
                        <td style="text-align: center;">{{seq_col2}}</td>
                    </tr>
                    <tr>
                        <td><span>Number of sequences</span></td>
                        <td style="text-align: center;">{{number_of_sequences1}}</td>
                        <td style="text-align: center;">{{number_of_sequences2}}</td>
                    </tr>
                    <tr>
                        <td><span>Unique sequences</span></td>
                        <td style="text-align: center;">{{dedup_sequences1}}</td>
                        <td style="text-align: center;">{{dedup_sequences2}}</td>
                    </tr>
                    <tr>
                        <td><span>Minimum length</span></td>
                        <td style="text-align: center;">{{min_length1}}</td>
                        <td style="text-align: center;">{{min_length2}}</td>
                    </tr>
                    <tr>
                        <td><span>Mean length</span></td>
                        <td style="text-align: center;">{{mean_length1}}</td>
                        <td style="text-align: center;">{{mean_length2}}</td>
                    </tr>
                    <tr>
                        <td><span>Maximum length</span></td>
                        <td style="text-align: center;">{{max_length1}}</td>
                        <td style="text-align: center;">{{max_length2}}</td>
                    </tr>
                    <tr>
                        <td><span>Number of bases</span></td>
                        <td style="text-align: center;">{{number_of_bases1}}</td>
                        <td style="text-align: center;">{{number_of_bases2}}</td>
                    </tr>
                    <tr>
                        <td><span>%GC content</span></td>
                        <td style="text-align: center;">{{gc_content1}}</td>
                        <td style="text-align: center;">{{gc_content2}}</td>
                    </tr>
                </table>
            </section>

            <section id="unique-bases" class="table-section">
                <div class="sidebar-item">
                    {{icon_unique_bases}}
                    <div class="section-header">
                        <h2>Unique Bases</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('unique-bases-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="unique-bases-explanation" class="explanation-text">
                    <strong>Unique Bases</strong> shows which nucleotides are present in each label. Differences in unique bases between labels may indicate data quality issues or biological differences that could bias machine learning models.
                </div>
                <div class="nucleotide-flags-container" id="unique-bases-flags">
                    {{unique_bases_flags}}
                </div>
                <table style="width: 100%; border-collapse: collapse; margin: 20px 0;">
                    <tr>
                        <td><span>Label</span></td>
                        <td style="text-align: center;">{{label1}}</td>
                        <td style="text-align: center;">{{label2}}</td>
                    </tr>
                    <tr>
                        <td><span>Unique bases</span></td>
                        <td style="text-align: center;">{{unique_bases1}}</td>
                        <td style="text-align: center;">{{unique_bases2}}</td>
                    </tr>
                </table>
            </section>

            <section id="sequence-duplications-within-classes">
                <div class="sidebar-item">
                    {{icon_sequence_duplications_within_classes}}
                    <div class="section-header">
                        <h2>Sequence Duplications within Labels</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('within-dup-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="within-dup-explanation" class="explanation-text">
                    <strong>Sequence Duplications within Labels</strong> shows how many sequences appear multiple times within each label. High duplication rates may indicate PCR artifacts, sequencing bias, or legitimate biological repeats. However, for machine learning and deep learning models, duplicate sequences can introduce bias during training. If the same sequence appears multiple times, the model may learn to overweight these repeated sequences, leading to skewed predictions and poor generalization.
                </div>
                <!-- This will be populated either with png plot showing duplication or with a message saying no duplications found -->
                {{sequence_duplications_within_classes}}
            </section>

            <section id="sequence-duplication-levels">
                <div class="sidebar-item">
                    {{icon_sequence_duplication_levels}}
                    <div class="section-header">
                        <h2>Duplicate Sequences between Labels</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('between-dup-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="between-dup-explanation" class="explanation-text">
                    <strong>Duplicate Sequences between Labels</strong> identifies sequences that appear in both labels. The same sequence annotated with both labels is indicating conflicting labels for identical sequence content.
                </div>
                <!-- This will be populated either with a table showing duplicate sequences or a message saying no duplications found -->
                {{sequence_duplication_levels}}
            </section>

            <section id="sequence-lengths">
                <div class="sidebar-item">
                    {{icon_sequence_lengths}}
                    <div class="section-header">
                        <h2>Sequence Lengths</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('lengths-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="lengths-explanation" class="explanation-text">
                    <strong>Sequence Lengths</strong> displays the distribution of sequence lengths in each label. The plot shows how lengths vary across your dataset. Significant differences in length distributions between labels may indicate bias or differences in the underlying biological processes. For ML models, length differences can sometimes be exploited as shortcuts.
                </div>

                <!-- This will be populated with png plot --->
                <img src="data:image/png;base64, {{sequence_length_plot_base64}}" alt="Sequence Lengths Plot" class="plot-half">
            </section>

            <section id="per-sequence-gc-content">
                <div class="sidebar-item">
                    {{icon_per_sequence_gc_content}}
                    <div class="section-header">
                        <h2>Per Sequence GC Content</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('gc-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="gc-explanation" class="explanation-text">
                    <strong>Per Sequence GC Content</strong> shows the distribution of GC% (percentage of G and C bases) across all sequences in each label. GC content affects DNA structure and stability. Significant differences in GC distribution between labels may indicate sequence composition bias that could impact model training.
                </div>
                <img src="data:image/png;base64, {{per-sequence-gc-content_base64}}" alt="Per Sequence GC Content" class="plot-half">
            </section>

            <section id="per-sequence-nucleotide-content">
                <div class="sidebar-item">
                    {{icon_per_sequence_nucleotide_content}}
                    <div class="section-header">
                        <h2>Per Sequence Nucleotide Content</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('nucleotide-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="nucleotide-explanation" class="explanation-text">
                    <strong>Per Sequence Nucleotide Content</strong> displays the distribution of individual nucleotide frequencies (A, C, G, T, N) across sequences. Each subplot shows how often that nucleotide appears in sequences from each label. Differences between labels may indicate composition bias or biological differences in your dataset.
                </div>
                {{per-sequence-nucleotide-content}}
            </section>

            <section id="per-sequence-dinucleotide-content">
                <div class="sidebar-item">
                    {{icon_per_sequence_dinucleotide_content}}
                    <div class="section-header">
                        <h2>Per Sequence Dinucleotide Content</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('dinucleotide-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="dinucleotide-explanation" class="explanation-text">
                    <strong>Per Sequence Dinucleotide Content</strong> shows the frequency of two-base combinations (e.g., AA, AC, AG, AT) across sequences. Dinucleotide frequencies can reveal sequence patterns and motifs. Each row shows all dinucleotides starting with a specific base. Significant differences between labels may indicate compositional bias.
                </div>
                {{per-sequence-dinucleotide-content}}
            </section>

            <section id="per-position-nucleotide-content">
                <div class="sidebar-item">
                    {{icon_per_position_nucleotide_content}}
                    <div class="section-header">
                        <h2>Per Position Nucleotide Content</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('per-position-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="per-position-explanation" class="explanation-text">
                    <strong>Per Position Nucleotide Content</strong> tracks nucleotide frequencies at each position along the sequence (5' to 3' direction). Each line shows one nucleotide's frequency across positions. Position-specific patterns can reveal adapter contamination, sequencing artifacts, or biological motifs. The bottom panel shows what proportion of sequences extend to each position.
                    {{position_window_note}}
                </div>
                {{per-position-nucleotide-content}}
            </section>

            <section id="per-position-reversed-nucleotide-content">
                <div class="sidebar-item">
                    {{icon_per_position_reversed_nucleotide_content}}
                    <div class="section-header">
                        <h2>Per Position Reversed Nucleotide Content</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('per-position-rev-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="per-position-rev-explanation" class="explanation-text">
                    <strong>Per Position Reversed Nucleotide Content</strong> is similar to the forward position plot, but reads sequences from 3' to 5' (reverse direction). This view helps identify patterns at sequence ends, which is particularly useful for detecting 3' adapter contamination or poly-A tails in RNA-seq data.
                    {{position_window_note_reversed}}
                </div>
                {{per-position-reversed-nucleotide-content}}
            </section>
        </div>
    </div>

{{report_scripts}}

</body>
</html>
"""

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
            f'reach it, and no analysed position has {MIN_SEQUENCES_PER_POSITION} sequences in both '
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
    """Describe the analysed per-position window and how well it is covered.

    The per-position plots and flags stop at a position chosen from the sequence
    lengths, and each position is scored only on the sequences that reach it. A
    reader who cannot see where the window ends, or how much of each class still
    stands behind its far end, cannot tell a thinly supported difference from a
    well supported one -- so both are stated. It goes inside the section's ?
    explanation rather than above the figure: it is background for reading the
    plot, not a finding, and on its own above the figure it read as a second
    explanation of the same check.
    """
    end_position = min(stats1.end_position, stats2.end_position)

    if end_position < 1:
        return ('<p>Sequences are too short for per-position '
                'statistics, so no positions were analysed.</p>')

    parts = []
    for stats in (stats1, stats2):
        label = stats.label if stats.label is not None else stats.filename
        coverage = stats.coverage_at(end_position)
        count = int(round(coverage * stats.stats['Number of sequences']))
        parts.append(f"{label}: {coverage:.1%} ({count:,} sequences)")

    thin = min(stats1.coverage_at(end_position), stats2.coverage_at(end_position)) < DEFAULT_MIN_COVERAGE
    caveat = ''
    if thin:
        caveat = (' Fewer sequences reach the end of this window than the '
                  f'{DEFAULT_MIN_COVERAGE:.0%} default, so its later positions rest on less data.')

    return (f'<p>Positions 1&ndash;{end_position} were analysed. '
            f'Each position is compared on the sequences that reach it; at position {end_position} '
            f'that is {parts[0]} and {parts[1]}.'
            f' A position is scored only where at least {MIN_SEQUENCES_PER_POSITION} sequences in '
            'each class reach it, and the difference it needs to show before it is flagged widens '
            'as its cohort shrinks, so that the worst case over the whole window is not set by '
            f'sampling noise.{caveat}</p>')


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
        tool_description = \
            """
            Toolkit for automated quality control of genomic datasets used in machine learning. 
            """
    if stats1.filepath == stats2.filepath:
        input_paths = stats1.filepath
    else:
        input_paths = f"{stats1.filepath}, {stats2.filepath}"
    generated_on = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_template = put_data(html_template, "{{tool_description}}", tool_description)
    html_template = put_data(html_template, "{{generated_on}}", generated_on)
    html_template = put_data(html_template, "{{input_paths}}", input_paths)
    html_template = put_data(html_template, "{{version}}", importlib.metadata.version("genomic-benchmarks-qc"))

    html_template = put_data(html_template, "{{filename1}}", stats1.filename)
    html_template = put_data(html_template, "{{filename2}}", stats2.filename)
    html_template = put_data(html_template, "{{label1}}", str(stats1.label) if stats1.label is not None else "N/A")
    html_template = put_data(html_template, "{{label2}}", str(stats2.label) if stats2.label is not None else "N/A")
    html_template = put_data(html_template, "{{seq_col1}}", str(stats1.seq_column) if stats1.seq_column is not None else "N/A")
    html_template = put_data(html_template, "{{seq_col2}}", str(stats2.seq_column) if stats2.seq_column is not None else "N/A")
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
    html_template = put_data(html_template, "{{unique_bases1}}", ', '.join(x for x in stats1.stats['Unique bases']))
    html_template = put_data(html_template, "{{unique_bases2}}", ', '.join(x for x in stats2.stats['Unique bases']))
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
    # 400-position window, which is what the zoom exists for. The PNG is still
    # written to the plots directory, it is just not what the page shows.
    per_position_payloads = per_position_payloads or {}
    for placeholder, direction, dom_id, name in (
        ("{{per-position-nucleotide-content}}", 'forward', 'ppv-fwd',
         'Per Position Nucleotide Content'),
        ("{{per-position-reversed-nucleotide-content}}", 'reversed', 'ppv-rev',
         'Per Position Reversed Nucleotide Content'),
    ):
        payload = per_position_payloads.get(direction)
        if payload is None:
            html_template = put_data(html_template, placeholder,
                                     f'<p class="no-plot-message">{disjoint_bases_message}</p>')
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
            html_template = put_data(html_template, "{{sequence_duplication_levels_file}}", f"All {len(duplicate_seqs)} duplicate sequences saved to {duplicate_seqs_file}.")
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
    html_template = put_data(html_template, "{{unique_bases_flags}}", unique_bases_flags)

    return html_template