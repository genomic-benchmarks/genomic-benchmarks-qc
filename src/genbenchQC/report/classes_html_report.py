from datetime import datetime
from genbenchQC.report.utils import put_data, encode_image_to_base64, escape_str, icon_html, COMMON_CSS, REPORT_HEADER_HTML, LOGO_BASE64
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

        nt = flag.split(" - ")[1]
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
                status_class = 'status-pass'
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
                <img src="{{logo_base64}}" alt="GenBenchQC Logo" style="max-width: 150px; height: auto;">
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
                <img src="data:image/png;base64, {{sequence_length_plot_base64}}" alt="Sequence Lengths Plot" style="max-width: 50%; height: auto; display: block; margin: 0 auto;">
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
                <img src="data:image/png;base64, {{per-sequence-gc-content_base64}}" alt="Per Sequence GC Content" style="max-width: 50%; height: auto; display: block; margin: 0 auto;">
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
                <img src="data:image/png;base64, {{per-sequence-nucleotide-content_base64}}" alt="Per Sequence Nucleotide Content" style="max-width: 100%; height: auto;">
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
                <img src="data:image/png;base64, {{per-sequence-dinucleotide-content_base64}}" alt="Per Sequence Dinucleotide Content" style="max-width: 100%; height: auto;">
            </section>

            <section id="per-position-nucleotide-content">
                <div class="sidebar-item">
                    {{icon_per_position_nucleotide_content}}
                    <div class="section-header">
                        <h2>Per Position Nucleotide Content</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('per-position-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div class="nucleotide-flags-container" id="per-position-nucleotide-flags">
                    {{per_position_nucleotide_content_flags}}
                </div>
                <div id="per-position-explanation" class="explanation-text">
                    <strong>Per Position Nucleotide Content</strong> tracks nucleotide frequencies at each position along the sequence (5' to 3' direction). Each line shows one nucleotide's frequency across positions. Position-specific patterns can reveal adapter contamination, sequencing artifacts, or biological motifs. The bottom panel shows what proportion of sequences extend to each position.
                </div>
                <img src="data:image/png;base64, {{per-position-nucleotide-content_base64}}" alt="Per Position Nucleotide Content" style="max-width: 108%; height: auto;">
            </section>

            <section id="per-position-reversed-nucleotide-content">
                <div class="sidebar-item">
                    {{icon_per_position_reversed_nucleotide_content}}
                    <div class="section-header">
                        <h2>Per Position Reversed Nucleotide Content</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('per-position-rev-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div class="nucleotide-flags-container" id="per-position-reversed-nucleotide-flags">
                    {{per_position_reversed_nucleotide_content_flags}}
                </div>
                <div id="per-position-rev-explanation" class="explanation-text">
                    <strong>Per Position Reversed Nucleotide Content</strong> is similar to the forward position plot, but reads sequences from 3' to 5' (reverse direction). This view helps identify patterns at sequence ends, which is particularly useful for detecting 3' adapter contamination or poly-A tails in RNA-seq data.
                </div>
                <img src="data:image/png;base64, {{per-position-reversed-nucleotide-content_base64}}" alt="Per Position Reversed Nucleotide Content" style="max-width: 108%; height: auto;">
            </section>
        </div>
    </div>

    <script>
        var sequenceDuplicationLevels = {{sequence_duplication_levels_seqs}};

        // Populate table for sequence duplication levels
        var tableBody = document.querySelector("#sequence-duplication-levels tbody");
        for (var i = 0; i < sequenceDuplicationLevels.length; i++) {
            var sequence = sequenceDuplicationLevels[i];

            var row = document.createElement("tr");
            var sequenceCell = document.createElement("td");

            sequenceCell.textContent = sequence;
            sequenceCell.className = "sequence_column";

            row.appendChild(sequenceCell);
            tableBody.appendChild(row);
        }
        
        // Toggle explanation visibility
        function toggleExplanation(elementId) {
            var element = document.getElementById(elementId);
            if (element.classList.contains('visible')) {
                element.classList.remove('visible');
            } else {
                element.classList.add('visible');
            }
        }
    </script>

</body>
</html>
"""

def get_dataset_html_template(stats1, stats2, plots_path, summary_statuses, duplicate_seqs, duplicate_seqs_file=None,
                             tool_description=None):
    """
    Returns the HTML template for the report.

    Args:
        stats1, stats2: objects containing dataset statistics (unchanged API).
        plots_path: dict with plot image paths.
        summary_statuses: dict with summary statuses for various checks.
        duplicate_seqs: list of duplicate sequences.
        duplicate_seqs_file: path to file with duplicate sequences (optional).
        tool_description: short description of the tool to include in the header (optional).

    Backwards compatible: if summary_statuses is None, placeholders are left empty.
    """
    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment
    html_template = put_data(html_template, "{{common_css}}", COMMON_CSS)
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
    html_template = put_data(html_template, "{{version}}", importlib.metadata.version("genbenchQC"))

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
                                 f'<img src="data:image/png;base64, {encode_image_to_base64(plots_path["Sequence Duplications within Labels"])}" alt="Sequence Duplications within Labels Plot" style="max-width: 100%; height: auto; display: block; margin: 0 auto;">')
    html_template = put_data(html_template, "{{sequence_length_plot_base64}}", 
                             encode_image_to_base64(plots_path['Sequence lengths']))
    html_template = put_data(html_template, "{{per-sequence-gc-content_base64}}", 
                             encode_image_to_base64(plots_path['Per sequence GC content']))
    html_template = put_data(html_template, "{{per-sequence-nucleotide-content_base64}}", 
                             encode_image_to_base64(plots_path['Per sequence nucleotide content']))
    html_template = put_data(html_template, "{{per-sequence-dinucleotide-content_base64}}", 
                             encode_image_to_base64(plots_path['Per sequence dinucleotide content']))
    html_template = put_data(html_template, "{{per-position-nucleotide-content_base64}}", 
                             encode_image_to_base64(plots_path['Per position nucleotide content']))
    html_template = put_data(html_template, "{{per-position-reversed-nucleotide-content_base64}}", encode_image_to_base64(
                             plots_path['Per position reversed nucleotide content']))

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
    html_template = put_data(html_template, "{{icon_per_position_reversed_nucleotide_content}}", icon_html(summary_statuses, 'Per reverse position nucleotide content'))
    html_template = put_data(html_template, "{{icon_per_sequence_gc_content}}", icon_html(summary_statuses, 'Per sequence GC content'))

    if duplicate_seqs == []:
        # no duplicate sequences found
        duplication_message = """
        <p>No duplicate sequences were found between classes.</p>
        """
        html_template = put_data(html_template, "{{sequence_duplication_levels}}", duplication_message)
        # set empty list for JS population
        html_template = put_data(html_template, "{{sequence_duplication_levels_seqs}}", "[]")
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
        html_template = put_data(html_template, "{{sequence_duplication_levels}}", duplication_table)
        # pass list of duplicate sequences for JS population
        escaped_seqs = [escape_str(seq) for seq in duplicate_seqs[:10]]
        html_template = put_data(html_template, "{{sequence_duplication_levels_seqs}}",
                                 '[' + ', '.join(f'{seq}' for seq in escaped_seqs) + ']')
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

    # Generate nucleotide-level flags for per-position sections
    per_position_nucleotide_flags = generate_nucleotide_flags_html(summary_statuses, 'Per position nucleotide content')
    html_template = put_data(html_template, "{{per_position_nucleotide_content_flags}}", per_position_nucleotide_flags)
    
    per_position_reversed_nucleotide_flags = generate_nucleotide_flags_html(summary_statuses, 'Per reverse position nucleotide content')
    html_template = put_data(html_template, "{{per_position_reversed_nucleotide_content_flags}}", per_position_reversed_nucleotide_flags)

    # Generate unique bases flags
    unique_bases_flags = generate_nucleotide_flags_html(summary_statuses, 'Unique bases')
    html_template = put_data(html_template, "{{unique_bases_flags}}", unique_bases_flags)

    return html_template