from datetime import datetime
from genbenchQC.report.utils import put_data, put_file_details, COMMON_CSS, REPORT_HEADER_HTML, LOGO_BASE64
from genbenchQC.utils.data_leakage_utils import add_alignments_to_results, build_alignment_string
import importlib.metadata

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Similar Sequences Report</title>
    <style>
    {{common_css}}

    /* Local styles for alignment display only */
    .alignment-block {
        font-family: monospace;
        font-size: 12px;
        background: #f7f7f7;
        padding: 10px;
        white-space: pre;
        overflow-x: auto;
        border: 1px solid #ddd;
    }
    /* DNA base coloring */
    .base-A { color: #2ca02c; font-weight: bold; } /* green */
    .base-C { color: #1f77b4; font-weight: bold; } /* blue */
    .base-G { color: #ff7f0e; font-weight: bold; } /* orange */
    .base-T { color: #d62728; font-weight: bold; } /* red */

    .base-gap { color: #555; } /* darker grey */
    .base-other { color: #999; }  /* N, ambiguous - grey */
    </style>
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
                <span style="display: inline-block; width: 40px"></span>
                <a href="#basic-descriptive-statistics">Basic Descriptive Statistics</a>
            </div>
            <div class="sidebar-item">
                {{icon_data_leakage}}
                <a href="#similarity-section">Data Leakage</a>
            </div>
        </div>

        <div class="content">
            {{report_header}}

            <section id="basic-descriptive-statistics" class="table-section">
                <div class="section-header">
                    <h2>Basic Descriptive Statistics</h2>
                </div>            
                <table style="width: 100%; border-collapse: collapse; margin: 20px 0;">
                    <tr id="filename">
                        <td><span>Filename</span></td>
                        <td style="text-align: center;">{{test_filename}}</td>
                        <td style="text-align: center;">{{train_filename}}</td>
                    </tr>
                    <tr id="num-sequences">
                        <td><span>Number of sequences</span></td>
                        <td style="text-align: center;">{{number_of_sequences_test}}</td>
                        <td style="text-align: center;">{{number_of_sequences_train}}</td>
                    </tr>
                    <tr id="min-length">
                        <td><span>Minimum length</span></td>
                        <td style="text-align: center;">{{min_length_test}}</td>
                        <td style="text-align: center;">{{min_length_train}}</td>
                    </tr>
                    <tr id="mean-length">
                        <td><span>Mean length</span></td>
                        <td style="text-align: center;">{{mean_length_test}}</td>
                        <td style="text-align: center;">{{mean_length_train}}</td>
                    </tr>
                    <tr id="max-length">
                        <td><span>Maximum length</span></td>
                        <td style="text-align: center;">{{max_length_test}}</td>
                        <td style="text-align: center;">{{max_length_train}}</td>
                    </tr>
                </table>            
            </section>

            <section id="similarity-section">
                <div class="sidebar-item">
                    {{icon_data_leakage}}
                    <div class="section-header">
                        <h2>Data Leakage</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('data-leakage-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="data-leakage-explanation" class="explanation-text">
                    Using MMSeqs2, test set sequencies (queries) are aligned against train set sequences (target database), to identify similar sequences across splits.
                    Query and target sequences are 0-based indexed, in the order they appear in the input files. Refer to fasta files for index-to-sequence mapping. 
                    Data leakage is defined as the percentage of sequences in the train/test sets that have exceeded a set coverage threshold.
                    Coverage is the sequence length overlap. The alignment covers at least this threshold of the query sequence and of the target sequence.
                </div>
                <table style="width: 49%; border-collapse: collapse; margin: 20px 0;">
                    <tr>
                        <td><span>Coverage threshold</span></td>
                        <td style="text-align: right;">{{coverage_threshold}}</td>
                    </tr>
                    <tr style="height: 15px;">
                        <td></td>
                        <td></td>
                    </tr>
                    <tr>
                        <td><span>Data leakage in train set (%)</span></td>
                        <td style="text-align: right;">{{perc_targets_above_thr}}</td>
                    </tr>                    
                    <tr>
                        <td><span>Data leakage in test set (%)</span></td>
                        <td style="text-align: right;">{{perc_queries_above_thr}}</td>
                    </tr>
                </table>
                <img src={{histogram_coverage}} alt="Histogram of Coverage" style="max-width: 50%; height: auto; display: block; margin: 0 auto;">  
            </section>

            <section id="results-section">
                <h3>Alignment of Leaked Sequences (Top 100)</h3>
                <table style="width: 100%; border-collapse: collapse;">
                    <thead>
                        <tr>
                            <th>Query (Q)</th>
                            <th>Target (T)</th>
                            <th>Q Cov.</th>
                            <th>T Cov.</th>
                            <th>Perc. Identity</th>
                            <th>E-value</th>
                            <th>Alignment</th>
                        </tr>
                    </thead>
                    <tbody>
                        {{results_rows}}
                    </tbody>
                </table>
            </section>
            
        </div>
    </div>

    <script>
        function toggleAlignment(id, btn) {
            const row = document.getElementById(id);
            if (!row) return;

            if (row.style.display === "none") {
                row.style.display = "table-row";
                btn.textContent = "Hide";
            } else {
                row.style.display = "none";
                btn.textContent = "Show";
            }
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

def get_splits_html_template(basic_stats, threshold_stats, results_filt, plots_paths_dict, tool_description=None):
    """Build Train-Test Split Check HTML using shared helpers."""

    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment
    html_template = put_data(html_template, "{{common_css}}", COMMON_CSS)
    html_template = put_data(html_template, "{{report_header}}", REPORT_HEADER_HTML)

    # insert logo
    html_template = put_data(html_template, "{{logo_base64}}", LOGO_BASE64)

    # header info
    if tool_description is None:
        tool_description = "Toolkit for automated quality control of genomic datasets used in machine learning."
    input_paths = f"{basic_stats['train_filename']}, {basic_stats['test_filename']}" if basic_stats['train_filename'] != basic_stats['test_filename'] else basic_stats['train_filename']
    generated_on = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_template = put_data(html_template, "{{tool_description}}", tool_description)
    html_template = put_data(html_template, "{{generated_on}}", generated_on)
    html_template = put_data(html_template, "{{input_paths}}", input_paths)
    html_template = put_data(html_template, "{{version}}", importlib.metadata.version("genbenchQC"))

    html_template = put_data(html_template, "{{train_filename}}", str(basic_stats['train_filename']))
    html_template = put_data(html_template, "{{test_filename}}", str(basic_stats['test_filename']))
    html_template = put_data(html_template, "{{number_of_sequences_train}}", str(basic_stats['number_of_sequences_train']))
    html_template = put_data(html_template, "{{number_of_sequences_test}}", str(basic_stats['number_of_sequences_test']))
    html_template = put_data(html_template, "{{min_length_train}}", str(basic_stats['min_length_train']))
    html_template = put_data(html_template, "{{mean_length_train}}", f"{basic_stats['mean_length_train']:.2f}")
    html_template = put_data(html_template, "{{max_length_train}}", str(basic_stats['max_length_train']))
    html_template = put_data(html_template, "{{min_length_test}}", str(basic_stats['min_length_test']))
    html_template = put_data(html_template, "{{mean_length_test}}", f"{basic_stats['mean_length_test']:.2f}")
    html_template = put_data(html_template, "{{max_length_test}}", str(basic_stats['max_length_test']))
    html_template = put_data(html_template, "{{coverage_threshold}}", f"{threshold_stats['coverage_threshold']:.2f}")
    html_template = put_data(html_template, "{{perc_queries_above_thr}}", f"{threshold_stats['perc_queries_above_thr']:.1f}%")
    html_template = put_data(html_template, "{{perc_targets_above_thr}}", f"{threshold_stats['perc_targets_above_thr']:.1f}%")
    html_template = put_data(html_template, "{{histogram_coverage}}", str(plots_paths_dict['Coverage histograms']))

    if len(results_filt) == 0:
        html_template = put_data(
            html_template,
            "{{results_rows}}",
            "<tr><td colspan='7'>No leaked sequences found across train/test splits.</td></tr>"
        )
        html_template = put_data(html_template, "{{icon_data_leakage}}", '<span class="status-icon status-pass">✔</span>')
        return html_template

    html_template = put_data(html_template, "{{icon_data_leakage}}", '<span class="status-icon status-fail">✘</span>')
    results_display = results_filt.head(100).copy() # limit to top 100 hits for display in HTML report
    results_display = add_alignments_to_results(results_display)

    rows = []

    for i, row in results_display.iterrows():
        alignment_str_color = build_alignment_string(row, color=True)

        rows.append(f"""
    <tr>
        <td>{row['query']}</td>
        <td>{row['target']}</td>
        <td style="text-align: center; ">{row['qcov']:.2f}</td>
        <td style="text-align: center; ">{row['tcov']:.2f}</td>
        <td style="text-align: right; padding-right: 20px;">{row['pident']:.1f}%</td>
        <td style="text-align: left; ">{row['evalue']}</td>
        <td style="text-align: center;">
            <button onclick="toggleAlignment('aln-{i}', this)">
                Show
            </button>
        </td>
    </tr>

    <tr id="aln-{i}" style="display:none;">
        <td colspan="7">
            <pre class="alignment-block">{alignment_str_color}</pre>
        </td>
    </tr>
    """)
        
    html_template = put_data(html_template, "{{results_rows}}", "\n".join(rows))
    return html_template