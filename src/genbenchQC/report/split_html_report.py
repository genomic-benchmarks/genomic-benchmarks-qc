from datetime import datetime
import html
from genbenchQC.report.utils import encode_image_to_base64, put_data, icon_html, COMMON_CSS, REPORT_HEADER_HTML, LOGO_BASE64
from genbenchQC.report.alignment_rendering import build_alignment_string
from genbenchQC.utils.split_stats import flag_split_data_leakage
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
                    <tr id="leakage-filename">
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

            <section id="similarity-section" class="table-section">
                <div class="sidebar-item">
                    {{icon_data_leakage}}
                    <div class="section-header">
                        <h2>Data Leakage</h2>
                        <button class="toggle-btn" onclick="toggleExplanation('data-leakage-explanation')" title="Show explanation">?</button>
                    </div>
                </div>
                <div id="data-leakage-explanation" class="explanation-text">
                    <p>
                        GenBenchQC evaluate-splits uses 
                        <a href="https://github.com/soedinglab/MMseqs2" target="_blank">MMseqs2</a> 
                        to perform an ultra fast and sensitive test sequence search against a train set database and compute alignment-based 
                        metrics for detecting data leakage across train–test splits.
                    </p>
                    <p>
                        Data leakage is the percentage of train/test sequences whose best alignment exceeds the configured similarity threshold. 
                        Here, similarity equals min(query coverage, target coverage) × percent identity, 
                        where coverage is the fraction of each sequence spanned by the alignment, 
                        and percent identity is the proportion of identical aligned positions.
                    </p>
                    <p>
                        Status is based on the percentage of test sequences with similarity at or above the threshold: 
                        Pass means 0% leaked test sequences, Warning means greater than 0% and less than 2%, 
                        and Fail means 2% or more.
                    </p>
                </div>
                <table style="width: 100%; border-collapse: collapse; margin: 20px 0;">
                    <tr id="filename">
                        <td><span>Filename</span></td>
                        <td style="text-align: center;">{{test_filename}}</td>
                        <td style="text-align: center;">{{train_filename}}</td>
                    </tr>
                    <tr id="leakage-percentage">
                        <td><span>Data Leakage (percentage)</span></td>
                        <td style="text-align: center;">{{perc_queries_above_thr}}%</td>
                        <td style="text-align: center;">{{perc_targets_above_thr}}%</td>
                    </tr>
                    <tr id="leakage-count">
                        <td><span>Data Leakage (count)</span></td>
                        <td style="text-align: center;">{{num_queries_above_thr}}</td>
                        <td style="text-align: center;">{{num_targets_above_thr}}</td>
                    </tr>
                </table> 

                <img src="data:image/png;base64, {{histogram_similarity_base64}}" alt="Similarity Histogram" style="max-width: 100%; width: 100%; height: auto; display: block; margin: 0 auto;">
            </section>

            <section id="results-section">
                <div class="section-header">
                    <h3 style="margin: 0;">Alignment of High-Similarity Sequences (Top 100)</h3>
                    <button class="toggle-btn"
                            onclick="toggleExplanation('results-explanation')"
                            title="Show explanation">?</button>
                </div>

                <div id="results-explanation" class="explanation-text">
                    <p>
                        This table lists the top 100 alignments between test (query) and train (target)
                        sequences that exceed the defined similarity threshold.
                    </p>

                    <div style="
                        display: grid;
                        grid-template-columns: 180px 1fr;
                        gap: 4px 14px;
                        margin-top: 10px;
                        margin-left: 20px;
                    ">
                        <div><strong>Query</strong></div>
                        <div>0-based index of the test sequence in the input test file order.</div>

                        <div><strong>Target</strong></div>
                        <div>0-based index of the train sequence in the input train file order.</div>

                        <div><strong>Q Cov.</strong></div>
                        <div>Fraction of the query sequence covered by the alignment (0–1).</div>

                        <div><strong>T Cov.</strong></div>
                        <div>Fraction of the target sequence covered by the alignment (0–1).</div>

                        <div><strong>% Identity</strong></div>
                        <div>Percent identical aligned positions in the aligned region.</div>

                        <div><strong>% Similarity</strong></div>
                        <div>Similarity score used for leakage detection, calculated as min(Q Cov., T Cov.) × % Identity.</div>

                        <div><strong>E-value</strong></div>
                        <div>MMSeqs2 E-value (lower is more "significant").</div>

                        <div><strong>Alignment</strong></div>
                        <div>Click “Show” to expand the alignment visualisation for a row.</div>
                    </div>
                </div>

                <table style="width: 100%; border-collapse: collapse;">
                    <thead>
                        <tr>
                            <th style="text-align: left;">Query (Q)</th>
                            <th style="text-align: left;">Target (T)</th>
                            <th>Q Cov.</th>
                            <th>T Cov.</th>
                            <th>% Identity</th>
                            <th>% Similarity</th>
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
    html_template = put_data(html_template, "{{perc_queries_above_thr}}", f"{threshold_stats['perc_queries_above_thr']:.2f}")
    html_template = put_data(html_template, "{{num_queries_above_thr}}", f"{threshold_stats['num_queries_above_thr']}")
    html_template = put_data(html_template, "{{perc_targets_above_thr}}", f"{threshold_stats['perc_targets_above_thr']:.2f}")
    html_template = put_data(html_template, "{{num_targets_above_thr}}", f"{threshold_stats['num_targets_above_thr']}")

    # insert plots as base64-encoded images
    html_template = put_data(html_template, "{{histogram_similarity_base64}}", 
                             encode_image_to_base64(plots_paths_dict['Similarity histograms']))

    leakage_flag = flag_split_data_leakage(threshold_stats["perc_queries_above_thr"])
    html_template = put_data(
        html_template,
        "{{icon_data_leakage}}",
        icon_html({"Data Leakage": leakage_flag}, "Data Leakage")
    )

    if len(results_filt) == 0:
        html_template = put_data(
            html_template,
            "{{results_rows}}",
            "<tr><td colspan='8'>No high-similarity sequences found across train/test splits.</td></tr>"
        )
        return html_template

    results_display = results_filt.head(100).copy() # limit to top 100 hits for display in HTML report

    rows = []

    for i, row in results_display.iterrows():
        try: 
            alignment_str_color = build_alignment_string(row, color=True)
        except Exception as e:
            alignment_str_color = (
                '<span style="color:red; font-weight:bold;">'
                'ALIGNMENT VISUALISATION ERROR</span><br>'
                f'{html.escape(type(e).__name__)}: {html.escape(str(e))}'
            )
        rows.append(f"""
    <tr>
        <td>{row['query']}</td>
        <td>{row['target']}</td>
        <td style="text-align: center; ">{row['qcov']:.2f}</td>
        <td style="text-align: center; ">{row['tcov']:.2f}</td>
        <td style="text-align: center; ">{row['pident']:.1f}</td>
        <td style="text-align: center; ">{row['min_cov*pident']:.2f}</td>
        <td style="text-align: center; ">{row['evalue']:.2e}</td>
        <td style="text-align: center;">
            <button onclick="toggleAlignment('aln-{i}', this)">
                Show
            </button>
        </td>
    </tr>

    <tr id="aln-{i}" style="display:none;">
        <td colspan="8">
            <pre class="alignment-block">{alignment_str_color}</pre>
        </td>
    </tr>
    """)
        
    html_template = put_data(html_template, "{{results_rows}}", "\n".join(rows))
    return html_template
