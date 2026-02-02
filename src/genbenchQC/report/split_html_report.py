from datetime import datetime
from genbenchQC.report.utils import put_data, put_file_details, COMMON_CSS, REPORT_HEADER_HTML, LOGO_BASE64
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
            <div class="sidebar-item"><a href="#basic-descriptive-statistics">Basic Descriptive Statistics</a></div>
            <div class="sidebar-item"><a href="#similarity-section">Train/Test Similarity</a></div>
        </div>

        <div class="content">
            {{report_header}}

            <section id="basic-descriptive-statistics">
                <h2>Basic Descriptive Statistics</h2>             
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
                <h2>Train/Test Similarity</h2>
                <p>Info about run with MMseqs2. </p>
                <p>Coverage threshold: {{coverage_threshold}} </p>
                <p>Percentage of hits above thresholds: {{perc_above_threshold}} </p>
                <p>Percentage of hits below thresholds: {{perc_below_threshold}} </p>
                <img src={{histogram_coverage}} alt="Histogram of Coverage" style="max-width: 50%; height: auto; display: block; margin: 0 auto;">  
            </section>

            <section id="results-section">
                <h3>Alignment Results</h3>
                <table style="width: 100%; border-collapse: collapse;">
                    <thead>
                        <tr>
                            <th>Query (Q)</th>
                            <th>Target (T)</th>
                            <th>Q Coverage</th>
                            <th>T Coverage</th>
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
</script>

</body>
</html>
"""

def get_splits_html_template(basic_stats, results_filt_aln, coverage_threshold, plots_paths_dict, tool_description=None):
    """Build train/test similarity HTML using shared helpers."""

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

    html_template = put_data(html_template, "{{train_filename}}", str(basic_stats["train_filename"]))
    html_template = put_data(html_template, "{{test_filename}}", str(basic_stats["test_filename"]))
    html_template = put_data(html_template, "{{number_of_sequences_train}}", str(basic_stats["number_of_sequences_train"]))
    html_template = put_data(html_template, "{{number_of_sequences_test}}", str(basic_stats["number_of_sequences_test"]))
    html_template = put_data(html_template, "{{min_length_train}}", str(basic_stats["min_length_train"]))
    html_template = put_data(html_template, "{{mean_length_train}}", str(basic_stats["mean_length_train"]))
    html_template = put_data(html_template, "{{max_length_train}}", str(basic_stats["max_length_train"]))
    html_template = put_data(html_template, "{{min_length_test}}", str(basic_stats["min_length_test"]))
    html_template = put_data(html_template, "{{mean_length_test}}", str(basic_stats["mean_length_test"]))
    html_template = put_data(html_template, "{{max_length_test}}", str(basic_stats["max_length_test"]))
    html_template = put_data(html_template, "{{coverage_threshold}}", coverage_threshold)
    html_template = put_data(html_template, "{{histogram_coverage}}", str(plots_paths_dict['Histogram of coverage']))

    if len(results_filt_aln) == 0:
        html_template = put_data(
            html_template,
            "{{results_rows}}",
            "<tr><td colspan='5'>No similar sequences found.</td></tr>"
        )
        return html_template

    rows = []

    for i, row in results_filt_aln.iterrows():
        rows.append(f"""
    <tr>
        <td>{row['query']}</td>
        <td>{row['target']}</td>
        <td style="text-align: center;">{row['qcov']:.3f}</td>
        <td style="text-align: center;">{row['tcov']:.3f}</td>

        <td style="text-align: center;">
            <button onclick="toggleAlignment('aln-{i}', this)">
                Show
            </button>
        </td>
    </tr>

    <tr id="aln-{i}" style="display:none;">
        <td colspan="5">
            <pre class="alignment-block">{row['alignment_str']}</pre>
        </td>
    </tr>
    """)
        
    html_template = put_data(html_template, "{{results_rows}}", "\n".join(rows))
    return html_template