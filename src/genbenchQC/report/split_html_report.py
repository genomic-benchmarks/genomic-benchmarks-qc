from datetime import datetime
from genbenchQC import __version__
from genbenchQC.report.report_common import put_data, COMMON_CSS, REPORT_HEADER_HTML

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Similar Sequences Report</title>
    <style>{{common_css}}</style>
</head>
<body>
    <div class="container">
        <div class="sidebar">
            <div class="logo" style="text-align: center; margin-bottom: 20px;">
                <img src="https://raw.githubusercontent.com/katarinagresova/GenBenchQC/main/assets/logo_with_text_transparent_small.png" alt="GenBenchQC Logo" style="max-width: 150px; height: auto;">
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
                        <td style="text-align: center;">{{train_filename}}</td>
                        <td style="text-align: center;">{{test_filename}}</td>
                    </tr>
                    <tr id="num-sequences">
                        <td><span>Number of sequences</span></td>
                        <td style="text-align: center;">{{number_of_sequences_train}}</td>
                        <td style="text-align: center;">{{number_of_sequences_test}}</td>
                    </tr>
                    <tr id="min-length">
                        <td><span>Minimum length</span></td>
                        <td style="text-align: center;">{{min_length_train}}</td>
                        <td style="text-align: center;">{{min_length_test}}</td>
                    </tr>
                    <tr id="mean-length">
                        <td><span>Mean length</span></td>
                        <td style="text-align: center;">{{mean_length_train}}</td>
                        <td style="text-align: center;">{{mean_length_test}}</td>
                    </tr>
                    <tr id="max-length">
                        <td><span>Maximum length</span></td>
                        <td style="text-align: center;">{{max_length_train}}</td>
                        <td style="text-align: center;">{{max_length_test}}</td>
                    </tr>
                </table>            
            </section>

            <section id="similarity-section">
                <h2>Train/Test Similarity</h2>
                <p>The test set was stratified using hashFrag. The score is the maximum pairwise SW alignment score of a test sequence when queried against all train set sequences. More information about the score is provided by <a href="https://github.com/de-Boer-Lab/hashFrag?tab=readme-ov-file#a-minor-note-on-the-use-of-blast-alignment-scores-in-hashfrag" target="_blank" rel="noopener noreferrer">hashFrag</a>. </p>
                <img src={{max_corr_alignment_scores_plot}} alt="Max Corrected Alignment Scores Plot" style="max-width: 50%; height: auto; display: block; margin: 0 auto;">          
                <p>Threshold: {{threshold}} </p>
                <p>Percentage of alignments above threshold: {{perc_above_threshold}} </p>
                <p>Percentage of alignments below threshold: {{perc_below_threshold}} </p>
                <img src={{threshold_selection_plot}} alt="Threshold Selection Plot" style="max-width: 50%; height: auto; display: block; margin: 0 auto;">    
        </section>

        </div>
    </div>

</body>
</html>
"""

def get_split_html_template(basic_stats, threshold_stats, plots_paths, tool_description=None):
    """Build train/test similarity HTML using shared helpers."""

    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment
    html_template = put_data(html_template, "{{common_css}}", COMMON_CSS)
    html_template = put_data(html_template, "{{report_header}}", REPORT_HEADER_HTML)

    # header info
    if tool_description is None:
        tool_description = "Toolkit for automated quality control of genomic datasets used in machine learning."
    input_paths = f"{basic_stats['train_filename']}, {basic_stats['test_filename']}" if basic_stats['train_filename'] != basic_stats['test_filename'] else basic_stats['train_filename']
    generated_on = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_template = put_data(html_template, "{{tool_description}}", tool_description)
    html_template = put_data(html_template, "{{generated_on}}", generated_on)
    html_template = put_data(html_template, "{{input_paths}}", input_paths)
    html_template = put_data(html_template, "{{version}}", __version__)

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

    html_template = put_data(html_template, "{{max_corr_alignment_scores_plot}}", str(plots_paths['Distribution of test set maximum corrected alignment scores']))
    html_template = put_data(html_template, "{{threshold}}", str(threshold_stats["threshold"]))
    html_template = put_data(html_template, "{{perc_above_threshold}}", str(round(threshold_stats["num_above_threshold"] / threshold_stats["total_alignments"] * 100, 2)))
    html_template = put_data(html_template, "{{perc_below_threshold}}", str(round(threshold_stats["num_below_threshold"] / threshold_stats["total_alignments"] * 100, 2)))

    html_template = put_data(html_template, "{{threshold_selection_plot}}", str(plots_paths['Threshold selection plot']))

    return html_template