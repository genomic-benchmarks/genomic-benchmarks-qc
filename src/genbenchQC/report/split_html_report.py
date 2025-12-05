from datetime import datetime
from genbenchQC.report.report_common import put_data, put_file_details, COMMON_CSS, REPORT_HEADER_HTML
import importlib.metadata

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
            <div class="sidebar-item"><a href="#clusters-section">Clusters</a></div>
        </div>

        <div class="content">
            {{report_header}}

            <section id="basic-descriptive-statistics">
                <h2>Basic Descriptive Statistics</h2>
                <div class="data-item"><span>Train set filename:</span> {{train_filename}}</div>
                <div class="data-item"><span>Number of sequences in train set:</span> {{number_of_sequences_train}}</div>
                <div class="data-item"><span>Number of train sequences overlapping with test set:</span> {{train_overlap}}</div>
                <div class="data-item"><span>Test set filename:</span> {{test_filename}}</div>
                <div class="data-item"><span>Number of sequences in test set:</span> {{number_of_sequences_test}}</div>
                <div class="data-item"><span>Number of test sequences overlapping with train set:</span> {{test_overlap}}</div>
            </section>

            <section id="clusters-section">
                <h2>Clusters of Similar Sequences</h2>
                <p>Clustering was done using cd-hit est-2d with identity threshold of {{identity_threshold}} and sequence alignment coverage of {{alignment_coverage}}.</p>
                {{clusters}}
            </section>

        </div>
    </div>

</body>
</html>
"""

def get_train_test_html_template(clusters, filename_train, sequences_train, filename_test, sequences_test, identity_threshold, alignment_coverage, tool_description=None):
    """Build train/test similarity HTML using shared helpers."""

    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment
    html_template = put_data(html_template, "{{common_css}}", COMMON_CSS)
    html_template = put_data(html_template, "{{report_header}}", REPORT_HEADER_HTML)

    # header info
    if tool_description is None:
        tool_description = "Toolkit for automated quality control of genomic datasets used in machine learning."
    input_paths = f"{filename_train}, {filename_test}" if filename_train != filename_test else filename_train
    generated_on = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_template = put_data(html_template, "{{tool_description}}", tool_description)
    html_template = put_data(html_template, "{{generated_on}}", generated_on)
    html_template = put_data(html_template, "{{input_paths}}", input_paths)
    html_template = put_data(html_template, "{{version}}", importlib.metadata.version("genbenchQC"))

    html_template = put_data(html_template, "{{train_filename}}", str(filename_train))
    html_template = put_data(html_template, "{{test_filename}}", str(filename_test))
    html_template = put_data(html_template, "{{number_of_sequences_train}}", str(len(sequences_train)))
    html_template = put_data(html_template, "{{number_of_sequences_test}}", str(len(sequences_test)))
    train_overlap = sum(len(cluster.get('train', [])) for cluster in clusters)
    test_overlap = sum(len(cluster.get('test', [])) for cluster in clusters)
    html_template = put_data(html_template, "{{train_overlap}}", str(train_overlap))
    html_template = put_data(html_template, "{{test_overlap}}", str(test_overlap))
    html_template = put_data(html_template, "{{identity_threshold}}", str(identity_threshold))
    html_template = put_data(html_template, "{{alignment_coverage}}", str(alignment_coverage))

    if not clusters:
        return html_template.replace("{{clusters}}", "<h2>No similar sequences found.</h2>")

    cluster_blocks = []
    max_seq_display = 1000
    n_sequences = 0

    for cluster in clusters:
        train_sequences = cluster.get('train', [])
        test_sequences = cluster.get('test', [])
        if n_sequences + len(train_sequences) > max_seq_display:
            n_sequences += len(train_sequences)
            train_sequences = train_sequences[:2] + ["..."] if len(train_sequences) > 2 else train_sequences
            test_sequences = test_sequences[:2] + ["..."] if len(test_sequences) > 2 else test_sequences
        elif n_sequences + len(train_sequences) + len(test_sequences) > max_seq_display:
            n_sequences += len(train_sequences) + len(test_sequences)
            test_sequences = test_sequences[:2] + ["..."] if len(test_sequences) > 2 else test_sequences
        else:
            n_sequences += len(train_sequences) + len(test_sequences)

        cluster_html = f"""
        <div class="cluster">
            <h3>Cluster #{cluster['cluster']}</h3>
            <div class="section-title">Train Sequences:</div>
            <pre>{chr(10).join(train_sequences)}</pre>
            <div class="section-title">Test Sequences:</div>
            <pre>{chr(10).join(test_sequences)}</pre>
        </div>
        """
        cluster_blocks.append(cluster_html)
        if n_sequences >= max_seq_display:
            cluster_blocks = [f"<b>Note:</b> There are too many clusters to display ({len(clusters)} clusters). Showing only part of the sequences. If you want to access all the clusters, toggle 'json' format and refer to the json report.  <br><br/>"] + cluster_blocks
            break

    return html_template.replace("{{clusters}}", "\n".join(cluster_blocks))