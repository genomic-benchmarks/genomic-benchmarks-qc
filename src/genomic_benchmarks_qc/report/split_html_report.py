"""The HTML template for the train-test split report.

Self-contained in the same way as the class report: embedded plots, inlined
styling, and the top alignments rendered inline as text.

The page is built from the same components as the class report - the cards, the
key/value tables, and the collapsible findings panel with a listing table inside
it (.qc-panel / .qc-listing in report_design.css). The high-similarity pairs are
evidence for the Data Leakage check in the way the flagged positions are
evidence for a per-position check, so they sit in that check's card in the same
kind of panel, and not in a card of their own with no flag and no nav entry.
"""

import html
import logging
from datetime import datetime

from genomic_benchmarks_qc import __version__
from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.alignment_rendering import build_alignment_string
from genomic_benchmarks_qc.report.utils import (
    COMMON_CSS,
    LOGO_BASE64,
    REPORT_HEADER_HTML,
    encode_image_to_base64,
    icon_html,
    put_data,
)
from genomic_benchmarks_qc.utils.split_stats import flag_split_data_leakage

# Rows of the alignment listing that the page carries. The listing is evidence a
# reader spot-checks, not a data file - every hit is exported to
# mmseqs/mmseqs2_search_result.tsv - and each row here embeds two coloured
# sequences, which is what the page weighs.
ROW_CAP = 100

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Similar Sequences Report</title>
    <style>
    {{common_css}}
    </style>
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
                <span class="sidebar-spacer" style="display: inline-block; width: 40px"></span>
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
                <table>
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
                        Genomic Benchmarks QC evaluate-splits uses
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
                    <p>
                        The panel under the histogram lists the alignments at or above the threshold, up to the
                        first {{row_cap}}; every hit is exported to
                        <code>mmseqs/mmseqs2_search_result.tsv</code> beside this report. Its columns are:
                    </p>
                    <dl class="explanation-defs">
                        <dt>Query (Q)</dt>
                        <dd>Identifier of the test sequence: <code>seq_&lt;i&gt;_test</code>, where i is its 0-based position in the input test file.</dd>

                        <dt>Target (T)</dt>
                        <dd>Identifier of the train sequence: <code>seq_&lt;i&gt;_train</code>, where i is its 0-based position in the input train file.</dd>

                        <dt>Q Cov.</dt>
                        <dd>Fraction of the query sequence covered by the alignment (0–1).</dd>

                        <dt>T Cov.</dt>
                        <dd>Fraction of the target sequence covered by the alignment (0–1).</dd>

                        <dt>% Identity</dt>
                        <dd>Percent identical aligned positions in the aligned region.</dd>

                        <dt>% Similarity</dt>
                        <dd>Similarity score used for leakage detection, calculated as min(Q Cov., T Cov.) × % Identity.</dd>

                        <dt>E-value</dt>
                        <dd>MMseqs2 E-value (lower is more "significant").</dd>

                        <dt>Alignment</dt>
                        <dd>Click “Show” to expand the alignment visualisation for a row.</dd>
                    </dl>
                </div>
                <table>
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

                <img src="data:image/png;base64, {{histogram_similarity_base64}}" alt="Similarity Histogram" class="plot-wide">

                <details class="qc-panel" id="results-section">
                    <summary>{{alignments_count}}</summary>
                    <div class="qc-panel-body">
                        {{results_body}}
                    </div>
                </details>
            </section>

        </div>
    </div>

{{report_scripts}}

</body>
</html>
"""

RESULTS_TABLE = """<table class="qc-listing">
                            <thead>
                                <tr>
                                    <th>Query (Q)</th>
                                    <th>Target (T)</th>
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
                        </table>"""


def alignments_count_text(total, shown):
    """The panel's summary line: how many pairs leaked, and how many are listed.

    The same job the flagged-position panel's count does in the class report -
    say what is inside before it is opened - and it has to distinguish a clean
    split from a listing that was capped.
    """
    if total == 0:
        return 'No high-similarity alignments'
    plural = '' if total == 1 else 's'
    capped = f' (first {shown} shown)' if shown < total else ''
    return f'{total} high-similarity alignment{plural}{capped}'


def get_splits_html_template(basic_stats, threshold_stats, results_filt, plots_paths_dict,
                             tool_description=None, total_hits=None):
    """Build Train-Test Split Check HTML using shared helpers.

    Args:
        results_filt: The hits to list, already capped by the caller or not.
        total_hits: How many hits there were before any cap, so the panel can say
            what it is not showing. Defaults to the number of rows given.
    """

    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment; split_report.css is layered on top,
    # the way per_position_viewer.css is in the class report
    html_template = put_data(html_template, "{{common_css}}",
                             COMMON_CSS + assets.stylesheet('split_report.css'))
    # The shared behaviour, plus this report's own alignment toggles. At the end
    # of the body so it runs against a complete page.
    html_template = put_data(html_template, "{{report_scripts}}",
                             assets.script('report_ui.js', 'split_report.js'))
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
    html_template = put_data(html_template, "{{version}}", __version__)
    html_template = put_data(html_template, "{{row_cap}}", str(ROW_CAP))

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

    results_display = results_filt.head(ROW_CAP).copy()  # the page lists at most this many hits
    total = len(results_filt) if total_hits is None else total_hits
    html_template = put_data(html_template, "{{alignments_count}}",
                             alignments_count_text(total, len(results_display)))

    if len(results_display) == 0:
        # An empty panel says the same thing the flagged-position panel does when
        # a check found nothing: the listing is empty because there is nothing in
        # it, not because it failed to build.
        return put_data(
            html_template,
            "{{results_body}}",
            '<p class="qc-empty">No test sequence aligns to a train sequence at or above '
            'the similarity threshold.</p>'
        )

    rows = []
    alignment_error_logged = False

    for i, row in results_display.iterrows():
        try:
            alignment_str_color = build_alignment_string(row, color=True)
        except Exception as e:
            if not alignment_error_logged:
                alignment_error_logged = True
                logging.warning(
                    "Alignment visualisation failed for query=%s target=%s: %s: %s. "
                    "This is often caused by a conda-installed mmseqs2 build. See the "
                    "mmseqs2 installation notes in README.md — installing precompiled "
                    "binaries is recommended over compiling from source. "
                    "Further alignment errors in this report are not logged.",
                    row['query'], row['target'], type(e).__name__, e,
                )
            alignment_str_color = (
                '<span class="aln-error">ALIGNMENT VISUALISATION ERROR</span><br>'
                f'{html.escape(type(e).__name__)}: {html.escape(str(e))}<br>'
                '<span class="aln-error-detail">This can happen with a '
                'conda-installed mmseqs2 — see the installation notes in README.md.</span>'
            )
        rows.append(f"""
    <tr>
        <td class="qc-mono">{row['query']}</td>
        <td class="qc-mono">{row['target']}</td>
        <td class="qc-num">{row['qcov']:.2f}</td>
        <td class="qc-num">{row['tcov']:.2f}</td>
        <td class="qc-num">{row['pident']:.1f}</td>
        <td class="qc-key">{row['min_cov*pident']:.2f}</td>
        <td class="qc-num">{row['evalue']:.2e}</td>
        <td>
            <button type="button" class="qc-btn" onclick="toggleAlignment('aln-{i}', this)">
                Show
            </button>
        </td>
    </tr>

    <tr id="aln-{i}" style="display:none;">
        <td class="aln-cell" colspan="8">
            <pre class="alignment-block">{alignment_str_color}</pre>
        </td>
    </tr>
    """)

    results_table = put_data(RESULTS_TABLE, "{{results_rows}}", "\n".join(rows))
    return put_data(html_template, "{{results_body}}", results_table)
