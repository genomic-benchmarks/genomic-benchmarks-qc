from datetime import datetime
from genbenchQC.report.report_common import put_data, put_file_details, COMMON_CSS, REPORT_HEADER_HTML
import importlib.metadata

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
                <img src="https://raw.githubusercontent.com/katarinagresova/GenBenchQC/main/assets/logo_with_text_transparent_small.png" alt="GenBenchQC Logo" style="max-width: 150px; height: auto;">
                <span style="display: block; font-size: 14px; color: #555;">v{{version}}</span>
            </div>
            <h2>Summary</h2>
            <div class="sidebar-item"><a href="#basic-descriptive-statistics">Basic Descriptive Statistics</a></div>
            <div class="sidebar-item"><a href="#sequence-lengths">Sequence lengths</a></div>
            <div class="sidebar-item"><a href="#sequence-duplication-levels">Sequence duplication levels</a></div>
            <div class="sidebar-item"><a href="#per-sequence-nucleotide-content">Per Sequence Nucleotide Content</a></div>
            <div class="sidebar-item"><a href="#per-sequence-dinucleotide-content">Per Sequence Dinucleotide Content</a></div>
            <div class="sidebar-item"><a href="#per-position-nucleotide-content">Per Position Nucleotide Content</a></div>
            <div class="sidebar-item"><a href="#per-position-reversed-nucleotide-content">Per Position Reversed Nucleotide Content</a></div>
            <a href="#per-sequence-gc-content">Per Sequence GC Content</a>
        </div>

        <div class="content">
            <!-- Report header: logo, short description and generated-on/data source info -->
            {{report_header}}

            <section id="basic-descriptive-statistics">
                <h2>Basic Descriptive Statistics</h2>
                <div class="data-item" id="filename">
                    <span>Filename:</span> {{filename}} <!-- Filename will be displayed here -->
                </div>
                <div class="data-item" id="label">
                    <span>Label:</span> {{label}} <!-- Label will be displayed here -->
                </div>
                <div class="data-item" id="seq_column">
                    <span>Sequence column:</span> {{seq_column}} <!-- Sequence column will be displayed here -->
                </div>
                <div class="data-item" id="num-sequences">
                    <span>Number of sequences:</span> {{number_of_sequences}} <!-- Number of sequences will be displayed here -->
                </div>
                <div class="data-item" id="dedup-sequences">
                    <span>Unique sequences:</span> {{dedup_sequences}} <!-- Number of sequences left after deduplication will be displayed here -->
                </div>
                <div class="data-item" id="min-length">
                    <span>Minimum length:</span> {{min_length}} <!-- Minimum length will be displayed here -->
                </div>
                <div class="data-item" id="mean-length">
                    <span>Mean length:</span> {{mean_length}} <!-- Mean length will be displayed here -->
                </div>
                <div class="data-item" id="max-length">
                    <span>Maximum length:</span> {{max_length}} <!-- Maximum length will be displayed here -->
                </div>
                <div class="data-item" id="num-bases">
                    <span>Number of bases:</span> {{number_of_bases}} <!-- Number of bases will be displayed here -->
                </div>
                <div class="data-item" id="unique-bases">
                    <span>Unique bases:</span> {{unique_bases}} <!-- Unique bases will be displayed here -->
                </div>
                <div class="data-item" id="gc-content">
                    <span>%GC content:</span> {{gc_content}} <!-- %GC content will be displayed here -->
                </div>
            </section>

            <section id="sequence-lengths">
                <h2>Sequence lengths</h2>

                <!-- This will be populated with png plot --->
                <img src={{sequence_length_plot}} alt="Sequence Lengths Plot" style="max-width: 100%; height: auto;">
            </section>

            <section id="sequence-duplication-levels">
                <h2>Sequence duplication levels</h2>
                <table>
                    <thead>
                        <tr>
                            <th class="count_column">Count</th>
                            <th class="sequence_column">Sequence</th>
                        </tr>
                    </thead>
                    <tbody>
                        <!-- Table rows will be dynamically populated -->
                    </tbody>
                </table>
                <div id="sequence-duplication-levels-info">
                    <p>And {{sequence_duplication_levels_rest}} more</p>
                </div>
            </section>

            <section id="per-sequence-nucleotide-content">
                <h2>Per Sequence Nucleotide Content</h2>
                <img src={{per-sequence-nucleotide-content}} alt="Per Sequence Nucleotide Content" style="max-width: 100%; height: auto;">
            </section>
            <section id="per-sequence-dinucleotide-content">
                <h2>Per Sequence Dinucleotide Content</h2>
                <img src={{per-sequence-dinucleotide-content}} alt="Per Sequence Dinucleotide Content" style="max-width: 100%; height: auto;">
            </section>
            <section id="per-position-nucleotide-content">
                <h2>Per Position Nucleotide Content</h2>
                <img src={{per-position-nucleotide-content}} alt="Per Position Nucleotide Content" style="max-width: 100%; height: auto;">
            </section>
            <section id="per-position-reversed-nucleotide-content">
                <h2>Per Position Reversed Nucleotide Content</h2>
                <img src={{per-position-reversed-nucleotide-content}} alt="Per Position Reversed Nucleotide Content" style="max-width: 100%; height: auto;">
            </section>
            <section id="per-sequence-gc-content">
                <h2>Per Sequence GC Content</h2>
                <img src={{per-sequence-gc-content}} alt="Per Sequence GC Content" style="max-width: 100%; height: auto;">
            </section>
        </div>
    </div>

    <script>
        var sequenceDuplicationLevels = {{sequence_duplication_levels}};

        // Populate table for sequence duplication levels
        var tableBody = document.querySelector("#sequence-duplication-levels tbody");
        for (var sequence in sequenceDuplicationLevels) {
            var row = document.createElement("tr");
            var countCell = document.createElement("td");
            var sequenceCell = document.createElement("td");

            countCell.textContent = sequenceDuplicationLevels[sequence];
            countCell.className = "count_column";
            sequenceCell.textContent = sequence;
            sequenceCell.className = "sequence_column";

            row.appendChild(countCell);
            row.appendChild(sequenceCell);
            tableBody.appendChild(row);
        }
    </script>

</body>
</html>
"""

def get_sequence_html_template(stats, plots_path, tool_description=None):
    """
    Returns the HTML template for the report.
    """
    html_template = HTML_TEMPLATE

    # insert shared CSS and header fragment
    html_template = put_data(html_template, "{{common_css}}", COMMON_CSS)
    html_template = put_data(html_template, "{{report_header}}", REPORT_HEADER_HTML)

    # populate header placeholders: tool description, generated timestamp and input paths
    # Provide sensible defaults when values are not supplied
    if tool_description is None:
        tool_description = \
            """
            Toolkit for automated quality control of genomic datasets used in machine learning. 
            """
    input_path = stats['Filepath']
    generated_on = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html_template = put_data(html_template, "{{tool_description}}", tool_description)
    html_template = put_data(html_template, "{{generated_on}}", generated_on)
    html_template = put_data(html_template, "{{input_paths}}", input_path)
    html_template = put_data(html_template, "{{version}}", importlib.metadata.version("genbenchQC"))

    html_template = put_file_details(html_template, stats['Filename'])
    html_template = put_data(html_template, "{{label}}", stats['Label'] if stats['Label'] else "N/A")
    html_template = put_data(html_template, "{{seq_column}}", stats['Sequence column'] if stats['Sequence column'] else "N/A")
    html_template = put_data(html_template, "{{number_of_sequences}}", str(stats['Number of sequences']))
    html_template = put_data(html_template, "{{dedup_sequences}}", str(stats['Number of sequences left after deduplication']))
    html_template = put_data(html_template, "{{min_length}}", str(int(stats['Sequence lengths']['Sequence lengths'].min())))
    html_template = put_data(html_template, "{{max_length}}", str(int(stats['Sequence lengths']['Sequence lengths'].max())))
    html_template = put_data(html_template, "{{mean_length}}", f"{stats['Sequence lengths']['Sequence lengths'].mean():.2f}")
    html_template = put_data(html_template, "{{number_of_bases}}", str(stats['Number of bases']))
    html_template = put_data(html_template, "{{unique_bases}}", ', '.join(x for x in stats['Unique bases']))
    html_template = put_data(html_template, "{{gc_content}}", f"{(stats['%GC content']*100):.2f}")  

    html_template = put_data(html_template, "{{sequence_length_plot}}", str(plots_path['Sequence lengths']))
    html_template = put_data(html_template, "{{per-sequence-nucleotide-content}}", str(plots_path['Per sequence nucleotide content']))
    html_template = put_data(html_template, "{{per-sequence-dinucleotide-content}}", str(plots_path['Per sequence dinucleotide content']))
    html_template = put_data(html_template, "{{per-position-nucleotide-content}}", str(plots_path['Per position nucleotide content']))
    html_template = put_data(html_template, "{{per-position-reversed-nucleotide-content}}", str(plots_path['Per position reversed nucleotide content']))
    html_template = put_data(html_template, "{{per-sequence-gc-content}}", str(plots_path['Per sequence GC content']))

    # take max 10 sequences for sequence duplication levels - stats['Sequence duplication levels'] is a dictionary
    # with sequence as key and count as value, we convert it to a list of tuples
    sequence_duplication_levels = list(stats['Sequence duplication levels'].items())
    # Sort by count in descending order and take the top 10
    sequence_duplication_levels.sort(key=lambda x: x[1], reverse=True)
    # Limit to the first 10 sequences if there are more than 10
    if len(sequence_duplication_levels) > 10:
        # Take the top 10 sequences
        sequence_duplication_levels = sequence_duplication_levels[:10]
        html_template = put_data(html_template, "{{sequence_duplication_levels_rest}}", str(len(stats['Sequence duplication levels']) - 10))
    else:
        # If there are 10 or fewer sequences, set the rest to 0
        html_template = put_data(html_template, "{{sequence_duplication_levels_rest}}", "0")
    # Convert sequence duplication levels back to a dictionary-like structure for JSON-like string with sequence and count
    sequence_duplication_levels_dict = {str(seq): count for seq, count in sequence_duplication_levels}
    # Convert to JSON-like string for JavaScript
    sequence_duplication_levels_str = str(sequence_duplication_levels_dict).replace("'", '"').replace(", ", ",\n")
    # Replace the placeholder with the JSON-like string
    html_template = put_data(html_template, "{{sequence_duplication_levels}}", sequence_duplication_levels_str)
    
    return html_template