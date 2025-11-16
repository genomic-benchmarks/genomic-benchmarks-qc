from datetime import datetime

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HTML Report Output</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 0;
            padding: 0;
            display: flex;
            max-width: 100%; /* Prevents content from exceeding the viewport width */
        }
        .sidebar {
            width: 250px;
            background: #f4f4f4;
            padding: 20px;
            box-shadow: 2px 0 5px rgba(0, 0, 0, 0.1);
            height: 100vh; /* Full height */
            position: fixed; /* Fixed position */
            overflow-y: auto;
            z-index: 1000; /* Ensures it stays above other elements */
        }
        .sidebar a {
            display: block;
            margin: 10px 0;
            text-decoration: none;
            color: #333;
        }
        .content {
            margin-left: 300px;
            padding: 20px;
            width: calc(100% - 350px); /* Adjust width to avoid overlap */
            overflow-x: hidden;
        }
        #sequence-duplication-levels table {
            table-layout: fixed; /* Ensures columns respect defined widths */
            width: 100%; /* Makes the table take up the full width of its container */
            border-collapse: collapse; /* Optional: Makes the table look cleaner */
        }

        #sequence-duplication-levels table th,
        #sequence-duplication-levels table td {
            padding: 8px; /* Adds padding for better readability */
            border: 1px solid #ddd; /* Adds a border for clarity */
        }

        #sequence-duplication-levels table td {
            white-space: nowrap; /* Prevents text from wrapping */
            overflow: hidden; /* Hides overflowing text */
            text-overflow: ellipsis; /* Adds "..." to indicate clipped text */
        }

        #sequence-duplication-levels table th {
            text-align: left; /* Aligns header text to the left */
        }

        .count_column {
            width: 45px; /* Fixed width for the first column */
        }

        .sequence_column {
            width: 95%; /* Fixed width for the second column */
        }

        h1 {
            text-align: center;
            margin-bottom: 50px;
        }

        section {
            margin-bottom: 50px;
        }

        h2 {
            color: #333;
            border-bottom: 2px solid #ddd;
            padding-bottom: 5px;
        }

        h3 {
            color: #555;
            margin-bottom: 15px;
        }

        .chart-container {
            margin-bottom: 30px;
        }

        .data-item {
            font-size: 1.2em;
            margin-bottom: 10px;
        }

        .data-item span {
            font-weight: bold;
            font-size: 1em;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="sidebar">
            <h2>Navigation</h2>
            <a href="#basic-descriptive-statistics">Basic Descriptive Statistics</a>
            <div style="margin-left: 15px;">
                <a href="#filename">Filename</a>
                <a href="#label">Label</a>
                <a href="#seq_column">Sequence column</a>
                <a href="#num-sequences">Number of sequences</a>
                <a href="#dedup-sequences">Unique sequences</a>
                <a href="#min-length">Minimum length</a>
                <a href="#mean-length">Mean length</a>
                <a href="#max-length">Maximum length</a>
                <a href="#num-bases">Number of bases</a>
                <a href="#unique-bases">Unique bases</a>
                <a href="#gc-content">%GC content</a>
            </div>
            <a href="#general-descriptive-statistics">General Descriptive Statistics</a>
            <div style="margin-left: 15px;">
                <a href="#sequence-lengths">Sequence lengths</a>
                <a href="#sequence-duplication-levels">Sequence duplication levels</a>
            </div>
            <a href="#per-sequence-descriptive-stats">Per Sequence Descriptive Stats</a>
            <div style="margin-left: 15px;">
                <a href="#per-sequence-nucleotide-content">Per Sequence Nucleotide Content</a>
                <a href="#per-sequence-dinucleotide-content">Per Sequence Dinucleotide Content</a>
                <a href="#per-position-nucleotide-content">Per Position Nucleotide Content</a>
                <a href="#per-position-reversed-nucleotide-content">Per Position Reversed Nucleotide Content</a>
                <a href="#per-sequence-gc-content">Per Sequence GC Content</a>
            </div>
        </div>

        <div class="content">
        <h1>HTML Report Output</h1>

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

            <section id="general-descriptive-statistics">
                <h2>General Descriptive Statistics</h2>

                <h3 id="sequence-lengths">Sequence lengths</h3>
                <!-- This will be populated with png plot --->
                <img src={{sequence_length_plot}} alt="Sequence Lengths Plot" style="max-width: 100%; height: auto;">

                <div id="sequence-duplication-levels">
                    <h3>Sequence duplication levels</h3>
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
                </div>
            </section>

            <section id="per-sequence-descriptive-stats">
                <h2>Per Sequence Descriptive Stats</h2>

                <h3 id="per-sequence-nucleotide-content">Per Sequence Nucleotide Content</h3>
                <img src={{per-sequence-nucleotide-content}} alt="Per Sequence Nucleotide Content" style="max-width: 100%; height: auto;">

                <h3 id="per-sequence-dinucleotide-content">Per Sequence Dinucleotide Content</h3>
                <img src={{per-sequence-dinucleotide-content}} alt="Per Sequence Dinucleotide Content" style="max-width: 100%; height: auto;">

                <h3 id="per-position-nucleotide-content">Per Position Nucleotide Content</h3>
                <img src={{per-position-nucleotide-content}} alt="Per Position Nucleotide Content" style="max-width: 100%; height: auto;">

                <h3 id="per-position-reversed-nucleotide-content">Per Position Reversed Nucleotide Content</h3>
                <img src={{per-position-reversed-nucleotide-content}} alt="Per Position Reversed Nucleotide Content" style="max-width: 100%; height: auto;">

                <h3 id="per-sequence-gc-content">Per Sequence GC Content</h3>
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

def put_file_details(html_template, filename):
    """
    Populates the placeholders {{filename}} and {{date}} in the HTML template.

    Args:
        html_template (str): The HTML template as a string.
        filename (str): The name of the file to insert into the template.

    Returns:
        str: The updated HTML template with placeholders replaced.
    """
    # Replace {{filename}} with the stripped filename
    html_template = html_template.replace("{{filename}}", filename)

    # Replace {{date}} with the current date and time
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html_template = html_template.replace("{{date}}", current_time)

    return html_template

def put_data(html_template, placeholder, data):
    """
    Replaces all occurrences of a placeholder in the HTML template with the provided data.

    Args:
        html_template (str): The HTML template as a string.
        placeholder (str): The placeholder to replace (e.g., "{{placeholder}}").
        data (str): The data to replace the placeholder with.

    Returns:
        str: The updated HTML template with placeholders replaced.

    Raises:
        ValueError: If the placeholder is not found in the HTML template.
    """
    if placeholder not in html_template:
        raise ValueError(f"Placeholder not found: {placeholder}")

    # Replace all occurrences of the placeholder with the data
    return html_template.replace(placeholder, data)

def escape_str(s):
    """
    Add \" around the string to escape it for HTML.
    """
    return '"' + s + '"'

def get_sequence_html_template(stats, plots_path):
    """
    Populate the HTML report template with provided statistics and plot paths.
    
    Parameters:
        stats (dict): Mapping with keys required to fill the template, including:
            - 'Filename' (str): file name to display
            - 'Label' (str or None)
            - 'Sequence column' (str or None)
            - 'Number of sequences' (int)
            - 'Number of sequences left after deduplication' (int)
            - 'Sequence lengths' (object with a `Sequence lengths` numeric array supporting min(), max(), mean())
            - 'Number of bases' (int)
            - 'Unique bases' (iterable of str)
            - '%GC content' (numeric fraction, e.g., 0.42)
            - 'Sequence duplication levels' (dict mapping sequence -> count)
        plots_path (dict): Mapping of plot placeholder names to their rendered paths or HTML snippets. Expected keys include:
            'Sequence lengths', 'Per sequence nucleotide content', 'Per sequence dinucleotide content',
            'Per position nucleotide content', 'Per position reversed nucleotide content', 'Per sequence GC content'.
    
    Returns:
        str: HTML text with template placeholders replaced by the corresponding values from `stats` and `plots_path`.
    """
    html_template = HTML_TEMPLATE

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