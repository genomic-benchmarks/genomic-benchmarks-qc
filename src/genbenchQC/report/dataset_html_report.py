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
            width: 220px;
            background: #f4f4f4;
            padding: 15px;
            box-shadow: 2px 0 5px rgba(0, 0, 0, 0.1);
            height: 100vh; /* Full height */
            position: fixed; /* Fixed position */
            overflow-y: auto;
            z-index: 1000; /* Ensures it stays above other elements */
        }

        /* Each sidebar entry (icon + link) is a horizontal row */
        .sidebar-item {
            display: flex;
            align-items: center;
            margin-bottom: 10px;
        }

        /* Icon that appears before each sidebar link. Rendered as a colored circle. */
        .status-icon {
            display: inline-flex;
            width: 40px;
            height: 40px;
            line-height: 40px;
            border-radius: 50%;
            color: white;
            font-size: 25px;
            font-family:verdana;
            text-align: center;
            justify-content: center;
        }

        /* Visual classes for icon states */
        .status-pass {
            background-color: #2e7d32; /* green */
            color: #ffffff; /* white check */
        }

        .status-warn {
            background-color: #f57f17; /* orange */
            color: #ffffff; /* white */
        }

        .status-fail {
            background-color: #c62828; /* red */
            color: #ffffff; /* white */
        }
        /* keep the <a> inline inside the flex row so the icon sits beside it */
        .sidebar a {
            display: inline-block;
            text-decoration: none;
            color: #333;
            margin-left: 10px;
            width: 175px;

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
            margin-left: 10px;
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

        /* Style basic descriptive statistics table */
        #basic-descriptive-statistics table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }
        #basic-descriptive-statistics table td {
            padding: 10px;
            border: 1px solid #ddd;
            text-align: center;
        }
        #basic-descriptive-statistics table td:first-child {
            text-align: left;
            width: 200px;
        }
        #basic-descriptive-statistics table tr:nth-child(even) {
            background-color: #f9f9f9;
        }
        #basic-descriptive-statistics table tr:hover {
            background-color: #f1f1f1;
        }

        #basic-descriptive-statistics table span {
            font-weight: bold;
        }

    </style>
</head>
<body>
    <div class="container">
            <div class="sidebar">
            <h2>Summary</h2>
            <div class="sidebar-item">{{icon_basic_descriptive_statistics}}<a href="#basic-descriptive-statistics">Basic Descriptive Statistics</a></div>
            <div class="sidebar-item">{{icon_sequence_lengths}}<a href="#sequence-lengths">Sequence lengths</a></div>
            <div class="sidebar-item">{{icon_sequence_duplication_levels}}<a href="#sequence-duplication-levels">Duplicate sequences</a></div>
            <div class="sidebar-item">{{icon_per_sequence_nucleotide_content}}<a href="#per-sequence-nucleotide-content">Per Sequence Nucleotide Content</a></div>
            <div class="sidebar-item">{{icon_per_sequence_dinucleotide_content}}<a href="#per-sequence-dinucleotide-content">Per Sequence Dinucleotide Content</a></div>
            <div class="sidebar-item">{{icon_per_position_nucleotide_content}}<a href="#per-position-nucleotide-content">Per Position Nucleotide Content</a></div>
            <div class="sidebar-item">{{icon_per_position_reversed_nucleotide_content}}<a href="#per-position-reversed-nucleotide-content">Per Position Reversed Nucleotide Content</a></div>
            <div class="sidebar-item">{{icon_per_sequence_gc_content}}<a href="#per-sequence-gc-content">Per Sequence GC Content</a></div>
        </div>

        <div class="content">
        <h1>HTML Report Output</h1>

            <section id="basic-descriptive-statistics">
                <div class="sidebar-item">
                    {{icon_basic_descriptive_statistics}}
                    <h2>Basic Descriptive Statistics</h2>
                </div>
                <table style="width: 100%; border-collapse: collapse; margin: 20px 0;">
                    <tr id="filename">
                        <td><span>Filename</span></td>
                        <td style="text-align: center;">{{filename1}}</td>
                        <td style="text-align: center;">{{filename2}}</td>
                    </tr>
                    <tr id="label">
                        <td><span>Label</span></td>
                        <td style="text-align: center;">{{label1}}</td>
                        <td style="text-align: center;">{{label2}}</td>
                    </tr>
                    <tr id="seq_column">
                        <td><span>Sequence column</span></td>
                        <td style="text-align: center;">{{seq_col1}}</td>
                        <td style="text-align: center;">{{seq_col2}}</td>
                    </tr>
                    <tr id="num-sequences">
                        <td><span>Number of sequences</span></td>
                        <td style="text-align: center;">{{number_of_sequences1}}</td>
                        <td style="text-align: center;">{{number_of_sequences2}}</td>
                    </tr>
                    <tr id="dedup-sequences">
                        <td><span>Unique sequences</span></td>
                        <td style="text-align: center;">{{dedup_sequences1}}</td>
                        <td style="text-align: center;">{{dedup_sequences2}}</td>
                    </tr>
                    <tr id="min-length">
                        <td><span>Minimum length</span></td>
                        <td style="text-align: center;">{{min_length1}}</td>
                        <td style="text-align: center;">{{min_length2}}</td>
                    </tr>
                    <tr id="mean-length">
                        <td><span>Mean length</span></td>
                        <td style="text-align: center;">{{mean_length1}}</td>
                        <td style="text-align: center;">{{mean_length2}}</td>
                    </tr>
                    <tr id="max-length">
                        <td><span>Maximum length</span></td>
                        <td style="text-align: center;">{{max_length1}}</td>
                        <td style="text-align: center;">{{max_length2}}</td>
                    </tr>
                    <tr id="num-bases">
                        <td><span>Number of bases</span></td>
                        <td style="text-align: center;">{{number_of_bases1}}</td>
                        <td style="text-align: center;">{{number_of_bases2}}</td>
                    </tr>
                    <tr id="unique-bases">
                        <td><span>Unique bases</span></td>
                        <td style="text-align: center;">{{unique_bases1}}</td>
                        <td style="text-align: center;">{{unique_bases2}}</td>
                    </tr>
                    <tr id="gc-content">
                        <td><span>%GC content</span></td>
                        <td style="text-align: center;">{{gc_content1}}</td>
                        <td style="text-align: center;">{{gc_content2}}</td>
                    </tr>
                </table>
            </section>

            <section id="sequence-lengths">
                <div class="sidebar-item">
                    {{icon_sequence_lengths}}
                    <h2>Sequence Lengths</h2>
                </div>

                <!-- This will be populated with png plot --->
                <img src={{sequence_length_plot}} alt="Sequence Lengths Plot" style="max-width: 50%; height: auto; display: block; margin: 0 auto;">
            </section>
            <section id="sequence-duplication-levels">
                <div class="sidebar-item">
                    {{icon_sequence_duplication_levels}}
                    <h2>Duplicate sequences</h2>
                </div>
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
                    <p>And {{sequence_duplication_levels_rest}} more</p>
                </div>
            </section>

            <section id="per-sequence-nucleotide-content">
                <div class="sidebar-item">
                    {{icon_per_sequence_nucleotide_content}}
                    <h2>Per Sequence Nucleotide Content</h2>
                </div>
                <img src={{per-sequence-nucleotide-content}} alt="Per Sequence Nucleotide Content" style="max-width: 100%; height: auto;">
            </section>
            <section id="per-sequence-dinucleotide-content">
                <div class="sidebar-item">
                    {{icon_per_sequence_dinucleotide_content}}
                    <h2>Per Sequence Dinucleotide Content</h2>
                </div>
                <img src={{per-sequence-dinucleotide-content}} alt="Per Sequence Dinucleotide Content" style="max-width: 100%; height: auto;">
            </section>
            <section id="per-position-nucleotide-content">
                <div class="sidebar-item">
                    {{icon_per_position_nucleotide_content}}
                    <h2>Per Position Nucleotide Content</h2>
                </div>
                <img src={{per-position-nucleotide-content}} alt="Per Position Nucleotide Content" style="max-width: 100%; height: auto;">
            </section>
            <section id="per-position-reversed-nucleotide-content">
                <div class="sidebar-item">
                    {{icon_per_position_reversed_nucleotide_content}}
                    <h2>Per Position Reversed Nucleotide Content</h2>
                </div>
                <img src={{per-position-reversed-nucleotide-content}} alt="Per Position Reversed Nucleotide Content" style="max-width: 100%; height: auto;">
            </section>
            <section id="per-sequence-gc-content">
                <div class="sidebar-item">
                    {{icon_per_sequence_gc_content}}
                    <h2>Per Sequence GC Content</h2>
                </div>
                <img src={{per-sequence-gc-content}} alt="Per Sequence GC Content" style="max-width: 50%; height: auto; display: block; margin: 0 auto;">
            </section>
        </div>
    </div>

    <script>
        var sequenceDuplicationLevels = {{sequence_duplication_levels}};

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

def get_dataset_html_template(stats1, stats2, plots_path, duplicate_seqs, summary_statuses=None):
    """
    Returns the HTML template for the report.

    Args:
        stats1, stats2: objects containing dataset statistics (unchanged API).
        plots_path: dict with plot image paths.
        duplicate_seqs: list of duplicate sequences.
        summary_statuses: optional pd.DataFrame mapping status keys to HTML/text snippets for
                          the sidebar placeholders.

    Backwards compatible: if summary_statuses is None, placeholders are left empty.
    """
    html_template = HTML_TEMPLATE

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

    html_template = put_data(html_template, "{{sequence_length_plot}}", str(plots_path['Sequence lengths']))
    html_template = put_data(html_template, "{{per-sequence-nucleotide-content}}", str(plots_path['Per sequence nucleotide content']))
    html_template = put_data(html_template, "{{per-sequence-dinucleotide-content}}", str(plots_path['Per sequence dinucleotide content']))
    html_template = put_data(html_template, "{{per-position-nucleotide-content}}", str(plots_path['Per position nucleotide content']))
    html_template = put_data(html_template, "{{per-position-reversed-nucleotide-content}}", str(plots_path['Per position reversed nucleotide content']))
    html_template = put_data(html_template, "{{per-sequence-gc-content}}", str(plots_path['Per sequence GC content']))

    # Populate sidebar icon placeholders (if provided). summary_statuses may contain
    # simple status keywords ('pass', 'warn', 'fail') or an HTML snippet. This helper
    # returns the HTML for the small circular icon shown before each section link.
    def _icon_html(key):
        if summary_statuses is None:
            return ''
        val = summary_statuses.get(key, '')
        if not val:
            return ''
        s = str(val).strip()
        lv = s.lower()
        if lv in ('pass', 'ok', 'good', 'success'):
            return '<span class="status-icon status-pass">✔</span>'
        if lv in ('warn', 'warning'):
            return '<span class="status-icon status-warn">!</span>'
        if lv in ('fail', 'failed', 'error'):
            return '<span class="status-icon status-fail">✖</span>'
        # Otherwise assume the value is an HTML snippet or a custom symbol and return as-is
        return s

    html_template = put_data(html_template, "{{icon_basic_descriptive_statistics}}", _icon_html('Unique bases'))
    html_template = put_data(html_template, "{{icon_sequence_lengths}}", _icon_html('Sequence lengths'))
    html_template = put_data(html_template, "{{icon_sequence_duplication_levels}}", _icon_html('Duplication between labels'))
    html_template = put_data(html_template, "{{icon_per_sequence_nucleotide_content}}", _icon_html('Per sequence nucleotide content'))
    html_template = put_data(html_template, "{{icon_per_sequence_dinucleotide_content}}", _icon_html('Per sequence dinucleotide content'))
    html_template = put_data(html_template, "{{icon_per_position_nucleotide_content}}", _icon_html('Per position nucleotide content'))
    html_template = put_data(html_template, "{{icon_per_position_reversed_nucleotide_content}}", _icon_html('Per reverse position nucleotide content'))
    html_template = put_data(html_template, "{{icon_per_sequence_gc_content}}", _icon_html('Per sequence GC content'))

    # take max 10 sequences for sequence duplication levels
    sequence_duplication_levels = duplicate_seqs[:10]
    html_template = put_data(html_template, "{{sequence_duplication_levels}}",
                             '[' + ', '.join(escape_str(seq) for seq in sequence_duplication_levels) + ']')
    if len(duplicate_seqs) > 10:
        # If there are more than 10 sequences, we show how many more there are
        # and set the rest to a placeholder
        html_template = put_data(html_template, "{{sequence_duplication_levels_rest}}", str(len(duplicate_seqs) - 10))
    else:
        # If there are 10 or fewer sequences, we set the rest to 0
        html_template = put_data(html_template, "{{sequence_duplication_levels_rest}}", "0")

    return html_template