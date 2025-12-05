from datetime import datetime


def put_file_details(html_template, filename):
    """
    Populates the placeholders {{filename}} and {{date}} in the HTML template.
    """
    # Replace {{filename}} with the stripped filename
    html_template = html_template.replace("{{filename}}", str(filename))

    # Replace {{date}} with the current date and time
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html_template = html_template.replace("{{date}}", current_time)

    return html_template


def put_data(html_template, placeholder, data):
    """
    Replaces all occurrences of a placeholder in the HTML template with the provided data.
    """
    if placeholder not in html_template:
        raise ValueError(f"Placeholder not found: {placeholder}")
    return html_template.replace(placeholder, str(data))


def escape_str(s):
    """Return a JSON/JS quoted string for safe embedding in templates."""
    return '"' + str(s).replace('"', '\\"') + '"'


def icon_html(summary_statuses, key):
    """
    Return small circular icon HTML for a summary status mapping.

    summary_statuses: mapping or None
    key: lookup key
    """
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


# Shared CSS used by all report HTML templates. Keep this as plain CSS (no <style> tags).
COMMON_CSS = """
/* Make sure the root uses border-box sizing and doesn't allow the
   page to horizontally overflow. Individual <pre> blocks will
   still allow internal horizontal scrolling. */
html, body {
    box-sizing: border-box;
    overflow-x: hidden; /* prevent accidental page-wide horizontal scroll */
    margin: 0;
    padding: 0;
}
body {
    font-family: Arial, sans-serif;
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
.sidebar a {
    display: inline-block;
    margin-left: 10px;
    text-decoration: none;
    color: #333;
    width: 175px;
}
.sidebar-item {
    display: flex;
    align-items: center;
    margin-bottom: 10px;
}
.status-icon {
    display: inline-flex;
    width: 40px;
    height: 40px;
    line-height: 40px;
    border-radius: 50%;
    color: white;
    font-size: 25px;
    font-family: verdana;
    text-align: center;
    justify-content: center;
}
.status-pass { background-color: #2e7d32; color: #ffffff; }
.status-warn { background-color: #f57f17; color: #ffffff; }
.status-fail { background-color: #c62828; color: #ffffff; }
.container {
    /* in some layouts the flex parent can grow to fit unbreakable children;
       make sure the container can shrink so .content's min-width:0 works */
    min-width: 0;
    box-sizing: border-box;
}
.content {
    margin-left: 260px;
    padding: 10px;
    width: calc(100% - 260px); /* Adjust width to avoid overlap */
    overflow-x: hidden;
    /* allow flex child to shrink below its content width so long lines
       inside (e.g. long sequences in <pre>) don't expand the whole page */
    min-width: 0;
    box-sizing: border-box;
}
#sequence-duplication-levels table {
    table-layout: fixed; /* Ensures columns respect defined widths */
    width: 100%;
    border-collapse: collapse;
}
#sequence-duplication-levels table th,
#sequence-duplication-levels table td {
    padding: 8px;
    border: 1px solid #ddd;
}
#sequence-duplication-levels table td {
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
}
.count_column { width: 45px; }
.sequence_column { width: 95%; }
h1 { text-align: center; margin-bottom: 50px; }
section { margin-bottom: 50px; }
h2 { color: #333; border-bottom: 2px solid #ddd; padding-bottom: 5px; margin-left: 10px; }
h3 { color: #555; margin-bottom: 15px; }
.chart-container { margin-bottom: 30px; }
.data-item { font-size: 1.2em; margin-bottom: 10px; }
.data-item span { font-weight: bold; font-size: 1em; }

/* basic descriptive statistics table */
#basic-descriptive-statistics table { width: 100%; border-collapse: collapse; margin: 20px 0; }
#basic-descriptive-statistics table td { padding: 10px; border: 1px solid #ddd; text-align: center; }
#basic-descriptive-statistics table td:first-child { text-align: left; width: 200px; }
#basic-descriptive-statistics table tr:nth-child(even) { background-color: #f9f9f9; }
#basic-descriptive-statistics table tr:hover { background-color: #f1f1f1; }
#basic-descriptive-statistics table span { font-weight: bold; }
.cluster { 
    border: 1px solid #ccc; 
    margin-bottom: 20px; 
    padding: 15px; 
    border-radius: 5px; 
    background: #fff; 
}
pre {
     background: #f9f9f9;
     padding: 10px;
     /* prevent <pre> from forcing page width; allow internal horizontal
         scrolling when the sequence is longer than the content area */
     max-width: 100%;
     overflow-x: auto;
     /* keep original whitespace formatting but don't let the pre tag
         itself grow the parent container in flex layout */
     white-space: pre;
     box-sizing: border-box;
}
.report-header img { 
    max-width: 350px; 
    height: auto;
    margin: 0 0 10px; 
}
"""


# Standard report header HTML fragment (uses placeholders)
REPORT_HEADER_HTML = """
<div class="report-header" style="text-align: left; margin-bottom: 30px;">
    <img src="https://raw.githubusercontent.com/katarinagresova/GenBenchQC/main/assets/logo_with_text_transparent.png" alt="GenBenchQC Logo">
    <div style="font-size: 1.05em; color: #333; margin-bottom: 6px;">{{tool_description}}</div>
    <div style="font-size: 0.9em; color: #666;">Report generated on {{generated_on}} based on data: {{input_paths}}</div>
</div>
"""
