"""Tests for what ends up in the generated HTML.

The per-position panels are not embedded as a picture: the page carries the
numbers and draws them in a canvas, so a reader can zoom into a flagged position
instead of squinting at a 1px band. Two properties have to hold for that to be
worth anything - the data has to be in the page, and the page has to stay a
single standalone file that opens with no network and no sibling assets.
"""

import json
import logging
import re
from urllib.parse import urlparse

import numpy as np
import pandas as pd
import pytest
from helpers import mmseqs_hit, write_csv, write_mmseqs_output

from genomic_benchmarks_qc import evaluate_splits
from genomic_benchmarks_qc.report.alignment_rendering import has_reversed_coordinates
from genomic_benchmarks_qc.report.report_generator import generate_dataset_html_report
from genomic_benchmarks_qc.report.split_html_report import (
    alignment_error_html,
    alignments_count_text,
)
from genomic_benchmarks_qc.report.utils import escape_str
from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics
from genomic_benchmarks_qc.utils.testing import flag_significant_differences

PAYLOAD_IDS = ('ppv-fwd-data', 'ppv-rev-data')


def make_stats(sequences, label='cls'):
    stats = SequenceStatistics(sequences=sequences, filename=f'{label}.fa',
                               filepath=f'/{label}.fa', label=label)
    stats.compute()
    return stats


def random_sequences(count, length, seed, composition=None):
    rng = np.random.default_rng(seed)
    return [''.join(rng.choice(list('ACGT'), size=length, p=composition)) for _ in range(count)]


def render(tmp_path, stats1, stats2):
    """Run the class report end to end and return the page."""
    results, failed_by_feature = flag_significant_differences(stats1, stats2)
    generate_dataset_html_report(
        stats1, stats2, tmp_path / 'gb-qc-report.html', plots_path=tmp_path / 'plots',
        plot_type='boxen',
        results=pd.DataFrame.from_dict(results, orient='index'),
        failed_by_feature=failed_by_feature,
    )
    return (tmp_path / 'gb-qc-report.html').read_text()


def payload_from(page, dom_id):
    """Parse one embedded payload back out of the page."""
    raw = re.search(rf'id="{dom_id}">(.*?)</script>', page, re.S).group(1)
    unescaped = (raw.replace('\\u003c', '<').replace('\\u003e', '>').replace('\\u0026', '&'))
    return json.loads(unescaped)


@pytest.fixture(scope='module')
def page(tmp_path_factory):
    """A report with enough sequences for the per-position checks to be scored."""
    tmp_path = tmp_path_factory.mktemp('report')
    stats1 = make_stats(random_sequences(300, 30, seed=1), label='a')
    stats2 = make_stats(random_sequences(300, 30, seed=2, composition=[.7, .1, .1, .1]), label='b')
    return render(tmp_path, stats1, stats2)


@pytest.fixture(scope='module')
def disjoint_page(tmp_path_factory):
    """Two classes with no base in common: the comparison cannot be made."""
    tmp_path = tmp_path_factory.mktemp('disjoint')
    stats1 = make_stats(['AAAAAAAAAA'] * 10, label='a')
    stats2 = make_stats(['CCCCCCCCCC'] * 10, label='c')
    return render(tmp_path, stats1, stats2)


class TestInteractivePerPosition:
    def test_both_directions_arrive_as_data(self, page):
        for dom_id in PAYLOAD_IDS:
            payload = payload_from(page, dom_id)
            assert payload['endPosition'] > 0
            assert payload['labels'] == ['a', 'b']
            assert all(len(series) == payload['endPosition'] for series in payload['coverage'])
            assert set(payload['nucleotides']) == set('ACGT')

    def test_the_flagged_positions_reach_the_page(self, page):
        payload = payload_from(page, 'ppv-fwd-data')

        flagged = {position for base in payload['flags'] for position in payload['flags'][base]}
        assert flagged, 'this comparison differs at every position'
        assert all('Unknown' not in payload['flags'][base].values() for base in payload['flags'])

    def test_the_viewer_and_its_data_are_both_in_the_page(self, page):
        assert 'id="ppv-fwd"' in page and 'id="ppv-rev"' in page
        assert 'initPerPositionViewers' in page
        # the viewer itself, not just the call
        assert 'global.initPerPositionViewers = function' in page

    def test_the_per_position_panels_are_no_longer_a_picture(self, page):
        images = re.findall(r'<img[^>]*alt="([^"]*)"', page)
        assert images, 'the other plots are still images'
        assert not [alt for alt in images if 'Position' in alt]

    def test_the_static_plots_are_still_written_to_disk(self, tmp_path):
        """The report does not show them, but the PNGs remain part of the output."""
        stats1 = make_stats(random_sequences(300, 30, seed=3), label='a')
        stats2 = make_stats(random_sequences(300, 30, seed=4), label='b')

        render(tmp_path, stats1, stats2)

        assert (tmp_path / 'plots' / 'per_position_nucleotide_content.png').is_file()
        assert (tmp_path / 'plots' / 'per_position_reversed_nucleotide_content.png').is_file()

    def test_the_figures_that_must_line_up_share_a_class(self, page):
        """Sizing lives in the stylesheet: the per-position canvas has to match
        the static figures above it, which it cannot do from inline widths."""
        assert 'class="plot-wide"' in page
        assert 'max-width: 108%' not in page

    def test_the_coverage_note_sits_inside_the_sections_explanation(self, page):
        """It is background for reading the plot, so it belongs behind the ?
        button with the rest of the explanation - but inside a block that button
        actually toggles, or it could never be read at all."""
        explanation = re.search(
            r'<div id="per-position-explanation" class="explanation-text">(.*?)</div>',
            page, re.S)

        assert explanation, 'the toggled explanation block is still there'
        assert 'Positions 1' in explanation.group(1)
        assert "toggleExplanation('per-position-explanation')" in page


class TestTemplates:
    """The page templates are .html files in report.assets, so they can be
    edited with an editor's help. put_data raises when the code fills a
    placeholder the template does not have; nothing was catching the other
    direction - a slot added to the markup that no code ever fills."""

    def test_no_placeholder_survives_into_a_class_report(self, page):
        assert re.findall(r'\{\{\w+\}\}', page) == []

    def test_no_placeholder_survives_into_a_split_report(self, tmp_path, monkeypatch):
        rendered = split_page(tmp_path, monkeypatch,
                              [mmseqs_hit('seq_0_test', 'seq_0_train')])

        assert re.findall(r'\{\{\w+\}\}', rendered) == []


class TestStandalone:
    def test_nothing_is_loaded_from_outside_the_file(self, page):
        """A report opens complete with no network. Subresources only - an <a>
        pointing at the documentation is somewhere the reader can go, not
        something the page fetches to render itself."""
        subresources = [
            ref for ref in re.findall(r'\bsrc="([^"]*)"', page)
            + re.findall(r'<link[^>]+href="([^"]*)"', page)
            if not ref.strip().startswith(('#', 'data:'))
        ]

        assert subresources == []

    def test_the_links_out_go_to_the_project(self, page):
        """The only outbound links are to the tool: nowhere else is worth
        sending a reader who was handed this file."""
        hosts = {urlparse(href).netloc
                 for href in re.findall(r'<a[^>]+href="(http[^"]*)"', page)}

        assert hosts == {'genomic-benchmarks.github.io', 'github.com'}

    def test_the_styling_and_behaviour_are_inlined(self, page):
        assert '<style>' in page
        assert 'window.toggleExplanation' in page          # report_ui.js
        assert '.ppv-canvas' in page                       # per_position_viewer.css
        assert '.plot-wide' in page                        # report_design.css

    def test_there_is_no_theme_switching(self, page):
        """Light only: the plots are matplotlib rasters that cannot follow a theme."""
        assert 'data-theme' not in page
        assert 'prefers-color-scheme' not in page


class TestNothingToDraw:
    def test_the_message_replaces_the_figure(self, disjoint_page):
        assert 'no-plot-message' in disjoint_page
        assert 'id="ppv-fwd"' not in disjoint_page

    def test_the_viewer_is_left_out_entirely(self, disjoint_page):
        """No point shipping 36 KB of drawing code for nothing to draw."""
        assert 'initPerPositionViewers' not in disjoint_page


class TestDuplicateSequences:
    def test_shared_sequences_arrive_as_data_not_as_code(self, tmp_path):
        shared = 'ACGTACGTAC' * 3
        stats1 = make_stats([shared] + random_sequences(20, 30, seed=5), label='a')
        stats2 = make_stats([shared] + random_sequences(20, 30, seed=6), label='b')

        page = render(tmp_path, stats1, stats2)

        block = re.search(r'id="duplicate-sequences">(.*?)</script>', page, re.S).group(1)
        assert json.loads(block) == [shared]


def split_page(tmp_path, monkeypatch, hits):
    """Run the split report end to end against a fixed set of MMseqs2 hits."""
    train = write_csv(tmp_path / 'train.csv', ['0'], rows_per_label=5)
    test = write_csv(tmp_path / 'test.csv', ['0'], rows_per_label=5)

    def fake_run_search(query_fasta, target_fasta, output_path, tmp_dir, **kwargs):
        return write_mmseqs_output(output_path, hits)

    monkeypatch.setattr(evaluate_splits.mmseqs_runtime, 'run_search', fake_run_search)
    evaluate_splits.run(train_files=[train], test_files=[test], format='csv',
                        out_folder=str(tmp_path / 'out'), report_types=['html'])

    return (tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test'
            / 'gb-qc-report.html').read_text()


class TestSplitReport:
    def test_it_gets_the_shared_behaviour_and_its_own(self, tmp_path, monkeypatch):
        page = split_page(tmp_path, monkeypatch,
                          [mmseqs_hit('seq_0_test', 'seq_0_train')])

        assert 'window.toggleExplanation' in page      # report_ui.js
        assert 'window.toggleAlignment' in page        # split_report.js
        assert 'initPerPositionViewers' not in page    # no per-position figure here
        assert '{{report_scripts}}' not in page

    def test_the_leakage_listing_uses_the_shared_components(self, tmp_path, monkeypatch):
        """The high-similarity pairs are the same kind of thing as the class
        report's flagged positions - a check's findings - and are shown in the
        same collapsible panel and listing table, styled from the shared sheet."""
        page = split_page(tmp_path, monkeypatch,
                          [mmseqs_hit('seq_0_test', 'seq_0_train')])

        assert 'class="qc-panel"' in page
        assert 'class="qc-listing"' in page
        assert '1 high-similarity alignment' in page
        # the panel sits inside the check it belongs to, not in a card of its own
        section = re.search(r'<section id="similarity-section".*?</section>', page, re.S)
        assert section and 'qc-panel' in section.group(0)
        assert '.qc-listing' in page                   # report_design.css
        assert '.alignment-block' in page              # split_report.css

    def test_the_count_says_when_the_listing_was_capped(self):
        """The page carries at most ROW_CAP rows, and the hits it drops are in
        the exported TSV - so a capped listing has to say so, or 100 rows read as
        'there were exactly 100'."""
        assert alignments_count_text(0, 0) == 'No high-similarity alignments'
        assert alignments_count_text(1, 1) == '1 high-similarity alignment'
        assert alignments_count_text(4, 4) == '4 high-similarity alignments'
        assert alignments_count_text(250, 100) == (
            '250 high-similarity alignments (first 100 shown)')

    def test_the_filenames_are_text_not_markup(self, tmp_path, monkeypatch):
        """File names are as much input as the sequences are, and land in the
        same page someone else opens."""
        train = write_csv(tmp_path / '<img src=x onerror=alert(1)>.csv', ['0'], rows_per_label=5)
        test = write_csv(tmp_path / 'test.csv', ['0'], rows_per_label=5)

        def fake_run_search(query_fasta, target_fasta, output_path, tmp_dir, **kwargs):
            return write_mmseqs_output(output_path, [mmseqs_hit('seq_0_test', 'seq_0_train')])

        monkeypatch.setattr(evaluate_splits.mmseqs_runtime, 'run_search', fake_run_search)
        evaluate_splits.run(train_files=[train], test_files=[test], format='csv',
                            out_folder=str(tmp_path / 'out'), report_types=['html'])
        page = next((tmp_path / 'out').rglob('gb-qc-report.html')).read_text()

        assert '<img src=x onerror=alert(1)>' not in page
        assert '&lt;img src=x onerror=alert(1)&gt;' in page

    def test_a_clean_split_says_so_in_the_panel(self, tmp_path, monkeypatch):
        """An empty listing has to read as 'nothing leaked', not as a table that
        failed to build - the same distinction the flag panel draws."""
        page = split_page(tmp_path, monkeypatch, [])

        assert 'No high-similarity alignments' in page
        assert 'qc-empty' in page
        assert '{{results_body}}' not in page


class TestUserDataIsEscaped:
    """Nothing read out of the input files may reach the page as markup.

    A report is a single file that gets mailed on and opened by someone who did
    not run the tool, so a label or a base is untrusted text by the time it is
    rendered - and a stray `<` swallowing the rest of a cell is the quiet half
    of the same bug.
    """

    def test_a_label_carrying_markup_stays_text(self, tmp_path):
        payload = '<img src=x onerror=alert(1)>'
        page = render(tmp_path,
                      make_stats(random_sequences(30, 20, seed=5), label=payload),
                      make_stats(random_sequences(30, 20, seed=6), label='b'))

        assert payload not in page
        assert '&lt;img src=x onerror=alert(1)&gt;' in page

    def test_a_label_cannot_become_a_placeholder(self, tmp_path):
        """The page is filled one placeholder at a time, so a value can carry one.

        `{{filename1}}` is filled before `{{label1}}` is, so a file whose name
        is literally `{{label1}}.csv` used to land in the page while `{{label1}}`
        was still to come, and the filename cell ended up showing the class
        name. `{` is written as `&#123;`, which renders the same and matches
        nothing.
        """
        stats1 = make_stats(random_sequences(30, 20, seed=5), label='a')
        stats1.filename = '{{label1}}.csv'
        page = render(tmp_path, stats1,
                      make_stats(random_sequences(30, 20, seed=6), label='b'))

        assert '&#123;&#123;label1}}.csv' in page
        assert '{{label1}}.csv' not in page
        assert 'a.csv' not in page

    def test_a_label_that_looks_like_a_placeholder_is_shown_as_one(self, tmp_path):
        """The label reaches every place a label reaches, still as its own text."""
        page = render(tmp_path,
                      make_stats(random_sequences(30, 20, seed=5), label='{{label2}}'),
                      make_stats(random_sequences(30, 20, seed=6), label='b'))

        assert '&#123;&#123;label2}}' in page
        assert re.findall(r'\{\{\w+\}\}', page) == []

    def test_a_label_in_the_per_position_data_stays_data(self, tmp_path):
        """The payload is JSON, so it escapes to \\u007b rather than to an entity."""
        page = render(tmp_path,
                      make_stats(random_sequences(300, 30, seed=1), label='{{report_scripts}}'),
                      make_stats(random_sequences(300, 30, seed=2), label='b'))

        assert payload_from(page, 'ppv-fwd-data')['labels'][0] == '{{report_scripts}}'
        assert '"{{report_scripts}}"' not in page

    def test_a_sequence_going_into_a_script_element_is_escaped_the_same_way(self):
        """The duplicate-sequences listing is JSON too, and lands the same way."""
        assert escape_str('{{label1}}') == '"\\u007b\\u007blabel1}}"'
        assert json.loads(escape_str('{{label1}}')) == '{{label1}}'

    def test_a_base_that_is_a_metacharacter_stays_text(self, tmp_path):
        """Unusual bases are exactly what the Unique bases check exists to
        report, so the check must survive being handed one."""
        page = render(tmp_path,
                      make_stats(['AC<GT' * 4] * 30, label='a'),
                      make_stats(['ACGT' * 5] * 30, label='b'))

        # the cell listing the bases survives intact, with the odd one escaped
        assert '<td style="text-align: center;">&lt;, A, C, G, T</td>' in page
        assert '<td style="text-align: center;"><, A, C, G, T</td>' not in page


class TestAHitThatCannotBeDrawn:
    """What the page says when an alignment will not render.

    Two causes, and only one of them is known. A hit whose coordinates run
    backwards came from the MMseqs2 build - the search asked for the forward
    strand only - and its scores are still right, so the page can say what
    happened and what avoids it. Anything else is the validator refusing a hit
    that disagrees with its sequences, and guessing at a cause there was how
    this message came to blame conda for everything.
    """

    def test_backwards_coordinates_are_recognised(self):
        forward = mmseqs_hit('seq_0_test', 'seq_0_train')
        assert not has_reversed_coordinates(forward)
        assert has_reversed_coordinates({**forward, 'tstart': 60, 'tend': 1})
        assert has_reversed_coordinates({**forward, 'qstart': 60, 'qend': 1})

    def test_a_backwards_row_is_told_what_happened_and_what_to_do(self):
        cell = alignment_error_html(ValueError('Target end (1) < start (60)'), True)

        assert 'reverse strand' in cell
        assert 'scores in this row are unaffected' in cell
        assert '--threads 1' in cell
        assert 'precompiled' in cell

    def test_any_other_failure_does_not_blame_the_installation(self):
        cell = alignment_error_html(ValueError('Query alignment mismatch'), False)

        assert 'disagree' in cell
        assert 'threads' not in cell
        assert 'conda' not in cell

    def test_the_error_text_is_escaped_into_the_cell(self):
        cell = alignment_error_html(ValueError('<img src=x onerror=alert(1)>'), False)

        assert '<img src=x onerror=alert(1)>' not in cell
        assert '&lt;img src=x onerror=alert(1)&gt;' in cell

    def test_the_page_carries_the_explanation_and_keeps_the_scores(
            self, tmp_path, monkeypatch):
        """End to end: the pair stays in the listing with its numbers, because
        what a backwards row loses is the picture, not the measurement."""
        backwards = {**mmseqs_hit('seq_0_test', 'seq_0_train'), 'tstart': 60, 'tend': 1}
        page = split_page(tmp_path, monkeypatch, [backwards])

        assert 'ALIGNMENT VISUALISATION ERROR' in page
        assert 'reverse strand' in page
        assert '--threads 1' in page
        assert '1 high-similarity alignment' in page
        assert '96.9' in page                       # the pident is still on the page

    def test_the_old_guess_is_gone(self, tmp_path, monkeypatch):
        """The message used to say every failure was 'often caused by a
        conda-installed mmseqs2' and to recommend precompiled binaries over
        compiling from source - a contrast that was never the one that mattered."""
        backwards = {**mmseqs_hit('seq_0_test', 'seq_0_train'), 'tstart': 60, 'tend': 1}
        page = split_page(tmp_path, monkeypatch, [backwards])

        assert 'often caused by' not in page
        assert 'compiling from source' not in page

    def test_a_backwards_row_is_counted_and_reported_once(
            self, tmp_path, monkeypatch, caplog):
        """The summariser says it once, with the totals. The per-row branch must
        not say it again for every row, which is what a build emitting hundreds
        of them would otherwise produce."""
        hits = [
            {**mmseqs_hit(f'seq_{i}_test', f'seq_{i}_train'), 'tstart': 60, 'tend': 1}
            for i in range(3)
        ]
        with caplog.at_level(logging.WARNING, logger='genomic_benchmarks_qc'):
            split_page(tmp_path, monkeypatch, hits)

        backwards_warnings = [
            record for record in caplog.records
            if 'running backwards' in record.getMessage()
        ]
        assert len(backwards_warnings) == 1
        assert '3 of 3' in backwards_warnings[0].getMessage()
        assert 'Alignment visualisation failed' not in caplog.text
