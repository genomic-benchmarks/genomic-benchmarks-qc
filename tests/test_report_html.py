"""Tests for what ends up in the generated HTML.

The per-position panels are not embedded as a picture: the page carries the
numbers and draws them in a canvas, so a reader can zoom into a flagged position
instead of squinting at a 1px band. Two properties have to hold for that to be
worth anything - the data has to be in the page, and the page has to stay a
single standalone file that opens with no network and no sibling assets.
"""

import json
import re
from urllib.parse import urlparse

import numpy as np
import pandas as pd
import pytest
from helpers import mmseqs_hit, write_csv, write_mmseqs_output

from genomic_benchmarks_qc import evaluate_splits
from genomic_benchmarks_qc.report.report_generator import generate_dataset_html_report
from genomic_benchmarks_qc.report.split_html_report import alignments_count_text
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

    def test_a_base_that_is_a_metacharacter_stays_text(self, tmp_path):
        """Unusual bases are exactly what the Unique bases check exists to
        report, so the check must survive being handed one."""
        page = render(tmp_path,
                      make_stats(['AC<GT' * 4] * 30, label='a'),
                      make_stats(['ACGT' * 5] * 30, label='b'))

        # the cell listing the bases survives intact, with the odd one escaped
        assert '<td style="text-align: center;">&lt;, A, C, G, T</td>' in page
        assert '<td style="text-align: center;"><, A, C, G, T</td>' not in page
