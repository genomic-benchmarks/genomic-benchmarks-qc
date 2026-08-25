"""Tests for what a report carries, and what it leaves on disk beside it.

A report is one standalone file, so everything it shows is inlined into it and
every byte of a figure is a byte of the file. Three things keep that from
running away, and all three are invisible in the output: a flagged figure is
drawn once and saved twice rather than built twice; the copy the page embeds is
half the resolution of the one written to plots/, which is the copy to reuse
elsewhere; and the logo is written into the stylesheet once instead of into
each of the two elements that show it.
"""

import base64
import io
import re

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from genomic_benchmarks_qc.report.report_generator import generate_dataset_html_report
from genomic_benchmarks_qc.report.utils import DISPLAY_DPI, FIGURE_DPI, LOGO_BASE64
from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics
from genomic_benchmarks_qc.utils.testing import flag_significant_differences

# The widest a report ever shows the logo: the header's copy, sized in
# report_design.css. The nav's is 132px.
WIDEST_LOGO_DISPLAY = 250


def make_stats(sequences, label='cls'):
    stats = SequenceStatistics(sequences=sequences, filename=f'{label}.fa',
                               filepath=f'/{label}.fa', label=label)
    stats.compute()
    return stats


def random_sequences(count, length, seed, composition=None):
    rng = np.random.default_rng(seed)
    return [''.join(rng.choice(list('ACGT'), size=length, p=composition)) for _ in range(count)]


@pytest.fixture(scope='module')
def report(tmp_path_factory):
    """A report on a comparison that flags, so both versions of a figure exist."""
    out = tmp_path_factory.mktemp('weight')
    stats1 = make_stats(random_sequences(300, 30, seed=1), label='a')
    stats2 = make_stats(random_sequences(300, 30, seed=2, composition=[.7, .1, .1, .1]),
                        label='b')
    results, failed_by_feature = flag_significant_differences(stats1, stats2)
    generate_dataset_html_report(
        stats1, stats2, out / 'gb-qc-report.html', plots_path=out / 'plots',
        plot_type='boxen',
        results=pd.DataFrame.from_dict(results, orient='index'),
        failed_by_feature=failed_by_feature,
    )
    return out


def embedded_images(page):
    """Every inlined PNG in the page, decoded, largest first."""
    found = re.findall(r'data:image/png;base64,\s*([A-Za-z0-9+/=\s]+?)["\')]', page)
    images = [base64.b64decode(raw.strip()) for raw in found]
    return sorted(images, key=len, reverse=True)


class TestBothVersionsComeFromOneFigure:
    """The clean file and the flagged one are the same figure, saved twice."""

    @pytest.mark.parametrize('stem', ['per_sequence_nucleotide_content',
                                      'per_sequence_dinucleotide_content',
                                      'per_position_nucleotide_content'])
    def test_the_pair_differs_only_by_the_marks(self, report, stem):
        plain = report / 'plots' / f'{stem}.png'
        flagged = report / 'plots' / f'{stem}_with_flags.png'
        assert flagged.is_file(), 'this comparison flags every one of these'

        # Same crop, so the boxes and curves sit in the same places and the two
        # files can be read side by side. `reserve_flag_margin` is what holds
        # this: the marks are drawn outside the axes and would otherwise grow
        # the tight bounding box.
        assert Image.open(plain).size == Image.open(flagged).size
        assert plain.read_bytes() != flagged.read_bytes(), 'the marks are missing'

    def test_the_page_shows_the_flagged_one(self, report):
        page = (report / 'gb-qc-report.html').read_text()
        flagged = report / 'plots' / 'per_sequence_nucleotide_content_with_flags.png'

        # Not by bytes - the page carries the smaller copy - but the flagged
        # figure is the wider of the two, so it is the one with more ink.
        assert flagged.stat().st_size > (
            report / 'plots' / 'per_sequence_nucleotide_content.png').stat().st_size
        assert len(embedded_images(page)) >= 3


class TestTheEmbeddedCopyIsSmallerThanTheFile:
    def test_the_file_is_at_print_resolution_and_the_page_at_half(self, report):
        page = (report / 'gb-qc-report.html').read_text()
        on_disk = Image.open(report / 'plots' / 'per_sequence_dinucleotide_content_with_flags.png')
        # The dinucleotide figure is the tallest thing in the report, so it is
        # the largest inlined image.
        in_page = Image.open(io.BytesIO(embedded_images(page)[0]))

        assert on_disk.width / in_page.width == pytest.approx(FIGURE_DPI / DISPLAY_DPI, abs=0.02)
        assert on_disk.height / in_page.height == pytest.approx(FIGURE_DPI / DISPLAY_DPI, abs=0.02)

    def test_no_inlined_figure_is_larger_than_its_file(self, report):
        page = (report / 'gb-qc-report.html').read_text()
        largest_embedded = len(embedded_images(page)[0])
        largest_file = max(png.stat().st_size for png in (report / 'plots').glob('*.png'))

        assert largest_embedded < largest_file


class TestTheLogo:
    def test_it_is_written_into_the_page_once(self, report):
        page = (report / 'gb-qc-report.html').read_text()

        # Two elements show it, and both read the one copy in the stylesheet.
        # `put_data` fills every occurrence of a placeholder, so a data URI in
        # the markup was embedded once per element.
        assert page.count(LOGO_BASE64) == 1
        assert len(re.findall(r'class="gb-qc-logo[^"]*"', page)) == 2

    def test_it_is_stored_near_the_size_it_is_shown_at(self):
        logo = Image.open(io.BytesIO(base64.b64decode(LOGO_BASE64.split(',', 1)[1])))

        # Twice the widest display, which is what a 2x screen asks for. It was
        # nearly six times that, in a report that showed it twice.
        assert logo.width == 2 * WIDEST_LOGO_DISPLAY
