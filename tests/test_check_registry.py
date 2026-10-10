"""Tests for the check registries: what a check has to be for the rest to work.

Everything that used to name the checks one by one - the order of the results,
the rows of `gb-qc-report.csv`, the figures' shading, the report's navigation and
sections - now reads the registry instead. So the registry is where the rules a
check has to follow are pinned: names that the sub-check convention can tell
apart, sections whose ids cannot collide, and a results table whose shape does
not depend on which check happens to come first.

And the reason there is a registry at all: a check listed in it is in every
report, with nothing else edited. The last two classes add one to each registry
and look for it everywhere it should appear.
"""

import csv
import re

import pandas as pd
import pytest
from helpers import mmseqs_hit, sequences, write_csv, write_mmseqs_output

from genomic_benchmarks_qc import evaluate_classes, evaluate_splits
from genomic_benchmarks_qc.checks import Check, CheckResult
from genomic_benchmarks_qc.checks import classes as class_checks
from genomic_benchmarks_qc.checks import splits as split_checks
from genomic_benchmarks_qc.checks.scorers import scalar_feature
from genomic_benchmarks_qc.report import assets
from genomic_benchmarks_qc.report.sections import Section, SectionContent
from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics
from genomic_benchmarks_qc.utils.testing import flag_significant_differences

REGISTRIES = {'classes': class_checks.CLASS_CHECKS, 'splits': split_checks.SPLIT_CHECKS}

# Ids the page chrome uses itself, which no section may take.
CHROME_IDS = {'basic-descriptive-statistics'}


def make_stats(sequences, label):
    stats = SequenceStatistics(sequences=sequences, filename=f'{label}.fa',
                               filepath=f'/{label}.fa', label=label)
    stats.compute()
    return stats


def comparison(count=300, length=40):
    return (make_stats(sequences(count, length, seed=1), 'a'),
            make_stats(sequences(count, length, seed=2), 'b'))


@pytest.mark.parametrize('registry', REGISTRIES.values(), ids=REGISTRIES.keys())
class TestEveryRegistry:
    def test_every_check_has_a_name_of_its_own(self, registry):
        names = [check.name for check in registry]

        assert len(names) == len(set(names))

    def test_no_name_looks_like_a_sub_check(self, registry):
        """' - ' is what separates a sub-check from its headline, in the results
        and in every consumer of `gb-qc-report.csv` that filters on it."""
        for check in registry:
            assert ' - ' not in check.name, check.name

    def test_no_two_sections_share_an_id(self, registry):
        """The navigation links to a section by its anchor and the ? button finds
        its explanation by id; a duplicate would send both to the first match."""
        ids = [id_ for check in registry
               for id_ in (check.section.anchor, check.section.explanation_id)]

        assert len(ids) == len(set(ids))
        assert not set(ids) & CHROME_IDS

    def test_every_section_has_its_fragment(self, registry):
        for check in registry:
            assert assets.template(check.section.template).strip(), check.name


class TestClassRegistry:
    def test_every_section_is_named_as_the_guide_names_it(self):
        """A section's id is its check's heading anchor in docs/guide/checks.md,
        the page its first docs link points to."""
        for check in class_checks.CLASS_CHECKS:
            page, anchor, _ = check.section.docs[0]
            assert (page, anchor) == ('checks', check.section.anchor), check.name

    def test_the_duplicate_listing_is_found_by_its_section_id(self):
        """report_ui.js fills the table of shared sequences, and the stylesheets
        lay it out, through the section's id."""
        anchors = {check.name: check.section.anchor for check in class_checks.CLASS_CHECKS}
        anchor = anchors['Duplicate Sequences between Labels']

        assert f"'#{anchor} tbody'" in assets.read_asset('report_ui.js')
        for stylesheet in ('report.css', 'report_design.css'):
            assert f'#{anchor} table' in assets.read_asset(stylesheet), stylesheet

    def test_every_row_belongs_to_exactly_one_check(self):
        """A row is a headline or one of its sub-checks, and no check's rows can be
        mistaken for another's: the report groups them by that prefix."""
        results, _ = flag_significant_differences(*comparison())
        names = [check.name for check in class_checks.CLASS_CHECKS]

        for row in results:
            owners = [name for name in names if row == name or row.startswith(f'{name} - ')]
            assert len(owners) == 1, (row, owners)

    def test_the_headlines_come_first_in_registry_order(self):
        results, _ = flag_significant_differences(*comparison())
        names = [check.name for check in class_checks.CLASS_CHECKS]

        assert list(results)[:len(names)] == names

    def test_every_check_with_a_figure_to_shade_hands_it_over(self):
        _, failed_by_feature = flag_significant_differences(*comparison())

        assert list(failed_by_feature) == [
            'Per sequence nucleotide content',
            'Per sequence dinucleotide content',
            'Per position nucleotide content',
            'Per position reversed nucleotide content',
        ]


class TestSimpleReport:
    def test_the_flag_is_the_first_column(self, tmp_path):
        """The header used to come out this way only because the first check in the
        table happened to carry nothing but a flag."""
        rows = ['sequence,label'] + [f'{seq},{index % 2}'
                                     for index, seq in enumerate(sequences(600, 40))]
        data = tmp_path / 'in.csv'
        data.write_text('\n'.join(rows) + '\n')

        evaluate_classes.run(input=[str(data)], format='csv', out_folder=str(tmp_path / 'out'),
                             report_types=['simple'])
        report = next((tmp_path / 'out').rglob('gb-qc-report.csv'))

        with report.open(newline='') as handle:
            header = next(csv.reader(handle))
        assert header == ['Check', 'Flag', 'AU-ROC', 'AU-PR', 'Accuracy']

    def test_a_split_report_has_the_same_first_columns(self, tmp_path, monkeypatch):
        report = split_report(tmp_path, monkeypatch, [], report_types=['simple'])

        with report.with_name('gb-qc-report.csv').open(newline='') as handle:
            header = next(csv.reader(handle))
        assert header == ['Check', 'Flag', 'Percentage of leaked queries',
                          'Percentage of leaked targets']


def split_report(tmp_path, monkeypatch, hits, report_types=('html', 'simple')):
    """Run the split command against a fixed set of MMseqs2 hits; the HTML path."""
    train = write_csv(tmp_path / 'train.csv', ['0'], rows_per_label=5)
    test = write_csv(tmp_path / 'test.csv', ['0'], rows_per_label=5)

    def fake_run_search(query_fasta, target_fasta, output_path, tmp_dir, **kwargs):
        return write_mmseqs_output(output_path, hits)

    monkeypatch.setattr(evaluate_splits.mmseqs_runtime, 'run_search', fake_run_search)
    evaluate_splits.run(train_files=[train], test_files=[test], format='csv',
                        out_folder=str(tmp_path / 'out'), report_types=list(report_types))
    return tmp_path / 'out' / 'split' / 'sequence' / 'train_vs_test' / 'gb-qc-report.html'


# A check that is not part of the tool: the number of CpG steps in a sequence,
# scored like any other feature. It computes its own feature from the sequences,
# so `SequenceStatistics` knows nothing about it either.
CPG = 'Per sequence CpG count'
CPG_FRAGMENT = """<div id="{{explanation_id}}" class="explanation-text">
    <p>How many CG steps each sequence has. {{docs_link}}</p>
</div>
<p class="cpg-body">{{cpg_body}}</p>
"""


def cpg_steps(stats):
    return pd.DataFrame({CPG: [float(seq.count('CG')) for seq in stats.sequences]})


def render_cpg(context):
    return SectionContent({'{{cpg_body}}': 'CpG steps compared'})


CPG_CHECK = Check(
    name=CPG,
    score=scalar_feature(CPG, cpg_steps),
    floor='per_sequence',
    section=Section(
        title='Per Sequence CpG Count',
        anchor='per-sequence-cpg-count',
        explanation_id='cpg-explanation',
        template='check_cpg_count.html',
        render=render_cpg,
        docs=(('checks', None, 'What to do about it'),),
    ),
)


def serve_fragment(monkeypatch, name, markup):
    """Make one more fragment loadable, as a new check's .html file would be."""
    template = assets.template
    monkeypatch.setattr(assets, 'template',
                        lambda asset: markup if asset == name else template(asset))


class TestANewClassCheckNeedsNothingElse:
    @pytest.fixture
    def outputs(self, tmp_path, monkeypatch):
        monkeypatch.setattr(class_checks, 'CLASS_CHECKS', class_checks.CLASS_CHECKS + (CPG_CHECK,))
        serve_fragment(monkeypatch, 'check_cpg_count.html', CPG_FRAGMENT)

        # One class rich in CG steps, so the check has a difference to find.
        plain = sequences(300, 40, seed=1)
        rich = [seq[:20] + 'CGCGCGCG' + seq[28:] for seq in sequences(300, 40, seed=2)]
        rows = ['sequence,label'] + [f'{seq},plain' for seq in plain] + [f'{seq},rich' for seq in rich]
        data = tmp_path / 'in.csv'
        data.write_text('\n'.join(rows) + '\n')

        evaluate_classes.run(input=[str(data)], format='csv', out_folder=str(tmp_path / 'out'))
        report_dir = tmp_path / 'out' / 'class' / 'sequence' / 'plain_vs_rich'
        with (report_dir / 'gb-qc-report.csv').open(newline='') as handle:
            flags = {row['Check']: row['Flag'] for row in csv.DictReader(handle)}
        return flags, (report_dir / 'gb-qc-report.html').read_text()

    def test_it_is_scored_and_flagged_like_any_other_check(self, outputs):
        flags, _ = outputs

        assert flags[CPG] == 'Fail'

    def test_it_is_in_the_navigation_and_has_a_section(self, outputs):
        _, page = outputs

        assert '<a href="#per-sequence-cpg-count">Per Sequence CpG Count</a>' in page
        section = re.search(r'<section id="per-sequence-cpg-count">.*?</section>', page, re.S)
        assert section, 'no section for the new check'
        assert 'status-fail' in section.group(0)
        assert "toggleExplanation('cpg-explanation')" in section.group(0)
        assert 'CpG steps compared' in section.group(0)

    def test_it_is_counted_in_the_verdict(self, outputs):
        _, page = outputs

        assert '10 checks:' in page

    def test_nothing_is_left_unfilled(self, outputs):
        _, page = outputs

        assert re.findall(r'\{\{\w+\}\}', page) == []


class TestANewSplitCheckNeedsNothingElse:
    """The same holds for the split report: listed, and it is in the page and the CSV."""

    def test_it_is_in_the_report_and_the_csv(self, tmp_path, monkeypatch):
        name = 'Queries without a hit'

        def score(threshold_stats):
            share = threshold_stats['num_queries_without_hits'] / threshold_stats['num_all_queries']
            return CheckResult({name: {'Flag': 'Warning' if share < 1 else 'Pass'}})

        check = Check(name=name, score=score, section=Section(
            title=name, anchor='queries-without-hits', explanation_id='no-hit-explanation',
            template='check_no_hits.html', render=lambda context: SectionContent(),
            docs=(('leakage', None, 'Read more'),)))
        monkeypatch.setattr(split_checks, 'SPLIT_CHECKS', split_checks.SPLIT_CHECKS + (check,))
        serve_fragment(monkeypatch, 'check_no_hits.html',
                       '<div id="{{explanation_id}}" class="explanation-text">'
                       '<p>Test sequences MMseqs2 found nothing for. {{docs_link}}</p></div>')

        page = split_report(tmp_path, monkeypatch, [mmseqs_hit('seq_0_test', 'seq_0_train')])
        with page.with_name('gb-qc-report.csv').open(newline='') as handle:
            flags = {row['Check']: row['Flag'] for row in csv.DictReader(handle)}
        page = page.read_text()

        assert flags[name] == 'Warning'
        assert '<a href="#queries-without-hits">Queries without a hit</a>' in page
        assert '<section id="queries-without-hits">' in page
        assert '2 checks:' in page
