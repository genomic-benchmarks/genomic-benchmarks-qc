"""Tests for the check registries: what a check has to be for the rest to work.

Everything that used to name the checks one by one - the order of the results,
the rows of `gb-qc-report.csv`, the figures' shading, the report's navigation and
sections - now reads the registry instead. So the registry is where the rules a
check has to follow are pinned: names that the sub-check convention can tell
apart, and a results table whose shape does not depend on which check happens to
come first.
"""

import csv

from helpers import sequences

from genomic_benchmarks_qc import evaluate_classes
from genomic_benchmarks_qc.checks import classes as class_checks
from genomic_benchmarks_qc.utils.seq_stats import SequenceStatistics
from genomic_benchmarks_qc.utils.testing import flag_significant_differences


def make_stats(sequences, label):
    stats = SequenceStatistics(sequences=sequences, filename=f'{label}.fa',
                               filepath=f'/{label}.fa', label=label)
    stats.compute()
    return stats


def comparison(count=300, length=40):
    return (make_stats(sequences(count, length, seed=1), 'a'),
            make_stats(sequences(count, length, seed=2), 'b'))


class TestNames:
    def test_every_check_has_a_name_of_its_own(self):
        names = [check.name for check in class_checks.CLASS_CHECKS]

        assert len(names) == len(set(names))

    def test_no_name_looks_like_a_sub_check(self):
        """' - ' is what separates a sub-check from its headline, in the results
        and in every consumer of `gb-qc-report.csv` that filters on it."""
        for check in class_checks.CLASS_CHECKS:
            assert ' - ' not in check.name, check.name

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
