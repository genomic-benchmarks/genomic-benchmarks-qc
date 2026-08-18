"""Unit tests for genomic_benchmarks_qc.utils.naming."""

import pytest

from genomic_benchmarks_qc.utils.naming import (
    DEFAULT_COLUMN_DIR,
    FALLBACK_SLUG,
    MAX_SEGMENT_LENGTH,
    MERGED_COLUMN_DIR,
    column_dirname,
    comparison_dirname,
    per_column_dirnames,
    slugify,
    strip_extensions,
    unique_slugs,
)


class TestStripExtensions:
    @pytest.mark.parametrize(
        "file_path, expected",
        [
            ("pos.fa", "pos"),
            ("pos.fasta", "pos"),
            ("pos.csv", "pos"),
            ("pos.tsv", "pos"),
            ("pos.fa.gz", "pos"),
            ("pos.fasta.gz", "pos"),
            ("pos.csv.gz", "pos"),
            # Extensions are matched case-insensitively, the stem keeps its case.
            ("Pos.FASTA", "Pos"),
            ("Pos.Fa.GZ", "Pos"),
            # Only the file name matters, never the directories above it.
            ("train/pos.fa", "pos"),
            ("/abs/results.v2/pos.fa", "pos"),
            # A file with no known extension keeps its whole name.
            ("dataset", "dataset"),
        ],
    )
    def test_strips_known_extensions(self, file_path, expected):
        assert strip_extensions(file_path) == expected

    @pytest.mark.parametrize(
        "file_path, expected",
        [
            ("hg38.chr1.fa", "hg38.chr1"),
            ("hg38.chr2.fa", "hg38.chr2"),
            ("enhancers.v2.fa", "enhancers.v2"),
            ("hg38.chr1.train.csv.gz", "hg38.chr1.train"),
        ],
    )
    def test_keeps_dots_that_are_part_of_the_name(self, file_path, expected):
        """Only known extensions are removed, so dotted names stay distinct.

        Stripping every suffix would collapse 'hg38.chr1.fa' and 'hg38.chr2.fa'
        onto a shared 'hg38' and make unrelated files share a report path.
        """
        assert strip_extensions(file_path) == expected


class TestSlugify:
    @pytest.mark.parametrize(
        "value, expected",
        [
            ("pos", "pos"),
            ("high", "high"),
            ("0", "0"),
            # Characters that are unusable or awkward in a path are replaced.
            ("A/B", "a-b"),
            ("5' UTR", "5-utr"),
            ("non-coding RNA", "non-coding-rna"),
            ("a\\b:c*d", "a-b-c-d"),
            # Dots and underscores are common in dataset names and survive.
            ("hg38.chr1", "hg38.chr1"),
            ("gene_expression", "gene_expression"),
            # Runs of replaced characters collapse, and edges are trimmed.
            ("a   ///   b", "a-b"),
            ("  spaced  ", "spaced"),
            (".hidden.", "hidden"),
        ],
    )
    def test_produces_safe_segments(self, value, expected):
        assert slugify(value) == expected

    def test_folds_case_so_labels_cannot_collide_on_case_insensitive_filesystems(self):
        assert slugify("Positives") == slugify("positives") == "positives"

    @pytest.mark.parametrize("value", ["", "   ", "...", "///", "con", "LPT1"])
    def test_falls_back_when_nothing_usable_survives(self, value):
        assert slugify(value) == FALLBACK_SLUG

    def test_custom_fallback_is_returned_verbatim(self):
        assert slugify("///", fallback="") == ""

    def test_truncates_to_a_length_filesystems_accept(self):
        slug = slugify("x" * 500)
        assert slug == "x" * MAX_SEGMENT_LENGTH

    def test_truncation_does_not_leave_a_trailing_separator(self):
        # Truncating mid-way through replaced characters must not end in '-'.
        slug = slugify("y" * (MAX_SEGMENT_LENGTH - 1) + " tail")
        assert not slug.endswith(("-", ".", "_"))


class TestUniqueSlugs:
    def test_leaves_distinct_values_untouched(self):
        assert unique_slugs(["pos", "neg"]) == ["pos", "neg"]

    def test_disambiguates_shared_names_with_their_context(self):
        slugs = unique_slugs(["pos", "pos"], contexts=["train", "test"])
        assert slugs == ["train-pos", "test-pos"]

    def test_disambiguates_with_a_numeric_suffix_when_no_context_helps(self):
        assert unique_slugs(["pos", "pos"]) == ["pos", "pos-2"]

    def test_disambiguates_values_that_only_collide_after_slugifying(self):
        # Distinct labels in the data, one slug: they still need separate paths.
        assert unique_slugs(["5' UTR", "5'-UTR"]) == ["5-utr", "5-utr-2"]

    def test_falls_back_to_numeric_suffix_when_contexts_also_collide(self):
        slugs = unique_slugs(["pos", "pos"], contexts=["same", "same"])
        assert slugs == ["same-pos", "same-pos-2"]

    def test_context_is_ignored_when_it_is_not_usable(self):
        # An empty parent directory name (a bare relative path) offers no hint.
        assert unique_slugs(["pos", "pos"], contexts=["", ""]) == ["pos", "pos-2"]

    def test_disambiguated_name_does_not_steal_another_input_slug(self):
        # 'train-pos' is already taken by the third input, so the disambiguated
        # first input must move further rather than overwrite it.
        slugs = unique_slugs(["pos", "pos", "train-pos"], contexts=["train", "test", None])
        assert len(set(slugs)) == 3
        assert slugs[2] == "train-pos"

    def test_preserves_input_order(self):
        assert unique_slugs(["c", "a", "b"]) == ["c", "a", "b"]

    def test_all_results_are_unique_for_a_pathological_input(self):
        values = ["A/B", "a-b", "a b", "A B", "a__b"]
        slugs = unique_slugs(values)
        assert len(set(slugs)) == len(values)

    def test_accepts_no_contexts(self):
        assert unique_slugs([]) == []


class TestComparisonDirname:
    def test_joins_two_slugs(self):
        assert comparison_dirname("pos", "neg") == "pos_vs_neg"


class TestColumnDirname:
    """The single column directory used when all columns are read at once."""

    @pytest.mark.parametrize("input_format", ['fasta', 'fa', 'fasta.gz', 'fa.gz'])
    def test_fasta_inputs_use_the_default_column_name(self, input_format):
        # FASTA has no columns at all, but its reports still go one level deep.
        assert column_dirname(input_format, ['sequence']) == DEFAULT_COLUMN_DIR
        assert column_dirname(input_format, ['anything', 'else']) == DEFAULT_COLUMN_DIR

    @pytest.mark.parametrize("input_format", ['csv', 'tsv', 'csv.gz', 'tsv.gz'])
    def test_a_single_column_names_the_directory(self, input_format):
        assert column_dirname(input_format, ['gene']) == 'gene'

    def test_several_columns_are_reported_as_merged(self):
        # They are concatenated into one analysis, which is what 'merged' means.
        assert column_dirname('csv', ['gene', 'noncodingRNA']) == MERGED_COLUMN_DIR

    def test_a_column_name_unusable_in_a_path_is_slugged(self):
        assert column_dirname('csv', ["5' UTR"]) == '5-utr'

    def test_a_missing_column_falls_back_to_the_default(self):
        assert column_dirname('csv', []) == DEFAULT_COLUMN_DIR
        assert column_dirname('csv', None) == DEFAULT_COLUMN_DIR


class TestPerColumnDirnames:
    """One directory per column, plus the merged one when there are several."""

    def test_a_single_column_has_no_merged_analysis(self):
        assert per_column_dirnames(['sequence']) == (['sequence'], None)

    def test_several_columns_each_get_a_directory_and_a_merged_one(self):
        assert per_column_dirnames(['gene', 'noncodingRNA']) == (
            ['gene', 'noncodingrna'],
            MERGED_COLUMN_DIR,
        )

    def test_column_order_is_preserved(self):
        per_column, _ = per_column_dirnames(['b', 'a', 'c'])
        assert per_column == ['b', 'a', 'c']

    def test_a_column_named_merged_keeps_its_own_directory(self):
        # The column owns the name it is called; the merged analysis gives way.
        per_column, merged = per_column_dirnames(['gene', 'merged'])
        assert per_column == ['gene', 'merged']
        assert merged != 'merged'
        assert merged not in per_column

    def test_columns_that_slug_alike_get_separate_directories(self):
        per_column, _ = per_column_dirnames(["5' UTR", '5 UTR'])
        assert len(set(per_column)) == 2
