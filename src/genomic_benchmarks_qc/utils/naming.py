"""The report layout: directory and file names, and how they are derived.

Every name the tool writes below the caller-provided output folder is defined
here, so the layout can be read in one place instead of being reassembled from
string literals scattered over the command and report modules.

Two kinds of name live here. The fixed ones - ``class``, ``split``,
``report.html`` - are plain constants. The derived ones come from user data
(class labels, input file names, sequence column names) and go through
``slugify``/``unique_slugs`` first, because report paths and report *contents*
have conflicting requirements: the HTML tables, plot legends and JSON reports
must show a class label exactly as it appears in the data, while a path segment
has to be filesystem-safe and unique within its parent directory. Only the path
form is produced here; the display form stays untouched on
``SequenceStatistics.label``.
"""

import logging
import re
from collections import Counter
from pathlib import Path

# Input extensions the CLI accepts, longest first so '.fasta' wins over '.fa'.
SEQUENCE_EXTENSIONS = ('.fasta', '.fa', '.csv', '.tsv')
COMPRESSION_EXTENSIONS = ('.gz',)

# Everything outside this set is replaced, which keeps path segments portable
# across filesystems. Dots and underscores survive because they are common and
# harmless in dataset names ('hg38.chr1', 'gene_expression').
_UNSAFE_RE = re.compile(r'[^A-Za-z0-9._-]+')
_DASH_RUN_RE = re.compile(r'-{2,}')

# Leave room for the suffixes callers append (e.g. '_vs_<other>') inside the
# 255-byte limit that common filesystems impose on a single path component.
MAX_SEGMENT_LENGTH = 80

FALLBACK_SLUG = 'unnamed'

# Reserved device names on Windows; cheap to avoid everywhere.
_RESERVED_NAMES = frozenset(
    ['con', 'prn', 'aux', 'nul']
    + [f'com{i}' for i in range(1, 10)]
    + [f'lpt{i}' for i in range(1, 10)]
)

# Sub-directory each command owns inside the caller-provided output folder.
# Everything above it - collection, dataset, split - is the caller's business.
CLASS_SUBDIR = 'class'
SPLIT_SUBDIR = 'split'

# One level below that: the sequence column an analysis was run on.
# FASTA inputs have no sequence column, but their reports still go one level
# deep so that the layout is identical for every input format. The CLI default
# column name is reused, which makes a FASTA dataset and a plain CSV one resolve
# to the same report path.
DEFAULT_COLUMN_DIR = 'sequence'

# Directory holding the analysis of all sequence columns concatenated together.
MERGED_COLUMN_DIR = 'merged'

# Fixed names inside a comparison directory. Because every comparison has a
# directory of its own, each report type can have the same name in all of them.
SIMPLE_REPORT_FILE = 'report.csv'
HTML_REPORT_FILE = 'report.html'
DUPLICATES_FILE = 'duplicates.txt'
PLOTS_DIR = 'plots'
MMSEQS_DIR = 'mmseqs'

# Prefix of the scratch directory holding the MMseqs2 working files. It is a
# prefix rather than a fixed name because the directory is created with
# `tempfile.mkdtemp`, which appends a random suffix: the create is atomic, so
# the run that gets the name is provably the only owner and may delete the
# directory afterwards without risking someone else's data. A fixed name would
# have to be adopted if it already existed, and removing it would then destroy
# whatever was there - including the files a concurrent run is still using.
TMP_PREFIX = 'gb-qc-mmseqs-'

# Per-class statistics are not a comparison, so they sit beside them.
PER_CLASS_DIR = 'per-class'


def strip_extensions(file_path):
    """Return the file name with only known input extensions removed.

    Stripping every entry of ``Path.suffixes`` would also eat dots that are
    part of the name itself, collapsing 'hg38.chr1.fa' and 'hg38.chr2.fa' to a
    shared 'hg38'. Genomic file names carry such dots routinely, so only the
    extensions the tool actually accepts are removed: 'hg38.chr1.fa' keeps its
    'hg38.chr1' stem.
    """
    name = Path(file_path).name

    for compression in COMPRESSION_EXTENSIONS:
        if name.lower().endswith(compression):
            name = name[: -len(compression)]
            break

    for extension in SEQUENCE_EXTENSIONS:
        if name.lower().endswith(extension):
            return name[: -len(extension)]

    return name


def slugify(value, fallback=FALLBACK_SLUG):
    """Reduce arbitrary text to a lowercase, filesystem-safe path segment.

    Lowercasing is deliberate: on case-insensitive filesystems labels 'A' and
    'a' would otherwise resolve to the same directory, so they are folded here
    and disambiguated by ``unique_slugs`` instead.
    """
    text = _UNSAFE_RE.sub('-', str(value).strip().lower())
    text = _DASH_RUN_RE.sub('-', text).strip('-._')[:MAX_SEGMENT_LENGTH].strip('-._')

    if not text or text in _RESERVED_NAMES:
        return fallback

    return text


def _fit(text, reserved=0):
    """Trim a slug to `MAX_SEGMENT_LENGTH`, leaving room for `reserved` more chars.

    `slugify` bounds each piece on its own, so combining a context slug with a
    value slug - or appending a numeric suffix - can produce twice the limit.
    Two such names then meet in `comparison_dirname`, and the result exceeds the
    255-character path component that common filesystems allow, so creating the
    report directory fails. Trimming here keeps every slug within the budget the
    limit was chosen for.
    """
    return text[:MAX_SEGMENT_LENGTH - reserved].strip('-._')


def unique_slugs(values, contexts=None):
    """Slugify ``values``, disambiguating any slugs that end up identical.

    Distinct inputs can share a slug either because they genuinely share a name
    ('train/pos.fa' and 'test/pos.fa') or because slugifying folded them
    together ("5' UTR" and "5'-UTR"). Both cases would otherwise write their
    reports to the same directory.

    ``contexts`` optionally supplies a disambiguating hint per value - for
    input files, the parent directory name - which is preferred over a bare
    numeric suffix because it says *which* file a report came from. A numeric
    suffix is the fallback when no context is available or it does not help.

    Returns one slug per value, in the order given.
    """
    values = list(values)
    slugs = [slugify(value) for value in values]
    contexts = list(contexts) if contexts is not None else [None] * len(values)

    counts = Counter(slugs)

    # Values whose slug is already unique keep it, and it is reserved up front:
    # disambiguating some other value must not take a name an input already owns.
    unique = [slug if counts[slug] == 1 else None for slug in slugs]
    used = {slug for slug in unique if slug is not None}

    for index, (value, slug, context) in enumerate(zip(values, slugs, contexts)):
        if unique[index] is not None:
            continue

        candidate = slug
        if context is not None:
            context_slug = slugify(context, fallback='')
            if context_slug:
                candidate = _fit(f'{context_slug}-{slug}')

        base, suffix = candidate, 1
        while candidate in used:
            suffix += 1
            marker = f'-{suffix}'
            candidate = f'{_fit(base, len(marker))}{marker}'

        if candidate != slug:
            logging.warning(
                f"Name '{slug}' derived from '{value}' is not unique among the inputs; "
                f"using '{candidate}' for its report path instead."
            )

        used.add(candidate)
        unique[index] = candidate

    return unique


def comparison_dirname(slug1, slug2):
    """Return the directory name holding the comparison of two slugs."""
    return f'{slug1}_vs_{slug2}'


def column_dirname(input_format, sequence_column):
    """Return the column directory for an analysis reading all columns at once.

    Used where one report covers every sequence column together, as the split
    evaluation does by concatenating them into a single FASTA. FASTA inputs have
    no column at all and fall back to `DEFAULT_COLUMN_DIR`, which keeps the
    layout identical across input formats.
    """
    if str(input_format).startswith('fa') or not sequence_column:
        return DEFAULT_COLUMN_DIR

    if len(sequence_column) > 1:
        return MERGED_COLUMN_DIR

    return slugify(sequence_column[0])


def per_column_dirnames(sequence_column):
    """Return the column directories for one analysis per sequence column.

    Yields ``(per_column, merged)``: a directory name for each column in the
    order given, and the name for the extra analysis of all columns
    concatenated - or None when there is only one column and no such analysis.

    The merged name is slugged together with the column names so that a column
    literally called 'merged' keeps its own directory and the merged analysis
    gets a disambiguated one, rather than the two overwriting each other.
    """
    names = list(sequence_column)

    if len(names) > 1:
        slugs = unique_slugs(names + [MERGED_COLUMN_DIR])
        return slugs[:-1], slugs[-1]

    return unique_slugs(names), None
