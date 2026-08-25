"""Tests for the links a report carries out to the documentation.

A report is one self-contained file that gets shared by sending it, so these
links are the only route from a report back to the project - and they are the one
part of the page nothing else can check. A renamed heading in `docs/guide/`
leaves the report pointing at a URL that still loads and scrolls nowhere, which
no build failure and no broken-link checker would catch: the report is not part
of the site, and `mkdocs --strict` never sees it.

So the anchors are checked here, against the headings the guide actually has.
"""

import pathlib
import re

import pytest

from genomic_benchmarks_qc.report import utils
from genomic_benchmarks_qc.report.classes_html_report import EXPLANATION_LINKS

DOCS = pathlib.Path(__file__).resolve().parents[1] / 'docs'

# Every link a ? explanation carries, as (page key, anchor): the class report's
# map, plus the four the two reports build by hand.
LINKS = [(page, anchor) for page, anchor, _ in EXPLANATION_LINKS.values()] + [
    ('checks', 'per-position-nucleotide-content'),   # the forward position check
    ('viewer', None),                                # driving the same figure
    ('leakage', None),                               # the split report's check
    ('flags', None),                                 # from the verdict line
]


def _anchors(page_key):
    """The heading anchors of one guide page, slugified as mkdocs does.

    The site uses python-markdown's default `toc` slugify: lowercase, spaces to
    hyphens, everything else that is not a word character dropped.
    """
    relative = utils.DOCS_PAGES[page_key].rstrip('/') + '.md'
    source = DOCS / relative
    assert source.is_file(), f'{relative} is not a docs page'
    anchors = set()
    for line in source.read_text().splitlines():
        heading = re.match(r'#{2,6}\s+(.*?)\s*$', line)
        if not heading:
            continue
        text = re.sub(r'[`*_]', '', heading.group(1))
        anchors.add(re.sub(r'[^\w\s-]', '', text).strip().lower().replace(' ', '-'))
    return anchors


@pytest.mark.parametrize('page,anchor', LINKS)
def test_every_link_lands_on_a_heading_that_exists(page, anchor):
    if anchor is None:
        assert page in utils.DOCS_PAGES        # the page itself, no anchor
        return
    assert anchor in _anchors(page), (
        f'{utils.docs_url(page, anchor)} names a heading that '
        f'docs/{utils.DOCS_PAGES[page].rstrip("/")}.md does not have')


def test_the_urls_are_the_ones_the_project_publishes():
    """Hard-coded in one place, and the same strings the packaging declares."""
    pyproject = (DOCS.parent / 'pyproject.toml').read_text()

    assert f'Documentation = "{utils.DOCS_URL}"' in pyproject
    assert f'Repository = "{utils.REPO_URL}"' in pyproject


def test_a_docs_link_opens_in_a_new_tab_without_leaking_the_referrer():
    """A report is often read from a file:// path or a share; a link that
    replaced the page would lose the reader's place in a long report."""
    link = utils.docs_link('checks', 'unique-bases')

    assert 'target="_blank"' in link
    assert 'rel="noopener"' in link


def test_an_unknown_page_is_a_failure_not_a_bad_url():
    with pytest.raises(ValueError):
        utils.docs_url('nonexistent')
