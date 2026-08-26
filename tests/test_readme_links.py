"""Tests that the README's links still work where the README is read.

The README is not only a file in this repository. setuptools puts it into the
package metadata as the long description, so it is also the PyPI project page -
and there a relative link resolves against `pypi.org`, not against the
repository, and 404s. The same goes for an in-page anchor: PyPI's renderer
strips the `id` attributes off headings, so `#something` lands nowhere.

Nothing in a normal checkout notices, because in a checkout every one of those
links is correct. This is the check that does.
"""

import pathlib
import re

README = pathlib.Path(__file__).resolve().parents[1] / 'README.md'
TEXT = README.read_text()

# `[text](target)`, and the `href=`/`src=` of the raw HTML the README also uses
# for its logo, badges and the two recorded animations.
MARKDOWN_LINK = re.compile(r'\[[^\]]*\]\(\s*([^)\s]+)')
HTML_TARGET = re.compile(r'(?:href|src)="([^"]+)"')

ABSOLUTE = ('http://', 'https://', 'mailto:')


def _targets():
    return sorted(set(MARKDOWN_LINK.findall(TEXT)) | set(HTML_TARGET.findall(TEXT)))


def test_every_link_is_absolute():
    relative = [target for target in _targets() if not target.startswith(ABSOLUTE)]
    assert not relative, (
        "README.md is the PyPI project page as well as a file in this repository, "
        "and these targets only resolve in the repository: "
        f"{', '.join(relative)}. Write them as full https:// URLs into "
        "github.com/genomic-benchmarks/genomic-benchmarks-qc.")


def test_the_quick_start_does_not_assume_a_checkout():
    """`pip install` brings no `examples/` directory, and Quick Start follows it.

    The rest of the README may name repository paths - the Input Formats section
    does, to show what an option looks like - but the commands directly under
    the install instructions are the ones a reader runs first, and they have to
    work with nothing but the installed package and what the README itself says
    to fetch.
    """
    quick_start = TEXT.split('## Quick Start', 1)[1].split('\n## ', 1)[0]

    # A bare path, not the same segment inside the raw.githubusercontent.com URL
    # the section uses to fetch those very files - which is the fix, not the bug.
    bare_path = re.search(r'(?<![\w/])examples/', quick_start)

    assert bare_path is None, (
        "The Quick Start names a path from this repository. A reader arriving "
        "from `pip install genomic-benchmarks-qc` has no `examples/` directory, "
        "so the first command they run fails.")
