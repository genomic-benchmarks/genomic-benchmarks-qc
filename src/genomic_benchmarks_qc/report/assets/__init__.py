"""The report's CSS and JavaScript, and the loader that inlines them.

The styling and the client-side behaviour live in real .css/.js files in this
package rather than in Python string literals, so an editor can highlight
them, a linter can read them, and `node --check` can be run over the scripts in
CI. They are still *inlined* into the generated page: a report has to be one
self-contained HTML file that can be mailed, archived or opened from a share
without dragging assets along behind it.

Anything read here must also be declared as package data in pyproject.toml,
otherwise it works from a source checkout and is missing from the wheel. The
test suite covers that by loading every asset through this module.
"""

from functools import lru_cache
from importlib.resources import files

from genomic_benchmarks_qc.report.colors import (
    FAIL_COLOR,
    PASS_COLOR,
    UNKNOWN_COLOR,
    WARN_COLOR,
)

# Colors are defined once in Python because the flags and the plots have to
# agree on them; the stylesheets take them through these placeholders.
COLOR_PLACEHOLDERS = {
    '{{FAIL_COLOR}}': FAIL_COLOR,
    '{{WARN_COLOR}}': WARN_COLOR,
    '{{PASS_COLOR}}': PASS_COLOR,
    '{{UNKNOWN_COLOR}}': UNKNOWN_COLOR,
}


@lru_cache(maxsize=None)
def read_asset(name: str) -> str:
    """Return the text of one asset file, with the color placeholders filled in.

    Args:
        name: File name inside the assets directory, e.g. 'report.css'.

    Returns:
        The file's contents.

    Raises:
        FileNotFoundError: If the asset is missing, which in an installed
            package means it was not declared as package data.
    """
    text = files(__name__).joinpath(name).read_text(encoding='utf-8')
    for placeholder, color in COLOR_PLACEHOLDERS.items():
        text = text.replace(placeholder, color)
    return text


def stylesheet(*names: str) -> str:
    """Return the named stylesheets concatenated, ready to drop inside <style>."""
    return '\n'.join(read_asset(name) for name in names)


def script(*names: str) -> str:
    """Return the named scripts, each wrapped in its own <script> element.

    One element per file keeps a syntax error in one script from taking the
    others down with it, and makes the page readable when debugging.
    """
    return '\n'.join(f'<script>\n{read_asset(name)}</script>' for name in names)
