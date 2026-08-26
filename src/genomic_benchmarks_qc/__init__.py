"""Find the shortcut. Learn the biology.

Automated quality control for genomic machine learning datasets: scores the
biases, duplicates and data leakage a classifier could exploit before you train
on it.
"""

import logging

# The one place the version is written down. pyproject.toml reads it from here
# when building, so the package metadata cannot drift from what the reports show,
# and reading it back at runtime costs nothing and does not depend on the package
# having been installed.
__version__ = "1.0.0"


# A library configures nothing on import. The null handler is what keeps a
# record from an unconfigured `genomic_benchmarks_qc` logger quiet instead of
# reaching logging's last-resort handler on stderr; `utils.input_utils.
# setup_logger` is what puts real handlers on it, and only when asked.
logging.getLogger(__name__).addHandler(logging.NullHandler())
