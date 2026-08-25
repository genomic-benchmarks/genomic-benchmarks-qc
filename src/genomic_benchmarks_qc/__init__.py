"""Find the shortcut. Learn the biology.

Automated quality control for genomic machine learning datasets: scores the
biases, duplicates and data leakage a classifier could exploit before you train
on it.
"""

# The one place the version is written down. pyproject.toml reads it from here
# when building, so the package metadata cannot drift from what the reports show,
# and reading it back at runtime costs nothing and does not depend on the package
# having been installed.
__version__ = "0.9.0"
