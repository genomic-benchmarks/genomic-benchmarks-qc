"""genbenchQC has been renamed to genomic-benchmarks-qc.

1.2.0 is the final release under this name, and exists only to say so.
"""

import warnings

# Python's default filter shows a DeprecationWarning only when it is attributed
# to __main__, and a warning raised while a package is being imported is
# attributed to the import machinery instead - no stacklevel reaches past that,
# so this line is invisible unless the caller runs with -W default, or under
# pytest, which turns it on. That is what it is: the channel that reaches people
# running the command is the notice printed in cli.py, and this one is here for
# anyone importing genbenchQC from their own code with warnings switched on.
warnings.warn(
    'genbenchQC has been renamed to genomic-benchmarks-qc. This is the final '
    'release under the old name; install the new package with '
    '`pip install genomic-benchmarks-qc` and use `gb-qc` in place of the '
    '`genbenchQC` command. '
    'https://genomic-benchmarks.github.io/genomic-benchmarks-qc/',
    DeprecationWarning,
)
