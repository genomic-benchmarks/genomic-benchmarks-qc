from setuptools import find_packages, setup

requirements = [
    'numpy>=1.23',
    'pandas>=1.5',
    'matplotlib>=3.6',
    'seaborn>=0.12',
    'biopython>=1.8',
    'scikit-learn>=1.2',
    'statsmodels>=0.13',
    'typer>=0.20',
    # The rename, expressed as a dependency: installing this package on a
    # Python that can run the new one brings the new one with it.
    #
    # The marker is what makes that safe. genomic-benchmarks-qc requires Python
    # 3.12 and this package is used on 3.8 upwards, so an unconditional
    # dependency would be unsatisfiable on the older interpreters - and pip does
    # not stop at unsatisfiable, it backtracks and quietly installs 1.1.0
    # instead, which carries none of this notice. Marked, 1.2.0 stays
    # installable everywhere and still pulls the new tool in wherever it can
    # run. There is deliberately no python_requires= below for the same reason.
    'genomic-benchmarks-qc>=1.0.0; python_version >= "3.12"',
]

test_requirements = [
    'pytest>=3',
]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='genbenchQC',
    version='1.2.0',
    description='DEPRECATED - renamed to genomic-benchmarks-qc. Automated Quality Control for Genomic Machine Learning Datasets',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Katarina Gresova",
    author_email='gresova11@gmail.com',
    # This package had no URLs at all in its metadata, so its PyPI page was the
    # rendered README and nothing else. The point of a final release is to be a
    # signpost, and a signpost wants links in the sidebar as well as in the
    # prose.
    url='https://github.com/genomic-benchmarks/genomic-benchmarks-qc',
    project_urls={
        'New package': 'https://pypi.org/project/genomic-benchmarks-qc/',
        'Documentation': 'https://genomic-benchmarks.github.io/genomic-benchmarks-qc/',
        'Source': 'https://github.com/genomic-benchmarks/genomic-benchmarks-qc',
    },
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=requirements,
    extras_require={
        "develop": test_requirements,
    },
    tests_require=["pytest"],
    test_suite='tests',
    entry_points='''
      [console_scripts]
      genbenchQC=genbenchQC.cli:main
      ''',
    keywords=["genomic benchmarks", "deep learning", "machine learning",
      "computational biology", "bioinformatics", "genomics", "quality control"],
    classifiers=[
        # Inactive: this package is finished, and the classifier is one of the
        # few parts of the metadata a person reads without reading the README.
        "Development Status :: 7 - Inactive",
        # Define that your audience are developers
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: MIT License",
        # Specify which pyhton versions that you want to support
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
)
