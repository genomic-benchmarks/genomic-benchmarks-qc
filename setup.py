from setuptools import find_packages, setup

requirements = [
    'numpy>=1.23',
    'pandas>=1.5',
    'matplotlib>=3.6',
    'seaborn>=0.12',
    'biopython>=1.8',
    'scikit-learn>=1.2',
    'cdhit-reader==0.2.0',
    'statsmodels>=0.13',
]

test_requirements = [
    'pytest>=3',
]

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='genbenchQC',
    version='1.0.2',
    description='Genomic Benchmarks QC: Automated Quality Control for Genomic Machine Learning Datasets',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Katarina Gresova",
    author_email='gresova11@gmail.com',
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
      evaluate_sequences=genbenchQC.evaluate_sequences:main
      evaluate_dataset=genbenchQC.evaluate_dataset:main
      evaluate_split=genbenchQC.evaluate_split:main
      ''',
    keywords=["genomic benchmarks", "deep learning", "machine learning",
      "computational biology", "bioinformatics", "genomics", "quality control"],
    classifiers=[
        # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        "Development Status :: 3 - Alpha",
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