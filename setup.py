"""
Setup script for oxagorn - Pure Python package.

This package requires the oxagorn binary to be installed separately
or available in PATH.

For the bundled binary version, use: maturin build
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_path = Path(__file__).parent / "README.md"
long_description = ""
if readme_path.exists():
    long_description = readme_path.read_text(encoding="utf-8")

setup(
    name="oxagorn",
    version="0.1.0",
    author="Sunju Kim",
    author_email="n.e.coli.1822@gmail.com",
    description="Python bindings for Oxagorn - fast tRNA/tmRNA gene finder (PyARAGORN compatible)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/necoli1822/oxagorn",
    license="GPL-3.0",
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=["bioinformatics", "tRNA", "tmRNA", "gene-finding", "aragorn"],
    package_data={
        "oxagorn": ["py.typed"],
    },
    include_package_data=True,
)
