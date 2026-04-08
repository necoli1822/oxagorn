"""
Oxagorn - Python bindings for Oxagorn (ARAGORN reimplementation)

A drop-in replacement for PyARAGORN with the same API.
Uses the high-performance Rust implementation for 5.6x faster detection.

Example:
    >>> import oxagorn
    >>> finder = oxagorn.RNAFinder(translation_table=11)
    >>> genes = finder.find_rna(b"ATGCATGC...")
    >>> for gene in genes:
    ...     print(gene.type, gene.begin, gene.end)
"""

from .lib import Gene, TRNAGene, TMRNAGene, RNAFinder

__version__ = "0.1.0"
__author__ = "Sunju Kim"
__all__ = ["Gene", "TRNAGene", "TMRNAGene", "RNAFinder"]
