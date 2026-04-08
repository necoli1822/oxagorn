"""
Oxagorn - Core library implementing PyARAGORN-compatible API

This module provides Gene, TRNAGene, TMRNAGene, and RNAFinder classes
that are API-compatible with PyARAGORN.
"""

import os
import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import List, Literal, Optional, Tuple, Union
from pathlib import Path


# Find oxagorn binary
def _find_oxagorn_binary() -> str:
    """Find the oxagorn binary in common locations."""
    # Check environment variable first
    env_path = os.environ.get("OXAGORN_BIN")
    if env_path and os.path.isfile(env_path):
        return env_path

    # Check common locations
    candidates = [
        # Relative to this file (when installed with package)
        Path(__file__).parent.parent.parent / "target" / "release" / "oxagorn",
        # System PATH
        "oxagorn",
        # Development location
        Path(__file__).parent.parent.parent / "target" / "debug" / "oxagorn",
    ]

    for candidate in candidates:
        if isinstance(candidate, Path):
            if candidate.is_file():
                return str(candidate)
        else:
            # Check if in PATH
            try:
                result = subprocess.run(
                    ["which", candidate],
                    capture_output=True,
                    text=True
                )
                if result.returncode == 0:
                    return result.stdout.strip()
            except Exception:
                pass

    raise RuntimeError(
        "Could not find oxagorn binary. Set OXAGORN_BIN environment variable "
        "or ensure oxagorn is in PATH."
    )


@dataclass
class Gene:
    """Base class for RNA genes detected by Oxagorn.

    Attributes:
        type: Gene type ("tRNA" or "tmRNA")
        begin: Start position (1-based)
        end: End position (1-based)
        strand: Strand direction (1 for forward, -1 for reverse)
        energy: Detection energy score
        raw_energy: Raw energy before normalization
        _sequence: Internal sequence storage
    """
    _type: str
    begin: int
    end: int
    strand: Literal[1, -1]
    energy: float
    raw_energy: float
    _sequence: str = ""

    @property
    def type(self) -> str:
        return self._type

    @property
    def length(self) -> int:
        return abs(self.end - self.begin) + 1

    def sequence(self) -> str:
        """Return the gene sequence."""
        return self._sequence

    def __repr__(self) -> str:
        return f"<Gene type='{self._type}' begin={self.begin} end={self.end} strand={self.strand}>"


@dataclass
class TRNAGene(Gene):
    """tRNA gene with anticodon information.

    Attributes:
        amino_acid: Three-letter amino acid code (e.g., "Ala", "Gly")
        anticodon: Anticodon sequence (e.g., "tgc")
        anticodon_offset: Position of anticodon within the gene
        anticodon_length: Length of anticodon (usually 3)
    """
    amino_acid: str = ""
    _amino_acids: Tuple[str, ...] = field(default_factory=tuple)
    anticodon: str = ""
    anticodon_offset: int = 0
    anticodon_length: int = 3

    @property
    def type(self) -> str:
        return "tRNA"

    @property
    def amino_acids(self) -> Tuple[str, ...]:
        if self._amino_acids:
            return self._amino_acids
        return (self.amino_acid,)

    def __repr__(self) -> str:
        return (f"<TRNAGene amino_acid='{self.amino_acid}' anticodon='{self.anticodon}' "
                f"begin={self.begin} end={self.end} strand={self.strand}>")


@dataclass
class TMRNAGene(Gene):
    """tmRNA gene with ORF information.

    Attributes:
        permuted: Whether the tmRNA has a permuted structure
        orf_offset: Position of the ORF within the gene
        orf_length: Length of the ORF
        _orf_sequence: ORF nucleotide sequence
        _peptide_sequence: Translated peptide sequence
    """
    permuted: bool = False
    orf_offset: int = 0
    orf_length: int = 0
    _orf_sequence: str = ""
    _peptide_sequence: str = ""

    @property
    def type(self) -> str:
        return "tmRNA"

    def orf(self, include_stop: bool = True) -> str:
        """Return the ORF sequence.

        Args:
            include_stop: Whether to include the stop codon
        """
        if include_stop:
            return self._orf_sequence
        # Remove last 3 bases (stop codon)
        if len(self._orf_sequence) >= 3:
            return self._orf_sequence[:-3]
        return self._orf_sequence

    def peptide(self, include_stop: bool = True) -> str:
        """Return the translated peptide sequence.

        Args:
            include_stop: Whether to include the stop codon marker (*)
        """
        if include_stop:
            return self._peptide_sequence
        return self._peptide_sequence.rstrip("*")

    def __repr__(self) -> str:
        return (f"<TMRNAGene permuted={self.permuted} orf_length={self.orf_length} "
                f"begin={self.begin} end={self.end} strand={self.strand}>")


class RNAFinder:
    """RNA gene finder using Oxagorn (ARAGORN reimplementation).

    A drop-in replacement for PyARAGORN's RNAFinder class.

    Args:
        translation_table: NCBI genetic code table (1-25, default: 1 = Standard)
        trna: Whether to search for tRNA genes (default: True)
        tmrna: Whether to search for tmRNA genes (default: True)
        linear: Whether sequence is linear (default: False = circular)
        threshold_scale: Scale factor for detection thresholds (default: 1.0)

    Example:
        >>> finder = RNAFinder(translation_table=11)
        >>> genes = finder.find_rna(sequence)
        >>> for gene in genes:
        ...     if gene.type == "tRNA":
        ...         print(gene.amino_acid, gene.anticodon)
    """

    def __init__(
        self,
        translation_table: int = 1,
        trna: bool = True,
        tmrna: bool = True,
        linear: bool = False,
        threshold_scale: float = 1.0,
    ):
        if not 1 <= translation_table <= 25:
            raise ValueError(f"Invalid translation table: {translation_table}. Must be 1-25.")

        self._translation_table = translation_table
        self._trna = trna
        self._tmrna = tmrna
        self._linear = linear
        self._threshold_scale = threshold_scale

        # Find binary
        self._binary = _find_oxagorn_binary()

    @property
    def translation_table(self) -> int:
        return self._translation_table

    @property
    def trna(self) -> bool:
        return self._trna

    @property
    def tmrna(self) -> bool:
        return self._tmrna

    @property
    def linear(self) -> bool:
        return self._linear

    @property
    def threshold_scale(self) -> float:
        return self._threshold_scale

    @threshold_scale.setter
    def threshold_scale(self, value: float):
        self._threshold_scale = value

    def find_rna(
        self,
        sequence: Union[str, bytes, bytearray]
    ) -> List[Union[TRNAGene, TMRNAGene]]:
        """Find RNA genes in a sequence.

        Args:
            sequence: DNA sequence as str, bytes, or bytearray

        Returns:
            List of Gene objects (TRNAGene or TMRNAGene)
        """
        # Convert to bytes if needed
        if isinstance(sequence, str):
            seq_bytes = sequence.encode('ascii')
        elif isinstance(sequence, bytearray):
            seq_bytes = bytes(sequence)
        else:
            seq_bytes = sequence

        # Run oxagorn and parse output
        return self._run_oxagorn(seq_bytes)

    def _run_oxagorn(self, sequence: bytes) -> List[Union[TRNAGene, TMRNAGene]]:
        """Run oxagorn binary and parse output."""
        # Create temporary file with sequence
        with tempfile.NamedTemporaryFile(mode='wb', suffix='.fasta', delete=False) as f:
            f.write(b">seq\n")
            f.write(sequence)
            f.write(b"\n")
            temp_path = f.name

        try:
            # Build command
            cmd = [self._binary]

            # Add options
            if self._trna and not self._tmrna:
                cmd.append("-t")  # tRNA only
            elif self._tmrna and not self._trna:
                cmd.append("-m")  # tmRNA only
            # Default: both tRNA and tmRNA

            if self._linear:
                cmd.append("-l")  # Linear sequence

            # Genetic code (format: -gc11, not -gc 11)
            if self._translation_table != 1:
                cmd.append(f"-gc{self._translation_table}")

            # Add batch output for parsing
            cmd.append("-w")  # Machine-readable output
            cmd.append("-seq")  # Include sequences in output

            # Input file
            cmd.append(temp_path)

            # Run oxagorn
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )

            if result.returncode != 0:
                raise RuntimeError(f"oxagorn failed: {result.stderr}")

            # Parse output
            return self._parse_output(result.stdout)

        finally:
            # Clean up temp file
            try:
                os.unlink(temp_path)
            except Exception:
                pass

    def _parse_output(self, output: str) -> List[Union[TRNAGene, TMRNAGene]]:
        """Parse oxagorn output into Gene objects.

        Output format with -seq flag:
        >seq_name
        N genes found
        1   tRNA-Ile  [start,end]  score  (anticodon)
        sequence_line1
        sequence_line2
        2   tRNA-Ala  [start,end]  score  (anticodon)
        ...
        """
        genes = []
        current_gene = None
        sequence_lines = []

        for line in output.strip().split('\n'):
            if not line.strip():
                continue

            # Skip header lines
            if line.startswith('>') or line.startswith('#'):
                continue

            # Skip "N genes found" line
            if 'genes found' in line:
                continue

            # Check if this is a gene info line or sequence line
            # Gene info lines start with a number followed by spaces and tRNA/tmRNA
            gene_match = re.match(r'\s*\d+\s+(tRNA|tmRNA)', line)

            if gene_match:
                # Save previous gene with its sequence
                if current_gene is not None:
                    current_gene._sequence = ''.join(sequence_lines)
                    genes.append(current_gene)
                    sequence_lines = []

                # Parse new gene
                current_gene = self._parse_gene_line(line)
            elif current_gene is not None:
                # This is a sequence line
                # Only add if it looks like nucleotide sequence (lowercase letters)
                clean_line = line.strip()
                if clean_line and re.match(r'^[acgtuACGTU]+$', clean_line):
                    sequence_lines.append(clean_line)

        # Don't forget the last gene
        if current_gene is not None:
            current_gene._sequence = ''.join(sequence_lines)
            genes.append(current_gene)

        return genes

    def _parse_gene_line(self, line: str) -> Optional[Union[TRNAGene, TMRNAGene]]:
        """Parse a single gene line from oxagorn output.

        Format examples:
        1   tRNA-Ile               [225381,225457]	35  	(gat)
        8   tRNA-Gln              c[696430,696504]	33  	(ctg)
        89  tmRNA                 c[2753620,2754007]	38  	Tag:ANDENYALAA*
        """
        # Match tRNA pattern: number, tRNA-AA, optional 'c', [start,end], tab, score, tab, (anticodon)
        trna_match = re.match(
            r'\s*\d+\s+tRNA-(\w+)\s+(c?)\[(\d+),(\d+)\]\s+(\d+)\s+\((\w+)\)',
            line
        )
        if trna_match:
            aa = trna_match.group(1)
            comp = trna_match.group(2) == 'c'
            start = int(trna_match.group(3))
            end = int(trna_match.group(4))
            energy = float(trna_match.group(5))
            anticodon = trna_match.group(6)

            return TRNAGene(
                _type="tRNA",
                begin=start,
                end=end,
                strand=-1 if comp else 1,
                energy=energy,
                raw_energy=energy,
                amino_acid=aa,
                anticodon=anticodon,
                anticodon_length=len(anticodon),
            )

        # Match tmRNA pattern: 46  tmRNA                [2755593,2755955]	90,125	ANDENYALAA**
        tmrna_match = re.match(
            r'\s*\d+\s+tmRNA\s+(c?)\[(\d+),(\d+)\]\s+(\d+),?(\d*)\s+(\S+)',
            line
        )
        if tmrna_match:
            comp = tmrna_match.group(1) == 'c'
            start = int(tmrna_match.group(2))
            end = int(tmrna_match.group(3))
            energy1 = float(tmrna_match.group(4))
            energy2 = float(tmrna_match.group(5)) if tmrna_match.group(5) else energy1
            tag = tmrna_match.group(6).rstrip('*')  # Remove trailing asterisks

            return TMRNAGene(
                _type="tmRNA",
                begin=start,
                end=end,
                strand=-1 if comp else 1,
                energy=energy1,
                raw_energy=energy2,
                _peptide_sequence=tag + "*",  # Add single stop marker
                orf_length=len(tag) * 3 if tag else 0,
            )

        return None

    def __repr__(self) -> str:
        return (f"RNAFinder(translation_table={self._translation_table}, "
                f"trna={self._trna}, tmrna={self._tmrna}, "
                f"linear={self._linear}, threshold_scale={self._threshold_scale})")
