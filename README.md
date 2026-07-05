# Oxagorn

A high-performance Rust reimplementation of [ARAGORN](http://www.ansikte.se/ARAGORN/) for tRNA and tmRNA gene detection.

## Why Oxagorn?

[ARAGORN](http://www.ansikte.se/ARAGORN/) (Laslett & Canback, 2004) is a widely-used and well-established tool for detecting tRNA and tmRNA genes in nucleotide sequences. It has become an essential dependency for bacterial genome annotation pipelines such as [Bakta](https://github.com/oschwengers/bakta), [Prokka](https://github.com/tseemann/prokka), and others.

However, ARAGORN was written in C over 20 years ago and has some limitations in modern computing environments:

- **Single-threaded only** - Cannot utilize multiple CPU cores
- **No Python bindings** - Requires subprocess calls or external wrappers
- **Legacy codebase** - Difficult to maintain and extend

**Oxagorn** addresses these limitations:

- **Multi-core support** - Parallel processing for 5-6x faster analysis
- **Native Python bindings** - Drop-in replacement for [PyARAGORN](https://github.com/althonos/pyaragorn)
- **Modern Rust implementation** - Memory-safe, maintainable codebase
- **100% compatible** - Identical output to ARAGORN v1.2.41

## Performance

Benchmarked on *Arabidopsis thaliana* genome (120 Mb):

| Tool | Time | Speedup |
|------|------|---------|
| ARAGORN v1.2.41 (C) | 26.08s | 1.0x |
| Oxagorn (single-threaded) | 24.21s | 1.08x |
| Oxagorn (6 threads) | **4.58s** | **5.7x** |

## Installation

### Binary (Rust)

```bash
# From source
git clone https://github.com/necoli1822/oxagorn.git
cd oxagorn
cargo build --release

# The binary will be at target/release/oxagorn
```

### Python Package

```bash
# From PyPI (coming soon)
pip install oxagorn

# From source
pip install -e .
```

**Note:** The Python package requires the `oxagorn` binary. Set `OXAGORN_BIN` environment variable or ensure it's in your PATH.

## Usage

### Command Line (ARAGORN compatible)

```bash
# Same options as ARAGORN
oxagorn -t -m -gc11 -w genome.fasta

# With parallel processing (new feature)
oxagorn -@4 genome.fasta  # Use 4 threads
```

### Python (PyARAGORN compatible)

```python
import oxagorn

# Create finder with bacterial genetic code
finder = oxagorn.RNAFinder(translation_table=11)

# Find RNA genes
genes = finder.find_rna(sequence)

for gene in genes:
    if gene.type == "tRNA":
        print(f"tRNA-{gene.amino_acid}({gene.anticodon}) at {gene.begin}-{gene.end}")
    elif gene.type == "tmRNA":
        print(f"tmRNA at {gene.begin}-{gene.end}, tag peptide: {gene.peptide()}")
```

## Bakta Integration

Oxagorn can be used as a drop-in replacement for ARAGORN in Bakta:

```bash
# Option 1: Symlink
ln -sf /path/to/oxagorn ~/.local/bin/aragorn

# Option 2: PATH priority
export PATH="/path/to/oxagorn:$PATH"
```

## API Reference

### Python Classes

**RNAFinder:**
```python
RNAFinder(
    translation_table: int = 1,  # NCBI genetic code (1-25)
    trna: bool = True,           # Search for tRNA genes
    tmrna: bool = True,          # Search for tmRNA genes
    linear: bool = False,        # Linear topology (default: circular)
)
```

**TRNAGene:** `type`, `begin`, `end`, `strand`, `length`, `amino_acid`, `anticodon`, `sequence()`

**TMRNAGene:** `type`, `begin`, `end`, `strand`, `length`, `permuted`, `peptide()`, `sequence()`

## License

GPL-3.0 (same as ARAGORN)

## Citation

If you use Oxagorn, please cite the original ARAGORN paper:

> Laslett, D. and Canback, B. (2004) ARAGORN, a program for the detection of transfer RNA and transfer-messenger RNA genes in nucleotide sequences. *Nucleic Acids Research*, 32:11-16.

## Authors

- **Dean Laslett** - Original ARAGORN algorithm and implementation
- **Sunju Kim** ([@necoli1822](https://github.com/necoli1822)) - Rust reimplementation
