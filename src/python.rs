/*
 * Oxagorn v0.1.0 - Python bindings
 * python.rs - PyO3 bindings for Python compatibility with PyARAGORN
 *
 * This module provides a drop-in replacement for PyARAGORN
 * API: RNAFinder, Gene, TRNAGene, TMRNAGene
 */

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use std::sync::Arc;

use crate::types::{Gene as RustGene, Opts, GeneCode, NGENECODE};
use crate::tables::AAMAP;
use crate::output::parallel_fastafile_bytes;

// =============================================================================
// Gene base class - common properties for tRNA and tmRNA
// =============================================================================

#[pyclass(subclass)]
#[derive(Clone)]
pub struct Gene {
    pub(crate) gene_type: String,
    pub(crate) begin: i64,
    pub(crate) end: i64,
    pub(crate) strand: i32,  // 1 or -1
    pub(crate) energy: f64,
    pub(crate) raw_energy: f64,
    pub(crate) sequence_data: String,
}

#[pymethods]
impl Gene {
    #[getter]
    fn r#type(&self) -> &str {
        &self.gene_type
    }

    #[getter]
    fn begin(&self) -> i64 {
        self.begin
    }

    #[getter]
    fn end(&self) -> i64 {
        self.end
    }

    #[getter]
    fn length(&self) -> i64 {
        (self.end - self.begin).abs() + 1
    }

    #[getter]
    fn strand(&self) -> i32 {
        self.strand
    }

    #[getter]
    fn energy(&self) -> f64 {
        self.energy
    }

    #[getter]
    fn raw_energy(&self) -> f64 {
        self.raw_energy
    }

    fn sequence(&self) -> &str {
        &self.sequence_data
    }

    fn __repr__(&self) -> String {
        format!(
            "<Gene type='{}' begin={} end={} strand={}>",
            self.gene_type, self.begin, self.end, self.strand
        )
    }
}

// =============================================================================
// TRNAGene - tRNA-specific gene with anticodon information
// =============================================================================

#[pyclass(extends=Gene)]
#[derive(Clone)]
pub struct TRNAGene {
    amino_acid: String,
    amino_acids: Vec<String>,
    anticodon: String,
    anticodon_offset: i32,
    anticodon_length: i32,
}

#[pymethods]
impl TRNAGene {
    #[getter]
    fn r#type(&self) -> &str {
        "tRNA"
    }

    #[getter]
    fn amino_acid(&self) -> &str {
        &self.amino_acid
    }

    #[getter]
    fn amino_acids(&self) -> Vec<String> {
        self.amino_acids.clone()
    }

    #[getter]
    fn anticodon(&self) -> &str {
        &self.anticodon
    }

    #[getter]
    fn anticodon_offset(&self) -> i32 {
        self.anticodon_offset
    }

    #[getter]
    fn anticodon_length(&self) -> i32 {
        self.anticodon_length
    }

    fn __repr__(self_: PyRef<'_, Self>) -> String {
        let base = self_.as_ref();
        format!(
            "<TRNAGene amino_acid='{}' anticodon='{}' begin={} end={} strand={}>",
            self_.amino_acid, self_.anticodon, base.begin, base.end, base.strand
        )
    }
}

// =============================================================================
// TMRNAGene - tmRNA-specific gene with ORF information
// =============================================================================

#[pyclass(extends=Gene)]
#[derive(Clone)]
pub struct TMRNAGene {
    permuted: bool,
    orf_offset: i32,
    orf_length: i32,
    orf_sequence: String,
    peptide_sequence: String,
}

#[pymethods]
impl TMRNAGene {
    #[getter]
    fn r#type(&self) -> &str {
        "tmRNA"
    }

    #[getter]
    fn permuted(&self) -> bool {
        self.permuted
    }

    #[getter]
    fn orf_offset(&self) -> i32 {
        self.orf_offset
    }

    #[getter]
    fn orf_length(&self) -> i32 {
        self.orf_length
    }

    #[pyo3(signature = (include_stop = true))]
    fn orf(&self, include_stop: bool) -> String {
        if include_stop {
            self.orf_sequence.clone()
        } else {
            // Remove stop codon (last 3 bases)
            let len = self.orf_sequence.len();
            if len >= 3 {
                self.orf_sequence[..len-3].to_string()
            } else {
                self.orf_sequence.clone()
            }
        }
    }

    #[pyo3(signature = (include_stop = true))]
    fn peptide(&self, include_stop: bool) -> String {
        if include_stop {
            self.peptide_sequence.clone()
        } else {
            // Remove stop codon marker
            self.peptide_sequence.trim_end_matches('*').to_string()
        }
    }

    fn __repr__(self_: PyRef<'_, Self>) -> String {
        let base = self_.as_ref();
        format!(
            "<TMRNAGene permuted={} orf_length={} begin={} end={} strand={}>",
            self_.permuted, self_.orf_length, base.begin, base.end, base.strand
        )
    }
}

// =============================================================================
// RNAFinder - main class for finding tRNA/tmRNA genes
// =============================================================================

#[pyclass]
pub struct RNAFinder {
    translation_table: i32,
    trna: bool,
    tmrna: bool,
    linear: bool,
    threshold_scale: f64,
}

#[pymethods]
impl RNAFinder {
    #[new]
    #[pyo3(signature = (translation_table = 1, trna = true, tmrna = true, linear = false, threshold_scale = 1.0))]
    fn new(
        translation_table: i32,
        trna: bool,
        tmrna: bool,
        linear: bool,
        threshold_scale: f64,
    ) -> PyResult<Self> {
        // Validate translation table (1-25 for NCBI genetic codes)
        if translation_table < 1 || translation_table > 25 {
            return Err(PyValueError::new_err(format!(
                "Invalid translation table: {}. Must be 1-25.",
                translation_table
            )));
        }

        Ok(RNAFinder {
            translation_table,
            trna,
            tmrna,
            linear,
            threshold_scale,
        })
    }

    #[getter]
    fn translation_table(&self) -> i32 {
        self.translation_table
    }

    #[getter]
    fn trna(&self) -> bool {
        self.trna
    }

    #[getter]
    fn tmrna(&self) -> bool {
        self.tmrna
    }

    #[getter]
    fn linear(&self) -> bool {
        self.linear
    }

    #[getter]
    fn threshold_scale(&self) -> f64 {
        self.threshold_scale
    }

    #[setter]
    fn set_threshold_scale(&mut self, value: f64) {
        self.threshold_scale = value;
    }

    /// Find RNA genes in a sequence
    ///
    /// Args:
    ///     sequence: DNA sequence as str, bytes, or bytearray
    ///
    /// Returns:
    ///     List of Gene objects (TRNAGene or TMRNAGene)
    fn find_rna(&self, py: Python<'_>, sequence: &Bound<'_, PyAny>) -> PyResult<Vec<PyObject>> {
        // Convert input to bytes
        let seq_bytes: Vec<u8> = if let Ok(s) = sequence.extract::<String>() {
            s.into_bytes()
        } else if let Ok(b) = sequence.extract::<Vec<u8>>() {
            b
        } else {
            return Err(PyValueError::new_err(
                "sequence must be str, bytes, or bytearray"
            ));
        };

        // Call the oxagorn detection
        let results = self.detect_genes(&seq_bytes)?;

        // Convert results to Python objects
        let mut py_results = Vec::new();
        for gene_result in results {
            let py_obj = self.rust_gene_to_python(py, gene_result)?;
            py_results.push(py_obj);
        }

        Ok(py_results)
    }

    fn __repr__(&self) -> String {
        format!(
            "RNAFinder(translation_table={}, trna={}, tmrna={}, linear={}, threshold_scale={})",
            self.translation_table, self.trna, self.tmrna, self.linear, self.threshold_scale
        )
    }
}

// Internal implementation for RNAFinder
impl RNAFinder {
    fn detect_genes(&self, sequence: &[u8]) -> PyResult<Vec<GeneResult>> {
        use crate::types::*;
        use crate::output::detect_genes_in_memory;

        // Set up options matching RNAFinder parameters
        let mut sw = Opts::default();
        sw.trna = if self.trna { 1 } else { 0 };
        sw.tmrna = if self.tmrna { 1 } else { 0 };
        sw.linear = if self.linear { 1 } else { 0 };
        sw.mtefold = if self.linear { 0 } else { 1 };  // Circular by default

        // Set genetic code
        sw.geneticcode = match self.translation_table {
            1 => STANDARD,
            2 => VERTEBRATE_MT,
            3..=25 => self.translation_table,
            _ => STANDARD,
        };

        // Set threshold based on scale
        // Default thresholds from ARAGORN
        sw.trnathresh = (95.0 * self.threshold_scale) as i32;
        sw.tmrnathresh = (90.0 * self.threshold_scale) as i32;

        // Detect genes
        let results = detect_genes_in_memory(sequence, &sw);

        Ok(results)
    }

    fn rust_gene_to_python(&self, py: Python<'_>, gene: GeneResult) -> PyResult<PyObject> {
        match gene.gene_type.as_str() {
            "tRNA" => {
                let base = Gene {
                    gene_type: "tRNA".to_string(),
                    begin: gene.begin,
                    end: gene.end,
                    strand: gene.strand,
                    energy: gene.energy,
                    raw_energy: gene.raw_energy,
                    sequence_data: gene.sequence,
                };
                let trna = TRNAGene {
                    amino_acid: gene.amino_acid,
                    amino_acids: gene.amino_acids,
                    anticodon: gene.anticodon,
                    anticodon_offset: gene.anticodon_offset,
                    anticodon_length: gene.anticodon_length,
                };
                Ok(Py::new(py, (trna, base))?.into_py(py))
            }
            "tmRNA" => {
                let base = Gene {
                    gene_type: "tmRNA".to_string(),
                    begin: gene.begin,
                    end: gene.end,
                    strand: gene.strand,
                    energy: gene.energy,
                    raw_energy: gene.raw_energy,
                    sequence_data: gene.sequence,
                };
                let tmrna = TMRNAGene {
                    permuted: gene.permuted,
                    orf_offset: gene.orf_offset,
                    orf_length: gene.orf_length,
                    orf_sequence: gene.orf_sequence,
                    peptide_sequence: gene.peptide_sequence,
                };
                Ok(Py::new(py, (tmrna, base))?.into_py(py))
            }
            _ => {
                // Generic gene
                let base = Gene {
                    gene_type: gene.gene_type,
                    begin: gene.begin,
                    end: gene.end,
                    strand: gene.strand,
                    energy: gene.energy,
                    raw_energy: gene.raw_energy,
                    sequence_data: gene.sequence,
                };
                Ok(Py::new(py, base)?.into_py(py))
            }
        }
    }
}

// =============================================================================
// GeneResult - Internal struct for passing gene data
// =============================================================================

#[derive(Clone, Default)]
pub struct GeneResult {
    pub gene_type: String,
    pub begin: i64,
    pub end: i64,
    pub strand: i32,
    pub energy: f64,
    pub raw_energy: f64,
    pub sequence: String,
    // tRNA specific
    pub amino_acid: String,
    pub amino_acids: Vec<String>,
    pub anticodon: String,
    pub anticodon_offset: i32,
    pub anticodon_length: i32,
    // tmRNA specific
    pub permuted: bool,
    pub orf_offset: i32,
    pub orf_length: i32,
    pub orf_sequence: String,
    pub peptide_sequence: String,
}

// =============================================================================
// Python module definition
// =============================================================================

#[pymodule]
#[pyo3(name = "pyoxagorn")]
fn pyoxagorn(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Gene>()?;
    m.add_class::<TRNAGene>()?;
    m.add_class::<TMRNAGene>()?;
    m.add_class::<RNAFinder>()?;

    // Add version info
    m.add("__version__", "0.1.0")?;
    m.add("__author__", "Sunju Kim")?;

    Ok(())
}
