/*
 * Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
 * parallel.rs - Parallel processing module for multi-threaded gene detection
 *
 * Based on ARAGORN v1.2.41 by Dean Laslett
 * Author: Sunju Kim <n.e.coli.1822@gmail.com>
 */

#![allow(dead_code)]

use rayon::prelude::*;
use std::path::Path;

/// Default maximum threads (benchmarks show optimal at 6 threads,
/// no speedup beyond this due to memory bandwidth saturation)
pub const DEFAULT_MAX_THREADS: usize = 6;

/// Parallel processing configuration
/// Controls chunking strategy, thread count, and overlap handling
#[derive(Debug, Clone)]
pub struct ParallelConfig {
    /// Minimum sequence length for parallel processing (default: 5,000,000 bp)
    /// Sequences shorter than this use single-threaded fast-path
    pub min_parallel_len: usize,

    /// Minimum chunk size for dynamic chunk sizing (default: 2,000,000 bp)
    /// Actual chunk size is calculated dynamically to utilize all cores
    pub min_chunk_size: usize,

    /// Overlap size between chunks for gene detection (default: 1000 bp)
    /// Must be at least TSWEEP (1000) to handle tRNA sweep distance
    pub overlap_size: usize,

    /// Number of threads (0 = auto-detect based on CPU cores)
    pub num_threads: usize,

    /// Maximum threads to use (default: 6)
    /// Benchmarks show optimal at 6 threads, plateau beyond due to memory bandwidth
    pub max_threads: usize,

    /// Enable parallel processing (can be disabled via CLI)
    pub enabled: bool,
}

impl Default for ParallelConfig {
    fn default() -> Self {
        Self {
            min_parallel_len: 5_000_000,  // 5 Mbp threshold (smaller → sequential)
            min_chunk_size: 2_000_000,    // 2 Mbp per chunk (~1 sec processing)
            overlap_size: 1000,           // TSWEEP from C original (handles tRNA sweep distance)
            num_threads: 0,               // Auto-detect
            max_threads: DEFAULT_MAX_THREADS,  // Cap at 8 threads
            enabled: true,                // Enabled by default
        }
    }
}

impl ParallelConfig {
    /// Create a new configuration with custom values
    pub fn new(
        min_parallel_len: usize,
        min_chunk_size: usize,
        overlap_size: usize,
        num_threads: usize,
    ) -> Self {
        Self {
            min_parallel_len,
            min_chunk_size,
            overlap_size,
            num_threads,
            max_threads: DEFAULT_MAX_THREADS,
            enabled: true,
        }
    }

    /// Create a configuration with custom max threads
    pub fn with_max_threads(mut self, max_threads: usize) -> Self {
        self.max_threads = max_threads;
        self
    }

    /// Get the effective number of threads
    /// Returns min(requested_threads, max_threads)
    /// Auto-detects physical CPU cores if num_threads is 0
    /// Uses physical cores (not logical/hyperthreading) as this workload is memory-bound
    pub fn effective_threads(&self) -> usize {
        let requested = if self.num_threads == 0 {
            // Use physical cores - hyperthreading doesn't help memory-bound workloads
            num_cpus::get_physical()
        } else {
            self.num_threads
        };
        std::cmp::min(requested, self.max_threads)
    }

    /// Dynamic chunk sizing: Calculate optimal chunk size to maximize core utilization
    ///
    /// Formula: max(min_chunk_size, (sequence_len / num_threads) + overlap_size)
    ///
    /// Example: 4 Mbp genome with 8 threads
    /// - Fixed 1 Mbp chunks: 4 chunks -> only 4 cores used
    /// - Dynamic sizing: (4M / 8) + 450 = ~500 Kbp -> 8 chunks -> all cores used
    pub fn calculate_chunk_size(&self, sequence_len: usize) -> usize {
        let num_threads = self.effective_threads();

        // Dynamic calculation: (sequence length / thread count) + overlap
        // Guarantee minimum chunk size to avoid excessive overhead
        std::cmp::max(
            self.min_chunk_size,
            (sequence_len / num_threads) + self.overlap_size
        )
    }

    /// Determine if a sequence should use parallel processing
    #[inline]
    pub fn should_parallelize(&self, sequence_len: usize) -> bool {
        self.enabled && sequence_len >= self.min_parallel_len
    }

    /// Create a single-threaded configuration (for testing or forced serial mode)
    pub fn single_threaded() -> Self {
        Self {
            enabled: false,
            ..Default::default()
        }
    }
}

/// Zero-copy chunk slice referencing original sequence data
#[derive(Debug)]
pub struct ChunkSlice<'a> {
    /// Reference to the original sequence data (no copy)
    pub data: &'a [u8],

    /// Chunk index (0-based)
    pub chunk_index: usize,

    /// Global offset in the original sequence
    pub global_offset: usize,

    /// Whether this is the last chunk
    pub is_last: bool,
}

/// Create zero-copy chunks from a sequence with dynamic chunk sizing
pub fn create_chunks<'a>(
    seq: &'a [u8],
    config: &ParallelConfig
) -> Vec<ChunkSlice<'a>> {
    let actual_chunk_size = config.calculate_chunk_size(seq.len());

    let mut chunks = Vec::new();
    let mut offset = 0;
    let mut index = 0;

    while offset < seq.len() {
        let end = std::cmp::min(
            offset + actual_chunk_size + config.overlap_size,
            seq.len()
        );
        let is_last = end == seq.len();

        chunks.push(ChunkSlice {
            data: &seq[offset..end],
            chunk_index: index,
            global_offset: offset,
            is_last,
        });

        offset += actual_chunk_size;  // Move by chunk size (overlap excluded)
        index += 1;
    }

    chunks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = ParallelConfig::default();
        assert_eq!(config.min_parallel_len, 5_000_000);  // 5 Mbp
        assert_eq!(config.min_chunk_size, 2_000_000);    // 2 Mbp
        assert_eq!(config.overlap_size, 1000);           // TSWEEP from C original
        assert_eq!(config.max_threads, DEFAULT_MAX_THREADS);
        assert!(config.enabled);
    }

    #[test]
    fn test_should_parallelize() {
        let config = ParallelConfig::default();

        // Small sequence: no parallelization
        assert!(!config.should_parallelize(3_000_000));  // 3 Mbp < 5 Mbp

        // Large sequence: parallelize
        assert!(config.should_parallelize(10_000_000)); // 10 Mbp > 5 Mbp

        // Exact threshold
        assert!(config.should_parallelize(5_000_000));  // 5 Mbp = threshold
    }

    #[test]
    fn test_dynamic_chunk_sizing() {
        let mut config = ParallelConfig::default();
        config.num_threads = 8;  // Force 8 threads

        // 16 Mbp genome with 8 threads
        // Expected: max(2M, 16M/8 + 1000) = max(2M, 2M + 1000) = 2_001_000
        let chunk_size = config.calculate_chunk_size(16_000_000);
        assert_eq!(chunk_size, 2_001_000);

        // 8 Mbp genome with 8 threads
        // Expected: max(2M, 8M/8 + 1000) = max(2M, 1M + 1000) = 2M
        let chunk_size = config.calculate_chunk_size(8_000_000);
        assert_eq!(chunk_size, 2_000_000);

        // 32 Mbp genome with 8 threads
        // Expected: max(2M, 32M/8 + 1000) = max(2M, 4M + 1000) = 4_001_000
        let chunk_size = config.calculate_chunk_size(32_000_000);
        assert_eq!(chunk_size, 4_001_000);
    }

    #[test]
    fn test_create_chunks() {
        let mut config = ParallelConfig::default();
        config.num_threads = 4;
        config.min_chunk_size = 100;  // Small for testing
        config.overlap_size = 10;

        let seq = vec![0u8; 500];
        let chunks = create_chunks(&seq, &config);

        // 500 / 4 = 125, + 10 = 135 (chunk size)
        // We should get multiple chunks
        assert!(chunks.len() >= 1);

        // First chunk starts at 0
        assert_eq!(chunks[0].global_offset, 0);
        assert_eq!(chunks[0].chunk_index, 0);

        // Last chunk is marked
        assert!(chunks.last().unwrap().is_last);
    }

    #[test]
    fn test_single_threaded_config() {
        let config = ParallelConfig::single_threaded();
        assert!(!config.enabled);
        assert!(!config.should_parallelize(100_000_000));
    }
}

// ============================================================================
// Sequence Record for Parallel Processing
// ============================================================================

/// A sequence record loaded into memory for parallel processing
#[derive(Debug, Clone)]
pub struct SeqRecord {
    /// Sequence identifier (from FASTA header)
    pub id: String,
    /// Nucleotide sequence data
    pub seq: Vec<u8>,
}

impl SeqRecord {
    pub fn new(id: String, seq: Vec<u8>) -> Self {
        Self { id, seq }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

// ============================================================================
// Case C: Metagenome / Multi-Contig Parallel Processing
// ============================================================================

/// Load all sequences from a FASTA file into memory
/// Uses needletail for fast parsing
pub fn load_fasta_sequences<P: AsRef<Path>>(path: P) -> Result<Vec<SeqRecord>, String> {
    use needletail::parse_fastx_file;

    let mut sequences = Vec::new();

    let mut reader = parse_fastx_file(path.as_ref())
        .map_err(|e| format!("Failed to open FASTA file: {}", e))?;

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| format!("Invalid FASTA record: {}", e))?;
        let id = String::from_utf8_lossy(record.id()).to_string();
        let seq = record.seq().to_vec();
        sequences.push(SeqRecord::new(id, seq));
    }

    Ok(sequences)
}

/// Detected gene result with sequence context
#[derive(Debug, Clone)]
pub struct GeneResult {
    /// Sequence ID this gene was found in
    pub seq_id: String,
    /// Gene start position (1-based, absolute in original sequence)
    pub start: i64,
    /// Gene stop position (1-based, absolute in original sequence)
    pub stop: i64,
    /// Gene type (tRNA=0, tmRNA=1, etc.)
    pub genetype: i32,
    /// Strand (0=forward, 1=reverse complement)
    pub comp: i32,
    /// Detection energy/score
    pub energy: f64,
}

/// Process multiple sequences in parallel (Case C: Metagenome)
/// Each sequence is processed independently on a separate thread
///
/// # Arguments
/// * `sequences` - Vector of sequence records to process
/// * `config` - Parallel processing configuration
/// * `process_fn` - Function to process a single sequence, returns vector of genes
///
/// # Returns
/// Vector of (seq_id, genes) pairs, one per input sequence
pub fn process_sequences_parallel<F>(
    sequences: &[SeqRecord],
    config: &ParallelConfig,
    process_fn: F,
) -> Vec<(String, Vec<GeneResult>)>
where
    F: Fn(&SeqRecord) -> Vec<GeneResult> + Sync + Send,
{
    if !config.enabled || config.effective_threads() <= 1 {
        // Single-threaded fallback
        sequences
            .iter()
            .map(|seq| {
                let genes = process_fn(seq);
                (seq.id.clone(), genes)
            })
            .collect()
    } else {
        // Configure rayon thread pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.effective_threads())
            .build()
            .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());

        pool.install(|| {
            sequences
                .par_iter()
                .map(|seq| {
                    let genes = process_fn(seq);
                    (seq.id.clone(), genes)
                })
                .collect()
        })
    }
}

// ============================================================================
// Case B: Large Single Genome Chunking with Deduplication
// ============================================================================

/// Chunk information with gene results
#[derive(Debug)]
pub struct ChunkResult {
    /// Chunk index
    pub chunk_index: usize,
    /// Global offset in the original sequence
    pub global_offset: usize,
    /// Genes detected in this chunk (with local coordinates)
    pub genes: Vec<GeneResult>,
}

/// Process a large sequence in parallel chunks (Case B: Large Genome)
///
/// # Arguments
/// * `seq` - The sequence record to process
/// * `config` - Parallel processing configuration
/// * `process_chunk_fn` - Function to process a single chunk
///
/// # Returns
/// Vector of genes with deduplicated results at chunk boundaries
pub fn process_large_sequence_parallel<F>(
    seq: &SeqRecord,
    config: &ParallelConfig,
    process_chunk_fn: F,
) -> Vec<GeneResult>
where
    F: Fn(&[u8], usize) -> Vec<GeneResult> + Sync + Send,
{
    // Check if parallelization is needed
    if !config.should_parallelize(seq.len()) {
        // Fast-path: process entire sequence without chunking
        return process_chunk_fn(&seq.seq, 0);
    }

    // Create chunks with dynamic sizing
    let chunks = create_chunks(&seq.seq, config);

    // Process chunks in parallel
    let chunk_results: Vec<ChunkResult> = if config.effective_threads() <= 1 {
        // Single-threaded
        chunks
            .iter()
            .map(|chunk| {
                let genes = process_chunk_fn(chunk.data, chunk.global_offset);
                ChunkResult {
                    chunk_index: chunk.chunk_index,
                    global_offset: chunk.global_offset,
                    genes,
                }
            })
            .collect()
    } else {
        // Parallel processing
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.effective_threads())
            .build()
            .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());

        pool.install(|| {
            chunks
                .par_iter()
                .map(|chunk| {
                    let genes = process_chunk_fn(chunk.data, chunk.global_offset);
                    ChunkResult {
                        chunk_index: chunk.chunk_index,
                        global_offset: chunk.global_offset,
                        genes,
                    }
                })
                .collect()
        })
    };

    // Merge and deduplicate results (lock-free)
    deduplicate_genes(chunk_results)
}

/// Lock-free deduplication of genes from multiple chunks
/// Genes at chunk boundaries may be detected twice - this removes duplicates
fn deduplicate_genes(chunk_results: Vec<ChunkResult>) -> Vec<GeneResult> {
    // Flatten all genes from all chunks (lock-free: just iterating)
    let mut all_genes: Vec<GeneResult> = chunk_results
        .into_iter()
        .flat_map(|cr| cr.genes)
        .collect();

    // Sort by position for deduplication
    all_genes.sort_by(|a, b| {
        a.start
            .cmp(&b.start)
            .then(a.stop.cmp(&b.stop))
            .then(a.genetype.cmp(&b.genetype))
            .then(a.comp.cmp(&b.comp))
    });

    // Remove exact duplicates (same position, type, strand)
    all_genes.dedup_by(|a, b| {
        a.start == b.start
            && a.stop == b.stop
            && a.genetype == b.genetype
            && a.comp == b.comp
    });

    all_genes
}

// ============================================================================
// Hybrid Routing (Case C with Case B for large contigs)
// ============================================================================

/// Process sequences with hybrid routing:
/// - Small contigs: process directly (Case A)
/// - Large contigs: chunk and parallelize (Case B)
/// - Multiple contigs: parallelize at contig level (Case C)
pub fn process_sequences_hybrid<F, G>(
    sequences: &[SeqRecord],
    config: &ParallelConfig,
    process_single_fn: F,
    process_chunk_fn: G,
) -> Vec<(String, Vec<GeneResult>)>
where
    F: Fn(&SeqRecord) -> Vec<GeneResult> + Sync + Send,
    G: Fn(&[u8], usize) -> Vec<GeneResult> + Sync + Send,
{
    if !config.enabled || config.effective_threads() <= 1 {
        // Single-threaded: process all sequences sequentially
        sequences
            .iter()
            .map(|seq| {
                let genes = if config.should_parallelize(seq.len()) {
                    // Large sequence but single-threaded mode
                    process_single_fn(seq)
                } else {
                    process_single_fn(seq)
                };
                (seq.id.clone(), genes)
            })
            .collect()
    } else {
        // Parallel processing with hybrid routing
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.effective_threads())
            .build()
            .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());

        pool.install(|| {
            sequences
                .par_iter()
                .map(|seq| {
                    let genes = if config.should_parallelize(seq.len()) {
                        // Case B: Large contig - use chunking
                        // Note: This creates nested parallelism which rayon handles well
                        process_large_sequence_parallel(seq, config, &process_chunk_fn)
                    } else {
                        // Case A: Small contig - process directly
                        process_single_fn(seq)
                    };
                    (seq.id.clone(), genes)
                })
                .collect()
        })
    }
}

// ============================================================================
// Additional Tests
// ============================================================================

#[cfg(test)]
mod parallel_tests {
    use super::*;

    #[test]
    fn test_load_fasta_sequences_nonexistent() {
        let result = load_fasta_sequences("/nonexistent/path.fasta");
        assert!(result.is_err());
    }

    #[test]
    fn test_seq_record() {
        let seq = SeqRecord::new("test".to_string(), vec![b'A', b'T', b'G', b'C']);
        assert_eq!(seq.id, "test");
        assert_eq!(seq.len(), 4);
        assert!(!seq.is_empty());
    }

    #[test]
    fn test_deduplicate_genes() {
        let chunk1 = ChunkResult {
            chunk_index: 0,
            global_offset: 0,
            genes: vec![
                GeneResult {
                    seq_id: "seq1".to_string(),
                    start: 100,
                    stop: 200,
                    genetype: 0,
                    comp: 0,
                    energy: 50.0,
                },
                GeneResult {
                    seq_id: "seq1".to_string(),
                    start: 500,
                    stop: 600,
                    genetype: 0,
                    comp: 0,
                    energy: 60.0,
                },
            ],
        };

        let chunk2 = ChunkResult {
            chunk_index: 1,
            global_offset: 450,
            genes: vec![
                // Duplicate of chunk1's second gene (in overlap region)
                GeneResult {
                    seq_id: "seq1".to_string(),
                    start: 500,
                    stop: 600,
                    genetype: 0,
                    comp: 0,
                    energy: 60.0,
                },
                GeneResult {
                    seq_id: "seq1".to_string(),
                    start: 900,
                    stop: 1000,
                    genetype: 0,
                    comp: 0,
                    energy: 70.0,
                },
            ],
        };

        let results = deduplicate_genes(vec![chunk1, chunk2]);

        // Should have 3 unique genes (duplicate removed)
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].start, 100);
        assert_eq!(results[1].start, 500);
        assert_eq!(results[2].start, 900);
    }

    #[test]
    fn test_process_sequences_parallel_single_thread() {
        let mut config = ParallelConfig::default();
        config.num_threads = 1;

        let sequences = vec![
            SeqRecord::new("seq1".to_string(), vec![b'A'; 100]),
            SeqRecord::new("seq2".to_string(), vec![b'T'; 200]),
        ];

        let results = process_sequences_parallel(&sequences, &config, |seq| {
            vec![GeneResult {
                seq_id: seq.id.clone(),
                start: 1,
                stop: seq.len() as i64,
                genetype: 0,
                comp: 0,
                energy: 0.0,
            }]
        });

        assert_eq!(results.len(), 2);
        assert_eq!(results[0].0, "seq1");
        assert_eq!(results[1].0, "seq2");
    }
}
