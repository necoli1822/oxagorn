/*
 * Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
 * thread_state.rs - Thread-local storage for parallel processing
 *
 * Based on ARAGORN v1.2.41 by Dean Laslett
 * Author: Sunju Kim <n.e.coli.1822@gmail.com>
 *
 * This module provides thread-local storage for the large working arrays
 * used by the gene detection functions. This enables safe parallel
 * processing of multiple sequences/chunks concurrently.
 */

#![allow(dead_code)]

use std::cell::UnsafeCell;
use crate::types::{
    TrnaLoop, TrnaDloop, MtTrnaLoop, MtTrnaTloop, MtTrnaAstem, MtCds, Gene,
    NTH, NC, NA, ND, NH,
    mtND, mtNTH, mtNA, mtNCDS,
    TRNA,
};

/// Thread-local state for tmioptimise function
/// These arrays are used as working buffers during gene detection
pub struct TmioptState {
    /// T-stem hits buffer
    pub thit: [TrnaLoop; NTH],
    /// C-stem hits buffer
    pub chit: [TrnaLoop; NC],
    /// Acceptor stem hits buffer
    pub ahit_arr: [TrnaLoop; NA],
    /// D-loop hits buffer
    pub dhit: [TrnaDloop; ND],
    /// C-stem encoding table (mutable during computation)
    pub ct: [u32; 6],
    /// Candidate gene buffer - reused to avoid 13KB memset per iteration
    /// Only te.energy needs to be reset; other fields are overwritten when candidate found
    pub te: Gene,
    /// Template gene - initialized once with default tRNA values
    /// Specific fields (tstem, tloop) updated per iteration without full reinit
    pub t: Gene,
}

impl Default for TmioptState {
    fn default() -> Self {
        // Initialize template gene t with default tRNA values (like C's static gene t)
        // These values match the C original: tmrna.c:334-336
        let mut t = Gene::default();
        t.astem1 = 7;
        t.astem2 = 7;
        t.spacer1 = 1;
        t.cstem = 2;
        t.cloop = 1;
        t.dstem = 3;
        t.dloop = 9;
        t.spacer2 = 5;
        t.var = 15;
        t.varbp = 0;
        t.intron = 0;
        t.nintron = 0;
        t.anticodon = 0;
        t.asst = 0;
        t.genetype = TRNA;
        t.energy = 0.0;

        Self {
            thit: [TrnaLoop::default(); NTH],
            chit: [TrnaLoop::default(); NC],
            ahit_arr: [TrnaLoop::default(); NA],
            dhit: [TrnaDloop::default(); ND],
            ct: [0u32; 6],
            te: Gene::default(),
            t,
        }
    }
}

/// Thread-local state for tmopt function (regular tmRNA)
pub struct TmoptState {
    /// R-hit buffer for tmRNA detection
    pub rhit: [TrnaLoop; NH],
}

impl Default for TmoptState {
    fn default() -> Self {
        Self {
            rhit: [TrnaLoop::default(); NH],
        }
    }
}

/// Thread-local state for find_mt_trna function (mitochondrial tRNA)
pub struct MtTrnaState {
    /// D-stem hits for mitochondrial tRNA
    pub dhit: [MtTrnaLoop; mtND + 1],
    /// T-stem hits for mitochondrial tRNA
    pub thit: [MtTrnaTloop; mtNTH + 1],
    /// Acceptor stem hits for mitochondrial tRNA
    pub ahit: [MtTrnaAstem; mtNA + 1],
    /// CDS hits for mitochondrial
    pub cdshit: [MtCds; mtNCDS],
    /// Template gene for optimization
    pub te: Gene,
    /// Encoding template
    pub tem: [u32; 6],
}

impl Default for MtTrnaState {
    fn default() -> Self {
        Self {
            dhit: [MtTrnaLoop::default(); mtND + 1],
            thit: [MtTrnaTloop::default(); mtNTH + 1],
            ahit: [MtTrnaAstem::default(); mtNA + 1],
            cdshit: [MtCds::default(); mtNCDS],
            te: Gene::default(),
            tem: [0u32; 6],
        }
    }
}

/// Thread-local state for find_astem5 function
pub struct Astem5State {
    /// Encoding template (modified during computation)
    pub tem: [u32; 6],
}

impl Default for Astem5State {
    fn default() -> Self {
        Self {
            tem: [0u32; 6],
        }
    }
}

/// Maximum capacity for thread-local gene storage
/// A. thaliana chromosome 1 has ~208 genes (highest observed)
/// Use 500 for safety margin while reducing memory footprint
pub const THREAD_LOCAL_GENE_CAPACITY: usize = 500;

/// Initial capacity for lazy allocation - most sequences have < 100 genes
const INITIAL_GENE_CAPACITY: usize = 100;

/// Thread-local gene storage for parallel processing
/// Each thread has its own gene array to avoid conflicts with the global TS
/// Uses lazy allocation: starts small and grows on demand to reduce memory usage
pub struct GeneStorage {
    /// Gene array - lazily allocated, grows on demand
    pub genes: Vec<Gene>,
    /// Maximum allowed capacity
    pub max_capacity: usize,
}

impl GeneStorage {
    /// Create new gene storage with lazy allocation
    /// Starts with small initial capacity, grows on demand
    pub fn new() -> Self {
        Self {
            genes: Vec::new(),  // Start empty - truly lazy
            max_capacity: THREAD_LOCAL_GENE_CAPACITY,
        }
    }

    /// Get the maximum capacity of this storage
    pub fn capacity(&self) -> usize {
        self.max_capacity
    }

    /// Get current allocated size
    pub fn current_size(&self) -> usize {
        self.genes.len()
    }

    /// Ensure the storage has at least `count` genes allocated and initialized
    /// This is the key function for lazy allocation
    pub fn ensure_capacity(&mut self, count: usize) {
        let target = count.min(self.max_capacity);
        if self.genes.len() < target {
            // Reserve space to avoid multiple reallocations
            self.genes.reserve(target - self.genes.len());
            // Push default genes until we reach target
            while self.genes.len() < target {
                self.genes.push(Gene::default());
            }
        }
    }

    /// Get a raw pointer to the gene array (for C-style access)
    /// IMPORTANT: Call ensure_capacity() before using this pointer
    pub fn as_mut_ptr(&mut self) -> *mut Gene {
        self.genes.as_mut_ptr()
    }

    /// Reset all allocated genes to default state
    pub fn reset(&mut self) {
        for gene in &mut self.genes {
            *gene = Gene::default();
        }
    }
}

// ============================================================================
// Thread-Local Storage using std::thread_local!
// Uses Box for heap allocation to avoid stack overflow with large structs
// Uses Option for lazy initialization to reduce memory usage
// ============================================================================

thread_local! {
    /// Thread-local state for tmioptimise (heap-allocated, lazily initialized)
    static TMIOPT_STATE: UnsafeCell<Option<Box<TmioptState>>> = const { UnsafeCell::new(None) };

    /// Thread-local state for tmopt (heap-allocated, lazily initialized)
    static TMOPT_STATE: UnsafeCell<Option<Box<TmoptState>>> = const { UnsafeCell::new(None) };

    /// Thread-local state for find_mt_trna (heap-allocated, lazily initialized)
    static MT_TRNA_STATE: UnsafeCell<Option<Box<MtTrnaState>>> = const { UnsafeCell::new(None) };

    /// Thread-local state for find_astem5 (heap-allocated, lazily initialized)
    static ASTEM5_STATE: UnsafeCell<Option<Box<Astem5State>>> = const { UnsafeCell::new(None) };

    /// Thread-local gene storage for parallel processing
    static GENE_STORAGE: std::cell::RefCell<Option<GeneStorage>> = const { std::cell::RefCell::new(None) };
}

// ============================================================================
// Accessor Functions for Thread-Local State
// Uses lazy initialization - allocates on first access
// ============================================================================

/// Get mutable reference to tmioptimise state for current thread
/// Lazily initializes on first access (heap-allocated)
///
/// # Safety
/// This function returns a raw pointer to thread-local storage.
/// The caller must ensure exclusive access within the current thread.
#[inline]
pub fn get_tmiopt_state() -> *mut TmioptState {
    TMIOPT_STATE.with(|cell| {
        let ptr = cell.get();
        unsafe {
            if (*ptr).is_none() {
                *ptr = Some(Box::new(TmioptState::default()));
            }
            (*ptr).as_mut().unwrap().as_mut() as *mut TmioptState
        }
    })
}

/// Get mutable reference to tmopt state for current thread
/// Lazily initializes on first access (heap-allocated)
#[inline]
pub fn get_tmopt_state() -> *mut TmoptState {
    TMOPT_STATE.with(|cell| {
        let ptr = cell.get();
        unsafe {
            if (*ptr).is_none() {
                *ptr = Some(Box::new(TmoptState::default()));
            }
            (*ptr).as_mut().unwrap().as_mut() as *mut TmoptState
        }
    })
}

/// Get mutable reference to mt_trna state for current thread
/// Lazily initializes on first access (heap-allocated)
#[inline]
pub fn get_mt_trna_state() -> *mut MtTrnaState {
    MT_TRNA_STATE.with(|cell| {
        let ptr = cell.get();
        unsafe {
            if (*ptr).is_none() {
                *ptr = Some(Box::new(MtTrnaState::default()));
            }
            (*ptr).as_mut().unwrap().as_mut() as *mut MtTrnaState
        }
    })
}

/// Get mutable reference to astem5 state for current thread
/// Lazily initializes on first access (heap-allocated)
#[inline]
pub fn get_astem5_state() -> *mut Astem5State {
    ASTEM5_STATE.with(|cell| {
        let ptr = cell.get();
        unsafe {
            if (*ptr).is_none() {
                *ptr = Some(Box::new(Astem5State::default()));
            }
            (*ptr).as_mut().unwrap().as_mut() as *mut Astem5State
        }
    })
}

/// Initialize and get thread-local gene storage with specified capacity
/// Returns a pointer to the gene array for C-style access
/// The storage is lazily allocated - only requested capacity is initialized
///
/// # Arguments
/// * `required_count` - Minimum number of genes that must be accessible
pub fn get_thread_local_ts_with_capacity(required_count: usize) -> *mut Gene {
    GENE_STORAGE.with(|cell| {
        let mut storage = cell.borrow_mut();
        if storage.is_none() {
            *storage = Some(GeneStorage::new());
        }
        let gs = storage.as_mut().unwrap();
        gs.ensure_capacity(required_count);
        gs.as_mut_ptr()
    })
}

/// Initialize and get thread-local gene storage (convenience wrapper)
/// Returns a pointer to the gene array for C-style access
/// Uses full capacity for backward compatibility with existing code
pub fn get_thread_local_ts() -> *mut Gene {
    get_thread_local_ts_with_capacity(THREAD_LOCAL_GENE_CAPACITY)
}

/// Reset thread-local gene storage to default state
pub fn reset_thread_local_ts() {
    GENE_STORAGE.with(|cell| {
        let mut storage = cell.borrow_mut();
        if let Some(ref mut gs) = *storage {
            gs.reset();
        }
    });
}

/// Get the maximum capacity of thread-local gene storage
pub fn get_thread_local_ts_capacity() -> i32 {
    THREAD_LOCAL_GENE_CAPACITY as i32
}

/// Get current allocated size of thread-local gene storage
/// Useful for debugging memory usage
pub fn get_thread_local_ts_current_size() -> usize {
    GENE_STORAGE.with(|cell| {
        let storage = cell.borrow();
        match storage.as_ref() {
            Some(gs) => gs.current_size(),
            None => 0,
        }
    })
}

// ============================================================================
// Convenience Macros (optional, for cleaner code)
// ============================================================================

/// Access tmiopt state arrays directly
/// Usage: with_tmiopt_state!(|state| { state.thit[0].pos = ...; })
#[macro_export]
macro_rules! with_tmiopt_state {
    ($body:expr) => {
        unsafe {
            let state = &mut *crate::thread_state::get_tmiopt_state();
            $body(state)
        }
    };
}

/// Access tmopt state arrays directly
#[macro_export]
macro_rules! with_tmopt_state {
    ($body:expr) => {
        unsafe {
            let state = &mut *crate::thread_state::get_tmopt_state();
            $body(state)
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;

    #[test]
    fn test_thread_local_isolation() {
        // Main thread: set a value
        unsafe {
            let state = &mut *get_tmiopt_state();
            state.ct[0] = 42;
        }

        // Spawn a new thread: verify it has its own state
        let handle = thread::spawn(|| {
            unsafe {
                let state = &mut *get_tmiopt_state();
                // Should be default (0), not 42
                assert_eq!(state.ct[0], 0);
                // Modify in child thread
                state.ct[0] = 100;
                assert_eq!(state.ct[0], 100);
            }
        });

        handle.join().unwrap();

        // Main thread: verify our value is unchanged
        unsafe {
            let state = &mut *get_tmiopt_state();
            assert_eq!(state.ct[0], 42);
        }
    }

    #[test]
    fn test_state_sizes() {
        // Log actual sizes for reference
        let tmiopt_size = std::mem::size_of::<TmioptState>();
        let tmopt_size = std::mem::size_of::<TmoptState>();
        let mt_size = std::mem::size_of::<MtTrnaState>();

        eprintln!("TmioptState size: {} KB", tmiopt_size / 1024);
        eprintln!("TmoptState size: {} KB", tmopt_size / 1024);
        eprintln!("MtTrnaState size: {} KB", mt_size / 1024);

        // Verify struct sizes are under 1 MB each (thread-local heap allocated)
        assert!(tmiopt_size < 1024 * 1024, "TmioptState > 1 MB");
        assert!(tmopt_size < 1024 * 1024, "TmoptState > 1 MB");
        assert!(mt_size < 1024 * 1024, "MtTrnaState > 1 MB");
    }
}
