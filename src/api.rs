//! Safe, high-level library API for oxagorn (tRNA / tmRNA detection).
//!
//! Wraps the raw `csw`-based detection ([`crate::output::detect_fastafile`]) so
//! consumers never touch unsafe pointers. [`detect`] builds a fully-initialized
//! search-config identical to the `oxagorn` binary's static setup (see
//! `main.rs`), runs the parallel detection in-process, and returns owned
//! [`DetectedGene`] records.

use crate::output::{detect_fastafile, DetectedGene};
use crate::types::*;

/// Options for a tRNA / tmRNA search.
#[derive(Clone, Copy, Debug)]
pub struct SearchOptions {
    /// Detect tRNA genes.
    pub trna: bool,
    /// Detect tmRNA genes.
    pub tmrna: bool,
    /// Mitochondrial tRNA mode (forces `trna` off, like the binary).
    pub mtrna: bool,
    /// Search both strands (`true`) or forward only (`false`).
    pub both_strands: bool,
    /// Linear sequence (`true`) or circular topology (`false`, the default).
    pub linear: bool,
    /// Genetic code id ([`STANDARD`] = 1 for bacteria).
    pub genetic_code: i32,
}

impl Default for SearchOptions {
    /// Bacterial default: tRNA + tmRNA, both strands, circular, standard code.
    fn default() -> Self {
        SearchOptions {
            trna: true,
            tmrna: true,
            mtrna: false,
            both_strands: true,
            linear: false,
            genetic_code: STANDARD,
        }
    }
}

/// Build a fully-initialized `Csw` matching the `oxagorn` binary's static csw
/// initialization (`main.rs`), then apply the mode flags from `opts`.
fn build_csw(opts: &SearchOptions) -> Csw {
    // SAFETY: Csw is plain-old-data (integers/floats/fixed arrays + a null FILE
    // pointer); zero-init then overwrite, exactly as the binary does.
    let mut sw: Csw = unsafe { std::mem::zeroed() };

    sw.f = std::ptr::null_mut();
    sw.both = if opts.both_strands { 2 } else { 0 };
    sw.geneticcode = opts.genetic_code;
    sw.discrim = METAZOAN_MT;
    sw.sp1max = 3;
    sw.sp2min = 0;
    sw.sp2max = 2;
    sw.mtxdetect = 1;
    sw.mtcdsscan = 1;
    sw.tmstrict = 1;
    sw.tvloop = 1;
    sw.extastem = 1;
    sw.showconfig = 1;
    sw.reportpsthresh = 100.0;
    sw.threshlevel = 1.0;
    sw.trnathresh = tRNAthresh;
    sw.ttscanthresh = 4.0;
    sw.ttarmthresh = 29.0;
    sw.tdarmthresh = 26.0;
    sw.tastemthresh = 7.5;
    sw.tascanthresh = 8.0;
    sw.mttthresh = mtRNAtthresh;
    sw.mtdthresh = mtRNAdthresh;
    sw.mtdtthresh = mtRNAdtthresh;
    sw.mttarmthresh = -7.9;
    sw.mtdarmthresh = -6.0;
    sw.tmrnathresh = tmRNAthresh;
    sw.tmathresh = 14.0;
    sw.tmcthresh = 10.0;
    sw.tmcathresh = 25.0;
    sw.tmrthresh = 9.0;
    sw.srpthresh = srpRNAthresh;
    sw.cdsthresh = CDSthresh;
    sw.tagthresh = NTAG as i32;

    // ARAGORN tmRNA signature block (main.rs).
    let sig: [i32; 100] = [
        45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
        45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
        45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
        10, 65, 82, 65, 71, 79, 82, 78, 32, 118,
        49, 46, 50, 46, 52, 49, 32, 32, 32, 68,
        101, 97, 110, 32, 76, 97, 115, 108, 101, 116,
        116, 10, 45, 45, 45, 45, 45, 45, 45, 45,
        45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
        45, 45, 45, 45, 45, 45, 45, 45, 45, 45,
        45, 45, 10, TERM, 0, 0, 0, 0, 0, 0,
    ];
    for i in 0..100 {
        sw.tmrna_struct[i] = sig[i];
    }

    // Modes.
    sw.trna = opts.trna as i32;
    sw.tmrna = opts.tmrna as i32;
    sw.mtrna = opts.mtrna as i32;
    if sw.mtrna != 0 {
        sw.trna = 0; // binary behaviour: mt mode disables nuclear tRNA
    }
    sw.linear = opts.linear as i32;
    sw.genespace = NT as i32;

    sw
}

/// Detect tRNA / tmRNA genes in a FASTA file, in-process, returning owned
/// records. `num_threads = 0` uses all available cores.
pub fn detect(
    filepath: &str,
    opts: &SearchOptions,
    num_threads: usize,
) -> Result<Vec<DetectedGene>, String> {
    let mut sw = build_csw(opts);
    // SAFETY: `sw` is fully initialized above; detect_fastafile reads config and
    // the FASTA file, never dereferencing the null output FILE pointer.
    unsafe { detect_fastafile(filepath, &mut sw as *mut Csw, num_threads) }
}
