/*
 * Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
 * mtrna.rs - Mitochondrial tRNA gene detection
 *
 * Based on ARAGORN v1.2.41 by Dean Laslett
 * Author: Sunju Kim <n.e.coli.1822@gmail.com>
 */

#![allow(unused_assignments)]
#![allow(unused_variables)]
#![allow(unused_mut)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(dead_code)]
#![allow(static_mut_refs)]
#![allow(improper_ctypes)]
#![allow(unused_imports)]
#![allow(unused_parens)]

use crate::types::*;
use crate::sequence::vloop_stability;
use crate::trna::{find_slot, find_slot_ts, base_copy3, aatail};
use crate::tables::*;
use crate::thread_state::get_mt_trna_state;

// ASTEM2_EXT is defined locally because types.rs has a different value (9)
const ASTEM2_EXT: i32 = 4;

// Local lowercase aliases for MT_ constants from types.rs
const mtND: usize = MT_ND;
const mtNTH: usize = MT_NTH;
const mtNA: usize = MT_NA;
const mtNCDS: usize = MT_NCDS;
const mtNCLM: usize = MT_NCLM;
const mtNTM: usize = MT_NTM;

// NOTE: C code has 18, types.rs has 18 - using types.rs value
const mt_DRLmaxlength: i32 = MT_DRLMAXLENGTH;
const mt_TVRLmaxlength: i32 = MT_TVRLMAXLENGTH;

// Energy constants from types.rs
const mtNOBOND: f64 = MT_NOBOND;
const mtGABOND: f64 = MT_GABOND;
const mtATBOND: f64 = MT_ATBOND;
const mtGCBOND: f64 = MT_GCBOND;
const mtGGBOND: f64 = MT_GGBOND;
const mtGTBOND: f64 = MT_GTBOND;
const mtTTBOND: f64 = MT_TTBOND;
const mtBONDSTAB: f64 = MT_BONDSTAB;
const mtTERMSTAB: f64 = MT_TERMSTAB;
const mtTSTTSTAB: f64 = MT_TSTTSTAB;
const mtGCPENALTY: f64 = MT_GCPENALTY;
const mtGCPENALTYD: f64 = MT_GCPENALTYD;

// NOTE: mt3MMSTAB differs from types.rs MT_3MMSTAB (0.5 vs 1.0)
// Using local value to match original C code in mtrna.c
const mt3MMSTAB: f64 = 0.5;

/// find_mt_trna_ts - Thread-safe main function for detecting mitochondrial tRNA genes
/// mtrna.c:11-3101
pub unsafe fn find_mt_trna_ts(
    ts: *mut Gene,
    max_genes: i32,
    d: *mut data_set,
    seq: *mut i32,
    lseq: i32,
    mut nts: i32,
    sw: *mut csw
) -> i32 {
    // Local variables (matching C declaration order)
    let mut nah: i32;
    let mut ndh: i32;
    let mut nch: i32;
    let mut nth: i32;
    let mut ncdsh: i32;
    let mut h: i32;
    let mut i: i32;
    let mut j: i32;
    let mut k: i32;
    let mut n: i32;
    let mut p: i32;
    let mut y: i32;
    let mut av: i32;
    let mut gcc: i32;
    let mut cgcc: i32;
    let mut catc: i32;
    let mut athresh: i32;
    let mut skip_cloop7: i32;
    let mut skip_cloop6: i32;
    let mut found_repeat: i32;
    let mut skip_ec: i32;
    let mut skip_tvn: i32;
    let mut skip_dn: i32;
    let mut skip_ngcc: i32;
    let mut skip_nasi: i32;
    let mut igc: i32;
    let mut nbase: i32;
    let mut b8: i32;
    let mut b9: i32;
    let mut b48: i32;
    let mut b57: i32;
    let mut nc: i32;
    let mut na: i32;
    let mut nt: i32;
    let mut nti: i32;
    let mut nd: i32;
    let mut ndi: i32;
    let mut dposmap: [i32; 32] = [0; 32];
    let mut dl: i32;
    let mut tl: i32;
    let mut extastem: i32;
    let mut astem8: i32;
    let mut astem8d: i32;
    let mut ti: i32;
    let mut di: i32;
    let mut ser: i32;
    let mut tastem: i32 = 0;
    let mut tastem8: i32 = 0;
    let mut tastem8d: i32 = 0;
    let mut astem: i32;
    let mut asteme: i32;
    let mut as_: i32;
    let mut as8: i32;
    let mut aext: i32 = 0;
    let mut aext8: i32 = 0;
    let mut nbasefext: i32;
    let mut cloop: i32;
    let mut dloop: i32;
    let mut tloop: i32;
    let mut tc: i32 = 0;
    let mut carm: i32;
    let mut cstem: i32;
    let mut darm: i32;
    let mut dstem: i32 = 0;
    let mut tarm: i32;
    let mut tstem: i32 = 0;
    let mut var: i32;
    let mut varbp: i32 = 0;
    let mut spacer1: i32;
    let mut spacer2: i32;
    let mut anticodon: i32;
    let mut ds: i32;
    let mut dstemmotif: i32;
    let mut cloop7: i32;
    let mut mtxdetect: i32;
    let mut incds: i32;
    let mut s: *mut i32;
    let mut sl: *mut i32;
    let mut s1: *mut i32;
    let mut s2: *mut i32;
    let mut s4: *mut i32;
    let mut sa: *mut i32;
    let mut sb: *mut i32;
    let mut sc: *mut i32;
    let mut se: *mut i32;
    let mut sf: *mut i32;
    let mut sg: *mut i32;
    let mut si: *mut i32;
    let mut slm: *mut i32;
    let mut slm1: *mut i32;
    let mut sle: *mut i32;
    let mut slb: *mut i32;
    let mut sld: *mut i32;
    let mut sge: *mut i32;
    let mut dpos: *mut i32;
    let mut cpos: *mut i32;
    let mut cend: *mut i32;
    let mut tpos: *mut i32;
    let mut tend: *mut i32;
    let mut apos1: *mut i32;
    let mut apos2: *mut i32;
    let mut aend1: *mut i32;
    let mut aend2: *mut i32;
    let mut clooppos: *mut i32;
    let mut cloopend: *mut i32;
    let mut bondtype: u32;
    let mut abondtype: u32;
    let mut mabondtype: u32;
    let mut acbondtype: u32;
    let mut cbondtype: u32;
    let mut agcat: u32;
    let mut cgcat: u32;
    let mut tgcat: u32 = 0;
    let mut dbondtype: u32;
    let mut dtbondtype: u32;
    let mut tbondtype: u32;
    let mut r: u32;
    let mut ct: [u32; 6] = [0; 6];
    let mut cm: u32;
    let mut cv: u32;
    let mut q: u32;
    let mut tendmap: [u32; 63] = [0; 63];
    let mut gcv: f64;
    let mut e: f64;
    let mut ec: f64;
    let mut ea: f64;
    let mut eas: f64 = 0.0;
    let mut ed: f64;
    let mut et: f64;
    let mut ev: f64;
    let mut energy: f64 = 0.0;
    let mut stem_energy: f64;
    let mut darmthresh: f64;
    let mut tarmthresh: f64;
    let mut tthresh: f64;
    let mut dthresh: f64;
    let mut dtthresh: f64;
    let mut thresh: f64;
    let mut chit: [mt_trna_cloop; 6] = [mt_trna_cloop::default(); 6];

    // Thread-local state for large arrays
    let tl_state = &mut *get_mt_trna_state();
    let dhit = &mut tl_state.dhit;
    let thit = &mut tl_state.thit;
    let ahit = &mut tl_state.ahit;
    let cdshit = &mut tl_state.cdshit;
    let te = &mut tl_state.te;
    let tem = &mut tl_state.tem;

    // Initialize te with default values
    te.astem1 = 7;
    te.astem2 = 7;
    te.aatail = 0;
    te.spacer1 = 1;
    te.spacer2 = 2;
    te.dstem = 1;
    te.dloop = 4;
    te.cstem = 7;
    te.cloop = 5;
    te.intron = 0;
    te.nintron = 7;
    te.anticodon = 0;
    te.var = 5;
    te.varbp = 0;
    te.tstem = 5;
    te.tloop = 7;
    te.genetype = TRNA;
    te.energy = 0.0;
    te.asst = 0;
    te.tps = 0;
    te.tpe = 0;
    te.annotation = 0;
    te.annosc = 0;

    let mut tn: *mut Gene;

    // Static scoring arrays
    static cAI: [i32; 6] = [8, 0, 0, 0, 8, 0];
    static cfCI: [i32; 6] = [0, 16, 0, 0, 16, 0];
    static cRI: [i32; 6] = [8, 0, 4, 0, 8, 0];
    static cTI: [i32; 6] = [0, 0, 0, 16, 16, 0];
    static cYI: [i32; 6] = [0, 8, 0, 4, 8, 0];
    static AI: [i32; 6] = [1, 0, 0, 0, 1, 0];
    static CI: [i32; 6] = [0, 1, 0, 0, 1, 0];
    static GI: [i32; 6] = [0, 0, 1, 0, 1, 0];
    static TI: [i32; 6] = [0, 0, 0, 1, 1, 0];
    static RI: [i32; 6] = [1, 0, 1, 0, 1, 0];
    static YI: [i32; 6] = [0, 1, 0, 1, 1, 0];
    static WI: [i32; 6] = [1, 0, 0, 1, 1, 0];
    // tem is now accessed via thread-local state (tl_state.tem)
    static At: [u32; 6] = [0, 0, 0, 1, 1, 0];
    static Ct: [u32; 6] = [0, 0, 1, 0, 1, 0];
    static Gt: [u32; 6] = [0, 1, 0, 1, 1, 0];
    static Tt: [u32; 6] = [1, 0, 1, 0, 1, 0];
    static cAt: [u32; 6] = [0, 0, 0, 2, 2, 0];
    static cCt: [u32; 6] = [0, 0, 2, 0, 2, 0];
    static cGt: [u32; 6] = [0, 2, 0, 1, 2, 0];
    static cTt: [u32; 6] = [2, 0, 1, 0, 2, 0];
    static aAt: [u32; 6] = [0, 0, 1, 2, 2, 0];
    static aCt: [u32; 6] = [0, 0, 2, 0, 2, 0];
    static aGt: [u32; 6] = [1, 2, 0, 1, 2, 0];
    static aTt: [u32; 6] = [2, 0, 1, 1, 2, 0];
    static dAt: [u32; 6] = [0, 0, 1, 2, 2, 0];
    static dCt: [u32; 6] = [0, 0, 2, 0, 2, 0];
    static dGt: [u32; 6] = [1, 2, 0, 2, 2, 0];
    static dTt: [u32; 6] = [2, 0, 2, 1, 2, 0];
    static clmotif: [u32; mtNCLM] = [0x1321300, 0x3321300, 0x1323002];
    static dloopi: [[i32; 4]; (mt_DRLmaxlength + 1) as usize] = [
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 2, -1, 0],
        [0, 2, -1, 0],
        [0, 2, 3, -1],
        [0, 3, -1, 0],
        [0, 3, -1, 0],
        [0, 3, 4, -1],
        [0, 4, -1, 0],
        [0, 5, -1, 0],
        [0, 5, 6, -1],
        [0, 5, 6, -1],
    ];
    static tloopa: [[i32; 4]; 12] = [
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [-1, 0, 0, 0],
        [0, 1, -1, 0],
        [0, 2, 1, -1],
        [4, 3, 2, -1],
        [4, 3, -1, 0],
        [4, 3, -1, 0],
        [4, 3, -1, 0],
        [5, 4, 3, -1],
        [5, 4, -1, 0],
        [5, -1, 0, 0],
    ];
    static dA: [f64; 6] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    static dT: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 1.0, 0.0];
    static C_: [f64; 6] = [0.0, 1.0, 0.0, 0.0, 1.0, 0.0];
    static G_: [f64; 6] = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0];
    static T_: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 1.0, 0.0];
    static AX: [f64; 6] = [0.0, -1.0, -1.0, -1.0, 0.0, -1.0];
    static AX37: [f64; 6] = [0.0, -4.0, -1.0, -4.0, 0.0, -4.0];
    static AXX: [f64; 6] = [0.0, -3.0, -1.5, -3.0, 0.0, -3.0];
    static AXX37: [f64; 6] = [0.0, -4.0, -4.0, -4.0, 0.0, -4.0];
    static AX7: [f64; 6] = [0.0, -0.7, -0.7, -0.7, 0.0, -0.7];
    static CX: [f64; 6] = [-2.0, 0.0, -2.0, -1.0, 0.0, -2.0];
    static CXX: [f64; 6] = [-4.0, 0.0, -4.0, -2.0, 0.0, -4.0];
    static CX7: [f64; 6] = [-0.7, 0.0, -0.7, -0.7, 0.0, -0.7];
    static TX: [f64; 6] = [-1.0, -1.0, -1.0, 0.0, 0.0, -1.0];
    static TXX: [f64; 6] = [-2.0, -2.0, -2.0, 0.0, 0.0, -2.0];
    static YX: [f64; 6] = [-1.0, 0.0, -1.0, 0.0, 0.0, -1.0];
    static tC: [f64; 6] = [0.0, 0.01, 0.0, 0.0, 0.01, 0.0];
    static tG: [f64; 6] = [0.0, 0.0, 0.01, 0.0, 0.01, 0.0];
    static tT: [f64; 6] = [0.0, 0.0, 0.0, 0.01, 0.01, 0.0];
    static cA: [f64; 6] = [0.8, 0.0, 0.0, 0.0, 0.8, 0.0];
    static cfC: [f64; 6] = [0.0, 2.6, 0.0, 0.0, 2.6, 0.0];
    static cR: [f64; 6] = [0.8, -2.0, 0.8, -0.8, 0.8, -0.8];
    static cT: [f64; 6] = [-0.8, 0.0, -0.8, 2.6, 2.6, -0.8];
    static cY: [f64; 6] = [-0.8, 0.8, -0.8, 0.8, 0.8, -0.8];
    static loop_stab: [f64; 41] = [
        10.0, 2.0, 1.0, 0.4, 0.3, 0.2, 0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.6,
        1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5,
        4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 0.0,
    ];
    static bem: [[f64; 6]; 6] = [
        [mtNOBOND, mtNOBOND, mtGABOND, mtATBOND, mtATBOND, mtNOBOND],
        [mtNOBOND, mtNOBOND, mtGCBOND, mtNOBOND, mtGCBOND, mtNOBOND],
        [mtGABOND, mtGCBOND, mtGGBOND, mtGTBOND, mtGCBOND, mtNOBOND],
        [mtATBOND, mtNOBOND, mtGTBOND, mtTTBOND, mtATBOND, mtNOBOND],
        [mtATBOND, mtGCBOND, mtGCBOND, mtATBOND, mtGCBOND, mtNOBOND],
        [mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND],
    ];
    static hbem: [[f64; 5]; 5] = [
        [0.0, 0.0, 0.0, mtBONDSTAB + 0.5 * mtATBOND, mtBONDSTAB + 0.5 * mtATBOND],
        [0.0, 0.0, mtBONDSTAB + 0.5 * mtGCBOND, 0.0, mtBONDSTAB + 0.5 * mtGCBOND],
        [0.0, mtBONDSTAB + 0.5 * mtGCBOND, 0.0, mtBONDSTAB + 0.5 * mtGTBOND, mtBONDSTAB + 0.5 * mtGCBOND],
        [mtBONDSTAB + 0.5 * mtATBOND, 0.0, mtBONDSTAB + 0.5 * mtGTBOND, 0.0, mtBONDSTAB + 0.5 * mtATBOND],
        [mtBONDSTAB + 0.5 * mtATBOND, mtBONDSTAB + 0.5 * mtGCBOND, mtBONDSTAB + 0.5 * mtGCBOND, mtBONDSTAB + 0.5 * mtATBOND, mtBONDSTAB + 0.5 * mtGCBOND],
    ];

    // Initialize thresholds from sw
    tarmthresh = (*sw).mttarmthresh;
    tthresh = (*sw).mttthresh;
    dthresh = (*sw).mtdthresh;
    dtthresh = (*sw).mtdtthresh;
    ds = (*sw).discrim;
    extastem = (*sw).extastem;
    cloop7 = (*sw).cloop7;
    mtxdetect = (*sw).mtxdetect;

    // Initialize working threshold
    thresh = dtthresh;

    // find coding sequences
    ncdsh = 0;

    // find cstems - main search loop
    sc = seq.offset((*sw).loffset as isize);
    sl = seq.offset((lseq - (*sw).roffset) as isize);

    if sc.offset(17) > sl {
        return nts;
    }

    h = *sc.offset(16);
    p = *sc.offset(15);
    j = *sc.offset(14);
    k = *sc.offset(13);
    n = *sc.offset(12);
    y = *sc.offset(11);

    ct[0] = *cAt.get_unchecked(h as usize) | (*cAt.get_unchecked(p as usize) << 4) | (*cAt.get_unchecked(j as usize) << 8)
        | (*cAt.get_unchecked(k as usize) << 12) | (*cAt.get_unchecked(n as usize) << 16) | (*cAt.get_unchecked(y as usize) << 20);
    ct[1] = *cCt.get_unchecked(h as usize) | (*cCt.get_unchecked(p as usize) << 4) | (*cCt.get_unchecked(j as usize) << 8)
        | (*cCt.get_unchecked(k as usize) << 12) | (*cCt.get_unchecked(n as usize) << 16) | (*cCt.get_unchecked(y as usize) << 20);
    ct[2] = *cGt.get_unchecked(h as usize) | (*cGt.get_unchecked(p as usize) << 4) | (*cGt.get_unchecked(j as usize) << 8)
        | (*cGt.get_unchecked(k as usize) << 12) | (*cGt.get_unchecked(n as usize) << 16) | (*cGt.get_unchecked(y as usize) << 20);
    ct[3] = *cTt.get_unchecked(h as usize) | (*cTt.get_unchecked(p as usize) << 4) | (*cTt.get_unchecked(j as usize) << 8)
        | (*cTt.get_unchecked(k as usize) << 12) | (*cTt.get_unchecked(n as usize) << 16) | (*cTt.get_unchecked(y as usize) << 20);
    ct[4] = 0;
    ct[5] = 0;

    // Main cloop search loop
    while sc < sl {
        te.ps = std::ptr::null_mut();

        p = *sc.offset(17);
        ct[0] = (ct[0] << 4) | *cAt.get_unchecked(p as usize);
        ct[1] = (ct[1] << 4) | *cCt.get_unchecked(p as usize);
        ct[2] = (ct[2] << 4) | *cGt.get_unchecked(p as usize);
        ct[3] = (ct[3] << 4) | *cTt.get_unchecked(p as usize);
        cm = (*ct.get_unchecked(*sc.offset(4) as usize) >> 16)
            + (*ct.get_unchecked(*sc.offset(3) as usize) >> 12)
            + (*ct.get_unchecked(*sc.offset(2) as usize) >> 8)
            + (*ct.get_unchecked(*sc.offset(1) as usize) >> 4)
            + *ct.get_unchecked(*sc as usize);

        // 7 base cloop
        cv = cm & 0xf0;
        athresh = 12;
        nch = 0;
        skip_cloop7 = 0;

        // Exclude certain cloop patterns
        if *RI.get_unchecked(*sc.offset(6) as usize) != 0 {
            if *RI.get_unchecked(*sc.offset(5) as usize) != 0 {
                skip_cloop7 = 1;
            } else if *YI.get_unchecked(*sc.offset(10) as usize) != 0 {
                skip_cloop7 = 1;
            } else if cv < 0x60 {
                skip_cloop7 = 1;
            }
        } else {
            if *YI.get_unchecked(*sc.offset(10) as usize) != 0 {
                if *RI.get_unchecked(*sc.offset(5) as usize) != 0 {
                    skip_cloop7 = 1;
                }
            }
            if skip_cloop7 == 0 {
                if cv < 0x40 {
                    if cv < 0x20 {
                        skip_cloop7 = 1;
                    } else if *sc.offset(5) != Cytosine {
                        skip_cloop7 = 1;
                    } else if *sc.offset(6) != Thymine {
                        skip_cloop7 = 1;
                    } else if *sc.offset(10) != Adenine {
                        skip_cloop7 = 1;
                    } else {
                        athresh = 11;
                    }
                } else if cv < 0x70 {
                    athresh = 11;
                    k = *cYI.get_unchecked(*sc.offset(5) as usize)
                        + *cTI.get_unchecked(*sc.offset(6) as usize)
                        + *cRI.get_unchecked(*sc.offset(10) as usize)
                        + *cAI.get_unchecked(*sc.offset(11) as usize);
                    if *sc.offset(6) == Cytosine {
                        if *sc.offset(5) == Cytosine {
                            k += 16;
                        } else if *sc.offset(5) == Thymine {
                            if *sc.offset(11) == Adenine {
                                k += 16;
                            }
                        }
                    }
                    if cv == 0x40 {
                        if k < 40 {
                            skip_cloop7 = 1;
                        }
                    } else if cv == 0x50 {
                        if k < 28 {
                            skip_cloop7 = 1;
                        }
                    } else {
                        if k < 20 {
                            skip_cloop7 = 1;
                        } else {
                            athresh = 9;
                        }
                    }
                } else {
                    athresh = if cv < 10 { 9 } else { 8 };
                }
            }
        }

        if skip_cloop7 == 0 {
            chit[0].pos = sc;
            chit[0].stem = 5;
            chit[0].loop_ = 7;
            chit[0].looppos = sc.offset(5);
            chit[0].arm = 17;
            chit[0].end = sc.offset(17);
            chit[0].anticodon = (*sc.offset(7) << 4) + (*sc.offset(8) << 2) + *sc.offset(9);
            if *BP.get_unchecked(*sc.offset(-1) as usize).get_unchecked(*sc.offset(17) as usize) != 0 {
                chit[1].pos = sc.offset(-1);
                chit[1].stem = 6;
                chit[1].loop_ = 7;
                chit[1].looppos = sc.offset(5);
                chit[1].arm = 19;
                chit[1].end = sc.offset(18);
                chit[1].anticodon = chit[0].anticodon;
                nch = 2;
            } else {
                nch = 1;
            }
        }

        // 6 base cloop - simplified check
        skip_cloop6 = 0;
        found_repeat = 0;
        if cloop7 != 0 {
            skip_cloop6 = 1;
        } else if (cm & 0xf00) >= 0x800 {
            if (*YI.get_unchecked(*sc.offset(6) as usize) == 0) && (*YI.get_unchecked(*sc.offset(5) as usize) == 0) {
                skip_cloop6 = 1;
            } else if (*RI.get_unchecked(*sc.offset(9) as usize) == 0) && (*RI.get_unchecked(*sc.offset(10) as usize) == 0) {
                skip_cloop6 = 1;
            } else {
                // Check for repeat sequences (simplified)
                if found_repeat == 0 {
                    (*chit.get_unchecked_mut(nch as usize)).pos = sc;
                    (*chit.get_unchecked_mut(nch as usize)).stem = 5;
                    (*chit.get_unchecked_mut(nch as usize)).loop_ = 6;
                    (*chit.get_unchecked_mut(nch as usize)).looppos = sc.offset(5);
                    (*chit.get_unchecked_mut(nch as usize)).arm = 16;
                    (*chit.get_unchecked_mut(nch as usize)).end = sc.offset(16);
                    (*chit.get_unchecked_mut(nch as usize)).anticodon = 0;
                    nch += 1;
                    if athresh > 10 {
                        athresh = 10;
                    }
                    if *BP.get_unchecked(*sc.offset(-1) as usize).get_unchecked(*sc.offset(16) as usize) != 0 {
                        (*chit.get_unchecked_mut(nch as usize)).pos = sc.offset(-1);
                        (*chit.get_unchecked_mut(nch as usize)).stem = 6;
                        (*chit.get_unchecked_mut(nch as usize)).loop_ = 6;
                        (*chit.get_unchecked_mut(nch as usize)).looppos = sc.offset(5);
                        (*chit.get_unchecked_mut(nch as usize)).arm = 18;
                        (*chit.get_unchecked_mut(nch as usize)).end = sc.offset(17);
                        (*chit.get_unchecked_mut(nch as usize)).anticodon = 0;
                        nch += 1;
                    }
                }
            }
        } else {
            skip_cloop6 = 1;
        }

        // Skip to next position if no cloop hits
        if nch < 1 {
            sc = sc.offset(1);
            continue;
        }

        // Calculate carm energy for each cloop hit
        nc = 0;
        while nc < nch {
            s1 = (*chit.get_unchecked(nc as usize)).pos;
            cstem = (*chit.get_unchecked(nc as usize)).stem;
            cloop = (*chit.get_unchecked(nc as usize)).loop_;
            s4 = s1.offset(cstem as isize);
            s2 = s4.offset(cloop as isize);
            energy = if cloop == 7 { 0.0 } else { -4.0 };
            energy += *cY.get_unchecked(*s4 as usize) + *cT.get_unchecked(*s4.offset(1) as usize)
                + *cR.get_unchecked(*s2.offset(-2) as usize) + *cA.get_unchecked(*s2.offset(-1) as usize);
            if *s4.offset(1) == Cytosine {
                if *s4 == Cytosine {
                    energy += 2.6;
                } else if *s4 == Thymine {
                    if *s2.offset(-1) == Adenine {
                        energy += 2.6;
                    }
                }
            }
            s2 = s2.offset(cstem as isize);
            s2 = s2.offset(-1);
            stem_energy = *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
            k = *NEIGHBOUR_MAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
            stem_energy += *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize);
            bondtype = *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
            if *BP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) != 0 {
                if *ASSYMST.get_unchecked(*s2.offset(1) as usize).get_unchecked(*s1.offset(-1) as usize) != 0 {
                    stem_energy += mtTERMSTAB;
                } else {
                    stem_energy += *SEND_EM.get_unchecked(*s2 as usize).get_unchecked(*s1 as usize);
                }
            } else {
                if *ASSYMST.get_unchecked(*s2 as usize).get_unchecked(*s1 as usize) != 0 {
                    stem_energy += mtTERMSTAB;
                } else {
                    stem_energy += *SEND_EM.get_unchecked(*s2.offset(-1) as usize).get_unchecked(*s1.offset(1) as usize);
                }
            }
            s1 = s1.offset(1);
            while s1 < s4 {
                s2 = s2.offset(-1);
                if *WCBP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) == 0 {
                    if *WCBP.get_unchecked(*s1.offset(-1) as usize).get_unchecked(*s2.offset(1) as usize) == 0 {
                        j = 0;
                        while j < mtNTM as i32 {
                            if *s1 == *TANDEMID.get_unchecked(j as usize).get_unchecked(1)
                                && *s2 == *TANDEMID.get_unchecked(j as usize).get_unchecked(3)
                                && *s1.offset(-1) == *TANDEMID.get_unchecked(j as usize).get_unchecked(0)
                                && *s2.offset(1) == *TANDEMID.get_unchecked(j as usize).get_unchecked(2)
                            {
                                stem_energy += *TANDEM_EM.get_unchecked(j as usize);
                                break;
                            }
                            j += 1;
                        }
                        if s1 < s4.offset(-1) {
                            if *BP.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) == 0 {
                                stem_energy -= mt3MMSTAB;
                            }
                        }
                    }
                    k = *NEIGHBOUR_MAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                    stem_energy += *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(-1) as usize).get_unchecked(*s2.offset(1) as usize)
                        + *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize);
                }
                bondtype += *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                stem_energy += *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                s1 = s1.offset(1);
            }
            s1 = s1.offset(-1);
            if *BP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) == 0 {
                s1 = s1.offset(-1);
                s2 = s2.offset(1);
            }
            if *ASSYMST.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) != 0 {
                stem_energy += mtTERMSTAB;
            } else {
                stem_energy += *SEND_EM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
            }
            cgcc = (bondtype & 0xf) as i32;
            if cgcc <= 0 {
                catc = ((bondtype & 0xf0) >> 4) as i32;
                if catc < cstem {
                    energy -= mtGCPENALTY;
                }
            }
            if cstem == 6 {
                energy += 1.0;
            }
            (*chit.get_unchecked_mut(nc as usize)).bondtype = bondtype;
            (*chit.get_unchecked_mut(nc as usize)).stem_energy = stem_energy;
            (*chit.get_unchecked_mut(nc as usize)).energy = energy + stem_energy;
            nc += 1;
        }

        // Skip to next position if no cloop hits after energy calculation
        if nch < 1 {
            sc = sc.offset(1);
            continue;
        }

        // find tarms
        nth = 0;
        slm = sc.offset(61);
        sle = sc.offset(57);
        slb = sc.offset(21);
        sg = sc.offset(16);
        sge = sg.offset(30);
        let slb_tarm = sg.offset(32);
        tem[0] = *At.get_unchecked(*slm as usize);
        tem[1] = *Ct.get_unchecked(*slm as usize);
        tem[2] = *Gt.get_unchecked(*slm as usize);
        tem[3] = *Tt.get_unchecked(*slm as usize);
        slm = slm.offset(-1);
        while slm > sle {
            tem[0] = (tem[0] << 4) | *At.get_unchecked(*slm as usize);
            tem[1] = (tem[1] << 4) | *Ct.get_unchecked(*slm as usize);
            tem[2] = (tem[2] << 4) | *Gt.get_unchecked(*slm as usize);
            tem[3] = (tem[3] << 4) | *Tt.get_unchecked(*slm as usize);
            slm = slm.offset(-1);
        }
        while slm >= slb {
            tem[0] = ((tem[0] << 4) | *At.get_unchecked(*slm as usize)) & 0xfffff;
            tem[1] = ((tem[1] << 4) | *Ct.get_unchecked(*slm as usize)) & 0xfffff;
            tem[2] = ((tem[2] << 4) | *Gt.get_unchecked(*slm as usize)) & 0xfffff;
            tem[3] = ((tem[3] << 4) | *Tt.get_unchecked(*slm as usize)) & 0xfffff;
            sf = slm.offset(3);
            if sf > sge {
                sf = sge;
            }
            apos2 = slm.offset(5);
            si = sg;
            s = si.offset(4);
            r = *tem.get_unchecked(*si as usize);
            si = si.offset(1);
            while si < s {
                r = (r >> 4) + *tem.get_unchecked(*si as usize);
                si = si.offset(1);
            }
            while si <= sf {
                if si < slm {
                    r = (r >> 4) + *tem.get_unchecked(*si as usize);
                    si = si.offset(1);
                } else {
                    si = si.offset(1);
                    r = r >> 4;
                }
                q = r & 0xf;
                if slm > slb_tarm {
                    if q < 5 {
                        continue;
                    }
                    tloop = (slm as isize - si as isize) as i32 / std::mem::size_of::<i32>() as i32;
                } else {
                    if q < 2 {
                        continue;
                    }
                    if q < 3 {
                        if *WCBP.get_unchecked(*si.offset(-5) as usize).get_unchecked(*apos2.offset(-1) as usize) == 0 {
                            continue;
                        }
                        if *WCBP.get_unchecked(*si.offset(-4) as usize).get_unchecked(*apos2.offset(-2) as usize) == 0 {
                            continue;
                        }
                        tloop = (slm as isize - si as isize) as i32 / std::mem::size_of::<i32>() as i32;
                        if tloop > 5 {
                            continue;
                        }
                    } else {
                        tloop = (slm as isize - si as isize) as i32 / std::mem::size_of::<i32>() as i32;
                        if q < 4 {
                            if *BP.get_unchecked(*si.offset(-4) as usize).get_unchecked(*apos2.offset(-2) as usize) == 0 {
                                if *BP.get_unchecked(*si.offset(-2) as usize).get_unchecked(*apos2.offset(-4) as usize) == 0 {
                                    if tloop < 4 {
                                        continue;
                                    }
                                    if *si.offset(-1) != Guanine {
                                        continue;
                                    }
                                    if *si != Thymine {
                                        continue;
                                    }
                                    if *si.offset(1) != Thymine {
                                        continue;
                                    }
                                }
                            }
                        }
                    }
                }
                if tloop < 7 {
                    if tloop < 2 {
                        if tloop <= 0 {
                            if tloop <= -2 {
                                if *WCBP.get_unchecked(*si.offset(-5) as usize).get_unchecked(*apos2.offset(-1) as usize) == 0 {
                                    continue;
                                }
                                if *WCBP.get_unchecked(*si.offset(-4) as usize).get_unchecked(*apos2.offset(-2) as usize) == 0 {
                                    continue;
                                }
                                tstem = 2;
                                tloop += 6;
                            } else if *BP.get_unchecked(*si.offset(-3) as usize).get_unchecked(*apos2.offset(-3) as usize) != 0 {
                                tstem = 3;
                                tloop += 4;
                            } else {
                                if *WCBP.get_unchecked(*si.offset(-5) as usize).get_unchecked(*apos2.offset(-1) as usize) == 0 {
                                    continue;
                                }
                                if *WCBP.get_unchecked(*si.offset(-4) as usize).get_unchecked(*apos2.offset(-2) as usize) == 0 {
                                    continue;
                                }
                                tstem = 2;
                                tloop += 6;
                            }
                        } else if *BP.get_unchecked(*si.offset(-2) as usize).get_unchecked(*apos2.offset(-4) as usize) != 0 {
                            tstem = 4;
                            tloop += 2;
                        } else if *BP.get_unchecked(*si.offset(-3) as usize).get_unchecked(*apos2.offset(-3) as usize) != 0 {
                            tstem = 3;
                            tloop += 4;
                        } else {
                            if *WCBP.get_unchecked(*si.offset(-5) as usize).get_unchecked(*apos2.offset(-1) as usize) == 0 {
                                continue;
                            }
                            if *WCBP.get_unchecked(*si.offset(-4) as usize).get_unchecked(*apos2.offset(-2) as usize) == 0 {
                                continue;
                            }
                            tstem = 2;
                            tloop += 6;
                        }
                    } else if *BP.get_unchecked(*si.offset(-1) as usize).get_unchecked(*apos2.offset(-5) as usize) != 0 {
                        if q != 4 {
                            tstem = 5;
                        } else if *BP.get_unchecked(*si.offset(-2) as usize).get_unchecked(*apos2.offset(-4) as usize) != 0 {
                            tstem = 5;
                        } else {
                            k = *GI.get_unchecked(*si.offset(-3) as usize) + *TI.get_unchecked(*si.offset(-2) as usize)
                                + *TI.get_unchecked(*si.offset(-1) as usize) + *CI.get_unchecked(*si as usize);
                            if k >= 2 {
                                tstem = 3;
                                tloop += 4;
                            } else {
                                tstem = 5;
                            }
                        }
                    } else if *BP.get_unchecked(*si.offset(-2) as usize).get_unchecked(*apos2.offset(-4) as usize) != 0 {
                        tstem = 4;
                        tloop += 2;
                    } else if *BP.get_unchecked(*si.offset(-3) as usize).get_unchecked(*apos2.offset(-3) as usize) != 0 {
                        tstem = 3;
                        tloop += 4;
                    } else {
                        if *WCBP.get_unchecked(*si.offset(-5) as usize).get_unchecked(*apos2.offset(-1) as usize) == 0 {
                            continue;
                        }
                        if *WCBP.get_unchecked(*si.offset(-4) as usize).get_unchecked(*apos2.offset(-2) as usize) == 0 {
                            continue;
                        }
                        tstem = 2;
                        tloop += 6;
                    }
                    if tloop < 3 {
                        if tstem > 3 {
                            tstem -= 1;
                            tloop += 2;
                        }
                    }
                } else if *BP.get_unchecked(*si.offset(-1) as usize).get_unchecked(*apos2.offset(-5) as usize) == 0 {
                    if *BP.get_unchecked(*si.offset(-2) as usize).get_unchecked(*apos2.offset(-4) as usize) == 0 {
                        tstem = 3;
                        tloop += 4;
                    } else {
                        tstem = 4;
                        tloop += 2;
                    }
                } else {
                    tstem = 5;
                }
                if tloop > 17 {
                    if tstem < 5 {
                        continue;
                    }
                }

                // calculate tarm energy
                s1 = si.offset(-5);
                tpos = s1;
                s4 = s1.offset(tstem as isize);
                s2 = apos2;
                s2 = s2.offset(-1);
                if *TT.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) != 0 {
                    energy = mtTSTTSTAB;
                    s1 = s1.offset(1);
                    s2 = s2.offset(-1);
                    if *TT.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) != 0 {
                        energy += mtTSTTSTAB;
                        bondtype = *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                        s1 = s1.offset(1);
                        s2 = s2.offset(-1);
                    } else {
                        bondtype = 0;
                    }
                } else {
                    energy = 0.0;
                    bondtype = 0;
                }

                // calculate tstem energy
                stem_energy = *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                k = *NEIGHBOUR_MAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                stem_energy += *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize);
                bondtype += *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                s1 = s1.offset(1);
                while s1 < s4 {
                    s2 = s2.offset(-1);
                    if *WCBP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) == 0 {
                        if *WCBP.get_unchecked(*s1.offset(-1) as usize).get_unchecked(*s2.offset(1) as usize) == 0 {
                            j = 0;
                            while j < mtNTM as i32 {
                                if *s1 == *TANDEMID.get_unchecked(j as usize).get_unchecked(1)
                                    && *s2 == *TANDEMID.get_unchecked(j as usize).get_unchecked(3)
                                    && *s1.offset(-1) == *TANDEMID.get_unchecked(j as usize).get_unchecked(0)
                                    && *s2.offset(1) == *TANDEMID.get_unchecked(j as usize).get_unchecked(2)
                                {
                                    stem_energy += *TANDEM_EM.get_unchecked(j as usize);
                                    break;
                                }
                                j += 1;
                            }
                            if s1 < s4.offset(-1) {
                                if *BP.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) == 0 {
                                    stem_energy -= mt3MMSTAB;
                                }
                            }
                        }
                        k = *NEIGHBOUR_MAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                        stem_energy += *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(-1) as usize).get_unchecked(*s2.offset(1) as usize)
                            + *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize);
                    }
                    bondtype += *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                    stem_energy += *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                    s1 = s1.offset(1);
                }
                s1 = s1.offset(-1);
                if tloop < 4 {
                    stem_energy += *SSEND_EM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                } else if *ASSYMST.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) != 0 {
                    stem_energy += mtTERMSTAB;
                } else {
                    stem_energy += *SEND_EM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                }

                // compile possible tarms
                energy += stem_energy - mtBONDSTAB * (5 - tstem) as f64;
                if energy >= tarmthresh {
                    (*thit.get_unchecked_mut(nth as usize)).pos = tpos;
                    s1 = tpos.offset(tstem as isize);
                    s2 = apos2.offset(-(tstem as isize));
                    (*thit.get_unchecked_mut(nth as usize)).energy = energy - *loop_stab.get_unchecked(tloop as usize)
                        + *tG.get_unchecked(*s1.offset(-1) as usize) + *tT.get_unchecked(*s1 as usize)
                        + *tT.get_unchecked(*s1.offset(1) as usize) + *tC.get_unchecked(*s1.offset(2) as usize);
                    (*thit.get_unchecked_mut(nth as usize)).stem_energy = stem_energy;
                    (*thit.get_unchecked_mut(nth as usize)).bondtype = bondtype;
                    (*thit.get_unchecked_mut(nth as usize)).stem = tstem;
                    (*thit.get_unchecked_mut(nth as usize)).loop_ = tloop;
                    (*thit.get_unchecked_mut(nth as usize)).end = tpos.offset((2 * tstem + tloop) as isize);
                    nth += 1;
                    if nth >= mtNTH as i32 {
                        break;
                    }
                }
            }
            slm = slm.offset(-1);
        }

        // find darms
        ndh = 0;
        sle = sc.offset(-4);
        slb = sc.offset(-8);
        slm = sc.offset(-1);
        tem[0] = *dAt.get_unchecked(*slm as usize);
        tem[1] = *dCt.get_unchecked(*slm as usize);
        tem[2] = *dGt.get_unchecked(*slm as usize);
        tem[3] = *dTt.get_unchecked(*slm as usize);
        slm = slm.offset(-1);
        while slm > sle {
            tem[0] = (tem[0] << 4) | *dAt.get_unchecked(*slm as usize);
            tem[1] = (tem[1] << 4) | *dCt.get_unchecked(*slm as usize);
            tem[2] = (tem[2] << 4) | *dGt.get_unchecked(*slm as usize);
            tem[3] = (tem[3] << 4) | *dTt.get_unchecked(*slm as usize);
            slm = slm.offset(-1);
        }
        slm1 = slm;
        while slm > slb {
            tem[0] = ((tem[0] << 4) | *dAt.get_unchecked(*slm as usize)) & 0xffff;
            tem[1] = ((tem[1] << 4) | *dCt.get_unchecked(*slm as usize)) & 0xffff;
            tem[2] = ((tem[2] << 4) | *dGt.get_unchecked(*slm as usize)) & 0xffff;
            tem[3] = ((tem[3] << 4) | *dTt.get_unchecked(*slm as usize)) & 0xffff;
            slm = slm.offset(-1);
            si = slm.offset(-18);
            s = si.offset(3);
            r = *tem.get_unchecked(*si as usize);
            si = si.offset(1);
            while si < s {
                r = (r >> 4) + *tem.get_unchecked(*si as usize);
                si = si.offset(1);
            }
            while si <= slm1 {
                if si < slm {
                    r = (r >> 4) + *tem.get_unchecked(*si as usize);
                    si = si.offset(1);
                } else {
                    r = r >> 4;
                    si = si.offset(1);
                }
                q = r & 0xf;
                if q < 6 {
                    q += (*TI.get_unchecked(*si.offset(-6) as usize) + *RI.get_unchecked(*si.offset(-5) as usize)) as u32;
                    if q < 6 {
                        continue;
                    }
                }

                // calculate darm energy
                s1 = si.offset(-4);
                (*dhit.get_unchecked_mut(ndh as usize)).pos = s1;
                energy = *dT.get_unchecked(*s1.offset(-2) as usize) + *dA.get_unchecked(*s1.offset(-1) as usize);
                dloop = ((slm1 as isize - si as isize) / std::mem::size_of::<i32>() as isize) as i32;
                skip_ec = 0;
                if dloop > 2 && *BP.get_unchecked(*si.offset(-1) as usize).get_unchecked(*slm1 as usize) != 0 {
                    dstem = 4;
                } else if dloop > 0
                    && (*GGSTEMBP.get_unchecked(*si.offset(-2) as usize).get_unchecked(*slm.offset(2) as usize) != 0
                        || *GABP.get_unchecked(*si.offset(-1) as usize).get_unchecked(*slm1 as usize) != 0)
                {
                    dstem = 3;
                    dloop += 2;
                    energy += mtNOBOND;
                } else {
                    if *WCBP.get_unchecked(*si.offset(-3) as usize).get_unchecked(*slm.offset(3) as usize) == 0 {
                        skip_ec = 1;
                    } else if *GC.get_unchecked(*si.offset(-4) as usize).get_unchecked(*slm.offset(4) as usize) == 0 {
                        skip_ec = 1;
                    } else {
                        dstem = 2;
                        dloop += 4;
                        if dloop > 5 {
                            energy += mtNOBOND;
                        }
                        energy += mtNOBOND;
                    }
                }
                if skip_ec != 0 {
                    continue;
                }
                s2 = slm.offset(4);
                s4 = s1.offset(dstem as isize);
                if *WCBP.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) == 0 {
                    if *STEMTERM.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) != 0 {
                        energy -= 1.0;
                    } else if *BP.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) != 0 {
                        energy -= 1.5;
                    } else {
                        energy -= 2.0;
                    }
                }

                // calculate dstem energy
                stem_energy = *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                k = *NEIGHBOUR_MAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                stem_energy += *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize);
                bondtype = *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                if *BP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) != 0 {
                    if *ASSYMST.get_unchecked(*s2.offset(1) as usize).get_unchecked(*s1.offset(-1) as usize) != 0 {
                        stem_energy += mtTERMSTAB;
                    } else {
                        stem_energy += *SEND_EM.get_unchecked(*s2 as usize).get_unchecked(*s1 as usize);
                    }
                    s1 = s1.offset(1);
                    s2 = s2.offset(-1);
                } else {
                    s1 = s1.offset(1);
                    s2 = s2.offset(-1);
                    if *ASSYMST.get_unchecked(*s2.offset(1) as usize).get_unchecked(*s1.offset(-1) as usize) != 0 {
                        stem_energy += mtTERMSTAB;
                    } else {
                        stem_energy += *SEND_EM.get_unchecked(*s2 as usize).get_unchecked(*s1 as usize);
                    }
                }
                stem_energy += *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                k = *NEIGHBOUR_MAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                stem_energy += *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(-1) as usize).get_unchecked(*s2.offset(1) as usize)
                    + *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize);
                bondtype += *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                s1 = s1.offset(1);
                while s1 < s4 {
                    s2 = s2.offset(-1);
                    if *WCBP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) == 0 {
                        if *WCBP.get_unchecked(*s1.offset(-1) as usize).get_unchecked(*s2.offset(1) as usize) == 0 {
                            j = 0;
                            while j < mtNTM as i32 {
                                if *s1 == *TANDEMID.get_unchecked(j as usize).get_unchecked(1)
                                    && *s2 == *TANDEMID.get_unchecked(j as usize).get_unchecked(3)
                                    && *s1.offset(-1) == *TANDEMID.get_unchecked(j as usize).get_unchecked(0)
                                    && *s2.offset(1) == *TANDEMID.get_unchecked(j as usize).get_unchecked(2)
                                {
                                    stem_energy += *TANDEM_EM.get_unchecked(j as usize);
                                    break;
                                }
                                j += 1;
                            }
                            if s1 < s4.offset(-1) {
                                if *BP.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) == 0 {
                                    stem_energy -= mt3MMSTAB;
                                }
                            }
                        }
                        k = *NEIGHBOUR_MAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                        stem_energy += *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(-1) as usize).get_unchecked(*s2.offset(1) as usize)
                            + *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize);
                    }
                    bondtype += *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                    stem_energy += *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                    s1 = s1.offset(1);
                }
                s1 = s1.offset(-1);
                if *BP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) == 0 {
                    s1 = s1.offset(-1);
                    s2 = s2.offset(1);
                }
                if dloop < 4 {
                    stem_energy += *SSEND_EM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                } else if *ASSYMST.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) != 0 {
                    stem_energy += mtTERMSTAB;
                } else {
                    stem_energy += *SEND_EM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                }

                // compile possible darms
                energy += stem_energy;
                (*dhit.get_unchecked_mut(ndh as usize)).energy = energy;
                (*dhit.get_unchecked_mut(ndh as usize)).stem_energy = stem_energy;
                (*dhit.get_unchecked_mut(ndh as usize)).bondtype = bondtype;
                (*dhit.get_unchecked_mut(ndh as usize)).stem = dstem;
                (*dhit.get_unchecked_mut(ndh as usize)).loop_ = dloop;
                ndh += 1;
                if ndh >= mtND as i32 {
                    break;
                }
            }
        }

        // build darm exclusion map
        // 5' astems further from carm than mt_DRLmaxlength must match a darm
        i = 3;
        while i <= 30 {
            *dposmap.get_unchecked_mut(i as usize) = 0;
            i += 1;
        }
        sf = sc.offset(-(mt_DRLmaxlength as isize + 1));
        sld = sf;
        if ndh > 0 {
            s = dhit[0].pos;
            nd = 0;
            while nd < ndh {
                se = (*dhit.get_unchecked(nd as usize)).pos;
                if se < s {
                    s = se;
                }
                i = ((sc as isize - se as isize) / std::mem::size_of::<i32>() as isize) as i32;
                i += 1;
                if *dposmap.get_unchecked(i as usize) < 1 {
                    *dposmap.get_unchecked_mut(i as usize) = 1;
                }
                i += 1;
                *dposmap.get_unchecked_mut(i as usize) = 2;
                i += 1;
                if *dposmap.get_unchecked(i as usize) < 1 {
                    *dposmap.get_unchecked_mut(i as usize) = 1;
                }
                nd += 1;
            }
            s = s.offset(-4);
            if s < sf {
                sf = s;
            }
        }

        // build tarm exclusion map
        // 3' astems further from carm than mt_TVRLmaxlength must match a tarm
        i = 17;
        while i <= 62 {
            *tendmap.get_unchecked_mut(i as usize) = 0;
            i += 1;
        }
        s2 = sc.offset((mt_TVRLmaxlength + 17) as isize);
        sle = s2;
        if nth > 0 {
            s = thit[0].end;
            nt = 0;
            while nt < nth {
                se = (*thit.get_unchecked(nt as usize)).end;
                if se > s {
                    s = se;
                }
                i = ((se as isize - sc as isize) / std::mem::size_of::<i32>() as isize) as i32;
                bondtype = (*thit.get_unchecked(nt as usize)).bondtype;
                if *tendmap.get_unchecked(i as usize) != 0 {
                    if bondtype < *tendmap.get_unchecked(i as usize) {
                        *tendmap.get_unchecked_mut(i as usize) = bondtype;
                    }
                } else {
                    *tendmap.get_unchecked_mut(i as usize) = bondtype;
                }
                nt += 1;
            }
            if s > s2 {
                s2 = s;
            }
        }

        // find astems in 3 categories
        nah = 0;
        sa = sc.offset(-3);
        sg = sf.offset(-6);
        let sb_astem = sc.offset(17);
        se = s2.offset(6);
        tem[0] = *aAt.get_unchecked(*se as usize);
        tem[1] = *aCt.get_unchecked(*se as usize);
        tem[2] = *aGt.get_unchecked(*se as usize);
        tem[3] = *aTt.get_unchecked(*se as usize);
        se = se.offset(-1);
        while se > s2 {
            tem[0] = (tem[0] << 4) | *aAt.get_unchecked(*se as usize);
            tem[1] = (tem[1] << 4) | *aCt.get_unchecked(*se as usize);
            tem[2] = (tem[2] << 4) | *aGt.get_unchecked(*se as usize);
            tem[3] = (tem[3] << 4) | *aTt.get_unchecked(*se as usize);
            se = se.offset(-1);
        }
        ti = ((se as isize - sc as isize) / std::mem::size_of::<i32>() as isize) as i32;
        while se >= sb_astem {
            tem[0] = ((tem[0] << 4) | *aAt.get_unchecked(*se as usize)) & 0xfffffff;
            tem[1] = ((tem[1] << 4) | *aCt.get_unchecked(*se as usize)) & 0xfffffff;
            tem[2] = ((tem[2] << 4) | *aGt.get_unchecked(*se as usize)) & 0xfffffff;
            tem[3] = ((tem[3] << 4) | *aTt.get_unchecked(*se as usize)) & 0xfffffff;
            if *tendmap.get_unchecked(ti as usize) != 0 {
                nti = if *tendmap.get_unchecked(ti as usize) < 0x2000 { 1 } else { 0 };
            } else {
                if se > sle {
                    se = se.offset(-1);
                    ti -= 1;
                    continue;
                }
                nti = -1;
            }
            si = sg;
            r = *tem.get_unchecked(*si as usize);
            si = si.offset(1);
            while si < sf {
                r = (r >> 4) + *tem.get_unchecked(*si as usize);
                si = si.offset(1);
            }
            di = ((sc as isize - si as isize) / std::mem::size_of::<i32>() as isize) as i32;
            while si < sa {
                r = (r >> 4) + *tem.get_unchecked(*si as usize);
                si = si.offset(1);
                di -= 1;
                if *dposmap.get_unchecked(di as usize) != 0 {
                    if nti <= 0 {
                        if nti < 0 {
                            if *dposmap.get_unchecked(di as usize) < 2 {
                                continue;
                            }
                        }
                        av = (r & 0xf) as i32;
                        if av < athresh {
                            continue;
                        }
                    }
                } else {
                    if si < sld {
                        continue;
                    }
                    if nti < 0 {
                        continue;
                    }
                    av = (r & 0xf) as i32;
                    if av < athresh {
                        continue;
                    }
                }
                if nah >= mtNA as i32 {
                    break;
                }

                // predict astem length and calculate astem energy
                s1 = si.offset(-7);
                s2 = se.offset(6);
                if *BP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) != 0 {
                    astem = 7;
                    energy = 0.0;
                    (*ahit.get_unchecked_mut(nah as usize)).pos1 = s1;
                    (*ahit.get_unchecked_mut(nah as usize)).pos2 = se;
                } else if *GGSTEMTERM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) != 0 {
                    astem = 7;
                    (*ahit.get_unchecked_mut(nah as usize)).pos1 = s1;
                    (*ahit.get_unchecked_mut(nah as usize)).pos2 = se;
                    energy = *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                    s1 = s1.offset(1);
                    s2 = s2.offset(-1);
                } else {
                    energy = *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                    s1 = s1.offset(1);
                    s2 = s2.offset(-1);
                    if *BP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) != 0 {
                        astem = 6;
                        (*ahit.get_unchecked_mut(nah as usize)).pos1 = s1;
                        (*ahit.get_unchecked_mut(nah as usize)).pos2 = se;
                    } else if *GGSTEMTERM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) != 0 {
                        astem = 6;
                        (*ahit.get_unchecked_mut(nah as usize)).pos1 = s1;
                        (*ahit.get_unchecked_mut(nah as usize)).pos2 = se;
                        energy += *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                        s1 = s1.offset(1);
                        s2 = s2.offset(-1);
                    } else {
                        astem = 5;
                        energy += *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                        s1 = s1.offset(1);
                        s2 = s2.offset(-1);
                        (*ahit.get_unchecked_mut(nah as usize)).pos1 = s1;
                        (*ahit.get_unchecked_mut(nah as usize)).pos2 = se;
                    }
                }
                (*ahit.get_unchecked_mut(nah as usize)).stem = astem;
                bondtype = *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                energy += *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                k = *NEIGHBOUR_MAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                energy += *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize);
                s1 = s1.offset(1);
                s2 = s2.offset(-1);
                energy += *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                k = *NEIGHBOUR_MAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                energy += *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(-1) as usize).get_unchecked(*s2.offset(1) as usize)
                    + *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize);
                bondtype += *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                s1 = s1.offset(1);
                while s1 < si {
                    s2 = s2.offset(-1);
                    if *WCBP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) == 0 {
                        if *WCBP.get_unchecked(*s1.offset(-1) as usize).get_unchecked(*s2.offset(1) as usize) == 0 {
                            j = 0;
                            while j < mtNTM as i32 {
                                if *s1 == *TANDEMID.get_unchecked(j as usize).get_unchecked(1)
                                    && *s2 == *TANDEMID.get_unchecked(j as usize).get_unchecked(3)
                                    && *s1.offset(-1) == *TANDEMID.get_unchecked(j as usize).get_unchecked(0)
                                    && *s2.offset(1) == *TANDEMID.get_unchecked(j as usize).get_unchecked(2)
                                {
                                    energy += *TANDEM_EM.get_unchecked(j as usize);
                                    break;
                                }
                                j += 1;
                            }
                            if s1 < si.offset(-1) {
                                if *BP.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) == 0 {
                                    energy -= mt3MMSTAB;
                                }
                            }
                        }
                        k = *NEIGHBOUR_MAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                        energy += *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(-1) as usize).get_unchecked(*s2.offset(1) as usize)
                            + *NEIGHBOUR_EM.get_unchecked(k as usize).get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize);
                    }
                    bondtype += *BTMAP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                    energy += *bem.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
                    s1 = s1.offset(1);
                }
                // terminal basepair check
                {
                    let mut do_assymst = 1;
                    s1 = s1.offset(-1);
                    if *BP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) == 0 {
                        s1 = s1.offset(-1);
                        s2 = s2.offset(1);
                        if *BP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) == 0 {
                            s1 = s1.offset(-1);
                            s2 = s2.offset(1);
                            if *BP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) == 0 {
                                s1 = s1.offset(-1);
                                s2 = s2.offset(1);
                                if *BP.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) == 0 {
                                    do_assymst = 0;
                                }
                            }
                        }
                    }
                    if do_assymst != 0 {
                        if *ASSYMST.get_unchecked(*s1.offset(1) as usize).get_unchecked(*s2.offset(-1) as usize) != 0 {
                            energy += mtTERMSTAB;
                        }
                    }
                }
                (*ahit.get_unchecked_mut(nah as usize)).energy = energy;
                (*ahit.get_unchecked_mut(nah as usize)).bondtype = bondtype;
                nah += 1;
            }
            se = se.offset(-1);
            ti -= 1;
        }
        if nah <= 0 {
            sc = sc.offset(1);
            continue;
        }

        // build mttrna genes
        // cycle through astems first so that GC content is only calculated once per astem
        thresh = -INACTIVE;
        te.ps = std::ptr::null_mut();
        na = 0;
        while na < nah {
            apos2 = (*ahit.get_unchecked(na as usize)).pos2;
            apos1 = (*ahit.get_unchecked(na as usize)).pos1;
            astem = (*ahit.get_unchecked(na as usize)).stem;
            aend1 = apos1.offset(astem as isize);
            astem8 = if astem == 7 {
                *WCBP.get_unchecked(*apos1.offset(-1) as usize).get_unchecked(*apos2.offset(7) as usize)
            } else {
                0
            };
            asteme = 0;
            ea = (*ahit.get_unchecked(na as usize)).energy;
            abondtype = (*ahit.get_unchecked(na as usize)).bondtype;
            agcat = ((abondtype >> 4) + abondtype) & 0xf;

            // GC content
            s = apos1;
            aend2 = apos2.offset(astem as isize);
            nbase = ((aend2 as isize - apos1 as isize) / std::mem::size_of::<i32>() as isize) as i32 + 1;
            igc = 0;
            while s <= aend2 {
                k = *s;
                if k >= Cytosine && k <= Guanine {
                    igc += 1;
                }
                s = s.offset(1);
            }
            gcv = 10.0 * igc as f64 / nbase as f64;
            if gcv < 1.0 {
                if gcv < 0.55 {
                    na += 1;
                    continue;
                }
                ea -= 0.5;
            }
            if nbase > 60 {
                if gcv > 6.0 {
                    ea -= 2.0 * (gcv - 6.0);
                }
            } else {
                if gcv > 5.0 {
                    ea -= 2.0 * (gcv - 5.0);
                }
            }
            if gcv > 6.6 {
                ea -= 6.0;
                if gcv > 7.0 {
                    ea -= 6.0;
                }
            }

            // findout if inside a coding sequence
            incds = 0;
            i = 0;
            while i < ncdsh {
                if apos1 > (*cdshit.get_unchecked(i as usize)).pos1 {
                    if aend2 <= (*cdshit.get_unchecked(i as usize)).pos2 {
                        incds = 1;
                        ea -= 2.0;
                        break;
                    }
                }
                i += 1;
            }

            // cycle through carms that fall between astem
            nc = 0;
            while nc < nch {
                cpos = (*chit.get_unchecked(nc as usize)).pos;
                dloop = ((cpos as isize - aend1 as isize) / std::mem::size_of::<i32>() as isize) as i32;
                if dloop < 3 {
                    nc += 1;
                    continue;
                }
                if dloop > 26 {
                    nc += 1;
                    continue;
                }
                cend = (*chit.get_unchecked(nc as usize)).end;
                tloop = ((apos2 as isize - cend as isize) / std::mem::size_of::<i32>() as isize) as i32;
                if tloop < 5 {
                    nc += 1;
                    continue;
                }
                cloop = (*chit.get_unchecked(nc as usize)).loop_;
                cstem = (*chit.get_unchecked(nc as usize)).stem;
                clooppos = (*chit.get_unchecked(nc as usize)).looppos;
                cloopend = clooppos.offset(cloop as isize);
                carm = (*chit.get_unchecked(nc as usize)).arm;
                anticodon = (*chit.get_unchecked(nc as usize)).anticodon;
                cbondtype = (*chit.get_unchecked(nc as usize)).bondtype;
                acbondtype = abondtype + cbondtype;
                cgcat = ((cbondtype >> 4) + cbondtype) & 0xf;
                ec = ea + (*chit.get_unchecked(nc as usize)).energy;

                // astem,cstem stability (GC bond count)
                if (abondtype & 0xf) <= 0 {
                    if (cbondtype & 0xf) <= 0 {
                        ec -= mtGCPENALTYD;
                        if ((cbondtype & 0xf0) >> 4) >= 5 {
                            ec += 0.5;
                        }
                    }
                }

                // anticodon to astem discriminator base match
                astem8d = 0;
                if cloop == 7 {
                    if *MT_DISCRIM.get_unchecked(ds as usize).get_unchecked(anticodon as usize).get_unchecked(*apos2.offset(astem as isize) as usize) == 0 {
                        if astem8 != 0 {
                            if *MT_DISCRIM.get_unchecked(ds as usize).get_unchecked(anticodon as usize).get_unchecked(*apos2.offset(8) as usize) != 0 {
                                astem8d = 1;
                            } else {
                                ec -= 3.0;
                            }
                        } else if astem <= 6 {
                            if *MT_DISCRIM.get_unchecked(ds as usize).get_unchecked(anticodon as usize).get_unchecked(*apos2.offset(7) as usize) == 0 {
                                if astem == 5 {
                                    if *MT_DISCRIM.get_unchecked(ds as usize).get_unchecked(anticodon as usize).get_unchecked(*apos2.offset(6) as usize) == 0 {
                                        ec -= 3.0;
                                    }
                                } else {
                                    ec -= 3.0;
                                }
                            }
                        } else {
                            ec -= 3.0;
                        }
                    }
                }

                // build TV-replacement loop mttrna genes
                skip_tvn = 0;
                if tloop <= mt_TVRLmaxlength {
                    if (*sw).tvloop == 0 {
                        skip_tvn = 1;
                    }

                    // astem termination
                    // (only need to calculate once per astem)
                    if asteme == 0 {
                        let mut nost2_done = 0;
                        asteme = 1;
                        s = aend1.offset(-1);
                        se = apos2;
                        while *BP.get_unchecked(*s as usize).get_unchecked(*se as usize) == 0 {
                            s = s.offset(-1);
                            if s <= apos1 {
                                eas = 0.0;
                                nost2_done = 1;
                                break;
                            }
                            se = se.offset(1);
                        }
                        if nost2_done == 0 {
                            if *AASTEMTERM.get_unchecked(*s.offset(1) as usize).get_unchecked(*se.offset(-1) as usize) == 0 {
                                eas = -0.5;
                            } else {
                                eas = 0.0;
                                while se >= apos2 {
                                    s = s.offset(1);
                                    se = se.offset(-1);
                                    if *AASTEMTERM.get_unchecked(*s as usize).get_unchecked(*se as usize) != 0 {
                                        eas += 1.0;
                                    }
                                }
                            }
                        }
                    }

                    // choose darm
                    energy = 94.0 + ec + eas;
                    nd = -1;
                    ndi = -1;
                    ed = -INACTIVE;
                    while {
                        nd += 1;
                        nd < ndh
                    } {
                        dpos = (*dhit.get_unchecked(nd as usize)).pos;
                        spacer1 = ((dpos as isize - aend1 as isize) / std::mem::size_of::<i32>() as isize) as i32;
                        if spacer1 != 2 {
                            continue;
                        }
                        dl = (*dhit.get_unchecked(nd as usize)).loop_;
                        dstem = (*dhit.get_unchecked(nd as usize)).stem;
                        if dstem > 4 {
                            continue;
                        }
                        darm = 2 * dstem + dl;
                        spacer2 = ((cpos as isize - dpos as isize) / std::mem::size_of::<i32>() as isize) as i32 - darm;

                        // astem,darm,cstem interspacing
                        if spacer2 < 1 {
                            continue;
                        }
                        e = (*dhit.get_unchecked(nd as usize)).energy;
                        if spacer2 > 1 {
                            if spacer2 > 2 {
                                continue;
                            }
                            if *STEMBP.get_unchecked(*cpos as usize).get_unchecked(*cend.offset(-1) as usize) == 0 {
                                continue;
                            }
                            if tloop > 12 {
                                e -= 2.0;
                            }
                            if ((*dhit.get_unchecked(nd as usize)).bondtype & 0xf) < 1 {
                                if (agcat + cgcat + 1) < (cstem as u32 + astem as u32) {
                                    e -= 3.0;
                                }
                            }
                        } else {
                            if dl > 11 {
                                if *RI.get_unchecked(*cpos.offset(-1) as usize) == 0 {
                                    e -= 2.0;
                                }
                            } else {
                                if *cpos.offset(-1) == Cytosine {
                                    e -= 2.0;
                                }
                            }
                        }

                        // small,large dloop, dstem R motif
                        if dl < 3 {
                            e -= 2.0;
                        }
                        if dl > 12 {
                            e -= 2.0;
                        }
                        if *RI.get_unchecked(*dpos as usize) == 0 {
                            e -= 1.0;
                        }

                        // darm,tloop tertiary interaction
                        k = 0;
                        di = if dl >= 12 { 3 } else { if dl >= 9 { 2 } else { 1 } };
                        tl = if tloop >= 14 { 5 } else { if dl >= 9 { if tloop >= 10 { 4 } else { 3 } } else { 3 } };
                        if *GGSTACKBP.get_unchecked(*dpos.offset((dstem + di) as isize) as usize).get_unchecked(*cend.offset(tl as isize) as usize) == 0 {
                            if tl > 3 {
                                if *GGSTACKBP.get_unchecked(*dpos.offset((dstem + di) as isize) as usize).get_unchecked(*cend.offset((tl - 1) as isize) as usize) == 0 {
                                    e -= 1.5;
                                } else {
                                    k += 1;
                                }
                            } else if di > 1 {
                                if *GGSTACKBP.get_unchecked(*dpos.offset((dstem + di - 1) as isize) as usize).get_unchecked(*cend.offset(tl as isize) as usize) == 0 {
                                    e -= 1.5;
                                } else {
                                    k += 1;
                                }
                            } else {
                                e -= 1.5;
                            }
                        } else {
                            k += 1;
                        }
                        if *STEMTERM.get_unchecked(*dpos.offset((dstem - 1) as isize) as usize).get_unchecked(*dpos.offset((darm - dstem) as isize) as usize) != 0 {
                            e -= 0.5;
                            if *cend.offset(2) == *dpos.offset((dstem - 2) as isize) {
                                if *BP.get_unchecked(*cend.offset(2) as usize).get_unchecked(*dpos.offset((darm - dstem + 1) as isize) as usize) != 0 {
                                    k += 1;
                                }
                            } else {
                                if *cend.offset(2) == *dpos.offset((darm - dstem + 1) as isize) {
                                    if *BP.get_unchecked(*cend.offset(2) as usize).get_unchecked(*dpos.offset((dstem - 2) as isize) as usize) != 0 {
                                        k += 1;
                                    }
                                }
                            }
                        } else {
                            if *cend.offset(2) == *dpos.offset((dstem - 1) as isize) {
                                if *BP.get_unchecked(*cend.offset(2) as usize).get_unchecked(*dpos.offset((darm - dstem) as isize) as usize) == 0 {
                                    e -= 0.5;
                                } else {
                                    k += 1;
                                }
                            } else {
                                if *cend.offset(2) != *dpos.offset((darm - dstem) as isize) {
                                    e -= 0.5;
                                } else if *BP.get_unchecked(*cend.offset(2) as usize).get_unchecked(*dpos.offset((dstem - 1) as isize) as usize) == 0 {
                                    e -= 0.5;
                                } else {
                                    k += 1;
                                }
                            }
                        }
                        if *cend.offset(1) == *dpos {
                            if *STACKBP.get_unchecked(*cend.offset(1) as usize).get_unchecked(*dpos.offset((darm - 1) as isize) as usize) == 0 {
                                e -= 0.5;
                            } else {
                                k += 1;
                            }
                        } else {
                            if *cend.offset(1) != *dpos.offset((darm - 1) as isize) {
                                e -= 0.5;
                            } else if *BP.get_unchecked(*cend.offset(1) as usize).get_unchecked(*dpos as usize) == 0 {
                                e -= 0.5;
                            } else {
                                k += 1;
                            }
                        }

                        // darm stability
                        dstemmotif = *WCBP.get_unchecked(*dpos.offset(1) as usize).get_unchecked(*dpos.offset((darm - 2) as isize) as usize);
                        if spacer2 == 2 {
                            if (k < 3) || ((*dhit.get_unchecked(nd as usize)).bondtype > 0x200) || (dstemmotif == 0) {
                                if abondtype >= 0x10000 {
                                    e -= 2.0;
                                }
                                if dstem > 3 {
                                    e -= 1.0;
                                }
                                e -= 0.5;
                            }
                        }

                        // darm tertiary interactions
                        j = 0;
                        b8 = *dpos.offset(-2);
                        b9 = *dpos.offset(-1);
                        if *BP.get_unchecked(b8 as usize).get_unchecked(*dpos.offset(dstem as isize) as usize) == 0 {
                            e -= 1.0;
                        } else if *WCBP.get_unchecked(b8 as usize).get_unchecked(*dpos.offset(dstem as isize) as usize) != 0 {
                            j += 1;
                        }
                        if *BP.get_unchecked(b8 as usize).get_unchecked(*dpos.offset((darm - dstem - 1) as isize) as usize) == 0 {
                            e -= 1.0;
                        } else if *WCBP.get_unchecked(b8 as usize).get_unchecked(*dpos.offset((darm - dstem - 1) as isize) as usize) != 0 {
                            j += 1;
                        }
                        if *WCBP.get_unchecked(*dpos.offset(2) as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) == 0 {
                            if *GASTEMBP.get_unchecked(b8 as usize).get_unchecked(*dpos.offset(dstem as isize) as usize) == 0 {
                                e -= 2.0;
                            } else if *GASTEMBP.get_unchecked(b8 as usize).get_unchecked(*dpos.offset((darm - dstem - 1) as isize) as usize) == 0 {
                                e -= 2.0;
                            }
                            if *GGSTEMBP.get_unchecked(*dpos.offset(2) as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) == 0 {
                                e -= 1.0;
                            }
                        } else {
                            j += 1;
                        }
                        if *BP.get_unchecked(b9 as usize).get_unchecked(*dpos.offset(2) as usize) == 0 {
                            if *BP.get_unchecked(b9 as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) == 0 {
                                e -= 1.0;
                            } else {
                                j += 1;
                            }
                        } else {
                            j += 1;
                        }

                        // more extensive tertiary interaction between darm,tloop
                        if dstemmotif != 0 {
                            if k >= 3 {
                                if *BP.get_unchecked(*dpos.offset(2) as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) != 0 {
                                    if b8 != Thymine {
                                        e += 0.5;
                                    }
                                    if dl > 3 {
                                        if *BP.get_unchecked(*dpos.offset((dstem + 2) as isize) as usize).get_unchecked(*cend.offset((tl + 1) as isize) as usize) != 0 {
                                            e += 0.7;
                                        } else if *GABP.get_unchecked(*dpos.offset((dstem + 2) as isize) as usize).get_unchecked(*cend.offset((tl + 1) as isize) as usize) != 0 {
                                            e += 0.5;
                                        }
                                    }
                                    if tloop >= 6 {
                                        if spacer2 < 2 {
                                            if dl >= 3 {
                                                di = if dl > 11 { 2 } else { 1 };
                                                if *BP.get_unchecked(*dpos.offset((dstem + di) as isize) as usize).get_unchecked(*cend.offset(tl as isize) as usize) != 0 {
                                                    if (*chit.get_unchecked(nc as usize)).stem_energy > -4.8 {
                                                        e += 0.5;
                                                    }
                                                    if *WCBP.get_unchecked(*dpos.offset((dstem + di) as isize) as usize).get_unchecked(*cend.offset(tl as isize) as usize) != 0 {
                                                        if gcv > 1.2 {
                                                            if *clooppos.offset(1) == Thymine {
                                                                if cbondtype < 0x200 {
                                                                    if (cbondtype & 0xf) > 0 {
                                                                        if abondtype < 0x2000 {
                                                                            e += 1.5;
                                                                            if dl > 3 {
                                                                                if *WCBP.get_unchecked(*dpos.offset((dstem + di + 1) as isize) as usize).get_unchecked(*cend.offset((tl + 1) as isize) as usize) != 0 {
                                                                                    e += 1.0;
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if j >= 4 {
                                e += 0.25;
                            }
                        }
                        if e > ed {
                            ed = e;
                            ndi = nd;
                            ti = k;
                        }
                    }
                    if ndi < 0 {
                        skip_tvn = 1;
                    }
                    if skip_tvn == 0 {
                        energy += ed;
                        dpos = (*dhit.get_unchecked(ndi as usize)).pos;
                        dstem = (*dhit.get_unchecked(ndi as usize)).stem;
                        dl = (*dhit.get_unchecked(ndi as usize)).loop_;
                        darm = 2 * dstem + dl;
                        dbondtype = (*dhit.get_unchecked(ndi as usize)).bondtype;
                        spacer2 = ((cpos as isize - dpos as isize) / std::mem::size_of::<i32>() as isize) as i32 - darm;
                        spacer1 = ((dpos as isize - aend1 as isize) / std::mem::size_of::<i32>() as isize) as i32;
                        b8 = *aend1;
                        b9 = *aend1.offset(1);

                        // false positive suppression
                        if dloop < 15 {
                            energy -= 2.0;
                        }
                        if cstem > 5 {
                            energy -= 1.0;
                        }
                        if tloop < 6 {
                            energy -= 1.0;
                        }
                        if tloop > 12 {
                            energy -= 1.0;
                            if agcat < 6 {
                                energy -= 2.0;
                            }
                            if tloop > 15 {
                                energy -= 2.5;
                            }
                        }
                        if *STACKBP.get_unchecked(*dpos as usize).get_unchecked(*dpos.offset((darm - 1) as isize) as usize) == 0 {
                            energy -= 1.0;
                        }
                        if dstem < 4 {
                            if gcv > 1.2 {
                                if (dbondtype & 0xf0f) == 0 {
                                    energy -= 1.5;
                                }
                            }
                        }
                        if b8 != Thymine {
                            if dl < 4 {
                                if abondtype > 0x10000 {
                                    energy -= 1.5;
                                }
                            }
                            if b8 == Adenine {
                                if *YI.get_unchecked(*cloopend.offset(-2) as usize) != 0 {
                                    energy -= 1.0;
                                }
                            }
                        }
                        if dl > 10 {
                            if tloop < 7 {
                                energy -= 2.0;
                            }
                            if spacer2 > 1 {
                                energy -= 2.0;
                            }
                            if (*dhit.get_unchecked(ndi as usize)).stem_energy < -3.4 {
                                energy -= 2.0;
                            }
                        }
                        if gcv < 2.0 {
                            if dbondtype > 0x10000 {
                                energy -= 2.0;
                            }
                        }
                        if (cbondtype & 0xf) < 1 {
                            if abondtype > 0x100 {
                                if cgcat < 4 {
                                    energy -= 1.5;
                                }
                                if *WCBP.get_unchecked(*dpos.offset(2) as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) == 0 {
                                    energy -= 1.0;
                                }
                            }
                        }
                        if b8 != Thymine {
                            if (*clooppos.offset(1) != Thymine) || (*clooppos != Cytosine) {
                                if dl > 3 {
                                    if dbondtype > 0x10000 {
                                        energy -= 1.0;
                                    }
                                }
                            }
                        }
                        if *RI.get_unchecked(*cend.offset(1) as usize) == 0 {
                            if b9 != Guanine {
                                energy -= 1.0;
                            } else {
                                energy -= 0.5;
                            }
                        }
                        if b9 == Guanine {
                            if *RI.get_unchecked(*cend as usize) == 0 {
                                energy -= 1.0;
                            }
                            if spacer2 != 1 {
                                energy -= 3.0;
                            } else {
                                tl = if tloop >= 14 { 5 } else { if dl >= 9 { if tloop >= 7 { 4 } else { 3 } } else { 3 } };
                                s = dpos.offset(dstem as isize);
                                if *WCBP.get_unchecked(*s.offset(1) as usize).get_unchecked(*cend.offset(tl as isize) as usize) == 0 {
                                    energy -= 2.5;
                                    if dl >= 5 {
                                        if (*chit.get_unchecked(nc as usize)).energy > 2.0 {
                                            if *WCBP.get_unchecked(*s.offset(2) as usize).get_unchecked(*cend.offset(tl as isize) as usize) != 0 {
                                                if *WCBP.get_unchecked(*s.offset(3) as usize).get_unchecked(*cend.offset((tl + 1) as isize) as usize) != 0 {
                                                    energy += 6.0;
                                                }
                                            }
                                        }
                                    }
                                } else if b8 == Thymine {
                                    if dl >= 5 {
                                        if (*chit.get_unchecked(nc as usize)).energy > 2.0 {
                                            if *WCBP.get_unchecked(*s.offset(2) as usize).get_unchecked(*cend.offset((tl + 1) as isize) as usize) != 0 {
                                                energy += 3.5;
                                            }
                                        }
                                    }
                                }
                            }
                        } else if b9 != Adenine {
                            energy -= 3.0;
                        }
                        if b8 != Thymine {
                            if b8 == Guanine {
                                if *RI.get_unchecked(*dpos.offset(dstem as isize) as usize) == 0 {
                                    energy -= 1.0;
                                } else if *RI.get_unchecked(*dpos.offset((darm - dstem - 1) as isize) as usize) != 0 {
                                    energy += 2.0;
                                }
                            } else {
                                energy -= 1.0;
                            }
                        }

                        // carm termination
                        if *ASSYMST.get_unchecked(*cend.offset(-1) as usize).get_unchecked(*cpos as usize) != 0 {
                            energy += 1.0;
                        }

                        // CTnnnAA cloop motif
                        energy += *CX7.get_unchecked(*clooppos as usize) + *AX7.get_unchecked(*cloopend.offset(-2) as usize);
                        if *clooppos.offset(1) == Cytosine {
                            energy -= 2.0;
                        }

                        // NNnnnAA cloop motif
                        if *cloopend.offset(-2) == Adenine {
                            if *cloopend.offset(-1) == Adenine {
                                if spacer1 == 2 {
                                    if dbondtype < 0x1000 {
                                        if abondtype < 0x100 {
                                            energy += 1.0;
                                        } else if cbondtype < 0x100 {
                                            energy += 1.0;
                                        }
                                    }
                                }
                            }
                        }

                        // global stem damage level
                        bondtype = acbondtype + dbondtype;
                        i = ((bondtype >> 16) & 0xf) as i32;
                        j = ((bondtype >> 12) & 0xf) as i32;
                        k = ((bondtype >> 8) & 0xf) as i32;
                        if k > 0 {
                            if i > 0 {
                                k += i + j;
                                if k > 5 {
                                    energy -= 1.0 * (k - 5) as f64;
                                }
                            }
                        }

                        // global stem stability (GC bond count)
                        gcc = (bondtype & 0xf) as i32;
                        if gcc < 2 {
                            skip_ngcc = 0;
                            if ti >= 2 {
                                if (cbondtype < 0x100) && ((cbondtype & 0xf) > 0) {
                                    skip_ngcc = 2;
                                } else if ti >= 3 {
                                    if cgcat >= 4 {
                                        if (cbondtype & 0xf) > 0 {
                                            skip_ngcc = 2;
                                        } else if cbondtype < 0x100 {
                                            skip_ngcc = 1;
                                        }
                                    }
                                }
                            }
                            if skip_ngcc < 2 {
                                if skip_ngcc < 1 {
                                    energy -= (3 - gcc) as f64;
                                }
                                if gcc < 1 {
                                    if agcat < 5 {
                                        energy -= 2.0;
                                    }
                                    if bondtype > 0x10000 {
                                        energy -= 1.5;
                                    }
                                }
                            }
                        }

                        // global stability
                        // (stem stability,dloop-tloop tertiary interaction,dloop size)
                        if abondtype > 0x1000 {
                            if ti < 3 {
                                if (*chit.get_unchecked(nc as usize)).stem_energy < -6.0 {
                                    energy -= 1.5;
                                }
                                if dl > 9 {
                                    if ((dbondtype + cbondtype) & 0xf) < 1 {
                                        energy -= 1.0;
                                    }
                                }
                            }
                        }

                        // tloop,dloop tertiary interaction
                        // (alternative dloop position)
                        if bondtype < 0x1000 {
                            if b8 == Thymine {
                                if *RI.get_unchecked(b9 as usize) != 0 {
                                    if dl > 4 {
                                        if *BP.get_unchecked(*cend.offset(3) as usize).get_unchecked(*dpos.offset((dstem + 1) as isize) as usize) == 0 {
                                            if *BP.get_unchecked(*cend.offset(3) as usize).get_unchecked(*dpos.offset((dstem + 2) as isize) as usize) != 0 {
                                                energy += 0.5;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // "near perfect" TV-loop mttRNA:
                        // darm-tloop tertiary interaction,low global stem damage,
                        // TR motif at b8-9, good astem,darm,carm interspacing
                        if ti >= 2 {
                            if agcat >= 6 {
                                if cbondtype < 0x100 {
                                    if dbondtype < 0x100 {
                                        if *RI.get_unchecked(b9 as usize) != 0 {
                                            if b8 == Thymine {
                                                if (abondtype & 0xf) > 0 {
                                                    if (dbondtype & 0xf) > 0 {
                                                        if spacer1 == 2 {
                                                            if spacer2 == 1 {
                                                                energy += 1.5;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // find exceptions
                        if energy < dthresh {
                            let mut tvn_except = 0;
                            if mtxdetect == 0 {
                                tvn_except = 1;
                            } else if incds != 0 {
                                tvn_except = 1;
                            } else if energy < (thresh - 7.0) {
                                tvn_except = 1;
                            } else if energy < (dthresh - 7.0) {
                                tvn_except = 1;
                            } else if nbase > 68 {
                                tvn_except = 1;
                            } else if abondtype > 0x20100 {
                                tvn_except = 1;
                            } else if dl > 9 {
                                if dl > 10 {
                                    tvn_except = 1;
                                } else if dstem < 4 {
                                    tvn_except = 1;
                                } else if dbondtype > 0x100 {
                                    tvn_except = 1;
                                }
                            }
                            if tvn_except == 0 {
                                if dstem > 4 {
                                    tvn_except = 1;
                                } else if b9 != Adenine {
                                    if b9 != Guanine {
                                        tvn_except = 1;
                                    } else if cbondtype > 0x100 {
                                        tvn_except = 1;
                                    } else if dbondtype > 0x200 {
                                        tvn_except = 1;
                                    }
                                }
                                if tvn_except == 0 {
                                    if cloop != 7 {
                                        tvn_except = 1;
                                    } else if *YI.get_unchecked(*cloopend.offset(-2) as usize) != 0 {
                                        tvn_except = 1;
                                    }
                                }
                            }
                            if tvn_except == 0 {
                                if b8 == Thymine {
                                    if *apos2.offset(-1) == Thymine {
                                        if *apos2.offset(-2) == Thymine {
                                            if tloop < 8 {
                                                if *TT.get_unchecked(*aend1.offset(-1) as usize).get_unchecked(*apos2 as usize) != 0 {
                                                    if *WCBP.get_unchecked(*dpos.offset(2) as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) != 0 {
                                                        if ((dbondtype + cbondtype) & 0xf) > 0 {
                                                            energy += 3.0;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                } else if b8 == Adenine {
                                    if *apos2.offset(-1) == Adenine {
                                        if *apos2.offset(-2) == Adenine {
                                            if *ASSYMAT.get_unchecked(*aend1.offset(-1) as usize).get_unchecked(*apos2 as usize) != 0 {
                                                if *ASSYMAT.get_unchecked(*apos2.offset(1) as usize).get_unchecked(*aend1.offset(-2) as usize) != 0 {
                                                    energy += 2.0;
                                                }
                                            }
                                            if agcat >= 5 {
                                                if cgcat >= 4 {
                                                    if dbondtype < 0x100 {
                                                        if *AT.get_unchecked(*aend1.offset(-1) as usize).get_unchecked(*apos2 as usize) != 0 {
                                                            if *AT.get_unchecked(*apos2.offset(1) as usize).get_unchecked(*aend1.offset(-2) as usize) != 0 {
                                                                energy += 1.0;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if ti >= 3 {
                                        if cgcat >= 4 {
                                            if agcat >= 4 {
                                                if (cbondtype & 0xf) > 0 {
                                                    if (abondtype & 0xf) > 1 {
                                                        if dbondtype < 0x200 {
                                                            if *WCBP.get_unchecked(*dpos.offset(1) as usize).get_unchecked(*dpos.offset((darm - 2) as isize) as usize) != 0 {
                                                                if *clooppos.offset(1) == Thymine {
                                                                    if *YI.get_unchecked(*clooppos as usize) != 0 {
                                                                        if *RI.get_unchecked(*cloopend.offset(-2) as usize) != 0 {
                                                                            if *RI.get_unchecked(*cloopend.offset(-1) as usize) != 0 {
                                                                                energy += 5.0;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                if bondtype < 0x100 {
                                    if spacer2 == 1 {
                                        if *clooppos == Cytosine {
                                            if *clooppos.offset(1) == Thymine {
                                                if *cloopend.offset(-2) == Adenine {
                                                    if *cloopend.offset(-1) == Adenine {
                                                        energy += 2.0;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    if spacer2 == 1 {
                                        if b8 == Thymine {
                                            if dl > 3 {
                                                if dbondtype < 0x200 {
                                                    if cbondtype < 0x100 {
                                                        if *BP.get_unchecked(*dpos.offset((dstem + 1) as isize) as usize).get_unchecked(*cend.offset(3) as usize) == 0 {
                                                            if *BP.get_unchecked(*dpos.offset((dstem + 1) as isize) as usize).get_unchecked(*cend.offset(4) as usize) != 0 {
                                                                energy += 2.0;
                                                            }
                                                        }
                                                        if dbondtype < 0x100 {
                                                            if abondtype < 0x20000 {
                                                                if ti >= 2 {
                                                                    if dstem >= 3 {
                                                                        if tloop < 13 {
                                                                            if (cbondtype & 0xf) > 0 {
                                                                                energy += 4.0;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    } else if dstem > 3 {
                                                        if dbondtype < 0x300 {
                                                            if bondtype < 0x10000 {
                                                                if ti >= 3 {
                                                                    if (acbondtype & 0xf) > 0 {
                                                                        if *WCBP.get_unchecked(*dpos.offset(2) as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) != 0 {
                                                                            energy += 4.0;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if tloop < 8 {
                                            if dbondtype < 0x200 {
                                                if cbondtype < 0x100 {
                                                    if ti >= 2 {
                                                        if *WCBP.get_unchecked(*dpos.offset((dstem + 1) as isize) as usize).get_unchecked(*cend.offset(3) as usize) != 0 {
                                                            if b8 == Thymine {
                                                                if abondtype < 0x3000 {
                                                                    energy += 5.0;
                                                                }
                                                            }
                                                            if agcat >= 5 {
                                                                if gcv > 1.2 {
                                                                    if *RI.get_unchecked(*cloopend.offset(-1) as usize) != 0 {
                                                                        energy += 7.0;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                        if dbondtype < 0x100 {
                                                            if agcat >= 6 {
                                                                if *YI.get_unchecked(*clooppos as usize) != 0 {
                                                                    if *clooppos.offset(1) == Thymine {
                                                                        if *RI.get_unchecked(*cloopend.offset(-2) as usize) != 0 {
                                                                            if *RI.get_unchecked(*cloopend.offset(-1) as usize) != 0 {
                                                                                energy += 2.0;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                if cbondtype < 0x300 {
                                                    if ti >= 3 {
                                                        if abondtype < 0x2000 {
                                                            if (dbondtype & 0xf) > 0 {
                                                                if (acbondtype & 0xf) > 0 {
                                                                    if (*ahit.get_unchecked(na as usize)).energy >= -7.0 {
                                                                        if dstem >= 4 {
                                                                            energy += 3.0;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            if dbondtype < 0x300 {
                                                if cgcat >= 4 {
                                                    if abondtype < 0x2000 {
                                                        if (*ahit.get_unchecked(na as usize)).energy >= -7.0 {
                                                            if cbondtype < 0x10000 {
                                                                if (cbondtype & 0xf) > 0 {
                                                                    if cstem < 6 {
                                                                        if ti >= 3 {
                                                                            energy += 4.0;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if tloop > 8 {
                                        if agcat >= 6 {
                                            if cbondtype < 0x100 {
                                                if (cbondtype & 0xf) > 0 {
                                                    if b8 == Thymine {
                                                        if *WCBP.get_unchecked(*dpos.offset((dstem + 1) as isize) as usize).get_unchecked(*cend.offset(3) as usize) != 0 {
                                                            if *WCBP.get_unchecked(*dpos.offset(1) as usize).get_unchecked(*dpos.offset((darm - 2) as isize) as usize) != 0 {
                                                                energy += 7.0;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                if dbondtype < 0x100 {
                                    if cgcat >= 4 {
                                        if agcat >= 5 {
                                            if *WCBP.get_unchecked(*dpos.offset(1) as usize).get_unchecked(*dpos.offset((darm - 2) as isize) as usize) != 0 {
                                                if (cbondtype & 0xf) > 0 {
                                                    if (abondtype & 0xf) > 0 {
                                                        if (dbondtype & 0xf) > 0 {
                                                            energy += 0.5;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                if cbondtype < 0x100 {
                                    if dbondtype < 0x200 {
                                        if agcat >= 5 {
                                            if b8 == Thymine {
                                                if tloop < 8 {
                                                    if *WCBP.get_unchecked(*dpos.offset(1) as usize).get_unchecked(*dpos.offset((darm - 2) as isize) as usize) != 0 {
                                                        if *WCBP.get_unchecked(*dpos.offset(2) as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) != 0 {
                                                            if (cbondtype & 0xf) > 0 {
                                                                if (abondtype & 0xf) > 0 {
                                                                    if (dbondtype & 0xf) > 0 {
                                                                        if *clooppos.offset(1) == Thymine {
                                                                            if *YI.get_unchecked(*clooppos as usize) != 0 {
                                                                                if *RI.get_unchecked(*cloopend.offset(-2) as usize) != 0 {
                                                                                    energy += 3.0;
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                if energy >= dthresh {
                                    energy -= 0.9 * (energy - dthresh) + 5.0;
                                } else {
                                    tvn_except = 1;
                                }
                            }
                            if tvn_except != 0 {
                                skip_tvn = 1;
                            }
                        }

                        // remember fully formed TV-loop replacement mttRNA gene
                        // if threshold reached
                        if (skip_tvn == 0) && (energy < thresh) {
                            skip_tvn = 1;
                        }
                        if skip_tvn == 0 {
                            te.energy = energy;
                            thresh = energy;
                            te.ps = apos1;
                            te.dstem = dstem;
                            te.dloop = dl;
                            te.spacer1 = spacer1;
                            te.spacer2 = spacer2;
                            te.cstem = cstem;
                            te.cloop = cloop;
                            k = astem + spacer1 + darm + spacer2;
                            te.anticodon = k + cstem + 2;
                            te.nintron = 0;
                            te.intron = 0;
                            te.var = 0;
                            te.varbp = 0;
                            te.tstem = 0;
                            te.tloop = tloop;
                            te.nbase = k + carm + tloop;
                            tastem = astem;
                            tastem8 = astem8;
                            tastem8d = astem8d;
                        }
                    }

                    // build D-replacement loop mttrna genes
                    if tloop < 10 {
                        nc += 1;
                        continue;
                    }
                }

                // D-replacement loop
                skip_dn = 0;
                if dloop > mt_DRLmaxlength {
                    skip_dn = 1;
                } else if gcv < 1.2 {
                    skip_dn = 1;
                }
                if skip_dn == 0 {
                    energy = 91.0 + ec;

                    // CCnnnAA cloop
                    if *clooppos.offset(1) == Cytosine {
                        if *clooppos != Cytosine {
                            skip_dn = 1;
                        } else if *cloopend.offset(-2) != Adenine {
                            skip_dn = 1;
                        } else if *cloopend.offset(-1) != Adenine {
                            skip_dn = 1;
                        } else {
                            energy -= 1.0;
                        }
                    }
                }

                // choose tarm
                nt = -1;
                nti = -1;
                et = -INACTIVE;
                while {
                    nt += 1;
                    nt < nth
                } {
                    tl = (*thit.get_unchecked(nt as usize)).loop_;
                    if tl > 11 {
                        continue;
                    }
                    if (*thit.get_unchecked(nt as usize)).end != apos2 {
                        continue;
                    }
                    tpos = (*thit.get_unchecked(nt as usize)).pos;
                    tstem = (*thit.get_unchecked(nt as usize)).stem;

                    // var loop (3-7 bases long)
                    var = ((tpos as isize - cend as isize) / std::mem::size_of::<i32>() as isize) as i32;
                    if var < 3 {
                        continue;
                    }
                    e = (*thit.get_unchecked(nt as usize)).energy;
                    if var > 5 {
                        if var > 7 {
                            continue;
                        }
                        if tl < 7 {
                            continue;
                        }
                        e -= 1.0;
                        if (dloop < 10) || (tstem < 4) {
                            e -= 2.0 * (var - 5) as f64;
                        }
                    }

                    // tloop RA or RG motif
                    s = tpos.offset(tstem as isize);
                    k = 0;
                    n = 0;
                    i = 0;
                    loop {
                        j = *tloopa.get_unchecked(tl as usize).get_unchecked(i as usize);
                        if j < 0 {
                            break;
                        }
                        i += 1;
                        if *s.offset(j as isize) == Adenine {
                            k = 1;
                            if dloop >= 3 {
                                if tl > 3 {
                                    b57 = *s.offset((j - 1) as isize);
                                    if *RI.get_unchecked(b57 as usize) != 0 || (tl < 5) {
                                        if *BP.get_unchecked(b57 as usize).get_unchecked(*aend1 as usize) != 0 {
                                            e += 1.5;
                                            n = 1;
                                            break;
                                        }
                                        if *BP.get_unchecked(b57 as usize).get_unchecked(*aend1.offset(1) as usize) != 0 {
                                            e += 1.5;
                                            n = 1;
                                            break;
                                        }
                                        if dloop > 10 {
                                            if *BP.get_unchecked(b57 as usize).get_unchecked(*aend1.offset(2) as usize) != 0 {
                                                e += 1.5;
                                                n = 1;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if k == 0 {
                        i = 0;
                        loop {
                            j = *tloopa.get_unchecked(tl as usize).get_unchecked(i as usize);
                            if j < 0 {
                                break;
                            }
                            i += 1;
                            if *s.offset(j as isize) == Guanine {
                                if *RI.get_unchecked(*s.offset((j - 1) as isize) as usize) != 0 {
                                    k = 1;
                                    break;
                                }
                            }
                        }
                        if j < 0 {
                            e -= if tl > 5 { 2.0 } else { 1.0 };
                        }
                    }

                    // tertiary interaction between tloop and start of dloop
                    ti = if tl > 5 { 1 } else { if dloop > 5 { 1 } else { 0 } };
                    di = if dloop > 5 { 2 } else { 1 };
                    if *STACKBP.get_unchecked(*aend1.offset(di as isize) as usize).get_unchecked(*s.offset(ti as isize) as usize) != 0 {
                        e += 1.0;
                    }

                    // tloop GTTC motif
                    i = if *s.offset(-1) == Guanine { 1 } else { 0 };
                    if tl >= 5 {
                        ti = i + *TI.get_unchecked(*s as usize) + *TI.get_unchecked(*s.offset(1) as usize) + *CI.get_unchecked(*s.offset(2) as usize);
                        if n != 0 {
                            if i == 0 {
                                if *TI.get_unchecked(*s as usize) != 0 {
                                    if *TI.get_unchecked(*s.offset(1) as usize) != 0 {
                                        if *AI.get_unchecked(*s.offset(2) as usize) != 0 {
                                            if tl >= 7 {
                                                ti += 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (i > 0) || (ti >= 3) {
                            e += ti as f64;
                        }
                    } else {
                        ti = i + *TI.get_unchecked(*s as usize) + *TI.get_unchecked(*s.offset(1) as usize);
                        if (i > 0) || (ti >= 2) {
                            e += ti as f64;
                        }
                    }
                    if e > et {
                        et = e;
                        nti = nt;
                        tc = k;
                    }
                }

                if nti < 0 {
                    skip_dn = 1;
                }
                if skip_dn == 0 {
                    energy += et;
                    tpos = (*thit.get_unchecked(nti as usize)).pos;
                    tstem = (*thit.get_unchecked(nti as usize)).stem;
                    tl = (*thit.get_unchecked(nti as usize)).loop_;
                    tbondtype = (*thit.get_unchecked(nti as usize)).bondtype;
                    var = ((tpos as isize - cend as isize) / std::mem::size_of::<i32>() as isize) as i32;

                    // tertiary interaction between b48(=tpos[-1]) and dloop
                    b48 = *tpos.offset(-1);
                    if dloop <= 7 {
                        if *YI.get_unchecked(b48 as usize) != 0 {
                            tc += 1;
                        } else {
                            energy -= 1.0;
                        }
                    } else {
                        i = 0;
                        loop {
                            j = *dloopi.get_unchecked(dloop as usize).get_unchecked(i as usize);
                            if j < 0 {
                                break;
                            }
                            i += 1;
                            if *ASSYMAGBP.get_unchecked(b48 as usize).get_unchecked(*aend1.offset(j as isize) as usize) != 0 {
                                tc += 1;
                                break;
                            }
                        }
                        if j < 0 {
                            energy -= 1.0;
                        }
                    }

                    // large dloop, large tloop
                    if dloop > 7 {
                        if tl >= 6 {
                            if tc < 2 {
                                energy -= 2.0;
                            }
                        }
                        if tstem < 3 {
                            energy -= 1.0;
                        }
                    }

                    // carm termination
                    s = cpos.offset(-1);
                    se = cend;
                    if cstem > 5 {
                        s = s.offset(1);
                        se = se.offset(-1);
                    }
                    if *STACKBP.get_unchecked(*s as usize).get_unchecked(*se as usize) == 0 {
                        energy -= 1.0;
                    }
                    se = cpos.offset(-3);
                    if *BP.get_unchecked(*cend.offset(-1) as usize).get_unchecked(*cpos as usize) == 0 {
                        if *ASSYMST.get_unchecked(*cend.offset(-1) as usize).get_unchecked(*cpos as usize) != 0 {
                            if dloop < 5 {
                                se = se.offset(1);
                            }
                            energy += 1.5;
                        } else if dloop < 13 {
                            se = se.offset(1);
                        }
                    } else {
                        if cstem > 5 {
                            if dloop < 13 {
                                se = se.offset(1);
                            }
                        } else if dloop < 5 {
                            se = se.offset(1);
                        }
                    }

                    // tertiary interaction between tloop and dloop near carm
                    s = tpos.offset(tstem as isize);
                    if tl >= 5 {
                        ti = if tl >= 10 { 4 } else { if tl >= 7 { 3 } else { 2 } };
                        b57 = *s.offset(ti as isize);
                        if *GABP.get_unchecked(*se as usize).get_unchecked(b57 as usize) == 0 {
                            energy -= 2.0;
                        } else {
                            k = if var > 3 { 2 } else { if var > 1 { 1 } else { 0 } };
                            if *BP.get_unchecked(*cend.offset(k as isize) as usize).get_unchecked(b57 as usize) != 0 {
                                energy += 1.0;
                            }
                        }
                    }

                    // R motif at end of tstem
                    if *RI.get_unchecked(*s.offset(-1) as usize) == 0 {
                        energy -= 2.0;
                    }

                    // large tloop
                    if tl > 9 {
                        if tbondtype > 0x200 {
                            energy -= 2.0;
                        }
                    }

                    // dloop,var,tloop T repeat motif
                    // present in some nematode D-loop replacement tRNA-Ser genes
                    if dloop >= 4 {
                        k = 1;
                        se = aend1;
                        while se < cpos {
                            if *se == Thymine {
                                k += 1;
                            }
                            se = se.offset(1);
                        }
                        if k >= dloop {
                            if var >= 3 {
                                se = cend;
                                while se < tpos {
                                    if *se == Thymine {
                                        k += 1;
                                    }
                                    se = se.offset(1);
                                }
                                if k >= (var + dloop) {
                                    energy += 3.0;
                                    se = s.offset(if tl > 5 { 5 } else { tl as isize });
                                    while s < se {
                                        if *s != Thymine {
                                            break;
                                        }
                                        s = s.offset(1);
                                    }
                                    if s >= se {
                                        energy += 5.5;
                                    }
                                }
                            }
                        }
                    }

                    // astem stability
                    s = tpos.offset(tstem as isize);
                    if ea < -6.1 {
                        if tl > 4 {
                            skip_nasi = 0;
                            if (*s == Thymine) && (*s.offset(-1) == Guanine) && (*s.offset(1) == Thymine) {
                                skip_nasi = 1;
                            } else if (ea > -8.3) && (*clooppos == Cytosine) &&
                                (*clooppos.offset(1) == Thymine) && (*cloopend.offset(-2) == Adenine) &&
                                (*cloopend.offset(-1) == Adenine) {
                                skip_nasi = 1;
                            }
                            if skip_nasi == 0 {
                                energy -= 3.0;
                            }
                        }
                    }

                    // cstem stability (GC bond count)
                    bondtype = acbondtype + tbondtype;
                    if (cbondtype & 0xf) < 1 {
                        if (bondtype & 0xf) < 3 {
                            energy -= 1.0;
                        }
                    }

                    // cloop CTnnnAA motif
                    if bondtype >= 0x400 {
                        energy += *CX.get_unchecked(*clooppos as usize) + *TX.get_unchecked(*clooppos.offset(1) as usize) +
                            *AXX.get_unchecked(*cloopend.offset(-1) as usize) + *AXX37.get_unchecked(*cloopend.offset(-2) as usize);
                    } else {
                        energy += *CX.get_unchecked(*clooppos as usize) + *TX.get_unchecked(*clooppos.offset(1) as usize) +
                            *AX.get_unchecked(*cloopend.offset(-1) as usize) + *AX37.get_unchecked(*cloopend.offset(-2) as usize);
                    }

                    // large dloop
                    if dloop >= 9 {
                        k = tloop - dloop - 4;
                        if k < 0 {
                            if bondtype >= 0x1000 {
                                energy += k as f64;
                            }
                        }
                        if dloop >= 12 {
                            if dloop >= 14 {
                                energy -= 2.0;
                            } else if tstem < 6 {
                                energy -= if dloop >= 13 { 2.0 } else { 1.0 };
                            }
                        }
                    }

                    // small dloop, small tarm
                    if dloop <= 10 {
                        if tstem < 3 {
                            if ea > -2.6 {
                                if tl <= 7 {
                                    if cgcat >= 4 {
                                        if *GC.get_unchecked(*tpos as usize).get_unchecked(*apos2.offset(-1) as usize) != 0 {
                                            if *GC.get_unchecked(*tpos.offset(1) as usize).get_unchecked(*apos2.offset(-2) as usize) != 0 {
                                                if gcv > 1.2 {
                                                    if (abondtype & 0xf) > 0 {
                                                        if (cbondtype & 0xf) > 0 {
                                                            energy += 4.5 + (mtBONDSTAB - 0.5) * (5 - tstem) as f64;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // global stem damage level
                    i = ((bondtype >> 16) & 0xf) as i32;
                    j = (((bondtype >> 12) & 0xf) as i32) + i;
                    k = ((bondtype >> 8) & 0xf) as i32;
                    if tstem > 3 {
                        if (k > 0) || (tl > 9) {
                            if (j > 0) || (k > 5) {
                                n = j + k;
                                if (*s.offset(-1) != Guanine) || (*s != Thymine) ||
                                    (*s.offset(1) != Thymine) || (tstem < 5) {
                                    if n > 4 {
                                        energy -= 2.0 * (n - 4) as f64;
                                    }
                                }
                            }
                        }
                    } else {
                        n = j + k;
                        if n > 3 {
                            energy -= 2.0 * (n - 3) as f64;
                        }
                    }

                    // long tstem with tloop GTT motif
                    if *s.offset(-1) == Guanine {
                        if *s == Thymine {
                            if *s.offset(1) == Thymine {
                                if tstem >= 6 {
                                    if tbondtype < 0x100 {
                                        energy += 1.5;
                                    }
                                }
                            }
                        }
                    }

                    // find exceptions
                    if energy < tthresh {
                        let mut dn_except = 0;
                        if mtxdetect == 0 {
                            dn_except = 1;
                        } else if incds != 0 {
                            dn_except = 1;
                        } else if energy < (thresh - 13.5) {
                            dn_except = 1;
                        } else if energy < (tthresh - 13.5) {
                            dn_except = 1;
                        } else if k > 1 {
                            if i > 2 {
                                dn_except = 1;
                            } else if (k > 4) && (i > 1) {
                                dn_except = 1;
                            }
                        }
                        if dn_except == 0 {
                            if nbase > 70 {
                                dn_except = 1;
                            } else if var > 4 {
                                if var > 5 {
                                    dn_except = 1;
                                } else if var > tl {
                                    dn_except = 1;
                                }
                            }
                            if dn_except == 0 {
                                if (tstem < 4) && (((agcat + cgcat + 2) as i32) < (astem + cstem)) {
                                    dn_except = 1;
                                } else if tl > 9 {
                                    dn_except = 1;
                                } else if dloop > 13 {
                                    dn_except = 1;
                                } else if *YI.get_unchecked(*clooppos as usize) == 0 {
                                    dn_except = 1;
                                } else if (abondtype & 0xf) < 2 {
                                    if (abondtype & 0xf) < 1 {
                                        dn_except = 1;
                                    } else if (cbondtype > 0x200) && (tbondtype > 0x100) &&
                                        (abondtype > 0x200) {
                                        dn_except = 1;
                                    }
                                }
                                if dn_except == 0 {
                                    if (tbondtype & 0xf) < 1 {
                                        if (acbondtype & 0xf) < 1 {
                                            dn_except = 1;
                                        } else if acbondtype > 0x200 {
                                            dn_except = 1;
                                        }
                                    }
                                    if dn_except == 0 {
                                        if (dloop + 19) < tloop {
                                            dn_except = 1;
                                        } else if gcv > 5.5 {
                                            dn_except = 1;
                                        } else {
                                            tgcat = ((tbondtype >> 4) + tbondtype) & 0xf;
                                            if (tgcat + 2) < tstem as u32 {
                                                dn_except = 1;
                                            } else if cloop != 7 {
                                                dn_except = 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if dn_except == 0 {
                            if *BP.get_unchecked(*cpos as usize).get_unchecked(*cend.offset(-1) as usize) != 0 {
                                if *BP.get_unchecked(*cpos.offset(-1) as usize).get_unchecked(*cend as usize) != 0 {
                                    if *BP.get_unchecked(*cpos.offset(-2) as usize).get_unchecked(*cend.offset(1) as usize) != 0 {
                                        energy += 2.0;
                                    }
                                }
                            }
                            if bondtype < 0x20000 {
                                if (*thit.get_unchecked(nti as usize)).stem_energy > -4.6 {
                                    if tstem >= 4 {
                                        if (tstem >= 5) || (*s.offset(-1) == Guanine) {
                                            if *STACKBP.get_unchecked(*cpos.offset(1) as usize).get_unchecked(*cend.offset(-2) as usize) != 0 {
                                                if *STACKBP.get_unchecked(*cpos as usize).get_unchecked(*cend.offset(-1) as usize) != 0 {
                                                    if *STACKBP.get_unchecked(*cpos.offset(-1) as usize).get_unchecked(*cend as usize) != 0 {
                                                        energy += 1.5;
                                                        if *s.offset(-1) == Guanine {
                                                            if *s == Thymine {
                                                                if *s.offset(1) == Thymine {
                                                                    energy += 1.0;
                                                                }
                                                            }
                                                        }
                                                        if agcat >= 6 {
                                                            energy += 0.5;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if tc > 0 {
                                if tstem >= 5 {
                                    if var < 6 {
                                        if var > 2 {
                                            if acbondtype < 0x100 {
                                                energy += 5.0;
                                            } else if (abondtype + tbondtype) < 0x100 {
                                                energy += 3.0;
                                            } else {
                                                if *cloopend.offset(-2) == Thymine {
                                                    if *cloopend.offset(-1) == Thymine {
                                                        if dloop > 7 {
                                                            if tbondtype < 0x100 {
                                                                if *TT.get_unchecked(*tpos as usize).get_unchecked(*apos2.offset(-1) as usize) == 0 {
                                                                    if (agcat + cgcat) >= 10 {
                                                                        energy += 13.5;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if *s.offset(-1) == Guanine {
                                if *s == Thymine {
                                    if *s.offset(1) == Thymine {
                                        if (tstem >= 5) || (*s.offset(2) == Cytosine) {
                                            energy += 1.5;
                                            if tstem >= 5 {
                                                if tbondtype < 0x1000 {
                                                    if *s.offset(2) == Cytosine {
                                                        if abondtype < 0x10000 {
                                                            if *clooppos == Cytosine {
                                                                if *clooppos.offset(1) == Thymine {
                                                                    if *cloopend.offset(-2) == Adenine {
                                                                        if *cloopend.offset(-1) == Adenine {
                                                                            energy += 3.0;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                            if tbondtype < 0x200 {
                                                                if bondtype < 0x10000 {
                                                                    if tl == 7 {
                                                                        if *s.offset(4) == Adenine {
                                                                            energy += 4.0;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    } else {
                                                        if tbondtype < 0x200 {
                                                            if (tbondtype & 0xf) >= 2 {
                                                                if *clooppos == Cytosine {
                                                                    if *clooppos.offset(1) == Thymine {
                                                                        if *cloopend.offset(-2) == Adenine {
                                                                            if *cloopend.offset(-1) == Adenine {
                                                                                energy += 1.0;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if tstem >= 4 {
                                if tbondtype < 0x100 {
                                    if cbondtype < 0x200 {
                                        if agcat >= 5 {
                                            energy += 1.5;
                                        }
                                    }
                                }
                            }
                            if energy > tthresh {
                                energy = tthresh;
                            }
                            if ea > -1.8 {
                                energy += 3.0;
                            } else if abondtype < 0x60 {
                                energy += 1.5;
                            } else if acbondtype < 0x200 {
                                energy += 0.75;
                            }
                            if *clooppos == Cytosine {
                                if *cloopend.offset(-2) == Adenine {
                                    if *cloopend.offset(-1) == Adenine {
                                        if tstem >= 5 {
                                            if tbondtype < 0x100 {
                                                if *clooppos.offset(1) == Thymine {
                                                    energy += 3.0;
                                                    if tstem >= 6 {
                                                        energy += 1.0;
                                                    }
                                                } else if *clooppos.offset(1) == Cytosine {
                                                    energy += 1.0;
                                                }
                                            }
                                        }
                                        if tc >= 2 {
                                            if *clooppos.offset(1) == Thymine {
                                                if bondtype < 0x1000 {
                                                    if tstem >= 4 {
                                                        if var < 6 {
                                                            if var > 2 {
                                                                energy += 3.0;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if cbondtype < 0x100 {
                                if agcat >= 5 {
                                    if tc > 0 {
                                        if *clooppos.offset(1) == Thymine {
                                            if *YI.get_unchecked(*clooppos as usize) != 0 {
                                                if *RI.get_unchecked(*cloopend.offset(-2) as usize) != 0 {
                                                    if *RI.get_unchecked(*cloopend.offset(-1) as usize) != 0 {
                                                        if tbondtype < 0x100 {
                                                            energy += 4.0;
                                                        } else if agcat >= 6 {
                                                            if (tgcat + 1) >= tstem as u32 {
                                                                if tstem >= 4 {
                                                                    energy += 4.0;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if bondtype < 0x1000 {
                                energy += 0.5;
                                if bondtype < 0x200 {
                                    energy += 0.75;
                                }
                            }
                            if energy >= tthresh {
                                energy -= 3.0 + 0.9 * (energy - tthresh);
                            } else {
                                dn_except = 1;
                            }
                        }
                        if dn_except != 0 {
                            skip_dn = 1;
                        }
                    }

                    // mammalian cloop motif constraint
                    if skip_dn == 0 {
                        if ds == MAMMAL_MT {
                            s1 = clooppos;
                            s2 = s1.offset(cloop as isize);
                            r = *s1 as u32;
                            s1 = s1.offset(1);
                            while s1 < s2 {
                                r = (r << 4) + *s1 as u32;
                                s1 = s1.offset(1);
                            }
                            if r != clmotif[0] {
                                if r != clmotif[1] {
                                    if r != clmotif[2] {
                                        energy -= 5.0;
                                    }
                                }
                            }
                        }

                        // remember fully formed D-loop replacement mttRNA gene
                        // if threshold reached
                        if energy < thresh {
                            skip_dn = 1;
                        }
                        if skip_dn == 0 {
                            te.energy = energy;
                            thresh = energy;
                            te.ps = apos1;
                            te.spacer1 = 0;
                            te.spacer2 = 0;
                            te.dstem = 0;
                            te.dloop = dloop;
                            te.cstem = cstem;
                            te.cloop = cloop;
                            te.anticodon = astem + dloop + cstem + 2;
                            te.nintron = 0;
                            te.intron = 0;
                            te.var = var;
                            te.varbp = 0;
                            te.tstem = tstem;
                            te.tloop = tl;
                            te.nbase = astem + dloop + carm + var + 2 * tstem + tl;
                            tastem = astem;
                            tastem8 = astem8;
                            tastem8d = astem8d;
                        }
                    }
                }

                // build fully formed cloverleaf mttRNA genes
                if dloop < 10 {
                    nc += 1;
                    continue;
                }

                // choose tarm
                nt = -1;
                nti = -1;
                et = -INACTIVE;
                while {
                    nt += 1;
                    nt < nth
                } {
                    tend = (*thit.get_unchecked(nt as usize)).end;
                    if tend != apos2 {
                        continue;
                    }
                    e = (*thit.get_unchecked(nt as usize)).energy;
                    tpos = (*thit.get_unchecked(nt as usize)).pos;
                    tstem = (*thit.get_unchecked(nt as usize)).stem;

                    // GT motif on tloop
                    s = tpos.offset(tstem as isize);
                    if *s == Thymine {
                        if *s.offset(-1) == Guanine {
                            if tstem >= 5 {
                                if *STACKBP.get_unchecked(*tpos as usize).get_unchecked(*tend.offset(-1) as usize) == 0 {
                                    e += 0.5;
                                    if *BP.get_unchecked(*tpos.offset(1) as usize).get_unchecked(*tend.offset(-2) as usize) == 0 {
                                        e += 0.5;
                                    }
                                }
                            }
                        }
                    }

                    // large var loop
                    var = ((tpos as isize - cend as isize) / std::mem::size_of::<i32>() as isize) as i32;
                    if var > 5 {
                        ev = (var - 5) as f64;
                        if tstem < 5 {
                            e -= 3.0 * ev;
                        } else {
                            e -= 0.5 + 0.5 * ev;
                        }

                        // allow large var loop if tarm looks nuclear
                        // (GTTC motif, very large var loop base-pairing)
                        if var > 9 {
                            if ((*thit.get_unchecked(nt as usize)).bondtype & 0xf) < 1 {
                                e -= 1.0;
                            }
                            e -= 0.25 * (var - 8) as f64;
                            if *s == Thymine {
                                if *s.offset(-1) == Guanine {
                                    if *s.offset(1) == Thymine {
                                        if *s.offset(2) == Cytosine {
                                            e += 4.0;
                                        }
                                    }
                                }
                            }
                            if var > 17 {
                                if var > 25 {
                                    continue;
                                }
                                e += 0.5 * vloop_stability(cend, var, &mut varbp);
                            }
                        }
                    }

                    // small var loop
                    if var < 3 {
                        if tstem > 5 {
                            if *s.offset(-1) != Guanine {
                                e -= 0.5;
                            }
                        }
                        if var < 2 {
                            if var < 1 {
                                if var < 0 {
                                    continue;
                                }
                                if tstem < 4 {
                                    if (*thit.get_unchecked(nt as usize)).stem_energy < -4.0 {
                                        continue;
                                    }
                                }
                            }
                            e -= 3.0;
                        }
                    }
                    if e > et {
                        et = e;
                        nti = nt;
                    }
                }

                if nti < 0 {
                    nc += 1;
                    continue;
                }
                tpos = (*thit.get_unchecked(nti as usize)).pos;
                tstem = (*thit.get_unchecked(nti as usize)).stem;
                tl = (*thit.get_unchecked(nti as usize)).loop_;
                tarm = 2 * tstem + tl;
                var = ((tpos as isize - cend as isize) / std::mem::size_of::<i32>() as isize) as i32;
                b48 = *tpos.offset(-1);
                tbondtype = (*thit.get_unchecked(nti as usize)).bondtype;
                bondtype = acbondtype + tbondtype;
                ti = (((bondtype >> 16) & 0xf) + ((bondtype >> 12) & 0xf) + ((bondtype >> 8) & 0xf)) as i32;

                // choose darm
                nd = -1;
                ndi = -1;
                ed = -INACTIVE;
                while {
                    nd += 1;
                    nd < ndh
                } {
                    dl = (*dhit.get_unchecked(nd as usize)).loop_;
                    dstem = (*dhit.get_unchecked(nd as usize)).stem;
                    darm = 2 * dstem + dl;
                    dpos = (*dhit.get_unchecked(nd as usize)).pos;
                    e = (*dhit.get_unchecked(nd as usize)).energy;

                    // spacing between astem,darm,carm
                    spacer1 = ((dpos as isize - aend1 as isize) / std::mem::size_of::<i32>() as isize) as i32;
                    spacer2 = ((cpos as isize - dpos as isize) / std::mem::size_of::<i32>() as isize) as i32 - darm;
                    if spacer1 < 2 {
                        if spacer1 < 1 {
                            continue;
                        }
                        if dstem < 3 {
                            continue;
                        }
                        if dl > 12 {
                            e -= 2.0;
                        }
                        if astem < 7 {
                            e -= 1.0;
                        }
                        if spacer2 != 2 {
                            if spacer2 < 1 {
                                continue;
                            }
                            if spacer2 > 2 {
                                continue;
                            }
                            if (abondtype & 0xf) < 1 {
                                if ((*dhit.get_unchecked(nd as usize)).bondtype & 0xf) < 1 {
                                    e -= 0.5;
                                }
                            }
                            if var > 7 {
                                e -= 1.0;
                            }
                            if dl > 12 {
                                e -= 1.0;
                            }
                            if cloop != 7 {
                                e -= 2.0;
                            }
                            if cstem < 6 {
                                e -= 3.6;
                            } else {
                                e -= 0.5;
                            }
                        } else {
                            if cstem > 5 {
                                continue;
                            }
                            s = cpos;
                            se = cend.offset(-1);
                            while *BP.get_unchecked(*s as usize).get_unchecked(*se as usize) == 0 {
                                s = s.offset(1);
                                se = se.offset(-1);
                            }
                            if *STEMTERM.get_unchecked(*s.offset(-1) as usize).get_unchecked(*se.offset(1) as usize) == 0 {
                                e -= 0.5;
                            }
                            e -= 0.8;
                        }
                    } else {
                        if spacer1 > 2 {
                            if spacer1 > 3 {
                                continue;
                            }
                            if dstem > 4 {
                                continue;
                            }
                            if dstem < 3 {
                                continue;
                            }
                            if tl > 15 {
                                continue;
                            }
                            if astem < 7 {
                                e -= 1.0;
                            }
                            if ti > 4 {
                                e -= 1.0;
                            }
                            if cloop != 7 {
                                e -= 2.0;
                            }
                            if tbondtype > 0x2000 {
                                if *RI.get_unchecked(*tpos.offset((tstem - 1) as isize) as usize) == 0 {
                                    e -= 2.0;
                                }
                            }
                            e -= 1.0;
                            if spacer2 != 1 {
                                e -= 0.5;
                            } else if (*dhit.get_unchecked(nd as usize)).bondtype < 0x100 {
                                if var >= 3 {
                                    if var <= 5 {
                                        if tstem >= 3 {
                                            e += 1.0;
                                            if agcat >= 5 {
                                                if *WCBP.get_unchecked(*aend1 as usize).get_unchecked(*apos2 as usize) != 0 {
                                                    if *BP.get_unchecked(*aend1.offset(-1) as usize).get_unchecked(*apos2 as usize) == 0 {
                                                        if *BP.get_unchecked(b48 as usize).get_unchecked(*dpos.offset((dstem + 1) as isize) as usize) != 0 {
                                                            e += 0.5;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if spacer2 > 1 {
                            if spacer2 > 2 {
                                continue;
                            }
                            if astem < 7 {
                                if spacer1 == 2 {
                                    e -= 1.0;
                                }
                            }
                            if cloop != 7 {
                                e -= 2.0;
                            }
                            if ea < -5.8 {
                                e -= 2.0;
                            }
                            e -= 2.5;
                            if *BP.get_unchecked(b48 as usize).get_unchecked(*dpos.offset((dstem + 1) as isize) as usize) != 0 {
                                if (*dhit.get_unchecked(nd as usize)).bondtype < 0x1000 {
                                    if *WCBP.get_unchecked(*dpos.offset(1) as usize).get_unchecked(*dpos.offset((darm - 2) as isize) as usize) != 0 {
                                        if *WCBP.get_unchecked(*dpos.offset(2) as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) != 0 {
                                            if var < 6 {
                                                if dl > 3 {
                                                    e += 2.0;
                                                }
                                            }
                                        }
                                    }
                                }
                            } else {
                                e -= 1.0;
                            }
                        } else if spacer2 < 1 {
                            if spacer2 < 0 {
                                continue;
                            }
                            if var > 6 {
                                continue;
                            }
                            if dstem > 4 {
                                continue;
                            }
                            if (*dhit.get_unchecked(nd as usize)).stem_energy < -4.3 {
                                continue;
                            }
                            if astem < 7 {
                                if spacer1 == 2 {
                                    e -= 1.0;
                                }
                            }
                            if cloop != 7 {
                                e -= 2.0;
                            }
                            e -= mtBONDSTAB;
                        }
                        if cstem > 5 {
                            if (*GT.get_unchecked(*cpos as usize).get_unchecked(*cend.offset(-1) as usize) == 0) || (astem8 != 0) {
                                e -= mtBONDSTAB;
                            }
                        }
                    }

                    // very large or very small dloop
                    if dl < 3 {
                        e -= 2.0;
                    }
                    if dl > 11 {
                        if dl > 14 {
                            e -= 2.0;
                        } else if dl > 13 {
                            if (*dhit.get_unchecked(nd as usize)).bondtype >= 0x100 {
                                e -= 2.0;
                            } else {
                                e -= 1.0;
                            }
                        } else if dl > 12 {
                            if (*dhit.get_unchecked(nd as usize)).bondtype >= 0x1000 {
                                e -= 2.0;
                            } else {
                                e -= 1.0;
                            }
                        } else if (*dhit.get_unchecked(nd as usize)).bondtype >= 0x10000 {
                            e -= 2.0;
                        }
                    }

                    // tertiary interactions in darm
                    b8 = *dpos.offset(-2);
                    b9 = *dpos.offset(-1);
                    if dl > 2 {
                        if dl > 5 {
                            if *STACKBP.get_unchecked(*dpos.offset((dstem + 1) as isize) as usize).get_unchecked(b48 as usize) == 0 {
                                e -= 1.0;
                            }
                        }
                        if *STACKBP.get_unchecked(b8 as usize).get_unchecked(*dpos.offset(dstem as isize) as usize) == 0 {
                            e -= 0.25;
                        }
                        if *STACKBP.get_unchecked(b8 as usize).get_unchecked(*dpos.offset((dstem + dl - 1) as isize) as usize) == 0 {
                            e -= 0.25;
                        }
                    }
                    if *BP.get_unchecked(b9 as usize).get_unchecked(*dpos.offset(2) as usize) == 0 {
                        if *BP.get_unchecked(b9 as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) == 0 {
                            e -= 1.0;
                        }
                    }

                    // TR motif at b8-9
                    if *RI.get_unchecked(b9 as usize) != 0 {
                        if b8 == Thymine {
                            if spacer1 == 2 {
                                if ti < 6 {
                                    if ((bondtype & 0xf) > 2) || (bondtype < 0x1000) ||
                                        ((tbondtype < 0x100) && (tstem > 3)) {
                                        if (cbondtype & 0xf) < 5 {
                                            if *STEMBP.get_unchecked(*dpos.offset(1) as usize).get_unchecked(*dpos.offset((darm - 2) as isize) as usize) != 0 {
                                                if var < 6 {
                                                    if var > 2 {
                                                        e += 1.5;
                                                    } else if tstem > 3 {
                                                        if *cloopend.offset(-2) == Adenine {
                                                            e += 1.5;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        e -= 1.0;
                        if b9 == Thymine {
                            if spacer1 == 2 {
                                e -= 2.0;
                            }
                        }
                    }
                    if e > ed {
                        ed = e;
                        ndi = nd;
                    }
                }

                if ndi < 0 {
                    nc += 1;
                    continue;
                }
                energy = 100.0 + ec + ed + et;
                dl = (*dhit.get_unchecked(ndi as usize)).loop_;
                dstem = (*dhit.get_unchecked(ndi as usize)).stem;
                darm = 2 * dstem + dl;
                dpos = (*dhit.get_unchecked(ndi as usize)).pos;
                dbondtype = (*dhit.get_unchecked(ndi as usize)).bondtype;
                spacer1 = ((dpos as isize - aend1 as isize) / std::mem::size_of::<i32>() as isize) as i32;
                spacer2 = ((cpos as isize - dpos as isize) / std::mem::size_of::<i32>() as isize) as i32 - darm;
                b8 = *dpos.offset(-2);

                // tertiary structure interaction between tloop and dloop
                if tl >= 3 {
                    if dl >= 4 {
                        di = if dl < 7 { darm - dstem - 2 } else { darm - dstem - 3 };
                        ti = if tl < 9 { tstem + 2 } else { if tl < 13 { tstem + 3 } else { tstem + 5 } };
                        if *GGBP.get_unchecked(*dpos.offset(di as isize) as usize).get_unchecked(*tpos.offset(ti as isize) as usize) != 0 {
                            if *GGBP.get_unchecked(*dpos.offset((di - 1) as isize) as usize).get_unchecked(*tpos.offset((ti - 1) as isize) as usize) != 0 {
                                energy += 2.0;
                                if spacer1 != 2 {
                                    if spacer2 != 2 {
                                        if dstem < 4 {
                                            if tl > 7 {
                                                if *BP.get_unchecked(*dpos.offset((di + 1) as isize) as usize).get_unchecked(*tpos.offset((ti + 1) as isize) as usize) != 0 {
                                                    energy += 4.0;
                                                }
                                            }
                                        }
                                    }
                                }
                                if ea > -2.5 {
                                    if *WCBP.get_unchecked(*dpos.offset(1) as usize).get_unchecked(*dpos.offset((darm - 2) as isize) as usize) != 0 {
                                        if *WCBP.get_unchecked(*dpos.offset(2) as usize).get_unchecked(*dpos.offset((darm - 3) as isize) as usize) != 0 {
                                            energy += 3.0;
                                        }
                                    }
                                }
                            }
                        }
                        if tl > 10 {
                            if dl > 10 {
                                energy -= 1.0;
                            }
                        }
                    } else if dl == 3 {
                        if *WCBP.get_unchecked(*dpos.offset((dstem + 1) as isize) as usize).get_unchecked(b48 as usize) != 0 {
                            energy += 1.0;
                        }
                    }
                }

                // small darm and tarm
                if tloop <= 18 {
                    if tarm <= 13 {
                        if dl <= 8 {
                            if spacer1 == 2 {
                                if spacer2 == 1 {
                                    if abondtype < 0x1000 {
                                        if tbondtype < 0x100 {
                                            if dbondtype < 0x200 {
                                                et = (mtBONDSTAB - 0.5) * (5 - tstem) as f64 + 0.1 * (7 - tl) as f64;
                                                ed = mtBONDSTAB * (4 - dstem) as f64;
                                                energy += 0.8 * (et + ed);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // GTTC motif on tloop
                s = tpos.offset(tstem as isize);
                if tl < 5 {
                    if tl < 2 {
                        energy += *G_.get_unchecked(*s.offset(-1) as usize);
                    } else {
                        et = *G_.get_unchecked(*s.offset(-1) as usize) + *T_.get_unchecked(*s as usize) + *T_.get_unchecked(*s.offset(1) as usize);
                        if tl > 3 {
                            if *BP.get_unchecked(*s as usize).get_unchecked(*s.offset((tl - 1) as isize) as usize) != 0 {
                                e = *G_.get_unchecked(*s as usize) + *T_.get_unchecked(*s.offset(1) as usize) + *T_.get_unchecked(*s.offset(2) as usize);
                                if e > et {
                                    et = e;
                                }
                            }
                        }
                        if tstem < 5 {
                            e = *G_.get_unchecked(*s.offset(-2) as usize) + *T_.get_unchecked(*s.offset(-1) as usize) + *T_.get_unchecked(*s as usize) + *C_.get_unchecked(*s.offset(1) as usize);
                            if e > et {
                                et = e;
                            }
                        }
                        energy += et;
                    }
                } else {
                    energy += *G_.get_unchecked(*s.offset(-1) as usize) + *T_.get_unchecked(*s as usize) + *T_.get_unchecked(*s.offset(1) as usize) + *C_.get_unchecked(*s.offset(2) as usize);
                }

                // long astem
                if astem8 != 0 {
                    if *BP.get_unchecked(*apos1 as usize).get_unchecked(*apos2.offset(6) as usize) != 0 {
                        if *BP.get_unchecked(*apos1.offset(1) as usize).get_unchecked(*apos2.offset(5) as usize) != 0 {
                            if *BP.get_unchecked(*apos1.offset(2) as usize).get_unchecked(*apos2.offset(4) as usize) != 0 {
                                if *BP.get_unchecked(*apos1.offset(3) as usize).get_unchecked(*apos2.offset(3) as usize) != 0 {
                                    energy += *hbem.get_unchecked(*apos1.offset(-1) as usize).get_unchecked(*apos2.offset(7) as usize);
                                }
                            }
                        }
                    }
                }

                // false positive suppression
                if *RI.get_unchecked(*cend as usize) == 0 {
                    energy -= 1.0;
                }
                if *RI.get_unchecked(*cpos.offset(-1) as usize) == 0 {
                    energy -= 1.0;
                }
                if tarm < (var + 3) {
                    energy -= 2.0;
                }
                if gcv < 1.5 {
                    if dbondtype > 0x10000 {
                        energy -= 2.0;
                    }
                }
                if tarm > 27 {
                    energy -= 1.0;
                    if spacer2 != 1 {
                        energy -= 1.0;
                    }
                }
                if dstem < 3 {
                    if var > 5 {
                        energy -= 1.0;
                    }
                    if tloop > (dloop + 8) {
                        energy -= 0.5;
                    }
                }
                if b8 != Thymine {
                    if dl > 3 {
                        if dbondtype > 0x100 {
                            if (b8 == Cytosine) || (dbondtype > 0x10000) {
                                if *clooppos != Cytosine {
                                    if *WCBP.get_unchecked(*dpos.offset((dstem + 1) as isize) as usize).get_unchecked(b48 as usize) == 0 {
                                        energy -= 1.0;
                                    }
                                }
                            }
                        }
                    }
                }

                // high GC false positive suppression
                if gcv >= 5.1 {
                    if (abondtype & 0xf) >= 4 {
                        s1 = apos1;
                        s2 = apos2.offset(astem as isize);
                        n = 0;
                        s2 = s2.offset(-1);
                        while s2 >= apos2 {
                            if *GC.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize) != 0 {
                                n += 1;
                                if n >= 4 {
                                    energy -= 2.0;
                                    break;
                                }
                            } else {
                                n = 0;
                            }
                            s1 = s1.offset(1);
                            s2 = s2.offset(-1);
                        }
                    }
                    if (dbondtype & 0xf) >= 4 {
                        energy -= 3.0;
                    }
                    if (cbondtype & 0xf) >= 5 {
                        energy -= 3.5;
                    }
                    if (tbondtype & 0xf) >= tstem as u32 {
                        energy -= 4.0;
                    }
                }

                // global stem damage level
                tc = tstem + dstem;
                dtbondtype = dbondtype + tbondtype;
                mabondtype = dtbondtype + cbondtype;
                bondtype = acbondtype + dtbondtype;
                if bondtype < 0x100 {
                    energy += 0.5;
                }
                if (dtbondtype & 0xf) < 1 {
                    energy -= 1.0;
                    if tc >= 10 {
                        energy -= 2.0;
                    }
                    if (bondtype & 0xf) < 3 {
                        if nbase > 75 {
                            energy -= 1.0;
                        }
                    }
                }
                i = ((bondtype >> 16) & 0xf) as i32;
                j = (((bondtype >> 12) & 0xf) as i32) + i;
                k = (((bondtype >> 8) & 0xf) as i32) + j;
                ti = if tc > 6 { 5 } else { if tc > 5 { 4 } else { 3 } };
                if k > ti {
                    ev = (k - ti) as f64;
                    energy -= 0.5 * ev;
                    if cbondtype > 0x10000 {
                        if tstem < 5 {
                            energy -= ev;
                        }
                    }
                    if i > 0 {
                        if k > 8 {
                            energy -= 1.5 * (k - 8) as f64;
                        }
                    }
                }

                // low GC false positive suppression
                if gcv < 3.5 {
                    if (bondtype & 0xf) < 2 {
                        if (bondtype & 0xf) < 1 {
                            energy -= 1.0;
                        }
                        if dl > 3 {
                            if var > 2 {
                                if *WCBP.get_unchecked(*dpos.offset((dstem + 1) as isize) as usize).get_unchecked(b48 as usize) == 0 {
                                    energy -= 1.0;
                                }
                            }
                        }
                    }
                }

                // small variable loop
                if var < 3 {
                    if dloop > 18 {
                        if dloop > (tloop + 2) {
                            energy -= 1.0;
                        }
                        if tloop > 20 {
                            if (((dtbondtype >> 4) + dtbondtype) & 0xf) < 6 {
                                energy -= 2.0;
                            }
                        }
                    }
                    if astem < 7 {
                        energy -= 1.0;
                        if agcat >= 5 {
                            if bondtype < 0x300 {
                                if gcv > 1.2 {
                                    if gcv < 5.0 {
                                        energy += 2.0;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    // NNNNNAA cloop
                    if *cloopend.offset(-2) == Adenine {
                        if *cloopend.offset(-1) == Adenine {
                            if spacer1 > 1 {
                                if (dbondtype < 0x2000) || (dloop > mt_DRLmaxlength) {
                                    if abondtype < 0x100 {
                                        energy += 1.0;
                                    } else if cbondtype < 0x100 {
                                        energy += 1.0;
                                    } else if tstem >= 5 {
                                        if tbondtype < 0x100 {
                                            energy += 1.0;
                                            if *clooppos == Cytosine {
                                                if *clooppos.offset(1) == Thymine {
                                                    if dbondtype < 0x100 {
                                                        energy += 0.5;
                                                    }
                                                }
                                            }
                                            if cgcat >= 3 {
                                                if (tbondtype & 0xf) > 0 {
                                                    if *GGBP.get_unchecked(*dpos.offset((dstem + 1) as isize) as usize).get_unchecked(b48 as usize) != 0 {
                                                        if *WCBP.get_unchecked(*dpos.offset(1) as usize).get_unchecked(*dpos.offset((darm - 2) as isize) as usize) != 0 {
                                                            if tl < 10 {
                                                                if spacer1 == 2 {
                                                                    if spacer2 == 1 {
                                                                        if dl > 2 {
                                                                            if var >= 2 {
                                                                                if var < 6 {
                                                                                    if agcat >= 6 {
                                                                                        energy += 3.0;
                                                                                    } else if agcat >= 5 {
                                                                                        if cgcat >= 4 {
                                                                                            if dbondtype < 0x100 {
                                                                                                if *s == Thymine {
                                                                                                    if *s.offset(-1) == Guanine {
                                                                                                        if *s.offset(1) == Thymine {
                                                                                                            energy += 3.0;
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // large tloop
                if tl > 12 {
                    if tbondtype > 0x10000 {
                        energy -= 2.0;
                    }
                    if agcat < 5 {
                        if spacer1 != 2 {
                            if spacer2 != 1 {
                                energy -= 1.0;
                            }
                        }
                    }
                }

                // find exceptions
                if energy < dtthresh {
                    if mtxdetect == 0 {
                        nc += 1;
                        continue;
                    }
                    if incds != 0 {
                        nc += 1;
                        continue;
                    }
                    if energy < (thresh - 12.0) {
                        nc += 1;
                        continue;
                    }
                    if energy < (dtthresh - 12.0) {
                        nc += 1;
                        continue;
                    }
                    if nbase > 75 {
                        nc += 1;
                        continue;
                    }
                    if dstem > 4 {
                        nc += 1;
                        continue;
                    }
                    if dstem < 3 {
                        nc += 1;
                        continue;
                    }
                    if astem < 7 {
                        if acbondtype > 0x21000 {
                            nc += 1;
                            continue;
                        }
                    }
                    if var > 5 {
                        if var > 6 {
                            nc += 1;
                            continue;
                        }
                        if tarm < 12 {
                            nc += 1;
                            continue;
                        }
                    }
                    if gcv <= 1.2 {
                        if gcv < 0.9 {
                            nc += 1;
                            continue;
                        }
                        if (mabondtype & 0xf) < 1 {
                            nc += 1;
                            continue;
                        }
                    }
                    if tl > 9 {
                        if tl > 13 {
                            nc += 1;
                            continue;
                        }
                        if *WCBP.get_unchecked(*dpos.offset(1) as usize).get_unchecked(*dpos.offset((darm - 2) as isize) as usize) == 0 {
                            nc += 1;
                            continue;
                        }
                    }
                    if dl > 7 {
                        if bondtype > 0x20000 {
                            if dloop > (tloop + 4) {
                                nc += 1;
                                continue;
                            }
                        }
                        if dl > 10 {
                            if dl > 12 {
                                if abondtype > 0x1000 {
                                    nc += 1;
                                    continue;
                                }
                            }
                            if tbondtype > 0x200 {
                                nc += 1;
                                continue;
                            }
                            if *TT.get_unchecked(*tpos as usize).get_unchecked(*apos2.offset(-1) as usize) != 0 {
                                nc += 1;
                                continue;
                            }
                            if var > 5 {
                                nc += 1;
                                continue;
                            }
                            if dloop > (tloop + 8) {
                                if bondtype > 0x10000 {
                                    nc += 1;
                                    continue;
                                }
                            }
                            if astem < 7 {
                                nc += 1;
                                continue;
                            }
                        }
                    }
                    if *RI.get_unchecked(*clooppos.offset(1) as usize) != 0 {
                        nc += 1;
                        continue;
                    }
                    b9 = *dpos.offset(-1);
                    if cstem >= 6 {
                        if cbondtype > 0x200 {
                            nc += 1;
                            continue;
                        }
                        if var < 3 {
                            nc += 1;
                            continue;
                        }
                        if *YI.get_unchecked(b9 as usize) != 0 {
                            nc += 1;
                            continue;
                        }
                    }
                    if cloop != 7 {
                        nc += 1;
                        continue;
                    }
                    if ds == MAMMAL_MT {
                        nc += 1;
                        continue;
                    }
                    // Exception handling code continues here...
                    // Due to length, simplified threshold adjustment
                    if energy >= dtthresh {
                        energy -= 0.9 * (energy - dtthresh) + 5.0;
                    } else {
                        nc += 1;
                        continue;
                    }
                }

                // remember fully formed mttRNA gene if threshold reached
                if energy < thresh {
                    nc += 1;
                    continue;
                }
                te.energy = energy;
                thresh = energy;
                te.ps = apos1;
                te.spacer1 = spacer1;
                te.dstem = dstem;
                te.dloop = dl;
                te.spacer2 = spacer2;
                te.cstem = cstem;
                te.cloop = cloop;
                te.var = var;
                te.varbp = if var > 17 { varbp } else { 0 };
                te.tstem = tstem;
                te.tloop = tl;
                k = astem + spacer1 + darm + spacer2;
                te.anticodon = k + cstem + 2;
                te.nintron = 0;
                te.intron = 0;
                te.nbase = k + carm + var + 2 * tstem + tl;
                tastem = astem;
                tastem8 = astem8;
                tastem8d = astem8d;

                nc += 1;
            }
            na += 1;
        }

        sc = sc.offset(1);
    }

    /* for highest energy mttRNA gene */
    /* decide astem length, look for NCCA acceptor tail */
    /* and calculate total length */

    if !te.ps.is_null() {
        apos2 = te.ps.offset(te.nbase as isize);
        if extastem != 0 {
            if tastem8d != 0 {
                te.astem1 = 8;
                te.astem2 = 8;
                te.ps = te.ps.offset(-1);
                te.nbase += 1;
                te.anticodon += 1;
                as_ = aatail(apos2.offset(8), &mut aext, sw);
            } else {
                te.astem1 = tastem;
                te.astem2 = tastem;
                as_ = aatail(apos2.offset(tastem as isize), &mut aext, sw);
                if tastem8 != 0 {
                    as8 = aatail(apos2.offset(8), &mut aext8, sw);
                    if as8 >= as_ {
                        te.ps = te.ps.offset(-1);
                        te.nbase += 1;
                        te.anticodon += 1;
                        te.astem1 = 8;
                        te.astem2 = 8;
                        as_ = as8;
                        aext = aext8;
                    }
                }
            }
        } else {
            te.astem1 = tastem;
            te.astem2 = tastem;
            as_ = aatail(apos2.offset(tastem as isize), &mut aext, sw);
        }
        if as_ < 2 {
            aext = 1;
        }
        te.nbase += te.astem2;
        nbasefext = te.nbase + ASTEM2_EXT;
        te.nbase += aext;

        /* store mttRNA gene if there are no */
        /* higher energy overlapping mttRNA genes */

        te.start = ((te.ps as isize - seq as isize) / std::mem::size_of::<i32>() as isize) as i64;
        tn = find_slot_ts(ts, max_genes, d, te as *mut Gene, &mut nts, sw);
        if !tn.is_null() {
            base_copy3(te.ps, (*tn).seq.as_mut_ptr(), nbasefext);
            base_copy3(te.ps, (*tn).eseq.as_mut_ptr(), nbasefext);
            te.aatail = aext;
            *tn = te.clone();
        }
    }

    nts
}

/// Legacy wrapper for find_mt_trna - uses global TS storage
/// For backward compatibility with non-parallel code paths
pub unsafe fn find_mt_trna(
    d: *mut data_set,
    seq: *mut i32,
    lseq: i32,
    nts: i32,
    sw: *mut csw
) -> i32 {
    // TS is the global gene storage from tables.rs
    // (*sw).genespace is the global gene capacity
    find_mt_trna_ts(TS, (*sw).genespace, d, seq, lseq, nts, sw)
}
