/*
 * Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
 * tmrna.rs - tmRNA gene detection functions
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

use crate::types::*;
use crate::sequence::vloop_stability;
use crate::trna::{find_slot, find_slot_ts, base_copy3, remove_intron, aatail, find_tstems, find_astem5, find_resume_seq};
use crate::mtrna::{find_mt_trna, find_mt_trna_ts};
use crate::thread_state::{get_tmiopt_state, get_tmopt_state};

// Local aliases for gene types (using types.rs constants)
const tRNA_: i32 = TRNA;
const tmRNA_: i32 = TMRNA;

// Local i32 versions of usize constants for arithmetic
const MAXTRNALEN_I32: i32 = MAXTRNALEN as i32;
const MINTRNALEN_I32: i32 = MINTRNALEN;

// Gene storage and tables from library
use crate::tables::{TS, BP};

/// ti_genedetected_ts - Process detected tRNA intron gene (with explicit ts parameter)
pub unsafe fn ti_genedetected_ts(
    ts: *mut Gene,
    max_genes: i32,
    d: *mut data_set,
    mut nts: i32,
    seq: *mut i32,
    te: *mut Gene,
    sw: *mut csw
) -> i32 {
    let mut as_: i32;
    let mut aext: i32 = 0;
    let mut as8: i32;
    let mut aext8: i32 = 0;
    let mut nbasefext: i32;
    let mut s: *mut i32;
    let mut pseq: [i32; 2 * MAXETRNALEN + 1] = [0; 2 * MAXETRNALEN + 1];
    let mut tn: *mut Gene;

    (*te).nbase = (*te).astem1 + (*te).spacer1 + (*te).spacer2 + 2 * (*te).dstem
                + (*te).dloop + 2 * (*te).cstem + (*te).cloop
                + (*te).var + 2 * (*te).tstem + (*te).tloop + (*te).astem2;
    s = (*te).ps.offset(((*te).nbase + (*te).nintron) as isize);
    as_ = aatail(s, &mut aext, sw);

    if (*sw).extastem != 0 {
        if (*te).astem1 == 7 {
            if *(*BP.get_unchecked(*(*te).ps.offset(-1) as usize)).get_unchecked(*s as usize) != 0 {
                as8 = aatail(s.offset(1), &mut aext8, sw);
                if as8 >= as_ {
                    (*te).ps = (*te).ps.offset(-1);
                    (*te).nbase += 2;
                    (*te).anticodon += 1;
                    if (*te).nintron > 0 {
                        (*te).intron += 1;
                    }
                    (*te).astem1 = 8;
                    (*te).astem2 = 8;
                    as_ = as8;
                    aext = aext8;
                }
            }
        }
    }

    nbasefext = (*te).nbase + ASTEM2_EXT;
    (*te).nbase += aext;
    (*te).start = (*te).ps.offset_from(seq) as i64;
    tn = find_slot_ts(ts, max_genes, d, te, &mut nts, sw);

    if !tn.is_null() {
        if (*te).nintron == 0 {
            base_copy3((*te).ps, (*te).seq.as_mut_ptr() as *mut i32, nbasefext);
        } else {
            base_copy3((*te).ps, (*te).eseq.as_mut_ptr() as *mut i32, nbasefext + (*te).nintron);
            remove_intron((*te).ps, pseq.as_mut_ptr(), nbasefext, (*te).intron, (*te).nintron);
            base_copy3(pseq.as_mut_ptr(), (*te).seq.as_mut_ptr() as *mut i32, nbasefext);
        }
        (*te).aatail = aext;
        std::ptr::copy_nonoverlapping(te, tn, 1);
    }
    nts
}

/// ti_genedetected - Process detected tRNA intron gene (legacy wrapper using global TS)
/// tmrna.c:11-47
pub unsafe fn ti_genedetected(
    d: *mut data_set,
    mut nts: i32,
    seq: *mut i32,
    te: *mut Gene,
    sw: *mut csw
) -> i32 {
    ti_genedetected_ts(TS, (*sw).genespace, d, nts, seq, te, sw)
}

/// tmopt_ts - Optimize tmRNA detection (thread-safe version)
/// Uses ts parameter instead of global TS
pub unsafe fn tmopt_ts(
    ts: *mut Gene,
    max_genes: i32,
    d: *mut data_set,
    th: *mut trna_loop,
    tarm: i32,
    the: f64,
    ahit: *mut trna_loop,
    nah: i32,
    mut nts: i32,
    seq: *mut i32,
    sw: *mut csw
) -> i32 {
    let mut r: i32;
    let mut na: i32;
    let mut nr: i32;
    let mut nrh: i32;
    let mut ibase: i32;
    let mut flag: i32;
    let mut as_: i32;
    let mut aext: i32 = 0;
    let mut nbasefext: i32;
    let mut s: *mut i32;
    let mut v: *mut i32;
    let mut s1: *mut i32;
    let mut s2: *mut i32;
    let mut sa: *mut i32;
    let mut sb: *mut i32;
    let mut se: *mut i32;
    let mut sf: *mut i32;
    let mut ps: *mut i32;
    let mut tpos: *mut i32;
    let mut pseq: [i32; MAXETRNALEN + 1] = [0; MAXETRNALEN + 1];

    static gtem: [i32; 6] = [0x00, 0x00, 0x11, 0x00, 0x00, 0x00];
    static A: [f64; 6] = [6.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    static Ar: [f64; 6] = [10.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    static Cr: [f64; 6] = [0.0, 10.0, 0.0, 0.0, 0.0, 0.0];
    static G: [f64; 6] = [0.0, 0.0, 6.0, 0.0, 0.0, 0.0];
    static Ga: [f64; 6] = [0.0, 0.0, 7.0, 0.0, 0.0, 0.0];
    static K: [f64; 6] = [0.0, 0.0, 6.0, 6.0, 0.0, 0.0];
    static Tr: [f64; 6] = [0.0, 0.0, 0.0, 10.0, 0.0, 0.0];

    let mut e: f64;
    let mut energy: f64;
    let mut penergy: f64;
    let mut tenergy: f64;
    let mut aenergy: f64;
    let athresh: f64;
    let cthresh: f64;
    let cathresh: f64;

    static bem: [[f64; 6]; 6] = [
        [-1.072, -0.214, -1.072, ATBOND, 0.000, 0.000],
        [-0.214, -1.072,  3.000, -1.072, 0.000, 0.000],
        [-1.072,  3.000, -1.072,  1.286, 0.000, 0.000],
        [ATBOND, -1.072,  1.286, -0.214, 0.000, 0.000],
        [ 0.000,  0.000,  0.000,  0.000, 0.000, 0.000],
        [ 0.000,  0.000,  0.000,  0.000, 0.000, 0.000]
    ];

    // Thread-local state for rhit array
    let tl_state = &mut *get_tmopt_state();
    let rhit = &mut tl_state.rhit;

    let mut te: Gene = Gene::default();
    let mut tn: *mut Gene;

    // Static gene template
    let mut t: Gene = Gene::default();
    t.astem1 = 7;
    t.astem2 = 7;
    t.aatail = 1;
    t.spacer1 = 0;
    t.spacer2 = 0;
    t.dstem = 0;
    t.dloop = 13;
    t.cstem = 8;
    t.cloop = 0;
    t.intron = 28;
    t.nintron = 0;
    t.anticodon = 0;
    t.var = 3;
    t.varbp = 0;
    t.tstem = 5;
    t.tloop = 7;
    t.genetype = tmRNA_;
    t.energy = 0.0;
    t.asst = 0;
    t.tps = 0;
    t.tpe = 0;

    tpos = (*th).pos;
    flag = 0;
    te.energy = (*sw).tmrnathresh;
    athresh = (*sw).tmathresh;
    cthresh = (*sw).tmcthresh;
    cathresh = (*sw).tmcathresh;
    s = tpos.offset((tarm + 4) as isize);
    v = tpos.offset(((*th).stem - 10) as isize);
    energy = *K.get_unchecked(*v as usize) + *G.get_unchecked(*v.offset(1) as usize) + *A.get_unchecked(*v.offset(2) as usize);
    e = *K.get_unchecked(*v.offset(1) as usize) + *G.get_unchecked(*v.offset(2) as usize) + *A.get_unchecked(*v.offset(3) as usize);
    if e > energy {
        energy = e;
    }
    if energy < 18.0 {
        energy = 0.0;
    }
    tenergy = *Tr.get_unchecked(*s as usize) + *Cr.get_unchecked(*s.offset(1) as usize) + *Cr.get_unchecked(*s.offset(2) as usize)
            + *Ar.get_unchecked(*s.offset(3) as usize) + energy + 1.59 * the;

    nrh = find_resume_seq(tpos.offset(-(MAXTPTSDIST as isize)), TPWINDOW, rhit.as_mut_ptr(), NH as i32, sw);
    nr = -1;
    nr += 1;
    while nr < nrh {
        ps = (*rhit.get_unchecked(nr as usize)).pos;
        penergy = tenergy + (*rhit.get_unchecked(nr as usize)).energy - 0.001 * ((tpos as isize - ps as isize) as f64 / std::mem::size_of::<i32>() as f64);
        if (*rhit.get_unchecked(nr as usize)).stem < 24 {
            penergy -= 15.0;
        }
        na = -1;
        na += 1;
        while na < nah {
            aenergy = (*ahit.offset(na as isize)).energy;
            if aenergy < athresh {
                na += 1;
                continue;
            }
            t.ps = (*ahit.offset(na as isize)).pos;
            if t.ps < ps.offset(-(MAXTPDIST as isize)) {
                na += 1;
                continue;
            }
            if t.ps > ps.offset(-(MINTPDIST as isize)) {
                break;
            }
            energy = -INACTIVE;
            sa = t.ps.offset(t.astem1 as isize);
            sb = sa.offset(9);
            while sb <= sa.offset(16) {
                sf = tpos.offset(-3);
                while sf >= tpos.offset(-7) {
                    s1 = sb;
                    s2 = sf;
                    se = s1.offset(t.cstem as isize);
                    s2 = s2.offset(-1);
                    e = *(*bem.get_unchecked(*s1 as usize)).get_unchecked(*s2 as usize);
                    s1 = s1.offset(1);
                    while s1 < se {
                        s2 = s2.offset(-1);
                        e += *(*bem.get_unchecked(*s1 as usize)).get_unchecked(*s2 as usize);
                        s1 = s1.offset(1);
                    }
                    if e > energy {
                        energy = e;
                        t.var = (tpos as isize - sf as isize) as i32 / std::mem::size_of::<i32>() as i32;
                        t.dloop = (sb as isize - sa as isize) as i32 / std::mem::size_of::<i32>() as i32;
                    }
                    sf = sf.offset(-1);
                }
                sb = sb.offset(1);
            }
            if energy < cthresh {
                na += 1;
                continue;
            }
            energy += aenergy;
            if energy < cathresh {
                na += 1;
                continue;
            }
            sb = sa.offset(3);
            sf = sa.offset(7);
            r = *gtem.get_unchecked(*sb as usize);
            sb = sb.offset(1);
            while sb < sf {
                r = (r >> 4) + *gtem.get_unchecked(*sb as usize);
                sb = sb.offset(1);
                if (r & 3) == 2 {
                    energy += 14.0;
                    break;
                }
            }
            t.energy = penergy + *Ga.get_unchecked(*t.ps.offset(1) as usize) + *Ga.get_unchecked(*t.ps.offset(2) as usize) + energy;
            if t.energy > te.energy {
                flag = 1;
                t.tstem = (*th).stem;
                t.tloop = (*th).loop_;
                t.tps = (ps as isize - t.ps as isize) as i32 / std::mem::size_of::<i32>() as i32;
                t.tpe = t.tps + (*rhit.get_unchecked(nr as usize)).stem;
                ibase = (tpos as isize - t.ps as isize) as i32 / std::mem::size_of::<i32>() as i32;
                t.nintron = ibase - t.var - 2 * t.cstem - t.dloop - t.astem1;
                t.nbase = ibase + tarm + t.astem2 - t.nintron;
                std::ptr::copy_nonoverlapping(&t as *const Gene, &mut te as *mut Gene, 1);
            }
            na += 1;
        }
        nr += 1;
    }

    if flag != 0 {
        te.start = (te.ps as isize - seq as isize) as i64 / std::mem::size_of::<i32>() as i64;
        s = te.ps.offset((te.nbase + te.nintron) as isize);
        as_ = aatail(s, &mut aext, sw);
        nbasefext = te.nbase + ASTEM2_EXT;
        te.nbase += aext;
        tn = find_slot_ts(ts, max_genes, d, &mut te, &mut nts, sw);
        if !tn.is_null() {
            te.intron = te.astem1 + te.dloop + te.cstem;
            te.asst = 0;
            base_copy3(te.ps, te.eseq.as_mut_ptr() as *mut i32, nbasefext + te.nintron);
            remove_intron(te.ps, pseq.as_mut_ptr(), nbasefext, te.intron, te.nintron);
            base_copy3(pseq.as_mut_ptr(), te.seq.as_mut_ptr() as *mut i32, te.nbase);
            te.aatail = aext;
            *tn = te;
        }
    }
    nts
}

/// tmopt - Optimize tmRNA detection (legacy wrapper using global TS)
/// tmrna.c:50-154
pub unsafe fn tmopt(
    d: *mut data_set,
    th: *mut trna_loop,
    tarm: i32,
    the: f64,
    ahit: *mut trna_loop,
    nah: i32,
    nts: i32,
    seq: *mut i32,
    sw: *mut csw
) -> i32 {
    tmopt_ts(TS, (*sw).genespace, d, th, tarm, the, ahit, nah, nts, seq, sw)
}

/// tmopt_perm_ts - Optimize permuted tmRNA detection (thread-safe version)
/// Uses ts parameter instead of global TS
pub unsafe fn tmopt_perm_ts(
    ts: *mut Gene,
    max_genes: i32,
    d: *mut data_set,
    th: *mut trna_loop,
    tarm: i32,
    the: f64,
    ahit: *mut trna_loop,
    nah: i32,
    mut nts: i32,
    seq: *mut i32,
    sw: *mut csw
) -> i32 {
    let mut r: i32;
    let mut na: i32;
    let mut nr: i32;
    let mut nrh: i32;
    let mut flag: i32;
    let mut as_: i32;
    let mut aext: i32 = 0;
    let mut s: *mut i32;
    let mut v: *mut i32;
    let mut s1: *mut i32;
    let mut s2: *mut i32;
    let mut sa: *mut i32;
    let mut sb: *mut i32;
    let mut se: *mut i32;
    let mut sf: *mut i32;
    let mut ps: *mut i32;
    let mut apos: *mut i32;
    let mut tpos: *mut i32;

    static gtem: [i32; 6] = [0x00, 0x00, 0x11, 0x00, 0x00, 0x00];
    static A: [f64; 6] = [6.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    static Ar: [f64; 6] = [10.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    static Cr: [f64; 6] = [0.0, 10.0, 0.0, 0.0, 0.0, 0.0];
    static G: [f64; 6] = [0.0, 0.0, 6.0, 0.0, 0.0, 0.0];
    static Ga: [f64; 6] = [0.0, 0.0, 7.0, 0.0, 0.0, 0.0];
    static K: [f64; 6] = [0.0, 0.0, 6.0, 6.0, 0.0, 0.0];
    static Tr: [f64; 6] = [0.0, 0.0, 0.0, 10.0, 0.0, 0.0];

    let mut e: f64;
    let mut energy: f64;
    let mut penergy: f64;
    let mut tenergy: f64;
    let mut aenergy: f64;
    let athresh: f64;
    let cthresh: f64;
    let cathresh: f64;

    static bem: [[f64; 6]; 6] = [
        [-1.072, -0.214, -1.072, ATBOND, 0.000, 0.000],
        [-0.214, -1.072, 3.000, -1.072, 0.000, 0.000],
        [-1.072, 3.000, -1.072, 1.286, 0.000, 0.000],
        [ATBOND, -1.072, 1.286, -0.214, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    ];

    // Thread-local state for rhit array
    let tl_state = &mut *get_tmopt_state();
    let rhit = &mut tl_state.rhit;

    let mut te: Gene = Gene::default();
    let mut tn: *mut Gene;

    // Static gene template
    let mut t: Gene = Gene::default();
    t.astem1 = 7;
    t.astem2 = 7;
    t.aatail = 1;
    t.spacer1 = 0;
    t.spacer2 = 0;
    t.dstem = 0;
    t.dloop = 13;
    t.cstem = 8;
    t.cloop = 0;
    t.intron = 28;
    t.nintron = 0;
    t.anticodon = 0;
    t.var = 3;
    t.varbp = 0;
    t.tstem = 5;
    t.tloop = 7;
    t.genetype = tmRNA_;
    t.energy = 0.0;
    t.asst = 0;
    t.tps = 0;
    t.tpe = 0;

    tpos = (*th).pos;
    flag = 0;
    te.energy = (*sw).tmrnathresh;
    athresh = (*sw).tmathresh;
    cthresh = (*sw).tmcthresh;
    cathresh = (*sw).tmcathresh;
    s = tpos.offset((tarm + 4) as isize);
    v = tpos.offset(((*th).stem - 10) as isize);
    energy = *K.get_unchecked(*v as usize) + *G.get_unchecked(*v.offset(1) as usize) + *A.get_unchecked(*v.offset(2) as usize);
    e = *K.get_unchecked(*v.offset(1) as usize) + *G.get_unchecked(*v.offset(2) as usize) + *A.get_unchecked(*v.offset(3) as usize);
    if e > energy {
        energy = e;
    }
    if energy < 18.0 {
        energy = 0.0;
    }
    tenergy = *Tr.get_unchecked(*s as usize) + *Cr.get_unchecked(*s.offset(1) as usize) + *Cr.get_unchecked(*s.offset(2) as usize)
            + *Ar.get_unchecked(*s.offset(3) as usize) + energy + 1.59 * the;

    na = -1;
    na += 1;
    while na < nah {
        aenergy = (*ahit.offset(na as isize)).energy;
        if aenergy < athresh {
            na += 1;
            continue;
        }
        apos = (*ahit.offset(na as isize)).pos;
        if apos < tpos.offset(MINTSTEM_DIST as isize) {
            na += 1;
            continue;
        }
        if apos > tpos.offset((MAXTSTEM_DIST + MAXPPINTRONDIST) as isize) {
            break;
        }
        energy = -INACTIVE;
        sa = apos.offset(t.astem1 as isize);
        sb = sa.offset(9);
        se = sb.offset(t.cstem as isize);
        while sb <= sa.offset(16) {
            sf = tpos.offset(-3);
            while sf >= tpos.offset(-7) {
                s1 = sb;
                s2 = sf;
                s2 = s2.offset(-1);
                e = *(*bem.get_unchecked(*s1 as usize)).get_unchecked(*s2 as usize);
                s1 = s1.offset(1);
                while s1 < se {
                    s2 = s2.offset(-1);
                    e += *(*bem.get_unchecked(*s1 as usize)).get_unchecked(*s2 as usize);
                    s1 = s1.offset(1);
                }
                if e > energy {
                    energy = e;
                    t.var = ((tpos as isize - sf as isize) / std::mem::size_of::<i32>() as isize) as i32;
                    t.dloop = ((sb as isize - sa as isize) / std::mem::size_of::<i32>() as isize) as i32;
                }
                sf = sf.offset(-1);
            }
            sb = sb.offset(1);
            se = se.offset(1);
        }
        if energy < cthresh {
            na += 1;
            continue;
        }
        energy += aenergy;
        if energy < cathresh {
            na += 1;
            continue;
        }
        sb = sa.offset(3);
        sf = sa.offset(7);
        r = *gtem.get_unchecked(*sb as usize);
        sb = sb.offset(1);
        while sb < sf {
            r = (r >> 4) + *gtem.get_unchecked(*sb as usize);
            sb = sb.offset(1);
            if (r & 3) == 2 {
                energy += 14.0;
                break;
            }
        }
        penergy = tenergy + *Ga.get_unchecked(*apos.offset(1) as usize) + *Ga.get_unchecked(*apos.offset(2) as usize) + energy;
        nrh = find_resume_seq(apos.offset(MINTPDIST as isize), TPWINDOW, rhit.as_mut_ptr(), NH as i32, sw);
        nr = -1;
        nr += 1;
        while nr < nrh {
            ps = (*rhit.get_unchecked(nr as usize)).pos;
            t.energy = penergy + (*rhit.get_unchecked(nr as usize)).energy;
            if (*rhit.get_unchecked(nr as usize)).stem < 24 {
                t.energy -= 15.0;
            }
            if t.energy > te.energy {
                flag = 1;
                t.tstem = (*th).stem;
                t.tloop = (*th).loop_;
                t.asst = ((apos as isize - tpos as isize) / std::mem::size_of::<i32>() as isize) as i32 + t.var + t.cstem;
                t.ps = tpos.offset(-t.var as isize).offset(-(t.cstem as isize));
                t.tps = ((ps as isize - t.ps as isize) / std::mem::size_of::<i32>() as isize) as i32;
                t.tpe = t.tps + (*rhit.get_unchecked(nr as usize)).stem;
                std::ptr::copy_nonoverlapping(&t, &mut te, 1);
            }
            nr += 1;
        }
        na += 1;
    }

    if flag != 0 {
        te.start = ((te.ps as isize - seq as isize) / std::mem::size_of::<i32>() as isize) as i64 - 54;
        te.intron = te.cstem + te.var + 2 * te.tstem + te.tloop + te.astem2;
        as_ = aatail(te.ps.offset(te.intron as isize), &mut aext, sw);
        te.aatail = aext;
        base_copy3(te.ps.offset(-54), te.eseq.as_mut_ptr() as *mut i32, te.tpe + 1 + TMPTRAILER);
        te.nbase = te.astem1 + te.dloop + te.cstem;
        base_copy3(te.ps.offset(te.asst as isize), te.seq.as_mut_ptr() as *mut i32, te.nbase);
        base_copy3(te.ps, te.seq.as_mut_ptr().offset(te.nbase as isize) as *mut i32, te.intron + ASTEM2_EXT);
        te.intron += aext;
        te.nbase += te.intron;
        te.nintron = te.tpe - te.nbase + 1 + TMPTRAILER;
        te.intron += 54;
        te.tps += 54;
        te.tpe += 54;
        te.asst += 54;
        tn = find_slot_ts(ts, max_genes, d, &mut te, &mut nts, sw);
        if !tn.is_null() {
            *tn = te;
        }
    }
    nts
}

/// tmopt_perm - Optimize permuted tmRNA detection (legacy wrapper using global TS)
/// tmrna.c:157-262
pub unsafe fn tmopt_perm(
    d: *mut data_set,
    th: *mut trna_loop,
    tarm: i32,
    the: f64,
    ahit: *mut trna_loop,
    nah: i32,
    nts: i32,
    seq: *mut i32,
    sw: *mut csw
) -> i32 {
    tmopt_perm_ts(TS, (*sw).genespace, d, th, tarm, the, ahit, nah, nts, seq, sw)
}

/// tmioptimise_ts - Main tmRNA/tRNA intron detection optimizer (thread-safe version)
/// Uses ts parameter instead of global TS
pub unsafe fn tmioptimise_ts(
    ts: *mut Gene,
    max_genes: i32,
    d: *mut data_set,
    seq: *mut i32,
    lseq: i32,
    mut nts: i32,
    sw: *mut csw
) -> i32 {
    // This is a large function with many local variables
    // Implementing the full logic here

    let mut i: i32;
    let mut j: i32;
    let mut k: i32;
    let mut intron: i32;
    let mut nt: i32;
    let mut nth: i32;
    let mut nd1: i32;
    let mut nd2: i32;
    let mut ndx: i32;
    let mut ndh: i32;
    let mut na: i32;
    let mut nah: i32;
    let mut nppah: i32;
    let mut nc: i32;
    let mut nch: i32;
    let mut tfold: i32;
    let mut tarm: i32;
    let mut dfl_overflow: i32;
    let mut dstem: i32;
    let mut dloop: i32;
    let mut flag: i32;
    let mindist: i32;
    let maxdist: i32;
    let tmindist: i32;
    let tmaxdist: i32;
    let tmmindist: i32;
    let tmmaxdist: i32;
    let tarmthresh: i32;
    let tmstrict: i32;
    let sp2min: i32;
    let sp2max: i32;

    let ethresh: f64;
    let mut he: f64;
    let mut the: f64;
    let mut thet: f64;
    let tdarmthresh: f64;
    let mut energy: f64;
    let mut e: f64;

    // Get thread-local state for mutable arrays
    let tl_state = unsafe { &mut *get_tmiopt_state() };

    static G: [f64; 6] = [0.0, 0.0, 6.0, 0.0, 0.0, 0.0];

    if (*sw).mtrna != 0 {
        nts = find_mt_trna_ts(ts, max_genes, d, seq, lseq, nts, sw);
        if (*sw).tmrna == 0 {
            return nts;
        }
    }

    ethresh = (*sw).trnathresh;
    tmmindist = MINTPTSDIST + MINTPDIST;
    tmmaxdist = MAXTPTSDIST + MAXTPDIST;
    tmindist = MINTRNALEN_I32 + (*sw).minintronlen - MAXTSTEM_DIST;
    tmaxdist = MAXTRNALEN_I32 + (*sw).maxintronlen - MINTSTEM_DIST;

    if (*sw).trna != 0 {
        if (*sw).tmrna != 0 {
            mindist = if tmindist < tmmindist { tmindist } else { tmmindist };
            maxdist = if tmaxdist > tmmaxdist { tmaxdist } else { tmmaxdist };
        } else {
            mindist = tmindist;
            maxdist = tmaxdist;
        }
    } else {
        mindist = tmmindist;
        maxdist = tmmaxdist;
    }

    tarmthresh = (*sw).ttarmthresh as i32;
    tdarmthresh = (*sw).tdarmthresh;
    tmstrict = (*sw).tmstrict;
    sp2min = (*sw).sp2min;
    sp2max = (*sw).sp2max;

    nth = find_tstems(seq, lseq, tl_state.thit.as_mut_ptr(), NTH as i32, sw);

    nt = -1;
    nt += 1;
    while nt < nth {
        let tpos = tl_state.thit[nt as usize].pos;
        let t_tloop = tl_state.thit[nt as usize].loop_;
        let t_tstem = tl_state.thit[nt as usize].stem;
        tfold = *tpos.offset(-1);
        tarm = 2 * t_tstem + t_tloop;
        let tend = tpos.offset(tarm as isize);
        let tmv = tpos.offset(-(VARMIN as isize));
        flag = 0;
        // Optimization: Reuse thread-local te instead of creating 13KB Gene each iteration
        // Only te.energy is read before copy from t; other fields are overwritten
        tl_state.te.energy = ethresh;
        the = tl_state.thit[nt as usize].energy;

        nah = find_astem5(
            tpos.offset(-(maxdist as isize)),
            tpos.offset(-(mindist as isize)),
            tend,
            7,
            tl_state.ahit_arr.as_mut_ptr(),
            NA as i32,
            sw
        );

        if (*sw).tmrna != 0 {
            thet = the - *G.get_unchecked(*tpos.offset(t_tstem as isize) as usize)
                      - *G.get_unchecked(*tpos.offset((t_tstem + 1) as isize) as usize);
            if tmstrict != 0 {
                if thet >= tarmthresh as f64 {
                    nts = tmopt_ts(ts, max_genes, d, tl_state.thit.as_mut_ptr().offset(nt as isize), tarm, thet,
                               tl_state.ahit_arr.as_mut_ptr(), nah, nts, seq, sw);
                }
            } else {
                nts = tmopt_ts(ts, max_genes, d, tl_state.thit.as_mut_ptr().offset(nt as isize), tarm, the,
                           tl_state.ahit_arr.as_mut_ptr(), nah, nts, seq, sw);
            }
            nppah = find_astem5(
                tpos.offset(MINPPASDIST as isize),
                tpos.offset(MAXPPASDIST as isize),
                tend,
                7,
                tl_state.ahit_arr.as_mut_ptr().offset(nah as isize),
                (NA as i32) - nah,
                sw
            );
            nts = tmopt_perm_ts(ts, max_genes, d, tl_state.thit.as_mut_ptr().offset(nt as isize), tarm, the,
                            tl_state.ahit_arr.as_mut_ptr().offset(nah as isize), nppah, nts, seq, sw);
            if thet < tarmthresh as f64 {
                nt += 1;
                continue;
            }
            the = thet;
        } else if (*sw).threshlevel < 1.0 {
            the -= *G.get_unchecked(*tpos.offset(t_tstem as isize) as usize)
                 + *G.get_unchecked(*tpos.offset((t_tstem + 1) as isize) as usize);
            if the < tarmthresh as f64 {
                nt += 1;
                continue;
            }
        }

        if (*sw).trna == 0 {
            nt += 1;
            continue;
        }

        // tmrna.c:391-723: D-stem, C-stem finding and gene assembly
        // Static arrays for D-stem and C-stem detection
        static TT: [u32; 6] = [0x00, 0x00, 0x00, 0x11, 0x00, 0x00];
        static GG: [u32; 6] = [0x00, 0x00, 0x11, 0x00, 0x00, 0x00];
        static cA: [u32; 6] = [0, 0, 0, 2, 0, 0];
        static cC: [u32; 6] = [0, 0, 2, 0, 0, 0];
        static cG: [u32; 6] = [0, 2, 0, 1, 0, 0];
        static cT: [u32; 6] = [2, 0, 1, 0, 0, 0];
        static yic: [i32; 9] = [1, 0, 0, 0, 0, 0, 0, 0, 0];
        static tic: [i32; 9] = [1, 1, 0, 0, 0, 0, 0, 0, 0];
        static a1ic: [i32; 9] = [1, 1, 1, 0, 0, 0, 0, 0, 0];
        static a2ic: [i32; 9] = [1, 1, 1, 1, 0, 0, 0, 0, 0];
        static a3ic: [i32; 9] = [1, 1, 1, 1, 1, 0, 0, 0, 0];
        static ric: [i32; 9] = [1, 1, 1, 1, 1, 1, 0, 0, 0];
        static goffb: [i32; 13] = [0, 0, 0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 2];
        static goffe: [i32; 13] = [0, 0, 0, 0, 2, 3, 4, 4, 5, 6, 6, 6, 6];
        static cY: [i32; 6] = [0, 1, 0, 1, 0, 0];
        static cR: [i32; 6] = [1, 0, 1, 0, 0, 0];
        static ilw: f64 = 0.002;
        static T: [f64; 6] = [0.0, 0.0, 0.0, 7.0, 0.0, 0.0];
        static Y: [f64; 6] = [0.0, 3.0, 0.0, 3.0, 0.0, 0.0];
        static R: [f64; 6] = [2.0, 0.0, 2.0, 0.0, 0.0, 0.0];
        static YP: [f64; 6] = [0.0, 3.0, 0.0, 3.0, 0.0, 0.0];
        static RP: [f64; 6] = [2.0, 0.0, 2.0, 0.0, 0.0, 0.0];
        static RI: [f64; 6] = [0.1, 0.0, 0.05, 0.0, 0.0, 0.0];
        static GI: [f64; 6] = [0.0, 0.0, 0.1, 0.0, 0.0, 0.0];
        static YI: [f64; 6] = [0.0, 0.1, 0.0, 0.1, 0.0, 0.0];
        static AI: [f64; 6] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        static GC: [f64; 6] = [0.0, 1.5, 6.0, 0.0, 0.0, 0.0];
        static G3: [f64; 6] = [0.0, 6.0, 12.0, 12.0, 0.0, 0.0];
        static dR: [f64; 6] = [6.0, 0.0, 6.0, 0.0, 0.0, 0.0];
        static RH: [f64; 6] = [3.0, 0.0, 3.0, 0.0, 0.0, 0.0];
        static AGT: [f64; 6] = [6.0, 0.0, 6.0, 6.0, 0.0, 0.0];
        static dT: [f64; 6] = [0.0, 0.0, 0.0, 6.0, 0.0, 0.0];
        // Flattened 1D arrays for better cache performance (6x6 = 36 elements)
        // Access: arr_flat[row * 6 + col] instead of arr[row][col]
        static dbem_flat: [f64; 36] = [
            -2.144, -0.428, -2.144, ATBOND, 0.000, 0.000,
            -0.428, -2.144, 3.000, -2.144, 0.000, 0.000,
            -2.144, 3.000, -2.144, 1.286, 0.000, 0.000,
            ATBOND, -2.144, 1.286, -0.428, 0.000, 0.000,
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
        ];
        static dfem_flat: [f64; 36] = [
            -4.000, -4.000, -4.000, ATBOND, 0.000, 0.000,
            -4.000, -4.000, 3.000, -4.000, 0.000, 0.000,
            -4.000, 3.000, -4.000, 1.286, 0.000, 0.000,
            ATBOND, -4.000, 1.286, -4.000, 0.000, 0.000,
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
        ];
        static cbem_flat: [f64; 36] = [
            -1.072, -0.214, -1.072, 2.0 * ATBOND, 0.000, 0.000,
            -0.214, -1.072, 6.000, -1.072, 0.000, 0.000,
            -1.072, 6.000, -1.072, 3.400, 0.000, 0.000,
            2.0 * ATBOND, -1.072, 3.400, -0.214, 0.000, 0.000,
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
        ];

        // Gene template for tRNA - use thread-local template
        // Optimization: Template t is pre-initialized in TmioptState::default()
        // Only tstem and tloop change per iteration (C version: static gene t)
        tl_state.t.tstem = t_tstem;
        tl_state.t.tloop = t_tloop;

        let tloopfold = tpos.offset((t_tstem + 1) as isize);

        // tmrna.c:391-392
        na = -1;
        na += 1;
        while na < nah {
            let apos = tl_state.ahit_arr[na as usize].pos;
            if apos < tpos.offset(-(tmaxdist as isize)) {
                na += 1;
                continue;
            }
            if apos > tpos.offset(-(tmindist as isize)) {
                break;
            }
            he = the + tl_state.ahit_arr[na as usize].energy;

            // tmrna.c:400-508: find dstems
            ndh = 0;
            dfl_overflow = 0;
            let mut sc = apos.offset(8);
            let mut energyf = *dfem_flat.get_unchecked((*sc.offset(5) as usize) * 6 + (tfold as usize));
            let sl = sc.offset((*sw).sp1max as isize);

            while sc < sl {
                let energy2 = *dT.get_unchecked(*sc.offset(-2) as usize) + *RH.get_unchecked(*sc.offset(-1) as usize)
                    + *GC.get_unchecked(*sc as usize) + *dfem_flat.get_unchecked((*sc.offset(-2) as usize) * 6 + (*sc.offset(4) as usize));
                let energyf6 = *dfem_flat.get_unchecked((*sc.offset(6) as usize) * 6 + (tfold as usize));

                dstem = 3;
                while dstem <= 4 {
                    let mut sd = sc.offset(dstem as isize);
                    dloop = 3;
                    let mut se = sd.offset(dloop as isize);
                    energy = energy2 + 6.0 + *dR.get_unchecked(*se.offset(-1) as usize) + energyf;
                    if dstem == 3 {
                        if energyf < 0.0 {
                            energyf = energyf6;
                        }
                    }
                    se = se.offset(dstem as isize);
                    let mut s1 = sc;
                    let mut s2 = se;
                    let sf = s1.offset(dstem as isize);
                    while s1 < sf {
                        s2 = s2.offset(-1);
                        energy += *dbem_flat.get_unchecked((*s1 as usize) * 6 + (*s2 as usize));
                        s1 = s1.offset(1);
                    }
                    if energy >= tdarmthresh {
                        if ndh >= ND as i32 {
                            dfl_overflow = 1;
                            break;
                        }
                        tl_state.dhit[ndh as usize].pos = sc;
                        tl_state.dhit[ndh as usize].end = se;
                        tl_state.dhit[ndh as usize].loop_ = dloop;
                        tl_state.dhit[ndh as usize].stem = dstem;
                        tl_state.dhit[ndh as usize].energy = energy;
                        ndh += 1;
                    }
                    if dfl_overflow != 0 {
                        break;
                    }

                    // tmrna.c:429-457: dloop 4-11
                    let mut sg1 = sd.offset(1);
                    let sg2 = sd.offset(6);
                    let mut q = *GG.get_unchecked(*sg1 as usize);
                    sg1 = sg1.offset(1);
                    let mut ige: [i32; 7] = [0; 7];
                    *ige.get_unchecked_mut(1) = (q & 3) as i32;
                    j = 2;
                    while sg1 <= sg2 {
                        q = (q >> 4) + *GG.get_unchecked(*sg1 as usize);
                        sg1 = sg1.offset(1);
                        *ige.get_unchecked_mut(j as usize) = (q & 3) as i32;
                        j += 1;
                    }
                    dloop = 4;
                    while dloop <= 11 {
                        j = *goffb.get_unchecked(dloop as usize);
                        k = *goffe.get_unchecked(dloop as usize);
                        let mut c: u32 = *ige.get_unchecked(j as usize) as u32;
                        j += 1;
                        while j <= k {
                            c = c | (*ige.get_unchecked(j as usize) as u32);
                            j += 1;
                        }
                        let genergy = *G3.get_unchecked(c as usize);
                        se = sd.offset(dloop as isize);
                        energy = energy2 + genergy + *dR.get_unchecked(*se.offset(-1) as usize) + energyf;
                        se = se.offset(dstem as isize);
                        s1 = sc;
                        s2 = se;
                        let sf = s1.offset(dstem as isize);
                        while s1 < sf {
                            s2 = s2.offset(-1);
                            energy += *dbem_flat.get_unchecked((*s1 as usize) * 6 + (*s2 as usize));
                            s1 = s1.offset(1);
                        }
                        if energy >= tdarmthresh {
                            if ndh >= ND as i32 {
                                dfl_overflow = 1;
                                break;
                            }
                            tl_state.dhit[ndh as usize].pos = sc;
                            tl_state.dhit[ndh as usize].end = se;
                            tl_state.dhit[ndh as usize].loop_ = dloop;
                            tl_state.dhit[ndh as usize].stem = dstem;
                            tl_state.dhit[ndh as usize].energy = energy;
                            ndh += 1;
                        }
                        dloop += 1;
                    }
                    if dfl_overflow != 0 {
                        break;
                    }
                    dstem += 1;
                }
                if dfl_overflow != 0 {
                    break;
                }

                // tmrna.c:460-482: 6-stem check
                let mut s1 = sc;
                let mut s2 = sc.offset(16);
                let mut sd = sc.offset(6);
                j = *(*BP.get_unchecked(*s1 as usize)).get_unchecked(*s2.offset(-1) as usize);
                s2 = s2.offset(-1);
                s1 = s1.offset(1);
                while s1 < sd {
                    s2 = s2.offset(-1);
                    j += *(*BP.get_unchecked(*s1 as usize)).get_unchecked(*s2 as usize);
                    s1 = s1.offset(1);
                }
                if j >= 6 {
                    energy = *dT.get_unchecked(*sc.offset(-1) as usize) + *RH.get_unchecked(*sc as usize)
                        + *GC.get_unchecked(*sc.offset(1) as usize) + energyf6;
                    sd = sd.offset(1);
                    energy += *G.get_unchecked(*sd as usize);
                    sd = sd.offset(1);
                    energy += *G.get_unchecked(*sd as usize);
                    sd = sd.offset(1);
                    energy += *AGT.get_unchecked(*sd as usize) + *dfem_flat.get_unchecked((*sc.offset(-1) as usize) * 6 + (*sc.offset(4) as usize));
                    sd = sd.offset(7);
                    s1 = sc;
                    s2 = sd;
                    let sf = s1.offset(6);
                    while s1 < sf {
                        s2 = s2.offset(-1);
                        energy += *dbem_flat.get_unchecked((*s1 as usize) * 6 + (*s2 as usize));
                        s1 = s1.offset(1);
                    }
                    if energy >= tdarmthresh {
                        if ndh >= ND as i32 {
                            dfl_overflow = 1;
                        } else {
                            tl_state.dhit[ndh as usize].pos = sc;
                            tl_state.dhit[ndh as usize].end = sd;
                            tl_state.dhit[ndh as usize].loop_ = 4;
                            tl_state.dhit[ndh as usize].stem = 6;
                            tl_state.dhit[ndh as usize].energy = energy;
                            ndh += 1;
                        }
                    }
                }
                if dfl_overflow != 0 {
                    break;
                }

                // tmrna.c:484-506: 7-stem check
                s1 = sc;
                s2 = sc.offset(18);
                sd = sc.offset(7);
                j = *(*BP.get_unchecked(*s1 as usize)).get_unchecked(*s2.offset(-1) as usize);
                s2 = s2.offset(-1);
                s1 = s1.offset(1);
                while s1 < sd {
                    s2 = s2.offset(-1);
                    j += *(*BP.get_unchecked(*s1 as usize)).get_unchecked(*s2 as usize);
                    s1 = s1.offset(1);
                }
                if j >= 7 {
                    energy = energy2 + *dfem_flat.get_unchecked((*sc.offset(7) as usize) * 6 + (tfold as usize));
                    sd = sd.offset(1);
                    energy += *G.get_unchecked(*sd as usize);
                    sd = sd.offset(1);
                    energy += *G.get_unchecked(*sd as usize);
                    sd = sd.offset(1);
                    energy += *AGT.get_unchecked(*sd as usize);
                    sd = sd.offset(8);
                    s1 = sc;
                    s2 = sd;
                    let sf = s1.offset(7);
                    while s1 < sf {
                        s2 = s2.offset(-1);
                        energy += *dbem_flat.get_unchecked((*s1 as usize) * 6 + (*s2 as usize));
                        s1 = s1.offset(1);
                    }
                    if energy >= tdarmthresh {
                        if ndh >= ND as i32 {
                            dfl_overflow = 1;
                        } else {
                            tl_state.dhit[ndh as usize].pos = sc;
                            tl_state.dhit[ndh as usize].end = sd;
                            tl_state.dhit[ndh as usize].loop_ = 4;
                            tl_state.dhit[ndh as usize].stem = 7;
                            tl_state.dhit[ndh as usize].energy = energy;
                            ndh += 1;
                        }
                    }
                }
                energyf = energyf6;
                sc = sc.offset(1);
            }
            if dfl_overflow != 0 {
                eprintln!("Too many D-stem hits");
            }
            // End of find dstems routine

            // tmrna.c:513-534: D-stem energy adjustments
            nd1 = ndh;
            nd1 -= 1;
            while nd1 >= 0 {
                dstem = tl_state.dhit[nd1 as usize].stem;
                let dpos = tl_state.dhit[nd1 as usize].pos;
                if ((dpos as isize - apos as isize) / (std::mem::size_of::<i32>() as isize)) < 9 {
                    tl_state.dhit[nd1 as usize].energy -= 3.0;
                }
                if *tloopfold == Guanine {
                    let mut sb = dpos.offset((dstem + 2) as isize);
                    let sc_inner = sb;
                    let se = sb.offset((tl_state.dhit[nd1 as usize].loop_ - 3) as isize);
                    let mut r: u32 = *TT.get_unchecked(*sb as usize);
                    sb = sb.offset(1);
                    while sb < se {
                        r = (r >> 4) + *TT.get_unchecked(*sb as usize);
                        sb = sb.offset(1);
                        if (r & 2) != 0 {
                            tl_state.dhit[nd1 as usize].energy += 10.0;
                            break;
                        }
                    }
                    let mut sc_gg = sc_inner;
                    let mut r: u32 = *GG.get_unchecked(*sc_gg as usize);
                    sc_gg = sc_gg.offset(1);
                    while sc_gg < se {
                        r = (r >> 4) + *GG.get_unchecked(*sc_gg as usize);
                        sc_gg = sc_gg.offset(1);
                        if (r & 2) != 0 {
                            tl_state.dhit[nd1 as usize].energy -= 12.0;
                            break;
                        }
                    }
                }
                nd1 -= 1;
            }

            // tmrna.c:535-548: Remove duplicate D-stems with same end position
            nd1 = ndh;
            nd1 -= 1;
            while nd1 >= 0 {
                if tl_state.dhit[nd1 as usize].end.is_null() {
                    nd1 -= 1;
                    continue;
                }
                let mut cpos = tl_state.dhit[nd1 as usize].end;
                let mut denergy = tl_state.dhit[nd1 as usize].energy;
                ndx = nd1;
                nd2 = nd1;
                nd2 -= 1;
                while nd2 >= 0 {
                    if tl_state.dhit[nd2 as usize].end != cpos {
                        nd2 -= 1;
                        continue;
                    }
                    e = tl_state.dhit[nd2 as usize].energy;
                    if e > denergy {
                        denergy = e;
                        tl_state.dhit[ndx as usize].end = std::ptr::null_mut();
                        ndx = nd2;
                    }
                    nd2 -= 1;
                }
                nd1 -= 1;
            }

            // tmrna.c:549-561: Find cposmin/cposmax
            let mut cposmin: *mut i32 = std::ptr::null_mut();
            let mut cposmax: *mut i32 = std::ptr::null_mut();
            nd1 = ndh;
            nd1 -= 1;
            while nd1 >= 0 {
                if tl_state.dhit[nd1 as usize].end.is_null() {
                    nd1 -= 1;
                    continue;
                }
                cposmin = tl_state.dhit[nd1 as usize].end;
                cposmax = cposmin;
                break;
            }
            nd2 = nd1;
            nd2 -= 1;
            while nd2 >= 0 {
                let cpos = tl_state.dhit[nd2 as usize].end;
                if cpos.is_null() {
                    nd2 -= 1;
                    continue;
                }
                if cpos < cposmin {
                    cposmin = cpos;
                }
                if cpos > cposmax {
                    cposmax = cpos;
                }
                nd2 -= 1;
            }

            // tmrna.c:562-623: C-stem finding and gene assembly
            if !cposmin.is_null() {
                let mut cpos = cposmin.offset(sp2min as isize);
                while cpos <= cposmax.offset(sp2max as isize) {
                    let mut denergy = -INACTIVE;
                    ndx = -1;
                    nd1 = ndh;
                    nd1 -= 1;
                    while nd1 >= 0 {
                        if tl_state.dhit[nd1 as usize].end.is_null() {
                            nd1 -= 1;
                            continue;
                        }
                        if tl_state.dhit[nd1 as usize].end.offset(sp2max as isize) < cpos {
                            nd1 -= 1;
                            continue;
                        }
                        if tl_state.dhit[nd1 as usize].end.offset(sp2min as isize) > cpos {
                            nd1 -= 1;
                            continue;
                        }
                        e = tl_state.dhit[nd1 as usize].energy;
                        if e > denergy {
                            denergy = e;
                            ndx = nd1;
                        }
                        nd1 -= 1;
                    }
                    if ndx < 0 {
                        cpos = cpos.offset(1);
                        continue;
                    }
                    denergy += he;
                    if denergy < (tl_state.te.energy - 49.0) {
                        cpos = cpos.offset(1);
                        continue;
                    }

                    // tmrna.c:578-622: find cstems
                    nch = 0;
                    let mut si = cpos;
                    let mut sc = cpos.offset(5);
                    let mut se = cpos.offset(4);
                    *tl_state.ct.get_unchecked_mut(0) = *cA.get_unchecked(*se as usize);
                    *tl_state.ct.get_unchecked_mut(1) = *cC.get_unchecked(*se as usize);
                    *tl_state.ct.get_unchecked_mut(2) = *cG.get_unchecked(*se as usize);
                    *tl_state.ct.get_unchecked_mut(3) = *cT.get_unchecked(*se as usize);
                    se = se.offset(-1);
                    while se >= cpos {
                        *tl_state.ct.get_unchecked_mut(0) = (*tl_state.ct.get_unchecked(0) << 4) + *cA.get_unchecked(*se as usize);
                        *tl_state.ct.get_unchecked_mut(1) = (*tl_state.ct.get_unchecked(1) << 4) + *cC.get_unchecked(*se as usize);
                        *tl_state.ct.get_unchecked_mut(2) = (*tl_state.ct.get_unchecked(2) << 4) + *cG.get_unchecked(*se as usize);
                        *tl_state.ct.get_unchecked_mut(3) = (*tl_state.ct.get_unchecked(3) << 4) + *cT.get_unchecked(*se as usize);
                        se = se.offset(-1);
                    }
                    si = si.offset(11);
                    se = tmv.offset(-(VARDIFF as isize)).offset(-5);
                    if si < se {
                        si = se;
                    }
                    let mut r: u32 = *tl_state.ct.get_unchecked(*si as usize);
                    si = si.offset(1);
                    r = (r >> 4) + *tl_state.ct.get_unchecked(*si as usize);
                    si = si.offset(1);
                    r = (r >> 4) + *tl_state.ct.get_unchecked(*si as usize);
                    si = si.offset(1);
                    r = (r >> 4) + *tl_state.ct.get_unchecked(*si as usize);
                    si = si.offset(1);
                    while si < tmv {
                        r = (r >> 4) + *tl_state.ct.get_unchecked(*si as usize);
                        si = si.offset(1);
                        if (r & 0xf) >= 5 {
                            if nch >= NC as i32 {
                                eprintln!("Too many cstem hits");
                                break;
                            }
                            tl_state.chit[nch as usize].pos = si;
                            tl_state.chit[nch as usize].stem = 5;
                            tl_state.chit[nch as usize].loop_ = ((si as isize - sc as isize) / std::mem::size_of::<i32>() as isize - 5) as i32;
                            if tl_state.chit[nch as usize].loop_ == 9 {
                                if *(*BP.get_unchecked(*sc as usize)).get_unchecked(*si.offset(-6) as usize) != 0 {
                                    if *cY.get_unchecked(*sc.offset(2) as usize) != 0 {
                                        if *cR.get_unchecked(*sc.offset(6) as usize) != 0 {
                                            if *cY.get_unchecked(*sc.offset(1) as usize) != 0 {
                                                tl_state.chit[nch as usize].stem = 6;
                                                tl_state.chit[nch as usize].loop_ = 7;
                                            }
                                        }
                                    }
                                }
                            }
                            let mut s1 = cpos;
                            let mut s2 = si;
                            se = s1.offset(tl_state.chit[nch as usize].stem as isize);
                            s2 = s2.offset(-1);
                            tl_state.chit[nch as usize].energy = *cbem_flat.get_unchecked((*s1 as usize) * 6 + (*s2 as usize));
                            s1 = s1.offset(1);
                            while s1 < se {
                                s2 = s2.offset(-1);
                                tl_state.chit[nch as usize].energy += *cbem_flat.get_unchecked((*s1 as usize) * 6 + (*s2 as usize));
                                s1 = s1.offset(1);
                            }
                            nch += 1;
                        }
                    }
                    // end of find cstems routine

                    // tmrna.c:626-723: gene assembly
                    nc = -1;
                    nc += 1;
                    while nc < nch {
                        energy = denergy + tl_state.chit[nc as usize].energy;
                        if energy < (tl_state.te.energy - 19.0) {
                            nc += 1;
                            continue;
                        }
                        let cend = tl_state.chit[nc as usize].pos;
                        tl_state.t.var = ((tpos as isize - cend as isize) / std::mem::size_of::<i32>() as isize) as i32;
                        tl_state.t.cloop = tl_state.chit[nc as usize].loop_;
                        tl_state.t.cstem = tl_state.chit[nc as usize].stem;
                        intron = 0;
                        if tl_state.t.cloop < 9 {
                            if (*sw).minintronlen > 0 {
                                nc += 1;
                                continue;
                            }
                            if (*sw).cloop7 != 0 {
                                if tl_state.t.cloop != 7 {
                                    nc += 1;
                                    continue;
                                }
                            }
                            tl_state.t.nintron = 0;
                            if tl_state.t.var > 17 {
                                energy += vloop_stability(cend, tl_state.t.var, &mut tl_state.t.varbp);
                            }
                            let sb = cpos.offset(tl_state.t.cstem as isize);
                            energy += *T.get_unchecked(*sb.offset(1) as usize) + *Y.get_unchecked(*sb as usize)
                                + *R.get_unchecked(*sb.offset(5) as usize)
                                - 0.05 * (tl_state.t.var as f64)
                                - (if tl_state.t.cloop == 7 { 0.0 } else { 6.0 });
                        } else {
                            tl_state.t.nintron = tl_state.t.cloop - 7;
                            if tl_state.t.nintron > (*sw).maxintronlen {
                                nc += 1;
                                continue;
                            }
                            if tl_state.t.nintron < (*sw).minintronlen {
                                nc += 1;
                                continue;
                            }
                            if tl_state.t.var > 17 {
                                energy += vloop_stability(cend, tl_state.t.var, &mut tl_state.t.varbp);
                            }
                            if energy < (tl_state.te.energy - 9.0) {
                                nc += 1;
                                continue;
                            }
                            tl_state.t.cloop = 7;
                            let sb = cpos.offset(tl_state.t.cstem as isize);
                            let se = sb.offset(tl_state.t.nintron as isize);
                            let mut cenergy: f64;
                            if (*sw).ifixedpos != 0 {
                                intron = 6;
                                cenergy = *YP.get_unchecked(*sb as usize) + *T.get_unchecked(*sb.offset(1) as usize) + *RP.get_unchecked(*sb.offset(5) as usize);
                            } else {
                                cenergy = *YP.get_unchecked(*se as usize) + *T.get_unchecked(*se.offset(1) as usize) + *RP.get_unchecked(*se.offset(5) as usize);
                                let mut ienergy = cenergy + *RI.get_unchecked(*sb as usize) + *GI.get_unchecked(*se.offset(-1) as usize)
                                    + *AI.get_unchecked(*se.offset(-2) as usize) * *YI.get_unchecked(*se.offset(-1) as usize);
                                j = 1;
                                while j <= 7 {
                                    let si = se.offset((j - 1) as isize);
                                    let ec = *YP.get_unchecked(*sb.offset((*yic.get_unchecked(j as usize) * tl_state.t.nintron) as isize) as usize)
                                        + *T.get_unchecked(*sb.offset((*tic.get_unchecked(j as usize) * tl_state.t.nintron + 1) as isize) as usize)
                                        + *RP.get_unchecked(*sb.offset((*ric.get_unchecked(j as usize) * tl_state.t.nintron + 5) as isize) as usize);
                                    let mut e_tmp = ec + *RI.get_unchecked(*sb.offset(j as isize) as usize)
                                        + *GI.get_unchecked(*si as usize)
                                        + *AI.get_unchecked(*si.offset(-1) as usize) * *YI.get_unchecked(*si as usize);
                                    if j == 6 {
                                        e_tmp += 0.01;
                                    }
                                    if e_tmp > ienergy {
                                        ienergy = e_tmp;
                                        cenergy = ec;
                                        intron = j;
                                    }
                                    j += 1;
                                }
                            }
                            energy += cenergy - 10.0 - ilw * (tl_state.t.nintron as f64 + 1.1 * tl_state.t.var as f64);
                            if tl_state.t.nintron >= 130 {
                                let si = se.offset(intron as isize);
                                j = *si.offset(-1);
                                if j != Guanine {
                                    if *si.offset(-2) != Adenine {
                                        energy -= 4.0;
                                    }
                                    if j != Cytosine {
                                        if j != Thymine {
                                            energy -= 8.0;
                                        }
                                    }
                                }
                            }
                        }
                        dstem = tl_state.dhit[ndx as usize].stem;
                        let dpos = tl_state.dhit[ndx as usize].pos;
                        if dstem >= 6 {
                            let sb = cpos.offset(tl_state.t.cstem as isize);
                            if *sb.offset((2 + *a1ic.get_unchecked(intron as usize) * tl_state.t.nintron) as isize) != Thymine {
                                nc += 1;
                                continue;
                            }
                            if *sb.offset((3 + *a2ic.get_unchecked(intron as usize) * tl_state.t.nintron) as isize) != Cytosine {
                                nc += 1;
                                continue;
                            }
                            if *sb.offset((4 + *a3ic.get_unchecked(intron as usize) * tl_state.t.nintron) as isize) != Adenine {
                                nc += 1;
                                continue;
                            }
                            energy += 3.0;
                        } else {
                            if (*dpos.offset(-1) & 5) == 0 {
                                i = 0;
                                let mut si = cend;
                                let se = cend.offset(4);
                                while si < se {
                                    if (*si & 5) == 0 {
                                        i += 1;
                                        if i >= 2 {
                                            energy += 3.0;
                                            break;
                                        }
                                    } else {
                                        i = 0;
                                    }
                                    si = si.offset(1);
                                }
                            }
                        }
                        if tl_state.t.cstem >= 6 {
                            let sb = cpos.offset(tl_state.t.cstem as isize);
                            if *sb.offset((2 + *a1ic.get_unchecked(intron as usize) * tl_state.t.nintron) as isize) == Cytosine {
                                if *sb.offset((3 + *a2ic.get_unchecked(intron as usize) * tl_state.t.nintron) as isize) == Thymine {
                                    if *sb.offset((4 + *a3ic.get_unchecked(intron as usize) * tl_state.t.nintron) as isize) == Adenine {
                                        energy += 4.0;
                                    }
                                }
                            }
                        }
                        if energy < ethresh {
                            nc += 1;
                            continue;
                        }
                        tl_state.t.energy = energy;
                        tl_state.t.dstem = dstem;
                        tl_state.t.astem1 = if tl_state.t.dstem < 6 { 7 } else { if tl_state.t.tstem < 5 { 9 } else { 8 } };
                        tl_state.t.astem2 = tl_state.t.astem1;
                        tl_state.t.ps = apos.offset((7 - tl_state.t.astem1) as isize);
                        tl_state.t.nbase = ((tend as isize - tl_state.t.ps as isize) / std::mem::size_of::<i32>() as isize) as i32 + tl_state.t.astem2;
                        tl_state.t.dloop = tl_state.dhit[ndx as usize].loop_;
                        tl_state.t.spacer1 = ((dpos as isize - apos as isize) / std::mem::size_of::<i32>() as isize) as i32 - 7;
                        tl_state.t.spacer2 = ((cpos as isize - tl_state.dhit[ndx as usize].end as isize) / std::mem::size_of::<i32>() as isize) as i32;
                        j = ((cpos as isize - tl_state.t.ps as isize) / std::mem::size_of::<i32>() as isize) as i32 + tl_state.t.cstem;
                        tl_state.t.anticodon = j + 2;
                        if tl_state.t.nintron > 0 {
                            tl_state.t.intron = j + intron;
                            if (tl_state.t.nbase + tl_state.t.nintron) > MAXTRNALEN_I32 {
                                nts = ti_genedetected_ts(ts, max_genes, d, nts, seq, &mut tl_state.t, sw);
                                nc += 1;
                                continue;
                            }
                        }
                        if energy < tl_state.te.energy {
                            nc += 1;
                            continue;
                        }
                        flag = 1;
                        std::ptr::copy_nonoverlapping(&tl_state.t, &mut tl_state.te, 1);
                        nc += 1;
                    }
                    cpos = cpos.offset(1);
                }
            }
            na += 1;
        }
        if flag != 0 {
            nts = ti_genedetected_ts(ts, max_genes, d, nts, seq, &mut tl_state.te, sw);
        }
        nt += 1;
    }

    nts
}

/// tmioptimise - Main tmRNA/tRNA intron detection optimizer (legacy wrapper using global TS)
/// tmrna.c:265-725
pub unsafe fn tmioptimise(
    d: *mut data_set,
    seq: *mut i32,
    lseq: i32,
    nts: i32,
    sw: *mut csw
) -> i32 {
    tmioptimise_ts(TS, (*sw).genespace, d, seq, lseq, nts, sw)
}
