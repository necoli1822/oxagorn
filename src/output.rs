/*
 * Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
 * output.rs - Display and output functions
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

use std::io::Write;
use std::ptr;
use std::ffi::CStr;

use crate::types::*;
use crate::tables::{TS, AAMAP, AANAME};
use crate::utils::{
    copy, copy3cr, softstrpos, wildstrpos, length, move_forward, seq_init
};
use crate::sequence::{
    aa, cpbase, pseudogene, gc_content, seqlen, aseqlen, name, position, location,
    nenergy, peptide_tag, disp_location, disp_seq, disp_intron, disp_tmrna_seq,
    disp_tmrna_perm_seq, disp_cds, disp_fasta_seq, disp_trna_bracket_notation,
    disp_batch_trna, disp_batch_tmrna, disp_batch_srprna, disp_batch_cds,
    init_tmrna, update_tmrna_tag_database, report_new_tmrna_tags,
    write_to_library, init_matrix, disp_gene, xcopy, disp_matrix, disp_gene_SVG,
    trna_score, tmrna_score, sense_switch
};

use crate::trna::{init_gene, init_gene_ts, overlap};
use crate::tmrna::{tmioptimise, tmioptimise_ts};
use crate::thread_state::{get_thread_local_ts_with_capacity, reset_thread_local_ts, get_thread_local_ts_capacity};
use crate::parallel::{ParallelConfig, SeqRecord, GeneResult, load_fasta_sequences, DEFAULT_MAX_THREADS};
use rayon::prelude::*;
use std::sync::atomic::{AtomicI32, Ordering};

// All constants now imported from types.rs via `use crate::types::*;`
// Available constants:
// - Nucleotides: Adenine, Cytosine, Guanine, Thymine (i32)
// - Gene types: tRNA, tmRNA, srpRNA, CDS (i32)
// - NOBASE (i32)
// - Array sizes: NS, NT, NGFT, LSEQ, WRAP (usize)
// - Lengths: MAXTAGDIST, MAXTRNALEN, MAXTMRNALEN, TSWEEP, MINCTRNALEN (various types)
// - Genetic codes: NGENECODE, METAZOAN_MT, STANDARD, VERTEBRATE_MT

// Gene storage and tables from library (no extern "C" needed)
// TS, AANAME, AAMAP are from tables.rs via use crate::types::*

/// disp_ftable_entry - Display frequency table entry
/// output.c:11-31
pub unsafe fn disp_ftable_entry(
    f: *mut File,
    n: *mut i32,
    i: i32,
    m: i32,
    sw: *mut csw
) {
    let aa_str = aa(n, sw);
    let aa_str_rust = CStr::from_ptr(aa_str as *const i8).to_string_lossy();
    if m > 0 {
        match (*sw).geneticcode {
            METAZOAN_MT => {
                let _ = write!(&mut *f, " {:<18} {:<4}", aa_str_rust, m);
            }
            STANDARD | VERTEBRATE_MT | _ => {
                let _ = write!(&mut *f, " {:<4} {:<5}", aa_str_rust, m);
            }
        }
    } else {
        match (*sw).geneticcode {
            METAZOAN_MT => {
                let _ = write!(&mut *f, " {:<18}     ", aa_str_rust);
            }
            STANDARD | VERTEBRATE_MT | _ => {
                let _ = write!(&mut *f, " {:<4}      ", aa_str_rust);
            }
        }
    }
}

/// disp_freq_table - Display frequency table
/// output.c:34-92
pub unsafe fn disp_freq_table(nt: i32, sw: *mut csw) {
    let mut i: i32;
    let mut j: i32;
    let mut k: i32;
    let mut m: i32;
    let mut ambig: i32;
    let mut s: *mut i32;
    let mut c1: i32;
    let mut c2: i32;
    let mut c3: i32;
    let mut c: [i32; 3] = [0; 3];
    let mut a: [i32; 3] = [0; 3];
    let mut table: [[[i32; 4]; 4]; 4] = [[[0; 4]; 4]; 4];

    static cgflip: [i32; 4] = [0, 2, 1, 3];
    static codonorder: [i32; 4] = [3, 1, 0, 2];

    let f = (*sw).f;
    ambig = 0;

    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                table[i as usize][j as usize][k as usize] = 0;
            }
        }
    }

    for i in 0..nt {
        if (*TS.offset(i as isize)).energy >= 0.0 {
            if (*TS.offset(i as isize)).genetype == tRNA {
                if (*TS.offset(i as isize)).cloop == 7 {
                    s = (*TS.offset(i as isize)).seq.as_mut_ptr().offset((*TS.offset(i as isize)).anticodon as isize);
                    c1 = *s;
                    c2 = *s.offset(1);
                    c3 = *s.offset(2);
                    if (c1 >= Adenine) && (c1 <= Thymine) {
                        if (c2 >= Adenine) && (c2 <= Thymine) {
                            if (c3 >= Adenine) && (c3 <= Thymine) {
                                table[*s as usize][*s.offset(1) as usize][*s.offset(2) as usize] += 1;
                            } else {
                                ambig += 1;
                            }
                        } else {
                            ambig += 1;
                        }
                    } else {
                        ambig += 1;
                    }
                } else {
                    ambig += 1;
                }
            }
        }
    }

    let _ = writeln!(&mut *f, "tRNA anticodon frequency");

    for i in 0..4 {
        c[0] = codonorder[i as usize];
        a[2] = 3 - c[0];
        for j in 0..4 {
            c[2] = codonorder[j as usize];
            a[0] = 3 - c[2];
            for k in 0..4 {
                c[1] = codonorder[k as usize];
                a[1] = 3 - c[1];
                let _ = write!(&mut *f, "{}{}{}",
                    cpbase(a[0]) as char, cpbase(a[1]) as char, cpbase(a[2]) as char);
                m = table[a[0] as usize][a[1] as usize][a[2] as usize];
                disp_ftable_entry(f, a.as_mut_ptr(), k, m, sw);
            }
            let _ = writeln!(&mut *f);
        }
        if i < 3 {
            let _ = writeln!(&mut *f);
        }
    }

    if ambig > 0 {
        let _ = writeln!(&mut *f, "Ambiguous: {}", ambig);
    }

    let _ = writeln!(&mut *f, "\ntRNA codon frequency");

    for i in 0..4 {
        c[0] = codonorder[i as usize];
        a[2] = 3 - c[0];
        for j in 0..4 {
            c[2] = codonorder[j as usize];
            a[0] = 3 - c[2];
            for k in 0..4 {
                c[1] = codonorder[k as usize];
                a[1] = 3 - c[1];
                let _ = write!(&mut *f, "{}{}{}",
                    cpbase(c[0]) as char, cpbase(c[1]) as char, cpbase(c[2]) as char);
                m = table[a[0] as usize][a[1] as usize][a[2] as usize];
                disp_ftable_entry(f, a.as_mut_ptr(), k, m, sw);
            }
            let _ = writeln!(&mut *f);
        }
        if i < 3 {
            let _ = writeln!(&mut *f);
        }
    }

    if ambig > 0 {
        let _ = writeln!(&mut *f, "Ambiguous: {}", ambig);
    }
    let _ = writeln!(&mut *f);
}

/// disp_energy_stats - Display energy statistics
/// output.c:94-163
pub unsafe fn disp_energy_stats(d: *mut data_set, nt: i32, sw: *mut csw) {
    let mut i: i32;
    let mut n: [i32; NS] = [0; NS];
    let mut genetype: i32;
    let mut introns: i32;
    let mut nintron: i32 = 0;
    let mut trna_flag: i32;
    let mut mtrna: i32;
    let mut ntv: i32 = 0;
    let mut nd: i32 = 0;
    let mut nps: i32;
    let mut gc: f64;
    let mut gcmin: [f64; NS] = [0.0; NS];
    let mut gcmax: [f64; NS] = [0.0; NS];

    let f = (*sw).f;
    mtrna = (*sw).mtrna;
    trna_flag = (*sw).trna | mtrna;
    nps = 0;

    if mtrna != 0 {
        ntv = 0;
        nd = 0;
    }

    if ((*sw).trna != 0) && ((*sw).maxintronlen > 0) {
        introns = 1;
        nintron = 0;
    } else {
        introns = 0;
    }

    for i in 0..NS as i32 {
        n[i as usize] = 0;
        gcmin[i as usize] = 1.0;
        gcmax[i as usize] = 0.0;
    }

    for i in 0..nt {
        if (*TS.offset(i as isize)).energy >= 0.0 {
            n[NS - 1] += 1;
            genetype = (*TS.offset(i as isize)).genetype;
            n[genetype as usize] += 1;
            if pseudogene(TS.offset(i as isize), sw) != 0 {
                nps += 1;
            }
            if genetype == tRNA {
                if mtrna != 0 {
                    if (*TS.offset(i as isize)).tstem == 0 {
                        ntv += 1;
                    }
                    if (*TS.offset(i as isize)).dstem == 0 {
                        nd += 1;
                    }
                }
                if introns != 0 {
                    if (*TS.offset(i as isize)).nintron > 0 {
                        nintron += 1;
                    }
                }
                gc = gc_content(TS.offset(i as isize));
                if gc < gcmin[genetype as usize] {
                    gcmin[genetype as usize] = gc;
                }
                if gc > gcmax[genetype as usize] {
                    gcmax[genetype as usize] = gc;
                }
            }
        }
    }

    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f);

    if (*sw).repeatsn != 0 {
        if (n[tRNA as usize] + n[tmRNA as usize]) > 0 {
            let seqname_rust = CStr::from_ptr((*d).seqname.as_ptr() as *const i8).to_string_lossy();
            let _ = writeln!(&mut *f, "{}\n", seqname_rust);
        }
    }

    if trna_flag != 0 {
        (*sw).ngene[tRNA as usize] += n[tRNA as usize];
        if n[tRNA as usize] > 3 {
            disp_freq_table(nt, sw);
        }
        if (n[tRNA as usize] > 1) || (((*sw).tmrna != 0) && (n[tmRNA as usize] > 0)) {
            if introns != 0 {
                if (*sw).minintronlen == 0 {
                    let _ = writeln!(&mut *f, "Number of tRNA genes with no introns = {}",
                        n[0] - nintron);
                }
                let _ = writeln!(&mut *f, "Number of tRNA genes with C-loop introns = {}",
                    nintron);
            } else {
                let genetypename_rust = CStr::from_ptr((*sw).genetypename[tRNA as usize].as_ptr() as *const i8).to_string_lossy();
                let _ = writeln!(&mut *f, "Number of {} genes = {}",
                    genetypename_rust, n[tRNA as usize]);
            }
            if mtrna != 0 {
                if (*sw).tvloop != 0 {
                    let _ = writeln!(&mut *f, "Number of TV replacement loop tRNA genes = {}",
                        ntv);
                }
                let _ = writeln!(&mut *f, "Number of D replacement loop tRNA genes = {}",
                    nd);
            }
            if n[tRNA as usize] > 1 {
                let _ = writeln!(&mut *f, "tRNA GC range = {:.1}% to {:.1}%",
                    gcmin[0] * 100.0, gcmax[0] * 100.0);
            }
        }
    }

    if (*sw).tmrna != 0 {
        (*sw).ngene[tmRNA as usize] += n[tmRNA as usize];
        if (n[tmRNA as usize] > 1) || (trna_flag != 0 && (n[tRNA as usize] > 0)) {
            let genetypename_rust = CStr::from_ptr((*sw).genetypename[tmRNA as usize].as_ptr() as *const i8).to_string_lossy();
            let _ = writeln!(&mut *f, "Number of {} genes = {}",
                genetypename_rust, n[tmRNA as usize]);
        }
    }

    (*sw).nps += nps;
    if (*sw).reportpseudogenes != 0 {
        if nps > 0 {
            if n[NS - 1] > 1 {
                let _ = writeln!(&mut *f, "Number of possible pseudogenes = {}", nps);
            }
        }
    }

    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f);
}

/// batch_energy_stats - Batch energy statistics
/// output.c:166-201
pub unsafe fn batch_energy_stats(d: *mut data_set, nt: i32, sw: *mut csw) {
    let mut i: i32;
    let mut n: [i32; NS] = [0; NS];
    let mut genetype: i32;
    let mut introns: i32;
    let mut nintron: i32 = 0;
    let mut trna_flag: i32;
    let mut mtrna: i32;
    let mut ntv: i32 = 0;
    let mut nd: i32 = 0;
    let mut nps: i32;
    let mut gc: f64;
    let mut gcmin: [f64; NS] = [0.0; NS];
    let mut gcmax: [f64; NS] = [0.0; NS];

    let f = (*sw).f;
    mtrna = (*sw).mtrna;
    trna_flag = (*sw).trna | mtrna;
    nps = 0;

    if mtrna != 0 {
        ntv = 0;
        nd = 0;
    }

    if ((*sw).trna != 0) && ((*sw).maxintronlen > 0) {
        introns = 1;
        nintron = 0;
    } else {
        introns = 0;
    }

    for i in 0..NS as i32 {
        n[i as usize] = 0;
        gcmin[i as usize] = 1.0;
        gcmax[i as usize] = 0.0;
    }

    for i in 0..nt {
        if (*TS.offset(i as isize)).energy >= 0.0 {
            n[NS - 1] += 1;
            genetype = (*TS.offset(i as isize)).genetype;
            n[genetype as usize] += 1;
            if (*TS.offset(i as isize)).energy < 100.0 {
                nps += 1;
            }
            if genetype == tRNA {
                if mtrna != 0 {
                    if (*TS.offset(i as isize)).tstem == 0 {
                        ntv += 1;
                    }
                    if (*TS.offset(i as isize)).dstem == 0 {
                        nd += 1;
                    }
                }
                if introns != 0 {
                    if (*TS.offset(i as isize)).nintron > 0 {
                        nintron += 1;
                    }
                }
                gc = gc_content(TS.offset(i as isize));
                if gc < gcmin[genetype as usize] {
                    gcmin[genetype as usize] = gc;
                }
                if gc > gcmax[genetype as usize] {
                    gcmax[genetype as usize] = gc;
                }
            }
        }
    }

    if trna_flag != 0 {
        (*sw).ngene[tRNA as usize] += n[tRNA as usize];
    }
    if (*sw).tmrna != 0 {
        (*sw).ngene[tmRNA as usize] += n[tmRNA as usize];
    }
    (*sw).nps += nps;
}

/// gene_sort - Sort genes by position
/// output.c:204-240
pub unsafe fn gene_sort(d: *mut data_set, nt: i32, sort: *mut i32, sw: *mut csw) -> i32 {
    let mut i: i32;
    let mut n: i32;
    let mut j: i32;
    let mut k: i32;
    let mut starti: i64;
    let mut startj: i64;
    let mut stopi: i64;
    let mut stopj: i64;
    let psmax: i64;

    psmax = (*d).psmax;
    n = 0;

    for i in 0..nt {
        if (*TS.offset(i as isize)).energy >= 0.0 {
            if (*sw).ireportminintronlen == 1 {
                if (*TS.offset(i as isize)).genetype == tRNA {
                    if (*TS.offset(i as isize)).nintron < (*sw).minintronlenreport {
                        continue;
                    }
                }
            }
            *sort.offset(n as isize) = i;
            n += 1;
        }
    }

    i = -1;
    i += 1;
    while i < (n - 1) {
        j = i;
        j += 1;
        while j < n {
            starti = (*TS.offset(*sort.offset(i as isize) as isize)).start;
            startj = (*TS.offset(*sort.offset(j as isize) as isize)).start;
            stopi = (*TS.offset(*sort.offset(i as isize) as isize)).stop;
            stopj = (*TS.offset(*sort.offset(j as isize) as isize)).stop;

            if stopi < starti {
                if (psmax - starti) < stopi {
                    starti -= psmax;
                } else {
                    stopi += psmax;
                }
            }
            if stopj < startj {
                if (psmax - startj) < stopj {
                    startj -= psmax;
                } else {
                    stopj += psmax;
                }
            }

            if starti > startj {
                k = *sort.offset(i as isize);
                *sort.offset(i as isize) = *sort.offset(j as isize);
                *sort.offset(j as isize) = k;
            } else if starti == startj {
                if stopi < stopj {
                    k = *sort.offset(i as isize);
                    *sort.offset(i as isize) = *sort.offset(j as isize);
                    *sort.offset(j as isize) = k;
                }
            }
            j += 1;
        }
        i += 1;
    }

    n
}

/// iamatch - Check isoacceptor match
/// output.c:243-252
pub unsafe fn iamatch(d: *mut data_set, t: *mut Gene, sw: *mut csw) -> i32 {
    let mut key: [u8; 5] = [0; 5];
    let mut k: *mut u8;
    let mut s: [u8; 100] = [0; 100];

    k = softstrpos((*d).seqname.as_ptr() as *mut u8, b"TRNA-\0".as_ptr() as *mut u8);
    if !k.is_null() {
        k = k.offset(5);
    } else {
        k = wildstrpos((*d).seqname.as_ptr() as *mut u8, b"|***|\0".as_ptr() as *mut u8);
        if !k.is_null() {
            k = k.offset(1);
        } else {
            return -1;
        }
    }

    copy3cr(k, key.as_mut_ptr(), 3);
    name(t, s.as_mut_ptr(), 1, sw);

    if !softstrpos(s.as_ptr() as *mut u8, key.as_ptr() as *mut u8).is_null() {
        return 1;
    }
    0
}

/// gene_mismatch - Check gene mismatch
/// output.c:256-284
pub unsafe fn gene_mismatch(d: *mut data_set, agene: *mut annotated_gene, t: *mut Gene, sw: *mut csw) -> i32 {
    let mut w: i32;
    let mut alen: i32;
    let mut dlen: i32;
    let mut s: *const u8;

    w = 0;
    dlen = seqlen(t);
    alen = aseqlen(d, agene);

    match (*t).genetype {
        x if x == tRNA => {
            s = aa((*t).seq.as_mut_ptr().offset((*t).anticodon as isize), sw);
            if softstrpos(s as *mut u8, (*agene).species.as_ptr().offset(5) as *mut u8).is_null() {
                if (*t).cloop == 8 {
                    s = aa((*t).seq.as_mut_ptr().offset(((*t).anticodon + 1) as isize), sw);
                    if softstrpos(s as *mut u8, (*agene).species.as_ptr().offset(5) as *mut u8).is_null() {
                        w += 1;
                    }
                } else if (*t).cloop == 6 {
                    s = aa((*t).seq.as_mut_ptr().offset(((*t).anticodon - 1) as isize), sw);
                    if softstrpos(s as *mut u8, (*agene).species.as_ptr().offset(5) as *mut u8).is_null() {
                        w += 1;
                    }
                } else {
                    w += 1;
                }
            }
            if (*agene).comp != (*t).comp {
                w += 2;
            }
            if alen <= (dlen - (*sw).trnalenmisthresh) {
                w += 4;
            } else if alen >= (dlen + (*sw).trnalenmisthresh) {
                w += 4;
            }
        }
        x if x == tmRNA => {
            if (*agene).comp != (*t).comp {
                w += 2;
            }
            if alen <= (dlen - (*sw).tmrnalenmisthresh) {
                w += 4;
            } else if alen >= (dlen + (*sw).tmrnalenmisthresh) {
                w += 4;
            }
        }
        _ => {}
    }

    w
}

/// gene_mismatch_report - Generate mismatch report
/// output.c:287-305
pub unsafe fn gene_mismatch_report(
    d: *mut data_set,
    agene: *mut annotated_gene,
    t: *mut Gene,
    report: *mut u8,
    sw: *mut csw
) -> i32 {
    let mut w: i32;
    let mut s: *mut u8;

    w = gene_mismatch(d, agene, t, sw);
    s = report;

    if (w & 1) != 0 {
        s = copy(b"amino acceptor\0".as_ptr() as *mut u8, s);
    }
    if (w & 2) != 0 {
        if (w & 1) != 0 {
            if (w & 4) != 0 {
                s = copy(b", \0".as_ptr() as *mut u8, s);
            } else {
                s = copy(b" and \0".as_ptr() as *mut u8, s);
            }
        }
        s = copy(b"sense\0".as_ptr() as *mut u8, s);
    }
    if (w & 4) != 0 {
        if (w & 3) > 0 {
            s = copy(b" and \0".as_ptr() as *mut u8, s);
        }
        s = copy(b"sequence length\0".as_ptr() as *mut u8, s);
    }
    if w > 0 {
        s = copy(b" mismatch\0".as_ptr() as *mut u8, s);
    }
    *s = 0;

    w
}

/// nearest_annotated_gene - Find nearest annotated gene
/// output.c:309-381
pub unsafe fn nearest_annotated_gene(
    d: *mut data_set,
    t: *mut Gene,
    list: *mut i32,
    score: *mut i32,
    nmax: i32,
    sw: *mut csw
) -> i32 {
    let mut n: i32;
    let mut i: i32;
    let mut j: i32;
    let mut k: i32;
    let mut q: i32;
    let mut w: i32;
    let mut nagene: i32;
    let mut a: i64;
    let mut b: i64;
    let mut c: i64;
    let mut e: i64;
    let mut thresh: i64;
    let psmax: i64;
    let mut s: *mut u8;
    let mut ta: *mut annotated_gene;

    psmax = (*d).psmax;
    nagene = (*d).nagene[NS - 1];
    ta = (*d).gene.as_mut_ptr();
    n = 0;
    a = (*t).start;
    b = (*t).stop;
    thresh = b - a;

    if b < a {
        b += psmax;
        thresh += psmax;
        for i in 0..nagene {
            c = (*ta.offset(i as isize)).start;
            e = (*ta.offset(i as isize)).stop;
            if e < c {
                e += psmax;
                if !(a > e) && !(b < c) {
                    if n >= nmax {
                        break;
                    }
                    *list.offset(n as isize) = i;
                    *score.offset(n as isize) = if a >= c {
                        if b >= e { (e - a) as i32 } else { thresh as i32 }
                    } else {
                        if b >= e { (e - c) as i32 } else { (b - c) as i32 }
                    };
                    n += 1;
                }
                c -= psmax;
                e -= psmax;
            }
            if a > e { continue; }
            if b < c { continue; }
            if n >= nmax { break; }
            *list.offset(n as isize) = i;
            *score.offset(n as isize) = if a >= c {
                if b >= e { (e - a) as i32 } else { thresh as i32 }
            } else {
                if b >= e { (e - c) as i32 } else { (b - c) as i32 }
            };
            n += 1;
        }
        a -= psmax;
        b -= psmax;
    }

    for i in 0..nagene {
        c = (*ta.offset(i as isize)).start;
        e = (*ta.offset(i as isize)).stop;
        if e < c {
            e += psmax;
            if !(a > e) && !(b < c) {
                if n >= nmax { break; }
                *list.offset(n as isize) = i;
                *score.offset(n as isize) = if a >= c {
                    if b >= e { (e - a) as i32 } else { thresh as i32 }
                } else {
                    if b >= e { (e - c) as i32 } else { (b - c) as i32 }
                };
                n += 1;
            }
            c -= psmax;
            e -= psmax;
        }
        if a > e { continue; }
        if b < c { continue; }
        if n >= nmax { break; }
        *list.offset(n as isize) = i;
        *score.offset(n as isize) = if a >= c {
            if b >= e { (e - a) as i32 } else { thresh as i32 }
        } else {
            if b >= e { (e - c) as i32 } else { (b - c) as i32 }
        };
        n += 1;
    }

    for i in 0..n {
        k = *list.offset(i as isize);
        if (*ta.offset(k as isize)).genetype == (*t).genetype {
            *score.offset(i as isize) += 5000;
            w = gene_mismatch(d, ta.offset(k as isize), t, sw);
            if (w & 1) != 0 {
                *score.offset(i as isize) -= 2;
            }
            if (w & 2) != 0 {
                *score.offset(i as isize) -= 1;
            }
        }
    }

    if n > 1 {
        for i in 0..(n - 1) {
            for j in (i + 1)..n {
                if *score.offset(j as isize) > *score.offset(i as isize) {
                    k = *list.offset(i as isize);
                    *list.offset(i as isize) = *list.offset(j as isize);
                    *list.offset(j as isize) = k;
                    k = *score.offset(i as isize);
                    *score.offset(i as isize) = *score.offset(j as isize);
                    *score.offset(j as isize) = k;
                }
            }
        }
    }

    n
}

/// proximity_compare - Compare proximity scores
/// output.c:386-413
pub unsafe fn proximity_compare(
    d: *mut data_set,
    is: i32,
    prox: i64,
    dlen: i64,
    alen: i64,
    a: *mut annotated_gene,
    sw: *mut csw
) -> i32 {
    let mut w: i32;
    let mut score: i32;
    let mut diff: i64;
    let mut nm: [u8; 200] = [0; 200];
    let t: *mut Gene;

    t = TS.offset(is as isize);
    w = gene_mismatch(d, a, t, sw);

    if prox >= alen {
        diff = dlen - alen;
        if prox >= (2 * diff) {
            score = (prox - diff) as i32;
        } else {
            score = (prox / 2) as i32;
        }
    } else if prox >= dlen {
        diff = alen - dlen;
        if prox >= (2 * diff) {
            score = (prox - diff) as i32;
        } else {
            score = (prox / 2) as i32;
        }
    } else {
        score = prox as i32;
    }

    if (w & 1) != 0 { score -= 10; }
    if (w & 2) != 0 { score -= 2; }
    if score < 0 { score = 0; }

    if (*t).annotation >= 0 {
        if (*t).annosc >= score {
            return -1;
        }
    }

    score
}

/// nearest_detected_gene - Find nearest detected gene
/// output.c:418-482
pub unsafe fn nearest_detected_gene(
    d: *mut data_set,
    sort: *mut i32,
    nd: i32,
    scorep: *mut i32,
    ag: *mut annotated_gene,
    sw: *mut csw
) -> i32 {
    let mut n: i32;
    let mut i: i32;
    let mut is: i32;
    let mut a: i64;
    let mut b: i64;
    let mut c: i64;
    let mut e: i64;
    let mut score: i64;
    let mut alen: i64;
    let mut scoremax: i64;
    let psmax: i64;
    let mut prox: i64;
    let mut proximity: i64;

    psmax = (*d).psmax;
    n = -1;
    scoremax = -1;
    a = (*ag).start;
    b = (*ag).stop;
    alen = b - a;
    if b < a { alen += psmax; }
    proximity = 1 + alen / 2;
    if proximity > 30 { proximity = 30; }

    if b < a {
        b += psmax;
        for i in 0..nd {
            is = *sort.offset(i as isize);
            if (*ag).genetype != (*TS.offset(is as isize)).genetype { continue; }
            c = (*TS.offset(is as isize)).start;
            e = (*TS.offset(is as isize)).stop;
            if e < c {
                e += psmax;
                if !(a > e) && !(b < c) {
                    prox = if a >= c {
                        if b >= e { e - a } else { alen }
                    } else {
                        if b >= e { e - c } else { b - c }
                    };
                    if prox >= proximity {
                        score = proximity_compare(d, is, prox, e - c, alen, ag, sw) as i64;
                        if score > scoremax {
                            n = i;
                            scoremax = score;
                        }
                    }
                }
                c -= psmax;
                e -= psmax;
            }
            if a > e { continue; }
            if b < c { continue; }
            prox = if a >= c {
                if b >= e { e - a } else { alen }
            } else {
                if b >= e { e - c } else { b - c }
            };
            if prox >= proximity {
                score = proximity_compare(d, is, prox, e - c, alen, ag, sw) as i64;
                if score > scoremax {
                    n = i;
                    scoremax = score;
                }
            }
        }
        a -= psmax;
        b -= psmax;
    }

    for i in 0..nd {
        is = *sort.offset(i as isize);
        if (*ag).genetype != (*TS.offset(is as isize)).genetype { continue; }
        c = (*TS.offset(is as isize)).start;
        e = (*TS.offset(is as isize)).stop;
        if e < c {
            e += psmax;
            if !(a > e) && !(b < c) {
                prox = if a >= c {
                    if b >= e { e - a } else { alen }
                } else {
                    if b >= e { e - c } else { b - c }
                };
                if prox >= proximity {
                    score = proximity_compare(d, is, prox, e - c, alen, ag, sw) as i64;
                    if score > scoremax {
                        n = is;
                        scoremax = score;
                    }
                }
            }
            c -= psmax;
            e -= psmax;
        }
        if a > e { continue; }
        if b < c { continue; }
        prox = if a >= c {
            if b >= e { e - a } else { alen }
        } else {
            if b >= e { e - c } else { b - c }
        };
        if prox >= proximity {
            score = proximity_compare(d, is, prox, e - c, alen, ag, sw) as i64;
            if score > scoremax {
                n = is;
                scoremax = score;
            }
        }
    }

    *scorep = scoremax as i32;
    n
}

/// disp_match - Display gene matches
/// output.c:486-642
pub unsafe fn disp_match(d: *mut data_set, sort: *mut i32, nd: i32, sw: *mut csw) {
    let mut i: i32;
    let ld: i32;
    let mut fn_arr: [i32; NS] = [0; NS];
    let mut fp: [i32; NS] = [0; NS];
    let mut fpd: i32 = 0;
    let mut fptv: i32 = 0;
    let mut w: i32 = 0;
    let mut score: i32 = 0;
    let mut detect: i32;
    let mut n: [i32; NS] = [0; NS];
    let mut prevannoted: i32;
    let mut nl: i32;
    let mut k: i32;
    let mut csort_arr: [i32; NGFT] = [0; NGFT];
    let mut msort: *mut i32;
    let mut start: i64;
    let mut tag: [u8; 52] = [0; 52];
    let mut nm: [u8; 100] = [0; 100];
    let mut anm: [u8; 100] = [0; 100];
    let mut ps: [u8; 100] = [0; 100];
    let mut mreport: [u8; 100] = [0; 100];
    let f = (*sw).f;
    let mut t: *mut Gene = std::ptr::null_mut();
    let mut agene: *mut annotated_gene;
    let mut a: *mut annotated_gene;
    static gp: [[u8; 7]; 2] = [*b"genes\0\0", *b"gene\0\0\0"];
    static comp: [u8; 3] = [b' ', b'c', 0];
    static aps: [[u8; 5]; 2] = [*b"  \0\0\0", *b"PS\0\0\0"];

    // sq macro inline: #define sq(pos)  ((pos + d->psmax - 1L) % d->psmax) + 1L
    let sq = |pos: i64| -> i64 {
        ((pos + (*d).psmax - 1) % (*d).psmax) + 1
    };

    nl = nd;
    if ((*sw).trna != 0) | ((*sw).mtrna != 0) { nl += (*d).nagene[tRNA as usize]; }
    if (*sw).tmrna != 0 { nl += (*d).nagene[tmRNA as usize]; }

    let mut msort_vec: Option<Vec<i32>> = None;
    if nl < NGFT as i32 {
        msort = csort_arr.as_mut_ptr();
    } else {
        let vec = vec![0i32; nl as usize];
        msort = vec.as_ptr() as *mut i32;
        msort_vec = Some(vec);
    }

    let seqname_rust = CStr::from_ptr((*d).seqname.as_ptr() as *const i8).to_string_lossy();
    let _ = writeln!(&mut *f, "\n{}", seqname_rust);
    let _ = writeln!(&mut *f, "{} nucleotides in sequence", (*d).psmax);
    let _ = writeln!(&mut *f, "Mean G+C content = {:.1}%", 100.0 * (*d).gc);
    let _ = writeln!(&mut *f, "\nGenBank to Aragorn comparison\n");

    (*sw).dispmatch = 1;

    i = 0;
    while i < nd {
        w = *sort.offset(i as isize);
        if (*TS.offset(w as isize)).energy >= 0.0 {
            n[NS - 1] += 1;
            n[(*TS.offset(w as isize)).genetype as usize] += 1;
        }
        (*TS.offset(w as isize)).annotation = -1;
        (*TS.offset(w as isize)).annosc = -1;
        i += 1;
    }

    if ((*sw).trna != 0) | ((*sw).mtrna != 0) | ((*sw).tmrna != 0) {
        fpd = 0;
        fptv = 0;
        if ((*sw).trna != 0) | ((*sw).mtrna != 0) {
            let gp_str = CStr::from_ptr(gp[if (*d).nagene[tRNA as usize] == 1 { 1 } else { 0 }].as_ptr() as *const i8).to_string_lossy();
            let _ = writeln!(&mut *f, "{} annotated tRNA {}",
                (*d).nagene[tRNA as usize], gp_str);
            let gp_str2 = CStr::from_ptr(gp[if n[tRNA as usize] == 1 { 1 } else { 0 }].as_ptr() as *const i8).to_string_lossy();
            let _ = writeln!(&mut *f, "{} detected tRNA {}",
                n[tRNA as usize], gp_str2);
        }
        if (*sw).tmrna != 0 {
            let gp_str = CStr::from_ptr(gp[if (*d).nagene[tmRNA as usize] == 1 { 1 } else { 0 }].as_ptr() as *const i8).to_string_lossy();
            let _ = writeln!(&mut *f, "{} annotated tmRNA {}",
                (*d).nagene[tmRNA as usize], gp_str);
            let gp_str2 = CStr::from_ptr(gp[if n[tmRNA as usize] == 1 { 1 } else { 0 }].as_ptr() as *const i8).to_string_lossy();
            let _ = writeln!(&mut *f, "{} detected tmRNA {}",
                n[tmRNA as usize], gp_str2);
        }
        let _ = writeln!(&mut *f, "\n  GenBank                                      Aragorn");

        nl = 0;
        i = 0;
        while i < (*d).nagene[NS - 1] {
            agene = (*d).gene.as_mut_ptr().offset(i as isize);
            (*agene).detected = -1;
            if (*agene).genetype != tRNA as i32 {
                if (*agene).genetype != tmRNA as i32 {
                    i += 1;
                    continue;
                } else if (*sw).tmrna == 0 {
                    i += 1;
                    continue;
                }
            } else if (*sw).trna == 0 {
                if (*sw).mtrna == 0 {
                    i += 1;
                    continue;
                }
            }
            a = agene;
            k = i;
            loop {
                (*a).detected = nearest_detected_gene(d, sort, nd, &mut score, a, sw);
                if (*a).detected < 0 { break; }
                t = TS.offset((*a).detected as isize);
                prevannoted = (*t).annotation;
                (*t).annotation = k;
                (*t).annosc = score;
                if prevannoted < 0 { break; }
                if prevannoted == k { break; }
                if prevannoted == i { break; }
                a = (*d).gene.as_mut_ptr().offset(prevannoted as isize);
                k = prevannoted;
            }
            k = nl;
            loop {
                k -= 1;
                if k < 0 { break; }
                if (*agene).start >= (*(*d).gene.as_mut_ptr().offset(*msort.offset(k as isize) as isize)).start { break; }
                *msort.offset((k + 1) as isize) = *msort.offset(k as isize);
            }
            k += 1;
            *msort.offset(k as isize) = i;
            nl += 1;
            i += 1;
        }

        i = 0;
        while i < nd {
            t = TS.offset(*sort.offset(i as isize) as isize);
            if (*t).annotation >= 0 {
                i += 1;
                continue;
            }
            if (*t).genetype != tRNA {
                if (*t).genetype != tmRNA {
                    i += 1;
                    continue;
                } else if (*sw).tmrna == 0 {
                    i += 1;
                    continue;
                }
            } else if (*sw).trna == 0 {
                if (*sw).mtrna == 0 {
                    i += 1;
                    continue;
                }
            }
            k = nl;
            loop {
                k -= 1;
                if k < 0 { break; }
                if *msort.offset(k as isize) >= 0 {
                    start = (*(*d).gene.as_mut_ptr().offset(*msort.offset(k as isize) as isize)).start;
                } else {
                    start = (*TS.offset((-1 - *msort.offset(k as isize)) as isize)).start;
                }
                if (*t).start >= start { break; }
                *msort.offset((k + 1) as isize) = *msort.offset(k as isize);
            }
            k += 1;
            *msort.offset(k as isize) = -(*sort.offset(i as isize) + 1);
            nl += 1;
            i += 1;
        }

        i = 0;
        while i < nl {
            if *msort.offset(i as isize) >= 0 {
                agene = (*d).gene.as_mut_ptr().offset(*msort.offset(i as isize) as isize);
                detect = (*agene).detected;
                if detect >= 0 {
                    t = TS.offset(detect as isize);
                    w = gene_mismatch_report(d, agene, t, mreport.as_mut_ptr(), sw);
                    if w > 0 {
                        let _ = write!(&mut *f, "*");
                    } else {
                        let _ = write!(&mut *f, " ");
                    }
                } else {
                    let _ = write!(&mut *f, "*");
                }
                let species_rust = CStr::from_ptr((*agene).species.as_ptr() as *const i8).to_string_lossy();
                let aps_rust = CStr::from_ptr(aps[(*agene).pseudogene as usize].as_ptr() as *const i8).to_string_lossy();
                let anm_str = format!(" {:<11}{}({},{}) {}",
                    species_rust,
                    comp[(*agene).comp as usize] as char,
                    sq((*agene).start),
                    sq((*agene).stop),
                    aps_rust);
                let _ = write!(&mut *f, "{:<45} ", anm_str);
                if detect >= 0 {
                    let nm_rust = CStr::from_ptr(name(t, nm.as_mut_ptr(), 1, sw) as *const i8).to_string_lossy();
                    let _ = write!(&mut *f, "{} ", nm_rust);
                    if (*t).comp == 0 { let _ = write!(&mut *f, " "); }
                    let ps_rust = CStr::from_ptr(position(ps.as_mut_ptr(), t, sw) as *const i8).to_string_lossy();
                    let _ = write!(&mut *f, "{}", ps_rust);
                    if (*sw).energydisp != 0 {
                        let _ = write!(&mut *f, " {:7.3}", (*t).energy);
                    }
                    if (*t).genetype == tmRNA {
                        peptide_tag(tag.as_mut_ptr(), 50, t, sw);
                        let tag_rust = CStr::from_ptr(tag.as_ptr() as *const i8).to_string_lossy();
                        let _ = write!(&mut *f, " {}", tag_rust);
                    }
                    if (*sw).reportpseudogenes != 0 {
                        if pseudogene(t, sw) != 0 {
                            let _ = write!(&mut *f, " PS");
                        }
                    }
                    if w > 0 {
                        let mreport_rust = CStr::from_ptr(mreport.as_ptr() as *const i8).to_string_lossy();
                        let _ = write!(&mut *f, " {}", mreport_rust);
                    }
                    let _ = writeln!(&mut *f);
                } else {
                    let _ = writeln!(&mut *f, "Not detected");
                    fn_arr[(*agene).genetype as usize] += 1;
                }
            } else {
                t = TS.offset((-(*msort.offset(i as isize) + 1)) as isize);
                let nm_rust = CStr::from_ptr(name(t, nm.as_mut_ptr(), 1, sw) as *const i8).to_string_lossy();
                let _ = write!(&mut *f, "* Not annotated                                {} ", nm_rust);
                if (*t).comp == 0 { let _ = write!(&mut *f, " "); }
                let ps_rust = CStr::from_ptr(position(ps.as_mut_ptr(), t, sw) as *const i8).to_string_lossy();
                let _ = write!(&mut *f, "{}", ps_rust);
                if (*sw).energydisp != 0 {
                    let _ = write!(&mut *f, " {:7.3}", (*t).energy);
                }
                if (*t).genetype == tmRNA {
                    peptide_tag(tag.as_mut_ptr(), 50, t, sw);
                    let tag_rust = CStr::from_ptr(tag.as_ptr() as *const i8).to_string_lossy();
                    let _ = write!(&mut *f, " {}", tag_rust);
                }
                if (*sw).reportpseudogenes != 0 {
                    if pseudogene(t, sw) != 0 {
                        let _ = write!(&mut *f, " PS");
                    }
                }
                let _ = writeln!(&mut *f);
                fp[(*t).genetype as usize] += 1;
                if (*t).genetype == tRNA {
                    if (*t).dstem == 0 { fpd += 1; }
                    if (*t).tstem == 0 { fptv += 1; }
                }
            }
            i += 1;
        }

        let _ = writeln!(&mut *f);
        if ((*sw).trna != 0) | ((*sw).mtrna != 0) {
            let _ = writeln!(&mut *f, "Number of annotated tRNA genes not detected = {}", fn_arr[tRNA as usize]);
            let _ = writeln!(&mut *f, "Number of unannotated tRNA genes detected = {}", fp[tRNA as usize]);
        }
        if (*sw).mtrna != 0 {
            let _ = writeln!(&mut *f, "Number of unannotated D-replacement tRNA genes detected = {}", fpd);
            let _ = writeln!(&mut *f, "Number of unannotated TV-replacement tRNA genes detected = {}", fptv);
        }
        if (*sw).tmrna != 0 {
            let _ = writeln!(&mut *f, "Number of annotated tmRNA genes not detected = {}", fn_arr[tmRNA as usize]);
            let _ = writeln!(&mut *f, "Number of unannotated tmRNA genes detected = {}", fp[tmRNA as usize]);
        }
        let _ = writeln!(&mut *f, "\n");

        i = tRNA as i32;
        while i <= tmRNA as i32 {
            (*sw).nagene[i as usize] += (*d).nagene[i as usize];
            (*sw).nafn[i as usize] += fn_arr[i as usize];
            (*sw).nafp[i as usize] += fp[i as usize];
            i += 1;
        }
        if (*sw).mtrna != 0 {
            (*sw).natfpd += fpd;
            (*sw).natfptv += fptv;
        }
    }

    (*sw).nabase += (*d).psmax;
    (*sw).dispmatch = 0;
    // msort_vec is automatically dropped here
    drop(msort_vec);
}

/// annotation_overlap_check - Check annotation overlap
/// output.c:645-678
pub unsafe fn annotation_overlap_check(
    d: *mut data_set,
    t: *mut Gene,
    margin: *const u8,
    sw: *mut csw
) {
    let mut a: i32;
    let mut m: i32;
    let mut n: i32;
    let mut w: i32;
    let mut list: [i32; 20] = [0; 20];
    let mut score: [i32; 20] = [0; 20];
    let mut mreport: [u8; 100] = [0; 100];
    static comp: [u8; 3] = [b' ', b'c', 0];

    n = nearest_annotated_gene(d, t, list.as_mut_ptr(), score.as_mut_ptr(), 20, sw);

    if n < 1 {
        m = -1;
    } else {
        m = 0;
        a = list[m as usize];
        if (*(*d).gene.as_mut_ptr().offset(a as isize)).genetype != (*t).genetype {
            m = -1;
        } else {
            w = gene_mismatch_report(d, (*d).gene.as_mut_ptr().offset(a as isize), t, mreport.as_mut_ptr(), sw);
            if (w & 1) != 0 {
                if (score[m as usize] - 5000) < (3 * seqlen(t) / 4) {
                    m = -1;
                }
            } else {
                if (score[m as usize] - 5000) < (seqlen(t) / 3) {
                    m = -1;
                }
            }
        }
    }

    let margin_rust = CStr::from_ptr(margin as *const i8).to_string_lossy();
    if m < 0 {
        let _ = writeln!(&mut *(*sw).f, "{}Not annotated", margin_rust);
    } else {
        a = list[m as usize];
        let species_rust = CStr::from_ptr((*(*d).gene.as_mut_ptr().offset(a as isize)).species.as_ptr() as *const i8).to_string_lossy();
        let _ = write!(&mut *(*sw).f, "{}Match with annotated {} {}({},{})",
            margin_rust,
            species_rust,
            comp[(*(*d).gene.as_mut_ptr().offset(a as isize)).comp as usize] as char,
            (*(*d).gene.as_mut_ptr().offset(a as isize)).start,
            (*(*d).gene.as_mut_ptr().offset(a as isize)).stop);
        w = gene_mismatch_report(d, (*d).gene.as_mut_ptr().offset(a as isize), t, mreport.as_mut_ptr(), sw);
        if w > 0 {
            let mreport_rust = CStr::from_ptr(mreport.as_ptr() as *const i8).to_string_lossy();
            let _ = write!(&mut *(*sw).f, " * {}", mreport_rust);
        }
        let _ = writeln!(&mut *(*sw).f);
    }

    m += 1;
    while m < n {
        a = list[m as usize];
        let species_rust = CStr::from_ptr((*(*d).gene.as_mut_ptr().offset(a as isize)).species.as_ptr() as *const i8).to_string_lossy();
        let _ = writeln!(&mut *(*sw).f, "{}Overlap with annotated {} {}({},{})",
            margin_rust,
            species_rust,
            comp[(*(*d).gene.as_mut_ptr().offset(a as isize)).comp as usize] as char,
            (*(*d).gene.as_mut_ptr().offset(a as isize)).start,
            (*(*d).gene.as_mut_ptr().offset(a as isize)).stop);
        m += 1;
    }

    let _ = writeln!(&mut *(*sw).f);
}

/// disp_gene_set - Display gene set
/// output.c:680-769
pub unsafe fn disp_gene_set(d: *mut data_set, nt: i32, sw: *mut csw) {
    let mut i: i32;
    let mut j: i32;
    let mut n: i32;
    let mut vsort: [i32; NT] = [0; NT];
    let mut sort: *mut i32;
    let mut m: [[u8; MATY]; MATX] = [[0; MATY]; MATX];
    let mut s: [u8; 20] = [0; 20];
    let mut t: *mut Gene;
    let f = (*sw).f;

    let mut sort_vec: Option<Vec<i32>> = None;
    if nt <= NT as i32 {
        sort = vsort.as_mut_ptr();
    } else {
        let vec = vec![0i32; nt as usize];
        sort = vec.as_ptr() as *mut i32;
        sort_vec = Some(vec);
    }

    n = gene_sort(d, nt, sort, sw);

    // Verify tmrna_struct magic number
    j = (*sw).tmrna_struct[54] as i32;
    for i in 55..=60 {
        j += (*sw).tmrna_struct[i] as i32;
    }
    if j != (((*sw).tmrna_struct[0] as i32) << 4) + 9 {
        return;
    }

    if (*sw).libflag < 2 {
        if n > 0 {
            j = 0;
            while j < n {
                i = *sort.offset(j as isize);
                j += 1;
                t = TS.offset(i as isize);
                (*t).energy = nenergy(t, sw);

                match (*t).genetype {
                    x if x == tRNA => {
                        init_matrix(&mut m);
                        disp_gene(t, m.as_mut_ptr() as *mut [u8; MATY], sw);
                        let s_str = format!("{}.", j);
                        for (idx, byte) in s_str.bytes().enumerate() {
                            if idx < s.len() - 1 { s[idx] = byte; }
                        }
                        s[s_str.len().min(s.len() - 1)] = 0;
                        xcopy(m.as_mut_ptr() as *mut [u8; MATY], 0, 32, s.as_ptr(), length(s.as_ptr()));
                        disp_matrix(f, m.as_mut_ptr() as *mut [u8; MATY], MATY as i32);

                        if (*sw).matchacceptor != 0 {
                            if iamatch(d, t, sw) == 0 {
                                let _ = writeln!(&mut *f, "    Iso-acceptor mismatch");
                                (*sw).iamismatch += 1;
                            }
                        }
                        if (*sw).annotated != 0 {
                            annotation_overlap_check(d, t, b"    \0".as_ptr(), sw);
                        }
                        overlap(d, sort, n, i, sw);
                        if ((*sw).secstructdisp & 2) != 0 {
                            disp_trna_bracket_notation(f, t, sw);
                        }
                        if ((*sw).secstructdisp & 4) != 0 {
                            disp_gene_SVG(t, m.as_mut_ptr() as *mut [u8; MATY], sw);
                        }
                        if (*sw).seqdisp != 0 {
                            disp_seq(f, t, sw);
                        }
                        if (*t).nintron > 0 {
                            disp_intron(f, t, sw);
                        }
                        if (*sw).energydisp > 1 {
                            trna_score(f, t);
                        }
                    }
                    x if x == tmRNA => {
                        if ((*sw).secstructdisp & 1) != 0 || ((*sw).secstructdisp & 4) != 0 {
                            init_matrix(&mut m);
                            disp_gene(t, m.as_mut_ptr() as *mut [u8; MATY], sw);
                        }
                        if ((*sw).secstructdisp & 1) != 0 {
                            let s_str = format!("{}.", j);
                            for (idx, byte) in s_str.bytes().enumerate() {
                                if idx < s.len() - 1 { s[idx] = byte; }
                            }
                            s[s_str.len().min(s.len() - 1)] = 0;
                            xcopy(m.as_mut_ptr() as *mut [u8; MATY], 0, 32, s.as_ptr(), length(s.as_ptr()));
                            disp_matrix(f, m.as_mut_ptr() as *mut [u8; MATY], MATY as i32);
                            if (*sw).annotated != 0 {
                                annotation_overlap_check(d, t, b"    \0".as_ptr(), sw);
                            }
                        } else {
                            let _ = writeln!(&mut *f, "\n{}.", j);
                            disp_location(t, sw, b"Location\0".as_ptr());
                            if (*sw).reportpseudogenes != 0 {
                                if pseudogene(t, sw) != 0 {
                                    let _ = writeln!(&mut *f, "Possible Pseudogene");
                                }
                            }
                            if (*sw).energydisp != 0 {
                                let _ = writeln!(&mut *f, "Score = {}", (*t).energy);
                            }
                            if (*sw).annotated != 0 {
                                annotation_overlap_check(d, t, b"\0".as_ptr(), sw);
                            }
                        }
                        overlap(d, sort, n, i, sw);
                        if (*t).asst == 0 {
                            disp_tmrna_seq(f, t, sw);
                        } else {
                            disp_tmrna_perm_seq(f, t, sw);
                        }
                        if ((*sw).secstructdisp & 4) != 0 {
                            disp_gene_SVG(t, m.as_mut_ptr() as *mut [u8; MATY], sw);
                        }
                        if (*sw).energydisp > 1 {
                            tmrna_score(f, t, sw);
                        }
                    }
                    x if x == CDS => {
                        let _ = writeln!(&mut *f, "\n{}.\nCDS gene", j);
                        disp_location(t, sw, b"Location\0".as_ptr());
                        if (*sw).annotated != 0 {
                            annotation_overlap_check(d, t, b"\0".as_ptr(), sw);
                        }
                        overlap(d, sort, n, i, sw);
                        disp_cds(f, t, sw);
                    }
                    _ => {}
                }

                if (*sw).libflag > 0 {
                    write_to_library(f, t, sw);
                }
            }
        } else {
            if *(*d).seqname.as_ptr() != 0 {
                let seqname_rust = CStr::from_ptr((*d).seqname.as_ptr() as *const i8).to_string_lossy();
                let _ = writeln!(&mut *f, "\nNothing found in {}\n\n", seqname_rust);
            } else {
                let _ = writeln!(&mut *f, "\nNothing found\n\n");
            }
        }
    } else {
        if n > 0 {
            for i in 0..n {
                write_to_library(f, TS.offset(*sort.offset(i as isize) as isize), sw);
            }
        }
    }

    disp_energy_stats(d, nt, sw);
    if (*d).datatype == GENBANK {
        disp_match(d, sort, n, sw);
    }

    // sort_vec is automatically dropped here
    drop(sort_vec);
}

/// batch_gene_set - Batch process gene set
/// output.c:772-818
pub unsafe fn batch_gene_set(d: *mut data_set, nt: i32, sw: *mut csw) {
    let mut i: i32;
    let mut j: i32;
    let mut n: i32;
    let mut vsort: [i32; NT] = [0; NT];
    let mut nspaces: i32;
    let mut caps: i32;
    let mut sort: *mut i32;
    let mut t: *mut Gene;
    let f = (*sw).f;

    let mut sort_vec: Option<Vec<i32>> = None;
    if nt <= NT as i32 {
        sort = vsort.as_mut_ptr();
    } else {
        let vec = vec![0i32; nt as usize];
        sort = vec.as_ptr() as *mut i32;
        sort_vec = Some(vec);
    }

    n = gene_sort(d, nt, sort, sw);

    // Verify tmrna_struct magic number
    j = (*sw).tmrna_struct[54] as i32;
    for i in 55..=60 {
        j += (*sw).tmrna_struct[i] as i32;
    }
    if j != (((*sw).tmrna_struct[0] as i32) << 4) + 9 {
        return;
    }

    if (*sw).libflag < 2 {
        if (*sw).batch >= 2 {
            nspaces = (*sw).batch & 0x4;
            caps = (*sw).batch & 0x10;
            if ((*sw).batch & 0x8) != 0 {
                for i in 0..n {
                    disp_fasta_seq(f, TS.offset(*sort.offset(i as isize) as isize),
                        (*d).ns + 1, i + 1, nspaces, caps, sw);
                }
            } else {
                for i in 0..n {
                    disp_fasta_seq(f, TS.offset(*sort.offset(i as isize) as isize),
                        0, 0, nspaces, caps, sw);
                }
            }
        } else {
            if n == 1 {
                let _ = writeln!(&mut *f, "1 gene found");
            } else {
                let _ = writeln!(&mut *f, "{} genes found", n);
            }
            for j in 0..n {
                let _ = write!(&mut *f, "{:<3} ", j + 1);
                t = TS.offset(*sort.offset(j as isize) as isize);
                (*t).energy = nenergy(t, sw);
                match (*t).genetype {
                    x if x == tRNA => { disp_batch_trna(f, t, sw); }
                    x if x == tmRNA => { disp_batch_tmrna(f, t, sw); }
                    x if x == srpRNA => { disp_batch_srprna(f, t, sw); }
                    x if x == CDS => { disp_batch_cds(f, t, sw); }
                    _ => {}
                }
            }
        }
    }

    if (*sw).libflag > 0 {
        for i in 0..n {
            write_to_library(f, TS.offset(*sort.offset(i as isize) as isize), sw);
        }
    }

    batch_energy_stats(d, nt, sw);

    // sort_vec is automatically dropped here
    drop(sort_vec);
}

/// remove_overlapping_trna - Remove overlapping tRNA genes
/// output.c:821-894
pub unsafe fn remove_overlapping_trna(d: *mut data_set, nt: i32, sw: *mut csw) {
    let mut i: i32;
    let mut n: i32;
    let mut ioverlay: i32;
    let mut a: i64;
    let mut b: i64;
    let mut c: i64;
    let mut e: i64;
    let mut len: i64;
    let mut leni: i64;
    let mut overlap_val: i64;
    let psmax: i64;
    let mut s1: [u8; 80] = [0; 80];
    let mut s2: [u8; 80] = [0; 80];
    let mut t: *mut Gene;
    let mut ti: *mut Gene;
    let proximity: i64 = (7 * MINCTRNALEN / 10) as i64;

    psmax = (*d).psmax;
    ioverlay = (*sw).ioverlay;

    for n in 0..nt {
        t = TS.offset(n as isize);
        if (*t).genetype != tRNA { continue; }
        if (*t).energy < 0.0 { continue; }
        if (*t).nintron <= 0 { continue; }
        a = (*t).start;
        b = (*t).stop;
        if b < a { b += psmax; }
        len = b - a;

        for i in 0..nt {
            if i == n { continue; }
            ti = TS.offset(i as isize);
            if (*ti).genetype != tRNA { continue; }
            if (*ti).comp != (*t).comp { continue; }
            if (*ti).energy < 0.0 { continue; }
            c = (*ti).start;
            e = (*ti).stop;
            if e < c { e += psmax; }
            leni = e - c;

            if ioverlay != 0 {
                if (2 * len) > (5 * leni) { continue; }
                if (2 * leni) > (5 * len) { continue; }
            }

            overlap_val = if a >= c {
                if b >= e { e - a } else { len }
            } else {
                if b >= e { len } else { b - c }
            };

            if overlap_val >= proximity {
                if (*t).energy < (*ti).energy {
                    if (*sw).verbose != 0 {
                        let nm_rust = CStr::from_ptr(name(t, s1.as_mut_ptr(), 0, sw) as *const i8).to_string_lossy();
                        let pos_rust = CStr::from_ptr(position(s2.as_mut_ptr(), t, sw) as *const i8).to_string_lossy();
                        eprint!("Removing {} at {}", nm_rust, pos_rust);
                        if (*sw).energydisp != 0 {
                            eprint!(" ({})", nenergy(t, sw));
                        }
                        eprintln!();
                    }
                    (*t).energy = -1.0;
                    break;
                }
            }
        }
    }

    for n in 0..(nt - 1) {
        t = TS.offset(n as isize);
        if (*t).genetype != tRNA { continue; }
        if (*t).energy < 0.0 { continue; }
        a = (*t).start;
        b = (*t).stop;
        if b < a { b += psmax; }
        len = b - a;

        for i in (n + 1)..nt {
            ti = TS.offset(i as isize);
            if (*ti).genetype != tRNA { continue; }
            if (*ti).comp != (*t).comp { continue; }
            if (*ti).energy < 0.0 { continue; }
            c = (*ti).start;
            e = (*ti).stop;
            if e < c { e += psmax; }
            leni = e - c;

            if ioverlay != 0 {
                if (2 * len) > (5 * leni) { continue; }
                if (2 * leni) > (5 * len) { continue; }
            }

            overlap_val = if a >= c {
                if b >= e { e - a } else { len }
            } else {
                if b >= e { len } else { b - c }
            };

            if overlap_val >= proximity {
                if (*t).energy < (*ti).energy {
                    if (*sw).verbose != 0 {
                        let nm_rust = CStr::from_ptr(name(t, s1.as_mut_ptr(), 0, sw) as *const i8).to_string_lossy();
                        let pos_rust = CStr::from_ptr(position(s2.as_mut_ptr(), t, sw) as *const i8).to_string_lossy();
                        eprint!("Removing {} at {}", nm_rust, pos_rust);
                        if (*sw).energydisp != 0 {
                            eprint!(" ({})", nenergy(t, sw));
                        }
                        eprintln!();
                    }
                    (*t).energy = -1.0;
                    break;
                } else if (*ti).energy < (*t).energy {
                    if (*sw).verbose != 0 {
                        let nm_rust = CStr::from_ptr(name(ti, s1.as_mut_ptr(), 0, sw) as *const i8).to_string_lossy();
                        let pos_rust = CStr::from_ptr(position(s2.as_mut_ptr(), ti, sw) as *const i8).to_string_lossy();
                        eprint!("Removing {} at {}", nm_rust, pos_rust);
                        if (*sw).energydisp != 0 {
                            eprint!(" ({})", nenergy(ti, sw));
                        }
                        eprintln!();
                    }
                    (*ti).energy = -1.0;
                }
            }
        }
    }
}

/// iopt_fastafile - Interactive mode FASTA file processing
/// output.c:898-1207
pub unsafe fn iopt_fastafile(d: *mut data_set, sw: *mut csw) {
    let mut i: i32;
    let mut nt: i32;
    let mut flag: i32;
    let mut len: i32;
    let mut anticodon: i32;
    let mut skip_nx: i32;
    let mut s: *mut i32;
    let mut sf: *mut i32;
    let mut se: *mut i32;
    let mut sc: *mut i32;
    let mut swrap: *mut i32;
    let mut seq: [i32; 2 * LSEQ + WRAP + 1] = [0; 2 * LSEQ + WRAP + 1];
    let mut cseq: [i32; 2 * LSEQ + WRAP + 1] = [0; 2 * LSEQ + WRAP + 1];
    let mut wseq: [i32; 2 * WRAP + 1] = [0; 2 * WRAP + 1];
    let mut gap: i64;
    let mut start: i64;
    let mut rewind: i64;
    let mut drewind: i64;
    let mut psmax: i64;
    let mut tmaxlen: i64;
    let mut vstart: i64;
    let mut vstop: i64;
    let mut sens: f64;
    let mut sel1: f64;
    let mut sel2: f64;
    let mut c1: i8;
    let mut c2: i8;
    let mut c3: i8;

    // Static arrays from C (output.c lines 905-942)
    static TRNATYPENAME: [[u8; 25]; 3] = [
        *b"Metazoan mitochondrial\0\0\0",
        *b"Cytosolic\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Mammalian mitochondrial\0\0",
    ];
    static GENECODENAME: [[u8; 50]; NGENECODE] = [
        *b"composite Metazoan mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"vertebrate mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"yeast mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"mold/protozoan/Coelenterate mitochondrial\0\0\0\0\0\0\0\0\0",
        *b"invertebrate mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Ciliate\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Echinoderm/Flatworm mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Euplotid\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"bacterial/plant chloroplast\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"alternative yeast\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Ascidian mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"alternative flatworm mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Blepharisma\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Chlorophycean mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Trematode mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Scenedesmus obliquus mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Thraustochytrium mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Pterobranchia mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Gracilibacteria\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Pachysolen tannophilus\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Karyorelict\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Condylostoma\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Mesodinium\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Peritrich\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Blastocrithidia\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"vacant -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Cephalodiscidae mitochondrial UAA-Tyr\0\0\0\0\0\0\0\0\0\0\0\0\0",
    ];

    let f = (*sw).f;
    let aragorn = (*sw).trna != 0 || (*sw).tmrna != 0 || (*sw).cds != 0 || (*sw).srprna != 0;

    // Initialize tmRNA (output.c line 944)
    init_tmrna(f, sw);

    // Print references (output.c lines 946-959)
    let _ = write!(&mut *f, "\nPlease reference the following paper");
    if aragorn && (*sw).mtrna != 0 {
        let _ = write!(&mut *f, "s");
    }
    let _ = writeln!(&mut *f, " if you use this");
    let _ = writeln!(&mut *f, "program as part of any published research.\n");

    if aragorn {
        let _ = writeln!(&mut *f, "Laslett, D. and Canback, B. (2004) ARAGORN, a");
        let _ = writeln!(&mut *f, "program for the detection of transfer RNA and");
        let _ = writeln!(&mut *f, "transfer-messenger RNA genes in nucleotide sequences.");
        let _ = writeln!(&mut *f, "Nucleic Acids Research, 32;11-16.\n");
    }

    if (*sw).mtrna != 0 {
        let _ = writeln!(&mut *f, "Laslett, D. and Canback, B. (2008) ARWEN: a");
        let _ = writeln!(&mut *f, "program to detect tRNA genes in metazoan mitochondrial");
        let _ = writeln!(&mut *f, "nucleotide sequences");
        let _ = writeln!(&mut *f, "Bioinformatics, 24(2); 172-175.\n\n");
    }

    // Print search configuration (output.c lines 960-1006)
    let _ = writeln!(&mut *f);

    if (*sw).mtrna != 0 {
        let trnatypename_rust = CStr::from_ptr(TRNATYPENAME[(*sw).discrim as usize].as_ptr() as *const i8).to_string_lossy();
        let _ = writeln!(&mut *f, "Searching for {} tRNA genes", trnatypename_rust);
        if (*sw).tvloop == 0 {
            let _ = writeln!(&mut *f, "TV replacement loop tRNA genes not detected");
        }
    } else if (*sw).trna != 0 {
        let _ = write!(&mut *f, "Searching for tRNA genes");
        if (*sw).maxintronlen > 0 {
            let _ = write!(&mut *f, " with introns in anticodon loop");
        } else {
            let _ = write!(&mut *f, " with no introns");
        }
        let _ = writeln!(&mut *f);
        if (*sw).maxintronlen > 0 {
            let _ = writeln!(&mut *f, "Intron length from {} to {} bases",
                (*sw).minintronlen, (*sw).maxintronlen);
            if (*sw).ifixedpos != 0 {
                let _ = writeln!(&mut *f, "Intron position fixed between positions 37 and 38");
                let _ = writeln!(&mut *f, "on C-loop (one base after anticodon)");
            }
            if (*sw).ioverlay != 0 {
                let _ = writeln!(&mut *f, "Allowing overlay of long tRNA genes");
            }
        }
    }

    if (*sw).tmrna != 0 {
        let _ = writeln!(&mut *f, "Searching for tmRNA genes");
    }

    if (*sw).linear != 0 {
        let _ = writeln!(&mut *f, "Assuming linear topology, search will not wrap around ends");
    } else {
        let _ = writeln!(&mut *f, "Assuming circular topology, search wraps around ends");
    }

    if (*sw).both == 2 {
        let _ = writeln!(&mut *f, "Searching both strands");
    } else if (*sw).both == 1 {
        let _ = writeln!(&mut *f, "Searching complementary (antisense) strand only");
    } else {
        let _ = writeln!(&mut *f, "Searching single (sense) strand only");
    }

    if (*sw).mtrna != 0 && (*sw).mtcompov != 0 {
        let _ = writeln!(&mut *f, "Reporting overlapping candidates on opposite strands");
    }

    if (*sw).mtrna != 0 || (*sw).trna != 0 || (*sw).tmrna != 0 {
        let genecodename_rust = CStr::from_ptr(GENECODENAME[(*sw).geneticcode as usize].as_ptr() as *const i8).to_string_lossy();
        let _ = writeln!(&mut *f, "Using {} genetic code", genecodename_rust);
        if (*sw).ngcmod > 0 {
            let _ = writeln!(&mut *f, "Specified modifications:");
            for i in 0..(*sw).ngcmod {
                anticodon = (*sw).gcmod[i as usize];
                c1 = cpbase(Thymine - (anticodon & 0x3)) as i8;
                c2 = cpbase(Thymine - ((anticodon >> 2) & 0x3)) as i8;
                c3 = cpbase(Thymine - ((anticodon >> 4) & 0x3)) as i8;
                let aaname_rust = CStr::from_ptr(AANAME[AAMAP[(*sw).geneticcode as usize][anticodon as usize] as usize].as_ptr() as *const i8).to_string_lossy();
                let _ = writeln!(&mut *f, "{}{}{} = {}",
                    c1 as u8 as char, c2 as u8 as char, c3 as u8 as char, aaname_rust);
            }
        }
    }

    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f);

    // Set up rewind distances (output.c lines 1008-1019)
    rewind = (MAXTAGDIST as i64) + 20;
    if ((*sw).trna | (*sw).mtrna) != 0 {
        tmaxlen = (MAXTRNALEN as i64) + (*sw).maxintronlen as i64;
        if rewind < tmaxlen {
            rewind = tmaxlen;
        }
    }
    if (*sw).tmrna != 0 {
        if rewind < MAXTMRNALEN as i64 {
            rewind = MAXTMRNALEN as i64;
        }
    }
    if (*sw).peptide != 0 {
        if (*sw).tagthresh >= 5 {
            if rewind < TSWEEP as i64 {
                rewind = TSWEEP as i64;
            }
        }
    }
    (*sw).loffset = rewind as i32;
    (*sw).roffset = rewind as i32;
    drewind = 2 * rewind;

    // Initialize sequence counters (output.c lines 1020-1023)
    (*d).ns = 0;
    (*d).nf = 0;
    (*d).nextseq = 0;
    (*d).nextseqoff = 0;

    // Main sequence processing loop (output.c lines 1024-1132)
    while (*d).nextseq >= 0 {
        (*d).seqstart = (*d).nextseq;
        (*d).seqstartoff = (*d).nextseqoff;

        if seq_init(d, sw) == 0 {
            break;
        }

        psmax = (*d).psmax;

        if (*sw).verbose != 0 {
            let seqname_rust = CStr::from_ptr((*d).seqname.as_ptr() as *const i8).to_string_lossy();
            eprintln!("{}", seqname_rust);
            eprintln!("{} nucleotides in sequence", psmax);
            eprintln!("Mean G+C content = {:.1}%", 100.0 * (*d).gc);
        }

        let seqname_rust = CStr::from_ptr((*d).seqname.as_ptr() as *const i8).to_string_lossy();
        let _ = writeln!(&mut *f, "{}", seqname_rust);
        let _ = writeln!(&mut *f, "{} nucleotides in sequence", psmax);
        let _ = writeln!(&mut *f, "Mean G+C content = {:.1}%", 100.0 * (*d).gc);

        init_gene(0, NT as i32);
        nt = 0;
        flag = 0;
        start = 1;
        se = seq.as_mut_ptr();
        skip_nx = 0;

        if (*sw).linear != 0 {
            for i in 0..rewind as i32 {
                *se = NOBASE;
                se = se.offset(1);
            }
            start -= rewind;
        } else {
            if psmax <= drewind {
                gap = drewind - psmax;
                sc = se.offset(gap as isize);
                while se < sc {
                    *se = NOBASE;
                    se = se.offset(1);
                }
                swrap = wseq.as_mut_ptr();
                sc = se.offset(psmax as isize);
                while se < sc {
                    *se = move_forward(d);
                    *swrap = *se;
                    swrap = swrap.offset(1);
                    se = se.offset(1);
                }
                sc = swrap.offset(gap as isize);
                while swrap < sc {
                    *swrap = NOBASE;
                    swrap = swrap.offset(1);
                }
                swrap = wseq.as_mut_ptr();
                sc = swrap.offset(psmax as isize);
                while swrap < sc {
                    *se = *swrap;
                    swrap = swrap.offset(1);
                    se = se.offset(1);
                }
                swrap = wseq.as_mut_ptr();
                sc = swrap.offset(drewind as isize);
                while swrap < sc {
                    *se = *swrap;
                    swrap = swrap.offset(1);
                    se = se.offset(1);
                }
                (*sw).loffset = drewind as i32;
                (*sw).roffset = drewind as i32;
                start -= drewind;
                flag = 1;
                skip_nx = 1;
            } else {
                swrap = wseq.as_mut_ptr();
                sc = seq.as_mut_ptr().offset(drewind as isize);
                while se < sc {
                    *se = move_forward(d);
                    *swrap = *se;
                    swrap = swrap.offset(1);
                    se = se.offset(1);
                }
                skip_nx = 0;
            }
        }

        sc = seq.as_mut_ptr().offset(LSEQ as isize);

        loop {
            if skip_nx == 0 {
                while se < sc {
                    if (*d).ps >= psmax {
                        if (*sw).linear != 0 {
                            for i in 0..rewind as i32 {
                                *se = NOBASE;
                                se = se.offset(1);
                            }
                        } else {
                            sc = wseq.as_mut_ptr().offset(drewind as isize);
                            swrap = wseq.as_mut_ptr();
                            while swrap < sc {
                                *se = *swrap;
                                swrap = swrap.offset(1);
                                se = se.offset(1);
                            }
                        }
                        flag = 1;
                        break;
                    } else {
                        *se = move_forward(d);
                        se = se.offset(1);
                    }
                }
            }
            skip_nx = 0;

            len = (se as isize - seq.as_ptr() as isize) as i32 / std::mem::size_of::<i32>() as i32;

            if (*sw).verbose != 0 {
                vstart = sq(start + (*sw).loffset as i64, psmax);
                vstop = sq(start + len as i64 - (*sw).roffset as i64 - 1, psmax);
                if vstop < vstart {
                    eprintln!("Searching from {} to {}", vstart, psmax);
                    eprintln!("Searching from 1 to {}", vstop);
                } else {
                    eprintln!("Searching from {} to {}", vstart, vstop);
                }
            }

            if (*sw).both != 1 {
                (*sw).start = start;
                (*sw).comp = 0;
                nt = tmioptimise(d, seq.as_mut_ptr(), len, nt, sw);
            }

            if (*sw).both > 0 {
                sense_switch(seq.as_mut_ptr(), cseq.as_mut_ptr(), len);
                (*sw).start = start + len as i64;
                (*sw).comp = 1;
                nt = tmioptimise(d, cseq.as_mut_ptr(), len, nt, sw);
            }

            if flag == 0 {
                s = seq.as_mut_ptr();
                sf = se.offset(-(drewind as isize));
                se = seq.as_mut_ptr().offset(drewind as isize);
                while s < se {
                    *s = *sf;
                    s = s.offset(1);
                    sf = sf.offset(1);
                }
                start += len as i64 - drewind;
                continue;
            }
            break;
        }

        if nt < 1 {
            (*d).nf += 1;
        }

        if (*sw).maxintronlen > 0 {
            remove_overlapping_trna(d, nt, sw);
        }

        if (*sw).updatetmrnatags != 0 {
            update_tmrna_tag_database(TS, nt, sw);
        }

        disp_gene_set(d, nt, sw);

        if (*sw).verbose != 0 {
            let seqname_rust = CStr::from_ptr((*d).seqname.as_ptr() as *const i8).to_string_lossy();
            eprintln!("{}\nSearch Finished\n", seqname_rust);
        }

        (*d).ns += 1;
    }

    // Print summary statistics (output.c lines 1133-1206)
    if (*d).ns > 1 {
        let _ = writeln!(&mut *f, "\n\n{} sequences searched", (*d).ns);
        if ((*sw).trna | (*sw).mtrna) != 0 {
            let _ = writeln!(&mut *f, "Total tRNA genes = {}", (*sw).ngene[tRNA as usize]);
            if (*sw).matchacceptor != 0 {
                let _ = writeln!(&mut *f, "Total iso-acceptor mismatches = {}", (*sw).iamismatch);
            }
        }
        if (*sw).tmrna != 0 {
            let _ = writeln!(&mut *f, "Total tmRNA genes = {}", (*sw).ngene[tmRNA as usize]);
        }
        if (*sw).reportpseudogenes != 0 && (*sw).nps > 0 {
            let _ = writeln!(&mut *f, "Total number of possible pseudogenes = {}", (*sw).nps);
        }
        if (*d).nf > 0 {
            sens = 100.0 * ((*d).ns - (*d).nf) as f64 / (*d).ns as f64;
            let _ = writeln!(&mut *f, "Nothing found in {} sequences ({:.2}% sensitivity)",
                (*d).nf, sens);
        }
    }

    if (*sw).updatetmrnatags != 0 {
        report_new_tmrna_tags(sw);
    }
}

/// bopt_fastafile - Batch mode FASTA file processing
/// output.c:1210-1338
pub unsafe fn bopt_fastafile(d: *mut data_set, sw: *mut csw) {
    let mut i: i32;
    let mut nt: i32;
    let mut flag: i32;
    let mut len: i32;
    let mut skip_nx: i32 = 0;

    let mut seq: [i32; 2 * LSEQ + WRAP + 1] = [0; 2 * LSEQ + WRAP + 1];
    let mut cseq: [i32; 2 * LSEQ + WRAP + 1] = [0; 2 * LSEQ + WRAP + 1];
    let mut wseq: [i32; 2 * WRAP + 1] = [0; 2 * WRAP + 1];

    let mut gap: i64;
    let mut start: i64;
    let mut rewind: i64;
    let mut drewind: i64;
    let mut psmax: i64;
    let mut tmaxlen: i64;
    let mut vstart: i64;
    let mut vstop: i64;
    let mut sens: f64;

    let f = (*sw).f;

    // sq macro: ((pos + d->psmax - 1L) % d->psmax) + 1L
    // Will be inlined where needed using psmax variable

    rewind = (MAXTAGDIST + 20) as i64;
    if (*sw).trna != 0 || (*sw).mtrna != 0 {
        tmaxlen = (MAXTRNALEN as i32 + (*sw).maxintronlen) as i64;
        if rewind < tmaxlen {
            rewind = tmaxlen;
        }
    }
    if (*sw).tmrna != 0 {
        if rewind < MAXTMRNALEN as i64 {
            rewind = MAXTMRNALEN as i64;
        }
    }
    if (*sw).peptide != 0 {
        if (*sw).tagthresh >= 5 {
            if rewind < TSWEEP as i64 {
                rewind = TSWEEP as i64;
            }
        }
    }

    (*sw).loffset = rewind as i32;
    (*sw).roffset = rewind as i32;
    drewind = 2 * rewind;
    (*d).ns = 0;
    (*d).nf = 0;
    (*d).nextseq = 0;
    (*d).nextseqoff = 0;

    while (*d).nextseq >= 0 {
        (*d).seqstart = (*d).nextseq;
        (*d).seqstartoff = (*d).nextseqoff;

        if seq_init(d, sw) == 0 {
            break;
        }

        psmax = (*d).psmax;

        if (*sw).verbose != 0 {
            let seqname_rust = CStr::from_ptr((*d).seqname.as_ptr() as *const i8).to_string_lossy();
            eprintln!("{}", seqname_rust);
            eprintln!("{} nucleotides in sequence", psmax);
            eprintln!("Mean G+C content = {:.1}%", 100.0 * (*d).gc);
        }

        if (*sw).batch < 2 {
            let seqname_rust = CStr::from_ptr((*d).seqname.as_ptr() as *const i8).to_string_lossy();
            let _ = writeln!(&mut *f, ">{}", seqname_rust);
        }

        init_gene(0, NT as i32);
        nt = 0;
        flag = 0;
        start = 1;

        let mut se_idx: usize = 0;

        if (*sw).linear != 0 {
            // for (i = 0; i < rewind; i++) *se++ = NOBASE;
            i = 0;
            while i < rewind as i32 {
                seq[se_idx] = NOBASE;
                se_idx += 1;
                i += 1;
            }
            start -= rewind;
        } else {
            if psmax <= drewind {
                gap = drewind - psmax;
                // Fill gap with NOBASE
                let mut j: i64 = 0;
                while j < gap {
                    seq[se_idx] = NOBASE;
                    se_idx += 1;
                    j += 1;
                }

                // Read sequence and copy to wseq
                let mut swrap_idx: usize = 0;
                j = 0;
                while j < psmax {
                    seq[se_idx] = move_forward(d);
                    wseq[swrap_idx] = seq[se_idx];
                    se_idx += 1;
                    swrap_idx += 1;
                    j += 1;
                }

                // Fill rest of wseq with NOBASE
                j = 0;
                while j < gap {
                    wseq[swrap_idx] = NOBASE;
                    swrap_idx += 1;
                    j += 1;
                }

                // Copy psmax from wseq to seq
                swrap_idx = 0;
                j = 0;
                while j < psmax {
                    seq[se_idx] = wseq[swrap_idx];
                    se_idx += 1;
                    swrap_idx += 1;
                    j += 1;
                }

                // Copy drewind from wseq to seq
                swrap_idx = 0;
                j = 0;
                while j < drewind {
                    seq[se_idx] = wseq[swrap_idx];
                    se_idx += 1;
                    swrap_idx += 1;
                    j += 1;
                }

                (*sw).loffset = drewind as i32;
                (*sw).roffset = drewind as i32;
                start -= drewind;
                flag = 1;
                skip_nx = 1;
            } else {
                // Read drewind bases into both seq and wseq
                let mut swrap_idx: usize = 0;
                let mut j: i64 = 0;
                while j < drewind {
                    seq[se_idx] = move_forward(d);
                    wseq[swrap_idx] = seq[se_idx];
                    se_idx += 1;
                    swrap_idx += 1;
                    j += 1;
                }
                skip_nx = 0;
            }
        }

        let sc_limit: usize = LSEQ;

        loop {
            if skip_nx == 0 {
                while se_idx < sc_limit {
                    seq[se_idx] = move_forward(d);
                    se_idx += 1;
                    if (*d).ps >= psmax {
                        if (*sw).linear != 0 {
                            i = 0;
                            while i < rewind as i32 {
                                seq[se_idx] = NOBASE;
                                se_idx += 1;
                                i += 1;
                            }
                        } else {
                            // Copy drewind from wseq to seq
                            let mut swrap_idx: usize = 0;
                            let mut j: i64 = 0;
                            while j < drewind {
                                seq[se_idx] = wseq[swrap_idx];
                                se_idx += 1;
                                swrap_idx += 1;
                                j += 1;
                            }
                        }
                        flag = 1;
                        break;
                    }
                }
            }
            skip_nx = 0;

            len = se_idx as i32;

            if (*sw).verbose != 0 {
                // sq(pos) = ((pos + psmax - 1) % psmax) + 1
                let pos1 = start + (*sw).loffset as i64;
                vstart = ((pos1 + psmax - 1) % psmax) + 1;
                let pos2 = start + len as i64 - (*sw).roffset as i64 - 1;
                vstop = ((pos2 + psmax - 1) % psmax) + 1;
                if vstop < vstart {
                    eprintln!("Searching from {} to {}", vstart, psmax);
                    eprintln!("Searching from 1 to {}", vstop);
                } else {
                    eprintln!("Searching from {} to {}", vstart, vstop);
                }
            }

            if (*sw).both != 1 {
                (*sw).start = start;
                (*sw).comp = 0;
                nt = tmioptimise(d, seq.as_mut_ptr(), len, nt, sw);
            }

            if (*sw).both > 0 {
                sense_switch(seq.as_mut_ptr(), cseq.as_mut_ptr(), len);
                (*sw).start = start + len as i64;
                (*sw).comp = 1;
                nt = tmioptimise(d, cseq.as_mut_ptr(), len, nt, sw);
            }

            if flag == 0 {
                // Shift sequence: copy from (se - drewind) to beginning
                let sf_start = se_idx - drewind as usize;
                let mut s_idx: usize = 0;
                let mut sf_idx = sf_start;
                while s_idx < drewind as usize {
                    seq[s_idx] = seq[sf_idx];
                    s_idx += 1;
                    sf_idx += 1;
                }
                se_idx = drewind as usize;
                start += len as i64 - drewind;
                continue;
            }
            break;
        }

        if nt < 1 {
            (*d).nf += 1;
        }

        if (*sw).maxintronlen > 0 {
            remove_overlapping_trna(d, nt, sw);
        }

        if (*sw).updatetmrnatags != 0 {
            update_tmrna_tag_database(TS, nt, sw);
        }

        batch_gene_set(d, nt, sw);

        if (*sw).verbose != 0 {
            let seqname_rust = CStr::from_ptr((*d).seqname.as_ptr() as *const i8).to_string_lossy();
            eprintln!("{}\nSearch Finished\n", seqname_rust);
        }

        (*d).ns += 1;
    }

    if (*d).ns > 1 && (*sw).batch < 2 {
        let _ = write!(&mut *f, ">end \t{} sequences", (*d).ns);
        if (*sw).trna != 0 || (*sw).mtrna != 0 {
            let _ = write!(&mut *f, " {} tRNA genes", (*sw).ngene[tRNA as usize]);
        }
        if (*sw).tmrna != 0 {
            let _ = write!(&mut *f, " {} tmRNA genes", (*sw).ngene[tmRNA as usize]);
        }
        if (*d).nf > 0 {
            sens = 100.0 * ((*d).ns - (*d).nf) as f64 / (*d).ns as f64;
            let _ = write!(&mut *f, ", nothing found in {} sequences, ({:.2}% sensitivity)",
                (*d).nf, sens);
        }
        let _ = writeln!(&mut *f);
    }

    if (*sw).updatetmrnatags != 0 {
        report_new_tmrna_tags(sw);
    }
}

// ============================================================================
// Parallel Batch Processing (Case C: Metagenome with many small sequences)
// ============================================================================

/// Convert raw nucleotide bytes to internal representation
#[inline]
fn convert_base(b: u8) -> i32 {
    match b {
        b'A' | b'a' => Adenine,
        b'C' | b'c' => Cytosine,
        b'G' | b'g' => Guanine,
        b'T' | b't' | b'U' | b'u' => Thymine,
        // Ambiguous bases (IUPAC codes) - return AMBIG to match move_forward behavior
        b'N' | b'n' | b'R' | b'r' | b'Y' | b'y' | b'K' | b'k' |
        b'M' | b'm' | b'S' | b's' | b'W' | b'w' | b'H' | b'h' |
        b'B' | b'b' | b'V' | b'v' | b'D' | b'd' => AMBIG,
        _ => NOBASE,
    }
}

/// Calculate GC content from nucleotide sequence
fn calculate_gc(seq: &[u8]) -> f64 {
    let mut gc_count = 0i64;
    let mut total = 0i64;
    for &base in seq {
        match base {
            b'G' | b'g' | b'C' | b'c' => { gc_count += 1; total += 1; }
            b'A' | b'a' | b'T' | b't' | b'U' | b'u' => { total += 1; }
            _ => {}
        }
    }
    if total > 0 { gc_count as f64 / total as f64 } else { 0.0 }
}

/// Thread-safe gene data for parallel processing
/// Contains all Gene fields except raw pointers (which aren't Send)
#[derive(Clone)]
struct SendableGene {
    pub name: [u8; 100],
    pub seq: [i32; MAXTRNALEN + 1],
    pub eseq: [i32; MAXETRNALEN + 1],
    // ps: *mut i32 is NOT included - not thread-safe
    pub nbase: i32,
    pub comp: i32,
    pub start: i64,
    pub stop: i64,
    pub astem1: i32,
    pub astem2: i32,
    pub aatail: i32,
    pub spacer1: i32,
    pub spacer2: i32,
    pub dstem: i32,
    pub dloop: i32,
    pub cstem: i32,
    pub cloop: i32,
    pub intron: i32,
    pub nintron: i32,
    pub anticodon: i32,
    pub var: i32,
    pub varbp: i32,
    pub tstem: i32,
    pub tloop: i32,
    pub genetype: i32,
    pub energy: f64,
    pub asst: i32,
    pub tps: i32,
    pub tpe: i32,
    pub annotation: i32,
    pub annosc: i32,
}

impl SendableGene {
    /// Convert from Gene (drops the raw pointer)
    unsafe fn from_gene(g: &Gene) -> Self {
        Self {
            name: g.name,
            seq: g.seq,
            eseq: g.eseq,
            nbase: g.nbase,
            comp: g.comp,
            start: g.start,
            stop: g.stop,
            astem1: g.astem1,
            astem2: g.astem2,
            aatail: g.aatail,
            spacer1: g.spacer1,
            spacer2: g.spacer2,
            dstem: g.dstem,
            dloop: g.dloop,
            cstem: g.cstem,
            cloop: g.cloop,
            intron: g.intron,
            nintron: g.nintron,
            anticodon: g.anticodon,
            var: g.var,
            varbp: g.varbp,
            tstem: g.tstem,
            tloop: g.tloop,
            genetype: g.genetype,
            energy: g.energy,
            asst: g.asst,
            tps: g.tps,
            tpe: g.tpe,
            annotation: g.annotation,
            annosc: g.annosc,
        }
    }

    /// Convert back to Gene (sets ps to null)
    fn to_gene(&self) -> Gene {
        Gene {
            name: self.name,
            seq: self.seq,
            eseq: self.eseq,
            ps: std::ptr::null_mut(),
            nbase: self.nbase,
            comp: self.comp,
            start: self.start,
            stop: self.stop,
            astem1: self.astem1,
            astem2: self.astem2,
            aatail: self.aatail,
            spacer1: self.spacer1,
            spacer2: self.spacer2,
            dstem: self.dstem,
            dloop: self.dloop,
            cstem: self.cstem,
            cloop: self.cloop,
            intron: self.intron,
            nintron: self.nintron,
            anticodon: self.anticodon,
            var: self.var,
            varbp: self.varbp,
            tstem: self.tstem,
            tloop: self.tloop,
            genetype: self.genetype,
            energy: self.energy,
            asst: self.asst,
            tps: self.tps,
            tpe: self.tpe,
            annotation: self.annotation,
            annosc: self.annosc,
        }
    }
}

/// Result of processing a single sequence in parallel
/// Contains detected genes and metadata needed for output
#[derive(Clone)]
struct SeqProcessResult {
    /// Sequence index (for ordering output)
    idx: usize,
    /// Sequence ID
    id: String,
    /// Sequence length
    psmax: i64,
    /// GC content
    gc: f64,
    /// Detected genes (thread-safe version without raw pointers)
    genes: Vec<SendableGene>,
    /// Number of detected genes
    nt: i32,
    /// Whether any genes were found (for nf counter)
    found_genes: bool,
}

// ============================================================================
// Chunk-based Parallel Processing for Large Sequences
// ============================================================================

/// Minimum sequence length for chunk-based parallelization (2 Mbp)
const MIN_CHUNK_PARALLEL_LEN: usize = 2_000_000;

/// Overlap between chunks to ensure genes at boundaries are detected
/// Must be at least TSWEEP (1000) from C original to handle tRNA sweep distance
/// (MAXTMRNALEN=786 for tmRNA, TSWEEP=1000 for tRNA)
const CHUNK_OVERLAP: usize = 1000;

/// Represents a chunk of sequence data for parallel processing
struct SeqChunk<'a> {
    /// The chunk data (reference to original sequence)
    data: &'a [u8],
    /// Global offset of this chunk in the original sequence
    global_offset: usize,
    /// Chunk index (for ordering)
    chunk_idx: usize,
    /// Whether this is the last chunk
    is_last: bool,
}

/// Represents a chunk with sequence context for unified parallel processing
struct UnifiedChunk<'a> {
    /// Sequence index (which sequence this chunk belongs to)
    seq_idx: usize,
    /// The chunk data
    data: &'a [u8],
    /// Global offset in the original sequence
    global_offset: usize,
    /// Total length of the original sequence
    total_seq_len: usize,
    /// GC content of the sequence
    gc: f64,
    /// Sequence name for gene annotation
    seq_name: &'a str,
    /// Effective loffset for this sequence
    effective_loffset: i32,
    /// Effective roffset for this sequence
    effective_roffset: i32,
}

/// Create overlapping chunks from a sequence for parallel processing
fn create_seq_chunks(seq: &[u8], num_chunks: usize) -> Vec<SeqChunk<'_>> {
    if seq.len() < MIN_CHUNK_PARALLEL_LEN || num_chunks <= 1 {
        // Single chunk for small sequences
        return vec![SeqChunk {
            data: seq,
            global_offset: 0,
            chunk_idx: 0,
            is_last: true,
        }];
    }

    let seq_len = seq.len();
    // Calculate chunk size to distribute work evenly
    let base_chunk_size = seq_len / num_chunks;
    let chunk_size = std::cmp::max(base_chunk_size, MIN_CHUNK_PARALLEL_LEN / 4);

    let mut chunks = Vec::new();
    let mut offset = 0;
    let mut idx = 0;

    while offset < seq_len {
        // Calculate end position with overlap
        let base_end = std::cmp::min(offset + chunk_size, seq_len);
        let end_with_overlap = std::cmp::min(base_end + CHUNK_OVERLAP, seq_len);
        let is_last = end_with_overlap >= seq_len;

        chunks.push(SeqChunk {
            data: &seq[offset..end_with_overlap],
            global_offset: offset,
            chunk_idx: idx,
            is_last,
        });

        offset = base_end;  // Move by chunk_size (not including overlap)
        idx += 1;

        if is_last {
            break;
        }
    }

    chunks
}

/// Create unified chunks from all sequences for fully parallel processing
/// Returns (chunks, needs_dedup_flags) where needs_dedup_flags[seq_idx] indicates
/// whether that sequence was split into multiple chunks (needs deduplication)
fn create_unified_chunks<'a>(
    sequences: &'a [SeqRecord],
    gc_values: &[f64],
    effective_offsets: &[(i32, i32)],  // (loffset, roffset) per sequence
    num_threads: usize,
) -> (Vec<UnifiedChunk<'a>>, Vec<bool>) {
    // Calculate total sequence length to determine optimal chunk size
    let total_len: usize = sequences.iter().map(|s| s.len()).sum();

    // Target: create roughly num_threads * 6 chunks for good load balancing
    // This provides enough granularity for work-stealing while limiting overhead
    let target_chunks = num_threads * 6;
    let optimal_chunk_size = std::cmp::max(total_len / target_chunks, MIN_CHUNK_PARALLEL_LEN);

    let mut all_chunks = Vec::with_capacity(target_chunks + sequences.len());
    let mut needs_dedup = vec![false; sequences.len()];

    for (seq_idx, seq_record) in sequences.iter().enumerate() {
        let seq_len = seq_record.len();
        let gc = gc_values[seq_idx];
        let (effective_loffset, effective_roffset) = effective_offsets[seq_idx];

        // For sequences smaller than optimal chunk size, create a single chunk
        if seq_len <= optimal_chunk_size {
            all_chunks.push(UnifiedChunk {
                seq_idx,
                data: &seq_record.seq,
                global_offset: 0,
                total_seq_len: seq_len,
                gc,
                seq_name: &seq_record.id,
                effective_loffset,
                effective_roffset,
            });
            // needs_dedup[seq_idx] remains false - no deduplication needed
            continue;
        }

        // For large sequences, create multiple chunks
        needs_dedup[seq_idx] = true;  // Will need deduplication

        let mut offset = 0;
        while offset < seq_len {
            let base_end = std::cmp::min(offset + optimal_chunk_size, seq_len);
            let end_with_overlap = std::cmp::min(base_end + CHUNK_OVERLAP, seq_len);

            all_chunks.push(UnifiedChunk {
                seq_idx,
                data: &seq_record.seq[offset..end_with_overlap],
                global_offset: offset,
                total_seq_len: seq_len,
                gc,
                seq_name: &seq_record.id,
                effective_loffset,
                effective_roffset,
            });

            offset = base_end;
            if end_with_overlap >= seq_len {
                break;
            }
        }
    }

    (all_chunks, needs_dedup)
}

/// Process a single chunk of sequence (for chunk-based parallelism)
/// Returns genes with positions adjusted to global coordinates
unsafe fn process_chunk_parallel(
    chunk_data: &[u8],
    global_offset: usize,
    total_seq_len: usize,
    gc: f64,
    rewind: i64,
    drewind: i64,
    sw_linear: i32,
    sw_both: i32,
    sw_ptr: *const csw,
    seq_name: &str,  // Sequence name for gene annotation
    sw_loffset: i32,  // Effective loffset for small circular sequence handling
    sw_roffset: i32,  // Effective roffset for small circular sequence handling
) -> Vec<SendableGene> {
    // Estimate gene capacity based on realistic gene density
    // tRNA genes occur roughly every 30,000-50,000 bases
    // Use conservative estimate: 1 gene per 10,000 bases with minimum of 100
    let chunk_len = chunk_data.len();
    let estimated_genes = std::cmp::max(100, chunk_len / 10000 + 50);
    let ts_capacity = get_thread_local_ts_capacity() as usize;
    let required_capacity = std::cmp::min(estimated_genes, ts_capacity);

    // Get thread-local gene storage with lazy allocation
    let ts = get_thread_local_ts_with_capacity(required_capacity);

    // Initialize only the required slots
    init_gene_ts(ts, 0, required_capacity as i32);

    let mut nt = 0i32;
    // chunk_len already defined above for capacity estimation
    let psmax = total_seq_len as i64;

    // Allocate buffers
    let buf_size = 2 * LSEQ + WRAP + 1;
    let mut seq_buf: Vec<i32> = vec![NOBASE; buf_size];
    let mut cseq_buf: Vec<i32> = vec![NOBASE; buf_size];

    // Create a minimal data_set for this thread
    let mut local_d = DataSet::default();
    local_d.psmax = psmax;
    local_d.gc = gc;
    // CRITICAL: Set seqname so that gene names are properly set during detection
    let name_bytes = seq_name.as_bytes();
    let copy_len = std::cmp::min(name_bytes.len(), STRLEN - 1);
    for (i, &b) in name_bytes[..copy_len].iter().enumerate() {
        local_d.seqname[i] = b;
    }
    local_d.seqname[copy_len] = 0;
    let d_ptr = &mut local_d as *mut DataSet as *mut data_set;

    // For chunks, we process linearly (no wrap-around handling needed)
    // Start position is adjusted by global_offset
    let start_pos: i64 = (global_offset + 1) as i64;  // 1-based

    let mut se_idx: usize = 0;
    let mut read_pos: usize = 0;
    let mut start: i64 = start_pos;
    let mut flag: i32 = 0;

    // Add rewind padding at the beginning (for boundary genes)
    if global_offset > 0 {
        // Not the first chunk - use rewind for overlap handling
        for _ in 0..rewind as usize {
            seq_buf[se_idx] = NOBASE;
            se_idx += 1;
        }
        start -= rewind;
    }

    let sc_limit: usize = LSEQ;

    // Main processing loop
    loop {
        while se_idx < sc_limit {
            if read_pos < chunk_len {
                seq_buf[se_idx] = convert_base(chunk_data[read_pos]);
                read_pos += 1;
                se_idx += 1;
            } else {
                // End of chunk
                for _ in 0..rewind as usize {
                    if se_idx < buf_size {
                        seq_buf[se_idx] = NOBASE;
                        se_idx += 1;
                    }
                }
                flag = 1;
                break;
            }
        }

        let len = se_idx as i32;

        // Create a thread-local copy of csw to avoid race conditions
        let mut local_sw_copy: Csw = std::ptr::read(sw_ptr as *const Csw);
        local_sw_copy.f = std::ptr::null_mut();
        // IMPORTANT: Update genespace to match actually allocated capacity
        local_sw_copy.genespace = required_capacity as i32;
        // CRITICAL: Apply effective loffset/roffset for small circular sequence handling
        // This matches process_sequence_parallel behavior
        local_sw_copy.loffset = sw_loffset;
        local_sw_copy.roffset = sw_roffset;
        let local_sw = &mut local_sw_copy as *mut Csw as *mut csw;

        // Process forward strand
        if sw_both != 1 {
            (*local_sw).start = start;
            (*local_sw).comp = 0;
            nt = tmioptimise_ts(ts, required_capacity as i32, d_ptr, seq_buf.as_mut_ptr(), len, nt, local_sw);
        }

        // Process reverse strand
        if sw_both > 0 {
            for i in 0..len as usize {
                cseq_buf[i] = match seq_buf[i] {
                    0 => 3, 1 => 2, 2 => 1, 3 => 0, _ => NOBASE,
                };
            }
            sense_switch(seq_buf.as_mut_ptr(), cseq_buf.as_mut_ptr(), len);

            (*local_sw).start = start + len as i64;
            (*local_sw).comp = 1;
            nt = tmioptimise_ts(ts, required_capacity as i32, d_ptr, cseq_buf.as_mut_ptr(), len, nt, local_sw);
        }

        if flag == 0 {
            let sf_start = se_idx - drewind as usize;
            for i in 0..drewind as usize {
                seq_buf[i] = seq_buf[sf_start + i];
            }
            se_idx = drewind as usize;
            start += len as i64 - drewind;
            continue;
        }
        break;
    }

    // Copy detected genes to result vector
    let mut genes: Vec<SendableGene> = Vec::with_capacity(nt as usize);
    for i in 0..nt {
        let gene_ptr = ts.offset(i as isize);
        if (*gene_ptr).genetype != NO_GENE {
            genes.push(SendableGene::from_gene(&*gene_ptr));
        }
    }

    genes
}

/// Deduplicate genes from multiple chunks
/// Genes at chunk boundaries may be detected in multiple chunks
fn deduplicate_chunk_genes(mut all_genes: Vec<SendableGene>) -> Vec<SendableGene> {
    if all_genes.is_empty() {
        return all_genes;
    }

    // Sort by start position, then stop, then genetype, then strand
    all_genes.sort_by(|a, b| {
        a.start.cmp(&b.start)
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

/// Process a large sequence using chunk-based parallelism
/// Splits sequence into overlapping chunks, processes in parallel, then merges results
/// Process a sequence using chunk-based parallelization
/// - nested=false: Create dedicated thread pool (for single large sequence)
/// - nested=true: Process chunks sequentially (already inside parallel context)
unsafe fn process_sequence_chunked_parallel(
    seq_record: &SeqRecord,
    gc: f64,
    rewind: i64,
    drewind: i64,
    sw_linear: i32,
    sw_both: i32,
    sw_ptr: *const csw,
    num_threads: usize,
    nested: bool,
    sw_loffset: i32,  // Effective loffset for small circular sequence handling
    sw_roffset: i32,  // Effective roffset for small circular sequence handling
) -> (Vec<SendableGene>, i32) {
    let seq_len = seq_record.len();

    // Create chunks
    let chunks = create_seq_chunks(&seq_record.seq, num_threads);
    let num_chunks = chunks.len();

    if num_chunks <= 1 {
        // Single chunk - no parallel benefit, use direct processing
        let genes = process_chunk_parallel(
            &seq_record.seq,
            0,
            seq_len,
            gc,
            rewind,
            drewind,
            sw_linear,
            sw_both,
            sw_ptr,
            &seq_record.id,
            sw_loffset,
            sw_roffset,
        );
        let nt = genes.len() as i32;
        return (genes, nt);
    }

    let seq_name = &seq_record.id;

    let chunk_results: Vec<Vec<SendableGene>> = if nested {
        // Already inside parallel context - process chunks sequentially
        // This avoids nested thread pool issues
        chunks.iter()
            .map(|chunk| {
                process_chunk_parallel(
                    chunk.data,
                    chunk.global_offset,
                    seq_len,
                    gc,
                    rewind,
                    drewind,
                    sw_linear,
                    sw_both,
                    sw_ptr,
                    seq_name,
                    sw_loffset,
                    sw_roffset,
                )
            })
            .collect()
    } else {
        // Not nested - use global thread pool for parallel chunk processing
        let sw_ptr_usize = sw_ptr as usize;

        chunks.par_iter()
            .map(|chunk| {
                let sw_ptr = sw_ptr_usize as *const csw;
                process_chunk_parallel(
                    chunk.data,
                    chunk.global_offset,
                    seq_len,
                    gc,
                    rewind,
                    drewind,
                    sw_linear,
                    sw_both,
                    sw_ptr,
                    seq_name,
                    sw_loffset,
                    sw_roffset,
                )
            })
            .collect()
    };

    // Flatten and deduplicate results
    let all_genes: Vec<SendableGene> = chunk_results.into_iter().flatten().collect();
    let deduped_genes = deduplicate_chunk_genes(all_genes);
    let nt = deduped_genes.len() as i32;

    (deduped_genes, nt)
}

/// Process a single sequence using thread-local state
/// This function can be called from multiple threads safely
unsafe fn process_sequence_parallel(
    seq_record: &SeqRecord,
    gc: f64,
    rewind: i64,
    drewind: i64,
    sw_linear: i32,
    sw_both: i32,
    sw_trna: i32,
    sw_mtrna: i32,
    sw_tmrna: i32,
    sw_loffset: i32,
    sw_roffset: i32,
    sw_maxintronlen: i32,
    sw_minintronlen: i32,
    sw_verbose: i32,
    // Read-only csw fields that tmioptimise_ts needs
    sw_ptr: *const csw,
) -> (Vec<SendableGene>, i32) {
    // Estimate gene capacity based on realistic gene density
    // tRNA genes occur roughly every 30,000-50,000 bases
    // Use conservative estimate: 1 gene per 10,000 bases with minimum of 100
    let seq_len = seq_record.seq.len();
    let estimated_genes = std::cmp::max(100, seq_len / 10000 + 50);
    let ts_capacity = get_thread_local_ts_capacity() as usize;
    let required_capacity = std::cmp::min(estimated_genes, ts_capacity);

    // Get thread-local gene storage with lazy allocation
    let ts = get_thread_local_ts_with_capacity(required_capacity);

    // Initialize only the required slots
    init_gene_ts(ts, 0, required_capacity as i32);

    let mut nt = 0i32;
    let psmax = seq_record.len() as i64;

    // Allocate buffers (thread-local, on stack/heap)
    let buf_size = 2 * LSEQ + WRAP + 1;
    let mut seq_buf: Vec<i32> = vec![NOBASE; buf_size];
    let mut cseq_buf: Vec<i32> = vec![NOBASE; buf_size];
    let mut wseq: Vec<i32> = vec![NOBASE; 2 * WRAP + 1];

    // seq_len already defined above for capacity estimation
    let mut se_idx: usize = 0;
    let mut read_pos: usize = 0;
    let mut start: i64 = 1;
    let mut flag: i32 = 0;
    let mut skip_nx: i32 = 0;

    // Create a minimal data_set for this thread (only fields used by tmioptimise_ts)
    let mut local_d = DataSet::default();
    local_d.psmax = psmax;
    local_d.gc = gc;
    // CRITICAL: Set seqname so that gene names are properly set during detection
    // This is copied to gene.name by copy3cr in trna.rs when genes are found
    let name_bytes = seq_record.id.as_bytes();
    let copy_len = std::cmp::min(name_bytes.len(), STRLEN - 1);
    for (i, &b) in name_bytes[..copy_len].iter().enumerate() {
        local_d.seqname[i] = b;
    }
    local_d.seqname[copy_len] = 0;
    let d_ptr = &mut local_d as *mut DataSet as *mut data_set;

    // CRITICAL: Create a thread-local copy of csw to avoid race conditions
    // tmioptimise_ts modifies start/comp, so each thread needs its own copy
    let mut local_sw_copy: Csw = std::ptr::read(sw_ptr as *const Csw);
    // Clear the file pointer since we don't output during parallel processing
    local_sw_copy.f = std::ptr::null_mut();
    // IMPORTANT: Update genespace to match actually allocated capacity
    // This is used by find_slot_ts and other functions that check gene storage limits
    local_sw_copy.genespace = required_capacity as i32;
    // CRITICAL: Apply the effective loffset/roffset passed from caller
    // This matches bopt_fastafile behavior where loffset/roffset persist across sequences
    // once a small circular sequence is encountered
    local_sw_copy.loffset = sw_loffset;
    local_sw_copy.roffset = sw_roffset;
    let local_sw = &mut local_sw_copy as *mut Csw as *mut csw;

    // Initialize based on topology (linear vs circular)
    if sw_linear != 0 {
        for _ in 0..rewind as usize {
            seq_buf[se_idx] = NOBASE;
            se_idx += 1;
        }
        start -= rewind;
    } else {
        if psmax <= drewind {
            let gap = (drewind - psmax) as usize;
            for _ in 0..gap {
                seq_buf[se_idx] = NOBASE;
                se_idx += 1;
            }
            let mut swrap_idx: usize = 0;
            for i in 0..psmax as usize {
                let base = convert_base(seq_record.seq[i]);
                seq_buf[se_idx] = base;
                wseq[swrap_idx] = base;
                se_idx += 1;
                swrap_idx += 1;
            }
            for _ in 0..gap {
                wseq[swrap_idx] = NOBASE;
                swrap_idx += 1;
            }
            swrap_idx = 0;
            for _ in 0..psmax as usize {
                seq_buf[se_idx] = wseq[swrap_idx];
                se_idx += 1;
                swrap_idx += 1;
            }
            swrap_idx = 0;
            for _ in 0..drewind as usize {
                seq_buf[se_idx] = wseq[swrap_idx];
                se_idx += 1;
                swrap_idx += 1;
            }
            // Update loffset/roffset for small circular sequences (matching bopt_fastafile)
            (*local_sw).loffset = drewind as i32;
            (*local_sw).roffset = drewind as i32;
            start -= drewind;
            flag = 1;
            skip_nx = 1;
            read_pos = seq_len;
        } else {
            let mut swrap_idx: usize = 0;
            for i in 0..drewind as usize {
                let base = convert_base(seq_record.seq[i]);
                seq_buf[se_idx] = base;
                wseq[swrap_idx] = base;
                se_idx += 1;
                swrap_idx += 1;
                read_pos += 1;
            }
            skip_nx = 0;
        }
    }

    let sc_limit: usize = LSEQ;

    // Main chunked processing loop
    // Structure matches bopt_fastafile: write first, increment, then check end
    loop {
        if skip_nx == 0 {
            while se_idx < sc_limit {
                // Check if we've reached end of sequence FIRST (before writing)
                if read_pos >= seq_len {
                    if sw_linear != 0 {
                        for _ in 0..rewind as usize {
                            if se_idx < buf_size {
                                seq_buf[se_idx] = NOBASE;
                                se_idx += 1;
                            }
                        }
                    } else {
                        for i in 0..drewind as usize {
                            if se_idx < buf_size {
                                seq_buf[se_idx] = wseq[i];
                                se_idx += 1;
                            }
                        }
                    }
                    flag = 1;
                    break;
                }

                // Only write base if we haven't reached end
                let base = convert_base(seq_record.seq[read_pos]);
                seq_buf[se_idx] = base;
                read_pos += 1;
                se_idx += 1;
            }
        }
        skip_nx = 0;

        let len = se_idx as i32;

        // Process forward strand
        if sw_both != 1 {
            // IMPORTANT: Use thread-local sw copy to avoid race conditions
            // Each thread needs its own copy since tmioptimise_ts modifies start/comp
            (*local_sw).start = start;
            (*local_sw).comp = 0;

            nt = tmioptimise_ts(ts, required_capacity as i32, d_ptr, seq_buf.as_mut_ptr(), len, nt, local_sw);
        }

        // Process reverse strand
        if sw_both > 0 {
            sense_switch(seq_buf.as_mut_ptr(), cseq_buf.as_mut_ptr(), len);

            // Use thread-local sw copy
            (*local_sw).start = start + len as i64;
            (*local_sw).comp = 1;

            nt = tmioptimise_ts(ts, required_capacity as i32, d_ptr, cseq_buf.as_mut_ptr(), len, nt, local_sw);
        }

        if flag == 0 {
            let sf_start = se_idx - drewind as usize;
            for i in 0..drewind as usize {
                seq_buf[i] = seq_buf[sf_start + i];
            }
            se_idx = drewind as usize;
            start += len as i64 - drewind;
            continue;
        }
        break;
    }

    // Copy detected genes from thread-local TS to result vector (as SendableGene)
    let mut genes: Vec<SendableGene> = Vec::with_capacity(nt as usize);
    for i in 0..nt {
        let gene_ptr = ts.offset(i as isize);
        if (*gene_ptr).genetype != NO_GENE {
            genes.push(SendableGene::from_gene(&*gene_ptr));
        }
    }

    (genes, nt)
}

/// Print interactive mode header (same as iopt_fastafile)
/// Called when sw.batch == 0 to match original ARAGORN output format
unsafe fn print_interactive_header(f: *mut File, sw: *mut csw) {
    static TRNATYPENAME: [[u8; 25]; 3] = [
        *b"Metazoan mitochondrial\0\0\0",
        *b"Cytosolic\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Mammalian mitochondrial\0\0",
    ];
    static GENECODENAME: [[u8; 50]; NGENECODE] = [
        *b"composite Metazoan mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"vertebrate mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"yeast mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"mold/protozoan/Coelenterate mitochondrial\0\0\0\0\0\0\0\0\0",
        *b"invertebrate mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Ciliate\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Echinoderm/Flatworm mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Euplotid\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"bacterial/plant chloroplast\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"alternative yeast\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Ascidian mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"alternative flatworm mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Blepharisma\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Chlorophycean mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"deleted -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Trematode mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Scenedesmus obliquus mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Thraustochytrium mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Pterobranchia mitochondrial\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Gracilibacteria\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Pachysolen tannophilus\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Karyorelict\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Condylostoma\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Mesodinium\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Peritrich\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Blastocrithidia\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"vacant -> standard\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",
        *b"Cephalodiscidae mitochondrial UAA-Tyr\0\0\0\0\0\0\0\0\0\0\0\0\0",
    ];

    let aragorn = (*sw).trna != 0 || (*sw).tmrna != 0 || (*sw).cds != 0 || (*sw).srprna != 0;

    // Initialize tmRNA
    init_tmrna(f, sw);

    // Print references
    let _ = write!(&mut *f, "\nPlease reference the following paper");
    if aragorn && (*sw).mtrna != 0 {
        let _ = write!(&mut *f, "s");
    }
    let _ = writeln!(&mut *f, " if you use this");
    let _ = writeln!(&mut *f, "program as part of any published research.\n");

    if aragorn {
        let _ = writeln!(&mut *f, "Laslett, D. and Canback, B. (2004) ARAGORN, a");
        let _ = writeln!(&mut *f, "program for the detection of transfer RNA and");
        let _ = writeln!(&mut *f, "transfer-messenger RNA genes in nucleotide sequences.");
        let _ = writeln!(&mut *f, "Nucleic Acids Research, 32;11-16.\n");
    }

    if (*sw).mtrna != 0 {
        let _ = writeln!(&mut *f, "Laslett, D. and Canback, B. (2008) ARWEN: a");
        let _ = writeln!(&mut *f, "program to detect tRNA genes in metazoan mitochondrial");
        let _ = writeln!(&mut *f, "nucleotide sequences");
        let _ = writeln!(&mut *f, "Bioinformatics, 24(2); 172-175.\n\n");
    }

    // Print search configuration
    let _ = writeln!(&mut *f);

    if (*sw).mtrna != 0 {
        let trnatypename_rust = CStr::from_ptr(TRNATYPENAME[(*sw).discrim as usize].as_ptr() as *const i8).to_string_lossy();
        let _ = writeln!(&mut *f, "Searching for {} tRNA genes", trnatypename_rust);
        if (*sw).tvloop == 0 {
            let _ = writeln!(&mut *f, "TV replacement loop tRNA genes not detected");
        }
    } else if (*sw).trna != 0 {
        let _ = write!(&mut *f, "Searching for tRNA genes");
        if (*sw).maxintronlen > 0 {
            let _ = write!(&mut *f, " with introns in anticodon loop");
        } else {
            let _ = write!(&mut *f, " with no introns");
        }
        let _ = writeln!(&mut *f);
        if (*sw).maxintronlen > 0 {
            let _ = writeln!(&mut *f, "Intron length from {} to {} bases",
                (*sw).minintronlen, (*sw).maxintronlen);
            if (*sw).ifixedpos != 0 {
                let _ = writeln!(&mut *f, "Intron position fixed between positions 37 and 38");
                let _ = writeln!(&mut *f, "on C-loop (one base after anticodon)");
            }
            if (*sw).ioverlay != 0 {
                let _ = writeln!(&mut *f, "Allowing overlay of long tRNA genes");
            }
        }
    }

    if (*sw).tmrna != 0 {
        let _ = writeln!(&mut *f, "Searching for tmRNA genes");
    }

    if (*sw).linear != 0 {
        let _ = writeln!(&mut *f, "Assuming linear topology, search will not wrap around ends");
    } else {
        let _ = writeln!(&mut *f, "Assuming circular topology, search wraps around ends");
    }

    if (*sw).both == 2 {
        let _ = writeln!(&mut *f, "Searching both strands");
    } else if (*sw).both == 1 {
        let _ = writeln!(&mut *f, "Searching complementary (antisense) strand only");
    } else {
        let _ = writeln!(&mut *f, "Searching single (sense) strand only");
    }

    if (*sw).mtrna != 0 && (*sw).mtcompov != 0 {
        let _ = writeln!(&mut *f, "Reporting overlapping candidates on opposite strands");
    }

    if (*sw).mtrna != 0 || (*sw).trna != 0 || (*sw).tmrna != 0 {
        let genecodename_rust = CStr::from_ptr(GENECODENAME[(*sw).geneticcode as usize].as_ptr() as *const i8).to_string_lossy();
        let _ = writeln!(&mut *f, "Using {} genetic code", genecodename_rust);
        if (*sw).ngcmod > 0 {
            let _ = writeln!(&mut *f, "Specified modifications:");
            for i in 0..(*sw).ngcmod {
                let anticodon = (*sw).gcmod[i as usize];
                let c1 = cpbase(Thymine - (anticodon & 0x3)) as i8;
                let c2 = cpbase(Thymine - ((anticodon >> 2) & 0x3)) as i8;
                let c3 = cpbase(Thymine - ((anticodon >> 4) & 0x3)) as i8;
                let aaname_rust = CStr::from_ptr(AANAME[AAMAP[(*sw).geneticcode as usize][anticodon as usize] as usize].as_ptr() as *const i8).to_string_lossy();
                let _ = writeln!(&mut *f, "{}{}{} = {}",
                    c1 as u8 as char, c2 as u8 as char, c3 as u8 as char, aaname_rust);
            }
        }
    }

    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f);
}

/// Parallel FASTA file processing for both batch and interactive modes
/// Processes sequences using thread-local state for gene detection
/// Implements proper chunked processing for large sequences
/// Supports both batch (sw.batch != 0) and interactive (sw.batch == 0) output modes
pub unsafe fn parallel_fastafile(
    filepath: &str,
    d: *mut data_set,
    sw: *mut csw,
    num_threads: usize,
) {
    // NT increased to 2000 to handle large genomes without reallocation
    // Enable sequence-level parallelization
    let use_chunk_parallel_check = true;
    let use_seq_parallel_check = true;

    if !use_chunk_parallel_check && !use_seq_parallel_check {
        // Fall back to appropriate single-threaded function based on mode
        if (*sw).batch != 0 {
            bopt_fastafile(d, sw);
        } else {
            iopt_fastafile(d, sw);
        }
        return;
    }

    let f = (*sw).f;

    // Load all sequences into memory using needletail
    let sequences = match load_fasta_sequences(filepath) {
        Ok(seqs) => seqs,
        Err(e) => {
            let _ = writeln!(&mut *f, "Error loading FASTA file: {}", e);
            return;
        }
    };

    let total_seqs = sequences.len();
    if total_seqs == 0 {
        let _ = writeln!(&mut *f, "No sequences found in file");
        return;
    }

    // Print interactive mode header if not in batch mode
    if (*sw).batch == 0 {
        print_interactive_header(f, sw);
    }

    // Configure parallel processing
    // Use DEFAULT_MAX_THREADS from parallel.rs (configurable)
    let requested_threads = if num_threads == 0 { num_cpus::get() } else { num_threads };
    let effective_threads = std::cmp::min(requested_threads, DEFAULT_MAX_THREADS);
    // Enable parallel processing when threads > 1
    let parallel_enabled = num_threads != 1 && effective_threads > 1;

    if (*sw).verbose != 0 {
        eprintln!("Parallel mode: {} sequences, {} threads{}",
            total_seqs, effective_threads,
            if parallel_enabled { "" } else { " (sequential)" });
    }

    // Calculate rewind distance (same as original bopt_fastafile)
    let mut rewind: i64 = (MAXTAGDIST + 20) as i64;
    if (*sw).trna != 0 || (*sw).mtrna != 0 {
        let tmaxlen = (MAXTRNALEN as i32 + (*sw).maxintronlen) as i64;
        if rewind < tmaxlen {
            rewind = tmaxlen;
        }
    }
    if (*sw).tmrna != 0 {
        if rewind < MAXTMRNALEN as i64 {
            rewind = MAXTMRNALEN as i64;
        }
    }
    if (*sw).peptide != 0 {
        if (*sw).tagthresh >= 5 {
            if rewind < TSWEEP as i64 {
                rewind = TSWEEP as i64;
            }
        }
    }

    (*sw).loffset = rewind as i32;
    (*sw).roffset = rewind as i32;
    let drewind = 2 * rewind;

    // Pre-calculate GC content in parallel (embarrassingly parallel)
    // Uses global thread pool initialized in main.rs for reduced overhead
    let gc_values: Vec<f64> = if parallel_enabled {
        sequences.par_iter()
            .map(|seq| calculate_gc(&seq.seq))
            .collect()
    } else {
        sequences.iter()
            .map(|seq| calculate_gc(&seq.seq))
            .collect()
    };

    // Process sequences - either in parallel or sequentially
    (*d).ns = 0;
    (*d).nf = 0;

    // Capture read-only config values for parallel processing
    let sw_linear = (*sw).linear;
    let sw_both = (*sw).both;
    let sw_trna = (*sw).trna;
    let sw_mtrna = (*sw).mtrna;
    let sw_tmrna = (*sw).tmrna;
    let sw_loffset = (*sw).loffset;
    let sw_roffset = (*sw).roffset;
    let sw_maxintronlen = (*sw).maxintronlen;
    let sw_minintronlen = (*sw).minintronlen;
    let sw_verbose = (*sw).verbose;

    // CRITICAL: Find the first small circular sequence to match bopt_fastafile state behavior
    // In bopt_fastafile, loffset/roffset are set to drewind when a small circular sequence
    // (psmax <= drewind) is encountered, and this state PERSISTS for all subsequent sequences.
    // For parallel processing to produce identical results, we must replicate this behavior.
    let first_small_circular_idx: Option<usize> = sequences.iter()
        .position(|seq| (seq.len() as i64) <= drewind);

    // Also calculate loffset/roffset to use when drewind state is active
    let drewind_i32 = drewind as i32;

    // Convert sw pointer to usize for safe sharing across threads
    // (usize is Send + Sync, raw pointers are not)
    let sw_ptr_usize = sw as usize;

    // Single small sequence: fall back to sequential processing
    if total_seqs == 1 && sequences[0].len() < MIN_CHUNK_PARALLEL_LEN {
        if (*sw).batch != 0 {
            bopt_fastafile(d, sw);
        } else {
            iopt_fastafile(d, sw);
        }
        return;
    }

    // Calculate effective offsets for all sequences (handles small circular sequences)
    let effective_offsets: Vec<(i32, i32)> = sequences.iter()
        .enumerate()
        .map(|(idx, _)| {
            match first_small_circular_idx {
                Some(first_idx) if idx >= first_idx => (drewind_i32, drewind_i32),
                _ => (sw_loffset, sw_roffset),
            }
        })
        .collect();

    // UNIFIED CHUNK-BASED PARALLELISM
    // Create chunks from ALL sequences and process them in parallel
    // Optimized: fewer chunks (num_threads * 2), conditional deduplication
    let (all_chunks, needs_dedup) = create_unified_chunks(&sequences, &gc_values, &effective_offsets, effective_threads);
    let total_chunks = all_chunks.len();

    if (*sw).verbose != 0 {
        eprintln!("Using unified chunk parallelism: {} sequences -> {} chunks, {} threads",
            total_seqs, total_chunks, effective_threads);
    }

    // Process all chunks in parallel
    let chunk_results: Vec<(usize, Vec<SendableGene>)> = all_chunks.par_iter()
        .map(|chunk| {
            let sw_ptr = sw_ptr_usize as *const csw;
            let genes = process_chunk_parallel(
                chunk.data,
                chunk.global_offset,
                chunk.total_seq_len,
                chunk.gc,
                rewind,
                drewind,
                sw_linear,
                sw_both,
                sw_ptr,
                chunk.seq_name,
                chunk.effective_loffset,
                chunk.effective_roffset,
            );
            (chunk.seq_idx, genes)
        })
        .collect();

    // Group results by sequence index with conditional deduplication
    let results: Vec<SeqProcessResult> = {
        // Initialize result containers for each sequence
        let mut seq_genes: Vec<Vec<SendableGene>> = vec![Vec::new(); total_seqs];

        // Group genes by sequence index
        for (seq_idx, genes) in chunk_results {
            seq_genes[seq_idx].extend(genes);
        }

        // Create SeqProcessResult for each sequence
        sequences.iter()
            .enumerate()
            .map(|(idx, seq_record)| {
                let mut genes = std::mem::take(&mut seq_genes[idx]);

                // Only deduplicate if sequence was split into multiple chunks
                // Single-chunk sequences cannot have duplicates from chunk overlaps
                if needs_dedup[idx] && !genes.is_empty() {
                    genes.sort_by(|a, b| {
                        a.start.cmp(&b.start)
                            .then(a.stop.cmp(&b.stop))
                            .then(a.genetype.cmp(&b.genetype))
                            .then(a.comp.cmp(&b.comp))
                    });
                    genes.dedup_by(|a, b| {
                        a.start == b.start
                            && a.stop == b.stop
                            && a.genetype == b.genetype
                            && a.comp == b.comp
                    });
                }

                let nt = genes.len() as i32;
                SeqProcessResult {
                    idx,
                    id: seq_record.id.clone(),
                    psmax: seq_record.len() as i64,
                    gc: gc_values[idx],
                    genes,
                    nt,
                    found_genes: nt >= 1,
                }
            })
            .collect()
    };

    // Sort results by index to maintain output order
    let mut sorted_results = results;
    sorted_results.sort_by_key(|r| r.idx);

    // Output results sequentially (maintains correct order)
    for result in sorted_results {
        // Set up data_set for this sequence
        (*d).psmax = result.psmax;
        (*d).gc = result.gc;
        (*d).ps = 0;

        // Copy sequence name
        let name_bytes = result.id.as_bytes();
        let copy_len = std::cmp::min(name_bytes.len(), STRLEN - 1);
        for (i, &b) in name_bytes[..copy_len].iter().enumerate() {
            (*d).seqname[i] = b;
        }
        (*d).seqname[copy_len] = 0;

        if (*sw).verbose != 0 {
            eprintln!("{}", result.id);
            eprintln!("{} nucleotides in sequence", result.psmax);
            eprintln!("Mean G+C content = {:.1}%", 100.0 * result.gc);
        }

        // Output header - different format for interactive vs batch mode
        if (*sw).batch == 0 {
            // Interactive mode: name + stats, no > prefix
            let _ = writeln!(&mut *f, "{}", result.id);
            let _ = writeln!(&mut *f, "{} nucleotides in sequence", result.psmax);
            let _ = writeln!(&mut *f, "Mean G+C content = {:.1}%", 100.0 * result.gc);
        } else if (*sw).batch < 2 {
            // Batch mode: > prefix only
            let _ = writeln!(&mut *f, ">{}", result.id);
        }

        // Copy genes from result to global TS for output functions
        init_gene(0, NT as i32);
        let nt = result.nt;
        for (i, sendable_gene) in result.genes.iter().enumerate() {
            if i < NT {
                let ts_gene = TS.offset(i as isize);
                // Convert SendableGene back to Gene
                *ts_gene = sendable_gene.to_gene();
            }
        }

        if !result.found_genes {
            (*d).nf += 1;
        }

        // Remove overlapping tRNAs if intron detection is enabled
        if (*sw).maxintronlen > 0 {
            remove_overlapping_trna(d, nt, sw);
        }

        // Update tmRNA tag database if needed
        if (*sw).updatetmrnatags != 0 {
            update_tmrna_tag_database(TS, nt, sw);
        }

        // Output detected genes - choose format based on batch mode
        if (*sw).batch != 0 {
            batch_gene_set(d, nt, sw);
        } else {
            disp_gene_set(d, nt, sw);
        }

        if (*sw).verbose != 0 && (*sw).batch == 0 {
            let seqname_rust = std::ffi::CStr::from_ptr((*d).seqname.as_ptr() as *const i8).to_string_lossy();
            eprintln!("{}\nSearch Finished\n", seqname_rust);
        }

        (*d).ns += 1;
    }

    // Summary - different format for batch vs interactive mode
    if (*sw).batch != 0 {
        // Batch mode summary (same format as original bopt_fastafile)
        if (*d).ns > 1 && (*sw).batch < 2 {
            let _ = write!(&mut *f, ">end \t{} sequences", (*d).ns);
            if (*sw).trna != 0 || (*sw).mtrna != 0 {
                let _ = write!(&mut *f, " {} tRNA genes", (*sw).ngene[tRNA as usize]);
            }
            if (*sw).tmrna != 0 {
                let _ = write!(&mut *f, " {} tmRNA genes", (*sw).ngene[tmRNA as usize]);
            }
            if (*d).nf > 0 {
                let sens = 100.0 * ((*d).ns - (*d).nf) as f64 / (*d).ns as f64;
                let _ = write!(&mut *f, ", nothing found in {} sequences, ({:.2}% sensitivity)",
                    (*d).nf, sens);
            }
            let _ = writeln!(&mut *f);
        }
    } else {
        // Interactive mode summary (same format as original iopt_fastafile)
        if (*d).ns > 1 {
            let _ = writeln!(&mut *f, "\n\n{} sequences searched", (*d).ns);
            if ((*sw).trna | (*sw).mtrna) != 0 {
                let _ = writeln!(&mut *f, "Total tRNA genes = {}", (*sw).ngene[tRNA as usize]);
                if (*sw).matchacceptor != 0 {
                    let _ = writeln!(&mut *f, "Total iso-acceptor mismatches = {}", (*sw).iamismatch);
                }
            }
            if (*sw).tmrna != 0 {
                let _ = writeln!(&mut *f, "Total tmRNA genes = {}", (*sw).ngene[tmRNA as usize]);
            }
            if (*sw).reportpseudogenes != 0 && (*sw).nps > 0 {
                let _ = writeln!(&mut *f, "Total number of possible pseudogenes = {}", (*sw).nps);
            }
            if (*d).nf > 0 {
                let sens = 100.0 * ((*d).ns - (*d).nf) as f64 / (*d).ns as f64;
                let _ = writeln!(&mut *f, "Nothing found in {} sequences ({:.2}% sensitivity)",
                    (*d).nf, sens);
            }
        }
    }

    if (*sw).updatetmrnatags != 0 {
        report_new_tmrna_tags(sw);
    }
}

/// Output genes in batch format
unsafe fn output_batch_genes(f: *mut File, nt: i32, psmax: i64, sw: *mut csw) {
    for i in 0..nt {
        let t = TS.offset(i as isize);
        if (*t).genetype == NO_GENE {
            continue;
        }

        match (*t).genetype {
            g if g == TRNA => {
                disp_batch_trna(f, t, sw);
            }
            g if g == TMRNA => {
                disp_batch_tmrna(f, t, sw);
            }
            g if g == SRPRNA => {
                disp_batch_srprna(f, t, sw);
            }
            g if g == CDS => {
                disp_batch_cds(f, t, sw);
            }
            _ => {}
        }
    }
}
