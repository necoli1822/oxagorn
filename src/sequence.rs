/*
 * Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
 * sequence.rs - Sequence handling and conversion functions
 *
 * Based on ARAGORN v1.2.41 by Dean Laslett
 * Author: Sunju Kim <n.e.coli.1822@gmail.com>
 */

#![allow(unused_assignments)]
#![allow(unused_variables)]
#![allow(unused_mut)]
#![allow(non_snake_case)]

use std::ptr;
use std::io::Write;
use std::ffi::CStr;
use crate::types::*;
use crate::tables::*;
use crate::utils::*;
use crate::tmrna_tags::TAGDATABASE;

/// Format float like C's %g (6 significant figures, removes trailing zeros)
fn format_g(val: f64) -> String {
    // C's %g uses 6 significant figures by default
    let s = format!("{:.6}", val);
    // Trim trailing zeros but keep at least one decimal digit for consistency
    let s = s.trim_end_matches('0');
    let s = s.trim_end_matches('.');
    // Handle edge case where we need more precision for small numbers
    // C's %g uses up to 6 significant figures
    if val != 0.0 && val.abs() < 0.1 {
        // For small numbers, use more decimal places
        let s2 = format!("{:.7}", val);
        let s2 = s2.trim_end_matches('0');
        let s2 = s2.trim_end_matches('.');
        return s2.to_string();
    }
    s.to_string()
}

// =============================================================================
// src_flatten/sequence.c lines 12-16
// char cbase(int c)
// =============================================================================
pub fn cbase(c: i32) -> u8 {
    static BASE: [u8; 6] = *b"acgt..";
    if c < Adenine {
        return b'#';
    }
    if c > NOBASE {
        return c as u8;
    }
    BASE[c as usize]
}

// =============================================================================
// src_flatten/sequence.c lines 19-23
// char cpbase(int c)
// =============================================================================
pub fn cpbase(c: i32) -> u8 {
    static BASE: [u8; 6] = *b"ACGT..";
    if c < Adenine {
        return b'#';
    }
    if c > NOBASE {
        return c as u8;
    }
    BASE[c as usize]
}

// =============================================================================
// src_flatten/sequence.c lines 26-31
// char *aa(int *anticodon, csw *sw)
// =============================================================================
pub unsafe fn aa(anticodon: *mut i32, sw: *mut csw) -> *const u8 {
    let p1: i32;
    let p2: i32;
    let p3: i32;
    p1 = *anticodon;
    if p1 >= AMBIG {
        return AMBIG_AANAME.as_ptr();
    }
    p2 = *anticodon.add(1);
    if p2 >= AMBIG {
        return AMBIG_AANAME.as_ptr();
    }
    p3 = *anticodon.add(2);
    if p3 >= AMBIG {
        return AMBIG_AANAME.as_ptr();
    }
    let idx = ((p1 << 4) + (p2 << 2) + p3) as usize;
    let aa_idx = AAMAP[(*sw).geneticcode as usize][idx] as usize;
    AANAME[aa_idx].as_ptr()
}

// =============================================================================
// src_flatten/sequence.c lines 34-41
// char *translate(int *codon, csw *sw)
// =============================================================================
pub unsafe fn translate(codon: *mut i32, sw: *mut csw) -> *const u8 {
    let p1: i32;
    let p2: i32;
    let p3: i32;
    let mut aac: i32;
    p1 = *codon;
    if p1 >= AMBIG {
        return AMBIG_AANAME.as_ptr();
    }
    p2 = *codon.add(1);
    if p2 >= AMBIG {
        return AMBIG_AANAME.as_ptr();
    }
    p3 = *codon.add(2);
    if p3 >= AMBIG {
        return AMBIG_AANAME.as_ptr();
    }
    let idx = (((3 - p3) << 4) + ((3 - p2) << 2) + (3 - p1)) as usize;
    aac = AAMAP[(*sw).geneticcode as usize][idx];
    if aac == SEC || aac == PYL {
        aac = STOP;
    }
    AANAME[aac as usize].as_ptr()
}

// =============================================================================
// src_flatten/sequence.c lines 44-51
// char ltranslate(int *codon, gene *t, csw *sw)
// =============================================================================
pub unsafe fn ltranslate(codon: *mut i32, t: *mut gene, sw: *mut csw) -> u8 {
    let code: i32;
    let p1: i32;
    let p2: i32;
    let p3: i32;
    if (*t).genetype == CDS {
        code = (*t).asst;
    } else {
        code = (*sw).geneticcode;
    }
    p1 = *codon;
    if p1 >= AMBIG {
        return AMBIG_AANAME[0];
    }
    p2 = *codon.add(1);
    if p2 >= AMBIG {
        return AMBIG_AANAME[0];
    }
    p3 = *codon.add(2);
    if p3 >= AMBIG {
        return AMBIG_AANAME[0];
    }
    let idx = (((3 - p3) << 4) + ((3 - p2) << 2) + (3 - p1)) as usize;
    let aa_idx = AAMAP[code as usize][idx] as usize;
    AALETTER[aa_idx]
}

// =============================================================================
// src_flatten/sequence.c lines 54-59
// char ptranslate(int *codon, csw *sw)
// =============================================================================
pub unsafe fn ptranslate(codon: *mut i32, sw: *mut csw) -> u8 {
    let p1: i32;
    let p2: i32;
    let p3: i32;
    p1 = *codon;
    if p1 >= AMBIG {
        return AMBIG_AANAME[0];
    }
    p2 = *codon.add(1);
    if p2 >= AMBIG {
        return AMBIG_AANAME[0];
    }
    p3 = *codon.add(2);
    if p3 >= AMBIG {
        return AMBIG_AANAME[0];
    }
    let idx = (((3 - p3) << 4) + ((3 - p2) << 2) + (3 - p1)) as usize;
    let aa_idx = AAMAP[(*sw).geneticcode as usize][idx] as usize;
    AAPOLARITY[aa_idx]
}

// =============================================================================
// src_flatten/sequence.c lines 62-65
// int seqlen(gene *t)
// =============================================================================
pub unsafe fn seqlen(t: *mut gene) -> i32 {
    (*t).nbase + (*t).nintron
}

// =============================================================================
// src_flatten/sequence.c lines 68-77
// int aseqlen(data_set *d, annotated_gene *a)
// =============================================================================
pub unsafe fn aseqlen(d: *mut data_set, a: *mut annotated_gene) -> i32 {
    let alen: i32;
    let astart: i64;
    let mut astop: i64;
    astart = (*a).start;
    astop = (*a).stop;
    if astart > astop {
        astop += (*d).psmax;
    }
    alen = (astop - astart) as i32 + 1;
    alen
}

// =============================================================================
// src_flatten/sequence.c lines 80-96
// double gc_content(gene *t)
// =============================================================================
pub unsafe fn gc_content(t: *mut gene) -> f64 {
    let mut s: *mut i32;
    let mut se: *mut i32;
    let mut ngc: f64;
    static SCORE: [f64; 6] = [0.0, 1.0, 1.0, 0.0, 0.0, 0.0];
    ngc = 0.0;
    if (*t).nintron > 0 && (*t).asst == 0 {
        s = (*t).eseq.as_mut_ptr();
        se = s.add((*t).intron as usize);
        while s < se {
            ngc += SCORE[*s as usize];
            s = s.add(1);
        }
        s = se.add((*t).nintron as usize);
        se = (*t).eseq.as_mut_ptr().add(((*t).nbase + (*t).nintron) as usize);
        while s < se {
            ngc += SCORE[*s as usize];
            s = s.add(1);
        }
    } else {
        s = (*t).seq.as_mut_ptr();
        se = s.add((*t).nbase as usize);
        while s < se {
            ngc += SCORE[*s as usize];
            s = s.add(1);
        }
    }
    ngc / (*t).nbase as f64
}

// =============================================================================
// src_flatten/sequence.c lines 99-108
// void write_seq(FILE *f, int *seq, int newline)
// =============================================================================
pub unsafe fn write_seq(f: *mut File, mut seq: *mut i32, newline: i32) {
    let mut i: i32;
    let mut c: i32;
    i = 0;
    loop {
        c = *seq;
        seq = seq.add(1);
        if c < Adenine {
            break;
        }
        let _ = write!(&mut *f, "{}", cbase(c) as char);
        if newline != 0 {
            i += 1;
            if i >= 50 {
                let _ = writeln!(&mut *f);
                i = 0;
            }
        }
    }
    if i > 0 {
        let _ = writeln!(&mut *f);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 111-170
// int find_var_hairpin(gene *t)
// =============================================================================
pub unsafe fn find_var_hairpin(t: *mut gene) -> i32 {
    let mut e: i32;
    let mut stem: i32;
    let mut vstem: i32 = 0;
    let mut loop_: i32;
    let mut sn: *mut i32;
    let mut sen: *mut i32;
    let mut pos1: *mut i32 = ptr::null_mut();
    let mut pos2: *mut i32 = ptr::null_mut();
    let mut sb: *mut i32;
    let mut se: *mut i32;
    let mut sc: *mut i32;
    let mut sd: *mut i32;
    let mut sf: *mut i32;
    let mut s: *mut i32;
    let mut c: u32;
    let mut cn: u32;
    let mut m: u32;

    static A: [u32; 6] = [0, 0, 0x100, 0x400, 0, 0];
    static C: [u32; 6] = [0, 0, 0x400, 0, 0, 0];
    static G: [u32; 6] = [0x100, 0x400, 0, 0x200, 0, 0];
    static T: [u32; 6] = [0x400, 0, 0x200, 0, 0, 0];
    let mut te: [u32; 6] = [0, 0, 0, 0, 0, 0];

    if (*t).genetype != tRNA {
        return 0;
    }
    if (*t).var < 13 {
        return 0;
    }
    e = 0;
    sb = (*t).seq.as_mut_ptr().add(
        ((*t).astem1 + (*t).spacer1 + 2 * (*t).dstem + (*t).dloop
            + (*t).spacer2 + 2 * (*t).cstem + (*t).cloop + (*t).nintron) as usize,
    );
    sc = sb.add(3);
    se = sb.add(((*t).var - 2) as usize);
    sf = se.sub(2);

    te[0] = A[*se as usize];
    te[1] = C[*se as usize];
    te[2] = G[*se as usize];
    te[3] = T[*se as usize];

    se = se.sub(1);
    while se > sf {
        te[0] = (te[0] >> 4) | A[*se as usize];
        te[1] = (te[1] >> 4) | C[*se as usize];
        te[2] = (te[2] >> 4) | G[*se as usize];
        te[3] = (te[3] >> 4) | T[*se as usize];
        se = se.sub(1);
    }

    while se >= sc {
        te[0] = (te[0] >> 4) | A[*se as usize];
        te[1] = (te[1] >> 4) | C[*se as usize];
        te[2] = (te[2] >> 4) | G[*se as usize];
        te[3] = (te[3] >> 4) | T[*se as usize];

        s = se.sub(5);
        sd = se.sub(7);
        m = te[*s as usize];
        s = s.sub(1);
        while s > sd {
            m = (m >> 4) + te[*s as usize];
            s = s.sub(1);
        }
        while s >= sb {
            m = (m >> 4) + te[*s as usize];
            c = m & 0xf;
            if c >= 9 {
                stem = 3;
                loop_ = (se.offset_from(s) as i32) - 3;
                sen = se;
                sn = s.add(2);
                while loop_ >= 6 {
                    cn = VBP[*sen.sub(1) as usize][*sn.add(1) as usize] as u32;
                    if cn == 0 {
                        break;
                    }
                    c += cn;
                    stem += 1;
                    loop_ -= 2;
                    sen = sen.sub(1);
                    sn = sn.add(1);
                }
                if c > e as u32 {
                    e = c as i32;
                    pos1 = s;
                    pos2 = sen;
                    vstem = stem;
                }
            }
            s = s.sub(1);
        }
        se = se.sub(1);
    }

    if e > 0 {
        return ((pos1.offset_from(sb) as i32) << 10)
            + ((pos2.offset_from(sb) as i32) << 5)
            + vstem;
    }
    0
}

// =============================================================================
// src_flatten/sequence.c lines 687-696
// void remove_inserts(int *s1, int *s2)
// =============================================================================
pub unsafe fn remove_inserts(mut s1: *mut i32, mut s2: *mut i32) {
    let mut flag: i32;
    let mut c: i32;
    flag = 0;
    loop {
        c = *s1;
        s1 = s1.add(1);
        if c == TERM {
            break;
        }
        if c == INSERT {
            flag = 1 - flag;
            continue;
        }
        if flag != 0 {
            continue;
        }
        *s2 = c;
        s2 = s2.add(1);
    }
    *s2 = TERM;
}

// =============================================================================
// src_flatten/sequence.c lines 699-712
// void sense_switch(int *seq1, int *seq2, int lseq)
// =============================================================================
pub unsafe fn sense_switch(seq1: *mut i32, seq2: *mut i32, lseq: i32) {
    let b: i32;
    let mut sseq: *mut i32;
    let mut cseq: *mut i32;
    sseq = seq1;
    cseq = seq2.add(lseq as usize);
    cseq = cseq.sub(1);
    while cseq >= seq2 {
        let b = *sseq;
        sseq = sseq.add(1);
        if b >= Adenine {
            if b <= Thymine {
                *cseq = Thymine - b;
            } else {
                if b <= NOBASE {
                    *cseq = b;
                } else {
                    *cseq = NOBASE;
                }
            }
        } else {
            *cseq = NOBASE;
        }
        cseq = cseq.sub(1);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 350-361
// int base_match(char b1, char b2)
// =============================================================================
pub fn base_match(b1: u8, b2: u8) -> i32 {
    static BASE1: [u8; 11] = *b"acgtgtagtg\0";
    static BASE2: [u8; 11] = *b"tgcatggatg\0";
    static SCORE: [i32; 11] = [2, 2, 2, 2, 1, 1, 3, 3, 3, 3, 0];
    let mut s: i32 = 0;
    for i in 0..10 {
        if b1 == BASE1[i] && b2 == BASE2[i] {
            s = SCORE[i];
            break;
        }
    }
    s
}

// =============================================================================
// src_flatten/sequence.c lines 713-721
// char *position(char *s, gene *t, csw *sw)
// =============================================================================
pub unsafe fn position(s: *mut u8, t: *mut gene, sw: *mut csw) -> *mut u8 {
    let mut start: i64 = (*t).start;
    if (*sw).linear != 0 {
        if start <= 0 {
            start -= 1;
        }
    }
    let result = if (*t).comp != 0 {
        format!("c[{},{}]", start, (*t).stop)
    } else {
        format!("[{},{}]", start, (*t).stop)
    };
    let bytes = result.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        *s.add(i) = b;
    }
    *s.add(bytes.len()) = 0;
    s.add(bytes.len())
}

// =============================================================================
// src_flatten/sequence.c lines 724-726
// void location(char *s, gene *t, csw *sw, char *m)
// =============================================================================
pub unsafe fn location(s: *mut u8, t: *mut gene, sw: *mut csw, m: *const u8) {
    let mut sp: [u8; 80] = [0; 80];
    position(sp.as_mut_ptr(), t, sw);
    // sprintf(s, "%s %s", m, position(sp, t, sw));
    let mut dst = s;
    let mut src = m;
    while *src != 0 {
        *dst = *src;
        dst = dst.add(1);
        src = src.add(1);
    }
    *dst = b' ';
    dst = dst.add(1);
    src = sp.as_ptr();
    while *src != 0 {
        *dst = *src;
        dst = dst.add(1);
        src = src.add(1);
    }
    *dst = 0;
}

// =============================================================================
// src_flatten/sequence.c lines 1277-1284
// int pseudogene(gene *t, csw *sw)
// =============================================================================
pub unsafe fn pseudogene(t: *mut gene, sw: *mut csw) -> i32 {
    if (*t).energy < (*sw).reportpsthresh {
        return 1;
    }
    if (*t).genetype == tRNA {
        if (*t).cloop != 7 {
            return 1;
        }
    }
    0
}

// =============================================================================
// src_flatten/sequence.c lines 1581-1595
// double stem_energy(int *s1, int *s2, int stem)
// =============================================================================
/// Optimized stem energy calculation - hot path function
/// SAFETY: Assumes s1 and s2 point to valid nucleotide values (0-5 range)
#[inline]
pub unsafe fn stem_energy(mut s1: *mut i32, mut s2: *mut i32, stem: i32) -> f64 {
    let se: *mut i32;
    let mut energy: f64;
    static BEM: [[f64; 6]; 6] = [
        [-1.072, -0.214, -1.072, ATBOND, 0.000, 0.000],
        [-0.214, -1.072, 3.000, -1.072, 0.000, 0.000],
        [-1.072, 3.000, -1.072, 1.286, 0.000, 0.000],
        [ATBOND, -1.072, 1.286, -0.214, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
    ];
    se = s1.add(stem as usize);
    s2 = s2.sub(1);
    // SAFETY: Nucleotide values are always in 0-5 range (validated at input)
    energy = *BEM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
    s1 = s1.add(1);
    while s1 < se {
        s2 = s2.sub(1);
        energy += *BEM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
        s1 = s1.add(1);
    }
    energy
}

// =============================================================================
// src_flatten/sequence.c lines 1597-1611
// double astem_energy(int *s1, int *s2, int stem)
// =============================================================================
/// Optimized acceptor stem energy calculation - hot path function
#[inline]
pub unsafe fn astem_energy(mut s1: *mut i32, mut s2: *mut i32, stem: i32) -> f64 {
    let se: *mut i32;
    let mut energy: f64;
    static ABEM: [[f64; 6]; 6] = [
        [-2.144, -0.428, -2.144, ATBOND, 0.000, 0.000],
        [-0.428, -2.144, 3.000, -2.144, 0.000, 0.000],
        [-2.144, 3.000, -2.144, 1.286, 0.000, 0.000],
        [ATBOND, -2.144, 1.286, -0.428, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
    ];
    se = s1.add(stem as usize);
    s2 = s2.sub(1);
    // SAFETY: Nucleotide values are always in 0-5 range
    energy = *ABEM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
    s1 = s1.add(1);
    while s1 < se {
        s2 = s2.sub(1);
        energy += *ABEM.get_unchecked(*s1 as usize).get_unchecked(*s2 as usize);
        s1 = s1.add(1);
    }
    energy
}

// =============================================================================
// src_flatten/sequence.c lines 1436-1493
// double vloop_stability(int *sb, int var, int *varbp)
// =============================================================================
pub unsafe fn vloop_stability(sb: *mut i32, var: i32, varbp: *mut i32) -> f64 {
    let mut e: i32;
    let mut stem: i32;
    let mut vstem: i32 = 0;
    let mut loop_: i32;
    let mut sn: *mut i32;
    let mut sen: *mut i32;
    let mut pos1: *mut i32 = ptr::null_mut();
    let mut pos2: *mut i32 = ptr::null_mut();
    let mut se: *mut i32;
    let sc: *mut i32;
    let mut sd: *mut i32;
    let sf: *mut i32;
    let mut s: *mut i32;
    let mut c: u32;
    let mut cn: u32;
    let mut m: u32;

    static A: [u32; 6] = [0, 0, 0x100, 0x400, 0, 0];
    static C: [u32; 6] = [0, 0, 0x400, 0, 0, 0];
    static G: [u32; 6] = [0x100, 0x400, 0, 0x200, 0, 0];
    static T: [u32; 6] = [0x400, 0, 0x200, 0, 0, 0];
    let mut te: [u32; 6] = [0, 0, 0, 0, 0, 0];

    e = 0;
    sc = sb.add(3);
    se = sb.add((var - 2) as usize);
    sf = se.sub(2);

    te[0] = A[*se as usize];
    te[1] = C[*se as usize];
    te[2] = G[*se as usize];
    te[3] = T[*se as usize];

    se = se.sub(1);
    while se > sf {
        te[0] = (te[0] >> 4) | A[*se as usize];
        te[1] = (te[1] >> 4) | C[*se as usize];
        te[2] = (te[2] >> 4) | G[*se as usize];
        te[3] = (te[3] >> 4) | T[*se as usize];
        se = se.sub(1);
    }

    while se >= sc {
        te[0] = (te[0] >> 4) | A[*se as usize];
        te[1] = (te[1] >> 4) | C[*se as usize];
        te[2] = (te[2] >> 4) | G[*se as usize];
        te[3] = (te[3] >> 4) | T[*se as usize];

        s = se.sub(5);
        sd = se.sub(7);
        m = te[*s as usize];
        s = s.sub(1);
        while s > sd {
            m = (m >> 4) + te[*s as usize];
            s = s.sub(1);
        }
        while s >= sb {
            m = (m >> 4) + te[*s as usize];
            c = m & 0xf;
            if c >= 9 {
                stem = 3;
                loop_ = (se.offset_from(s) as i32) - 3;
                sen = se;
                sn = s.add(2);
                while loop_ >= 6 {
                    cn = VBP[*sen.sub(1) as usize][*sn.add(1) as usize] as u32;
                    if cn == 0 {
                        break;
                    }
                    c += cn;
                    stem += 1;
                    loop_ -= 2;
                    sen = sen.sub(1);
                    sn = sn.add(1);
                }
                if c > e as u32 {
                    e = c as i32;
                    pos1 = s;
                    pos2 = sen;
                    vstem = stem;
                }
            }
            s = s.sub(1);
        }
        se = se.sub(1);
    }

    if e > 0 {
        *varbp = ((pos1.offset_from(sb) as i32) << 10)
            + ((pos2.offset_from(sb) as i32) << 5)
            + vstem;
        return (3 * (vstem - 4)) as f64;
    }
    *varbp = 0;
    -12.0
}

// =============================================================================
// src_flatten/sequence.c lines 1987-1997
// double nenergy(gene *t, csw *sw)
// =============================================================================
pub unsafe fn nenergy(t: *mut gene, sw: *mut csw) -> f64 {
    let eref: f64;
    if (*t).genetype != tRNA {
        eref = (*sw).eref[(*t).genetype as usize];
    } else {
        if (*sw).mtrna != 0 {
            if (*t).dstem == 0 {
                eref = MTRNA_T_THRESH;
            } else if (*t).tstem == 0 {
                eref = MTRNA_D_THRESH;
            } else {
                eref = MTRNA_DT_THRESH;
            }
        } else {
            eref = (*sw).eref[tRNA as usize];
        }
    }
    100.0 * (*t).energy / eref
}

// =============================================================================
// src_flatten/sequence.c lines 2122-2125
// void init_matrix(char matrix[][MATY])
// =============================================================================
pub fn init_matrix(matrix: &mut [[u8; MATY]; MATX]) {
    for i in 0..MATY {
        for j in 0..MATX {
            matrix[j][i] = b' ';
        }
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1794-1802
// void xcopy(char m[][MATY], int x, int y, char *s, int l)
// =============================================================================
pub unsafe fn xcopy(m: *mut [u8; MATY], mut x: i32, y: i32, mut s: *const u8, l: i32) {
    let mut i: i32 = 0;
    while i < l {
        if x >= MATX as i32 {
            break;
        }
        let c = *s;
        if c == 0 {
            break;
        }
        s = s.add(1);
        (*m.add(x as usize))[y as usize] = c;
        x += 1;
        i += 1;
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1905-1915
// int string_compare(char *s1, char *s2)
// =============================================================================
pub unsafe fn string_compare(mut s1: *const u8, mut s2: *const u8) -> i32 {
    let mut r: i32 = 0;
    loop {
        let c1 = *s1;
        if c1 == 0 {
            break;
        }
        s1 = s1.add(1);
        let c2 = *s2;
        if c2 == 0 {
            break;
        }
        s2 = s2.add(1);
        r = upcasec(c1) as i32 - upcasec(c2) as i32;
        if r != 0 {
            break;
        }
    }
    r
}

// =============================================================================
// src_flatten/sequence.c lines 1840-1855
// int peptide_tag(char tag[], int maxlen, gene *t, csw *sw)
// =============================================================================
pub unsafe fn peptide_tag(tag: *mut u8, maxlen: i32, t: *mut gene, sw: *mut csw) -> i32 {
    let mut i: i32;
    let mut lx: i32;
    let mut se: *mut i32;
    se = (*t).eseq.as_mut_ptr().add((*t).tps as usize);
    lx = ((*t).tpe - (*t).tps + 1) as i32;
    if ltranslate(se.add(lx as usize), t, sw) == b'*' {
        lx += 3;
        if ltranslate(se.add(lx as usize), t, sw) == b'*' {
            lx += 3;
        }
    }
    lx /= 3;
    if lx > maxlen {
        lx = maxlen;
    }
    i = 0;
    while i < lx {
        *tag.add(i as usize) = ltranslate(se, t, sw);
        se = se.add(3);
        i += 1;
    }
    *tag.add(i as usize) = 0;
    lx
}

// =============================================================================
// src_flatten/sequence.c lines 254-347
// int *make_tv(int *seq, char matrix[][MATY], int *x, int *y, int orient, int tv)
// =============================================================================
pub unsafe fn make_tv(mut seq: *mut i32, matrix: *mut [u8; MATY], x: *mut i32, y: *mut i32, orient: i32, mut tv: i32) -> *mut i32 {
    let mut i: i32;
    let mut px: i32;
    let mut py: i32;
    let mut stem: i32;
    static UX: [i32; 4] = [1, 0, -1, 0];
    static UY: [i32; 4] = [0, 1, 0, -1];
    static VX: [i32; 4] = [0, -1, 0, 1];
    static VY: [i32; 4] = [1, 0, -1, 0];
    static LOOPU: [[i32; 26]; 26] = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    ];
    static LOOPV: [[i32; 26]; 26] = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ];

    px = *x;
    py = *y;
    stem = 0;
    if tv < 6 {
        px += UX[orient as usize];
        py += UY[orient as usize];
        i = 0;
        while i < tv {
            px += VX[orient as usize];
            py += VY[orient as usize];
            (*matrix.add(px as usize))[py as usize] = cbase(*seq);
            seq = seq.add(1);
            i += 1;
        }
        py += (6 - i) * VY[orient as usize];
    } else {
        if tv > 25 {
            if tv % 2 != 0 {
                stem = (tv - 25) / 2;
            } else {
                stem = (tv - 24) / 2;
            }
            tv = tv - 2 * stem;
        }
        i = 0;
        while i < stem {
            px += UX[orient as usize];
            py += UY[orient as usize];
            (*matrix.add(px as usize))[py as usize] = cbase(*seq);
            seq = seq.add(1);
            i += 1;
        }
        i = 0;
        while i < tv {
            px += UX[orient as usize] * LOOPU[tv as usize][i as usize] + VX[orient as usize] * LOOPV[tv as usize][i as usize];
            py += UY[orient as usize] * LOOPU[tv as usize][i as usize] + VY[orient as usize] * LOOPV[tv as usize][i as usize];
            (*matrix.add(px as usize))[py as usize] = cbase(*seq);
            seq = seq.add(1);
            i += 1;
        }
        px += UX[orient as usize] * LOOPU[tv as usize][i as usize] + VX[orient as usize] * LOOPV[tv as usize][i as usize];
        py += UY[orient as usize] * LOOPU[tv as usize][i as usize] + VY[orient as usize] * LOOPV[tv as usize][i as usize];
        i = 0;
        while i < stem {
            (*matrix.add(px as usize))[py as usize] = cbase(*seq);
            seq = seq.add(1);
            px -= UX[orient as usize];
            py -= UY[orient as usize];
            i += 1;
        }
    }
    *x = px;
    *y = py;
    seq
}

// =============================================================================
// src_flatten/sequence.c lines 364-515
// int *make_clover(int *seq, int b, int e, int stemlength,
//                   char matrix[][MATY], int *x, int *y, int orient)
// =============================================================================
pub unsafe fn make_clover(mut seq: *mut i32, b: i32, e: i32, stemlength: i32, matrix: *mut [u8; MATY], x: *mut i32, y: *mut i32, orient: i32) -> *mut i32 {
    let mut i: i32;
    let mut px: i32;
    let mut py: i32;
    let mut pxb: i32;
    let mut pyb: i32;
    let mut pxe: i32;
    let mut pye: i32;
    let mut l: i32;
    let mut xlg: i32;
    let mut xlgd: i32;
    let mut ylgh: i32;
    let mut ylg: i32;
    let mut noloop: i32;
    let mut s: *mut i32;
    let mut se: *mut i32;

    static UX: [i32; 9] = [1, 0, -1, 0, 0, 1, 1, -1, -1];
    static UY: [i32; 9] = [0, 1, 0, -1, 1, -1, 1, 1, -1];
    static VX: [i32; 9] = [0, -1, 0, 1, 1, 1, 1, -1, -1];
    static VY: [i32; 9] = [1, 0, -1, 0, 0, 0, 0, 0, 0];
    static LOOPU: [[i32; 18]; 18] = [
        [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, -1, -1, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, 0],
        [0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1],
    ];
    static LOOPV: [[i32; 18]; 18] = [
        [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 2, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 1, 1, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 1, 1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 1, 1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 1, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 1, 1, 1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 1, 1, 1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 1, 1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, -1, 0, 0],
        [-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, -1, 0],
        [-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, -1],
    ];
    static DLOOPU: [[i32; 18]; 18] = [
        [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 0, 0, -1, -1, -1, -1, -1, -1, -1, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, 0],
        [0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1],
    ];
    static DLOOPV: [[i32; 18]; 18] = [
        [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 2, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 1, 1, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 1, 1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 1, 1, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 1, 1, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 1, 1, 1, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 1, 1, 1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 1, 1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, -1, 0, 0],
        [-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, -1, 0],
        [-1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, -1],
    ];
    static BOND1: [u8; 4] = *b" +!:";
    static BOND2: [u8; 4] = *b" +-.";

    px = *x;
    py = *y;
    s = seq.add(b as usize);
    se = s.add(stemlength as usize);
    while s < se {
        (*matrix.add(px as usize))[py as usize] = cbase(*s);
        s = s.add(1);
        px += UX[orient as usize];
        py += UY[orient as usize];
    }
    l = e - b - 2 * stemlength;
    if l < 0 {
        l = 0;
    }
    noloop = 0;
    if l < 18 {
        i = 0;
        if orient == DOWN {
            while i < l {
                px += UX[orient as usize] * DLOOPU[l as usize][i as usize] + VX[orient as usize] * DLOOPV[l as usize][i as usize];
                py += UY[orient as usize] * DLOOPU[l as usize][i as usize] + VY[orient as usize] * DLOOPV[l as usize][i as usize];
                (*matrix.add(px as usize))[py as usize] = cbase(*s);
                s = s.add(1);
                i += 1;
            }
            px += UX[orient as usize] * DLOOPU[l as usize][i as usize] + VX[orient as usize] * DLOOPV[l as usize][i as usize];
            py += UY[orient as usize] * DLOOPU[l as usize][i as usize] + VY[orient as usize] * DLOOPV[l as usize][i as usize];
        } else {
            while i < l {
                px += UX[orient as usize] * LOOPU[l as usize][i as usize] + VX[orient as usize] * LOOPV[l as usize][i as usize];
                py += UY[orient as usize] * LOOPU[l as usize][i as usize] + VY[orient as usize] * LOOPV[l as usize][i as usize];
                (*matrix.add(px as usize))[py as usize] = cbase(*s);
                s = s.add(1);
                i += 1;
            }
            px += UX[orient as usize] * LOOPU[l as usize][i as usize] + VX[orient as usize] * LOOPV[l as usize][i as usize];
            py += UY[orient as usize] * LOOPU[l as usize][i as usize] + VY[orient as usize] * LOOPV[l as usize][i as usize];
        }
    } else {
        ylgh = ((l >> 2) - 2) >> 1;
        ylg = (ylgh << 1) + 2;
        xlgd = l - ylg - 2 * ylgh + 1;
        xlg = (xlgd + 1) >> 1;
        pxb = px - ylgh * VX[orient as usize];
        pyb = py - ylgh * VY[orient as usize];
        pxe = px + xlg * UX[orient as usize] + (ylg - ylgh + 1) * VX[orient as usize];
        pye = py + xlg * UY[orient as usize] + (ylg - ylgh + 1) * VY[orient as usize];
        if pxb < 0 || pxb >= MATX as i32 || pyb < 0 || pyb >= MATY as i32 ||
           pxe < 0 || pxe >= MATX as i32 || pye < 0 || pye >= MATY as i32 {
            noloop = 1;
        } else {
            for _ii in 0..ylgh {
                px -= VX[orient as usize];
                py -= VY[orient as usize];
                (*matrix.add(px as usize))[py as usize] = cbase(*s);
                s = s.add(1);
            }
            for _ii in 0..xlg {
                px += UX[orient as usize];
                py += UY[orient as usize];
                (*matrix.add(px as usize))[py as usize] = cbase(*s);
                s = s.add(1);
            }
            for _ii in 1..ylg {
                px += VX[orient as usize];
                py += VY[orient as usize];
                (*matrix.add(px as usize))[py as usize] = cbase(*s);
                s = s.add(1);
            }
            px += VX[orient as usize];
            py += VY[orient as usize];
            if (xlgd & 1) == 0 {
                (*matrix.add(px as usize))[py as usize] = cbase(*s);
                s = s.add(1);
            }
            for _ii in 0..xlg {
                px -= UX[orient as usize];
                py -= UY[orient as usize];
                (*matrix.add(px as usize))[py as usize] = cbase(*s);
                s = s.add(1);
            }
            for _ii in 1..ylgh {
                px -= VX[orient as usize];
                py -= VY[orient as usize];
                (*matrix.add(px as usize))[py as usize] = cbase(*s);
                s = s.add(1);
            }
            px -= UX[orient as usize] + VX[orient as usize];
            py -= UY[orient as usize] + VY[orient as usize];
        }
    }
    if noloop != 0 {
        px += UX[orient as usize] * LOOPU[0][0] + VX[orient as usize] * LOOPV[0][0];
        py += UY[orient as usize] * LOOPU[0][0] + VY[orient as usize] * LOOPV[0][0];
    }
    se = seq.add(e as usize);
    s = se.sub(stemlength as usize);
    while s < se {
        (*matrix.add(px as usize))[py as usize] = cbase(*s);
        s = s.add(1);
        let bi = base_match(
            (*matrix.add(px as usize))[py as usize],
            (*matrix.add((px - 2 * VX[orient as usize]) as usize))[(py - 2 * VY[orient as usize]) as usize]
        );
        match orient {
            RIGHT | LEFT => {
                (*matrix.add((px - VX[orient as usize]) as usize))[(py - VY[orient as usize]) as usize] = BOND1[bi as usize];
            }
            _ => {
                (*matrix.add((px - VX[orient as usize]) as usize))[(py - VY[orient as usize]) as usize] = BOND2[bi as usize];
            }
        }
        px -= UX[orient as usize];
        py -= UY[orient as usize];
    }
    *x = px;
    *y = py;
    se
}

// =============================================================================
// src_flatten/sequence.c lines 518-586
// int *make_dv(int *seq, char matrix[][MATY], int dloop,
//                   int orient, int *xp, int *yp)
// =============================================================================
pub unsafe fn make_dv(mut seq: *mut i32, matrix: *mut [u8; MATY], dloop: i32, orient: i32, xp: *mut i32, yp: *mut i32) -> *mut i32 {
    let mut i: i32;
    let mut x: i32;
    let mut y: i32;
    static UX: [i32; 5] = [1, 0, -1, 0, 0];
    static UY: [i32; 5] = [0, 1, 0, -1, 1];
    static VX: [i32; 5] = [0, -1, 0, 1, 1];
    static VY: [i32; 5] = [1, 0, -1, 0, 0];
    static LOOPU: [[i32; 22]; 22] = [
        [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, 0, -1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    ];
    static LOOPV: [[i32; 22]; 22] = [
        [-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-2, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-2, -1, -1, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -2, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, -1, 0, -1, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, 0, -1, -1, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, 0, -1, -1, -1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, -1, 0, -1, -1, -1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, -1, 0, -1, -1, -1, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, -1, 0, -1, -1, -1, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, -1, 0, -1, -1, -1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, -1, 0, -1, -1, -1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, -1, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, -1, 0, -1, -1, -1, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, -1, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, -1, 0, -1, -1, -1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, -1, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, -1, -1, 0, -1, 0, 0, 0, 0, 0, 0, -1],
    ];

    x = *xp;
    y = *yp;
    if dloop < 2 || dloop > 21 {
        x -= 1;
        y -= 6;
        seq = seq.add(dloop as usize);
    } else {
        i = 0;
        while i < dloop {
            x += UX[orient as usize] * LOOPU[dloop as usize][i as usize] + VX[orient as usize] * LOOPV[dloop as usize][i as usize];
            y += UY[orient as usize] * LOOPU[dloop as usize][i as usize] + VY[orient as usize] * LOOPV[dloop as usize][i as usize];
            (*matrix.add(x as usize))[y as usize] = cbase(*seq);
            seq = seq.add(1);
            i += 1;
        }
        x += UX[orient as usize] * LOOPU[dloop as usize][i as usize] + VX[orient as usize] * LOOPV[dloop as usize][i as usize];
        y += UY[orient as usize] * LOOPU[dloop as usize][i as usize] + VY[orient as usize] * LOOPV[dloop as usize][i as usize];
    }
    *xp = x;
    *yp = y;
    seq
}

// =============================================================================
// src_flatten/sequence.c lines 589-684
// int *make_var(int *seq, char matrix[][MATY],
//                int *x, int *y, int orient, int var, int varbp)
// =============================================================================
pub unsafe fn make_var(mut seq: *mut i32, matrix: *mut [u8; MATY], x: *mut i32, y: *mut i32, orient: i32, mut var: i32, varbp: i32) -> *mut i32 {
    let mut i: i32;
    let mut b: i32;
    let e: i32;
    let p: i32;
    let mut px: i32;
    let mut py: i32;
    let pxf: i32;
    let pyf: i32;
    let mut l: i32;
    let mut stem: i32;
    let mut skip_nbp: i32;
    static UX: [i32; 4] = [1, 0, -1, 0];
    static UY: [i32; 4] = [0, 1, 0, -1];
    static VX: [i32; 4] = [0, -1, 0, 1];
    static VY: [i32; 4] = [1, 0, -1, 0];
    static PREU: [[[i32; 4]; 5]; 4] = [
        [[0, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]],
        [[0, 0, 0, 0], [1, 1, 0, 0], [1, 1, 0, 0], [1, 1, 0, 0], [1, 1, 0, 0]],
        [[0, 0, 0, 0], [0, 1, 1, 0], [0, 1, 1, 0], [1, 1, 1, 0], [1, 1, 1, 0]],
        [[0, 0, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 1], [0, 1, 1, 1]],
    ];
    static PREV: [[[i32; 4]; 5]; 4] = [
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [-1, 0, 0, 0]],
        [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [-1, 0, 0, 0], [-1, 0, 0, 0]],
        [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [-1, 0, 0, 0], [-1, 0, -1, 0]],
        [[0, 0, 0, 0], [0, 1, 0, 0], [0, 1, 0, -1], [0, 0, 0, -1], [0, 0, -1, -1]],
    ];
    static POSTU: [[[i32; 4]; 5]; 4] = [
        [[0, 0, 0, 0], [0, 0, 0, 0], [1, -1, 0, 0], [0, -1, 0, 0], [1, -1, -1, 0]],
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, -1, 0, 0], [0, -1, 0, 0], [0, -1, -1, 0]],
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, -1, 0, 0], [0, -1, -1, 0], [1, -1, -1, -1]],
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, -1, 0, 0], [0, -1, -1, 0], [1, -1, -1, -1]],
    ];
    static POSTV: [[[i32; 4]; 5]; 4] = [
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 1, 1]],
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 1, 1]],
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 1, 1]],
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 1, 0, 0], [0, 1, 1, 0], [0, 1, 1, 1]],
    ];
    static LOOPU: [[i32; 10]; 10] = [
        [2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, -1, -1, 1, 0, 0, 0, 0],
        [1, 1, 1, 0, 0, -1, 0, 0, 0, 0],
        [1, 1, 1, 1, 0, -1, -1, 0, 0, 0],
        [1, 1, 1, 1, 0, 0, -1, -1, 0, 0],
        [1, 1, 1, 1, 1, 0, -1, -1, -1, 0],
    ];
    static LOOPV: [[i32; 10]; 10] = [
        [3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 2, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 2, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        [-1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
        [-1, 0, 1, 1, 1, 1, 0, 0, 0, 0],
        [-1, -1, 1, 1, 1, 1, 1, 0, 0, 0],
        [-1, -1, 0, 1, 1, 1, 1, 1, 0, 0],
        [-1, -1, -1, 1, 1, 1, 1, 1, 1, 0],
        [-1, -1, -1, 0, 1, 1, 1, 1, 1, 1],
    ];

    px = *x;
    py = *y;
    if var < 0 {
        var = 0;
    }
    if var > 30 {
        var = 30;
    }
    skip_nbp = 0;
    if varbp > 0 {
        b = (varbp >> 10) & 0x1f;
        if b <= 3 {
            stem = varbp & 0x1f;
            e = stem + ((varbp >> 5) & 0x1f);
            p = var - e;
            if p >= 1 && p <= 4 {
                pxf = px + 2 * UX[orient as usize] + 3 * VX[orient as usize];
                pyf = py + 2 * UY[orient as usize] + 3 * VY[orient as usize];
                i = 0;
                while i < b {
                    px += UX[orient as usize] * PREU[b as usize][p as usize][i as usize] + VX[orient as usize] * PREV[b as usize][p as usize][i as usize];
                    py += UY[orient as usize] * PREU[b as usize][p as usize][i as usize] + VY[orient as usize] * PREV[b as usize][p as usize][i as usize];
                    (*matrix.add(px as usize))[py as usize] = cbase(*seq);
                    seq = seq.add(1);
                    i += 1;
                }
                px += UX[orient as usize] * PREU[b as usize][p as usize][b as usize] + VX[orient as usize] * PREV[b as usize][p as usize][b as usize];
                py += UY[orient as usize] * PREU[b as usize][p as usize][b as usize] + VY[orient as usize] * PREV[b as usize][p as usize][b as usize];
                seq = make_clover(seq, 0, e - b, stem, matrix, &mut px, &mut py, orient + SLANT);
                i = 0;
                while i < p {
                    px += UX[orient as usize] * POSTU[b as usize][p as usize][i as usize] + VX[orient as usize] * POSTV[b as usize][p as usize][i as usize];
                    py += UY[orient as usize] * POSTU[b as usize][p as usize][i as usize] + VY[orient as usize] * POSTV[b as usize][p as usize][i as usize];
                    (*matrix.add(px as usize))[py as usize] = cbase(*seq);
                    seq = seq.add(1);
                    i += 1;
                }
                *x = pxf;
                *y = pyf;
                skip_nbp = 1;
            }
        }
    }
    if skip_nbp == 0 {
        if var > 9 {
            if var % 2 != 0 {
                stem = (var - 7) / 2;
            } else {
                stem = (var - 6) / 2;
            }
        } else {
            stem = 0;
        }
        l = var - 2 * stem;
        i = 0;
        while i < stem {
            px += UX[orient as usize] - VX[orient as usize];
            py += UY[orient as usize] - VY[orient as usize];
            (*matrix.add(px as usize))[py as usize] = cbase(*seq);
            seq = seq.add(1);
            i += 1;
        }
        i = 0;
        while i < l {
            px += UX[orient as usize] * LOOPU[l as usize][i as usize] + VX[orient as usize] * LOOPV[l as usize][i as usize];
            py += UY[orient as usize] * LOOPU[l as usize][i as usize] + VY[orient as usize] * LOOPV[l as usize][i as usize];
            (*matrix.add(px as usize))[py as usize] = cbase(*seq);
            seq = seq.add(1);
            i += 1;
        }
        px += UX[orient as usize] * LOOPU[l as usize][i as usize] + VX[orient as usize] * LOOPV[l as usize][i as usize];
        py += UY[orient as usize] * LOOPU[l as usize][i as usize] + VY[orient as usize] * LOOPV[l as usize][i as usize];
        i = 0;
        while i < stem {
            (*matrix.add(px as usize))[py as usize] = cbase(*seq);
            seq = seq.add(1);
            px -= UX[orient as usize] - VX[orient as usize];
            py -= UY[orient as usize] - VY[orient as usize];
            i += 1;
        }
        *x = px;
        *y = py;
    }
    seq
}

// =============================================================================
// src_flatten/sequence.c lines 173-245
// void write_to_library(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn write_to_library(f: *mut File, t: *mut gene, sw: *mut csw) {
    let mut s: *mut i32;
    static TRNATYPE: [[u8; 6]; 2] = [*b"tRNA\0\0", *b"mtRNA\0"];
    s = (*t).seq.as_mut_ptr().add((*t).anticodon as usize);
    let _ = write!(&mut *f, ">{}", CStr::from_ptr((*t).name.as_ptr() as *const i8).to_string_lossy());
    if softstrpos((*t).name.as_mut_ptr(), b"RNA\0".as_ptr() as *mut u8).is_null() {
        match (*t).genetype {
            CDS => {
                let _ = write!(&mut *f, " CDS");
            }
            SRPRNA => {
                let _ = write!(&mut *f, " srpRNA");
            }
            TMRNA => {
                if (*t).asst > 0 {
                    let _ = write!(&mut *f, " Permuted");
                }
                let _ = write!(&mut *f, " tmRNA");
            }
            TRNA | _ => {
                (*t).varbp = find_var_hairpin(t);
                if (*t).tstem == 0 {
                    let _ = write!(&mut *f, " TV-loop");
                } else if (*t).dstem == 0 {
                    let _ = write!(&mut *f, " D-loop");
                }
                match (*t).cloop {
                    6 => {
                        let _ = write!(&mut *f, " {}-???({}{})",
                            CStr::from_ptr(TRNATYPE[(*sw).mtrna as usize].as_ptr() as *const i8).to_string_lossy(),
                            cbase(*s) as char,
                            cbase(*s.add(1)) as char);
                    }
                    8 => {
                        let _ = write!(&mut *f, " {}-???({}{}{}{})",
                            CStr::from_ptr(TRNATYPE[(*sw).mtrna as usize].as_ptr() as *const i8).to_string_lossy(),
                            cbase(*s) as char,
                            cbase(*s.add(1)) as char,
                            cbase(*s.add(2)) as char,
                            cbase(*s.add(3)) as char);
                    }
                    7 | _ => {
                        let _ = write!(&mut *f, " {}-{}({}{}{})",
                            CStr::from_ptr(TRNATYPE[(*sw).mtrna as usize].as_ptr() as *const i8).to_string_lossy(),
                            CStr::from_ptr(aa(s, sw) as *const i8).to_string_lossy(),
                            cbase(*s) as char,
                            cbase(*s.add(1)) as char,
                            cbase(*s.add(2)) as char);
                    }
                }
            }
        }
    }
    if !strpos((*t).name.as_mut_ptr(), b"bases)\0".as_ptr() as *mut u8).is_null() {
        let _ = writeln!(&mut *f);
    } else {
        let _ = writeln!(&mut *f, " ({} bases)", (*t).nbase);
    }
    let _ = writeln!(&mut *f, "sequence =");
    write_seq(f, (*t).seq.as_mut_ptr(), 1);
    if *(*t).eseq.as_ptr() >= Adenine {
        let _ = writeln!(&mut *f, "extended sequence =");
        write_seq(f, (*t).eseq.as_mut_ptr(), 1);
    }
    let _ = writeln!(&mut *f, "nbase = {}", (*t).nbase);
    let _ = writeln!(&mut *f, "sense = {}", (*t).comp);
    let _ = writeln!(&mut *f, "start = {}", (*t).start);
    let _ = writeln!(&mut *f, "stop = {}", (*t).stop);
    let _ = writeln!(&mut *f, "astem1 = {}", (*t).astem1);
    let _ = writeln!(&mut *f, "astem2 = {}", (*t).astem2);
    let _ = writeln!(&mut *f, "atail = {}", (*t).aatail);
    let _ = writeln!(&mut *f, "spacer1 = {}", (*t).spacer1);
    let _ = writeln!(&mut *f, "spacer2 = {}", (*t).spacer2);
    let _ = writeln!(&mut *f, "dstem = {}", (*t).dstem);
    let _ = writeln!(&mut *f, "dloop = {}", (*t).dloop);
    let _ = writeln!(&mut *f, "cstem = {}", (*t).cstem);
    let _ = writeln!(&mut *f, "cloop = {}", (*t).cloop);
    let _ = writeln!(&mut *f, "anticodon = {}", (*t).anticodon);
    let _ = writeln!(&mut *f, "nintron = {}", (*t).nintron);
    let _ = writeln!(&mut *f, "intron = {}", (*t).intron);
    let _ = write!(&mut *f, "asst = {}", (*t).asst);
    if (*t).genetype == tmRNA {
        if (*t).asst > 0 {
            let _ = write!(&mut *f, " permuted");
        }
    }
    let _ = writeln!(&mut *f, "\ntps = {}", (*t).tps);
    let _ = writeln!(&mut *f, "tpe = {}", (*t).tpe);
    let _ = writeln!(&mut *f, "var = {}", (*t).var);
    let _ = writeln!(&mut *f, "varbp = {},{},{}",
        ((*t).varbp >> 10) & 0x1f,
        ((*t).varbp >> 5) & 0x1f,
        (*t).varbp & 0x1f);
    let _ = writeln!(&mut *f, "tstem = {}", (*t).tstem);
    let _ = writeln!(&mut *f, "tloop = {}", (*t).tloop);
    let _ = writeln!(&mut *f, "gc = {}\n", format_g(gc_content(t)));
}

// =============================================================================
// src_flatten/sequence.c lines 733-797
// char *name(gene *t, char *si, int proc, csw *sw)
// =============================================================================
pub unsafe fn name(t: *mut gene, si: *mut u8, proc: i32, sw: *mut csw) -> *mut u8 {
    let mut s: [i32; 5] = [0; 5];
    let ss: *mut i32;
    let mut sin: *mut i32;
    let mut sm: *mut i32;
    let mut s0: *mut i32;
    let mut s1: *mut i32;
    let mut s2: *mut i32;
    let mut s3: *mut i32;
    let mut nintron: i32;
    let mut sb: *mut u8;
    static TRNATYPE: [[u8; 6]; 2] = [*b"tRNA\0\0", *b"mtRNA\0"];

    match (*t).genetype {
        CDS => {
            copy(b"CDS\0".as_ptr(), si);
        }
        SRPRNA => {
            copy(b"srpRNA\0".as_ptr(), si);
        }
        TMRNA => {
            if (*sw).dispmatch != 0 {
                if (*t).asst > 0 {
                    copy(b"tmRNA(Perm)  \0".as_ptr(), si);
                } else {
                    copy(b"tmRNA        \0".as_ptr(), si);
                }
            } else {
                if (*t).asst > 0 {
                    copy(b"tmRNA (Permuted)\0".as_ptr(), si);
                } else {
                    copy(b"tmRNA\0".as_ptr(), si);
                }
            }
        }
        TRNA | _ => {
            ss = if proc != 0 { (*t).seq.as_mut_ptr() } else { (*t).ps };
            sm = ss.add(((*t).anticodon - 1) as usize);
            s0 = sm.add(1);
            s1 = s0.add(1);
            s2 = s1.add(1);
            s3 = s2.add(1);
            nintron = (*t).nintron;
            if proc == 0 && nintron > 0 {
                sin = ss.add((*t).intron as usize);
                if sm >= sin { sm = sm.add(nintron as usize); }
                if s0 >= sin { s0 = s0.add(nintron as usize); }
                if s1 >= sin { s1 = s1.add(nintron as usize); }
                if s2 >= sin { s2 = s2.add(nintron as usize); }
                if s3 >= sin { s3 = s3.add(nintron as usize); }
            }
            s[0] = *sm;
            s[1] = *s0;
            s[2] = *s1;
            s[3] = *s2;
            s[4] = *s3;
            sb = si;
            if (*t).dstem == 0 {
                copy(b"D-loop \0".as_ptr(), sb);
                sb = sb.add(7);
            }
            if (*t).tstem == 0 {
                copy(b"TV-loop \0".as_ptr(), sb);
                sb = sb.add(8);
            }
            if (*t).cloop == 8 {
                let aa1 = aa(s.as_mut_ptr().add(1), sw);
                let aa2 = aa(s.as_mut_ptr().add(2), sw);
                let result = format!("{}-?({}|{})({}{}{}{})",
                    String::from_utf8_lossy(&TRNATYPE[(*sw).mtrna as usize]).trim_matches('\0'),
                    std::ffi::CStr::from_ptr(aa1 as *const i8).to_str().unwrap_or(""),
                    std::ffi::CStr::from_ptr(aa2 as *const i8).to_str().unwrap_or(""),
                    cbase(s[1]) as char,
                    cbase(s[2]) as char,
                    cbase(s[3]) as char,
                    cbase(s[4]) as char);
                let bytes = result.as_bytes();
                for (i, &b) in bytes.iter().enumerate() {
                    *sb.add(i) = b;
                }
                *sb.add(bytes.len()) = 0;
            } else if (*t).cloop == 6 {
                let aa1 = aa(s.as_mut_ptr(), sw);
                let aa2 = aa(s.as_mut_ptr().add(1), sw);
                let result = format!("{}-?({}|{})({}{})",
                    String::from_utf8_lossy(&TRNATYPE[(*sw).mtrna as usize]).trim_matches('\0'),
                    std::ffi::CStr::from_ptr(aa1 as *const i8).to_str().unwrap_or(""),
                    std::ffi::CStr::from_ptr(aa2 as *const i8).to_str().unwrap_or(""),
                    cbase(s[1]) as char,
                    cbase(s[2]) as char);
                let bytes = result.as_bytes();
                for (i, &b) in bytes.iter().enumerate() {
                    *sb.add(i) = b;
                }
                *sb.add(bytes.len()) = 0;
            } else {
                let aa1 = aa(s.as_mut_ptr().add(1), sw);
                let result = format!("{}-{}({}{}{})",
                    String::from_utf8_lossy(&TRNATYPE[(*sw).mtrna as usize]).trim_matches('\0'),
                    std::ffi::CStr::from_ptr(aa1 as *const i8).to_str().unwrap_or(""),
                    cbase(s[1]) as char,
                    cbase(s[2]) as char,
                    cbase(s[3]) as char);
                let bytes = result.as_bytes();
                for (i, &b) in bytes.iter().enumerate() {
                    *sb.add(i) = b;
                }
                *sb.add(bytes.len()) = 0;
            }
        }
    }
    si
}

// =============================================================================
// src_flatten/sequence.c lines 857-885
// void disp_seq(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_seq(f: *mut File, t: *mut gene, sw: *mut csw) {
    let mut i: i32;
    let mut s: *mut i32;
    let mut se: *mut i32;
    let mut genename: [u8; 100] = [0; 100];

    if (*sw).seqdisp >= 3 {
        if (*sw).batch == 0 {
            let _ = writeln!(&mut *f);
        }
        // disp_fasta_seq would be called here
        // For now, just display the basic sequence
    } else {
        if (*sw).batch == 0 {
            name(t, genename.as_mut_ptr(), 1, sw);
            let _ = writeln!(&mut *f, "\nPrimary sequence for {}",
                CStr::from_ptr(genename.as_ptr() as *const i8).to_string_lossy());
            let _ = writeln!(&mut *f, "1   .   10    .   20    .   30    .   40    .   50");
        }
        if (*t).nintron > 0 {
            s = (*t).eseq.as_mut_ptr();
            se = s.add(((*t).nbase + (*t).nintron) as usize);
        } else {
            s = (*t).seq.as_mut_ptr();
            se = s.add((*t).nbase as usize);
        }
        i = 0;
        while s < se {
            let _ = write!(&mut *f, "{}", cbase(*s) as char);
            s = s.add(1);
            i += 1;
            if i >= 50 {
                let _ = writeln!(&mut *f);
                i = 0;
            }
        }
        if i > 0 {
            let _ = writeln!(&mut *f);
        }
    }
    if (*sw).batch == 0 {
        let _ = writeln!(&mut *f);
        let _ = writeln!(&mut *f);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 2128-2136
// void disp_matrix(FILE *f, char matrix[][MATY], int ylines)
// =============================================================================
pub unsafe fn disp_matrix(f: *mut File, matrix: *const [u8; MATY], ylines: i32) {
    let mut i: i32;
    let mut j: i32;
    let mut k: i32;
    i = ylines;
    while { i -= 1; i >= 0 } {
        k = MATX as i32;
        while { k -= 1; k > 0 } {
            if (*matrix.add(k as usize))[i as usize] != b' ' {
                break;
            }
        }
        j = 0;
        while j <= k {
            let _ = write!(&mut *f, "{}", (*matrix.add(j as usize))[i as usize] as char);
            j += 1;
        }
        let _ = writeln!(&mut *f);
    }
    let _ = writeln!(&mut *f);
}

// =============================================================================
// src_flatten/sequence.c lines 2000-2069
// void build_trna(gene *t, char matrix[][MATY], int x, int y, csw *sw)
// =============================================================================
pub unsafe fn build_trna(t: *mut gene, matrix: *mut [u8; MATY], mut x: i32, mut y: i32, sw: *mut csw) {
    let mut i: i32;
    let mut j: i32;
    let mut e: i32;
    let mut c: i32;
    let mut seq: *mut i32;
    let mut rseq: [i32; 150] = [0; 150];
    static BOND2: [u8; 4] = *b" +-.";

    (*t).varbp = find_var_hairpin(t);
    remove_inserts((*t).seq.as_mut_ptr(), rseq.as_mut_ptr());
    seq = rseq.as_mut_ptr();

    i = 0;
    while i < (*t).astem1 {
        (*matrix.add(x as usize))[y as usize] = cbase(*seq);
        seq = seq.add(1);
        y -= 1;
        i += 1;
    }
    if (*t).spacer1 > 0 {
        x -= 1;
        if (*t).spacer1 >= 3 {
            (*matrix.add(x as usize))[(y + 1) as usize] = cbase(*seq);
            seq = seq.add(1);
        }
        (*matrix.add(x as usize))[y as usize] = cbase(*seq);
        seq = seq.add(1);
        y -= 1;
        x -= 1;
        if (*t).spacer1 >= 2 {
            (*matrix.add(x as usize))[y as usize] = cbase(*seq);
            seq = seq.add(1);
        }
        if (*t).spacer2 < 2 || (*t).spacer1 > 1 {
            x -= 1;
            y -= 1;
        }
    }
    if (*t).dstem > 0 {
        e = 2 * (*t).dstem + (*t).dloop;
        seq = make_clover(seq, 0, e, (*t).dstem, matrix, &mut x, &mut y, LEFT);
        if (*t).spacer2 > 1 {
            x -= 1;
        }
        y -= 1;
        if (*t).spacer2 > 0 {
            (*matrix.add(x as usize))[y as usize] = cbase(*seq);
            seq = seq.add(1);
        }
        y -= 1;
        if (*t).spacer2 > 1 {
            if (*t).spacer1 > 1 {
                x += 1;
            }
            (*matrix.add(x as usize))[y as usize] = cbase(*seq);
            seq = seq.add(1);
            if (*t).spacer1 < 2 {
                y -= 1;
            }
        }
        x += 1;
    } else {
        seq = make_dv(seq, matrix, (*t).dloop, RIGHT, &mut x, &mut y);
    }
    e = 2 * (*t).cstem + (*t).cloop;
    seq = make_clover(seq, 0, e, (*t).cstem, matrix, &mut x, &mut y, DOWN);
    if (*t).tstem > 0 {
        seq = make_var(seq, matrix, &mut x, &mut y, RIGHT, (*t).var, (*t).varbp);
        e = 2 * (*t).tstem + (*t).tloop;
        seq = make_clover(seq, 0, e, (*t).tstem, matrix, &mut x, &mut y, RIGHT);
        y += 1;
    } else {
        seq = make_tv(seq, matrix, &mut x, &mut y, RIGHT, (*t).tloop);
    }
    e = (*t).astem2;
    i = 0;
    while i < e {
        c = *seq;
        seq = seq.add(1);
        if c < Adenine {
            break;
        }
        (*matrix.add(x as usize))[y as usize] = cbase(c);
        j = base_match(
            (*matrix.add(x as usize))[y as usize],
            (*matrix.add((x - 2) as usize))[y as usize]
        );
        (*matrix.add((x - 1) as usize))[y as usize] = BOND2[j as usize];
        y += 1;
        i += 1;
    }
    i = 0;
    e = if (*sw).aataildisp != 0 { ASTEM2_EXTD } else { (*t).aatail };
    j = if e < 2 { e } else { 2 };
    while i < j {
        c = *seq;
        seq = seq.add(1);
        if c < Adenine {
            break;
        }
        (*matrix.add(x as usize))[y as usize] = cbase(c);
        x += 1;
        y += 1;
        i += 1;
    }
    e -= j;
    i = 0;
    while i < e {
        c = *seq;
        seq = seq.add(1);
        if c < Adenine {
            break;
        }
        (*matrix.add(x as usize))[y as usize] = cbase(c);
        x += 1;
        i += 1;
    }
}

// =============================================================================
// src_flatten/sequence.c lines 2075-2119
// void build_tmrna(gene *t, char matrix[][MATY], int x, int y, csw *sw)
// =============================================================================
pub unsafe fn build_tmrna(t: *mut gene, matrix: *mut [u8; MATY], mut x: i32, mut y: i32, sw: *mut csw) {
    let mut i: i32;
    let mut j: i32;
    let mut e: i32;
    let mut c: i32;
    let tarm: i32;
    let mut seq: *mut i32;
    let mut rseq: [i32; 2 * MAXTMRNALEN as usize + 1] = [0; 2 * MAXTMRNALEN as usize + 1];
    static BOND2: [u8; 4] = *b" +-.";

    remove_inserts((*t).eseq.as_mut_ptr(), rseq.as_mut_ptr());
    seq = rseq.as_mut_ptr().add((*t).asst as usize);

    i = 0;
    while i < (*t).astem1 {
        (*matrix.add(x as usize))[y as usize] = cbase(*seq);
        seq = seq.add(1);
        y -= 1;
        i += 1;
    }
    seq = make_dv(seq, matrix, (*t).dloop, RIGHT, &mut x, &mut y);
    tarm = 2 * (*t).tstem + (*t).tloop;
    e = if (*t).asst > 0 {
        (*t).cstem - (*t).dloop - (*t).astem1 - (*t).asst + 54
    } else {
        2 * (*t).cstem + (*t).cloop + (*t).nintron
    };
    seq = make_clover(seq, 0, e, (*t).cstem, matrix, &mut x, &mut y, DOWN);
    seq = make_var(seq, matrix, &mut x, &mut y, RIGHT, (*t).var, (*t).varbp);
    seq = make_clover(seq, 0, tarm, (*t).tstem, matrix, &mut x, &mut y, RIGHT);
    y += 1;
    e = (*t).astem2;
    i = 0;
    while i < e {
        c = *seq;
        seq = seq.add(1);
        if c == TERM {
            break;
        }
        (*matrix.add(x as usize))[y as usize] = cbase(c);
        j = base_match(
            (*matrix.add(x as usize))[y as usize],
            (*matrix.add((x - 2) as usize))[y as usize]
        );
        (*matrix.add((x - 1) as usize))[y as usize] = BOND2[j as usize];
        y += 1;
        i += 1;
    }
    e = if (*sw).aataildisp != 0 { ASTEM2_EXTD } else { (*t).aatail };
    j = if e < 2 { e } else { 2 };
    i = 0;
    while i < j {
        c = *seq;
        seq = seq.add(1);
        if c == TERM {
            break;
        }
        (*matrix.add(x as usize))[y as usize] = cbase(c);
        x += 1;
        y += 1;
        i += 1;
    }
    e -= j;
    i = 0;
    while i < e {
        c = *seq;
        seq = seq.add(1);
        if c == TERM {
            break;
        }
        (*matrix.add(x as usize))[y as usize] = cbase(c);
        x += 1;
        i += 1;
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1936-1984
// void disp_peptide_tag(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_peptide_tag(f: *mut File, t: *mut gene, sw: *mut csw) {
    let mut i: i32;
    let lx: i32;
    let c1: i32;
    let c2: i32;
    let c3: i32;
    let mut s: *mut i32;
    let se: *mut i32;
    let mut tag: [u8; 52] = [0; 52];

    let _ = write!(&mut *f, "Tag peptide at [{},{}]\nTag sequence: ", (*t).tps + 1, (*t).tpe + 1);
    lx = peptide_tag(tag.as_mut_ptr(), 50, t, sw);
    se = (*t).eseq.as_mut_ptr().add((*t).tps as usize);
    s = se;
    i = 0;
    while i < lx {
        if i > 0 {
            let _ = write!(&mut *f, "-");
        }
        let c1 = *s;
        s = s.add(1);
        if c1 >= AMBIG {
            i += 1;
            continue;
        }
        let c2 = *s;
        s = s.add(1);
        if c2 >= AMBIG {
            i += 1;
            continue;
        }
        let c3 = *s;
        s = s.add(1);
        if c3 >= AMBIG {
            i += 1;
            continue;
        }
        let _ = write!(&mut *f, "{}{}{}", cbase(c1) as char, cbase(c2) as char, cbase(c3) as char);
        i += 1;
    }
    s = se;
    let _ = write!(&mut *f, "\nTag peptide:  ");
    i = 0;
    while i < lx {
        let _ = write!(&mut *f, "{}", CStr::from_ptr(translate(s, sw) as *const i8).to_string_lossy());
        s = s.add(3);
        if i < lx - 1 {
            let _ = write!(&mut *f, "-");
        }
        i += 1;
    }
    let _ = write!(&mut *f, "\nTag peptide:  {}", CStr::from_ptr(tag.as_ptr() as *const i8).to_string_lossy());
    if (*sw).energydisp != 0 {
        s = se;
        let _ = write!(&mut *f, "\nTag Polarity: ");
        i = 0;
        while i < lx {
            let _ = write!(&mut *f, "{}", ptranslate(s, sw) as u8 as char);
            s = s.add(3);
            i += 1;
        }
    }
    let _ = writeln!(&mut *f);
    // Identify tag
    let mut thit: [[u8; 50]; 21] = [[0; 50]; 21];
    let nmh = identify_tag(tag.as_ptr(), lx, thit.as_mut_ptr(), 21);
    if nmh > 0 {
        if nmh > 1 {
            let _ = writeln!(&mut *f, "Match with tmRNA tags from:");
            let mut line_count = 0;
            for nm in 0..nmh {
                line_count += 1;
                if line_count > 3 {
                    let _ = writeln!(&mut *f);
                    line_count = 1;
                }
                if line_count > 1 {
                    let _ = write!(&mut *f, ", ");
                }
                let _ = write!(&mut *f, "{}", CStr::from_ptr(thit[nm as usize].as_ptr() as *const i8).to_string_lossy());
            }
            let _ = writeln!(&mut *f);
        } else {
            let _ = writeln!(&mut *f, "Match with {} tmRNA tag", CStr::from_ptr(thit[0].as_ptr() as *const i8).to_string_lossy());
        }
    } else if nmh == -1 {
        let _ = writeln!(&mut *f, "Match with many tmRNA tags");
    } else {
        let _ = writeln!(&mut *f, "Tag not identified");
    }
    let _ = writeln!(&mut *f);
}

// =============================================================================
// src_flatten/sequence.c lines 958-1044
// void disp_trna_bracket_notation(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_trna_bracket_notation(f: *mut File, t: *mut gene, sw: *mut csw) {
    let mut j: i32;
    let mut k: i32;
    let mut varbp: i32;
    let mut stem: i32;
    let ab: i32;
    let ae: i32;
    let hl: i32;
    let mut s: *mut i32;
    let se: *mut i32;
    let mut sl: *mut i32;
    let mut sb: *mut i32;
    let mut sr: *mut i32;
    let mut genename: [u8; 100] = [0; 100];
    static BPLB: [u8; 2] = [b'.', b'('];
    static BPRB: [u8; 2] = [b'.', b')'];

    if (*sw).batch == 0 {
        name(t, genename.as_mut_ptr(), 1, sw);
        let _ = writeln!(&mut *f, "\nSecondary structure (bracket notation) for {}", CStr::from_ptr(genename.as_ptr() as *const i8).to_string_lossy());
    }
    if (*t).nintron > 0 {
        s = (*t).eseq.as_mut_ptr();
        se = s.add(((*t).nbase + (*t).nintron) as usize);
    } else {
        s = (*t).seq.as_mut_ptr();
        se = s.add((*t).nbase as usize);
    }
    sl = s;
    while sl < se {
        let _ = write!(&mut *f, "{}", cbase(*sl) as char);
        sl = sl.add(1);
    }
    let _ = writeln!(&mut *f);
    sl = s;
    sr = se.sub(((*t).aatail + 1) as usize);
    for _i in 0..(*t).astem1 {
        let _ = write!(&mut *f, "{}", BPLB[BP[*sl as usize][*sr as usize] as usize] as char);
        sl = sl.add(1);
        sr = sr.sub(1);
    }
    for _i in 0..(*t).spacer1 {
        let _ = write!(&mut *f, "s");
    }
    sl = sl.add((*t).spacer1 as usize);
    sb = sl.add(((*t).dstem - 1) as usize);
    sr = sb.add(((*t).dstem + (*t).dloop) as usize);
    for _i in 0..(*t).dstem {
        let _ = write!(&mut *f, "{}", BPLB[BP[*sl as usize][*sr as usize] as usize] as char);
        sl = sl.add(1);
        sr = sr.sub(1);
    }
    for _i in 0..(*t).dloop {
        let _ = write!(&mut *f, "d");
    }
    sl = sl.add((*t).dloop as usize);
    for _i in 0..(*t).dstem {
        let _ = write!(&mut *f, "{}", BPRB[BP[*sl as usize][*sb as usize] as usize] as char);
        sl = sl.add(1);
        sb = sb.sub(1);
    }
    for _i in 0..(*t).spacer2 {
        let _ = write!(&mut *f, "s");
    }
    sl = sl.add((*t).spacer2 as usize);
    sb = sl.add(((*t).cstem - 1) as usize);
    sr = sb.add(((*t).cstem + (*t).cloop + (*t).nintron) as usize);
    for _i in 0..(*t).cstem {
        let _ = write!(&mut *f, "{}", BPLB[BP[*sl as usize][*sr as usize] as usize] as char);
        sl = sl.add(1);
        sr = sr.sub(1);
    }
    hl = (*t).astem1 + (*t).spacer1 + 2 * (*t).dstem + (*t).dloop + (*t).spacer2 + (*t).cstem;
    if (*t).nintron > 0 {
        j = (*t).intron - hl;
        ab = (*t).anticodon - hl;
        ae = ab + (*t).cloop - 5;
        for ii in 0..j {
            if ii <= ae && ii >= ab {
                let _ = write!(&mut *f, "A");
            } else {
                let _ = write!(&mut *f, "c");
            }
        }
        for _ii in 0..(*t).nintron {
            let _ = write!(&mut *f, "i");
        }
        for ii in j..(*t).cloop {
            if ii <= ae && ii >= ab {
                let _ = write!(&mut *f, "A");
            } else {
                let _ = write!(&mut *f, "c");
            }
        }
    } else {
        j = (*t).cloop - 4;
        ab = (*t).anticodon - hl;
        ae = (*t).cloop - ab - j;
        for _ii in 0..ab {
            let _ = write!(&mut *f, "c");
        }
        for _ii in 0..j {
            let _ = write!(&mut *f, "A");
        }
        for _ii in 0..ae {
            let _ = write!(&mut *f, "c");
        }
    }
    sl = sl.add(((*t).cloop + (*t).nintron) as usize);
    for _i in 0..(*t).cstem {
        let _ = write!(&mut *f, "{}", BPRB[BP[*sl as usize][*sb as usize] as usize] as char);
        sl = sl.add(1);
        sb = sb.sub(1);
    }
    varbp = find_var_hairpin(t);
    if varbp > 0 {
        j = varbp >> 10;
        k = (varbp >> 5) & 0x1f;
        stem = varbp & 0x1f;
        sr = sl.add((k + stem - 1) as usize);
        sl = sl.add(j as usize);
        sb = sl.add((stem - 1) as usize);
        for _ii in 0..j {
            let _ = write!(&mut *f, "s");
        }
        for _ii in 0..stem {
            let _ = write!(&mut *f, "{}", BPLB[BP[*sl as usize][*sr as usize] as usize] as char);
            sl = sl.add(1);
            sr = sr.sub(1);
        }
        for _ii in (j + stem)..k {
            let _ = write!(&mut *f, "v");
            sl = sl.add(1);
        }
        for _ii in 0..stem {
            let _ = write!(&mut *f, "{}", BPRB[BP[*sl as usize][*sb as usize] as usize] as char);
            sl = sl.add(1);
            sb = sb.sub(1);
        }
        for _ii in (k + stem)..(*t).var {
            let _ = write!(&mut *f, "s");
            sl = sl.add(1);
        }
    } else {
        for _ii in 0..(*t).var {
            let _ = write!(&mut *f, "v");
        }
        sl = sl.add((*t).var as usize);
    }
    sb = sl.add(((*t).tstem - 1) as usize);
    sr = sb.add(((*t).tstem + (*t).tloop) as usize);
    for _i in 0..(*t).tstem {
        let _ = write!(&mut *f, "{}", BPLB[BP[*sl as usize][*sr as usize] as usize] as char);
        sl = sl.add(1);
        sr = sr.sub(1);
    }
    for _i in 0..(*t).tloop {
        let _ = write!(&mut *f, "t");
    }
    sl = sl.add((*t).tloop as usize);
    for _i in 0..(*t).tstem {
        let _ = write!(&mut *f, "{}", BPRB[BP[*sl as usize][*sb as usize] as usize] as char);
        sl = sl.add(1);
        sb = sb.sub(1);
    }
    sb = s.add(((*t).astem1 - 1) as usize);
    for _i in 0..(*t).astem2 {
        let _ = write!(&mut *f, "{}", BPRB[BP[*sl as usize][*sb as usize] as usize] as char);
        sl = sl.add(1);
        sb = sb.sub(1);
    }
    let _ = writeln!(&mut *f);
    if (*sw).batch == 0 {
        let _ = writeln!(&mut *f);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1048-1120
// void disp_tmrna_trnadomain_bracket_notation(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_tmrna_trnadomain_bracket_notation(f: *mut File, t: *mut gene, sw: *mut csw) {
    let mut j: i32;
    let mut k: i32;
    let mut varbp: i32;
    let mut stem: i32;
    let s: *mut i32;
    let sa: *mut i32;
    let mut sb: *mut i32;
    let sc: *mut i32;
    let sd: *mut i32;
    let se: *mut i32;
    let sf: *mut i32;
    let mut sl: *mut i32;
    let mut sr: *mut i32;
    static BPLB: [u8; 2] = [b'.', b'('];
    static BPRB: [u8; 2] = [b'.', b')'];

    if (*t).nintron <= 0 {
        return;
    }
    if (*sw).batch == 0 {
        let _ = writeln!(&mut *f, "Secondary structure (bracket notation) for tRNA domain:");
    }
    s = (*t).eseq.as_mut_ptr();
    se = s.add(((*t).nbase + (*t).nintron) as usize);
    if (*t).asst > 0 {
        sc = s.add((*t).asst as usize);
        sa = sc.add(((*t).astem1 + (*t).dloop + (*t).cstem) as usize);
        sd = s.add(54);
        sf = s.add((*t).intron as usize);
    } else {
        sc = s;
        sa = s.add(((*t).astem1 + (*t).dloop + (*t).cstem) as usize);
        sd = se.sub(((*t).aatail + (*t).astem2 + 2 * (*t).tstem + (*t).tloop + (*t).var + (*t).cstem) as usize);
        sf = se;
    }
    sl = sc;
    while sl < sa {
        let _ = write!(&mut *f, "{}", cbase(*sl) as char);
        sl = sl.add(1);
    }
    let _ = write!(&mut *f, "|");
    sr = sd;
    while sr < sf {
        let _ = write!(&mut *f, "{}", cbase(*sr) as char);
        sr = sr.add(1);
    }
    let _ = writeln!(&mut *f);
    sl = sc;
    sr = sf.sub(((*t).aatail + 1) as usize);
    for _i in 0..(*t).astem1 {
        let _ = write!(&mut *f, "{}", BPLB[BP[*sl as usize][*sr as usize] as usize] as char);
        sl = sl.add(1);
        sr = sr.sub(1);
    }
    for _i in 0..(*t).spacer1 {
        let _ = write!(&mut *f, "s");
    }
    sl = sl.add((*t).spacer1 as usize);
    sb = sl.add(((*t).dstem - 1) as usize);
    sr = sb.add(((*t).dstem + (*t).dloop) as usize);
    for _i in 0..(*t).dstem {
        let _ = write!(&mut *f, "{}", BPLB[BP[*sl as usize][*sr as usize] as usize] as char);
        sl = sl.add(1);
        sr = sr.sub(1);
    }
    for _i in 0..(*t).dloop {
        let _ = write!(&mut *f, "d");
    }
    sl = sl.add((*t).dloop as usize);
    for _i in 0..(*t).dstem {
        let _ = write!(&mut *f, "{}", BPRB[BP[*sl as usize][*sb as usize] as usize] as char);
        sl = sl.add(1);
        sb = sb.sub(1);
    }
    for _i in 0..(*t).spacer2 {
        let _ = write!(&mut *f, "s");
    }
    sl = sl.add((*t).spacer2 as usize);
    sb = sl.add(((*t).cstem - 1) as usize);
    sr = sd.add(((*t).cstem - 1) as usize);
    for _i in 0..(*t).cstem {
        let _ = write!(&mut *f, "{}", BPLB[BP[*sl as usize][*sr as usize] as usize] as char);
        sl = sl.add(1);
        sr = sr.sub(1);
    }
    let _ = write!(&mut *f, "|");
    sl = sd;
    for _i in 0..(*t).cstem {
        let _ = write!(&mut *f, "{}", BPRB[BP[*sl as usize][*sb as usize] as usize] as char);
        sl = sl.add(1);
        sb = sb.sub(1);
    }
    varbp = find_var_hairpin(t);
    if varbp > 0 {
        j = varbp >> 10;
        k = (varbp >> 5) & 0x1f;
        stem = varbp & 0x1f;
        sr = sl.add((k + stem - 1) as usize);
        sl = sl.add(j as usize);
        sb = sl.add((stem - 1) as usize);
        for _ii in 0..j {
            let _ = write!(&mut *f, " ");
        }
        for _ii in 0..stem {
            let _ = write!(&mut *f, "{}", BPLB[BP[*sl as usize][*sr as usize] as usize] as char);
            sl = sl.add(1);
            sr = sr.sub(1);
        }
        for _ii in (j + stem)..k {
            let _ = write!(&mut *f, "v");
            sl = sl.add(1);
        }
        for _ii in 0..stem {
            let _ = write!(&mut *f, "{}", BPRB[BP[*sl as usize][*sb as usize] as usize] as char);
            sl = sl.add(1);
            sb = sb.sub(1);
        }
        for _ii in (k + stem)..(*t).var {
            let _ = write!(&mut *f, " ");
            sl = sl.add(1);
        }
    } else {
        for _ii in 0..(*t).var {
            let _ = write!(&mut *f, "v");
        }
        sl = sl.add((*t).var as usize);
    }
    sb = sl.add(((*t).tstem - 1) as usize);
    sr = sb.add(((*t).tstem + (*t).tloop) as usize);
    for _i in 0..(*t).tstem {
        let _ = write!(&mut *f, "{}", BPLB[BP[*sl as usize][*sr as usize] as usize] as char);
        sl = sl.add(1);
        sr = sr.sub(1);
    }
    for _i in 0..(*t).tloop {
        let _ = write!(&mut *f, "t");
    }
    sl = sl.add((*t).tloop as usize);
    for _i in 0..(*t).tstem {
        let _ = write!(&mut *f, "{}", BPRB[BP[*sl as usize][*sb as usize] as usize] as char);
        sl = sl.add(1);
        sb = sb.sub(1);
    }
    sb = sc.add(((*t).astem1 - 1) as usize);
    for _i in 0..(*t).astem2 {
        let _ = write!(&mut *f, "{}", BPRB[BP[*sl as usize][*sb as usize] as usize] as char);
        sl = sl.add(1);
        sb = sb.sub(1);
    }
    for _i in 0..(*t).aatail {
        let _ = write!(&mut *f, "a");
    }
    let _ = writeln!(&mut *f);
    if (*sw).batch == 0 {
        let _ = writeln!(&mut *f);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 888-955
// void disp_gene_SVG(gene *t, char m[][MATY], csw *sw)
// =============================================================================
pub unsafe fn disp_gene_SVG(t: *mut gene, m: *mut [u8; MATY], sw: *mut csw) {
    let mut i: i32;
    let mut xb: i32;
    let mut xe: i32;
    let mut yb: i32;
    let mut ye: i32;
    let mut xbm: i32;
    let mut xem: i32;
    let mut ybm: i32;
    let mut yem: i32;
    let mut xdiff: i32;
    let mut ydiff: i32;
    let mut xpos: f64;
    let mut ypos: f64;
    let xsc: f64;
    let ysc: f64;
    let fontsize: f64;
    let yv: f64;
    let mut genename: [u8; 100] = [0; 100];
    let f = (*sw).f;

    xe = 0;
    xb = MATX as i32;
    ye = 0;
    yb = MATY as i32;
    for yy in 5..=30 {
        for xx in 0..MATX as i32 {
            if (*m.add(xx as usize))[yy as usize] > b' ' && (*m.add(xx as usize))[yy as usize] <= b'~' {
                if xx < xb { xb = xx; }
                if xx > xe { xe = xx; }
                if yy < yb { yb = yy; }
                if yy > ye { ye = yy; }
            }
        }
    }
    xbm = xb - 5;
    if xbm < 0 { xbm = 0; }
    xem = xe + 5;
    if xem >= MATX as i32 { xem = MATX as i32 - 1; }
    if (xb - xbm) > (xem - xe) { xbm = xb - (xem - xe); }
    else { xem = xe + (xb - xbm); }
    xdiff = xem - xbm;
    ybm = yb - 5;
    if ybm < 5 { ybm = 5; }
    yem = ye + 5;
    if yem > 30 { yem = 30; }
    if (yb - ybm) > (yem - ye) { ybm = yb - (yem - ye); }
    else { yem = ye + (yb - ybm); }
    xdiff = xem - xbm;
    ydiff = yem - ybm + 4;
    ysc = 10.0;
    xsc = 0.1 * ((10.0 * ysc * (xdiff as f64 / ydiff as f64) + 0.5) as i32) as f64;
    fontsize = 1.4;
    yv = 0.01 * ((18.5 * ydiff as f64 * fontsize / ysc + 0.5) as i32) as f64;
    name(t, genename.as_mut_ptr(), 1, sw);
    if (*sw).batch == 0 {
        let _ = writeln!(&mut *f, "Scalable vector graphics (SVG) image:");
    }
    let _ = write!(&mut *f, "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' ");
    let _ = writeln!(&mut *f, "width='{}cm' height='{}cm' viewBox='0 0 {} {}'>", xsc, ysc, xdiff, ydiff);
    let _ = writeln!(&mut *f, "<title>{}</title>", CStr::from_ptr(genename.as_ptr() as *const i8).to_string_lossy());
    let _ = write!(&mut *f, "<g font-family='Courier New,Courier,monospace' font-size='{}' ", fontsize);
    let _ = writeln!(&mut *f, "text-anchor='middle' fill='black' stroke='none'>");
    i = 0;
    for yy in ybm..=yem {
        ypos = (ydiff - 2 - (yy - ybm)) as f64;
        for xx in xbm..=xem {
            if (*m.add(xx as usize))[yy as usize] > b' ' && (*m.add(xx as usize))[yy as usize] <= b'~' && (*m.add(xx as usize))[yy as usize] != b'!' {
                xpos = (xx - xbm) as f64;
                let _ = write!(&mut *f, "<text x='{}' y='{}'>{}</text>", xpos, ypos, (*m.add(xx as usize))[yy as usize] as char);
                i += 1;
                if i >= 4 {
                    let _ = writeln!(&mut *f);
                    i = 0;
                }
            }
        }
    }
    if i > 0 {
        let _ = writeln!(&mut *f);
    }
    let _ = writeln!(&mut *f, "</g><g fill='none' stroke='black' stroke-width='0.075'>");
    i = 0;
    for yy in ybm..=yem {
        ypos = (ydiff - 2 - (yy - ybm)) as f64;
        for xx in xbm..=xem {
            if (*m.add(xx as usize))[yy as usize] == b'!' {
                xpos = (xx - xbm) as f64;
                let _ = write!(&mut *f, "<line x1='{}' y1='{}' x2='{}' y2='{}'/>", xpos, ypos, xpos, ypos - yv);
                i += 1;
                if i >= 2 {
                    let _ = writeln!(&mut *f);
                    i = 0;
                }
            }
        }
    }
    if i > 0 {
        let _ = writeln!(&mut *f);
    }
    let _ = writeln!(&mut *f, "</g></svg>");
}

// =============================================================================
// src_flatten/sequence.c lines 248-251
// void init_tmrna(FILE *f, csw *sw)
// =============================================================================
pub unsafe fn init_tmrna(f: *mut File, sw: *mut csw) {
    let mut c: i32;
    let mut s: *mut i32;
    s = (*sw).tmrna_struct.as_mut_ptr();
    loop {
        c = *s;
        s = s.add(1);
        if c == TERM {
            break;
        }
        // itmparam(cbase(c), f) - would be called here
        let _ = write!(&mut *f, "{}", cbase(c) as char);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1497-1536
// double find_tag_upstream_hairpin(int *se)
// =============================================================================
pub unsafe fn find_tag_upstream_hairpin(mut se: *mut i32) -> f64 {
    let sb: *mut i32;
    let mut sd: *mut i32;
    let sf: *mut i32;
    let mut sh: *mut i32;
    let mut s: *mut i32;
    let mut found: i32;
    let mut c: u32;
    let mut m: u32;
    let mut mx: u32;
    static A: [u32; 6] = [0, 0, 0, 0x10000, 0, 0];
    static C: [u32; 6] = [0, 0, 0x10000, 0, 0, 0];
    static G: [u32; 6] = [0, 0x10000, 0, 0x10000, 0, 0];
    static T: [u32; 6] = [0x10000, 0, 0x10000, 0, 0, 0];
    let mut t: [u32; 6] = [0, 0, 0, 0, 0, 0];

    mx = 0;
    found = 0;
    sf = se.sub(4);
    sb = se.sub(20);
    t[0] = A[*se as usize];
    t[1] = C[*se as usize];
    t[2] = G[*se as usize];
    t[3] = T[*se as usize];
    se = se.sub(1);
    while se > sf {
        t[0] = (t[0] >> 4) | A[*se as usize];
        t[1] = (t[1] >> 4) | C[*se as usize];
        t[2] = (t[2] >> 4) | G[*se as usize];
        t[3] = (t[3] >> 4) | T[*se as usize];
        se = se.sub(1);
    }
    sh = se.sub(4);
    sd = se.sub(30);
    while se > sb && found == 0 {
        t[0] = (t[0] >> 4) | A[*se as usize];
        t[1] = (t[1] >> 4) | C[*se as usize];
        t[2] = (t[2] >> 4) | G[*se as usize];
        t[3] = (t[3] >> 4) | T[*se as usize];
        s = sh;
        m = t[*s as usize];
        s = s.sub(1);
        while s > sd {
            m = (m >> 4) + t[*s as usize];
            c = m & 0xf;
            if c > mx { mx = c; }
            if mx == 5 {
                found = 1;
                break;
            }
            s = s.sub(1);
        }
        if found == 0 {
            sd = sd.sub(1);
            sh = sh.sub(1);
            se = se.sub(1);
        }
    }
    if found != 0 { 15.0 } else { 0.0 }
}

// =============================================================================
// src_flatten/sequence.c lines 1540-1579
// double find_taghairpin(int *seq)
// =============================================================================
pub unsafe fn find_taghairpin(seq: *mut i32) -> f64 {
    let mut i: i32;
    let mut s: *mut i32;
    let mut sb: *mut i32;
    let mut se: *mut i32;
    let sf: *mut i32;
    let mut c: u32;
    let mut m: u32;
    let mut mx: u32;
    static A: [u32; 6] = [0, 0, 0, 1, 0, 0];
    static C: [u32; 6] = [0, 0, 1, 0, 0, 0];
    static G: [u32; 6] = [0, 1, 0, 1, 0, 0];
    static T: [u32; 6] = [1, 0, 1, 0, 0, 0];
    let mut t: [u32; 6] = [0, 0, 0, 0, 0, 0];

    mx = 0;
    sb = seq.sub(20);
    se = seq.sub(13);
    sf = seq.sub(4);
    t[0] = A[*sb as usize];
    t[1] = C[*sb as usize];
    t[2] = G[*sb as usize];
    t[3] = T[*sb as usize];
    sb = sb.add(1);
    while sb < se {
        t[0] = (t[0] << 4) | A[*sb as usize];
        t[1] = (t[1] << 4) | C[*sb as usize];
        t[2] = (t[2] << 4) | G[*sb as usize];
        t[3] = (t[3] << 4) | T[*sb as usize];
        sb = sb.add(1);
    }
    while sb < sf {
        t[0] = ((t[0] << 4) | A[*sb as usize]) & 0xffffffff;
        t[1] = ((t[1] << 4) | C[*sb as usize]) & 0xffffffff;
        t[2] = ((t[2] << 4) | G[*sb as usize]) & 0xffffffff;
        t[3] = ((t[3] << 4) | T[*sb as usize]) & 0xffffffff;
        sb = sb.add(1);
        s = seq.add(20);
        se = seq.add(2);
        m = t[*s as usize];
        s = s.sub(1);
        while s > se {
            m = (m >> 4) + t[*s as usize];
            c = m & 0xf;
            if c > mx { mx = c; }
            s = s.sub(1);
        }
        i = 7 - mx as i32;
        while { let r = i > 0; i -= 1; r } {
            m = m >> 4;
            c = m & 0xf;
            if c > mx { mx = c; }
        }
    }
    (mx << 1) as f64
}

// =============================================================================
// src_flatten/sequence.c lines 723-726
// void disp_location(gene *t, csw *sw, char *m)
// =============================================================================
pub unsafe fn disp_location(t: *mut gene, sw: *mut csw, m: *const u8) {
    let f = (*sw).f;
    let mut sp: [u8; 80] = [0; 80];
    position(sp.as_mut_ptr(), t, sw);
    let _ = writeln!(&mut *f, "{} {}", CStr::from_ptr(m as *const i8).to_string_lossy(), CStr::from_ptr(sp.as_ptr() as *const i8).to_string_lossy());
}

// =============================================================================
// src_flatten/sequence.c lines 800-827
// void disp_intron(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_intron(f: *mut File, t: *mut gene, sw: *mut csw) {
    let mut i: i32;
    let mut s: *mut i32;
    let sb: *mut i32;
    let se: *mut i32;
    let mut genename: [u8; 100] = [0; 100];

    if (*t).nintron <= 0 {
        return;
    }
    name(t, genename.as_mut_ptr(), 1, sw);
    let _ = writeln!(&mut *f, "Intron from {}", CStr::from_ptr(genename.as_ptr() as *const i8).to_string_lossy());
    let _ = writeln!(&mut *f, "1   .   10    .   20    .   30    .   40    .   50");
    sb = (*t).eseq.as_mut_ptr().add((*t).intron as usize);
    s = sb;
    se = sb.add((*t).nintron as usize);
    i = 0;
    while s < se {
        let c = *s;
        s = s.add(1);
        if c < Adenine {
            break;
        }
        let _ = write!(&mut *f, "{}", cbase(c) as char);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    if i > 0 {
        let _ = writeln!(&mut *f);
    }
    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f, "Intron Length: {}", (*t).nintron);
    let _ = write!(&mut *f, "Intron Insertion Position({}-{}): ", (*t).intron, (*t).intron + 1);
    s = sb.sub(5);
    for _ii in 0..5 {
        let _ = write!(&mut *f, "{}", cbase(*s) as char);
        s = s.add(1);
    }
    let _ = write!(&mut *f, "-Intron-");
    s = se;
    for _ii in 0..5 {
        let _ = write!(&mut *f, "{}", cbase(*s) as char);
        s = s.add(1);
    }
    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f);
}

// =============================================================================
// src_flatten/sequence.c lines 830-854
// void disp_fasta_seq(FILE *f, gene *t, int ns, int n, int nsp, int c, csw *sw)
// =============================================================================
pub unsafe fn disp_fasta_seq(f: *mut File, t: *mut gene, ns: i32, n: i32, nsp: i32, c: i32, sw: *mut csw) {
    let mut i: i32;
    let mut s: *mut i32;
    let se: *mut i32;
    let mut genename: [u8; 100] = [0; 100];
    let mut genepos: [u8; 100] = [0; 100];

    if (*t).nintron > 0 {
        s = (*t).eseq.as_mut_ptr();
        se = s.add(((*t).nbase + (*t).nintron) as usize);
    } else {
        s = (*t).seq.as_mut_ptr();
        se = s.add((*t).nbase as usize);
    }
    name(t, genename.as_mut_ptr(), 1, sw);
    position(genepos.as_mut_ptr(), t, sw);
    let genename_str = CStr::from_ptr(genename.as_ptr() as *const i8).to_string_lossy();
    let genepos_str = CStr::from_ptr(genepos.as_ptr() as *const i8).to_string_lossy();
    if nsp > 0 {
        if ns > 0 {
            let _ = writeln!(&mut *f, ">{}-{}{}{}", ns, n, genename_str, genepos_str);
        } else {
            let _ = writeln!(&mut *f, ">{}{}", genename_str, genepos_str);
        }
    } else {
        if ns > 0 {
            let _ = writeln!(&mut *f, ">{}-{} {} {}", ns, n, genename_str, genepos_str);
        } else {
            let _ = writeln!(&mut *f, ">{} {}", genename_str, genepos_str);
        }
    }
    i = 0;
    while s < se {
        if c != 0 {
            let _ = write!(&mut *f, "{}", cpbase(*s) as char);
        } else {
            let _ = write!(&mut *f, "{}", cbase(*s) as char);
        }
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    if i > 0 {
        let _ = writeln!(&mut *f);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1125-1176
// void disp_tmrna_seq(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_tmrna_seq(f: *mut File, t: *mut gene, sw: *mut csw) {
    let mut i: i32;
    let mut s: *mut i32;
    let sb: *mut i32;
    let mut se: *mut i32;

    if (*t).nintron <= 0 {
        return;
    }
    if (*t).name[0] == 0 {
        let _ = writeln!(&mut *f, "tmRNA sequence\n");
    } else {
        let _ = writeln!(&mut *f, "tmRNA sequence in {}\n", CStr::from_ptr((*t).name.as_ptr() as *const i8).to_string_lossy());
    }
    let _ = writeln!(&mut *f, "1   .   10    .   20    .   30    .   40    .   50");
    sb = (*t).eseq.as_mut_ptr();
    s = sb;
    se = sb.add((*t).intron as usize);
    i = 0;
    while s < se {
        let _ = write!(&mut *f, "{}", cbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    se = sb.add((*t).tps as usize);
    while s < se {
        let _ = write!(&mut *f, "{}", cpbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    se = sb.add(((*t).tpe + 1) as usize);
    while ltranslate(se, t, sw) == b'*' {
        se = se.add(3);
    }
    while s < se {
        let _ = write!(&mut *f, "{}", cbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    se = sb.add(((*t).intron + (*t).nintron) as usize);
    while s < se {
        let _ = write!(&mut *f, "{}", cpbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    se = sb.add(((*t).nbase + (*t).nintron) as usize);
    while s < se {
        let _ = write!(&mut *f, "{}", cbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    if i > 0 {
        let _ = writeln!(&mut *f);
    }
    let _ = writeln!(&mut *f, "\n5' tRNA domain at [1,{}]", (*t).intron);
    let _ = writeln!(&mut *f, "3' tRNA domain at [{},{}]",
        (*t).intron + (*t).nintron + 1, (*t).nbase + (*t).nintron);
    if ((*sw).secstructdisp & 2) != 0 {
        disp_tmrna_trnadomain_bracket_notation(f, t, sw);
    }
    let _ = write!(&mut *f, "Resume consensus sequence at [{},{}]: ", (*t).tps - 6, (*t).tps + 11);
    s = (*t).eseq.as_mut_ptr().add(((*t).tps - 7) as usize);
    for _ii in 0..18 {
        let _ = write!(&mut *f, "{}", cbase(*s) as char);
        s = s.add(1);
    }
    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f);
    disp_peptide_tag(f, t, sw);
}

// =============================================================================
// src_flatten/sequence.c lines 1180-1244
// void disp_tmrna_perm_seq(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_tmrna_perm_seq(f: *mut File, t: *mut gene, sw: *mut csw) {
    let mut i: i32;
    let mut s: *mut i32;
    let sb: *mut i32;
    let mut se: *mut i32;

    if (*t).nintron <= 0 {
        return;
    }
    if (*t).name[0] == 0 {
        let _ = writeln!(&mut *f, "tmRNA Sequence\n");
    } else {
        let _ = writeln!(&mut *f, "tmRNA sequence in {}\n", CStr::from_ptr((*t).name.as_ptr() as *const i8).to_string_lossy());
    }
    let _ = writeln!(&mut *f, "Permuted");
    let _ = writeln!(&mut *f, "1   .   10    .   20    .   30    .   40    .   50");
    sb = (*t).eseq.as_mut_ptr();
    s = sb;
    se = sb.add(54);
    i = 0;
    while s < se {
        let _ = write!(&mut *f, "{}", cpbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    se = sb.add((*t).intron as usize);
    while s < se {
        let _ = write!(&mut *f, "{}", cbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    se = sb.add((*t).asst as usize);
    while s < se {
        let _ = write!(&mut *f, "{}", cpbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    se = sb.add(((*t).asst + (*t).astem1 + (*t).dloop + (*t).cstem) as usize);
    while s < se {
        let _ = write!(&mut *f, "{}", cbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    se = sb.add((*t).tps as usize);
    while s < se {
        let _ = write!(&mut *f, "{}", cpbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    se = sb.add(((*t).tpe + 1) as usize);
    while ltranslate(se, t, sw) == b'*' {
        se = se.add(3);
    }
    while s < se {
        let _ = write!(&mut *f, "{}", cbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    se = sb.add(((*t).tpe + TMPTRAILER - 54) as usize);
    while s <= se {
        let _ = write!(&mut *f, "{}", cpbase(*s) as char);
        s = s.add(1);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
    }
    if i > 0 {
        let _ = writeln!(&mut *f);
    }
    let _ = writeln!(&mut *f, "\n5' tRNA domain at [{},{}]",
        (*t).asst + 1, (*t).asst + (*t).astem1 + (*t).dloop + (*t).cstem);
    let _ = writeln!(&mut *f, "3' tRNA domain at [55,{}]", (*t).intron);
    if ((*sw).secstructdisp & 2) != 0 {
        disp_tmrna_trnadomain_bracket_notation(f, t, sw);
    }
    let _ = write!(&mut *f, "Resume consensus sequence at [{},{}]: ", (*t).tps - 6, (*t).tps + 11);
    s = (*t).eseq.as_mut_ptr().add(((*t).tps - 7) as usize);
    for _ii in 0..18 {
        let _ = write!(&mut *f, "{}", cbase(*s) as char);
        s = s.add(1);
    }
    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f);
    disp_peptide_tag(f, t, sw);
}

// =============================================================================
// src_flatten/sequence.c lines 1247-1274
// void disp_cds(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_cds(f: *mut File, t: *mut gene, sw: *mut csw) {
    let mut i: i32;
    let mut c: u8;
    let mut s: *mut i32;

    let mut ncodon = (*t).nbase / 3;
    if (*t).tps == 0 {
        ncodon -= 1;
    }
    let _ = write!(&mut *f, "\n{} codons, start = {}{}{}, stop = ",
        ncodon,
        cbase((*t).seq[0]) as char,
        cbase((*t).seq[1]) as char,
        cbase((*t).seq[2]) as char);
    s = (*t).seq.as_mut_ptr().add(3);
    loop {
        i = *s;
        s = s.add(1);
        if i == TERM {
            break;
        }
        let _ = write!(&mut *f, "{}", cbase(i) as char);
    }
    if (*t).tps != 0 {
        let _ = write!(&mut *f, " incomplete");
    }
    let _ = writeln!(&mut *f, "\n1   .   10    .   20    .   30    .   40    .   50");
    s = (*t).eseq.as_mut_ptr();
    let mut se = s;
    while *se != TERM {
        se = se.add(1);
    }
    if (*t).tps != 0 {
        se = se.sub(3);
    }
    i = 0;
    while s < se {
        c = ltranslate(s, t, sw);
        let _ = write!(&mut *f, "{}", c as char);
        i += 1;
        if i >= 50 {
            let _ = writeln!(&mut *f);
            i = 0;
        }
        s = s.add(3);
    }
    if i > 0 {
        let _ = writeln!(&mut *f);
    }
    if (*sw).energydisp != 0 {
        let _ = writeln!(&mut *f, "Score = {}", (*t).energy);
    }
    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f);
}

// =============================================================================
// src_flatten/sequence.c lines 1287-1310
// void disp_gene(gene *t, char m[][MATY], csw *sw)
// =============================================================================
pub unsafe fn disp_gene(t: *mut gene, m: *mut [u8; MATY], sw: *mut csw) {
    let gc: f64;
    let mut stat: [u8; 80] = [0; 80];

    match (*t).genetype {
        TMRNA => {
            build_tmrna(t, m, 13, 27, sw);
            xcopy(m, 4, 3, b"tmRNA (tRNA domain)\0".as_ptr(), 19);
        }
        TRNA | _ => {
            build_trna(t, m, 13, 27, sw);
            name(t, stat.as_mut_ptr(), 1, sw);
            xcopy(m, 4, 3, stat.as_ptr(), length(stat.as_mut_ptr()));
        }
    }
    location(stat.as_mut_ptr(), t, sw, b"Sequence\0".as_ptr());
    xcopy(m, 4, 1, stat.as_ptr(), length(stat.as_mut_ptr()));
    gc = gc_content(t);
    let result = format!("{} bases, %GC = {:.1}\0", (*t).nbase, 100.0 * gc);
    let bytes = result.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        stat[i] = b;
    }
    xcopy(m, 4, 2, stat.as_ptr(), length(stat.as_mut_ptr()));
    if (*sw).reportpseudogenes != 0 {
        if pseudogene(t, sw) != 0 {
            xcopy(m, 4, 4, b"Possible Pseudogene\0".as_ptr(), 19);
        }
    }
    if (*sw).energydisp != 0 {
        let result = format!("Score = {}\n\0", (*t).energy);
        let bytes = result.as_bytes();
        for (i, &b) in bytes.iter().enumerate() {
            stat[i] = b;
        }
        xcopy(m, 4, 0, stat.as_ptr(), length(stat.as_mut_ptr()));
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1313-1378
// void disp_batch_trna(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_batch_trna(f: *mut File, t: *mut gene, sw: *mut csw) {
    let ls: i32;
    let mut ps: i32;
    let s: *mut i32;
    let mut anticodon: i32;
    let mut pos: [u8; 50] = [0; 50];
    let mut species: [u8; 50] = [0; 50];
    let mut m: [[u8; MATY]; MATX] = [[0; MATY]; MATX];
    static TYPE: [[u8; 6]; 2] = [*b"tRNA\0\0", *b"mtRNA\0"];
    static ASTERISK: [u8; 2] = [b' ', b'*'];

    s = (*t).seq.as_mut_ptr().add((*t).anticodon as usize);
    ps = if (*sw).reportpseudogenes != 0 {
        if pseudogene(t, sw) != 0 { 1 } else { 0 }
    } else {
        0
    };

    // Build species string
    match (*t).cloop {
        6 | 8 => {
            let type_str = CStr::from_ptr(TYPE[(*sw).mtrna as usize].as_ptr() as *const i8).to_string_lossy();
            let asterisk = ASTERISK[ps as usize] as char;
            let formatted = format!("{}-???{}", type_str, asterisk);
            let bytes = formatted.as_bytes();
            species[..bytes.len()].copy_from_slice(bytes);
            species[bytes.len()] = 0;
        }
        _ => {
            let type_str = CStr::from_ptr(TYPE[(*sw).mtrna as usize].as_ptr() as *const i8).to_string_lossy();
            let aa_str = CStr::from_ptr(aa(s, sw) as *const i8).to_string_lossy();
            let asterisk = ASTERISK[ps as usize] as char;
            let formatted = format!("{}-{}{}", type_str, aa_str, asterisk);
            let bytes = formatted.as_bytes();
            species[..bytes.len()].copy_from_slice(bytes);
            species[bytes.len()] = 0;
        }
    }

    position(pos.as_mut_ptr(), t, sw);
    ls = length(species.as_mut_ptr());

    let species_str = CStr::from_ptr(species.as_ptr() as *const i8).to_string_lossy();
    let pos_str = CStr::from_ptr(pos.as_ptr() as *const i8).to_string_lossy();
    if ls <= 10 {
        let _ = write!(&mut *f, "{:<10}{:>28}", species_str, pos_str);
    } else if ls <= 17 {
        let _ = write!(&mut *f, "{:<17}{:>21}", species_str, pos_str);
    } else {
        let _ = write!(&mut *f, "{:<25}{:>13}", species_str, pos_str);
    }

    if (*sw).energydisp != 0 {
        let _ = write!(&mut *f, "\t{:5.1}", (*t).energy);
    }

    anticodon = 1 + (*t).anticodon;
    if (*t).nintron > 0 {
        if (*t).intron <= (*t).anticodon {
            anticodon += (*t).nintron;
        }
    }
    let _ = write!(&mut *f, "\t{:<4}", anticodon);

    match (*t).cloop {
        6 => {
            let _ = write!(&mut *f, "\t({}{}) ", cbase(*s) as char, cbase(*s.add(1)) as char);
        }
        8 => {
            let _ = write!(&mut *f, "\t({}{}{}{}) ",
                cbase(*s) as char, cbase(*s.add(1)) as char,
                cbase(*s.add(2)) as char, cbase(*s.add(3)) as char);
        }
        _ => {
            let _ = write!(&mut *f, "\t({}{}{})",
                cbase(*s) as char, cbase(*s.add(1)) as char, cbase(*s.add(2)) as char);
        }
    }

    if (*t).nintron > 0 {
        let _ = write!(&mut *f, "i({},{})", (*t).intron + 1, (*t).nintron);
    }
    let _ = writeln!(&mut *f);

    if ((*sw).secstructdisp & 2) != 0 {
        disp_trna_bracket_notation(f, t, sw);
    }
    if ((*sw).secstructdisp & 4) != 0 {
        init_matrix(&mut m);
        build_trna(t, m.as_mut_ptr(), 13, 27, sw);
        disp_gene_SVG(t, m.as_mut_ptr(), sw);
    }
    if (*sw).seqdisp != 0 {
        disp_seq(f, t, sw);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1381-1409
// void disp_batch_tmrna(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_batch_tmrna(f: *mut File, t: *mut gene, sw: *mut csw) {
    let ps: i32;
    let mut tpe: i32;
    let sb: *mut i32;
    let mut se: *mut i32;
    let mut pos: [u8; 50] = [0; 50];
    let mut m: [[u8; MATY]; MATX] = [[0; MATY]; MATX];
    static PERMASK: [[[u8; 3]; 2]; 2] = [
        [[b' ', b' ', 0], [b'p', b' ', 0]],
        [[b'*', b' ', 0], [b'p', b'*', 0]],
    ];

    ps = if (*t).energy < 100.0 { 1 } else { 0 };
    position(pos.as_mut_ptr(), t, sw);
    let perm_idx = if (*t).asst == 0 { 0 } else { 1 };
    let perm_str = CStr::from_ptr(PERMASK[perm_idx][ps as usize].as_ptr() as *const i8).to_string_lossy();
    let pos_str = CStr::from_ptr(pos.as_ptr() as *const i8).to_string_lossy();
    let _ = write!(&mut *f, "tmRNA{:2}{:>31}", perm_str, pos_str);

    if (*sw).energydisp != 0 {
        let _ = write!(&mut *f, "\t{:5.1}\t", (*t).energy);
    }

    tpe = (*t).tpe;
    sb = (*t).eseq.as_mut_ptr().add((*t).tps as usize);
    se = (*t).eseq.as_mut_ptr().add((tpe + 1) as usize);
    while ltranslate(se, t, sw) == b'*' {
        se = se.add(3);
        tpe += 3;
    }
    let _ = write!(&mut *f, "\t{},{}\t", (*t).tps + 1, tpe + 1);
    let mut s = sb;
    while s < se {
        let _ = write!(&mut *f, "{}", ltranslate(s, t, sw) as char);
        s = s.add(3);
    }
    let _ = writeln!(&mut *f);

    if ((*sw).secstructdisp & 2) != 0 {
        disp_tmrna_trnadomain_bracket_notation(f, t, sw);
    }
    if ((*sw).secstructdisp & 4) != 0 {
        init_matrix(&mut m);
        build_tmrna(t, m.as_mut_ptr(), 13, 27, sw);
        disp_gene_SVG(t, m.as_mut_ptr(), sw);
    }
    if (*sw).seqdisp != 0 {
        disp_seq(f, t, sw);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1412-1422
// void disp_batch_srprna(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_batch_srprna(f: *mut File, t: *mut gene, sw: *mut csw) {
    let ps: i32;
    let mut pos: [u8; 50] = [0; 50];
    static ASTERISK: [u8; 2] = [b' ', b'*'];

    ps = if (*t).energy < 100.0 { 1 } else { 0 };
    position(pos.as_mut_ptr(), t, sw);
    let _ = write!(&mut *f, "srpRNA{}   {:>25}",
        ASTERISK[ps as usize] as char, CStr::from_ptr(pos.as_ptr() as *const i8).to_string_lossy());
    if (*sw).energydisp != 0 {
        let _ = write!(&mut *f, "\t{:5.1}", (*t).energy);
    }
    let _ = writeln!(&mut *f);
    if (*sw).seqdisp != 0 {
        disp_seq(f, t, sw);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1424-1434
// void disp_batch_cds(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn disp_batch_cds(f: *mut File, t: *mut gene, sw: *mut csw) {
    let ps: i32;
    let mut pos: [u8; 50] = [0; 50];
    static ASTERISK: [u8; 2] = [b' ', b'*'];

    ps = if (*t).energy < 100.0 { 1 } else { 0 };
    position(pos.as_mut_ptr(), t, sw);
    let _ = write!(&mut *f, "CDS{}      {:>25}",
        ASTERISK[ps as usize] as char, CStr::from_ptr(pos.as_ptr() as *const i8).to_string_lossy());
    if (*sw).energydisp != 0 {
        let _ = write!(&mut *f, "\t{:5.1}", (*t).energy);
    }
    let _ = writeln!(&mut *f);
    if (*sw).seqdisp != 0 {
        disp_seq(f, t, sw);
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1614-1650
// void trna_score(FILE *f, gene *t)
// =============================================================================
pub unsafe fn trna_score(f: *mut File, t: *mut gene) {
    let tpos: *mut i32;
    let tarm: i32;
    let mut varbp: i32 = 0;
    let ea: f64;
    let mut eta: f64;
    let evls: f64;
    static BEM: [[f64; 6]; 6] = [
        [-2.144, -0.428, -2.144, ATBOND, 0.000, 0.000],
        [-0.428, -2.144, 3.000, -2.144, 0.000, 0.000],
        [-2.144, 3.000, -2.144, 1.286, 0.000, 0.000],
        [ATBOND, -2.144, 1.286, -0.428, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
    ];
    static A_SCORE: [f64; 6] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    static C_SCORE: [f64; 6] = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    static G_SCORE: [f64; 6] = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    static T_SCORE: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];

    if (*t).genetype != tRNA {
        return;
    }
    tarm = 2 * (*t).tstem + (*t).tloop;
    tpos = (*t).seq.as_mut_ptr().add(
        ((*t).astem1 + (*t).spacer1 + (*t).dloop + 2 * (*t).dstem + 1 +
         2 * (*t).cstem + (*t).cloop + (*t).var) as usize);
    let mut s = tpos.add(((*t).tstem - 1) as usize);
    eta = 6.0 * (G_SCORE[*s as usize] + T_SCORE[*s.add(1) as usize] +
                 T_SCORE[*s.add(2) as usize] + C_SCORE[*s.add(3) as usize]) +
          3.0 * A_SCORE[*s.add(1) as usize];
    s = s.add(((*t).tloop - 3) as usize);
    eta += 2.0 * (G_SCORE[*s as usize] + A_SCORE[*s.add(1) as usize] +
                  T_SCORE[*s.add(3) as usize] + C_SCORE[*s.add(4) as usize] +
                  C_SCORE[*s.add(5) as usize]);
    eta += astem_energy(tpos, tpos.add(tarm as usize), (*t).tstem);
    eta += BEM[*tpos.add((*t).tstem as usize) as usize][*tpos.add(((*t).tstem + 4) as usize) as usize];
    eta -= 3.0 * (5 - (*t).tstem) as f64;
    if (*t).tloop > 7 {
        eta -= 3.0 * ((*t).tloop - 7) as f64;
    } else {
        eta -= 3.0 * (7 - (*t).tloop) as f64;
    }
    s = (*t).seq.as_mut_ptr();
    if (*t).astem1 > 7 {
        s = s.add(1);
    }
    ea = astem_energy(s, tpos.add((tarm + 7) as usize), 7);
    if (*t).var > 17 {
        evls = vloop_stability(tpos.sub((*t).var as usize), (*t).var, &mut varbp);
    } else {
        evls = 0.0;
    }
    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f, "               T-arm score: {}", eta);
    let _ = writeln!(&mut *f, "              A-stem score: {}", ea);
    let _ = writeln!(&mut *f, "          V-loop stability: {}", evls);
    let _ = writeln!(&mut *f);
}

// =============================================================================
// src_flatten/sequence.c lines 1653-1790
// void tmrna_score(FILE *f, gene *t, csw *sw)
// =============================================================================
pub unsafe fn tmrna_score(f: *mut File, t: *mut gene, sw: *mut csw) {
    let mut r: i32;
    let mut j: i32;
    let mut te: i32;
    let mut s: *mut i32;
    let mut sb: *mut i32;
    let mut se: *mut i32;
    let mut tpos: *mut i32;
    let tarm: i32;
    let mut e: f64;
    let er: f64;
    let et: f64;
    let eal: f64;
    let esp: f64;
    let ed: f64;
    let ec: f64;
    let ea: f64;
    let mut egga: f64;
    let etcca: f64;
    let egg: f64;
    let mut edgg: f64;
    let mut eta: f64;
    let ehairpin: f64;
    let euhairpin: f64;

    static GTEM: [i32; 6] = [0x00, 0x00, 0x11, 0x00, 0x00, 0x00];
    static TAGEND_SCORE: [f64; 4] = [36.0, 66.0, 62.0, 72.0];
    static NPS: [i32; 126] = [
        0,0,0,0,  // 0-3
        0,0,0,0,  // 4-7
        0,0,0,0,  // 8-11
        1,1,1,1,  // 12-15
        0,0,0,0,  // 16-19
        1,1,1,1,  // 20-23
        0,0,0,0,  // 24-27
        1,1,1,1,  // 28-31
        0,0,0,0,  // 32-35
        1,1,1,1,  // 36-39
        1,1,1,1,  // 40-43
        1,1,1,1,  // 44-47
        2,1,2,1,  // 48-51
        0,0,0,0,  // 52-55
        2,1,1,1,  // 56-59
        1,1,1,1,  // 60-63
        0,0,0,0,  // 64-67
        0,0,0,0,  // 68-71
        0,0,0,0,  // 72-75
        0,0,0,0,  // 76-79
        0,0,0,0,  // 80-83
        0,0,0,0,  // 84-87
        // remaining 38 elements (88-125) are 0
        0,0,0,0,0,0,0,0,0,0,  // 88-97
        0,0,0,0,0,0,0,0,0,0,  // 98-107
        0,0,0,0,0,0,0,0,0,0,  // 108-117
        0,0,0,0,0,0,0,0,      // 118-125
    ];
    static BEM: [[f64; 6]; 6] = [
        [-2.144, -0.428, -2.144, ATBOND, 0.000, 0.000],
        [-0.428, -2.144, 3.000, -2.144, 0.000, 0.000],
        [-2.144, 3.000, -2.144, 1.286, 0.000, 0.000],
        [ATBOND, -2.144, 1.286, -0.428, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
    ];
    static A: [f64; 6] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    static C: [f64; 6] = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    static G: [f64; 6] = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    static K: [f64; 6] = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0];
    static R: [f64; 6] = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    static T: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    static Y: [f64; 6] = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0];
    static NA: [f64; 6] = [0.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    static NV: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    static NM: [f64; 6] = [0.0, 0.0, 1.0, 1.0, 1.0, 1.0];

    if (*t).genetype != tmRNA {
        return;
    }

    tarm = 2 * (*t).tstem + (*t).tloop;
    s = (*t).eseq.as_mut_ptr().offset((*t).tps as isize - 7);
    er = (A[*s.offset(1) as usize] + 2.0 * T[*s.offset(2) as usize] + C[*s.offset(2) as usize]
        + 3.0 * A[*s.offset(3) as usize] + R[*s.offset(4) as usize] + Y[*s.offset(6) as usize]
        + 3.0 * G[*s.offset(7) as usize] + C[*s.offset(8) as usize]
        - if (*sw).tmstrict != 0 {
            NV[*s.offset(10) as usize] + NV[*s.offset(11) as usize]
            + NM[*s.offset(14) as usize] + NA[*s.offset(17) as usize]
        } else { 0.0 }) * 4.0;

    s = (*t).eseq.as_mut_ptr().offset((*t).tpe as isize - 8);
    te = ((NPS[((*s.offset(0) << 4) + (*s.offset(1) << 2) + *s.offset(2)) as usize] & 1) << 1)
        | (NPS[((*s.offset(3) << 4) + (*s.offset(4) << 2) + *s.offset(5)) as usize] & 1);
    et = TAGEND_SCORE[te as usize];

    if (*sw).tmstrict != 0 {
        let mut eal_val = 0.0;
        j = -3;
        while j < 6 {
            te = *s.offset(j as isize);
            j += 1;
            te = (te << 2) | *s.offset(j as isize);
            j += 1;
            if te == 9 {
                eal_val = (11 + 2 * ((j + 1) / 3)) as f64;
            }
            j += 1;
        }
        eal = eal_val;
        ehairpin = find_taghairpin(s.offset(8));
        euhairpin = find_tag_upstream_hairpin((*t).eseq.as_mut_ptr().offset((*t).tps as isize - 10));
    } else {
        eal = 15.0;
        ehairpin = 16.0;
        euhairpin = 15.0;
    }

    tpos = (*t).eseq.as_mut_ptr();
    if (*t).asst > 0 {
        tpos = tpos.offset(((*t).cstem + (*t).var + 54) as isize);
        ed = 0.0;
    } else {
        tpos = tpos.offset(((*t).astem1 + (*t).dloop + 2 * (*t).cstem + (*t).nintron + (*t).var) as isize);
        ed = 0.001 * ((*t).tps - (tpos as isize - (*t).eseq.as_mut_ptr() as isize) as i32) as f64;
    }

    s = tpos.offset((*t).tstem as isize - 10);
    e = K[*s.offset(0) as usize] + G[*s.offset(1) as usize] + A[*s.offset(2) as usize];
    egga = K[*s.offset(1) as usize] + G[*s.offset(2) as usize] + A[*s.offset(3) as usize];
    if e > egga { egga = e; }
    egga *= 6.0;
    if egga < 18.0 { egga = 0.0; }

    s = tpos.offset((tarm + 4) as isize);
    etcca = 10.0 * (T[*s.offset(0) as usize] + C[*s.offset(1) as usize]
        + C[*s.offset(2) as usize] + A[*s.offset(3) as usize]);

    s = (*t).eseq.as_mut_ptr().offset((*t).asst as isize);
    egg = 7.0 * (G[*s.offset(1) as usize] + G[*s.offset(2) as usize]);

    edgg = 0.0;
    s = (*t).eseq.as_mut_ptr().offset(((*t).asst + (*t).astem1) as isize);
    sb = s.offset(3);
    se = s.offset(7);
    r = GTEM[*sb as usize];
    sb = sb.offset(1);
    while sb < se {
        r = (r >> 4) + GTEM[*sb as usize];
        sb = sb.offset(1);
        if (r & 3) == 2 {
            edgg = 14.0;
            break;
        }
    }

    s = tpos.offset((*t).tstem as isize - 1);
    if (*sw).tmstrict != 0 && (*t).asst == 0 {
        eta = 6.0 * (G[*s.offset(0) as usize] + T[*s.offset(1) as usize]
            + T[*s.offset(2) as usize] + C[*s.offset(3) as usize])
            + 3.0 * A[*s.offset(1) as usize];
    } else {
        eta = 6.0 * (G[*s.offset(0) as usize]
            + (G[*s.offset(1) as usize] + T[*s.offset(1) as usize])
            + (G[*s.offset(2) as usize] + T[*s.offset(2) as usize])
            + C[*s.offset(3) as usize])
            + 3.0 * A[*s.offset(1) as usize];
    }
    s = s.offset(((*t).tloop - 3) as isize);
    eta += 2.0 * (G[*s as usize] + A[*s.offset(1) as usize] + T[*s.offset(3) as usize]
        + C[*s.offset(4) as usize] + C[*s.offset(5) as usize]);
    eta += astem_energy(tpos, tpos.offset(tarm as isize), (*t).tstem);
    eta += BEM[*tpos.offset((*t).tstem as isize) as usize][*tpos.offset(((*t).tstem + 4) as isize) as usize];
    eta -= 3.0 * (5 - (*t).tstem) as f64;
    if (*t).tloop > 7 {
        eta -= 3.0 * ((*t).tloop - 7) as f64;
    } else {
        eta -= 3.0 * (7 - (*t).tloop) as f64;
    }
    eta *= 1.59;

    s = (*t).eseq.as_mut_ptr().offset(((*t).asst + (*t).astem1 + (*t).dloop) as isize);
    ec = stem_energy(s, tpos.offset(-((*t).var as isize)), (*t).cstem);

    s = (*t).eseq.as_mut_ptr().offset((*t).asst as isize);
    ea = astem_energy(s, tpos.offset((tarm + (*t).astem1) as isize), (*t).astem1);

    esp = if ((*t).tpe - (*t).tps) < 24 { -15.0 } else { 0.0 };

    e = er + et + ed + eal + esp + egga + egg + etcca + eta + ec + ea + edgg + ehairpin + euhairpin;

    let _ = writeln!(&mut *f);
    let _ = writeln!(&mut *f, "     Resume sequence score: {}", er);
    let _ = writeln!(&mut *f, "Resume-Tarm distance score: {}", ed);
    let _ = writeln!(&mut *f, "         Tag peptide score: {}", et);
    let _ = writeln!(&mut *f, "     Tag end alanine score: {}", eal);
    let _ = writeln!(&mut *f, "         Short tag penalty: {}", esp);
    let _ = writeln!(&mut *f, "         Tag hairpin score: {}", ehairpin);
    let _ = writeln!(&mut *f, "Tag upstream hairpin score: {}", euhairpin);
    let _ = writeln!(&mut *f, "          V-loop GGA score: {}", egga);
    let _ = writeln!(&mut *f, "           A-stem GG score: {}", egg);
    let _ = writeln!(&mut *f, "         A-stem TCCA score: {}", etcca);
    let _ = writeln!(&mut *f, "           D-loop GG score: {}", edgg);
    let _ = writeln!(&mut *f, "               T-arm score: {}", eta);
    let _ = writeln!(&mut *f, "              C-stem score: {}", ec);
    let _ = writeln!(&mut *f, "              A-stem score: {}", ea);
    let _ = writeln!(&mut *f, "     C-stem + A-stem score: {}", ea + ec);
    let _ = writeln!(&mut *f, "               Total score: {}", e);
    let _ = writeln!(&mut *f, "          Normalised score: {}", nenergy(t, sw));
    let _ = writeln!(&mut *f);
}

// =============================================================================
// src_flatten/sequence.c lines 1807-1836
// int identify_tag(char tag[], int len, char (*thit)[50], int nt)
// =============================================================================
pub unsafe fn identify_tag(tag: *const u8, len: i32, thit: *mut [u8; 50], nt: i32) -> i32 {
    let mut n: i32 = 0;
    let mut too_many: i32 = 0;

    // st = tag + len; while (*--st == '*');
    let mut st = tag.offset(len as isize);
    loop {
        st = st.offset(-1);
        if *st != b'*' { break; }
    }

    // for (i = 0; i < NTAG && !too_many; i++)
    let mut i: i32 = 0;
    while i < NTAG as i32 && too_many == 0 {
        let mut s = st;
        let sb = TAGDATABASE[i as usize].tag.as_ptr();
        let mut sd = sb;

        // while (*++sd);  -- find end of tag string
        loop {
            sd = sd.offset(1);
            if *sd == 0 { break; }
        }

        // while (*s-- == *--sd)
        loop {
            sd = sd.offset(-1);
            if *s != *sd { break; }
            s = s.offset(-1);

            let mut do_partial: i32 = 0;

            if s < tag {
                if sd > sb {
                    do_partial = 1;
                } else {
                    if n >= nt {
                        too_many = 1;
                        break;
                    }
                    copy(TAGDATABASE[i as usize].name.as_ptr(), (*thit.offset(n as isize)).as_mut_ptr());
                    n += 1;
                    break;
                }
            } else if sd > sb {
                continue;
            } else {
                do_partial = 1;
            }

            if do_partial != 0 {
                if n >= nt {
                    too_many = 1;
                    break;
                }
                let s_out = copy(TAGDATABASE[i as usize].name.as_ptr(), (*thit.offset(n as isize)).as_mut_ptr());
                copy(b" (partial match)\0".as_ptr(), s_out);
                n += 1;
                break;
            }
        }
        i += 1;
    }

    if too_many != 0 { -1 } else { n }
}

// =============================================================================
// src_flatten/sequence.c lines 1858-1903
// void update_tmrna_tag_database(gene ts[], int nt, csw *sw)
// =============================================================================
pub unsafe fn update_tmrna_tag_database(ts: *mut gene, nt: i32, sw: *mut csw) {
    let mut species: [u8; STRLEN] = [0; STRLEN];
    let mut tag_buf: [u8; 100] = [0; 100];

    if (*sw).tagend >= NTAGMAX as i32 {
        return;
    }

    for i in 0..nt {
        let t = ts.offset(i as isize);
        if (*t).genetype != tmRNA {
            continue;
        }

        // Find last '|' in name
        let mut s = (*t).name.as_ptr();
        let mut se: *const u8 = std::ptr::null();
        while *s != 0 {
            if *s == b'|' {
                se = s;
            }
            s = s.offset(1);
        }

        if se.is_null() || *se == 0 {
            continue;
        }

        // while (++se) if (space(*se)) break;
        loop {
            se = se.offset(1);
            if *se == 0 || space(*se) {
                break;
            }
        }
        if *se == 0 {
            continue;
        }

        // while (++se) if (!space(*se)) break;
        loop {
            se = se.offset(1);
            if *se == 0 || !space(*se) {
                break;
            }
        }
        if *se == 0 {
            continue;
        }

        // Check for " sp. " pattern
        if !softstrpos(se as *mut u8, b" sp. \0".as_ptr() as *mut u8).is_null() {
            let mut sp = softstrpos(se as *mut u8, b"two-piece\0".as_ptr() as *mut u8);
            if sp.is_null() {
                sp = softstrpos(se as *mut u8, b"tmRNA\0".as_ptr() as *mut u8);
            }
            if sp.is_null() {
                continue;
            }
            // while (space(sp[-1])) sp--;
            while space(*sp.offset(-1)) {
                sp = sp.offset(-1);
            }
            copy2sp(se as *mut u8, sp, species.as_mut_ptr(), 49);
        } else {
            // Extract first two words
            let mut s_out = species.as_mut_ptr();
            let mut c: i32 = 2;
            let mut se_mut = se;
            while *se_mut != 0 {
                if space(*se_mut) {
                    c -= 1;
                    if c <= 0 {
                        break;
                    }
                }
                *s_out = *se_mut;
                s_out = s_out.offset(1);
                se_mut = se_mut.offset(1);
            }
            *s_out = 0;
        }

        // Check if species already in database
        let mut k: i32 = 0;
        while k < (*sw).tagend {
            if !softstrpos(TAGDATABASE[k as usize].name.as_ptr() as *mut u8, species.as_mut_ptr()).is_null() {
                break;
            }
            k += 1;
        }
        if k < (*sw).tagend {
            continue;
        }

        // Add new entry
        copy(species.as_ptr(), TAGDATABASE[(*sw).tagend as usize].name.as_mut_ptr());

        let lx = peptide_tag(tag_buf.as_mut_ptr(), 50, t, sw);
        let mut s_tag = tag_buf.as_mut_ptr().offset((lx - 1) as isize);
        while *s_tag == b'*' {
            s_tag = s_tag.offset(-1);
        }
        s_tag = s_tag.offset(1);
        *s_tag = 0;

        copy(tag_buf.as_ptr(), TAGDATABASE[(*sw).tagend as usize].tag.as_mut_ptr());

        (*sw).tagend += 1;
        if (*sw).tagend >= NTAGMAX as i32 {
            break;
        }
    }
}

// =============================================================================
// src_flatten/sequence.c lines 1917-1933
// void report_new_tmrna_tags(csw *sw)
// =============================================================================
pub unsafe fn report_new_tmrna_tags(sw: *mut csw) {
    let f = (*sw).f;
    let mut sort: [i32; NTAGMAX] = [0; NTAGMAX];

    // Insertion sort by name
    for n in 0..(*sw).tagend {
        let mut k = n;
        while k > 0 {
            k -= 1;
            if string_compare(
                TAGDATABASE[n as usize].name.as_ptr(),
                TAGDATABASE[sort[k as usize] as usize].name.as_ptr()
            ) >= 0 {
                k += 1;
                break;
            }
            sort[(k + 1) as usize] = sort[k as usize];
        }
        sort[k as usize] = n;
    }

    let _ = writeln!(&mut *f, "\ntmRNA tag database update:");

    for k in 0..(*sw).tagend {
        let n = sort[k as usize];
        let _ = writeln!(&mut *f, "     {{ \"{}\",\"{}\"}},",
            CStr::from_ptr(TAGDATABASE[n as usize].name.as_ptr() as *const i8).to_string_lossy(),
            CStr::from_ptr(TAGDATABASE[n as usize].tag.as_ptr() as *const i8).to_string_lossy()
        );
    }

    let _ = writeln!(&mut *f, "\n{} tmRNA peptide tags", (*sw).tagend);
    let _ = writeln!(&mut *f, "{} new tmRNA peptide tags\n", (*sw).tagend - NTAG as i32);
}
