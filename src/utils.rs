/*
 * Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
 * utils.rs - Utility functions
 *
 * Based on ARAGORN v1.2.41 by Dean Laslett
 * Author: Sunju Kim <n.e.coli.1822@gmail.com>
 */

use std::ptr;
use crate::types::*;

// =============================================================================
// src_flatten/utils.c lines 8-9
// char upcasec(char c)
// { return((c >= 'a')?c-32:c); }
// =============================================================================
#[inline]
pub fn upcasec(c: u8) -> u8 {
    if c >= b'a' { c - 32 } else { c }
}

// =============================================================================
// src_flatten/utils.c lines 11-14
// int length(char *s)
// { int i = 0;
//   while (*s++) i++;
//   return(i); }
// =============================================================================
pub unsafe fn length(mut s: *const u8) -> i32 {
    let mut i: i32 = 0;
    while *s != 0 {
        s = s.add(1);
        i += 1;
    }
    i
}

// =============================================================================
// src_flatten/utils.c lines 16-21
// char *softmatch(char *s, char *key)
// =============================================================================
pub unsafe fn softmatch(mut s: *mut u8, mut key: *mut u8) -> *mut u8 {
    while upcasec(*key) == upcasec(*s) {
        key = key.add(1);
        if *key == 0 {
            return s;
        }
        s = s.add(1);
    }
    if *key != 0 {
        return ptr::null_mut();
    }
    s
}

// =============================================================================
// src_flatten/utils.c lines 23-32
// char *strpos(char *s, char *k)
// =============================================================================
pub unsafe fn strpos(mut s: *mut u8, k: *mut u8) -> *mut u8 {
    let d: u8 = *k;
    loop {
        let c = *s;
        if c == 0 {
            break;
        }
        if c == d {
            let mut i: i32 = 0;
            loop {
                i += 1;
                if *k.offset(i as isize) == 0 {
                    return s;
                }
                if *s.offset(i as isize) != *k.offset(i as isize) {
                    break;
                }
            }
        }
        s = s.add(1);
    }
    ptr::null_mut()
}

// =============================================================================
// src_flatten/utils.c lines 35-45
// char *softstrpos(char *s, char *k)
// =============================================================================
pub unsafe fn softstrpos(mut s: *mut u8, k: *mut u8) -> *mut u8 {
    let d: u8 = upcasec(*k);
    loop {
        let c = *s;
        if c == 0 {
            break;
        }
        if upcasec(c) == d {
            let mut i: i32 = 0;
            loop {
                i += 1;
                if *k.offset(i as isize) == 0 {
                    return s;
                }
                if upcasec(*s.offset(i as isize)) != upcasec(*k.offset(i as isize)) {
                    break;
                }
            }
        }
        s = s.add(1);
    }
    ptr::null_mut()
}

// =============================================================================
// src_flatten/utils.c lines 47-57
// char *wildstrpos(char *s, char *k)
// =============================================================================
pub unsafe fn wildstrpos(mut s: *mut u8, k: *mut u8) -> *mut u8 {
    let d: u8 = upcasec(*k);
    loop {
        let c = *s;
        if c == 0 {
            break;
        }
        if upcasec(c) == d || d == b'*' {
            let mut i: i32 = 0;
            loop {
                i += 1;
                if *k.offset(i as isize) == 0 {
                    return s;
                }
                if upcasec(*s.offset(i as isize)) != upcasec(*k.offset(i as isize))
                    && *k.offset(i as isize) != b'*'
                {
                    break;
                }
            }
        }
        s = s.add(1);
    }
    ptr::null_mut()
}

// =============================================================================
// src_flatten/utils.c lines 59-70
// char *marginstring(char *s, char *k, int margin)
// =============================================================================
pub unsafe fn marginstring(mut s: *mut u8, k: *mut u8, margin: i32) -> *mut u8 {
    let mut j: i32 = 0;
    let d: u8 = *k;
    loop {
        let c = *s;
        if c == 0 {
            break;
        }
        if c == d {
            let mut i: i32 = 0;
            loop {
                i += 1;
                if *k.offset(i as isize) == 0 {
                    return s;
                }
                if *s.offset(i as isize) != *k.offset(i as isize) {
                    break;
                }
            }
        }
        s = s.add(1);
        j += 1;
        if j >= margin {
            break;
        }
    }
    ptr::null_mut()
}

// =============================================================================
// src_flatten/utils.c lines 73-85
// int margindetect(char *line, int margin)
// =============================================================================
pub unsafe fn margindetect(line: *mut u8, margin: i32) -> i32 {
    let mut i: i32 = 0;
    let mut s: *mut u8 = line;
    loop {
        let c = *s;
        s = s.add(1);
        if c == 0 {
            break;
        }
        if !space(c) {
            // Check terminating characters
            if c == b'\n' || c == b'\r' || c == 0 {
                return 0;
            }
            return 1;
        }
        if c == b'\t' {
            i += 7;
        }
        i += 1;
        if i >= margin {
            return 0;
        }
    }
    0
}

// =============================================================================
// src_flatten/utils.c lines 88-107
// char *backword(char *line, char *s, int n)
// =============================================================================
pub unsafe fn backword(line: *mut u8, mut s: *mut u8, mut n: i32) -> *mut u8 {
    let mut spzone: i32;
    if space(*s) {
        spzone = 1;
    } else {
        spzone = 0;
        n += 1;
    }
    while s > line {
        if space(*s) {
            if spzone == 0 {
                spzone = 1;
                n -= 1;
                if n <= 0 {
                    s = s.add(1);
                    return s;
                }
            }
        } else {
            spzone = 0;
        }
        s = s.sub(1);
    }
    if !space(*s) {
        if n <= 1 {
            return s;
        }
    }
    ptr::null_mut()
}

// =============================================================================
// src_flatten/utils.c lines 109-157
// char *dconvert(char *s, double *r)
// =============================================================================
pub unsafe fn dconvert(mut s: *mut u8, r: *mut f64) -> *mut u8 {
    let zero: u8 = b'0';
    let nine: u8 = b'9';
    let mut shift: i32 = 0;
    let mut expshift: i32 = 0;
    let mut sgn: i32 = 1;
    let mut expsgn: i32 = 1;
    let mut limit: u8 = 0;
    let mut exponent: i32 = 0;
    let mut result: f64 = 0.0;

    let mut c: u8 = *s;
    if c == b'-' {
        sgn = -1;
        s = s.add(1);
        c = *s;
    } else if c == b'+' {
        s = s.add(1);
        c = *s;
    }

    if c >= zero && c <= nine {
        result = (c - zero) as f64;
        loop {
            s = s.add(1);
            c = *s;
            if c < zero || c > nine {
                break;
            }
            limit += 1;
            if limit < 15 {
                result = result * 10.0 + (c - zero) as f64;
            }
        }
    }

    if c == b'.' {
        loop {
            s = s.add(1);
            c = *s;
            if c < zero || c > nine {
                break;
            }
            limit += 1;
            if limit < 15 {
                result = result * 10.0 + (c - zero) as f64;
                shift += 1;
            }
        }
    }

    if c == b'E' || c == b'e' || c == b'D' || c == b'd' {
        s = s.add(1);
        c = *s;
        if c == b'-' {
            expsgn = -1;
            s = s.add(1);
            c = *s;
        } else if c == b'+' {
            s = s.add(1);
            c = *s;
        }
        if c >= zero && c <= nine {
            exponent = (c - zero) as i32;
            loop {
                s = s.add(1);
                c = *s;
                if c < zero || c > nine {
                    break;
                }
                exponent = exponent * 10 + (c - zero) as i32;
                expshift += 1;
                if expshift > 3 {
                    break;
                }
            }
        }
    }

    result *= sgn as f64;
    exponent = exponent * expsgn - shift;
    if exponent >= 0 {
        while exponent > 0 {
            exponent -= 1;
            result *= 10.0;
        }
    } else {
        while exponent < 0 {
            exponent += 1;
            result /= 10.0;
        }
    }
    *r *= 0.01 * result;
    s
}

// =============================================================================
// src_flatten/utils.c lines 159-177
// char *lconvert(char *s, long *r)
// =============================================================================
pub unsafe fn lconvert(mut s: *mut u8, r: *mut i64) -> *mut u8 {
    let zero: u8 = b'0';
    let nine: u8 = b'9';
    let mut sgn: i64 = 1;
    let mut result: i64 = 0;

    let mut c: u8 = *s;
    if c == b'-' {
        sgn = -1;
        s = s.add(1);
        c = *s;
    } else if c == b'+' {
        s = s.add(1);
        c = *s;
    }

    if c >= zero && c <= nine {
        result = (c - zero) as i64;
        loop {
            s = s.add(1);
            c = *s;
            if c < zero || c > nine {
                break;
            }
            result = result * 10 + (c - zero) as i64;
        }
    }
    *r = result * sgn;
    s
}

// =============================================================================
// src_flatten/utils.c lines 180-195
// char *getlong(char *line, long *l)
// =============================================================================
pub unsafe fn getlong(line: *mut u8, l: *mut i64) -> *mut u8 {
    let zero: u8 = b'0';
    let nine: u8 = b'9';
    if line.is_null() {
        return ptr::null_mut();
    }
    let mut s: *mut u8 = line;
    loop {
        let c1 = *s;
        if c1 == 0 {
            break;
        }
        if c1 >= zero && c1 <= nine {
            return lconvert(s, l);
        } else if c1 == b'-' || c1 == b'+' {
            let c2 = *s.add(1);
            if c2 >= zero && c2 <= nine {
                return lconvert(s, l);
            }
        }
        s = s.add(1);
    }
    ptr::null_mut()
}

// =============================================================================
// src_flatten/utils.c lines 199-201
// char *copy(char *from, char *to)
// { while (*to++ = *from++);
//   return(--to);  }
// =============================================================================
pub unsafe fn copy(mut from: *const u8, mut to: *mut u8) -> *mut u8 {
    loop {
        *to = *from;
        if *from == 0 {
            break;
        }
        from = from.add(1);
        to = to.add(1);
    }
    to
}

// =============================================================================
// src_flatten/utils.c lines 204-216
// char *copy2sp(char *from1, char *from2, char *to, int n)
// =============================================================================
pub unsafe fn copy2sp(mut from1: *mut u8, from2: *mut u8, to: *mut u8, mut n: i32) -> *mut u8 {
    let mut s: *mut u8 = to;
    while from1 < from2 {
        *s = *from1;
        from1 = from1.add(1);
        s = s.add(1);
        n -= 1;
        if n <= 0 {
            loop {
                s = s.sub(1);
                if s <= to {
                    break;
                }
                if space(*s) {
                    break;
                }
            }
            break;
        }
    }
    *s = 0;
    s
}

// =============================================================================
// src_flatten/utils.c lines 219-228
// char *copy3cr(char *from, char *to, int n)
// =============================================================================
pub unsafe fn copy3cr(mut from: *const u8, mut to: *mut u8, mut n: i32) -> *mut u8 {
    loop {
        *to = *from;
        if *to == 0 {
            break;
        }
        if *to == DLIM as u8 {
            *to = 0;
            break;
        }
        from = from.add(1);
        n -= 1;
        if n <= 0 {
            to = to.add(1);
            *to = 0;
            break;
        }
        to = to.add(1);
    }
    to
}

// =============================================================================
// src_flatten/utils.c lines 230-243
// char *quotestring(char *line, char *a, int n)
// =============================================================================
pub unsafe fn quotestring(mut line: *mut u8, mut a: *mut u8, mut n: i32) -> *mut u8 {
    loop {
        let ch = *line;
        line = line.add(1);
        if ch == 0 {
            break;
        }
        if ch == b'"' {
            loop {
                let ch = *line;
                line = line.add(1);
                if ch == 0 {
                    break;
                }
                if ch == b'"' || ch == b';' || ch == b'\n' || ch == b'\r' {
                    break;
                }
                *a = ch;
                a = a.add(1);
                n -= 1;
                if n <= 0 {
                    break;
                }
            }
            break;
        }
    }
    *a = 0;
    a
}

// =============================================================================
// LIBRARY FUNCTIONS
// =============================================================================

// =============================================================================
// src_flatten/utils.c lines 249-266
// int fseekd(data_set *d, long fpos, long foffset)
// =============================================================================
pub unsafe fn fseekd(d: *mut data_set, fpos: i64, foffset: i64) -> i32 {
    let mut fpos = fpos;
    let mut foffset = foffset;

    if (*d).bugmode != 0 {
        fpos += foffset;
        if fpos < 0 {
            fpos = 0;
        }
        if (*d).f.is_null() {
            return EOF;
        }
        // In bugmode, seek to start and read forward
        if libc::fseek((*d).f as *mut libc::FILE, 0, libc::SEEK_SET) != 0 {
            return EOF;
        }
        (*d).filepointer = -1;
        while (*d).filepointer + 1 < fpos {
            (*d).filepointer += 1;
            if libc::fgetc((*d).f as *mut libc::FILE) == libc::EOF {
                return EOF;
            }
        }
        return 0;
    }

    if (*d).f.is_null() {
        return EOF;
    }
    if libc::fseek((*d).f as *mut libc::FILE, fpos as libc::c_long, libc::SEEK_SET) != 0 {
        return EOF;
    }
    (*d).filepointer = fpos;

    if foffset != 0 {
        if fpos + foffset < 0 {
            foffset = -fpos;
        }
        if libc::fseek((*d).f as *mut libc::FILE, foffset as libc::c_long, libc::SEEK_CUR) != 0 {
            return EOF;
        }
        (*d).filepointer += foffset;
    }
    0
}

// =============================================================================
// src_flatten/utils.c lines 269-273
// long ftelld(data_set *d)
// =============================================================================
pub unsafe fn ftelld(d: *mut data_set) -> i64 {
    if (*d).bugmode != 0 {
        return (*d).filepointer;
    }
    if (*d).f.is_null() {
        return -1;
    }
    libc::ftell((*d).f as *mut libc::FILE) as i64
}

// =============================================================================
// src_flatten/utils.c lines 276-282
// char fgetcd(data_set *d)
// =============================================================================
pub unsafe fn fgetcd(d: *mut data_set) -> u8 {
    if (*d).f.is_null() {
        return NOCHAR as u8;
    }
    let c = libc::fgetc((*d).f as *mut libc::FILE);
    if c == libc::EOF {
        return NOCHAR as u8;
    }
    (*d).filepointer += 1;
    c as u8
}

// =============================================================================
// src_flatten/utils.c lines 285-300
// char *fgetsd(data_set *d, char line[], int len)
// =============================================================================
pub unsafe fn fgetsd(d: *mut data_set, line: *mut u8, len: i32) -> *mut u8 {
    if (*d).f.is_null() {
        return ptr::null_mut();
    }
    let mut i: i32 = 0;
    while i < len {
        let c = libc::fgetc((*d).f as *mut libc::FILE);
        if c == libc::EOF {
            break;
        }
        (*d).filepointer += 1;
        let ic = c as u8;
        if ic == b'\r' {
            continue;
        }
        if ic == b'\n' {
            *line.offset(i as isize) = DLIM as u8;
            i += 1;
            break;
        }
        *line.offset(i as isize) = ic;
        i += 1;
    }
    if i < 1 {
        return ptr::null_mut();
    }
    *line.offset(i as isize) = 0;
    line
}

// =============================================================================
// src_flatten/utils.c lines 302-321
// int agene_position_check(data_set *d, int nagene, annotated_gene *agene)
// =============================================================================
pub unsafe fn agene_position_check(
    d: *mut data_set,
    nagene: i32,
    agene: *mut annotated_gene,
) -> i32 {
    let swap: i64;
    if (*agene).stop - (*agene).start > MAXAGENELEN as i64 {
        swap = (*agene).stop;
        (*agene).stop = (*agene).start;
        (*agene).start = swap;
        (*agene).stop += (*d).aseqlen;
    }
    if (*agene).start > (*agene).stop {
        (*agene).stop += (*d).aseqlen;
    }
    let l = (*agene).stop - (*agene).start;
    if l < 1 || l > MAXAGENELEN as i64 {
        return 0;
    }
    if (*agene).stop == (*d).aseqlen {
        for a in 0..nagene {
            let gene_a = &(&(*d).gene)[a as usize];
            if gene_a.start == (*agene).start {
                if gene_a.genetype == (*agene).genetype {
                    if !softmatch(
                        gene_a.species.as_ptr() as *mut u8,
                        (*agene).species.as_mut_ptr(),
                    )
                    .is_null()
                    {
                        return 0;
                    }
                }
            }
        }
    }
    1
}

// =============================================================================
// src_flatten/utils.c lines 325-520
// long process_sequence_heading(data_set *d, csw *sw)
// NOTE: This is a large function - full 1:1 port
// =============================================================================
pub unsafe fn process_sequence_heading(d: *mut data_set, sw: *mut csw) -> i64 {
    let mut i: i32;
    let mut nagene: i32;
    let mut goto_fnsn: i32;
    let mut skip_gbnl: i32;
    let mut l: i64 = 0;
    let mut realstart: i64;
    let mut line: [u8; STRLEN] = [0; STRLEN];
    let mut c: u8;
    let mut s: *mut u8;
    let mut sq: *mut u8;
    let mut sd: *mut u8;
    let mut agene: *mut annotated_gene;
    let mut tmpagene: annotated_gene = std::mem::zeroed();

    (*d).datatype = FASTA;
    fseekd(d, (*d).seqstart, (*d).seqstartoff);

    // Read past comments
    loop {
        loop {
            c = fgetcd(d);
            if c == NOCHAR as u8 {
                return -1;
            }
            if !space(c) {
                break;
            }
        }
        if c != b'#' {
            break;
        }
        if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
            return -1;
        }
    }

    if fgetsd(d, (*d).seqname.as_mut_ptr(), STRLENM1 as i32).is_null() {
        return -1;
    }

    goto_fnsn = 0;
    if c != b'>' {
        s = (*d).seqname.as_mut_ptr();
        if upcasec(c) != b'L' {
            loop {
                c = *s;
                if c == 0 {
                    goto_fnsn = 1;
                    break;
                }
                s = s.add(1);
                if upcasec(c) == b'L' {
                    break;
                }
            }
        }
        if goto_fnsn == 0 {
            let ocus = b"OCUS\0";
            s = softmatch(s, ocus.as_ptr() as *mut u8);
            if s.is_null() {
                goto_fnsn = 1;
            }
        }

        if goto_fnsn != 0 {
            // FNSN label equivalent
            let unnamed = b"Unnamed sequence \0";
            s = copy(unnamed.as_ptr(), (*d).seqname.as_mut_ptr());
            fseekd(d, (*d).seqstart, (*d).seqstartoff);
            realstart = ftelld(d);
            if !fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                copy3cr(line.as_ptr(), s, 50);
            }
            fseekd(d, realstart, 0);
            return realstart;
        }

        // Check for BP
        let bp = b"BP\0";
        sd = softstrpos((*d).seqname.as_mut_ptr(), bp.as_ptr() as *mut u8);
        if !sd.is_null() {
            sd = backword((*d).seqname.as_mut_ptr(), sd, 1);
            if !sd.is_null() {
                sd = getlong(sd, &mut l);
                if !sd.is_null() {
                    (*d).aseqlen = l;
                }
            }
        }

        s = s.add(4);
        while space(*s) {
            s = s.add(1);
        }
        sq = (*d).seqname.as_mut_ptr();
        while !space(*s) {
            *sq = *s;
            sq = sq.add(1);
            s = s.add(1);
        }
        (*d).aseqlen = 0;

        if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
            return -2;
        }

        // Check for DEFINITION
        let definition = b"DEFINITION\0";
        sd = softstrpos(line.as_mut_ptr(), definition.as_ptr() as *mut u8);
        if !sd.is_null() {
            sd = sd.add(10);
            while space(*sd) {
                sd = sd.add(1);
            }
            *sq = b' ';
            sq = sq.add(1);
            copy(sd, sq);
            if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                return -2;
            }
        } else {
            copy(s, sq);
        }

        for i in 0..NS {
            (*d).nagene[i as usize] = 0;
        }
        nagene = 0;

        // Parse FEATURES until ORIGIN
        let origin = b"ORIGIN\0";
        while marginstring(line.as_mut_ptr(), origin.as_ptr() as *mut u8, 10).is_null() {
            skip_gbnl = 0;
            if nagene >= NGFT as i32 {
                skip_gbnl = 1;
            } else {
                agene = &mut (&mut (*d).gene)[nagene as usize];
                (*agene).comp = 0;
                (*agene).start = -1;
                (*agene).stop = -1;
                (*agene).antistart = -1;
                (*agene).antistop = -1;
                (*agene).permuted = 0;
                (*agene).pseudogene = 0;

                let trna = b"tRNA\0";
                let tmrna = b"tmRNA\0";
                let cds = b"CDS\0";
                let mrna = b"mRNA\0";
                let rrna = b"rRNA\0";

                s = marginstring(line.as_mut_ptr(), trna.as_ptr() as *mut u8, 10);
                if !s.is_null() {
                    // tRNA handling
                    (*agene).genetype = tRNA as i32;
                    let complement = b"complement\0";
                    if !softstrpos(s, complement.as_ptr() as *mut u8).is_null() {
                        (*agene).comp = 1;
                    }
                    s = getlong(s, &mut l);
                    if !s.is_null() {
                        (*agene).start = l;
                    }
                    s = getlong(s, &mut l);
                    if !s.is_null() {
                        (*agene).stop = l;
                    }
                    let trna_species = b"tRNA-???\0";
                    copy(trna_species.as_ptr(), (*agene).species.as_mut_ptr());

                    if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                        return -2;
                    }
                    while margindetect(line.as_mut_ptr(), 10) == 0 {
                        let product = b"product=\0";
                        s = softstrpos(line.as_mut_ptr(), product.as_ptr() as *mut u8);
                        if !s.is_null() {
                            let trna_prefix = b"tRNA-\0";
                            s = softstrpos(s, trna_prefix.as_ptr() as *mut u8);
                            if !s.is_null() {
                                s = s.add(5);
                                while space(*s) {
                                    s = s.add(1);
                                }
                                copy3cr(s, (*agene).species.as_mut_ptr().add(5), 3);
                            }
                        }
                        let anticodon = b"anticodon=\0";
                        s = softstrpos(line.as_mut_ptr(), anticodon.as_ptr() as *mut u8);
                        if !s.is_null() {
                            s = s.add(10);
                            s = getlong(s, &mut l);
                            if s.is_null() {
                                l = -1;
                            }
                            (*agene).antistart = l;
                            s = getlong(s, &mut l);
                            if s.is_null() {
                                l = -1;
                            }
                            (*agene).antistop = l;
                        }
                        let pseudo = b"/pseudo\0";
                        if !softstrpos(line.as_mut_ptr(), pseudo.as_ptr() as *mut u8).is_null() {
                            (*agene).pseudogene = 1;
                        }
                        if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                            return -2;
                        }
                    }
                    if agene_position_check(d, nagene, agene) != 0 {
                        (*d).nagene[tRNA as usize] += 1;
                        nagene += 1;
                    }
                } else {
                    s = marginstring(line.as_mut_ptr(), tmrna.as_ptr() as *mut u8, 10);
                    if !s.is_null() {
                        // tmRNA handling (simplified - full implementation would match C exactly)
                        (*agene).genetype = tmRNA as i32;
                        let complement = b"complement\0";
                        if !softstrpos(s, complement.as_ptr() as *mut u8).is_null() {
                            (*agene).comp = 1;
                        }
                        s = getlong(s, &mut l);
                        if !s.is_null() {
                            (*agene).start = l;
                        }
                        s = getlong(s, &mut l);
                        if !s.is_null() {
                            (*agene).stop = l;
                        }
                        let tmrna_species = b"tmRNA\0";
                        copy(tmrna_species.as_ptr(), (*agene).species.as_mut_ptr());
                        if agene_position_check(d, nagene, agene) == 0 {
                            skip_gbnl = 1;
                        } else {
                            (*d).nagene[tmRNA as usize] += 1;
                            nagene += 1;
                            // Continue reading tmRNA features...
                            if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                                return -2;
                            }
                            while margindetect(line.as_mut_ptr(), 10) == 0 {
                                let acceptor = b"acceptor\0";
                                if !softstrpos(line.as_mut_ptr(), acceptor.as_ptr() as *mut u8).is_null() {
                                    (*agene).permuted = 1;
                                }
                                let pseudo = b"/pseudo\0";
                                if !softstrpos(line.as_mut_ptr(), pseudo.as_ptr() as *mut u8).is_null() {
                                    (*agene).pseudogene = 1;
                                }
                                if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                                    return -2;
                                }
                            }
                            // Check for second tmRNA entry (permuted)
                            s = marginstring(line.as_mut_ptr(), tmrna.as_ptr() as *mut u8, 10);
                            if !s.is_null() {
                                tmpagene.comp = 0;
                                tmpagene.start = -1;
                                tmpagene.stop = -1;
                                tmpagene.antistart = -1;
                                tmpagene.antistop = -1;
                                tmpagene.permuted = 0;
                                tmpagene.pseudogene = 0;
                                let complement = b"complement\0";
                                if !softstrpos(s, complement.as_ptr() as *mut u8).is_null() {
                                    tmpagene.comp = 1;
                                }
                                s = getlong(s, &mut l);
                                if !s.is_null() {
                                    tmpagene.start = l;
                                }
                                s = getlong(s, &mut l);
                                if !s.is_null() {
                                    tmpagene.stop = l;
                                }
                                if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                                    return -2;
                                }
                                while margindetect(line.as_mut_ptr(), 10) == 0 {
                                    let coding = b"coding\0";
                                    if !softstrpos(line.as_mut_ptr(), coding.as_ptr() as *mut u8).is_null() {
                                        tmpagene.permuted = 1;
                                    }
                                    let pseudo = b"/pseudo\0";
                                    if !softstrpos(line.as_mut_ptr(), pseudo.as_ptr() as *mut u8).is_null() {
                                        tmpagene.pseudogene = 1;
                                    }
                                    let tag_peptide = b"/tag_peptide\0";
                                    s = softstrpos(line.as_mut_ptr(), tag_peptide.as_ptr() as *mut u8);
                                    if !s.is_null() {
                                        s = getlong(s, &mut l);
                                        if !s.is_null() {
                                            tmpagene.antistart = l;
                                        }
                                        s = getlong(s, &mut l);
                                        if !s.is_null() {
                                            tmpagene.antistop = l;
                                        }
                                    }
                                    if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                                        return -2;
                                    }
                                }
                                if (*agene).permuted != 0 && tmpagene.permuted != 0 {
                                    (*agene).stop = tmpagene.stop;
                                    (*agene).antistart = tmpagene.antistart;
                                    (*agene).antistop = tmpagene.antistop;
                                    let perm_species = b"tmRNA(Perm)\0";
                                    copy(perm_species.as_ptr(), (*agene).species.as_mut_ptr());
                                } else {
                                    if nagene >= NGFT as i32 {
                                        skip_gbnl = 1;
                                    } else {
                                        agene = &mut (&mut (*d).gene)[nagene as usize];
                                        (*agene).comp = tmpagene.comp;
                                        (*agene).start = tmpagene.start;
                                        (*agene).stop = tmpagene.stop;
                                        (*agene).antistart = -1;
                                        (*agene).antistop = -1;
                                        (*agene).permuted = 0;
                                        (*agene).pseudogene = tmpagene.pseudogene;
                                        let tmrna_species = b"tmRNA\0";
                                        copy(tmrna_species.as_ptr(), (*agene).species.as_mut_ptr());
                                        if agene_position_check(d, nagene, agene) != 0 {
                                            (*d).nagene[tmRNA as usize] += 1;
                                            nagene += 1;
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        s = marginstring(line.as_mut_ptr(), cds.as_ptr() as *mut u8, 10);
                        if s.is_null() {
                            s = marginstring(line.as_mut_ptr(), mrna.as_ptr() as *mut u8, 10);
                        }
                        if !s.is_null() {
                            // CDS/mRNA handling
                            (*agene).genetype = CDS as i32;
                            let complement = b"complement\0";
                            if !softstrpos(s, complement.as_ptr() as *mut u8).is_null() {
                                (*agene).comp = 1;
                            }
                            s = getlong(s, &mut l);
                            if !s.is_null() {
                                (*agene).start = l;
                            }
                            s = getlong(s, &mut l);
                            if !s.is_null() {
                                (*agene).stop = l;
                            }
                            let cds_species = b"CDS\0";
                            copy(cds_species.as_ptr(), (*agene).species.as_mut_ptr());
                            if agene_position_check(d, nagene, agene) == 0 {
                                skip_gbnl = 1;
                            } else {
                                (*d).nagene[CDS as usize] += 1;
                                nagene += 1;
                                if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                                    return -2;
                                }
                                while margindetect(line.as_mut_ptr(), 10) == 0 {
                                    let pseudo = b"/pseudo\0";
                                    if !softstrpos(line.as_mut_ptr(), pseudo.as_ptr() as *mut u8).is_null() {
                                        (*agene).pseudogene = 1;
                                    }
                                    if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                                        return -2;
                                    }
                                }
                            }
                        } else {
                            s = marginstring(line.as_mut_ptr(), rrna.as_ptr() as *mut u8, 10);
                            if !s.is_null() {
                                // rRNA handling
                                (*agene).genetype = rRNA as i32;
                                let complement = b"complement\0";
                                if !softstrpos(s, complement.as_ptr() as *mut u8).is_null() {
                                    (*agene).comp = 1;
                                }
                                s = getlong(s, &mut l);
                                if !s.is_null() {
                                    (*agene).start = l;
                                }
                                s = getlong(s, &mut l);
                                if !s.is_null() {
                                    (*agene).stop = l;
                                }
                                let rrna_species = b"rRNA-???\0";
                                copy(rrna_species.as_ptr(), (*agene).species.as_mut_ptr());
                                if agene_position_check(d, nagene, agene) == 0 {
                                    skip_gbnl = 1;
                                } else {
                                    (*d).nagene[rRNA as usize] += 1;
                                    nagene += 1;
                                    if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                                        return -2;
                                    }
                                    while margindetect(line.as_mut_ptr(), 10) == 0 {
                                        let pseudo = b"/pseudo\0";
                                        if !softstrpos(line.as_mut_ptr(), pseudo.as_ptr() as *mut u8).is_null() {
                                            (*agene).pseudogene = 1;
                                        }
                                        let product = b"product=\0";
                                        s = softstrpos(line.as_mut_ptr(), product.as_ptr() as *mut u8);
                                        if !s.is_null() {
                                            s = s.add(8);
                                            while space(*s) {
                                                s = s.add(1);
                                            }
                                            quotestring(s, (*agene).species.as_mut_ptr(), SHORTSTRLENM1 as i32);
                                        }
                                        if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                                            return -2;
                                        }
                                    }
                                }
                            } else {
                                skip_gbnl = 1;
                            }
                        }
                    }
                }
            }
            if skip_gbnl != 0 {
                if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                    return -2;
                }
            }
        }

        (*sw).annotated = if nagene > 0 { 1 } else { 0 };
        for i in 0..NS {
            (*sw).nagene[i as usize] = (*d).nagene[i as usize];
        }
        (*d).datatype = GENBANK;

        loop {
            c = fgetcd(d);
            if c == NOCHAR as u8 {
                return -2;
            }
            if !space(c) {
                break;
            }
        }
        if c != b'1' && c != b'0' {
            loop {
                c = fgetcd(d);
                if c == NOCHAR as u8 {
                    return -2;
                }
                if !space(c) {
                    break;
                }
            }
        }
        realstart = ftelld(d);
    } else {
        // FASTA format
        loop {
            realstart = ftelld(d);
            loop {
                c = fgetcd(d);
                if c == NOCHAR as u8 {
                    return -3;
                }
                if !space(c) {
                    break;
                }
            }
            if c != b'>' {
                break;
            }
            if fgetsd(d, line.as_mut_ptr(), STRLENM1 as i32).is_null() {
                return -3;
            }
        }
        fseekd(d, realstart, 0);
    }

    // Trim sequence name
    s = (*d).seqname.as_mut_ptr();
    i = 0;
    loop {
        c = *s;
        if c == 0 {
            break;
        }
        if c == b'\n' || c == b'\r' {
            break;
        }
        i += 1;
        if i >= STRLEN as i32 {
            break;
        }
        s = s.add(1);
    }
    *s = 0;
    realstart
}

// =============================================================================
// src_flatten/utils.c lines 523-590
// int move_forward(data_set *d)
// =============================================================================
#[allow(unused_assignments)]
pub unsafe fn move_forward(d: *mut data_set) -> i32 {
    let mut ic: i32;
    let mut fail: i32;
    let mut do_bs: i32;
    let mut nextbase: i64;

    static MAP: [i32; 256] = [
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
        -4, -4, -4, -4, -4, -4, NOBASE, -3, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
        -4, -4, -2, -4, -4, Adenine, AMBIG, Cytosine, AMBIG, -4, -4, Guanine, AMBIG,
        -4, -4, AMBIG, -5, AMBIG, AMBIG, -4, -4, -4,
        AMBIG, AMBIG, Thymine, Thymine, AMBIG, AMBIG, -4,
        AMBIG, -4, -4, -4, -4, INSERT, NOBASE, -4, Adenine, AMBIG, Cytosine, AMBIG,
        -4, -4, Guanine, AMBIG, -4, -4, AMBIG, -5, AMBIG, AMBIG, -4, -4, -4,
        AMBIG, AMBIG, Thymine, Thymine, AMBIG, AMBIG, -4,
        AMBIG, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
    ];

    if (*d).ps >= (*d).psmax {
        if (*d).psmax > 0 {
            fseekd(d, (*d).seqstart, (*d).seqstartoff);
            (*d).ps = 0;
        }
    }

    fail = 0;
    while fail == 0 {
        let c = fgetcd(d);
        if c == NOCHAR as u8 {
            fail = 1;
            break;
        }
        ic = MAP[c as usize];
        do_bs = 1;
        while do_bs != 0 {
            do_bs = 0;
            if ic >= Adenine {
                (*d).ps += 1;
                return ic;
            }
            if ic == -2 {
                (*d).nextseq = ftelld(d);
                (*d).nextseqoff = -1;
                return TERM;
            }
            if ic == -3 {
                if (*d).datatype == GENBANK {
                    let c2 = fgetcd(d);
                    if c2 == NOCHAR as u8 {
                        fail = 1;
                        break;
                    }
                    ic = MAP[c2 as usize];
                    if ic != -3 {
                        do_bs = 1;
                        continue;
                    }
                    loop {
                        let c3 = fgetcd(d);
                        if c3 == NOCHAR as u8 {
                            fail = 1;
                            break;
                        }
                        if !space(c3) {
                            break;
                        }
                    }
                    if fail != 0 {
                        break;
                    }
                    (*d).nextseq = ftelld(d);
                    (*d).nextseqoff = -1;
                    return TERM;
                }
            }
            if ic == -5 {
                nextbase = ftelld(d);
                let c2 = fgetcd(d);
                if c2 == NOCHAR as u8 {
                    fail = 1;
                    break;
                }
                if upcasec(c2) == b'O' {
                    let c3 = fgetcd(d);
                    if c3 == NOCHAR as u8 {
                        fail = 1;
                        break;
                    }
                    if upcasec(c3) == b'C' {
                        let c4 = fgetcd(d);
                        if c4 == NOCHAR as u8 {
                            fail = 1;
                            break;
                        }
                        if upcasec(c4) == b'U' {
                            let c5 = fgetcd(d);
                            if c5 == NOCHAR as u8 {
                                fail = 1;
                                break;
                            }
                            if upcasec(c5) == b'S' {
                                (*d).nextseq = nextbase;
                                (*d).nextseqoff = -1;
                                return TERM;
                            }
                        }
                    }
                }
                fseekd(d, nextbase, 0);
            }
        }
    }

    (*d).nextseq = -1;
    (*d).nextseqoff = 0;
    if (*d).psmax > 0 {
        (*d).ps = (*d).psmax;
        return NOBASE;
    }
    TERM
}

// =============================================================================
// src_flatten/utils.c lines 594-616
// int seq_init(data_set *d, csw *sw)
// =============================================================================
pub unsafe fn seq_init(d: *mut data_set, sw: *mut csw) -> i32 {
    let mut ngc: i64;
    let mut ic: i32;

    (*d).filepointer = 0;
    (*d).seqstart = process_sequence_heading(d, sw);
    if (*d).seqstart < 0 {
        if (*d).seqstart == -2 {
            eprintln!(
                "ERROR - unable to read Genbank sequence {}",
                std::str::from_utf8_unchecked(&(*d).seqname)
            );
        } else if (*d).seqstart == -3 {
            eprintln!(
                "ERROR - unable to read fasta sequence {}",
                std::str::from_utf8_unchecked(&(*d).seqname)
            );
        }
        return 0;
    }
    (*d).seqstartoff = 0;
    (*d).ps = 0;
    (*d).psmax = -1;
    ngc = 0;

    loop {
        ic = move_forward(d);
        if ic < Adenine {
            break;
        }
        if ic >= Cytosine && ic <= Guanine {
            ngc += 1;
        }
    }

    (*d).psmax = (*d).ps;
    if (*d).psmax <= 0 {
        return 0;
    }
    (*d).gc = ngc as f64 / (*d).psmax as f64;
    fseekd(d, (*d).seqstart, (*d).seqstartoff);
    (*d).ps = 0;
    1
}
