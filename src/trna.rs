/*
 * Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
 * trna.rs - tRNA detection helper functions
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

use std::io::Write;

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::types::*;
use crate::sequence::{find_taghairpin, find_tag_upstream_hairpin, name, location, position, nenergy};
use crate::utils::{copy, copy3cr};
use crate::types::sq;
use crate::tables::TS;
use crate::thread_state::get_astem5_state;

// Local alias for NO_GENE constant (types.rs uses NO_GENE, C code uses noGENE)
const noGENE: i32 = NO_GENE;

/// find_tstems - Find T-stem candidates in sequence
/// trna.c:12-77
pub unsafe fn find_tstems(
    s: *mut i32,
    ls: i32,
    hit: *mut trna_loop,
    nh: i32,
    sw: *mut csw
) -> i32 {
    let mut i: i32;
    let mut r: i32;
    let mut c: i32;
    let mut tstem: i32;
    let mut tloop: i32;
    let mut ithresh1: i32;
    let mut too_many: i32;
    let mut s1: *mut i32;
    let mut s2: *mut i32;
    let mut se: *mut i32;
    let mut ss: *mut i32;
    let mut si: *mut i32;
    let mut sb: *mut i32;
    let mut sc: *mut i32;
    let mut sf: *mut i32;
    let mut sl: *mut i32;
    let mut sx: *mut i32;
    let mut tem: *mut i32;
    let mut ec: f64;
    let mut energy: f64;
    let mut penalty: f64;
    let mut thresh2: f64;

    // Flattened 1D array for consistency with other optimized arrays
    static bem_flat: [f64; 36] = [
        -2.144, -0.428, -2.144, ATBOND, 0.000, 0.000,
        -0.428, -2.144,  3.000, -2.144, 0.000, 0.000,
        -2.144,  3.000, -2.144,  1.286, 0.000, 0.000,
        ATBOND, -2.144,  1.286, -0.428, 0.000, 0.000,
         0.000,  0.000,  0.000,  0.000, 0.000, 0.000,
         0.000,  0.000,  0.000,  0.000, 0.000, 0.000,
    ];
    static mut A: [f64; 6] = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    static mut C: [f64; 6] = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0];
    static mut G: [f64; 6] = [0.0, 0.0, 2.0, 0.0, 0.0, 0.0];
    static mut T: [f64; 6] = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
    static mut tem_trna: [i32; 6] = [0x0100, 0x0002, 0x2000, 0x0220, 0x0000, 0x0000];
    static mut tem_tmrna: [i32; 6] = [0x0100, 0x0002, 0x2220, 0x0220, 0x0000, 0x0000];

    i = 0;
    too_many = 0;
    tem = if (*sw).tmrna != 0 || (*sw).threshlevel < 1.0 {
        tem_tmrna.as_mut_ptr()
    } else {
        tem_trna.as_mut_ptr()
    };
    ithresh1 = (*sw).ttscanthresh as i32;
    thresh2 = (*sw).ttarmthresh;
    ss = s.offset((*sw).loffset as isize);
    si = ss.offset(4 - 1);
    sl = s.offset((ls - (*sw).roffset + 5 + 3) as isize);
    r = *tem.offset(*si as isize);
    si = si.offset(1);
    r = (r >> 4) + *tem.offset(*si as isize);
    si = si.offset(1);
    r = (r >> 4) + *tem.offset(*si as isize);
    si = si.offset(1);

    while si < sl && too_many == 0 {
        r = (r >> 4) + *tem.offset(*si as isize);
        si = si.offset(1);
        c = r & 0xF;
        if c < ithresh1 {
            continue;
        }
        sb = si.offset(-7);
        sf = sb.offset(13);
        ec = (3 * c) as f64;

        tstem = 4;
        while tstem <= 5 && too_many == 0 {
            if sb < sl.offset(-8) {
                sc = sf;
                sx = si.offset(-2);
                tloop = 5;
                while tloop <= 9 && too_many == 0 {
                    if tloop > 7 {
                        penalty = 3.0 * (tloop - tstem - 2) as f64;
                    } else {
                        penalty = 3.0 * (12 - tloop - tstem) as f64;
                    }
                    s1 = sb;
                    s2 = sc;
                    se = s1.offset(tstem as isize);
                    s2 = s2.offset(-1);
                    energy = ec + *bem_flat.get_unchecked((*se as usize) * 6 + (*se.offset(4) as usize))
                           + *bem_flat.get_unchecked((*s1 as usize) * 6 + (*s2 as usize)) - penalty;
                    s1 = s1.offset(1);
                    while s1 < se {
                        s2 = s2.offset(-1);
                        energy += *bem_flat.get_unchecked((*s1 as usize) * 6 + (*s2 as usize));
                        s1 = s1.offset(1);
                    }
                    energy += *G.get_unchecked(*sx as usize) + *A.get_unchecked(*sx.offset(1) as usize)
                            + *T.get_unchecked(*sx.offset(3) as usize) + *C.get_unchecked(*sx.offset(4) as usize)
                            + *C.get_unchecked(*sx.offset(5) as usize);
                    if energy >= thresh2 {
                        if i >= nh {
                            eprintln!("Too many tstem hits");
                            too_many = 1;
                        } else {
                            (*hit.offset(i as isize)).pos = sb;
                            (*hit.offset(i as isize)).loop_ = tloop;
                            (*hit.offset(i as isize)).stem = tstem;
                            (*hit.offset(i as isize)).energy = energy;
                            i += 1;
                        }
                    }
                    sx = sx.offset(1);
                    sc = sc.offset(1);
                    tloop += 1;
                }
            }
            sb = sb.offset(-1);
            if sb < ss {
                break;
            }
            sf = sf.offset(1);
            tstem += 1;
        }
    }
    i
}

/// find_astem5 - Find 5' acceptor stem candidates
/// trna.c:80-133
#[inline]
pub unsafe fn find_astem5(
    mut si: *mut i32,
    mut sl: *mut i32,
    astem3: *mut i32,
    n3: i32,
    hit: *mut trna_loop,
    nh: i32,
    sw: *mut csw
) -> i32 {
    let mut i: i32;
    let mut k: i32;
    let mut too_many: i32;
    let mut s1: *mut i32;
    let mut s2: *mut i32;
    let mut se: *mut i32;
    let mut r: u32;
    let mut tascanthresh: u32;
    let mut tastemthresh: f64;
    let mut energy: f64;

    // Thread-local state for tem array
    let tl_state = &mut *get_astem5_state();
    let tem = &mut tl_state.tem;

    static A: [u32; 6] = [0, 0, 0, 2, 0, 0];
    static C: [u32; 6] = [0, 0, 2, 0, 0, 0];
    static G: [u32; 6] = [0, 2, 0, 1, 0, 0];
    static T: [u32; 6] = [2, 0, 1, 0, 0, 0];
    // Flattened 1D array for better cache performance (6x6 = 36 elements)
    // Access: abem_flat[row * 6 + col] instead of abem[row][col]
    static abem_flat: [f64; 36] = [
        -2.144, -0.428, -2.144, ATBOND, 0.000, 0.000,  // row 0
        -0.428, -2.144,  3.000, -2.144, 0.000, 0.000,  // row 1
        -2.144,  3.000, -2.144,  1.286, 0.000, 0.000,  // row 2
        ATBOND, -2.144,  1.286, -0.428, 0.000, 0.000,  // row 3
         0.000,  0.000,  0.000,  0.000, 0.000, 0.000,  // row 4
         0.000,  0.000,  0.000,  0.000, 0.000, 0.000,  // row 5
    ];

    tascanthresh = (*sw).tascanthresh as u32;
    tastemthresh = (*sw).tastemthresh;
    i = 0;
    too_many = 0;
    sl = sl.offset(n3 as isize);
    se = astem3.offset(n3 as isize - 1);

    *tem.get_unchecked_mut(0) = *A.get_unchecked(*se as usize);
    *tem.get_unchecked_mut(1) = *C.get_unchecked(*se as usize);
    *tem.get_unchecked_mut(2) = *G.get_unchecked(*se as usize);
    *tem.get_unchecked_mut(3) = *T.get_unchecked(*se as usize);
    se = se.offset(-1);
    while se >= astem3 {
        *tem.get_unchecked_mut(0) = (*tem.get_unchecked(0) << 4) + *A.get_unchecked(*se as usize);
        *tem.get_unchecked_mut(1) = (*tem.get_unchecked(1) << 4) + *C.get_unchecked(*se as usize);
        *tem.get_unchecked_mut(2) = (*tem.get_unchecked(2) << 4) + *G.get_unchecked(*se as usize);
        *tem.get_unchecked_mut(3) = (*tem.get_unchecked(3) << 4) + *T.get_unchecked(*se as usize);
        se = se.offset(-1);
    }

    r = *tem.get_unchecked(*si as usize);
    si = si.offset(1);
    k = 1;
    k += 1;
    while k < n3 {
        r = (r >> 4) + *tem.get_unchecked(*si as usize);
        si = si.offset(1);
        k += 1;
    }

    // Specialized unrolled version for n3=7 (most common case)
    if n3 == 7 {
        // Pre-compute astem3 row indices (multiplied by 6) - these don't change during scan
        let a0_row = (*astem3.offset(0) as i32) * 6;
        let a1_row = (*astem3.offset(1) as i32) * 6;
        let a2_row = (*astem3.offset(2) as i32) * 6;
        let a3_row = (*astem3.offset(3) as i32) * 6;
        let a4_row = (*astem3.offset(4) as i32) * 6;
        let a5_row = (*astem3.offset(5) as i32) * 6;
        let a6_row = (*astem3.offset(6) as i32) * 6;

        #[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
        {
            // AVX2 SIMD version using gather instructions
            // _mm256_set_epi32 order: (pos7, pos6, pos5, pos4, pos3, pos2, pos1, pos0)
            // _mm256_castsi256_si128 extracts pos0-pos3 (lower 128 bits)
            // So we put indices in positions 0-3: a0, a1, a2, a3
            let row_lo = _mm256_set_epi32(0, 0, 0, 0, a3_row, a2_row, a1_row, a0_row);
            let row_hi = _mm256_set_epi32(0, 0, 0, 0, 0, a6_row, a5_row, a4_row);

            while si < sl && too_many == 0 {
                r = (r >> 4) + *tem.get_unchecked(*si as usize);
                si = si.offset(1);
                if (r & 15) >= tascanthresh {
                    // Get column indices
                    let b0 = *si.offset(-1);
                    let b1 = *si.offset(-2);
                    let b2 = *si.offset(-3);
                    let b3 = *si.offset(-4);
                    let b4 = *si.offset(-5);
                    let b5 = *si.offset(-6);
                    let b6 = *si.offset(-7);

                    // Build indices for gather: row + col
                    let col_lo = _mm256_set_epi32(0, 0, 0, 0, b3, b2, b1, b0);
                    let col_hi = _mm256_set_epi32(0, 0, 0, 0, 0, b6, b5, b4);
                    let idx_lo = _mm256_add_epi32(row_lo, col_lo);
                    let idx_hi = _mm256_add_epi32(row_hi, col_hi);

                    // Gather: idx_lo positions 0-3 give us values for pairs 0-3
                    // Gather: idx_hi positions 0-2 give us values for pairs 4-6
                    let vals_lo = _mm256_i32gather_pd(
                        abem_flat.as_ptr(),
                        _mm256_castsi256_si128(idx_lo),
                        8  // scale=8 for f64
                    );
                    let vals_hi = _mm256_i32gather_pd(
                        abem_flat.as_ptr(),
                        _mm256_castsi256_si128(idx_hi),
                        8
                    );

                    // Sum all 7 values (vals_lo has 4 valid, vals_hi has 3 valid + 1 zero)
                    // Horizontal sum
                    let sum_vec = _mm256_add_pd(vals_lo, vals_hi);
                    // Extract and add: [a+e, b+f, c+g, d+0]
                    let hi128 = _mm256_extractf128_pd(sum_vec, 1);  // [c+g, d+0]
                    let lo128 = _mm256_castpd256_pd128(sum_vec);    // [a+e, b+f]
                    let sum2 = _mm_add_pd(lo128, hi128);            // [a+e+c+g, b+f+d+0]
                    let sum3 = _mm_hadd_pd(sum2, sum2);             // [sum, sum]
                    energy = _mm_cvtsd_f64(sum3);

                    if energy >= tastemthresh {
                        if i >= nh {
                            eprintln!("Too many astem5 hits");
                            too_many = 1;
                        } else {
                            (*hit.offset(i as isize)).pos = si.offset(-7);
                            (*hit.offset(i as isize)).energy = energy;
                            i += 1;
                        }
                    }
                }
            }
        }

        #[cfg(not(all(target_arch = "x86_64", target_feature = "avx2")))]
        {
            // Scalar fallback
            while si < sl && too_many == 0 {
                r = (r >> 4) + *tem.get_unchecked(*si as usize);
                si = si.offset(1);
                if (r & 15) >= tascanthresh {
                    let b0 = *si.offset(-1) as usize;
                    let b1 = *si.offset(-2) as usize;
                    let b2 = *si.offset(-3) as usize;
                    let b3 = *si.offset(-4) as usize;
                    let b4 = *si.offset(-5) as usize;
                    let b5 = *si.offset(-6) as usize;
                    let b6 = *si.offset(-7) as usize;

                    energy = *abem_flat.get_unchecked(a0_row as usize + b0)
                           + *abem_flat.get_unchecked(a1_row as usize + b1)
                           + *abem_flat.get_unchecked(a2_row as usize + b2)
                           + *abem_flat.get_unchecked(a3_row as usize + b3)
                           + *abem_flat.get_unchecked(a4_row as usize + b4)
                           + *abem_flat.get_unchecked(a5_row as usize + b5)
                           + *abem_flat.get_unchecked(a6_row as usize + b6);

                    if energy >= tastemthresh {
                        if i >= nh {
                            eprintln!("Too many astem5 hits");
                            too_many = 1;
                        } else {
                            (*hit.offset(i as isize)).pos = si.offset(-7);
                            (*hit.offset(i as isize)).energy = energy;
                            i += 1;
                        }
                    }
                }
            }
        }
    } else {
        // Generic version for other n3 values
        while si < sl && too_many == 0 {
            r = (r >> 4) + *tem.get_unchecked(*si as usize);
            si = si.offset(1);
            if (r & 15) >= tascanthresh {
                s1 = astem3;
                s2 = si;
                se = s1.offset(n3 as isize);
                s2 = s2.offset(-1);
                // Use flattened 1D array: abem_flat[row * 6 + col]
                energy = *abem_flat.get_unchecked((*s1 as usize) * 6 + (*s2 as usize));
                s1 = s1.offset(1);
                while s1 < se {
                    s2 = s2.offset(-1);
                    energy += *abem_flat.get_unchecked((*s1 as usize) * 6 + (*s2 as usize));
                    s1 = s1.offset(1);
                }
                if energy >= tastemthresh {
                    if i >= nh {
                        eprintln!("Too many astem5 hits");
                        too_many = 1;
                    } else {
                        (*hit.offset(i as isize)).pos = si.offset(-(n3 as isize));
                        (*hit.offset(i as isize)).energy = energy;
                        i += 1;
                    }
                }
            }
        }
    }
    i
}

/// find_resume_seq - Find tmRNA resume sequence
/// trna.c:152-257
pub unsafe fn find_resume_seq(
    mut s: *mut i32,
    ls: i32,
    hit: *mut trna_loop,
    nh: i32,
    sw: *mut csw
) -> i32 {
    let mut e: i32;
    let mut i: i32;
    let mut j: i32;
    let mut k: i32;
    let mut a: i32;
    let mut aa: [i32; 3] = [0; 3];
    let mut si: *mut i32;
    let mut sb: *mut i32;
    let mut sf: *mut i32;
    let mut st: *mut i32;
    let mut sl: *mut i32;
    let mut too_many: i32;
    let mut found: i32;
    let mut al: f64;
    let mut r: u32;
    let mut c: u32;
    let mut thresh: u32;

    static nps: [i32; 105] = [
        0,0,0,0, 0,0,0,0,
        0,0,0,0, 1,1,1,1,
        0,0,0,0, 1,1,1,1,
        0,0,0,0, 1,1,1,1,
        0,0,0,0, 1,1,1,1,
        1,1,1,1, 1,1,1,1,
        0,1,0,1, 0,0,0,0,
        0,1,1,1, 1,1,1,1,
        0,0,0,0, 0,0,0,0,
        0,0,0,0, 0,0,0,0,
        0,0,0,0, 0,0,0,0,
        0,0,0,0, 0,0,0,0,
        0,0,0,0, 0,0,0,0,0
    ];
    static score: [f64; 4] = [36.0, 66.0, 62.0, 72.0];
    static tem: [u32; 6] = [0x10310000, 0x01000101, 0x00010030, 0x02000100, 0x00000000, 0x00000000];
    static A: [i32; 6] = [0, 1, 1, 1, 1, 1];
    static V: [i32; 6] = [0, 0, 0, 1, 1, 1];
    static M: [i32; 6] = [0, 0, 1, 1, 1, 1];

    thresh = (*sw).tmrthresh as u32;
    i = 0;
    too_many = 0;
    sl = s.offset(ls as isize);
    r = *tem.get_unchecked(*s as usize);
    s = s.offset(1);
    r = (r >> 4) + *tem.get_unchecked(*s as usize);
    s = s.offset(1);
    r = (r >> 4) + *tem.get_unchecked(*s as usize);
    s = s.offset(1);
    r = (r >> 4) + *tem.get_unchecked(*s as usize);
    s = s.offset(1);
    r = (r >> 4) + *tem.get_unchecked(*s as usize);
    s = s.offset(1);
    r = (r >> 4) + *tem.get_unchecked(*s as usize);
    s = s.offset(1);
    r = (r >> 4) + *tem.get_unchecked(*s as usize);
    s = s.offset(1);

    if (*sw).tmstrict != 0 {
        while s < sl && too_many == 0 {
            r = (r >> 4) + *tem.get_unchecked(*s as usize);
            s = s.offset(1);
            c = r & 0xF;
            if c < thresh {
                continue;
            }
            c = c.wrapping_sub((*V.get_unchecked(*s.offset(1) as usize) + *V.get_unchecked(*s.offset(2) as usize)
                + *M.get_unchecked(*s.offset(5) as usize) + *A.get_unchecked(*s.offset(8) as usize)) as u32);
            if c < thresh {
                continue;
            }
            if i >= nh {
                too_many = 1;
                continue;
            }
            st = s.offset(-2);
            si = st;
            sb = st.offset((MINTAGDIST + 2) as isize);
            sf = st.offset(MAXTAGDIST as isize);
            found = 0;
            // C: while (si < sf) { if (*si++ != Thymine) si++; else if ... si++; }
            // The key insight: C uses *si++ which ALWAYS increments, even when checking Thymine
            while si < sf && found == 0 {
                if *si != Thymine {
                    // C: *si++ (increment) then si++ = total +2, then si++ at end = +3
                    si = si.offset(1);  // simulate *si++
                    si = si.offset(1);  // C: si++ in the if body
                } else if *si.offset(1) == Adenine {
                    // C: *si++ incremented, so now at position+1, check Adenine
                    // then *++si increments to position+2 before checking & 5
                    si = si.offset(2);  // C: *si++ then *++si = +2
                    if (*si & 5) == 0 {
                        found = 1;
                    }
                } else if *si.offset(1) == Guanine {
                    si = si.offset(2);  // C: *si++ then *++si = +2
                    if *si == Adenine {
                        found = 1;
                    }
                } else {
                    // C: *si++ incremented, Thymine found but not A or G
                    si = si.offset(2);  // C: *si++ then else si++
                }
                if found == 0 {
                    si = si.offset(1);  // C: si++ at end of while body
                }
            }
            if found == 0 {
                continue;
            }
            if si < sb {
                continue;
            }
            al = 0.0;
            k = 0;
            j = -11;
            while j < -2 {
                a = *si.offset(j as isize);
                j += 1;
                a = (a << 2) | *si.offset(j as isize);
                j += 1;
                if a == 9 {
                    al = (11 + 2 * ((j + 9) / 3)) as f64;
                }
                a = (a << 2) | *si.offset(j as isize);
                j += 1;
                *aa.get_unchecked_mut(k as usize) = a;
                k += 1;
            }
            (*hit.offset(i as isize)).pos = st;
            (*hit.offset(i as isize)).stem = (si as isize - st as isize) as i32 / std::mem::size_of::<i32>() as i32;
            e = (*nps.get_unchecked(*aa.get_unchecked(1) as usize) << 1) | *nps.get_unchecked(*aa.get_unchecked(2) as usize);
            (*hit.offset(i as isize)).energy = ((c << 2) as f64) + score[e as usize] + al
                + find_taghairpin(si)
                + find_tag_upstream_hairpin(st.offset(-10));
            i += 1;
        }
    } else {
        while s < sl && too_many == 0 {
            r = (r >> 4) + *tem.get_unchecked(*s as usize);
            s = s.offset(1);
            c = r & 0xF;
            if c < thresh {
                continue;
            }
            if i >= nh {
                too_many = 1;
                continue;
            }
            st = s.offset(-2);
            si = st.offset(MINTAGDIST as isize);
            sf = st.offset(MAXTAGDIST as isize);
            found = 0;
            // C: while (si < sf) { if (*si++ != Thymine) si++; else if ... si++; }
            while si < sf && found == 0 {
                if *si != Thymine {
                    si = si.offset(1);  // simulate *si++
                    si = si.offset(1);  // C: si++ in the if body
                } else if *si.offset(1) == Adenine {
                    si = si.offset(2);  // C: *si++ then *++si = +2
                    if (*si & 5) == 0 {
                        found = 1;
                    }
                } else if *si.offset(1) == Guanine {
                    si = si.offset(2);  // C: *si++ then *++si = +2
                    if *si == Adenine {
                        found = 1;
                    }
                } else {
                    si = si.offset(2);  // C: *si++ then else si++
                }
                if found == 0 {
                    si = si.offset(1);  // C: si++ at end of while body
                }
            }
            if found == 0 {
                continue;
            }
            (*hit.offset(i as isize)).pos = st;
            (*hit.offset(i as isize)).stem = (si as isize - st as isize) as i32 / std::mem::size_of::<i32>() as i32;
            e = (*nps.get_unchecked(((*si.offset(-8) << 4) | (*si.offset(-7) << 2) | *si.offset(-6)) as usize) << 1)
              | *nps.get_unchecked(((*si.offset(-5) << 4) | (*si.offset(-4) << 2) | *si.offset(-3)) as usize);
            (*hit.offset(i as isize)).energy = 46.0 + ((c << 2) as f64) + score[e as usize];
            i += 1;
        }
    }
    if too_many != 0 {
        eprintln!("Too many resume sequence hits");
    }
    i
}

/// base_copy3 - Copy bases and terminate
/// trna.c:260-263
#[inline(always)]
pub unsafe fn base_copy3(from: *mut i32, to: *mut i32, mut n: i32) -> *mut i32 {
    let mut from_ptr = from;
    let mut to_ptr = to;
    while n > 0 {
        *to_ptr = *from_ptr;
        to_ptr = to_ptr.offset(1);
        from_ptr = from_ptr.offset(1);
        n -= 1;
    }
    *to_ptr = TERM;
    to_ptr
}

/// remove_intron - Remove intron from sequence
/// trna.c:266-274
pub unsafe fn remove_intron(
    mut s1: *mut i32,
    mut s2: *mut i32,
    mut nbase: i32,
    intron: i32,
    nintron: i32
) {
    let mut s1e: *mut i32;
    s1e = s1.offset(intron as isize);
    nbase -= intron;
    while s1 < s1e {
        *s2 = *s1;
        s2 = s2.offset(1);
        s1 = s1.offset(1);
    }
    s1 = s1.offset(nintron as isize);
    s1e = s1.offset(nbase as isize);
    while s1 < s1e {
        *s2 = *s1;
        s2 = s2.offset(1);
        s1 = s1.offset(1);
    }
    *s2 = TERM;
}

/// nearest_trna_gene - Find nearest overlapping tRNA gene
/// trna.c:277-362
/// nearest_trna_gene_ts - Find nearest overlapping tRNA gene (with explicit ts parameter)
pub unsafe fn nearest_trna_gene_ts(
    ts: *mut Gene,
    d: *mut data_set,
    nt: i32,
    t: *mut Gene,
    sw: *mut csw
) -> *mut Gene {
    let mut n: i32;
    let mut i: i32;
    let mut comp: i32;
    let mut mtrna: i32;
    let mut mtcompov: i32;
    let mut maxintronlen: i32;
    let mut ilength: i64;
    let mut a: i64;
    let mut b: i64;
    let mut c: i64;
    let mut e: i64;
    let mut score: i64;
    let mut thresh: i64;
    let mut psmax: i64;
    static proximity: i64 = 7 * MINCTRNALEN as i64 / 10;
    let mut energy: f64;

    psmax = (*d).psmax;
    comp = (*t).comp;
    mtrna = (*sw).mtrna;
    mtcompov = (*sw).mtcompov;
    maxintronlen = (*sw).maxintronlen;
    n = -1;
    energy = INACTIVE;
    a = (*t).start;
    b = (*t).stop;
    thresh = b - a;

    if b < a {
        b += psmax;
        thresh += psmax;
        for i in 0..nt {
            c = (*ts.offset(i as isize)).start;
            e = (*ts.offset(i as isize)).stop;
            if e < c {
                e += psmax;
                if !(a > e) && !(b < c) {
                    if (*ts.offset(i as isize)).genetype == tRNA {
                        if (*ts.offset(i as isize)).comp == comp || (mtrna != 0 && mtcompov == 0) {
                            if maxintronlen <= 0
                               || ((2 * thresh) <= (5 * (e - c)) && (2 * (e - c)) <= (5 * thresh))
                            {
                                score = if a >= c {
                                    if b >= e { e - a } else { thresh }
                                } else {
                                    if b >= e { e - c } else { b - c }
                                };
                                if score >= proximity {
                                    if (*ts.offset(i as isize)).energy < energy {
                                        n = i;
                                        energy = (*ts.offset(i as isize)).energy;
                                    }
                                }
                            }
                        }
                    }
                }
                c -= psmax;
                e -= psmax;
            }
            if a > e { continue; }
            if b < c { continue; }
            if (*ts.offset(i as isize)).genetype != tRNA { continue; }
            if (*ts.offset(i as isize)).comp != comp {
                if mtrna == 0 { continue; }
                if mtcompov != 0 { continue; }
            }
            if maxintronlen > 0 {
                ilength = e - c;
                if (2 * thresh) > (5 * ilength) { continue; }
                if (2 * ilength) > (5 * thresh) { continue; }
            }
            score = if a >= c {
                if b >= e { e - a } else { thresh }
            } else {
                if b >= e { e - c } else { b - c }
            };
            if score >= proximity {
                if (*ts.offset(i as isize)).energy < energy {
                    n = i;
                    energy = (*ts.offset(i as isize)).energy;
                }
            }
        }
        a -= psmax;
        b -= psmax;
    }

    for i in 0..nt {
        c = (*ts.offset(i as isize)).start;
        e = (*ts.offset(i as isize)).stop;
        if e < c {
            e += psmax;
            if !(a > e) && !(b < c) {
                if (*ts.offset(i as isize)).genetype == tRNA {
                    if (*ts.offset(i as isize)).comp == comp || (mtrna != 0 && mtcompov == 0) {
                        if maxintronlen <= 0
                           || ((2 * thresh) <= (5 * (e - c)) && (2 * (e - c)) <= (5 * thresh))
                        {
                            score = if a >= c {
                                if b >= e { e - a } else { thresh }
                            } else {
                                if b >= e { e - c } else { b - c }
                            };
                            if score >= proximity {
                                if (*ts.offset(i as isize)).energy < energy {
                                    n = i;
                                    energy = (*ts.offset(i as isize)).energy;
                                }
                            }
                        }
                    }
                }
            }
            c -= psmax;
            e -= psmax;
        }
        if a > e { continue; }
        if b < c { continue; }
        if (*ts.offset(i as isize)).genetype != tRNA { continue; }
        if (*ts.offset(i as isize)).comp != comp {
            if mtrna == 0 { continue; }
            if mtcompov != 0 { continue; }
        }
        if maxintronlen > 0 {
            ilength = e - c;
            if (2 * thresh) > (5 * ilength) { continue; }
            if (2 * ilength) > (5 * thresh) { continue; }
        }
        score = if a >= c {
            if b >= e { e - a } else { thresh }
        } else {
            if b >= e { e - c } else { b - c }
        };
        if score >= proximity {
            if (*ts.offset(i as isize)).energy < energy {
                n = i;
                energy = (*ts.offset(i as isize)).energy;
            }
        }
    }

    if n >= 0 {
        return ts.offset(n as isize);
    }
    std::ptr::null_mut()
}

/// nearest_trna_gene - Find nearest overlapping tRNA gene (legacy wrapper using global TS)
pub unsafe fn nearest_trna_gene(
    d: *mut data_set,
    nt: i32,
    t: *mut Gene,
    sw: *mut csw
) -> *mut Gene {
    nearest_trna_gene_ts(TS, d, nt, t, sw)
}

/// nearest_tmrna_gene_ts - Find nearest overlapping tmRNA gene (with explicit ts parameter)
pub unsafe fn nearest_tmrna_gene_ts(
    ts: *mut Gene,
    d: *mut data_set,
    nt: i32,
    t: *mut Gene
) -> *mut Gene {
    let mut n: i32;
    let mut i: i32;
    let mut comp: i32;
    let mut a: i64;
    let mut b: i64;
    let mut c: i64;
    let mut e: i64;
    let mut score: i64;
    let mut smax: i64;
    let mut thresh: i64;
    let mut psmax: i64;

    psmax = (*d).psmax;
    comp = (*t).comp;
    smax = -1;
    n = -1;
    a = (*t).start;
    b = (*t).stop;
    thresh = b - a;

    if b < a {
        b += psmax;
        thresh += psmax;
        for i in 0..nt {
            c = (*ts.offset(i as isize)).start;
            e = (*ts.offset(i as isize)).stop;
            if e < c {
                e += psmax;
                if !(a > e) && !(b < c) {
                    if (*ts.offset(i as isize)).genetype == tmRNA && (*ts.offset(i as isize)).comp == comp {
                        score = if a >= c {
                            if b >= e { e - a } else { thresh }
                        } else {
                            if b >= e { e - c } else { b - c }
                        };
                        if score >= smax {
                            if score > smax {
                                n = i;
                                smax = score;
                            } else if (*ts.offset(i as isize)).energy < (*ts.offset(n as isize)).energy {
                                n = i;
                            }
                        }
                    }
                }
                c -= psmax;
                e -= psmax;
            }
            if a > e { continue; }
            if b < c { continue; }
            if (*ts.offset(i as isize)).genetype != tmRNA { continue; }
            if (*ts.offset(i as isize)).comp != comp { continue; }
            score = if a >= c {
                if b >= e { e - a } else { thresh }
            } else {
                if b >= e { e - c } else { b - c }
            };
            if score >= smax {
                if score > smax {
                    n = i;
                    smax = score;
                } else if (*ts.offset(i as isize)).energy < (*ts.offset(n as isize)).energy {
                    n = i;
                }
            }
        }
        a -= psmax;
        b -= psmax;
    }

    for i in 0..nt {
        c = (*ts.offset(i as isize)).start;
        e = (*ts.offset(i as isize)).stop;
        if e < c {
            e += psmax;
            if !(a > e) && !(b < c) {
                if (*ts.offset(i as isize)).genetype == tmRNA && (*ts.offset(i as isize)).comp == comp {
                    score = if a >= c {
                        if b >= e { e - a } else { thresh }
                    } else {
                        if b >= e { e - c } else { b - c }
                    };
                    if score >= smax {
                        if score > smax {
                            n = i;
                            smax = score;
                        } else if (*ts.offset(i as isize)).energy < (*ts.offset(n as isize)).energy {
                            n = i;
                        }
                    }
                }
            }
            c -= psmax;
            e -= psmax;
        }
        if a > e { continue; }
        if b < c { continue; }
        if (*ts.offset(i as isize)).genetype != tmRNA { continue; }
        if (*ts.offset(i as isize)).comp != comp { continue; }
        score = if a >= c {
            if b >= e { e - a } else { thresh }
        } else {
            if b >= e { e - c } else { b - c }
        };
        if score >= smax {
            if score > smax {
                n = i;
                smax = score;
            } else if (*ts.offset(i as isize)).energy < (*ts.offset(n as isize)).energy {
                n = i;
            }
        }
    }

    if (10 * smax) > (9 * thresh) {
        return ts.offset(n as isize);
    }
    std::ptr::null_mut()
}

/// nearest_tmrna_gene - Find nearest overlapping tmRNA gene (legacy wrapper using global TS)
pub unsafe fn nearest_tmrna_gene(
    d: *mut data_set,
    nt: i32,
    t: *mut Gene
) -> *mut Gene {
    nearest_tmrna_gene_ts(TS, d, nt, t)
}

/// overlap_ts - Check and report gene overlaps (with explicit ts parameter)
pub unsafe fn overlap_ts(
    ts: *mut Gene,
    d: *mut data_set,
    sort: *mut i32,
    n: i32,
    it: i32,
    sw: *mut csw
) {
    let mut i: i32;
    let mut j: i32;
    let mut flag: i32;
    let mut cross: i32;
    let mut crosstoo: i32;
    let mut ov_found: i32;
    let mut a: i64;
    let mut b: i64;
    let mut e: i64;
    let mut f: i64;
    let mut a2: i64 = 0;
    let mut b2: i64 = 0;
    let mut e2: i64;
    let mut f2: i64;
    let mut psmax: i64;
    let mut sname: [u8; 100] = [0; 100];
    let mut s: [u8; 100] = [0; 100];

    flag = 0;
    cross = 0;
    psmax = (*d).psmax;
    a = (*ts.offset(it as isize)).start;
    b = (*ts.offset(it as isize)).stop;
    if b < a {
        a2 = a - psmax;
        b2 = b;
        b += psmax;
        cross = 1;
    }

    j = -1;
    j += 1;
    while j < n {
        i = *sort.offset(j as isize);
        if i == it {
            j += 1;
            continue;
        }
        e = (*ts.offset(i as isize)).start;
        f = (*ts.offset(i as isize)).stop;
        crosstoo = 0;
        e2 = 0;
        f2 = 0;
        if f < e {
            e2 = e - psmax;
            f2 = f;
            f += psmax;
            crosstoo = 1;
        }
        ov_found = 0;
        if a <= f && b >= e {
            ov_found = 1;
        } else if crosstoo != 0 && a <= f2 && b >= e2 {
            ov_found = 1;
        } else if cross != 0 {
            if a2 <= f && b2 >= e {
                ov_found = 1;
            } else if crosstoo != 0 && a2 <= f2 && b2 >= e2 {
                ov_found = 1;
            }
        }
        if ov_found != 0 {
            if flag == 0 {
                let _ = write!((*sw).f.as_mut().unwrap(), "\n");
            }
            name(ts.offset(i as isize), sname.as_mut_ptr(), 1, sw);
            location(s.as_mut_ptr(), ts.offset(i as isize), sw, sname.as_ptr());
            let s_str = std::ffi::CStr::from_ptr(s.as_ptr() as *const i8).to_str().unwrap_or("");
            let _ = writeln!((*sw).f.as_mut().unwrap(), "Overlap with {}: {}", j + 1, s_str);
            flag = 1;
        }
        j += 1;
    }
    if flag != 0 {
        let _ = write!((*sw).f.as_mut().unwrap(), "\n");
    }
}

/// overlap - Check and report gene overlaps (legacy wrapper using global TS)
pub unsafe fn overlap(
    d: *mut data_set,
    sort: *mut i32,
    n: i32,
    it: i32,
    sw: *mut csw
) {
    overlap_ts(TS, d, sort, n, it, sw);
}

/// init_gene - Initialize gene array
/// trna.c:487-493
/// init_gene - Initialize gene slots in storage
/// ts: Gene storage array (global TS or thread-local)
pub unsafe fn init_gene_ts(ts: *mut Gene, nstart: i32, nstop: i32) {
    for i in nstart..nstop {
        (*ts.offset(i as isize)).energy = -1.0;
        (*ts.offset(i as isize)).genetype = noGENE;
        (*ts.offset(i as isize)).tps = 0;
        (*ts.offset(i as isize)).name[0] = 0;
        (*ts.offset(i as isize)).eseq[0] = TERM; // Initialize eseq to prevent false "extended sequence" output
    }
}

/// init_gene - Initialize gene slots (legacy wrapper using global TS)
pub unsafe fn init_gene(nstart: i32, nstop: i32) {
    init_gene_ts(TS, nstart, nstop);
}

/// find_slot_ts - Find or allocate slot for gene (with explicit ts parameter)
/// Note: Thread-local version does NOT reallocate - returns null if out of space
/// max_genes: maximum genes in ts array (for bounds checking)
pub unsafe fn find_slot_ts(
    ts: *mut Gene,
    max_genes: i32,
    d: *mut data_set,
    t: *mut Gene,
    nts: *mut i32,
    sw: *mut csw
) -> *mut Gene {
    let mut s1: [u8; 80] = [0; 80];
    let mut s2: [u8; 80] = [0; 80];
    let mut s3: [u8; 80] = [0; 80];
    let mut s4: [u8; 80] = [0; 80];
    let mut tn: *mut Gene;

    if (*sw).comp != 0 {
        (*t).stop = (*sw).start - (*t).start - 1;
        (*t).start = (*t).stop - (*t).nbase as i64 - (*t).nintron as i64 + 1;
        (*t).comp = 1;
    } else {
        (*t).start += (*sw).start;
        (*t).stop = (*t).start + (*t).nbase as i64 + (*t).nintron as i64 - 1;
        (*t).comp = 0;
    }

    if (*sw).linear == 0 {
        (*t).start = sq((*t).start, (*d).psmax);
        (*t).stop = sq((*t).stop, (*d).psmax);
    }

    if (*t).genetype == tRNA {
        tn = nearest_trna_gene_ts(ts, d, *nts, t, sw);
    } else if (*t).genetype == tmRNA {
        tn = nearest_tmrna_gene_ts(ts, d, *nts, t);
    } else {
        tn = std::ptr::null_mut();
    }

    if !tn.is_null() {
        if (*t).energy <= (*tn).energy {
            return std::ptr::null_mut();
        }
        copy((*tn).name.as_mut_ptr(), (*t).name.as_mut_ptr());
        if (*sw).verbose != 0 {
            eprint!("{} {} ",
                std::ffi::CStr::from_ptr(name(t, s1.as_mut_ptr(), 0, sw) as *const i8).to_str().unwrap_or(""),
                std::ffi::CStr::from_ptr(position(s3.as_mut_ptr(), t, sw) as *const i8).to_str().unwrap_or(""));
            if (*sw).energydisp != 0 {
                eprint!("({}) ", nenergy(t, sw));
            }
            eprint!("replacing {} {}",
                std::ffi::CStr::from_ptr(name(tn, s2.as_mut_ptr(), 1, sw) as *const i8).to_str().unwrap_or(""),
                std::ffi::CStr::from_ptr(position(s4.as_mut_ptr(), tn, sw) as *const i8).to_str().unwrap_or(""));
            if (*sw).energydisp != 0 {
                eprint!(" ({})", nenergy(tn, sw));
            }
            eprintln!();
        }
    } else {
        // Thread-local version: no reallocation, just bounds check
        if *nts >= max_genes {
            eprintln!("Thread-local gene storage full ({} genes)", max_genes);
            eprintln!("Gene lost");
            return std::ptr::null_mut();
        }
        copy3cr((*d).seqname.as_mut_ptr(), (*t).name.as_mut_ptr(), 99);
        tn = ts.offset(*nts as isize);
        *nts = *nts + 1;
        if (*sw).verbose != 0 {
            eprint!("{} at {}",
                std::ffi::CStr::from_ptr(name(t, s1.as_mut_ptr(), 0, sw) as *const i8).to_str().unwrap_or(""),
                std::ffi::CStr::from_ptr(position(s2.as_mut_ptr(), t, sw) as *const i8).to_str().unwrap_or(""));
            if (*sw).energydisp != 0 {
                eprint!(" ({})", nenergy(t, sw));
            }
            eprintln!();
        }
    }
    tn
}

/// find_slot - Find or allocate slot for gene (legacy wrapper using global TS)
/// trna.c:496-550
pub unsafe fn find_slot(
    d: *mut data_set,
    t: *mut Gene,
    nts: *mut i32,
    sw: *mut csw
) -> *mut Gene {
    let mut newspace: i32;
    let mut s1: [u8; 80] = [0; 80];
    let mut s2: [u8; 80] = [0; 80];
    let mut s3: [u8; 80] = [0; 80];
    let mut s4: [u8; 80] = [0; 80];
    let mut tn: *mut Gene;
    let mut tsn: *mut Gene;

    if (*sw).comp != 0 {
        (*t).stop = (*sw).start - (*t).start - 1;
        (*t).start = (*t).stop - (*t).nbase as i64 - (*t).nintron as i64 + 1;
        (*t).comp = 1;
    } else {
        (*t).start += (*sw).start;
        (*t).stop = (*t).start + (*t).nbase as i64 + (*t).nintron as i64 - 1;
        (*t).comp = 0;
    }

    if (*sw).linear == 0 {
        (*t).start = sq((*t).start, (*d).psmax);
        (*t).stop = sq((*t).stop, (*d).psmax);
    }

    if (*t).genetype == tRNA {
        tn = nearest_trna_gene(d, *nts, t, sw);
    } else if (*t).genetype == tmRNA {
        tn = nearest_tmrna_gene(d, *nts, t);
    } else {
        tn = std::ptr::null_mut();
    }

    if !tn.is_null() {
        if (*t).energy <= (*tn).energy {
            return std::ptr::null_mut();
        }
        copy((*tn).name.as_mut_ptr(), (*t).name.as_mut_ptr());
        if (*sw).verbose != 0 {
            eprint!("{} {} ",
                std::ffi::CStr::from_ptr(name(t, s1.as_mut_ptr(), 0, sw) as *const i8).to_str().unwrap_or(""),
                std::ffi::CStr::from_ptr(position(s3.as_mut_ptr(), t, sw) as *const i8).to_str().unwrap_or(""));
            if (*sw).energydisp != 0 {
                eprint!("({}) ", nenergy(t, sw));
            }
            eprint!("replacing {} {}",
                std::ffi::CStr::from_ptr(name(tn, s2.as_mut_ptr(), 1, sw) as *const i8).to_str().unwrap_or(""),
                std::ffi::CStr::from_ptr(position(s4.as_mut_ptr(), tn, sw) as *const i8).to_str().unwrap_or(""));
            if (*sw).energydisp != 0 {
                eprint!(" ({})", nenergy(tn, sw));
            }
            eprintln!();
        }
    } else {
        if *nts >= (*sw).genespace {
            newspace = if (*d).ps > 0 {
                (*sw).genespace * (1 + ((*d).psmax / (*d).ps) as i32)
            } else {
                (*sw).genespace + NT as i32
            };
            tsn = std::alloc::realloc(
                TS as *mut u8,
                std::alloc::Layout::array::<Gene>((*sw).genespace as usize).unwrap(),
                newspace as usize * std::mem::size_of::<Gene>()
            ) as *mut Gene;
            if tsn.is_null() {
                eprintln!("No more memory to store detected genes");
                eprintln!("Gene lost");
                return std::ptr::null_mut();
            }
            if (*sw).verbose != 0 {
                eprintln!("Expanding detected gene store from {} genes to {} genes",
                    (*sw).genespace, newspace);
            }
            TS = tsn;
            init_gene((*sw).genespace, newspace);
            (*sw).genespace = newspace;
        }
        copy3cr((*d).seqname.as_mut_ptr(), (*t).name.as_mut_ptr(), 99);
        tn = TS.offset(*nts as isize);
        *nts = *nts + 1;
        if (*sw).verbose != 0 {
            eprint!("{} at {}",
                std::ffi::CStr::from_ptr(name(t, s1.as_mut_ptr(), 0, sw) as *const i8).to_str().unwrap_or(""),
                std::ffi::CStr::from_ptr(position(s2.as_mut_ptr(), t, sw) as *const i8).to_str().unwrap_or(""));
            if (*sw).energydisp != 0 {
                eprint!(" ({})", nenergy(t, sw));
            }
            eprintln!();
        }
    }
    tn
}

/// aatail - Score acceptor stem CCA tail
/// trna.c:553-586
pub unsafe fn aatail(s: *mut i32, ext: *mut i32, sw: *mut csw) -> i32 {
    let mut score: i32;
    let mut e: i32;
    static A: [i32; 6] = [1, 0, 0, 0, 0, 0];
    static C: [i32; 6] = [0, 1, 0, 0, 0, 0];

    if (*sw).aataildiv != 0 {
        score = 0;
        e = 0;
        if *A.get_unchecked(*s.offset(3) as usize) != 0 {
            score += 1;
            e = 3;
        }
        if *C.get_unchecked(*s.offset(2) as usize) != 0 {
            score += 1;
            if e == 0 { e = 2; }
        }
        if *C.get_unchecked(*s.offset(1) as usize) != 0 {
            score += 1;
            if e == 0 { e = 1; }
        }
        if score < 2 {
            if *A.get_unchecked(*s as usize) != 0 {
                score += 1;
            }
        }
        e += 1;
        *ext = e;
        score
    } else {
        score = 1;
        e = 1;
        if *C.get_unchecked(*s.offset(1) as usize) != 0 {
            score += 1;
            e = 2;
            if *C.get_unchecked(*s.offset(2) as usize) != 0 {
                score += 1;
                e = 3;
                if *A.get_unchecked(*s.offset(3) as usize) != 0 {
                    score += 1;
                    e = 4;
                }
            }
        }
        *ext = e;
        score
    }
}

