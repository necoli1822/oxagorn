/*
 * Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
 * main.rs - Main program entry point and command-line processing
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

use std::env;
use std::ptr;
use std::io::Write;

use rayon;
use oxagorn::types::*;
use oxagorn::tables::*;
use oxagorn::utils::*;
use oxagorn::output::*;
use oxagorn::output::parallel_fastafile;

// Constants imported from types.rs via `use oxagorn::types::*`
// All constants (NHELPLINE, NGENECODE, MAXGCMOD, NAMINOACID, NT, STANDARD, METAZOAN_MT,
// VERTEBRATE_MT, MAMMAL_MT, LSEQ, MAXTRNALEN, MAXINTRONLEN, Adenine, Cytosine, Guanine,
// Thymine, TERM, tRNAthresh, mtRNAtthresh, mtRNAdthresh, mtRNAdtthresh, tmRNAthresh,
// srpRNAthresh, CDSthresh, PSEUDOGENElevel, NTAG) are defined in types.rs

// Use tables from library (no extern "C" needed)
// HELPMENU, AANAME, AAMAP are from tables.rs
// TS is the global gene storage from tables.rs

// Helper to get C stdout (platform-specific)
#[inline]
unsafe fn c_stdout() -> *mut libc::FILE {
    #[cfg(target_os = "macos")]
    {
        extern "C" {
            static __stdoutp: *mut libc::FILE;
        }
        __stdoutp
    }
    #[cfg(not(target_os = "macos"))]
    {
        extern "C" {
            static stdout: *mut libc::FILE;
        }
        stdout
    }
}

// Helper to get C stderr
#[inline]
unsafe fn c_stderr() -> *mut libc::FILE {
    extern "C" {
        static stderr: *mut libc::FILE;
    }
    stderr
}

/// oxagorn_help_menu - Display help menu
/// main.c:11-15
fn oxagorn_help_menu() {
    for h in 0..NHELPLINE {
        println!("{}", HELPMENU[h]);
    }
}

/// error_report - Report errors and exit
/// main.c:17-36
unsafe fn error_report(n: i32, s: *const u8) {
    match n {
        0 => {
            let msg = std::ffi::CStr::from_ptr(s as *const i8).to_string_lossy();
            eprintln!("-{} not recognised, type oxagorn -h for help", msg);
        }
        1 => {
            let msg = std::ffi::CStr::from_ptr(s as *const i8).to_string_lossy();
            eprintln!("-{} not understood, type oxagorn -h for help", msg);
        }
        2 => {
            let msg = std::ffi::CStr::from_ptr(s as *const i8).to_string_lossy();
            eprintln!("Could not open {}", msg);
        }
        3 => {
            eprintln!("No sequence file specified, type oxagorn -h for help");
        }
        4 => {
            let msg = std::ffi::CStr::from_ptr(s as *const i8).to_string_lossy();
            eprintln!("Don't know genetic code {}", msg);
        }
        5 => {
            eprintln!("Too many genetic code modifications (max={})", MAXGCMOD);
        }
        _ => {}
    }
    std::process::exit(0);
}

/// process_genecode_switch - Process genetic code command-line switch
/// main.c:39-97
unsafe fn process_genecode_switch(mut s: *mut u8, sw: *mut csw) {
    let mut i: i32;
    let mut m: i32;
    let mut lmax: i32;
    let mut len: [i32; NGENECODE] = [0; NGENECODE];
    let mut anticodon: i32;
    let mut b: [i32; 3] = [0; 3];
    let mut l: i64;
    let mut c: u8;
    let mut ss: *mut u8;
    let mut se: *mut u8;

    static genecodetag: [[u8; 10]; NGENECODE] = [
        *b"MET\0\0\0\0\0\0\0",
        *b"STD\0\0\0\0\0\0\0", *b"VERT\0\0\0\0\0\0", *b"YEAST\0\0\0\0\0", *b"PROT\0\0\0\0\0\0", *b"INVERT\0\0\0\0",
        *b"CILIATE\0\0\0", *b"DELETED\0\0\0", *b"DELETED\0\0\0", *b"FLATWORM\0\0", *b"EUPLOT\0\0\0\0",
        *b"BACT\0\0\0\0\0\0", *b"ALTYEAST\0\0", *b"ASCID\0\0\0\0\0", *b"ALTFLAT\0\0\0", *b"BLEP\0\0\0\0\0\0",
        *b"CHLOROPH\0\0", *b"DELETED\0\0\0", *b"DELETED\0\0\0", *b"DELETED\0\0\0", *b"DELETED\0\0\0",
        *b"TREM\0\0\0\0\0\0", *b"SCEN\0\0\0\0\0\0", *b"THRAUST\0\0\0", *b"PTERO\0\0\0\0\0", *b"GRAC\0\0\0\0\0\0",
        *b"PACH\0\0\0\0\0\0", *b"KARY\0\0\0\0\0\0", *b"COND\0\0\0\0\0\0", *b"MESO\0\0\0\0\0\0", *b"PERI\0\0\0\0\0\0",
        *b"BLAST\0\0\0\0\0", *b"VACANT\0\0\0\0", *b"CEPH\0\0\0\0\0\0"
    ];

    (*sw).geneticcode = STANDARD;
    (*sw).gcfix = 1;
    c = *s;

    if c >= b'0' && c <= b'9' {
        l = 0;
        lconvert(s, &mut l);
        i = l as i32;
        if (i >= 0) && (i < NGENECODE as i32) {
            (*sw).geneticcode = i;
        }
    } else {
        for i in 0..NGENECODE {
            len[i] = 0;
            ss = s;
            se = genecodetag[i].as_ptr() as *mut u8;
            c = *ss;
            ss = ss.offset(1);
            while c == *ss.offset(-1) {
                if upcasec(c) as u8 != *se {
                    break;
                }
                se = se.offset(1);
                len[i] += 1;
                c = *ss;
                ss = ss.offset(1);
            }
        }
        m = -1;
        lmax = 0;
        i = -1;
        i += 1;
        while i < NGENECODE as i32 {
            if len[i as usize] > lmax {
                m = i;
                lmax = len[i as usize];
            }
            i += 1;
        }
        if m >= 0 {
            (*sw).geneticcode = m;
        } else {
            error_report(4, s);
        }
    }

    (*sw).ngcmod = 0;
    ss = s;

    // Process genetic code modifications (comma-separated)
    loop {
        ss = strpos(ss, b",\0".as_ptr() as *mut u8);
        if ss.is_null() {
            break;
        }
        if (*sw).ngcmod >= MAXGCMOD as i32 {
            error_report(5, ptr::null());
        }
        ss = ss.offset(1);

        for i in 0..3 {
            b[i as usize] = Adenine;
            c = upcasec(*ss.offset(i as isize)) as u8;
            if c == b'C' { b[i as usize] = Cytosine; }
            if c == b'G' { b[i as usize] = Guanine; }
            if c == b'T' { b[i as usize] = Thymine; }
            if c == b'U' { b[i as usize] = Thymine; }
        }

        anticodon = ((Thymine - b[2]) << 4) + ((Thymine - b[1]) << 2) + (Thymine - b[0]);

        se = strpos(ss, b"=\0".as_ptr() as *mut u8);
        if se.is_null() {
            break;
        }
        se = se.offset(1);

        for i in 0..NAMINOACID as i32 {
            if upcasec(*se) == upcasec(AANAME[i as usize][0]) {
                if upcasec(*se.offset(1)) == upcasec(AANAME[i as usize][1]) {
                    if upcasec(*se.offset(2)) == upcasec(AANAME[i as usize][2]) {
                        AAMAP[(*sw).geneticcode as usize][anticodon as usize] = i;
                        (*sw).gcmod[(*sw).ngcmod as usize] = anticodon;
                        (*sw).ngcmod += 1;
                        break;
                    }
                }
            }
        }
    }
}

/// change_thresholds - Modify detection thresholds
/// main.c:101-111
unsafe fn change_thresholds(sw: *mut csw, psthresh: f64) {
    (*sw).threshlevel = psthresh;
    (*sw).cdsthresh *= psthresh;
    (*sw).srpthresh *= psthresh;
    (*sw).tmrnathresh *= psthresh;
    (*sw).mtdtthresh *= psthresh;
    (*sw).mttthresh *= psthresh;
    (*sw).mtdthresh *= psthresh;
    (*sw).trnathresh *= psthresh;
}

/// main - Program entry point
/// main.c:114-387
fn main() {
    unsafe {
        let args: Vec<String> = env::args().collect();
        let z = args.len() as i32;

        let mut i: i32;
        let mut lv: i32;
        let mut filecounter: i32;
        let mut l: i64 = 0;
        let mut psthresh: f64;
        let mut num_threads: usize = 0;  // 0 = auto-detect
        let mut input_file_path: String = String::new();  // For parallel processing
        let mut c1: u8;
        let mut c2: u8;
        let mut c3: u8;
        let mut c4: u8;
        let mut s: *mut u8;

        // Initialize data_set
        let mut d: data_set = Default::default();
        d.bugmode = 0;

        // Initialize control switches with defaults
        let mut sw: csw = std::mem::zeroed();

        // Set default gene type names (10 bytes each)
        sw.genetypename[0] = *b"tRNA\0\0\0\0\0\0";
        sw.genetypename[1] = *b"tmRNA\0\0\0\0\0";
        sw.genetypename[2] = *b"\0\0\0\0\0\0\0\0\0\0";
        sw.genetypename[3] = *b"\0\0\0\0\0\0\0\0\0\0";
        sw.genetypename[4] = *b"CDS\0\0\0\0\0\0\0";
        sw.genetypename[5] = *b"overall\0\0\0";

        // Initialize default values (from static csw initialization in C)
        sw.f = c_stdout() as *mut File;
        sw.both = 2;
        sw.geneticcode = STANDARD;
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

        // Initialize tmrna_struct (ARAGORN signature)
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
            45, 45, 10, TERM, 0, 0, 0, 0, 0, 0
        ];
        for i in 0..100 {
            sw.tmrna_struct[i] = sig[i];
        }

        filecounter = 0;
        i = 0;

        // Parse command-line arguments
        i += 1;
        while i < z {
            let arg = &args[i as usize];
            let arg_bytes = arg.as_bytes();

            if arg_bytes[0] == b'-' {
                lv = arg_bytes.len() as i32;
                if lv < 2 {
                    i += 1;
                    continue;
                }

                s = arg_bytes.as_ptr().offset(1) as *mut u8;
                c1 = upcasec(*s);
                c2 = if lv > 2 { upcasec(*s.offset(1)) } else { b' ' };
                c3 = if lv > 3 { upcasec(*s.offset(2)) } else { b' ' };
                c4 = if lv > 4 { upcasec(*s.offset(3)) } else { b' ' };

                match c1 {
                    b'E' => {
                        sw.energydisp = if c2 == b'S' { 2 } else { 1 };
                    }
                    b'A' => {
                        if c2 == b'7' {
                            sw.extastem = 0;
                        } else if c2 == b'A' {
                            sw.matchacceptor = 1;
                        } else if c2 == b'M' {
                            // Handle -AM[T|M]<n> options
                            l = 1;
                            if c3 == b'T' {
                                if lv > 4 {
                                    lconvert(s.offset(3), &mut l);
                                    if l < 1 { l = 1; }
                                    sw.trnalenmisthresh = l as i32;
                                } else {
                                    sw.trnalenmisthresh = 1;
                                }
                            } else if c3 == b'M' {
                                if lv > 4 {
                                    lconvert(s.offset(3), &mut l);
                                    if l < 1 { l = 1; }
                                    sw.tmrnalenmisthresh = l as i32;
                                } else {
                                    sw.tmrnalenmisthresh = 1;
                                }
                            } else if lv > 3 {
                                lconvert(s.offset(2), &mut l);
                                if l < 1 { l = 1; }
                                sw.trnalenmisthresh = l as i32;
                                sw.tmrnalenmisthresh = l as i32;
                            } else {
                                sw.trnalenmisthresh = 1;
                                sw.tmrnalenmisthresh = 1;
                            }
                        } else {
                            sw.secstructdisp |= 1;
                        }
                    }
                    b'B' => {
                        // -BR: secondary structure display
                        // -B: library flag (C original compatibility)
                        if c2 == b'R' {
                            sw.secstructdisp |= 2;
                        } else {
                            sw.libflag = 1;
                        }
                    }
                    b'X' => {
                        sw.libflag = 2;
                    }
                    b'W' => {
                        if c2 == b'U' && c3 == b'N' && c4 == b'I' {
                            d.bugmode = 1;
                        } else {
                            if sw.batch < 1 {
                                sw.batch = 1;
                            }
                            if c2 == b'A' {
                                sw.batchfullspecies = 1;
                            }
                        }
                    }
                    b'V' => {
                        sw.verbose = 1;
                    }
                    b'S' => {
                        if c2 == b'S' {
                            sw.sp1max = 2;
                            sw.sp2min = 1;
                            sw.sp2max = 1;
                        } else if c2 == b'E' {
                            if sw.seqdisp < 1 {
                                sw.seqdisp = 1;
                            }
                        } else if c2 == b'C' || c2 == b'-' {
                            sw.both = 1;
                        } else if c2 == b'V' && c3 == b'G' {
                            sw.secstructdisp |= 4;
                        } else {
                            sw.both = 0;
                        }
                    }
                    b'F' => {
                        if !softstrpos(s, b"O\0".as_ptr() as *mut u8).is_null() {
                            sw.batch = 2;
                            if !softstrpos(s, b"S\0".as_ptr() as *mut u8).is_null() {
                                sw.batch |= 0x4;
                            }
                            if !softstrpos(s, b"N\0".as_ptr() as *mut u8).is_null() {
                                sw.batch |= 0x8;
                            }
                            if !softstrpos(s, b"C\0".as_ptr() as *mut u8).is_null() {
                                sw.batch |= 0x10;
                            }
                        } else {
                            if !softstrpos(s, b"C\0".as_ptr() as *mut u8).is_null() {
                                sw.seqdisp = 4;
                            } else {
                                sw.seqdisp = 3;
                            }
                        }
                    }
                    b'D' => {
                        sw.both = 2;
                    }
                    b'L' => {
                        sw.linear = 1;
                    }
                    b'C' => {
                        if c2 == b'7' {
                            sw.cloop7 = 1;
                        } else {
                            sw.linear = 0;
                        }
                    }
                    b'J' => {
                        if lv > 2 {
                            if c2 == b'R' {
                                sw.aataildiv = 1;
                            }
                            if c3 == b'4' {
                                sw.aataildisp = 1;
                            }
                        } else {
                            sw.aataildisp = 1;
                        }
                    }
                    b'1' => {
                        sw.minintronlen = 10;
                    }
                    b'I' => {
                        // Handle -I[O|F|R]<min>,<max> options
                        let mut s_ptr = s;
                        let mut lv_local = lv;
                        if c2 == b'O' { sw.ioverlay = 1; s_ptr = s_ptr.offset(1); lv_local -= 1; }
                        else if c2 == b'F' { sw.ifixedpos = 1; s_ptr = s_ptr.offset(1); lv_local -= 1; }
                        else if c2 == b'R' { sw.ireportminintronlen = 1; s_ptr = s_ptr.offset(1); lv_local -= 1; }
                        let c3_local = if lv_local > 2 { upcasec(*s_ptr.offset(1)) } else { b' ' };
                        if c3_local == b'O' { sw.ioverlay = 1; s_ptr = s_ptr.offset(1); lv_local -= 1; }
                        else if c3_local == b'F' { sw.ifixedpos = 1; s_ptr = s_ptr.offset(1); lv_local -= 1; }
                        else if c3_local == b'R' { sw.ireportminintronlen = 1; s_ptr = s_ptr.offset(1); lv_local -= 1; }
                        let c4_local = if lv_local > 3 { upcasec(*s_ptr.offset(2)) } else { b' ' };
                        if c4_local == b'O' { sw.ioverlay = 1; s_ptr = s_ptr.offset(1); lv_local -= 1; }
                        else if c4_local == b'F' { sw.ifixedpos = 1; s_ptr = s_ptr.offset(1); lv_local -= 1; }
                        else if c4_local == b'R' { sw.ireportminintronlen = 1; s_ptr = s_ptr.offset(1); lv_local -= 1; }

                        if lv_local > 2 {
                            let mut s_num = lconvert(s_ptr.offset(1), &mut l);
                            if *s_num == b',' {
                                if sw.ireportminintronlen == 1 {
                                    sw.minintronlenreport = l as i32;
                                } else {
                                    sw.minintronlen = l as i32;
                                }
                                lconvert(s_num.offset(1), &mut l);
                                sw.maxintronlen = l as i32;
                            } else {
                                sw.maxintronlen = l as i32;
                            }
                            if sw.maxintronlen > (LSEQ as i32 - MAXTRNALEN as i32) {
                                sw.maxintronlen = LSEQ as i32 - MAXTRNALEN as i32;
                            }
                            if sw.maxintronlen > MAXINTRONLEN as i32 {
                                sw.maxintronlen = MAXINTRONLEN as i32;
                            }
                            if sw.minintronlen < 0 || sw.maxintronlen < sw.minintronlen {
                                error_report(1, arg_bytes.as_ptr());
                            }
                            if sw.minintronlenreport < 0 || sw.maxintronlen < sw.minintronlenreport {
                                error_report(1, arg_bytes.as_ptr());
                            }
                        } else {
                            sw.maxintronlen = MAXINTRONLEN as i32;
                        }
                    }
                    b'T' => {
                        if c2 == b'V' {
                            sw.tvloop = 0;
                        } else {
                            sw.trna = 1;
                            if lv > 2 {
                                let s_next = dconvert(s.offset(1), &mut sw.trnathresh);
                                if *s_next == b',' {
                                    dconvert(s_next.offset(1), &mut sw.ttarmthresh);
                                }
                            }
                        }
                    }
                    b'M' => {
                        if c2 == b'T' {
                            sw.mtrna = 1;
                            if sw.gcfix == 0 {
                                sw.geneticcode = METAZOAN_MT;
                            }
                            if lv > 3 {
                                let mut s_mt = s.offset(2);
                                let mut c3_mt = upcasec(*s_mt);
                                if c3_mt == b'M' {
                                    loop {
                                        s_mt = s_mt.offset(1);
                                        c3_mt = upcasec(*s_mt);
                                        if !((c3_mt == b'A') || (c3_mt == b'M') || (c3_mt == b'L')) {
                                            break;
                                        }
                                    }
                                    sw.tvloop = 0;
                                    sw.geneticcode = VERTEBRATE_MT;
                                    sw.discrim = MAMMAL_MT;
                                }
                                while c3_mt == b'X' || c3_mt == b'C' || c3_mt == b'D' {
                                    if c3_mt == b'X' { sw.mtxdetect = 0; }
                                    else if c3_mt == b'C' { sw.mtcdsscan = 0; }
                                    else if c3_mt == b'D' { sw.mtcompov = 1; }
                                    s_mt = s_mt.offset(1);
                                    c3_mt = upcasec(*s_mt);
                                }
                                if c3_mt != b'-' {
                                    if c3_mt != b'.' {
                                        if (c3_mt < b'0') || (c3_mt > b'9') {
                                            // break from M case - no threshold parsing
                                        } else {
                                            s_mt = dconvert(s_mt, &mut sw.mtdtthresh);
                                            if *s_mt == b',' { s_mt = dconvert(s_mt.offset(1), &mut sw.mttthresh); }
                                            if *s_mt == b',' { s_mt = dconvert(s_mt.offset(1), &mut sw.mtdthresh); }
                                            if *s_mt == b',' { s_mt = dconvert(s_mt.offset(1), &mut sw.mttarmthresh); }
                                            if *s_mt == b',' { dconvert(s_mt.offset(1), &mut sw.mtdarmthresh); }
                                        }
                                    } else {
                                        s_mt = dconvert(s_mt, &mut sw.mtdtthresh);
                                        if *s_mt == b',' { s_mt = dconvert(s_mt.offset(1), &mut sw.mttthresh); }
                                        if *s_mt == b',' { s_mt = dconvert(s_mt.offset(1), &mut sw.mtdthresh); }
                                        if *s_mt == b',' { s_mt = dconvert(s_mt.offset(1), &mut sw.mttarmthresh); }
                                        if *s_mt == b',' { dconvert(s_mt.offset(1), &mut sw.mtdarmthresh); }
                                    }
                                } else {
                                    s_mt = dconvert(s_mt, &mut sw.mtdtthresh);
                                    if *s_mt == b',' { s_mt = dconvert(s_mt.offset(1), &mut sw.mttthresh); }
                                    if *s_mt == b',' { s_mt = dconvert(s_mt.offset(1), &mut sw.mtdthresh); }
                                    if *s_mt == b',' { s_mt = dconvert(s_mt.offset(1), &mut sw.mttarmthresh); }
                                    if *s_mt == b',' { dconvert(s_mt.offset(1), &mut sw.mtdarmthresh); }
                                }
                            }
                        } else {
                            sw.tmrna = 1;
                            let mut s_tm = s;
                            let mut lv_tm = lv;
                            if c2 == b'U' {
                                if c3 == b'T' {
                                    sw.updatetmrnatags = 1;
                                    lv_tm -= 2;
                                    s_tm = s_tm.offset(2);
                                }
                            }
                            if lv_tm > 2 {
                                dconvert(s_tm.offset(1), &mut sw.tmrnathresh);
                            }
                        }
                    }
                    b'P' => {
                        if c2 == b'S' {
                            if c3 != b'-' {
                                if c3 != b'.' {
                                    if (c3 < b'0') || (c3 > b'9') {
                                        change_thresholds(&mut sw, PSEUDOGENElevel);
                                    } else {
                                        psthresh = 1.0;
                                        dconvert(s.offset(2), &mut psthresh);
                                        change_thresholds(&mut sw, psthresh);
                                    }
                                } else {
                                    psthresh = 1.0;
                                    dconvert(s.offset(2), &mut psthresh);
                                    change_thresholds(&mut sw, psthresh);
                                }
                            } else {
                                psthresh = 1.0;
                                dconvert(s.offset(2), &mut psthresh);
                                change_thresholds(&mut sw, psthresh);
                            }
                        }
                    }
                    b'G' => {
                        if c2 == b'C' {
                            process_genecode_switch(s.offset(2), &mut sw);
                        }
                    }
                    b'R' => {
                        if c2 == b'N' {
                            sw.repeatsn = 1;
                        } else if c2 == b'P' {
                            sw.reportpseudogenes = 1;
                            if lv > 3 {
                                dconvert(s.offset(2), &mut sw.reportpsthresh);
                            }
                        } else {
                            sw.tmstrict = 0;
                        }
                    }
                    b'Q' => {
                        sw.showconfig = 0;
                    }
                    b'@' => {
                        // -@<N> or -@ <N>: Set number of threads for parallel processing
                        // -@0 or -@ 0 or -@ : auto-detect (default)
                        // -@1 or -@ 1: single-threaded (disable parallelism)
                        // -@4 or -@ 4: use 4 threads
                        if lv > 2 {
                            // Number attached: -@1, -@4, etc.
                            lconvert(s.offset(1), &mut l);
                            num_threads = l as usize;
                        } else if i + 1 < z {
                            // Check if next arg is a number: -@ 1, -@ 4, etc.
                            let next_arg = &args[(i + 1) as usize];
                            if let Ok(n) = next_arg.parse::<usize>() {
                                num_threads = n;
                                i += 1;  // consume the number argument
                            } else {
                                num_threads = 0;  // auto-detect
                            }
                        } else {
                            num_threads = 0;  // auto-detect
                        }
                    }
                    b'H' => {
                        oxagorn_help_menu();
                        std::process::exit(0);
                    }
                    b'O' => {
                        // Output file handling
                        let output_file: &str;
                        if lv > 2 {
                            output_file = &arg[2..];
                        } else {
                            i += 1;
                            if i >= z {
                                i += 1;
                                continue;
                            }
                            output_file = &args[i as usize];
                        }
                        let c_str = std::ffi::CString::new(output_file).unwrap();
                        sw.f = libc::fopen(c_str.as_ptr(), b"w\0".as_ptr() as *const i8) as *mut File;
                        if sw.f.is_null() {
                            error_report(2, c_str.as_ptr() as *const u8);
                        }
                    }
                    _ => {
                        error_report(0, s);
                    }
                }
            } else {
                // Input/output file argument
                if filecounter < 1 {
                    input_file_path = arg.clone();  // Store for parallel processing
                    let c_str = std::ffi::CString::new(arg.as_str()).unwrap();
                    d.f = libc::fopen(c_str.as_ptr(), b"r\0".as_ptr() as *const i8) as *mut File;
                    if !d.f.is_null() {
                        filecounter += 1;
                    } else {
                        error_report(2, c_str.as_ptr() as *const u8);
                    }
                } else if filecounter < 2 {
                    let c_str = std::ffi::CString::new(arg.as_str()).unwrap();
                    sw.f = libc::fopen(c_str.as_ptr(), b"w\0".as_ptr() as *const i8) as *mut File;
                    if sw.f.is_null() {
                        error_report(2, c_str.as_ptr() as *const u8);
                    }
                    filecounter += 1;
                } else {
                    error_report(0, arg.as_ptr());
                }
            }
            i += 1;
        }

        // Check if input file was provided
        if filecounter < 1 {
            error_report(3, ptr::null());
        }

        // Set default modes if not specified
        if (sw.trna == 0) && (sw.tmrna == 0) {
            sw.trna = 1;
            sw.tmrna = 1;
        }
        if sw.mtrna != 0 {
            sw.trna = 0;
        }

        // Allocate gene storage
        TS = libc::malloc(NT * std::mem::size_of::<Gene>()) as *mut Gene;
        if TS.is_null() {
            eprintln!("Not enough memory available to store detected genes");
            std::process::exit(1);
        }
        sw.genespace = NT as i32;

        // Run detection
        if sw.libflag != 0 {
            let _ = writeln!(&mut *sw.f, "Library");
        }

        // Use parallel processing when threads != 1, regardless of batch mode
        if num_threads != 1 && !input_file_path.is_empty() {
            // Initialize global thread pool (only once, before parallel processing)
            let effective_threads = if num_threads == 0 {
                std::cmp::min(num_cpus::get_physical(), 6)  // Default max 6 threads
            } else {
                std::cmp::min(num_threads, 6)
            };
            let _ = rayon::ThreadPoolBuilder::new()
                .num_threads(effective_threads)
                .build_global();

            parallel_fastafile(&input_file_path, &mut d, &mut sw, num_threads);
        } else if sw.batch != 0 {
            bopt_fastafile(&mut d, &mut sw);
        } else {
            iopt_fastafile(&mut d, &mut sw);
        }

        // Cleanup
        libc::free(TS as *mut libc::c_void);
        if !d.f.is_null() {
            libc::fclose(d.f as *mut libc::FILE);
        }

        if sw.batch == 0 && sw.showconfig != 0 {
            let _ = write!(&mut *sw.f, "Configuration: ");
            for i in 0..z {
                let _ = write!(&mut *sw.f, "{} ", args[i as usize]);
            }
            let _ = writeln!(&mut *sw.f);
        }

        if sw.f != c_stdout() as *mut File {
            libc::fclose(sw.f as *mut libc::FILE);
        }
    }
}
