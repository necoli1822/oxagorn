// =============================================================================
// Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
// types.rs - Type definitions and constants
//
// Based on ARAGORN v1.2.41 by Dean Laslett
// Author: Sunju Kim <n.e.coli.1822@gmail.com>
// GPL License applies.
// =============================================================================
//
// Allow C-style naming for 1:1 compatibility with original C code
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]

// For C FFI: File is an opaque type representing C's FILE
// When using libc functions, cast *mut File to *mut libc::FILE
#[repr(C)]
pub struct File {
    _private: [u8; 0],  // Zero-sized, opaque
}

// Implement Write for File to support Rust's write!/writeln! macros
impl std::io::Write for File {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        unsafe {
            let ptr = self as *mut File as *mut libc::FILE;
            let written = libc::fwrite(
                buf.as_ptr() as *const libc::c_void,
                1,
                buf.len(),
                ptr
            );
            Ok(written)
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        unsafe {
            let ptr = self as *mut File as *mut libc::FILE;
            libc::fflush(ptr);
            Ok(())
        }
    }
}

// =============================================================================
// CONSTANTS - common.h lines 21-35
// =============================================================================
// C: #define NOCHAR '\0'
pub const NOCHAR: u8 = 0;
// C: #define DLIM '\n'
pub const DLIM: u8 = b'\n';
// C: #define STRLEN 4001
pub const STRLEN: usize = 4001;
// C: #define STRLENM1 4000
pub const STRLENM1: usize = 4000;
// C: #define SHORTSTRLEN 51
pub const SHORTSTRLEN: usize = 51;
// C: #define SHORTSTRLENM1 50
pub const SHORTSTRLENM1: usize = 50;
// C: #define KEYLEN 15
pub const KEYLEN: usize = 15;
// C: #define NHELPLINE 181 (original), 186 for Oxagorn (added -j, -@ options)
pub const NHELPLINE: usize = 186;
// C: #define INACTIVE 2.0e+35
pub const INACTIVE: f64 = 2.0e+35;
// C: #define IINACTIVE 2000000001L
pub const IINACTIVE: i64 = 2000000001;
// C: #define ITHRESHOLD 2000000000L
pub const ITHRESHOLD: i64 = 2000000000;

// =============================================================================
// FILE FORMAT TYPES - common.h lines 36-38
// =============================================================================
// C: #define FASTA 0
pub const FASTA: i32 = 0;
// C: #define GENBANK 1
pub const GENBANK: i32 = 1;

// =============================================================================
// GENE TYPES - common.h lines 39-45
// =============================================================================
// C: #define noGENE -1
pub const NO_GENE: i32 = -1;
// C: #define tRNA 0
pub const TRNA: i32 = 0;
// C: #define tmRNA 1
pub const TMRNA: i32 = 1;
// C: #define srpRNA 2
pub const SRPRNA: i32 = 2;
// C: #define rRNA 3
pub const RRNA: i32 = 3;
// C: #define CDS 4
pub const CDS: i32 = 4;
// C: #define NS 6
pub const NS: usize = 6;

// =============================================================================
// GENETIC CODE CONSTANTS - common.h lines 47-52
// =============================================================================
// C: #define MAXGCMOD 16
pub const MAXGCMOD: usize = 16;
// C: #define MAMMAL_MT 2
pub const MAMMAL_MT: i32 = 2;
// C: #define NGENECODE 34
pub const NGENECODE: usize = 34;
// C: #define METAZOAN_MT 0
pub const METAZOAN_MT: i32 = 0;
// C: #define STANDARD 1
pub const STANDARD: i32 = 1;
// C: #define VERTEBRATE_MT 2
pub const VERTEBRATE_MT: i32 = 2;

// =============================================================================
// AMINO ACID CONSTANTS - common.h lines 54-77
// =============================================================================
// C: #define NAMINOACID 27
pub const NAMINOACID: usize = 27;
// C: #define Phe 0 ... #define Pyl 22
pub const PHE: i32 = 0;
pub const VAL: i32 = 1;
pub const LEU: i32 = 2;
pub const ILE: i32 = 3;
pub const CYS: i32 = 4;
pub const GLY: i32 = 5;
pub const ARG: i32 = 6;
pub const SER: i32 = 7;
pub const ALA: i32 = 8;
pub const PRO: i32 = 9;
pub const THR: i32 = 10;
pub const TYR: i32 = 11;
pub const ASP: i32 = 12;
pub const HIS: i32 = 13;
pub const ASN: i32 = 14;
pub const MET: i32 = 15;
pub const TRP: i32 = 16;
pub const GLU: i32 = 17;
pub const GLN: i32 = 18;
pub const LYS: i32 = 19;
pub const STOP: i32 = 20;
// C: #define SeC 21
pub const SEC: i32 = 21;
// C: #define Pyl 22
pub const PYL: i32 = 22;

// =============================================================================
// NUCLEOTIDE CONSTANTS - common.h lines 79-86
// =============================================================================
// C: #define INSERT -2
pub const INSERT: i32 = -2;
// C: #define TERM -1
pub const TERM: i32 = -1;
// C: #define Adenine 0
pub const ADENINE: i32 = 0;
// C: #define Cytosine 1
pub const CYTOSINE: i32 = 1;
// C: #define Guanine 2
pub const GUANINE: i32 = 2;
// C: #define Thymine 3
pub const THYMINE: i32 = 3;
// C: #define AMBIG 4
pub const AMBIG: i32 = 4;
// C: #define NOBASE 5
pub const NOBASE: i32 = 5;

// =============================================================================
// THRESHOLD CONSTANTS - common.h lines 88-95
// =============================================================================
// C: #define tRNAthresh 132.0
pub const TRNA_THRESH: f64 = 132.0;
// C: #define mtRNAdtthresh 91.5
pub const MTRNA_DT_THRESH: f64 = 91.5;
// C: #define mtRNAtthresh 83.5
pub const MTRNA_T_THRESH: f64 = 83.5;
// C: #define mtRNAdthresh 85.0
pub const MTRNA_D_THRESH: f64 = 85.0;
// C: #define tmRNAthresh 325.0
pub const TMRNA_THRESH: f64 = 325.0;
// C: #define srpRNAthresh 175.0
pub const SRPRNA_THRESH: f64 = 175.0;
// C: #define CDSthresh 100.0
pub const CDS_THRESH: f64 = 100.0;
// C: #define PSEUDOGENElevel 0.95
pub const PSEUDOGENE_LEVEL: f64 = 0.95;

// =============================================================================
// DIRECTION CONSTANTS - common.h lines 97-106
// =============================================================================
// C: #define RIGHT 0 ... #define SLANT 5
pub const RIGHT: i32 = 0;
pub const UP: i32 = 1;
pub const LEFT: i32 = 2;
pub const DOWN: i32 = 3;
pub const UPRIGHT: i32 = 4;
pub const SLANTDR: i32 = 5;
pub const SLANTUR: i32 = 6;
pub const SLANTUL: i32 = 7;
pub const SLANTDL: i32 = 8;
pub const SLANT: i32 = 5;

// =============================================================================
// MATRIX DIMENSIONS - common.h lines 108-109
// =============================================================================
// C: #define MATX 42
pub const MATX: usize = 42;
// C: #define MATY 34
pub const MATY: usize = 34;

// =============================================================================
// STEM/LOOP CONSTANTS - common.h lines 111-149
// =============================================================================
// C: #define ASTEM2_EXT 9
pub const ASTEM2_EXT: i32 = 9;
// C: #define ASTEM2_EXTD 4
pub const ASTEM2_EXTD: i32 = 4;
// C: #define ASTEM2_EXTE 5
pub const ASTEM2_EXTE: i32 = 5;
// C: #define MINTSTEM_DIST (17 + ASTEM2_EXT)
pub const MINTSTEM_DIST: i32 = 17 + ASTEM2_EXT;
// C: #define MAXTSTEM_DIST (26 + ASTEM2_EXT)
pub const MAXTSTEM_DIST: i32 = 26 + ASTEM2_EXT;
// C: #define MAXDSTEM_DIST 9
pub const MAXDSTEM_DIST: i32 = 9;
// C: #define MINDSTEM_DIST 8
pub const MINDSTEM_DIST: i32 = 8;
// C: #define MININTRONLEN 0
pub const MININTRONLEN: i32 = 0;
// C: #define MAXINTRONLEN 3000
pub const MAXINTRONLEN: usize = 3000;
// C: #define MINCTRNALEN 62
pub const MINCTRNALEN: i32 = 62;
// C: #define MAXCTRNALEN 110
pub const MAXCTRNALEN: i32 = 110;
// C: #define MINTRNALEN (MINCTRNALEN + 1)
pub const MINTRNALEN: i32 = MINCTRNALEN + 1;
// C: #define MAXTRNALEN (MAXCTRNALEN + ASTEM2_EXT)
pub const MAXTRNALEN: usize = (MAXCTRNALEN + ASTEM2_EXT) as usize;
// C: #define MAXETRNALEN (MAXTRNALEN + MAXINTRONLEN)
pub const MAXETRNALEN: usize = MAXTRNALEN + MAXINTRONLEN;
// C: #define VARMAX 26
pub const VARMAX: i32 = 26;
// C: #define VARMIN 3
pub const VARMIN: i32 = 3;
// C: #define VARDIFF 23
pub const VARDIFF: i32 = 23;
// C: #define MINTPTSDIST 50
pub const MINTPTSDIST: i32 = 50;
// C: #define MAXTPTSDIST 321
pub const MAXTPTSDIST: i32 = 321;
// C: #define TPWINDOW (MAXTPTSDIST - MINTPTSDIST + 1)
pub const TPWINDOW: i32 = MAXTPTSDIST - MINTPTSDIST + 1;
// C: #define MINTPDIST 50
pub const MINTPDIST: i32 = 50;
// C: #define MAXTPDIST 250
pub const MAXTPDIST: i32 = 250;
// C: #define TPDISTWINDOW (MAXTPDIST - MINTPDIST + 1)
pub const TPDISTWINDOW: i32 = MAXTPDIST - MINTPDIST + 1;
// C: #define MINTAGDIST 12
pub const MINTAGDIST: i32 = 12;
// C: #define MAXTAGDIST 102
pub const MAXTAGDIST: i32 = 102;
// C: #define TAGWINDOW (MAXTAGDIST - MINTAGDIST)
pub const TAGWINDOW: i32 = MAXTAGDIST - MINTAGDIST;
// C: #define MINRNACDIST (MINTPDIST - 5)
pub const MINRNACDIST: i32 = MINTPDIST - 5;
// C: #define MAXRNACDIST (MAXTPDIST - 5)
pub const MAXRNACDIST: i32 = MAXTPDIST - 5;
// C: #define MAXPPINTRONDIST 250
pub const MAXPPINTRONDIST: i32 = 250;
// C: #define TMPTRAILER 145
pub const TMPTRAILER: i32 = 145;
// C: #define MINPPASDIST MINTSTEM_DIST
pub const MINPPASDIST: i32 = MINTSTEM_DIST;
// C: #define MAXPPASDIST (MAXTSTEM_DIST + MAXPPINTRONDIST)
pub const MAXPPASDIST: i32 = MAXTSTEM_DIST + MAXPPINTRONDIST;
// C: #define MINPPTSTPDIST (MINTSTEM_DIST + MINTPDIST)
pub const MINPPTSTPDIST: i32 = MINTSTEM_DIST + MINTPDIST;
// C: #define MAXPPTSTPDIST (MAXTSTEM_DIST+ASTEM2_EXT+MAXTPDIST+MAXPPINTRONDIST)
pub const MAXPPTSTPDIST: i32 = MAXTSTEM_DIST + ASTEM2_EXT + MAXTPDIST + MAXPPINTRONDIST;
// C: #define MAXTMRNALEN (4 + MAXPPASDIST + MAXTPDIST + MAXTAGDIST + TMPTRAILER)
pub const MAXTMRNALEN: i32 = 4 + MAXPPASDIST + MAXTPDIST + MAXTAGDIST + TMPTRAILER;
// C: #define TSWEEP 1000
pub const TSWEEP: i32 = 1000;
// C: #define WRAP (2*MAXETRNALEN)
pub const WRAP: usize = 2 * MAXETRNALEN;
// C: #define NPTAG 33
pub const NPTAG: i32 = 33;
// C: #define MAXAGENELEN (MAXETRNALEN + MAXTMRNALEN)
pub const MAXAGENELEN: usize = MAXETRNALEN + MAXTMRNALEN as usize;

// =============================================================================
// STEM TYPES - common.h lines 151-153
// =============================================================================
// C: #define BASE 0
pub const BASE: i32 = 0;
// C: #define FSTEM 1
pub const FSTEM: i32 = 1;
// C: #define BSTEM 2
pub const BSTEM: i32 = 2;

// =============================================================================
// LOOP IDENTIFICATION - common.h lines 155-159
// =============================================================================
// C: #define NOID 0
pub const NOID: i32 = 0;
// C: #define DLOOP 1
pub const DLOOP: i32 = 1;
// C: #define DSTEM 2
pub const DSTEM: i32 = 2;
// C: #define CLOOP 3
pub const CLOOP: i32 = 3;
// C: #define VAR 4
pub const VAR: i32 = 4;

// =============================================================================
// ARRAY SIZE CONSTANTS - common.h lines 161-196
// =============================================================================
// C: #define NA MAXINTRONLEN
pub const NA: usize = MAXINTRONLEN;
// C: #define ND 100
pub const ND: usize = 100;
// C: #define NT 200
// Increased from 200 to 2000 for large genomes (A. thaliana has 667 genes, C. elegans has ~650)
// The original C code dynamically reallocates, but find_slot_ts cannot realloc in thread-local storage
pub const NT: usize = 2000;
// C: #define NH 2000
pub const NH: usize = 2000;
// C: #define NTH 3000
pub const NTH: usize = 3000;
// C: #define NC 5000
pub const NC: usize = 5000;
// C: #define NGFT 5000
pub const NGFT: usize = 5000;
// C: #define NTAG 1273
pub const NTAG: usize = 1273;
// C: #define NTAGMAX 1300
pub const NTAGMAX: usize = 1300;
// C: #define LSEQ 20000
pub const LSEQ: usize = 20000;
// C: #define ATBOND 2.5
pub const ATBOND: f64 = 2.5;
// C: #define mtNA 1500
pub const MT_NA: usize = 1500;
// C: #define mtND 150
pub const MT_ND: usize = 150;
// C: #define mtNTH 3000
pub const MT_NTH: usize = 3000;
// C: #define mtNTM 3
pub const MT_NTM: usize = 3;
// C: #define mtNCDS 200
pub const MT_NCDS: usize = 200;
// C: #define mtNCDSCODON 6000
pub const MT_NCDSCODON: usize = 6000;
// C: #define mtGCBOND 0.0
pub const MT_GCBOND: f64 = 0.0;
// C: #define mtATBOND -0.5
pub const MT_ATBOND: f64 = -0.5;
// C: #define mtGTBOND -1.2
pub const MT_GTBOND: f64 = -1.2;
// C: #define mtTTBOND -2.9
pub const MT_TTBOND: f64 = -2.9;
// C: #define mtGGBOND -3.0
pub const MT_GGBOND: f64 = -3.0;
// C: #define mtGABOND -3.0
pub const MT_GABOND: f64 = -3.0;
// C: #define mtNOBOND -3.0
pub const MT_NOBOND: f64 = -3.0;
// C: #define mtBONDSTAB 1.5
pub const MT_BONDSTAB: f64 = 1.5;
// C: #define mtABONDSTAB 2.0
pub const MT_ABONDSTAB: f64 = 2.0;
// C: #define mtTSTTSTAB -2.5
pub const MT_TSTTSTAB: f64 = -2.5;
// C: #define mtTERMSTAB 0.01
pub const MT_TERMSTAB: f64 = 0.01;
// C: #define mtSENDSTAB 0.01
pub const MT_SENDSTAB: f64 = 0.01;
// C: #define mtNSTAB 0.1
pub const MT_NSTAB: f64 = 0.1;
// C: #define mt3MMSTAB 1.0
pub const MT_3MMSTAB: f64 = 1.0;
// C: #define mtGCPENALTY 0.8
pub const MT_GCPENALTY: f64 = 0.8;
// C: #define mtGCPENALTYD 2.0
pub const MT_GCPENALTYD: f64 = 2.0;
// C: #define mt_DRLmaxlength 16
pub const MT_DRLMAXLENGTH: i32 = 16;
// C: #define mt_TVRLmaxlength 18
pub const MT_TVRLMAXLENGTH: i32 = 18;
// C: #define mtNCLM 3
pub const MT_NCLM: usize = 3;

// =============================================================================
// rRNA CONSTANTS - common.h lines 198-201
// =============================================================================
// C: #define SRRNAMAXLEN 1500
pub const SRRNAMAXLEN: i32 = 1500;
// C: #define SRRNAMINLEN 600
pub const SRRNAMINLEN: i32 = 600;
// C: #define LRRNAMINLEN 1200
pub const LRRNAMINLEN: i32 = 1200;
// C: #define LRRNAMAXLEN 3000
pub const LRRNAMAXLEN: i32 = 3000;

// =============================================================================
// srpRNA CONSTANTS - common.h lines 203-216
// =============================================================================
// C: #define srpMAXLEN 650
pub const SRP_MAXLEN: i32 = 650;
// C: #define srpUMAXLEN 300
pub const SRP_UMAXLEN: i32 = 300;
// C: #define srpUMINLEN 100
pub const SRP_UMINLEN: i32 = 100;
// C: #define srpDMAXLEN 300
pub const SRP_DMAXLEN: i32 = 300;
// C: #define srpDMINLEN 100
pub const SRP_DMINLEN: i32 = 100;
// C: #define srpNH 200
pub const SRP_NH: usize = 200;
// C: #define srpNS 500
pub const SRP_NS: usize = 500;
// C: #define srpMAXHPL 14
pub const SRP_MAXHPL: i32 = 14;
// C: #define srpMAXSP 6
pub const SRP_MAXSP: i32 = 6;
// C: #define srpMAXSTEM 6500
pub const SRP_MAXSTEM: usize = 6500;
// C: #define srpDISPMAX (4*srpMAXLEN)
pub const SRP_DISPMAX: i32 = 4 * SRP_MAXLEN;
// C: #define srpMAXSPACER 12
pub const SRP_MAXSPACER: i32 = 12;
// C: #define srpMAXNISTEMS 10
pub const SRP_MAXNISTEMS: i32 = 10;
// C: #define srpNESTMAX 2
pub const SRP_NESTMAX: i32 = 2;

// =============================================================================
// CDS CONSTANTS - common.h lines 218-220
// =============================================================================
// C: #define cdsMAXLEN 3000
pub const CDS_MAXLEN: i32 = 3000;
// C: #define NCDS 200
pub const NCDS: usize = 200;
// C: #define NCDSCODON 1000
pub const NCDSCODON: usize = 1000;

// =============================================================================
// MACROS AS INLINE FUNCTIONS - common.h lines 32-34
// =============================================================================
// C: #define space(c) (c==' ')||(c=='\t')||(c=='\n')||(c=='\r')
#[inline]
pub fn space(c: u8) -> bool {
    c == b' ' || c == b'\t' || c == b'\n' || c == b'\r'
}

// C: #define sq(pos) ((pos + d->psmax - 1L) % d->psmax) + 1L
#[inline]
pub fn sq(pos: i64, psmax: i64) -> i64 {
    ((pos + psmax - 1) % psmax) + 1
}

// =============================================================================
// TYPE DEFINITIONS - common.h lines 222-457
// =============================================================================

// -----------------------------------------------------------------------------
// C: typedef struct { long start; ... char species[SHORTSTRLEN]; } annotated_gene;
// common.h lines 224-233
// -----------------------------------------------------------------------------
#[derive(Clone)]
pub struct AnnotatedGene {
    pub start: i64,          // C: long start
    pub stop: i64,           // C: long stop
    pub comp: i32,           // C: int comp
    pub antistart: i64,      // C: long antistart
    pub antistop: i64,       // C: long antistop
    pub genetype: i32,       // C: int genetype
    pub pseudogene: i32,     // C: int pseudogene
    pub permuted: i32,       // C: int permuted
    pub detected: i32,       // C: int detected
    pub species: [u8; SHORTSTRLEN],  // C: char species[SHORTSTRLEN]
}

impl Default for AnnotatedGene {
    fn default() -> Self {
        AnnotatedGene {
            start: 0,
            stop: 0,
            comp: 0,
            antistart: 0,
            antistop: 0,
            genetype: 0,
            pseudogene: 0,
            permuted: 0,
            detected: 0,
            species: [0; SHORTSTRLEN],
        }
    }
}

// -----------------------------------------------------------------------------
// C: typedef struct { char filename[80]; ... annotated_gene gene[NGFT]; } data_set;
// common.h lines 235-251
// -----------------------------------------------------------------------------
pub struct DataSet {
    pub filename: [u8; 80],       // C: char filename[80]
    pub f: *mut File,             // C: FILE *f (raw pointer for 1:1 mapping)
    pub seqname: [u8; STRLEN],    // C: char seqname[STRLEN]
    pub bugmode: i32,             // C: int bugmode
    pub datatype: i32,            // C: int datatype
    pub gc: f64,                  // C: double gc
    pub filepointer: i64,         // C: long filepointer
    pub ps: i64,                  // C: long ps
    pub psmax: i64,               // C: long psmax
    pub seqstart: i64,            // C: long seqstart
    pub seqstartoff: i64,         // C: long seqstartoff
    pub nextseq: i64,             // C: long nextseq
    pub nextseqoff: i64,          // C: long nextseqoff
    pub ns: i32,                  // C: int ns
    pub nf: i32,                  // C: int nf
    pub aseqlen: i64,             // C: long aseqlen
    pub nagene: [i32; NS],        // C: int nagene[NS]
    pub gene: Vec<AnnotatedGene>, // C: annotated_gene gene[NGFT]
}

impl Default for DataSet {
    fn default() -> Self {
        DataSet {
            filename: [0; 80],
            f: std::ptr::null_mut(),
            seqname: [0; STRLEN],
            bugmode: 0,
            datatype: 0,
            gc: 0.0,
            filepointer: 0,
            ps: 0,
            psmax: 0,
            seqstart: 0,
            seqstartoff: 0,
            nextseq: 0,
            nextseqoff: 0,
            ns: 0,
            nf: 0,
            aseqlen: 0,
            nagene: [0; NS],
            gene: Vec::with_capacity(NGFT),
        }
    }
}

// -----------------------------------------------------------------------------
// C: typedef struct { char name[100]; ... int annosc; } gene;
// common.h lines 253-283
// -----------------------------------------------------------------------------
#[derive(Clone)]
pub struct Gene {
    pub name: [u8; 100],          // C: char name[100]
    pub seq: [i32; MAXTRNALEN + 1],   // C: int seq[MAXTRNALEN+1]
    pub eseq: [i32; MAXETRNALEN + 1], // C: int eseq[MAXETRNALEN+1]
    pub ps: *mut i32,             // C: int *ps
    pub nbase: i32,               // C: int nbase
    pub comp: i32,                // C: int comp
    pub start: i64,               // C: long start
    pub stop: i64,                // C: long stop
    pub astem1: i32,              // C: int astem1
    pub astem2: i32,              // C: int astem2
    pub aatail: i32,              // C: int aatail
    pub spacer1: i32,             // C: int spacer1
    pub spacer2: i32,             // C: int spacer2
    pub dstem: i32,               // C: int dstem
    pub dloop: i32,               // C: int dloop
    pub cstem: i32,               // C: int cstem
    pub cloop: i32,               // C: int cloop
    pub intron: i32,              // C: int intron
    pub nintron: i32,             // C: int nintron
    pub anticodon: i32,           // C: int anticodon
    pub var: i32,                 // C: int var
    pub varbp: i32,               // C: int varbp
    pub tstem: i32,               // C: int tstem
    pub tloop: i32,               // C: int tloop
    pub genetype: i32,            // C: int genetype
    pub energy: f64,              // C: double energy
    pub asst: i32,                // C: int asst
    pub tps: i32,                 // C: int tps
    pub tpe: i32,                 // C: int tpe
    pub annotation: i32,          // C: int annotation
    pub annosc: i32,              // C: int annosc
}

impl Default for Gene {
    fn default() -> Self {
        // Initialize eseq with TERM (-1) at index 0 to prevent false "extended sequence" output
        let mut eseq = [0; MAXETRNALEN + 1];
        eseq[0] = TERM;
        Gene {
            name: [0; 100],
            seq: [0; MAXTRNALEN + 1],
            eseq,
            ps: std::ptr::null_mut(),
            nbase: 0,
            comp: 0,
            start: 0,
            stop: 0,
            astem1: 0,
            astem2: 0,
            aatail: 0,
            spacer1: 0,
            spacer2: 0,
            dstem: 0,
            dloop: 0,
            cstem: 0,
            cloop: 0,
            intron: 0,
            nintron: 0,
            anticodon: 0,
            var: 0,
            varbp: 0,
            tstem: 0,
            tloop: 0,
            genetype: 0,
            energy: 0.0,
            asst: 0,
            tps: 0,
            tpe: 0,
            annotation: 0,
            annosc: 0,
        }
    }
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos; int stem; int loop; double energy; } trna_loop;
// common.h lines 285-288
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct TrnaLoop {
    pub pos: *mut i32,    // C: int *pos
    pub stem: i32,        // C: int stem
    pub loop_: i32,       // C: int loop (loop is reserved in Rust)
    pub energy: f64,      // C: double energy
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos; ... double stem_energy; } mt_trna_loop;
// common.h lines 290-295
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct MtTrnaLoop {
    pub pos: *mut i32,        // C: int *pos
    pub stem: i32,            // C: int stem
    pub loop_: i32,           // C: int loop
    pub bondtype: u32,        // C: unsigned int bondtype
    pub energy: f64,          // C: double energy
    pub stem_energy: f64,     // C: double stem_energy
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos; ... double stem_energy; } mt_trna_cloop;
// common.h lines 297-306
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct MtTrnaCloop {
    pub pos: *mut i32,        // C: int *pos
    pub looppos: *mut i32,    // C: int *looppos
    pub end: *mut i32,        // C: int *end
    pub stem: i32,            // C: int stem
    pub loop_: i32,           // C: int loop
    pub arm: i32,             // C: int arm
    pub anticodon: i32,       // C: int anticodon
    pub bondtype: u32,        // C: unsigned int bondtype
    pub energy: f64,          // C: double energy
    pub stem_energy: f64,     // C: double stem_energy
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos; ... double stem_energy; } mt_trna_tloop;
// common.h lines 308-314
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct MtTrnaTloop {
    pub pos: *mut i32,        // C: int *pos
    pub stem: i32,            // C: int stem
    pub loop_: i32,           // C: int loop
    pub end: *mut i32,        // C: int *end
    pub bondtype: u32,        // C: unsigned int bondtype
    pub energy: f64,          // C: double energy
    pub stem_energy: f64,     // C: double stem_energy
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos; int *end; int stem; int loop; double energy; } trna_dloop;
// common.h lines 316-320
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct TrnaDloop {
    pub pos: *mut i32,    // C: int *pos
    pub end: *mut i32,    // C: int *end
    pub stem: i32,        // C: int stem
    pub loop_: i32,       // C: int loop
    pub energy: f64,      // C: double energy
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos1; int *pos2; int stem; double energy; } trna_astem;
// common.h lines 322-325
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct TrnaAstem {
    pub pos1: *mut i32,   // C: int *pos1
    pub pos2: *mut i32,   // C: int *pos2
    pub stem: i32,        // C: int stem
    pub energy: f64,      // C: double energy
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos1; ... double energy; } mt_trna_astem;
// common.h lines 327-331
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct MtTrnaAstem {
    pub pos1: *mut i32,   // C: int *pos1
    pub pos2: *mut i32,   // C: int *pos2
    pub stem: i32,        // C: int stem
    pub bondtype: u32,    // C: unsigned int bondtype
    pub energy: f64,      // C: double energy
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos; int comp; int frame; int codon; int win; } mt_cds_codon;
// common.h lines 333-337
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct MtCdsCodon {
    pub pos: *mut i32,    // C: int *pos
    pub comp: i32,        // C: int comp
    pub frame: i32,       // C: int frame
    pub codon: i32,       // C: int codon
    pub win: i32,         // C: int win
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos1; int *pos2; int comp; } mt_cds;
// common.h lines 339-341
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct MtCds {
    pub pos1: *mut i32,   // C: int *pos1
    pub pos2: *mut i32,   // C: int *pos2
    pub comp: i32,        // C: int comp
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos1; int *pos2; int comp; } mt_rrna;
// common.h lines 343-345
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct MtRrna {
    pub pos1: *mut i32,   // C: int *pos1
    pub pos2: *mut i32,   // C: int *pos2
    pub comp: i32,        // C: int comp
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos1; int *pos2; int stem; int loop; } rrna_hairpin;
// common.h lines 347-350
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct RrnaHairpin {
    pub pos1: *mut i32,   // C: int *pos1
    pub pos2: *mut i32,   // C: int *pos2
    pub stem: i32,        // C: int stem
    pub loop_: i32,       // C: int loop
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos1; int *pos2; int stem; } rrna_stem;
// common.h lines 352-354
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct RrnaStem {
    pub pos1: *mut i32,   // C: int *pos1
    pub pos2: *mut i32,   // C: int *pos2
    pub stem: i32,        // C: int stem
}

// -----------------------------------------------------------------------------
// C: typedef struct { int *pos; int comp; int frame; int codon; int win; } cds_codon;
// common.h lines 356-360
// -----------------------------------------------------------------------------
#[derive(Clone, Copy, Default)]
pub struct CdsCodon {
    pub pos: *mut i32,    // C: int *pos
    pub comp: i32,        // C: int comp
    pub frame: i32,       // C: int frame
    pub codon: i32,       // C: int codon
    pub win: i32,         // C: int win
}

// -----------------------------------------------------------------------------
// C: typedef struct { char name[50]; char tag[50]; } tmrna_tag_entry;
// common.h lines 362-363
// -----------------------------------------------------------------------------
#[derive(Clone, Copy)]
#[repr(C)]
pub struct TmrnaTagEntry {
    pub name: [u8; 50],   // C: char name[50]
    pub tag: [u8; 50],    // C: char tag[50]
}

impl TmrnaTagEntry {
    /// Create entry from string literals (const fn for static initialization)
    pub const fn new(name_str: &str, tag_str: &str) -> Self {
        let mut name = [0u8; 50];
        let mut tag = [0u8; 50];

        let name_bytes = name_str.as_bytes();
        let tag_bytes = tag_str.as_bytes();

        let mut i = 0;
        while i < name_bytes.len() && i < 49 {
            name[i] = name_bytes[i];
            i += 1;
        }

        i = 0;
        while i < tag_bytes.len() && i < 49 {
            tag[i] = tag_bytes[i];
            i += 1;
        }

        TmrnaTagEntry { name, tag }
    }

    pub const fn empty() -> Self {
        TmrnaTagEntry { name: [0u8; 50], tag: [0u8; 50] }
    }
}

impl Default for TmrnaTagEntry {
    fn default() -> Self {
        TmrnaTagEntry::empty()
    }
}

// -----------------------------------------------------------------------------
// C: typedef struct { char genetypename[NS][10]; ... int tmrna_struct[200]; } csw;
// common.h lines 365-457
// -----------------------------------------------------------------------------
pub struct Csw {
    pub genetypename: [[u8; 10]; NS],  // C: char genetypename[NS][10]
    pub f: *mut File,                   // C: FILE *f (raw pointer for 1:1 mapping)
    pub batch: i32,                     // C: int batch
    pub batchfullspecies: i32,          // C: int batchfullspecies
    pub repeatsn: i32,                  // C: int repeatsn
    pub trna: i32,                      // C: int trna
    pub tmrna: i32,                     // C: int tmrna
    pub srprna: i32,                    // C: int srprna
    pub cds: i32,                       // C: int cds
    pub mtrna: i32,                     // C: int mtrna
    pub tvloop: i32,                    // C: int tvloop
    pub cloop7: i32,                    // C: int cloop7
    pub peptide: i32,                   // C: int peptide
    pub geneticcode: i32,               // C: int geneticcode
    pub ngcmod: i32,                    // C: int ngcmod
    pub gcmod: [i32; MAXGCMOD],         // C: int gcmod[MAXGCMOD]
    pub gcfix: i32,                     // C: int gcfix
    pub discrim: i32,                   // C: int discrim
    pub extastem: i32,                  // C: int extastem
    pub tarm: i32,                      // C: int tarm
    pub tagthresh: i32,                 // C: int tagthresh
    pub tarmlength: i32,                // C: int tarmlength
    pub showconfig: i32,                // C: int showconfig
    pub libflag: i32,                   // C: int libflag
    pub verbose: i32,                   // C: int verbose
    pub linear: i32,                    // C: int linear
    pub both: i32,                      // C: int both
    pub reportpseudogenes: i32,         // C: int reportpseudogenes
    pub energydisp: i32,                // C: int energydisp
    pub secstructdisp: i32,             // C: int secstructdisp
    pub seqdisp: i32,                   // C: int seqdisp
    pub aataildisp: i32,                // C: int aataildisp
    pub aataildiv: i32,                 // C: int aataildiv
    pub sp1max: i32,                    // C: int sp1max
    pub sp2min: i32,                    // C: int sp2min
    pub sp2max: i32,                    // C: int sp2max
    pub mtxdetect: i32,                 // C: int mtxdetect
    pub mtcdsscan: i32,                 // C: int mtcdsscan
    pub mtcompov: i32,                  // C: int mtcompov
    pub matchacceptor: i32,             // C: int matchacceptor
    pub maxintronlen: i32,              // C: int maxintronlen
    pub minintronlen: i32,              // C: int minintronlen
    pub minintronlenreport: i32,        // C: int minintronlenreport
    pub ioverlay: i32,                  // C: int ioverlay
    pub ifixedpos: i32,                 // C: int ifixedpos
    pub ireportminintronlen: i32,       // C: int ireportminintronlen
    pub tmstrict: i32,                  // C: int tmstrict
    pub iamismatch: i32,                // C: int iamismatch
    pub loffset: i32,                   // C: int loffset
    pub roffset: i32,                   // C: int roffset
    pub start: i64,                     // C: long start
    pub comp: i32,                      // C: int comp
    pub genespace: i32,                 // C: int genespace
    pub srpspace: i32,                  // C: int srpspace
    pub ngene: [i32; NS],               // C: int ngene[NS]
    pub nps: i32,                       // C: int nps
    pub annotated: i32,                 // C: int annotated
    pub dispmatch: i32,                 // C: int dispmatch
    pub updatetmrnatags: i32,           // C: int updatetmrnatags
    pub tagend: i32,                    // C: int tagend
    pub trnalenmisthresh: i32,          // C: int trnalenmisthresh
    pub tmrnalenmisthresh: i32,         // C: int tmrnalenmisthresh
    pub nagene: [i32; NS],              // C: int nagene[NS]
    pub nafn: [i32; NS],                // C: int nafn[NS]
    pub nafp: [i32; NS],                // C: int nafp[NS]
    pub natfpd: i32,                    // C: int natfpd
    pub natfptv: i32,                   // C: int natfptv
    pub lacds: i32,                     // C: int lacds
    pub ldcds: i32,                     // C: int ldcds
    pub nabase: i64,                    // C: long nabase
    pub reportpsthresh: f64,            // C: double reportpsthresh
    pub threshlevel: f64,               // C: double threshlevel
    pub trnathresh: f64,                // C: double trnathresh
    pub ttscanthresh: f64,              // C: double ttscanthresh
    pub ttarmthresh: f64,               // C: double ttarmthresh
    pub tdarmthresh: f64,               // C: double tdarmthresh
    pub tastemthresh: f64,              // C: double tastemthresh
    pub tascanthresh: f64,              // C: double tascanthresh
    pub mttthresh: f64,                 // C: double mttthresh
    pub mtdthresh: f64,                 // C: double mtdthresh
    pub mtdtthresh: f64,                // C: double mtdtthresh
    pub mttarmthresh: f64,              // C: double mttarmthresh
    pub mtdarmthresh: f64,              // C: double mtdarmthresh
    pub tmrnathresh: f64,               // C: double tmrnathresh
    pub tmathresh: f64,                 // C: double tmathresh
    pub tmcthresh: f64,                 // C: double tmcthresh
    pub tmcathresh: f64,                // C: double tmcathresh
    pub tmrthresh: f64,                 // C: double tmrthresh
    pub srpthresh: f64,                 // C: double srpthresh
    pub cdsthresh: f64,                 // C: double cdsthresh
    pub eref: [f64; NS],                // C: double eref[NS]
    pub tmrna_struct: [i32; 200],       // C: int tmrna_struct[200]
}

impl Default for Csw {
    fn default() -> Self {
        Csw {
            genetypename: [[0; 10]; NS],
            f: std::ptr::null_mut(),
            batch: 0,
            batchfullspecies: 0,
            repeatsn: 0,
            trna: 0,
            tmrna: 0,
            srprna: 0,
            cds: 0,
            mtrna: 0,
            tvloop: 1,  // C default: 1
            cloop7: 0,
            peptide: 0,
            geneticcode: 1,  // STANDARD genetic code
            ngcmod: 0,
            gcmod: [0; MAXGCMOD],
            gcfix: 0,
            discrim: 0,
            extastem: 1,  // C default: 1
            tarm: 0,
            tagthresh: 0,
            tarmlength: 0,
            showconfig: 0,
            libflag: 0,
            verbose: 0,
            linear: 0,
            both: 2,  // C default: 2 (search both strands)
            reportpseudogenes: 0,
            energydisp: 0,
            secstructdisp: 0,
            seqdisp: 0,
            aataildisp: 0,
            aataildiv: 0,
            sp1max: 3,
            sp2min: 0,
            sp2max: 2,
            mtxdetect: 0,
            mtcdsscan: 0,
            mtcompov: 0,
            matchacceptor: 0,
            maxintronlen: 0,
            minintronlen: 0,
            minintronlenreport: 0,
            ioverlay: 0,
            ifixedpos: 0,
            ireportminintronlen: 0,
            tmstrict: 1,  // C default: 1
            iamismatch: 0,
            loffset: 0,
            roffset: 0,
            start: 0,
            comp: 0,
            genespace: 0,
            srpspace: 0,
            ngene: [0; NS],
            nps: 0,
            annotated: 0,
            dispmatch: 0,
            updatetmrnatags: 0,
            tagend: 0,
            trnalenmisthresh: 0,
            tmrnalenmisthresh: 0,
            nagene: [0; NS],
            nafn: [0; NS],
            nafp: [0; NS],
            natfpd: 0,
            natfptv: 0,
            lacds: 0,
            ldcds: 0,
            nabase: 0,
            reportpsthresh: 100.0,
            threshlevel: 1.0,
            trnathresh: 132.0,      // tRNAthresh
            ttscanthresh: 4.0,
            ttarmthresh: 29.0,
            tdarmthresh: 26.0,
            tastemthresh: 7.5,
            tascanthresh: 8.0,
            mttthresh: 83.5,        // mtRNAtthresh
            mtdthresh: 85.0,        // mtRNAdthresh
            mtdtthresh: 91.5,       // mtRNAdtthresh
            mttarmthresh: -7.9,
            mtdarmthresh: -6.0,
            tmrnathresh: 325.0,     // tmRNAthresh
            tmathresh: 14.0,
            tmcthresh: 10.0,
            tmcathresh: 25.0,
            tmrthresh: 9.0,
            srpthresh: 175.0,       // srpRNAthresh
            cdsthresh: 100.0,       // CDSthresh
            eref: [132.0, 325.0, 175.0, 0.0, 100.0, 0.0],  // tRNAthresh, tmRNAthresh, srpRNAthresh, 0.0, CDSthresh, 0.0
            tmrna_struct: [0; 200],
        }
    }
}

// =============================================================================
// C-STYLE ALIASES FOR 1:1 COMPATIBILITY
// These allow utils.rs and other modules to use C-style naming
// =============================================================================

// Gene type aliases (lowercase as in C)
pub const tRNA: i32 = TRNA;
pub const tmRNA: i32 = TMRNA;
pub const srpRNA: i32 = SRPRNA;
pub const rRNA: i32 = RRNA;

// Nucleotide aliases (PascalCase as in C)
pub const Adenine: i32 = ADENINE;
pub const Cytosine: i32 = CYTOSINE;
pub const Guanine: i32 = GUANINE;
pub const Thymine: i32 = THYMINE;

// Threshold aliases (camelCase as in C)
pub const tRNAthresh: f64 = TRNA_THRESH;
pub const mtRNAtthresh: f64 = MTRNA_T_THRESH;
pub const mtRNAdthresh: f64 = MTRNA_D_THRESH;
pub const mtRNAdtthresh: f64 = MTRNA_DT_THRESH;
pub const tmRNAthresh: f64 = TMRNA_THRESH;
pub const srpRNAthresh: f64 = SRPRNA_THRESH;
pub const CDSthresh: f64 = CDS_THRESH;
pub const PSEUDOGENElevel: f64 = PSEUDOGENE_LEVEL;

// MT constant aliases (lowercase as in C)
pub const mtNA: usize = MT_NA;
pub const mtND: usize = MT_ND;
pub const mtNTH: usize = MT_NTH;
pub const mtNTM: usize = MT_NTM;
pub const mtNCDS: usize = MT_NCDS;
pub const mtNCDSCODON: usize = MT_NCDSCODON;
pub const mtGCBOND: f64 = MT_GCBOND;
pub const mtATBOND: f64 = MT_ATBOND;
pub const mtGTBOND: f64 = MT_GTBOND;
pub const mtTTBOND: f64 = MT_TTBOND;
pub const mtGGBOND: f64 = MT_GGBOND;
pub const mtGABOND: f64 = MT_GABOND;
pub const mtNOBOND: f64 = MT_NOBOND;
pub const mtBONDSTAB: f64 = MT_BONDSTAB;
pub const mtABONDSTAB: f64 = MT_ABONDSTAB;
pub const mtTSTTSTAB: f64 = MT_TSTTSTAB;
pub const mtTERMSTAB: f64 = MT_TERMSTAB;
pub const mtSENDSTAB: f64 = MT_SENDSTAB;
pub const mtNSTAB: f64 = MT_NSTAB;
pub const mt3MMSTAB: f64 = MT_3MMSTAB;
pub const mtGCPENALTY: f64 = MT_GCPENALTY;
pub const mtGCPENALTYD: f64 = MT_GCPENALTYD;
pub const mt_DRLmaxlength: i32 = MT_DRLMAXLENGTH;
pub const mt_TVRLmaxlength: i32 = MT_TVRLMAXLENGTH;
pub const mtNCLM: usize = MT_NCLM;

// EOF constant (not in original common.h, but needed for file operations)
pub const EOF: i32 = -1;

// Type aliases for C-style struct names
pub type data_set = DataSet;
pub type csw = Csw;
pub type annotated_gene = AnnotatedGene;
pub type gene = Gene;
pub type trna_loop = TrnaLoop;
pub type mt_trna_loop = MtTrnaLoop;
pub type mt_trna_cloop = MtTrnaCloop;
pub type mt_trna_tloop = MtTrnaTloop;
pub type trna_dloop = TrnaDloop;
pub type trna_astem = TrnaAstem;
pub type mt_trna_astem = MtTrnaAstem;
pub type mt_cds_codon = MtCdsCodon;
pub type mt_cds = MtCds;
pub type mt_rrna = MtRrna;
pub type rrna_hairpin = RrnaHairpin;
pub type rrna_stem = RrnaStem;
pub type cds_codon = CdsCodon;
pub type tmrna_tag_entry = TmrnaTagEntry;
