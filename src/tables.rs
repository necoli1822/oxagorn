// =============================================================================
// tables.rs - Rust implementation of data_tables.c
// =============================================================================
//
// SOURCE FILE: src_flatten/data_tables.c
// This file is a 1:1 port of the C source file data_tables.c
//
// MAPPING GUIDE:
//   C int array[N][M]         →  pub static ARRAY: [[i32; M]; N]
//   C double array[N][M]      →  pub static ARRAY: [[f64; M]; N]
//   C unsigned int array[]    →  pub static ARRAY: [[u32; N]; M]
//   C char string[N]          →  pub static STRING: &[u8; N] or [[u8; N]; M]
//   C pointer *var            →  pub static mut VAR: *mut Type
//
// NOTES:
//   - All array dimensions and values are preserved exactly from C source
//   - Variable names are converted to UPPER_SNAKE_CASE for Rust constants
//   - C uses row-major order, Rust preserves this layout
//
// =============================================================================

/*
 * Oxagorn v0.1.0 - Rust reimplementation of ARAGORN
 * tables.rs - Global data arrays and tables
 *
 * Based on ARAGORN v1.2.41 data_tables.c by Dean Laslett
 * Author: Sunju Kim <n.e.coli.1822@gmail.com>
 * GPL License applies.
 */

use crate::types::*;

// =============================================================================
// BASEPAIR MATCHING MATRICES
// data_tables.c lines 8-186
// =============================================================================

// -----------------------------------------------------------------------------
// C: int lbp[3][6][6] = { ... };
// data_tables.c lines 10-28
// -----------------------------------------------------------------------------
pub static LBP: [[[i32; 6]; 6]; 3] = [
    [
        [0, 0, 1, 1, 1, 0],
        [0, 0, 1, 0, 1, 0],
        [1, 1, 0, 1, 1, 0],
        [1, 0, 1, 0, 1, 0],
        [1, 1, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0],
    ],
    [
        [0, 0, 0, 1, 1, 0],
        [0, 0, 1, 0, 1, 0],
        [0, 1, 0, 1, 1, 0],
        [1, 0, 1, 0, 1, 0],
        [1, 1, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0],
    ],
    [
        [0, 0, 0, 1, 1, 0],
        [0, 0, 1, 0, 1, 0],
        [0, 1, 0, 0, 1, 0],
        [1, 0, 0, 0, 1, 0],
        [1, 1, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0],
    ],
];

// -----------------------------------------------------------------------------
// C: int bp[6][6] = { ... };
// data_tables.c lines 30-35
// -----------------------------------------------------------------------------
pub static BP: [[i32; 6]; 6] = [
    [0, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 1, 0, 1, 1, 0],
    [1, 0, 1, 0, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int wbp[6][6] = { ... };
// data_tables.c lines 37-43
// -----------------------------------------------------------------------------
pub static WBP: [[i32; 6]; 6] = [
    [0, 0, 0, 2, 2, 0],
    [0, 0, 2, 0, 2, 0],
    [0, 2, 0, 1, 2, 0],
    [2, 0, 1, 0, 2, 0],
    [2, 2, 2, 2, 2, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int wcbp[6][6] = { ... };
// data_tables.c lines 45-50
// -----------------------------------------------------------------------------
pub static WCBP: [[i32; 6]; 6] = [
    [0, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 1, 0, 0, 1, 0],
    [1, 0, 0, 0, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int gc[6][6] = { ... };
// data_tables.c lines 52-57
// -----------------------------------------------------------------------------
pub static GC: [[i32; 6]; 6] = [
    [0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 1, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 1, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int gt[6][6] = { ... };
// data_tables.c lines 59-64
// -----------------------------------------------------------------------------
pub static GT: [[i32; 6]; 6] = [
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int at[6][6] = { ... };
// data_tables.c lines 66-71
// -----------------------------------------------------------------------------
pub static AT: [[i32; 6]; 6] = [
    [0, 0, 0, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int tt[6][6] = { ... };
// data_tables.c lines 73-78
// -----------------------------------------------------------------------------
pub static TT: [[i32; 6]; 6] = [
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0],
    [0, 0, 0, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int stemterm[6][6] = { ... };
// data_tables.c lines 80-85
// -----------------------------------------------------------------------------
pub static STEMTERM: [[i32; 6]; 6] = [
    [0, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 1, 0],
    [0, 0, 0, 1, 1, 0],
    [1, 0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int aastemterm[6][6] = { ... };
// data_tables.c lines 87-93
// -----------------------------------------------------------------------------
pub static AASTEMTERM: [[i32; 6]; 6] = [
    [1, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 1, 0],
    [0, 0, 0, 1, 1, 0],
    [1, 0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int ggstemterm[6][6] = { ... };
// data_tables.c lines 95-101
// -----------------------------------------------------------------------------
pub static GGSTEMTERM: [[i32; 6]; 6] = [
    [0, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 1, 0],
    [0, 0, 0, 1, 1, 0],
    [1, 0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int assymst[6][6] = { ... };
// data_tables.c lines 103-108
// -----------------------------------------------------------------------------
pub static ASSYMST: [[i32; 6]; 6] = [
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 1, 0],
    [0, 0, 0, 1, 1, 0],
    [1, 0, 0, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int assymat[6][6] = { ... };
// data_tables.c lines 110-115
// -----------------------------------------------------------------------------
pub static ASSYMAT: [[i32; 6]; 6] = [
    [0, 0, 0, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int stackbp[6][6] = { ... };
// data_tables.c lines 117-122
// -----------------------------------------------------------------------------
pub static STACKBP: [[i32; 6]; 6] = [
    [0, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 1, 0, 1, 1, 0],
    [1, 0, 1, 1, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int ggstackbp[6][6] = { ... };
// data_tables.c lines 124-130
// -----------------------------------------------------------------------------
pub static GGSTACKBP: [[i32; 6]; 6] = [
    [0, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 1, 1, 1, 1, 0],
    [1, 0, 1, 1, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int ggbp[6][6] = { ... };
// data_tables.c lines 132-138
// -----------------------------------------------------------------------------
pub static GGBP: [[i32; 6]; 6] = [
    [0, 0, 0, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 1, 1, 1, 1, 0],
    [1, 0, 1, 0, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int gabp[6][6] = { ... };
// data_tables.c lines 140-146
// -----------------------------------------------------------------------------
pub static GABP: [[i32; 6]; 6] = [
    [0, 0, 1, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [1, 1, 0, 1, 1, 0],
    [1, 0, 1, 0, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int assymagbp[6][6] = { ... };
// data_tables.c lines 148-154
// -----------------------------------------------------------------------------
pub static ASSYMAGBP: [[i32; 6]; 6] = [
    [0, 0, 1, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [0, 1, 0, 1, 1, 0],
    [1, 0, 1, 0, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int stembp[6][6] = { ... };
// data_tables.c lines 156-162
// -----------------------------------------------------------------------------
pub static STEMBP: [[i32; 6]; 6] = [
    [0, 0, 1, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [1, 1, 0, 1, 1, 0],
    [1, 0, 1, 1, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int ggstembp[6][6] = { ... };
// data_tables.c lines 164-170
// -----------------------------------------------------------------------------
pub static GGSTEMBP: [[i32; 6]; 6] = [
    [0, 0, 1, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [1, 0, 1, 1, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int gastembp[6][6] = { ... };
// data_tables.c lines 172-178
// -----------------------------------------------------------------------------
pub static GASTEMBP: [[i32; 6]; 6] = [
    [1, 0, 1, 1, 1, 0],
    [0, 0, 1, 0, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [1, 0, 1, 1, 1, 0],
    [1, 1, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: int vbp[6][6] = { ... };
// data_tables.c lines 180-186
// -----------------------------------------------------------------------------
pub static VBP: [[i32; 6]; 6] = [
    [0, 0, 1, 4, 4, 0],
    [0, 0, 4, 0, 4, 0],
    [1, 4, 0, 2, 4, 0],
    [4, 0, 2, 0, 4, 0],
    [4, 4, 4, 4, 4, 0],
    [0, 0, 0, 0, 0, 0],
];

// =============================================================================
// TANDEM ARRAYS
// data_tables.c lines 188-193
// =============================================================================

// -----------------------------------------------------------------------------
// C: int tandemid[mtNTM][4] = { ... };
// data_tables.c lines 188-191
// -----------------------------------------------------------------------------
pub static TANDEMID: [[i32; 4]; MT_NTM] = [
    [3, 2, 2, 3],
    [2, 3, 3, 2],
    [3, 3, 3, 3],
];

// -----------------------------------------------------------------------------
// C: double tandem_em[mtNTM] = { -0.5,-0.5,2.0 };
// data_tables.c line 193
// -----------------------------------------------------------------------------
pub static TANDEM_EM: [f64; MT_NTM] = [-0.5, -0.5, 2.0];

// =============================================================================
// ENERGY MATRICES
// data_tables.c lines 195-248
// =============================================================================

// -----------------------------------------------------------------------------
// C: double send_em[6][6] = { ... };
// data_tables.c lines 195-201
// Note: Uses mtSENDSTAB constant from common.h (MT_SENDSTAB in Rust)
// -----------------------------------------------------------------------------
pub static SEND_EM: [[f64; 6]; 6] = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.5 * MT_SENDSTAB, 0.0, 0.5 * MT_SENDSTAB, 0.0],
    [0.0, 0.5 * MT_SENDSTAB, 0.0, MT_SENDSTAB, MT_SENDSTAB, 0.0],
    [0.0, 0.0, MT_SENDSTAB, 0.0, MT_SENDSTAB, 0.0],
    [0.0, 0.5 * MT_SENDSTAB, MT_SENDSTAB, MT_SENDSTAB, MT_SENDSTAB, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
];

// -----------------------------------------------------------------------------
// C: double ssend_em[6][6] = { ... };
// data_tables.c lines 203-209
// Note: Uses mtSENDSTAB constant from common.h (MT_SENDSTAB in Rust)
// -----------------------------------------------------------------------------
pub static SSEND_EM: [[f64; 6]; 6] = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, MT_SENDSTAB, 0.0, MT_SENDSTAB, 0.0],
    [0.0, MT_SENDSTAB, 0.0, MT_SENDSTAB, MT_SENDSTAB, 0.0],
    [0.0, 0.0, MT_SENDSTAB, 0.0, MT_SENDSTAB, 0.0],
    [0.0, MT_SENDSTAB, MT_SENDSTAB, MT_SENDSTAB, MT_SENDSTAB, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
];

// -----------------------------------------------------------------------------
// C: int neighbour_map[6][6] = { ... };
// data_tables.c lines 211-217
// -----------------------------------------------------------------------------
pub static NEIGHBOUR_MAP: [[i32; 6]; 6] = [
    [0, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 1, 0],
    [0, 0, 0, 1, 1, 0],
    [1, 0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 0],
];

// -----------------------------------------------------------------------------
// C: double neighbour_em[2][6][6] = { ... };
// data_tables.c lines 219-232
// Note: Uses mtNSTAB constant from common.h (MT_NSTAB in Rust)
// -----------------------------------------------------------------------------
pub static NEIGHBOUR_EM: [[[f64; 6]; 6]; 2] = [
    [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ],
    [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, MT_NSTAB, 0.0, MT_NSTAB, 0.0],
        [0.0, MT_NSTAB, 0.0, 0.0, MT_NSTAB, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, MT_NSTAB, MT_NSTAB, 0.0, MT_NSTAB, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ],
];

// -----------------------------------------------------------------------------
// C: unsigned int btmap[6][6] = { ... };
// data_tables.c lines 234-240
// -----------------------------------------------------------------------------
pub static BTMAP: [[u32; 6]; 6] = [
    [0x10000, 0x10000, 0x1000, 0x10, 0x00000, 0x10000],
    [0x10000, 0x10000, 0x1, 0x10000, 0x00000, 0x10000],
    [0x1000, 0x1, 0x10000, 0x100, 0x00000, 0x10000],
    [0x10, 0x10000, 0x100, 0x1000, 0x00000, 0x10000],
    [0x00000, 0x00000, 0x00000, 0x00000, 0x00000, 0x10000],
    [0x10000, 0x10000, 0x10000, 0x10000, 0x10000, 0x10000],
];

// -----------------------------------------------------------------------------
// C: double bem[6][6] = { ... };
// data_tables.c lines 242-248
// Note: Uses ATBOND constant from common.h
// -----------------------------------------------------------------------------
pub static BEM: [[f64; 6]; 6] = [
    [-2.144, -0.428, -2.144, ATBOND, 0.000, 0.000],
    [-0.428, -2.144, 3.000, -2.144, 0.000, 0.000],
    [-2.144, 3.000, -2.144, 1.286, 0.000, 0.000],
    [ATBOND, -2.144, 1.286, -0.428, 0.000, 0.000],
    [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
    [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
];

// =============================================================================
// MT_DISCRIM - 3D array [3][64][6]
// data_tables.c lines 250-310
// =============================================================================

// -----------------------------------------------------------------------------
// C: int mt_discrim[3][64][6] = { ... };
// data_tables.c lines 250-310
// Index 0: metazoan mt
// Index 1: standard
// Index 2: mammal mt
// -----------------------------------------------------------------------------
pub static MT_DISCRIM: [[[i32; 6]; 64]; 3] = [
    // metazoan mt (index 0)
    // data_tables.c lines 252-270
    [
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [0,0,0,0,0,0], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [0,0,0,0,0,0], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
    ],
    // standard (index 1)
    // data_tables.c lines 272-290
    [
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1], [1,1,1,1,1,1],
    ],
    // mammal mt (index 2)
    // data_tables.c lines 292-310
    [
        [1,0,0,0,1,1], [1,0,0,0,1,1], [1,0,0,0,1,1], [1,0,0,0,1,1],
        [0,0,0,1,1,1], [1,0,0,0,1,1], [1,0,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,0,1,0,1,1], [1,0,0,0,1,1], [1,1,1,1,1,1],
        [1,0,0,0,1,1], [1,1,1,1,1,1], [0,1,0,0,1,1], [0,0,1,0,1,1],
        [1,0,0,0,1,1], [1,0,0,0,1,1], [1,0,0,0,1,1], [1,0,1,0,1,1],
        [0,0,1,0,1,1], [1,0,0,0,1,1], [1,0,1,1,1,1], [1,0,1,1,1,1],
        [1,1,1,1,1,1], [1,0,1,0,1,1], [1,0,0,0,1,1], [1,1,1,1,1,1],
        [0,0,0,0,0,0], [1,0,1,1,1,1], [0,0,1,0,1,1], [1,1,1,1,1,1],
        [1,0,0,0,1,1], [1,0,0,0,1,1], [1,0,0,0,1,1], [1,0,1,0,1,1],
        [0,1,0,1,1,1], [1,0,0,0,1,1], [1,0,1,1,1,1], [1,1,1,1,1,1],
        [1,1,1,1,1,1], [1,0,1,0,1,1], [1,0,0,0,1,1], [1,1,1,1,1,1],
        [1,1,0,0,1,1], [1,1,1,1,1,1], [0,1,0,0,1,1], [0,0,1,0,1,1],
        [1,1,0,1,1,1], [1,0,0,0,1,1], [1,0,1,1,1,1], [1,0,0,0,1,1],
        [1,0,1,0,1,1], [1,0,0,0,1,1], [1,1,1,1,1,1], [0,0,0,0,0,0],
        [1,1,1,1,1,1], [1,0,1,0,1,1], [1,0,1,0,1,1], [1,1,1,1,1,1],
        [1,0,0,0,1,1], [1,1,1,1,1,1], [1,0,1,0,1,1], [1,1,1,1,1,1],
    ],
];

// =============================================================================
// AMINO ACID DATA
// data_tables.c lines 312-329
// =============================================================================

// -----------------------------------------------------------------------------
// C: char aapolarity[NAMINOACID+1] = "NNNNPNPPNNPNPPPNNPPP***????";
// data_tables.c line 314
// Note: Added null terminator for C string compatibility
// -----------------------------------------------------------------------------
pub static AAPOLARITY: &[u8; NAMINOACID + 1] = b"NNNNPNPPNNPNPPPNNPPP***????\0";

// -----------------------------------------------------------------------------
// C: char aaletter[NAMINOACID+1] = "FVLICGRSAPTYDHNMWEQK***????";
// data_tables.c line 315
// Note: Added null terminator for C string compatibility
// -----------------------------------------------------------------------------
pub static AALETTER: &[u8; NAMINOACID + 1] = b"FVLICGRSAPTYDHNMWEQK***????\0";

// -----------------------------------------------------------------------------
// C: char aaname[NAMINOACID][20] = { "Phe", "Val", ... };
// data_tables.c lines 316-327
// Note: Each entry is padded to 20 bytes with null terminators
// -----------------------------------------------------------------------------
pub static AANAME: [[u8; 20]; NAMINOACID] = [
    *b"Phe\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 0: Phe
    *b"Val\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 1: Val
    *b"Leu\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 2: Leu
    *b"Ile\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 3: Ile
    *b"Cys\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 4: Cys
    *b"Gly\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 5: Gly
    *b"Arg\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 6: Arg
    *b"Ser\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 7: Ser
    *b"Ala\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 8: Ala
    *b"Pro\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 9: Pro
    *b"Thr\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 10: Thr
    *b"Tyr\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 11: Tyr
    *b"Asp\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 12: Asp
    *b"His\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 13: His
    *b"Asn\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 14: Asn
    *b"Met\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 15: Met
    *b"Trp\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 16: Trp
    *b"Glu\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 17: Glu
    *b"Gln\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 18: Gln
    *b"Lys\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 19: Lys
    *b"Stop\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",   // 20: Stop
    *b"SeC\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 21: SeC (selenocysteine)
    *b"Pyl\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0",  // 22: Pyl (pyrrolysine)
    *b"(Arg|Stop|Ser|Gly)\0\0",                 // 23: ambiguous
    *b"(Ile|Met)\0\0\0\0\0\0\0\0\0\0\0",        // 24: ambiguous
    *b"(Stop|Trp)\0\0\0\0\0\0\0\0\0\0",         // 25: ambiguous
    *b"(Lys|Asn)\0\0\0\0\0\0\0\0\0\0\0",        // 26: ambiguous
];

// -----------------------------------------------------------------------------
// C: char ambig_aaname[4] = "???";
// data_tables.c line 329
// -----------------------------------------------------------------------------
pub static AMBIG_AANAME: [u8; 4] = *b"???\0";

// =============================================================================
// GENETIC CODE TABLES (INDEXED BY ANTICODON)
// data_tables.c lines 331-910
// =============================================================================

// -----------------------------------------------------------------------------
// C: int aamap[NGENECODE][64] = { ... };
// data_tables.c lines 331-910
// Contains 34 genetic code tables (indices 0-33)
// Each table has 64 entries (for each codon, indexed by anticodon)
// Values are amino acid indices defined in common.h (types.rs)
// -----------------------------------------------------------------------------
pub static mut AAMAP: [[i32; 64]; NGENECODE] = [
    // 0. composite metazoan mt
    // data_tables.c lines 333-348
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, 23,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, 24, 25, GLY, ARG, 23,
        SER, ALA, PRO, THR, STOP, GLU, GLN, 26,
    ],
    // 1. standard
    // data_tables.c lines 350-365
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 2. vertebrate mt
    // data_tables.c lines 367-382
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, STOP,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, STOP,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 3. yeast mt
    // data_tables.c lines 384-399
    [
        PHE, VAL, THR, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, THR, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, THR, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, THR, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 4. mold, protozoan, and coelenterate mt
    // data_tables.c lines 401-416
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 5. invertebrate mt
    // data_tables.c lines 418-433
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, SER,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, SER,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 6. ciliate
    // data_tables.c lines 435-450
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLN, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLN, GLU, GLN, LYS,
    ],
    // 7. deleted -> standard
    // data_tables.c lines 452-467
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 8. deleted -> standard
    // data_tables.c lines 469-484
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 9. echinoderm and flatworm mt
    // data_tables.c lines 486-501
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, SER,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, TRP, GLY, ARG, SER,
        SER, ALA, PRO, THR, STOP, GLU, GLN, ASN,
    ],
    // 10. Euplotid
    // data_tables.c lines 503-518
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, CYS, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 11. bacterial and plant chloroplast
    // data_tables.c lines 520-535
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 12. alternate yeast
    // data_tables.c lines 537-552
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, SER, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 13. Ascidian mt
    // data_tables.c lines 554-569
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, GLY,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, GLY,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 14. alternate flatworm mt
    // data_tables.c lines 571-586
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, SER,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, TRP, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, GLU, GLN, ASN,
    ],
    // 15. Blepharisma
    // data_tables.c lines 588-603
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLN, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 16. Chlorophycean mt
    // data_tables.c lines 605-620
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, LEU, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 17. deleted -> standard
    // data_tables.c lines 622-637
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 18. deleted -> standard
    // data_tables.c lines 639-654
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 19. deleted -> standard
    // data_tables.c lines 656-671
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 20. deleted -> standard
    // data_tables.c lines 673-688
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 21. trematode mt
    // data_tables.c lines 690-705
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, SER,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, SER,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 22. Scenedesmus obliquus mt
    // data_tables.c lines 707-722
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, LEU, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        STOP, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 23. Thraustochytrium mt
    // data_tables.c lines 724-739
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, SER, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        STOP, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 24. Pterobranchia mt
    // data_tables.c lines 741-756
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, LYS,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, TRP, GLY, ARG, SER,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 25. Gracilibacteria
    // data_tables.c lines 758-773
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, GLY, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 26. Pachysolen tannophilus
    // data_tables.c lines 775-790
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, ALA, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 27. Karyorelict
    // data_tables.c lines 792-807
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLN, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLN, GLU, GLN, LYS,
    ],
    // 28. Condylostoma
    // data_tables.c lines 809-824
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLN, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLN, GLU, GLN, LYS,
    ],
    // 29. Mesodinium
    // data_tables.c lines 826-841
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, TYR, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, TYR, GLU, GLN, LYS,
    ],
    // 30. Peritrich
    // data_tables.c lines 843-858
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLU, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLU, GLU, GLN, LYS,
    ],
    // 31. Blastocrithidia
    // data_tables.c lines 860-875
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLU, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, GLU, GLU, GLN, LYS,
    ],
    // 32. vacant -> standard
    // data_tables.c lines 877-892
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, ARG,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, SEC, GLY, ARG, ARG,
        SER, ALA, PRO, THR, STOP, GLU, GLN, LYS,
    ],
    // 33. Cephalodiscidae Mitochondrial UAA-Tyr
    // data_tables.c lines 894-909
    [
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, MET, TRP, GLY, ARG, LYS,
        SER, ALA, PRO, THR, PYL, GLU, GLN, LYS,
        PHE, VAL, LEU, ILE, CYS, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, ASP, HIS, ASN,
        LEU, VAL, LEU, ILE, TRP, GLY, ARG, SER,
        SER, ALA, PRO, THR, TYR, GLU, GLN, LYS,
    ],
];

// =============================================================================
// GLOBAL MUTABLE STATE
// data_tables.c lines 912-913
// =============================================================================

// -----------------------------------------------------------------------------
// C: gene *ts;
// data_tables.c line 913
// Pointer to detected genes array (mutable global)
// -----------------------------------------------------------------------------
pub static mut TS: *mut Gene = std::ptr::null_mut();

// =============================================================================
// HELP MENU
// data_tables.c lines 915-1099
// =============================================================================

// -----------------------------------------------------------------------------
// C: char helpmenu[NHELPLINE][81] = { ... };
// data_tables.c lines 916-1099
// Note: Converted to array of string slices for Rust
// -----------------------------------------------------------------------------
pub static HELPMENU: [&str; NHELPLINE] = [
    "-------------------------------------------",
    "Oxagorn v0.1.0 (ARAGORN reimplementation)",
    "-------------------------------------------\n",
    "Please reference the following papers if you use this",
    "program as part of any published research.\n",
    "Laslett, D. and Canback, B. (2004) ARAGORN, a",
    "program for the detection of transfer RNA and transfer-messenger",
    "RNA genes in nucleotide sequences",
    "Nucleic Acids Research, 32;11-16\n",
    "Laslett, D. and Canback, B. (2008) ARWEN: a",
    "program to detect tRNA genes in metazoan mitochondrial",
    "nucleotide sequences",
    "Bioinformatics, 24(2); 172-175.\n",
    "Rust reimplementation by Sunju Kim <n.e.coli.1822@gmail.com>",
    "for bactars bacterial genome annotator.\n\n",
    "Oxagorn detects tRNA, mtRNA, and tmRNA genes.\n",
    "Usage:",
    "oxagorn -v -e -s -d -c -l -j -a -q -rn -w -ifro<min>,<max> -t -mt -m",
    "        -rp -ps -gc -tv -seq -br -fasta -fo -o <outfile> <filename>\n",
    "<filename> is assumed to contain one or more sequences",
    "in FASTA or GENBANK format. Results of the search are printed",
    "to STDOUT. All switches are optional and case-insensitive.",
    "Unless -i is specified, tRNA genes containing introns",
    "are not detected.\n",
    "    -m            Search for tmRNA genes.",
    "    -t            Search for tRNA genes.",
    "                  By default, all are detected. If one of",
    "                  -m or -t is specified, then the other",
    "                  is not detected unless specified as well.",
    "    -mt           Search for Metazoan mitochondrial tRNA genes.",
    "                  tRNA genes with introns not detected. -i,-sr switchs",
    "                  ignored. Composite Metazoan mitochondrial",
    "                  genetic code used.",
    "    -mtmam        Search for Mammalian mitochondrial tRNA",
    "                  genes. -i switch ignored. -tv switch set.",
    "                  Mammalian mitochondrial genetic code used.",
    "    -mtx          Same as -mt but low scoring tRNA genes are",
    "                  not reported.",
    "    -mtd          Overlapping metazoan mitochondrial tRNA genes",
    "                  on opposite strands are reported.",
    "    -gc<num>      Use the GenBank transl_table = <num> genetic code.",
    "    -gcstd        Use standard genetic code.",
    "    -gcmet        Use composite Metazoan mitochondrial genetic code.",
    "    -gcvert       Use Vertebrate mitochondrial genetic code.",
    "    -gcinvert     Use Invertebrate mitochondrial genetic code.",
    "    -gcyeast      Use Yeast mitochondrial genetic code.",
    "    -gcprot       Use Mold/Protozoan/Coelenterate mitochondrial genetic code.",
    "    -gcciliate    Use Ciliate genetic code.",
    "    -gcflatworm   Use Echinoderm/Flatworm mitochondrial genetic code",
    "    -gceuplot     Use Euplotid genetic code.",
    "    -gcbact       Use Bacterial/Plant chloroplast genetic code.",
    "    -gcaltyeast   Use alternative Yeast genetic code.",
    "    -gcascid      Use Ascidian mitochondrial genetic code.",
    "    -gcaltflat    Use alternative flatworm mitochondrial genetic code.",
    "    -gcblep       Use Blepharisma genetic code.",
    "    -gcchloroph   Use Chlorophycean mitochondrial genetic code.",
    "    -gctrem       Use Trematode mitochondrial genetic code.",
    "    -gcscen       Use Scenedesmus obliquus mitochondrial genetic code.",
    "    -gcthraust    Use Thraustochytrium mitochondrial genetic code.",
    "    -gcptero      Use Pterobranchia mitochondrial genetic code.",
    "    -gcgrac       Use Gracilibacteria genetic code.",
    "    -gcpach       Use Pachysolen tannophilus genetic code.",
    "    -gckary       Use Karyorelict genetic code.",
    "    -gccond       Use Condylostoma genetic code.",
    "    -gcmeso       Use Mesodinium genetic code.",
    "    -gcperi       Use Peritrich genetic code.",
    "    -gcblast      Use Blastocrithidia genetic code.",
    "    -gcceph       Use Cephalodiscidae mitochondrial UAA-Tyr genetic code.",
    "                  Individual modifications can be appended using",
    "    ,BBB=<aa>     B = A,C,G, or T. <aa> is the three letter",
    "                  code for an amino-acid. More than one modification",
    "                  can be specified. eg -gcvert,aga=Trp,agg=Trp uses",
    "                  the Vertebrate mitochondrial code and the codons",
    "                  AGA and AGG changed to Tryptophan.",
    "    -c            Assume that each sequence has a circular",
    "                  topology. Search wraps around each end.",
    "                  Default setting.",
    "    -l            Assume that each sequence has a linear",
    "                  topology. Search does not wrap.",
    "    -d            Double. Search both strands of each",
    "                  sequence. Default setting.",
    "    -s  or -s+    Single. Do not search the complementary",
    "                  (antisense) strand of each sequence.",
    "    -sc or -s-    Single complementary. Do not search the sense",
    "                  strand of each sequence.",
    "    -i            Search for tRNA genes with introns in",
    "                  anticodon loop with maximum length 3000",
    "                  bases. Minimum intron length is 0 bases.",
    "                  Ignored if -m is specified.",
    "    -i<max>       Search for tRNA genes with introns in",
    "                  anticodon loop with maximum length <max>",
    "                  bases. Minimum intron length is 0 bases.",
    "                  Ignored if -m is specified.",
    "    -i<min>,<max> Search for tRNA genes with introns in",
    "                  anticodon loop with maximum length <max>",
    "                  bases, and minimum length <min> bases.",
    "                  Ignored if -m is specified.",
    "    -io           Same as -i, but allow tRNA genes with long",
    "                  introns to overlap shorter tRNA genes.",
    "    -if           Same as -i, but fix intron between positions",
    "                  37 and 38 on C-loop (one base after anticodon).",
    "    -ifo          Same as -if and -io combined.",
    "    -ir           Same as -i, but report tRNA genes with minimum",
    "                  length <min> bases rather than search for",
    "                  tRNA genes with minimum length <min> bases.",
    "                  With this switch, <min> acts as an output filter,",
    "                  minimum intron length for searching is still 0 bases.",
    "    -tv           Do not search for mitochondrial TV replacement",
    "                  loop tRNA genes. Only relevant if -mt used.",
    "    -c7           Search for tRNA genes with 7 base C-loops only.",
    "    -ss           Use the stricter canonical 1-2 bp spacer1 and",
    "                  1 bp spacer2. Ignored if -mt set. Default is to",
    "                  allow 3 bp spacer1 and 0-2 bp spacer2, which may",
    "                  degrade selectivity.",
    "    -j            Display 4-base sequence on 3' end of astem",
    "                  regardless of predicted amino-acyl acceptor length.",
    "    -jr           Allow some divergence of 3' amino-acyl acceptor",
    "                  sequence from NCCA.",
    "    -jr4          Allow some divergence of 3' amino-acyl acceptor",
    "                  sequence from NCCA, and display 4 bases.",
    "    -e            Print out score for each reported gene.",
    "    -ps           Lower scoring thresholds to 95% of default levels.",
    "    -ps<num>      Change scoring thresholds to <num> percent of default levels.",
    "    -rp           Flag possible pseudogenes (score < 100 or tRNA anticodon",
    "                  loop <> 7 bases long). Note that genes with score < 100",
    "                  will not be detected or flagged if scoring thresholds are not",
    "                  also changed to below 100% (see -ps switch).",
    "    -rp<num>      Flag possible pseudogenes and change score threshold to <num>",
    "                  percent of default levels.",
    "    -seq          Print out primary sequence.",
    "    -br           Show secondary structure of tRNA gene primary sequence,",
    "                  or tRNA domain for tmRNA genes, using round brackets.",
    "    -svg          Generate SVG image file code for secondary structure.",
    "    -fasta        Print out primary sequence in fasta format.",
    "    -fo           Print out primary sequence in fasta format only",
    "                  (no secondary structure).",
    "    -fon          Same as -fo, with sequence and gene numbering in header.",
    "    -fos          Same as -fo, with no spaces in header.",
    "    -fons         Same as -fo, with sequence and gene numbering, but no spaces.",
    "    -v            Verbose. Prints out information during",
    "                  search to STDERR.",
    "    -a            Print out tRNA domain for tmRNA genes.",
    "    -a7           Restrict tRNA astem length to a maximum of 7 bases",
    "    -aa           Display message if predicted iso-acceptor species",
    "                  does not match species in sequence name (if present).",
    "    -amt<num>     Change annotated tRNA length mismatch reporting threshold to",
    "                  <num> bases when searching GENBANK files. Default is 10 bases.",
    "    -amm<num>     Change annotated tmRNA length mismatch reporting threshold to",
    "                  <num> bases when searching GENBANK files. Default is 30 bases.",
    "    -q            Dont print configuration line (which switches",
    "                  and files were used).",
    "    -rn           Repeat sequence name before summary information.",
    "    -o <outfile>  Print output to <outfile>. If <outfile>",
    "                  already exists, it is overwritten. By default",
    "                  all output goes to stdout.",
    "    -w            Print out in batch mode.",
    "    -wa           Same as -w, but for 6 or 8 base anticodon",
    "                  loops, print possible iso-acceptor species",
    "                  as ?(<species>|<species>) instead of ???",
    "                  For tRNA genes, batch mode output is in the form:\n",
    "                  Sequence name",
    "                  N genes found",
    "                  1 tRNA-<species> [locus 1] <Apos> (nnn)",
    "                  i(<intron position>,<intron length>)",
    "                            .          ",
    "                            .          ",
    "                  N tRNA-<species> [Locus N] <Apos> (nnn)",
    "                  i(<intron position>,<intron length>)\n",
    "                  N is the number of genes found",
    "                  <species> is the tRNA iso-acceptor species",
    "                  <Apos> is the tRNA anticodon relative position",
    "                  (nnn) is the tRNA anticodon base triplet",
    "                  i means the tRNA gene has a C-loop intron\n",
    "                  For tmRNA genes, output is in the form:\n",
    "                  n tmRNA(p) [Locus n] <tag offset>,<tag end offset>",
    "                  <tag peptide>\n",
    "                  p means the tmRNA gene is permuted",
    "    -@<N>         Set number of threads for parallel processing (1-6).",
    "                  -@ or -@0 auto-detect (uses physical cores, max 6).",
    "                  -@1 single-threaded (C original compatible).",
    "    -wunix        Get around problem with some windows gcc compilers",
    "                  (found so far in Strawberry Perl and Active Perl)",
    "                  when reading Unix files.",
    "                  Execution speed may be slower for large files.",
    "                  Execution speed will be a lot slower for files",
    "                  with many small sequences.",
];

// =============================================================================
// tmRNA tag database
// Note: tmRNA tag database is defined in tmrna_tags.rs
// data_tables.c line 1101
// =============================================================================
