/********************************************************************************
 * IBM Cell Miniature 6-station 1024-point correlator
 * Copyright (C) 2007 Jan Wagner, Metsahovi Radio Observatory, Finland
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 ********************************************************************************/
#ifndef _MINICORRELATOR_H
#define _MINICORRELATOR_H

#include <stdlib.h>
#include <math.h>

#ifdef __SPU__
    #include <simdmath.h>
    #include <libgmath.h>
    #include <fft_1d.h>
    #include <fft_1d_r2.h>
    #include <spu_mfcio.h>
    //#define SEPARATE_RE_IM 1
    //#include "fft_1d_r2_modified.h"
#else
    #include <libspe.h>
    // #include <libspe2.h>
    // #include <pthread.h>
    #include <stdio.h>
    #include <sched.h>
    #include <sys/wait.h>
    #include <malloc_align.h>
    #include <free_align.h>
    #include <sys/time.h>
    #include <string.h>
#endif

// Alignment helpers
#define _CLINE_ALIGN     __attribute__ ((aligned (128)))
#define _QUAD_ALIGN      __attribute__ ((aligned (16)))

// Vector and optimization helpers
#define UNROLL_BY_1(x)  { x }
#define UNROLL_BY_2(x)  { x }{ x }
#define UNROLL_BY_4(x)  { x }{ x }{ x }{ x }
#define UNROLL_BY_8(x)  UNROLL_BY_4(x) UNROLL_BY_4(x)
#define UNROLL_BY_16(x) UNROLL_BY_8(x) UNROLL_BY_8(x)

#define REPLICATE_2(x,y) { y=0; { x } y=1; { x } }
#define REPLICATE_4(x,y) { y=0; { x } y=1; { x } y=2; { x } y=3; { x } }

#define VEC_4f(x)       (vector float){x, x, x, x}

// Branch hints for the compiler
#ifndef unlikely
#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)
#endif

// Struct for 32 vs 64-bit address passing between PPU and SPU 
typedef union {
  unsigned long long ull;
  unsigned int ui[2];
} addr64;

// Compile modes
#define USE_COMBINED_3BASELINE_CALCS 1
#define USE_COUPLED_OSCILLATOR       0

// Executing time profiling
#define __freq__      3200     // CPU frequency in MHz (cat /proc/cpuinfo)
#define __timebase__  79800000 // Cell Blade: 14318000  PS3: 79800000 (cat /proc/cpuinfo)
#define start_timer(pts) { spu_write_decrementer(0xffffffff); *pts=spu_read_decrementer(); }
#define stop_timer(pts)  { *pts -= spu_read_decrementer(); }

// All things hard-coded :
#define NUM_SPE_THREADS 6
#define NUM_BASELINES   (((NUM_SPE_THREADS-1)*NUM_SPE_THREADS)/2)
#define BASELNS_ON_SPE  ((NUM_BASELINES/NUM_SPE_THREADS) + 1)

#define SPE_RAW_BYTES        16384                 // raw data buffer
#define SAMPLES_PER_BYTE     4                     // 4 floats for 1 byte of 2-bit data

#define SPE_FIXEDBUFSIZE     16384                 // bytes
#define SPE_FIXEDFLOATSIZE   (SPE_FIXEDBUFSIZE/4)  // floats     : 4096
#define SPE_FIXEDCMPLXSIZE   (SPE_FIXEDBUFSIZE/8)  // complexs   : 2048
#define SPE_FIXEDVECSIZE     (SPE_FIXEDBUFSIZE/16) // vectors    : 1024
#define N_FFT                1024                  // 1024pt FFT
#define N_FFT_RVECS          (N_FFT/4)             // 1024pt FFT : 256 vecs real, 256 vecs imag

#define INTEGRATE_NUM_FFTS  PPU_RAW_TOTAL_FFTS      // how many FFT:s to integrate over
#define PPU_CORRELBUF_BYTES (N_FFT*2*sizeof(float)) // number of bytes in single auto- and cross-correlation result
#define SPU_CORRELBUF_VECS  (2*N_FFT_RVECS)         // how many float vectors the auto/crosscorrelation bufs consist of
#define PPU_RAW_BUFCOUNT      256                   // numer of SPE_RAW_BYTES data buffers reserved on PPU side for each of 6 SPE's
#define PPU_RAW_1SPEBYTES   (PPU_RAW_BUFCOUNT*SPE_RAW_BYTES) // source raw data bytes in total, for one SPE
#define PPU_RAW_TOTAL_FFTS  ((PPU_RAW_1SPEBYTES * SAMPLES_PER_BYTE) / N_FFT)

// -- HARD-CODED 6 STATION EXPERIMENT
typedef struct _control_block {
    addr64 rawdata_src          _CLINE_ALIGN;
    addr64 autocorrelation_out;
    addr64 baselines_out[15];
    addr64 fmultipliers_src;
    addr64 syncline;
    int spe_num;
    int num_spes;
    unsigned int ffts_total;
    unsigned int ffts_to_integrate;
    addr64 spe_ls_listOnPPU;
    addr64 spe_sig_listOnPPU;
} control_block                 _CLINE_ALIGN;

#endif
