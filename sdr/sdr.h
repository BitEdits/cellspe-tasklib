/********************************************************************************
 * IBM Cell Software Defined Radio
 * Copyright (C) 2007 Jan Wagner
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 ********************************************************************************/
#ifndef _SDR_H
#define _SDR_H

#include <stdlib.h>
#include <math.h>
#ifdef __SPU__
    #include <simdmath.h>
    #include <libgmath.h>
    //#include <fft_1d.h>
    //#include <fft_1d_r2.h>
    #include <spu_mfcio.h>
    #include <stdio.h>
    #include <sin18_v.h>
    #include <float.h>
#else
    #include <libspe.h>
    #include <stdio.h>
    #include <sched.h>
    #include <sys/wait.h>
    #include <malloc_align.h>
    #include <free_align.h>
    #include <sys/time.h>
#endif


// ------------------------------------------------------------------------------------------------
//    D E F I N E S
// ------------------------------------------------------------------------------------------------

// Branch hints
#ifndef unlikely
#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)
#endif

// Executing time profiling
#define __freq__      3200     // CPU frequency in MHz (cat /proc/cpuinfo)
#define __timebase__  79800000 // Cell Blade: 14318000  PS3: 79800000 (cat /proc/cpuinfo)

// Unrolling helpers
#define REPLICATE_2(x,y) { y=0; { x } y=1; { x } }
#define REPLICATE_4(x,y) { y=0; { x } y=1; { x } y=2; { x } y=3; { x } }
#define UNROLL_BY_2(x)  { x }{ x }
#define UNROLL_BY_4(x)  { x }{ x }{ x }{ x }
#define UNROLL_BY_8(x)  UNROLL_BY_4(x) UNROLL_BY_4(x)
#define UNROLL_BY_16(x) UNROLL_BY_8(x) UNROLL_BY_8(x)

// Generic vector helpers
#define _CLINE_ALIGN     __attribute__ ((aligned (128)))
#define _QUAD_ALIGN      __attribute__ ((aligned (16)))
#define VEC_4f(x)          (vector float){x, x, x, x}
#define SHREP_VEC_ELEM(x)  (vector unsigned char) { \
                                0+4*(x),1+4*(x),2+4*(x),3+4*(x), \
                                0+4*(x),1+4*(x),2+4*(x),3+4*(x), \
                                0+4*(x),1+4*(x),2+4*(x),3+4*(x), \
                                0+4*(x),1+4*(x),2+4*(x),3+4*(x) }
#define vec_sld(A,B,bytes)  spu_shuffle(A, B, \
     ( (vector unsigned char){0+(bytes),  1+(bytes),  2+(bytes),  3+(bytes), \
                              4+(bytes),  5+(bytes),  6+(bytes),  7+(bytes), \
                              8+(bytes),  9+(bytes), 10+(bytes), 11+(bytes), \
                             12+(bytes), 13+(bytes), 14+(bytes), 15+(bytes)} ) )

// Alignment helpers
#define _CLINE_ALIGN     __attribute__ ((aligned (128)))
#define _QUAD_ALIGN      __attribute__ ((aligned (16)))

// Data buffers
#define IF_VECSPERBLOCK 1024                 // 1kV
#define IF_BYTEPERBLOCK (16*IF_VECSPERBLOCK) // 16kB, real input data

// SPE and PPU run time related
#define NUM_SPE_THREADS    1     // how many SPE threads the PPU prog should start
#define MODE_CONSOLE       1     // SPE thread was started from console
#define MODE_EMBEDDED      2     // SPE thread was started from PPU program

// ------------------------------------------------------------------------------------------------
//    S T R U C T S
// ------------------------------------------------------------------------------------------------

// Address conversion 32/64-bit
typedef union {
    unsigned long long ull;
    unsigned int ui[2];
    void* p;
} addr64;

// Task context
typedef struct _context_block {
    addr64 ifband_src   _CLINE_ALIGN;
    addr64 demod_dst;
    int    rank;
    int    blockscount;
} context_block         _CLINE_ALIGN;

// Complex nums
typedef struct _complex_t {
    float re, im;
} complex_t;


// ------------------------------------------------------------------------------------------------
//    "C L A S S E S"
// ------------------------------------------------------------------------------------------------

#ifdef __SPU__

// Not using C++ classes after all, some SDK2.x libraries and headers
// have trouble with C++. Templates could've been rather nice... :-/

// "Class" for sine or sine/cosine (IQ) oscillator
typedef struct _Oscillator_t {
    vector float sin4Y,  cos4Y;
    vector float sinXpY, cosXpY;
    vector float phase,  phaseinc;
    vector float oldsin[2], oldcos[2], twoCosOmega;
} Oscillator;
void Oscillator_initialize(Oscillator* osc, float phaseOffset, float freq, float samplerate);
void Oscillator_moreComplexes(Oscillator* osc, vector float* out_re, vector float* out_im, int sampsDiv16);
void Oscillator_moreReals(Oscillator* osc, vector float* out, int sampsDiv16);

// "Class" for 8th degree FIR/IIR filter
typedef struct _FltDeg8_t {
    // http://www.cl.cam.ac.uk/~al407/research/abstracts/acaces06.pdf
    // http://www.freescale.com/webapp/sps/site/taxonomy.jsp?nodeId=02VS0l81285Nf2F9DH "AltiVec Code Samples"
    // http://teleinfo.pb.edu.pl/~krashan/altivec/G_Kraszewski_XI_Symposium_AES_Bialystok_2006.pdf
    vector float a[8], b[8];
    vector float compactA[2], compactB[2];
    vector float xhist0, xhist1, yhist0, yhist1;
} FltDeg8;
void FltDeg8_initialize(FltDeg8* flt, vector float A[2], vector float B[2]);
void FltDeg8_doFIR_I(FltDeg8* flt, vector float* inout, int vectors);
void FltDeg8_doFIR  (FltDeg8* flt, vector float* in, vector float* out, int vectors);
void FltDeg8_doIIR_I(FltDeg8* flt, vector float* inout, int sampsDiv16);
void FltDeg8_doIIR  (FltDeg8* flt, vector float* in, vector float* out, int vectors);

// "Class" for signal properties analyzer
typedef struct _Analyzer_t {
    complex_t min_peak, max_peak, average, variance;
} Analyzer;
void Analyzer_reset(Analyzer* ana);
void Analyzer_processReal(Analyzer* ana, vector float* in, int vectors);
void Analyzer_processComplex(Analyzer* ana, vector float* re, vector float* im, int vectors);

#endif // __SPU__

#endif // _SDR_H
