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

#include "sdr.h"

#ifndef min
    #define min(x,y) ((x>y) ? y : x)
    #define max(x,y) ((x>y) ? x : y)
#endif

// -- CONTROL BLOCK AND SPE RANK
static context_block   cb              _CLINE_ALIGN;
static int             m_rank;
static int             m_mode;

// -- IF- AND BASEBAND DATA as floats
static vector float ifband[2][IF_VECSPERBLOCK]     _CLINE_ALIGN;
static vector float baseband[2][IF_VECSPERBLOCK]   _CLINE_ALIGN;


// ------------------------------------------------------------------------------------------------
//    D M A   H E L P E R S
// ------------------------------------------------------------------------------------------------

#define TAG_CB          1 // control block DMA access tag
#define TAG_IFBAND      2 // IF band data input tag
#define TAG_BASEBAND    3 // baseband data output tag
#define TAG_A           4 // generic tags
#define TAG_B           5
#define tag2mask(x) (1<<x)

inline void load_data(addr64 _PPE_basePtr, void* _field, int _blocknr, size_t blocksize, int _tag) {
    mfc_get((float*)_field,
            _PPE_basePtr.ull + (unsigned long long)(_blocknr * blocksize),
            blocksize, _tag, 0, 0);
}

inline void store_data(addr64 _PPE_basePtr, void* _feld, int _blocknr, size_t blocksize, int _tag) {
    mfc_putf((float*)_feld,
            _PPE_basePtr.ull + (unsigned long long)(_blocknr * blocksize),
            blocksize, _tag, 0, 0);
}

inline void wait_for_tag(int _tag) {
    mfc_write_tag_mask(tag2mask(_tag));
    mfc_read_tag_status_all();
    return;
}

inline void wait_for_mask(int _mask) {
    mfc_write_tag_mask(_mask);
    mfc_read_tag_status_all();
    return;
}


// ------------------------------------------------------------------------------------------------
//    O T H E R   H E L P E R S
// ------------------------------------------------------------------------------------------------

inline void print4f(vector float f) {
    printf("%f %f %f %f ", spu_extract(f,0), spu_extract(f,1),spu_extract(f,2),spu_extract(f,3));
}


// ------------------------------------------------------------------------------------------------
//    M A I N
// ------------------------------------------------------------------------------------------------

int main(unsigned long long speid, addr64 argp, addr64 envp) {

    unsigned int i;
    unsigned int t_start = 0, t_spu = 0;
    unsigned int count16kB = 0;
    unsigned int blocknr;
    unsigned long* argc = (unsigned long*)&speid;
    Oscillator oscLO;
    FltDeg8 lowpass, bandpass, highpass;
    Analyzer ana;

    vector float lp_iir_A[2] = {
       (vector float){ 3.369674E-4, 0.0023587719, 0.0070763156, 0.011793859 },
       (vector float){ 0.011793859, 0.0070763156, 0.0023587719, 3.369674E-4 }
    };
    vector float lp_iir_B[2] = {
       (vector float){ 1.0, -3.4815187, 5.7111516, -5.5168867 },
       (vector float){ 3.3460727, -1.2626357, 0.27289245, -0.02594379 }
    };
    vector float lp_fir_B[2] = {
       (vector float){ -0.07406144, -0.029266886, 0.11199167, 0.2721453 },
       (vector float){ 0.34239945, 0.2721453, 0.11199167, -0.029266886 }
    };

    /* check mode - SPE thread started by main PPU program, or from console */
    if (argc[0] != 1) {

        /* started by main PPU prog */
        m_mode = MODE_EMBEDDED;

        /* fetch the control block via DMA */
        mfc_get(&cb, argp.ui[1], sizeof(cb), TAG_CB, 0, 0);
        mfc_write_tag_mask(1<<TAG_CB);
        mfc_read_tag_status_all();

        /* initialize rank, locks/incrementers */
        m_rank = cb.rank;

        /* tell the PPE that the SPE thread is running */
        spu_write_out_mbox(m_rank);

        /* get the first raw data set */
        load_data(cb.ifband_src, (void*)ifband[0], 0, IF_BYTEPERBLOCK, TAG_IFBAND);

    } else {

        /* started by user */
        m_mode = MODE_CONSOLE;
        m_rank = 0;
        cb.blockscount = 4096;

        printf("Console mode, not really implemented yet...\n");

    }

    /* init some tests */
    Oscillator_initialize(&oscLO, 0.0, 1.0 /* f,MHz */ , 64.0 /* fs,MHz */);
    FltDeg8_initialize(&lowpass, lp_iir_A, lp_fir_B);
    Analyzer_reset(&ana);

    /* start to measure execution time */
    spu_write_decrementer(0x7fffffff);
    t_start = spu_read_decrementer();

    /* do whatever maths should be done */
    for (blocknr = 0; blocknr < cb.blockscount; blocknr++ )
    {
        Oscillator_moreComplexes(&oscLO, ifband[0], ifband[1], /*sampsDiv16:*/ (4*IF_VECSPERBLOCK/16)); /* 3.26 Gfloat/s */
        // FltDeg8_doFIR_I(&lowpass, ifband[0], IF_VECSPERBLOCK); /* 0.406887 Gfloat/s */
        // Oscillator_moreReals(&oscLO, ifband[0], /*sampsDiv16:*/ (4*IF_VECSPERBLOCK/16)); /*  1.75 Gfloat/s, not too well parallized */
        Analyzer_processReal(&ana, ifband[0], IF_VECSPERBLOCK);

        count16kB+=1;
    }

    /* stop the SPE-side measurement of the execution time */
    t_spu = t_start - spu_read_decrementer();

    /* results dump */
    if (MODE_CONSOLE == m_mode) {
        for (blocknr = 0; blocknr < IF_VECSPERBLOCK/32; blocknr++) {
            print4f(ifband[0][blocknr]);
        }
        printf("\nMin %f, max %f, average %f (variance %f)\n", 
               ana.min_peak.re, ana.max_peak.re, ana.average.re, ana.variance.re);
    }

    /* return or display statistics */
    if (MODE_EMBEDDED==m_mode) {
        spu_write_out_mbox(t_spu);
        spu_write_out_mbox(count16kB);
    } else {
        printf("Time delta: %fs\n", (double)t_spu/((double)__timebase__));
        printf("16kB count: %u\n", count16kB);
        printf("Throughput: %f Gfloat/s\n", 
                 (0.25 * 16384.0 / (1024.0*1024.0*1024.0)) * ((double)count16kB) 
                 * ((double)__timebase__) / ((double)t_spu));
    }

    return 0;
}



// ------------------------------------------------------------------------------------------------
//    C A L C U L A T I O N   H E L P E R S
// ------------------------------------------------------------------------------------------------

// Single 4-tap FIR block
inline vector float FIRtap4block(vector float oldx, vector float newx, vector float b) {
    // b: [b0 b1 b2 b3]
    vector float acc;
    vector float xshifted[4];
    xshifted[0] = oldx;
    xshifted[1] = vec_sld(oldx, newx, 4);
    xshifted[2] = vec_sld(oldx, newx, 8);
    xshifted[3] = vec_sld(oldx, newx, 12);
            // instead of extract+splat perhaps "spu_shuffle(b, b, SHREP_VEC_ELEM(0..3))"
    acc = spu_madd(spu_splats(spu_extract(b,0)), xshifted[0], VEC_4f(0.0));
    acc = spu_madd(spu_splats(spu_extract(b,1)), xshifted[1], acc);
    acc = spu_madd(spu_splats(spu_extract(b,2)), xshifted[2], acc);
    acc = spu_madd(spu_splats(spu_extract(b,3)), xshifted[3], acc);
    return acc;
}



// ------------------------------------------------------------------------------------------------
//    C L A S S E S
// ------------------------------------------------------------------------------------------------


// =======================
// FIR FILTER OF 8th ORDER
// =======================

// Initialize the 8th order filter with coefficients and reset the history
//   Coeffs layout is A[] = [a0 a1 a2 a3 | a4 a5 a6 a7], B[] likewise
//   Copies coeffs into a[0]=[a0 a0 a0 a0], a[1]=[a1 a1 a1 a1], ...
void FltDeg8_initialize(FltDeg8* flt, vector float A[2], vector float B[2]) {
    int i;
    flt->compactA[0] = A[0]; flt->compactA[1] = A[1];
    flt->compactB[0] = B[0]; flt->compactB[1] = B[1];
    for (i=0; i<4; i++) {
        flt->a[0+i] = spu_splats(spu_extract(A[0], i));
        flt->a[4+i] = spu_splats(spu_extract(A[1], i));
        flt->b[0+i] = spu_splats(spu_extract(B[0], i));
        flt->b[4+i] = spu_splats(spu_extract(B[1], i));
    }
    flt->xhist0 = VEC_4f(0.0); flt->xhist1 = VEC_4f(0.0);
    flt->yhist0 = VEC_4f(0.0); flt->yhist1 = VEC_4f(0.0);
    return;
}

// Filter the input data through an 8-tap FIR filter
// Has an inplace version and a separate input/output array version.
void FltDeg8_doFIR_I(FltDeg8* flt, vector float* inout, int vectors) {
    FltDeg8_doFIR(flt, inout, inout, vectors);
}
void FltDeg8_doFIR(FltDeg8* flt, vector float* in, vector float* out, int vectors) {
    int i;
#if 0
    for (i=0; i<vectors/(2 * /*unroll:*/ 4); i++) {
        // perf: 0.406887 Gfloat/s
        vector float xnew0 = *(in+0);
        vector float xnew1 = *(in+1);
        vector float acc0 = (VEC_4f(0.0));
        vector float acc1 = (VEC_4f(0.0));
        #undef FIR_DO_COEFF
        #define FIR_DO_COEFF(n) \
            acc0 = spu_madd(flt->b[n], flt->xhist0, acc0); \
            acc1 = spu_madd(flt->b[n], flt->xhist1, acc1); \
            flt->xhist0 = vec_sld(flt->xhist0, flt->xhist1, 4); \
            flt->xhist1 = vec_sld(flt->xhist1, xnew0, 4); \
            xnew0  = vec_sld(xnew0, xnew1, 4); \
            xnew1  = spu_slqwbyte(xnew1, 4)
        // process 2 input vectors (8 floats)
        UNROLL_BY_4 ( \
            FIR_DO_COEFF(0); \
            FIR_DO_COEFF(1); \
            FIR_DO_COEFF(2); \
            FIR_DO_COEFF(3); \
            FIR_DO_COEFF(4); \
            FIR_DO_COEFF(5); \
            FIR_DO_COEFF(6); \
            FIR_DO_COEFF(7); \
            *(out+0) = acc0; \
            *(out+1) = acc1; \
            in += 2; \
            out += 2; \
        );
    }
#else
    for (i=0; i<vectors/(2 * /*partitioning*/ 4); i++) {
        // TODO, perhaps faster since more parallelism
        // perf: 1.988649 Gfloat/s
        vector float acc00, acc01, acc10, acc11;
        acc00 = FIRtap4block(flt->xhist0, flt->xhist1, flt->compactB[0]);
        acc01 = FIRtap4block(*in, *(in+1), flt->compactB[0]);
        acc10 = FIRtap4block(flt->xhist0, flt->xhist1, flt->compactB[1]);
        acc11 = FIRtap4block(*in, *(in+1), flt->compactB[1]);
        *(out+0) = spu_add(acc00, acc10);
        *(out+1) = spu_add(acc01, acc11);
    }
#endif
    return;
}



// =======================
// IIR FILTER OF 8th ORDER
// =======================

// Filter the input data through an 8-tap IIR filter
// Has an inplace version and a separate input/output array version.
void FltDeg8_doIIR_I(FltDeg8* flt, vector float* inout, int vectors) {
    FltDeg8_doIIR(flt, inout, inout, vectors);
}
void FltDeg8_doIIR(FltDeg8* flt, vector float* in, vector float* out, int vectors) {
    int i;
// TODO: not quite correct yet
    for (i=0; i<vectors/(2 * /*unroll:*/ 4); i++) {
        vector float xnew0 = *(in+0);
        vector float xnew1 = *(in+1);
        vector float ynew0, ynew1;
        vector float acc0 = (VEC_4f(0.0));
        vector float acc1 = (VEC_4f(0.0));
        #undef FIR_DO_COEFF
        #define FIR_DO_COEFF(n) \
            acc0 = spu_madd(flt->b[n], flt->xhist0, acc0); \
            acc1 = spu_madd(flt->b[n], flt->xhist1, acc1); \
            flt->xhist0 = vec_sld(flt->xhist0, flt->xhist1, 4); \
            flt->xhist1 = vec_sld(flt->xhist1, xnew0, 4); \
            xnew0  = vec_sld(xnew0, xnew1, 4); \
            xnew1  = spu_slqwbyte(xnew1, 4)
        #define IIR_DO_YCOEFF(n) \
            acc0 = spu_madd(flt->a[n], flt->yhist0, acc0); \
            acc1 = spu_madd(flt->a[n], flt->yhist1, acc1); \
            flt->yhist0 = vec_sld(flt->yhist0, flt->yhist1, 4); \
            flt->yhist1 = vec_sld(flt->yhist1, ynew0, 4); \
            ynew0  = vec_sld(ynew0, ynew1, 4); \
            ynew1  = spu_slqwbyte(ynew1, 4)
        // process 2 input vectors (8 floats)
        UNROLL_BY_4 ( \
            FIR_DO_COEFF(0); \
            FIR_DO_COEFF(1); \
            FIR_DO_COEFF(2); \
            FIR_DO_COEFF(3); \
            FIR_DO_COEFF(4); \
            FIR_DO_COEFF(5); \
            FIR_DO_COEFF(6); \
            FIR_DO_COEFF(7); \
            *(out+0) = acc0; \
            *(out+1) = acc1; \
            ynew0 = acc0; \
            ynew1 = acc1; \
            IIR_DO_YCOEFF(0); \
            IIR_DO_YCOEFF(1); \
            IIR_DO_YCOEFF(2); \
            IIR_DO_YCOEFF(3); \
            IIR_DO_YCOEFF(4); \
            IIR_DO_YCOEFF(5); \
            IIR_DO_YCOEFF(6); \
            IIR_DO_YCOEFF(7); \
            *(out+0) = spu_add(*(out+0), acc0); \
            *(out+1) = spu_add(*(out+1), acc1); \
            in += 2; \
            out += 2; \
        );
    }
    return;
}



// =====================
// QUADRATURE OSCILLATOR
// =====================

// Oscillator initialization and seeding
// The oscillator is based on trigonometric identity of
//   sin(x+inc) = sin(x)*cos(inc) + cos(x)*sin(inc)
//   cos(x+inc) = cos(x)*cos(inc) - sin(x)*sin(inc)
//   TODO: Numeric drift is still a problem!
// Could use digital waveguide osc
//   http://ccrma.stanford.edu/~jos/pasp/Digital_Waveguide_Oscillator.html
void Oscillator_initialize(Oscillator* osc, float phaseOffset, float freq, float samplerate) {
    vector float ph_inc = spu_splats((float) ((2*3.1415926535897932384626) * freq/samplerate));
    osc->phase    = spu_madd(((vector float){ 0.0, 1.0, 2.0, 3.0 }),/*'*'*/ ph_inc,  /*'+'*/ spu_splats(phaseOffset));
    osc->phaseinc = spu_mul(VEC_4f(4.0), ph_inc);
    #ifdef OSC_COUPLED
    // Coupled oscillators
    osc->sinXpY = sin18_v(osc->phase);
    osc->cosXpY = cos18_v(osc->phase);
    osc->sin4Y  = sin18_v(osc->phaseinc);
    osc->cos4Y  = cos18_v(osc->phaseinc);
    #else
    // Direct form oscillator
    osc->oldsin[1] = sin18_v(osc->phase);
    osc->oldcos[1] = cos18_v(osc->phase);
    osc->oldsin[0] = sin18_v(spu_add(osc->phase, osc->phaseinc));
    osc->oldcos[0] = cos18_v(spu_add(osc->phase, osc->phaseinc));
    osc->twoCosOmega = spu_mul(VEC_4f(2.0), cos18_v(osc->phaseinc));
    #endif
    return;
}

// Get more complex samples out of the oscillator
void Oscillator_moreComplexes(Oscillator* osc, vector float* out_re, vector float* out_im, int sampsDiv16) {
    int i;
    #ifdef OSC_COUPLED
    // Coupled oscillators
    for (i=0; i<sampsDiv16; i++) {
        vector float nsin, ncos;
        UNROLL_BY_4( \
            *(out_re++) = osc->cosXpY; \
            *(out_im++) = osc->sinXpY; \
            nsin = spu_madd (osc->sinXpY, osc->cos4Y, spu_mul(osc->cosXpY, osc->sin4Y)); \
            ncos = spu_msub (osc->cosXpY, osc->cos4Y, spu_mul(osc->sinXpY, osc->sin4Y)); \
            osc->sinXpY = nsin; \
            osc->cosXpY = ncos; \
        );
    }
    /* osc->phases += 16.0 * sampsDiv16 * osc->phaseinc; NEEDS wrapping! */
    /* TODO recalculate correct values with sin18_v()/cos..() after some iterations to avoid drift */
    #else
    // Direct form oscillator
    for (i=0; i<sampsDiv16; i++) {
        vector float nsin, ncos;
        UNROLL_BY_4( \
            nsin = spu_msub (osc->twoCosOmega, osc->oldsin[0], osc->oldsin[1]); \
            ncos = spu_msub (osc->twoCosOmega, osc->oldcos[0], osc->oldcos[1]); \
            *(out_re++) = osc->oldcos[1]; \
            *(out_im++) = osc->oldsin[1]; \
            osc->oldcos[1] = osc->oldcos[0]; \
            osc->oldsin[1] = osc->oldsin[0]; \
            osc->oldcos[0] = ncos; \
            osc->oldsin[0] = nsin; \
        );
    }
    #endif
}

// Get more real-valued samples out of the oscillator
void Oscillator_moreReals(Oscillator* osc, vector float* out, int sampsDiv16) {
    int i;
    #ifdef OSC_COUPLED
    // Coupled oscillators
    for (i=0; i<sampsDiv16; i++) {
        vector float nsin, ncos;
        UNROLL_BY_4( \
            *(out++) = osc->sinXpY; \
            nsin = spu_madd (osc->sinXpY, osc->cos4Y, spu_mul(osc->cosXpY, osc->sin4Y)); \
            ncos = spu_msub (osc->cosXpY, osc->cos4Y, spu_mul(osc->sinXpY, osc->sin4Y)); \
            osc->sinXpY = nsin; \
            osc->cosXpY = ncos; \
        );
    }
    /* osc->phases += 16.0 * sampsDiv16 * osc->phaseinc; NEEDS wrapping! */
    /* TODO recalculate correct values with sin18_v()/cos..() after some iterations to avoid drift */
    #else
    // Direct form oscillator
    for (i=0; i<sampsDiv16; i++) {
        vector float nsin;
        UNROLL_BY_4( \
            nsin = spu_msub (osc->twoCosOmega, osc->oldsin[0], osc->oldsin[1]); \
            *(out++) = osc->oldsin[1]; \
            osc->oldsin[1] = osc->oldsin[0]; \
            osc->oldsin[0] = nsin; \
        );
    }
    #endif
}



// ======================
// SIMPLE SIGNAL ANALYZER
// ======================

void Analyzer_reset(Analyzer* ana) {
    vector float *overwrite = (vector float*)&ana->min_peak;
    *overwrite = VEC_4f(0.0);     // complex min_peak, max_peak
    *(overwrite+1) = VEC_4f(0.0); // complex average, variance
}

void Analyzer_processReal(Analyzer* ana, vector float* in, int vectors) {
    int chunks = vectors / 16;
    int leftovers = vectors - 16*chunks;
    vector float minpeak = VEC_4f(FLT_MAX);
    vector float maxpeak = VEC_4f(FLT_MIN);
    vector float ave = VEC_4f(0.0), var = VEC_4f(0.0);
    vector unsigned int c0, c1;
    vector float n = VEC_4f(0.0);
    vector float tmp, vec;
    int i;
    for (i=0; i<chunks; i++) {
       UNROLL_BY_16( \
       vec = *(in++); \
       /* min max */ \
       c0 = spu_cmpgt(vec, maxpeak); \
       maxpeak = spu_sel(maxpeak, vec, c0); \
       c1 = spu_cmpgt(minpeak, vec); \
       minpeak = spu_sel(minpeak, vec, c1); \
       /* ave : my[n] = (1/n) ( (n - 1)*my[n-1] + x ) -- not parallel enough yet */ \
       tmp = spu_madd(n, ave, vec); \
       n = spu_add(n, VEC_4f(4.0)); \
       ave = spu_mul(spu_re(n), tmp); \
       );
    }
    ana->average.re = spu_extract(ave,0) + spu_extract(ave,1) + spu_extract(ave,2) + spu_extract(ave,3);
    ana->min_peak.re = min( min(spu_extract(minpeak,0), spu_extract(minpeak,1)),
                         min(spu_extract(minpeak,2), spu_extract(minpeak,3)) );
    ana->max_peak.re = max( max(spu_extract(maxpeak,0), spu_extract(maxpeak,1)),
                         max(spu_extract(maxpeak,2), spu_extract(maxpeak,3)) );
    // TODO: recursive variance
    return;
}

void Analyzer_processComplex(Analyzer* ana, vector float* re, vector float* im, int vectors);
