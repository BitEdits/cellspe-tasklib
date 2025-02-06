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

//
// SPE debugging -- http://www.ibm.com/developerworks/power/library/pa-celldebug/
//   $ sudo bash
//   $ SPU_DEBUG_START=1 ./minicorrelator &
//   $ spu-gdb minicorrelator_spu -p [processID]
//

#include "minicorrelator.h"


// -- CONTROL BLOCK AND SPE RANK
control_block   cb               _CLINE_ALIGN;
addr64          spe_ls_Ptrs[16]  _CLINE_ALIGN;
addr64          spe_sig_Ptrs[16] _CLINE_ALIGN;
int             g_rank;

// -- STATISTICS
unsigned int totaldmascount;
unsigned int rawblocknr;

// -- HARD-CODED 2-bit SAMPLE LOOKUP FOR BYTES CONTAINING |Samp0[bit0 bit1] Samp1[bit0 bit1] Samp2[bit0 bit1] Samp3[bit0 bit1]|
vector float    byte_to_vecfloat[256]; // 2-bit samples to floats lookup (1 byte to 4 float samples), uninitialized for now...

// -- HARD-CODED 6 STATION EXPERIMENT, 1024-point FFT (re&im 2048 floats = 512 vectors = 8kB)
vector unsigned char raw_in[SPE_RAW_BYTES/16]          _CLINE_ALIGN; // was: raw_in[SPE_FIXEDVECSIZE] i.e. 16kB, now just 1kB
vector float    stationdata0[SPE_FIXEDVECSIZE]         _CLINE_ALIGN; // timeslot 0 unpacked and station-processed data, 2 x 1024pt complex
vector float    stationdata1[SPE_FIXEDVECSIZE]         _CLINE_ALIGN; // timeslot 1 unpacked and station-processed data, 2 x 1024pt complex
vector float    autocorr[SPU_CORRELBUF_VECS]           _CLINE_ALIGN; // time-integrated station autocorrelation, 2 x 1024pt complex
vector float    baselines[BASELNS_ON_SPE][SPU_CORRELBUF_VECS] _CLINE_ALIGN; // time-integrated baseline cross-correlation
vector float    stationdata[6][SPE_FIXEDVECSIZE]       _CLINE_ALIGN; // station data tbat all other SPEs have DMA'ed in

vector float    *bl0Pair, *bl1Pair, *bl2Pair;                        // pointers to stationdata[x] this SPE should MAC into
unsigned int    spreadout_targetSPEs[NUM_SPE_THREADS][BASELNS_ON_SPE] =
                { {2, 3, 5}, {0, 3, 4}, {1, 4, 5}, {0, 2, 5}, {0, 1, 3}, {1, 2, 4} };

volatile vector unsigned short stationdata_locks[8]    _CLINE_ALIGN; // only first ..locks[0] is used, rest is padding
vector float    fmultipliers[SPE_FIXEDVECSIZE/4]       _CLINE_ALIGN; // 4kB fits one 1024-point complex FFT, these are real multipliers
vector float    fftTwiddleF[N_FFT_RVECS / 2]           _CLINE_ALIGN; // FFT twiddle factors

// -- INLINE CALC FUNCS WITH HARD-CODED 1024-point FFT AND 16kB BLOCK SETTINGS
inline void sincos_contiguous(float startangle, float angleinc, vector float * in1real, vector float * in1imag, vector float * in2real, vector float * in2imag);
inline void sincos_fast(vector float* phase, vector float* mulwithsin, vector float* mulwithcos);
void        more_raws_to_complex(vector float* re_im);
inline void complex_multiply_accumulate(vector float* ReA, vector float* ImA, vector float* ReB, vector float* ImB, vector float* ReAcc, vector float* ImAcc);
inline void complex_multiply_accumulate_autocorrelation(vector float* A, vector float* Acc);
void        three_baseline_naive_MACauto(vector float* center, vector float* A, vector float* B, vector float* C);

// -- OTHER HELPER FUNCS, WITH 1024-point FFT 16kB BLOCK HARD-CODED
inline void check_baseline_integration_period(unsigned int timeperiod);
void zero_baselines_and_autocorr();


// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
//
//    D M A   H E L P E R S
//
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------

#define TAG_CB      1 // control block DMA access tag
#define TAG_RAW     2 // raw data input tag
#define TAG_LOCKS   3 // stationdata_locks[] access
#define TAG_SELFDMA 4 // inside SPE DMA
#define TAG_A       5 // generic tags
#define TAG_B       6
#define TAG_SIGNALS 7
#define TAG_PPUOUT  8
#define tag2mask(x) (1<<x)

inline void load_data(addr64 _PPE_basePtr, float *_field, int _blocknr, size_t blocksize, int _tag) {
    mfc_get( _field,
            _PPE_basePtr.ull + (unsigned long long)(_blocknr * blocksize), 
            blocksize, 
            _tag, 0, 0);
}

inline void store_data(addr64 _PPE_basePtr, float *_field, int _blocknr, size_t blocksize, int _tag) {
    mfc_putf( _field,
            _PPE_basePtr.ull + (unsigned long long)(_blocknr * blocksize), 
            blocksize, 
            _tag, 0, 0);
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

inline void send_signal(addr64 sigAddr) {
    unsigned int signal[4] _CLINE_ALIGN;
    unsigned long long ea = sigAddr.ull + 12;
    char *ls = ((char*)&signal[0]) + 12;
    signal[3] = (1 << g_rank);
    mfc_sndsig(ls, ea, TAG_SIGNALS, 0, 0);
    mfc_write_tag_mask(tag2mask(TAG_SIGNALS));
    mfc_read_tag_status_all();
}

inline int receive_signal(){
    while(!spu_stat_signal1());
    return spu_read_signal1();
}

void barrier(){
    int i, temp=1;
    unsigned int check;

    if(g_rank!=0) {
        // signal SPU0 that we're done, then wait for SPE0 to wake us
        send_signal(spe_sig_Ptrs[0]);
        receive_signal();
    } else {
        check = 0;
        while(check < ((1<<NUM_SPE_THREADS)-1)) {
            check |= spu_read_signal1();
        }
    }
    for(i=0;i<3; i++){
        if(g_rank<temp) {
            send_signal(spe_sig_Ptrs[g_rank+temp]);
        }
        temp=temp<<1;
    }
}

// Increment SPE's own entry in the state vector array
// inc=1 -> final=3, inc 2 -> final=6, inc 3 -> final=9, inc 4 -> final=12 : something odd going on
void atomic_statevec_inc(int _tag)
{
    unsigned int status, i;
    vector unsigned short vec_inc = spu_insert( (unsigned short)1,
                                                ((vector unsigned short){ 0, 0, 0, 0, 0, 0, 0, 0 }),
                                                g_rank );
    do {
        // Read the PPU-side cacheline containing the state vector
        mfc_getllar(stationdata_locks, cb.syncline.ull, 0, 0);
        status = mfc_read_atomic_status();
        // Increment it
        stationdata_locks[0] = spu_add(stationdata_locks[0], vec_inc);
        // Conditionally (lock not lost) put it back
        mfc_putllc(stationdata_locks, cb.syncline.ull, 0, 0);
        status = mfc_read_atomic_status();
        // Reattempt if the lock was lost
    } while (status & MFC_PUTLLC_STATUS);
    return;
}

// Write the stationdata0 processed data (of station 'g_rank') out to
// every SPE's stationdata[g_rank] buffer
void spread_out_stationdata(vector float* sourcebuf, int _tag) {
    int i;
    addr64 out_ls;
    // LS-addresses, same on all SPEs as it is the same binary after all
    unsigned long long out_flt_offset  = (unsigned long)&stationdata[g_rank];
    for (i=0; i<BASELNS_ON_SPE; i++) {
        int k = spreadout_targetSPEs[g_rank][i];
        mfc_put(sourcebuf, spe_ls_Ptrs[k].ull + out_flt_offset, SPE_FIXEDBUFSIZE, _tag, 0, 0);
    }
    // mark this iteration as 'done' on own SPU (no overflow check...)
    atomic_statevec_inc(TAG_LOCKS);
    return;
}

// Wait first for own DMAs from spead_out_stationdata0() to
// finish. Then wait for other SPUs to have incremented the 
// stationdata "lock" counter (=dma finished) to the value 
// of this SPEs lock
void join_otherSPUs_stationdata(int _tag)
{
    unsigned int oldmask, events, status;
    unsigned int ea_lo, ea_hi;
    volatile vector unsigned short* pl = &stationdata_locks[0];
    vector unsigned short expected = spu_splats(spu_extract(*pl,g_rank));
    vector unsigned short equal;
    wait_for_tag(_tag);
    // Setup for possible use of lock line reservation lost events.
    // Detect and discard phantom events.
    oldmask = spu_read_event_mask();
    spu_write_event_mask(0);
    if (spu_stat_event_status()) {
        spu_write_event_ack(spu_read_event_status());
    }
    spu_write_event_mask(MFC_LLR_LOST_EVENT);
    while (1) {
        if (spu_stat_event_status()) {
            spu_write_event_ack(spu_read_event_status());
        }
        // Reserve cache line, get the sync lock
        mfc_getllar(stationdata_locks, cb.syncline.ull, 0, 0);
        status = mfc_read_atomic_status();
        equal = spu_cmpeq(expected, stationdata_locks[0]);
        if(spu_extract(spu_gather(equal), 0) == 0xFF)
        {
            break;
        } else {
            // Wait until LLAR is lost before checking again.            
        break; // somehow "atomic" is broken or otherwise LLR_LOST_EVENT is never received!!
            events = spu_read_event_status(); // 'events' should be just MFC_LLR_LOST_EVENT
            spu_write_event_ack(events);
        }
    }
    // Clean up any remnant events, restore the event mask
    spu_write_event_mask(0);
    if (spu_stat_event_status()) {
        spu_write_event_ack(spu_read_event_status());
    }
    spu_write_event_mask(oldmask);
    return;
}


// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
//
//    M A I N 
//
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------

int main(unsigned long long speid, addr64 argp, addr64 envp) {

    int i;
    unsigned int timeperiod = 0;
    unsigned int t_start = 0, t_spu;


    /* tell the PPE that the SPE thread is running */
    spu_write_out_mbox(0);

    /* this mailbox call blocks until the PPE has created the control block */
    spu_read_in_mbox();

    /* start the SPE-side measurement of the execution time */
    start_timer(&t_spu);

    /* fetch the control block via DMA */
    mfc_get(&cb, argp.ull, sizeof(cb), TAG_CB, 0, 0);
    mfc_write_tag_mask(1<<TAG_CB);
    mfc_read_tag_status_all();

    /* pre-fetch the list of other SPU LS virtual EA's */
    mfc_get(spe_ls_Ptrs,  cb.spe_ls_listOnPPU.ull,  sizeof(spe_ls_Ptrs),  0, 0, TAG_CB);
    mfc_get(spe_sig_Ptrs, cb.spe_sig_listOnPPU.ull, sizeof(spe_sig_Ptrs), 0, 0, TAG_CB);

    /* initialize rank, locks/incrementers */
    g_rank = cb.spe_num;

#if 0
    wait_for_tag(TAG_CB);
    barrier();
    /*
    for (rawblocknr = 0, timeperiod = 0; rawblocknr < 8; ) {
        //atomic_statevec_inc(TAG_LOCKS);
        //join_otherSPUs_stationdata(TAG_LOCKS);
        rawblocknr += 1;
    }
    */
#else

    /* pre-fetch first raw data set and the delay variables */
    load_data(cb.rawdata_src, (float*)raw_in, 0, SPE_RAW_BYTES, TAG_RAW);
    rawblocknr = 1;
    load_data(cb.fmultipliers_src, (float*)fmultipliers, 0, SPE_FIXEDBUFSIZE/4, TAG_CB);

    /* precalculate some things */
    // -- assign 3 stations to xcorrelate together with this SPE's own station data (hard-coded 6 station setup...)
    // no DMA optimizations here...!
    bl0Pair = &stationdata[(g_rank+1)%NUM_SPE_THREADS][0];
    bl1Pair = &stationdata[(g_rank+3)%NUM_SPE_THREADS][0];
    bl2Pair = &stationdata[(g_rank+4)%NUM_SPE_THREADS][0];
    // -- reset baseline's data
    zero_baselines_and_autocorr();
    // -- FFT twiddle factors for 1024-point FFT
    {
        typedef struct { float re; float im; } fc32;
        fc32* fft_W = (fc32*)&fftTwiddleF[0];
        fft_W->re = 1.0;
        fft_W->im = 0.0;
        for (i=1; i<N_FFT/4; i++) {
            float arg;
            // (fft_W+i)->re = cos(i * 2*M_PI/N_FFT);
            arg = i * 2*M_PI/N_FFT;
            (fft_W+i)->re = spu_extract(cos14_v(spu_splats(arg)), 0);
            (fft_W + N_FFT/4 - i)->im = -(fft_W + i)->re;
        }
    }
    // -- TODO initialize 2-bit sample to float lookup table properly
    for (i=0; i<256; i++) {
        byte_to_vecfloat[i] = spu_mul(spu_splats((float)i), ((vector float){0.1,0.25,0.15,0.05}));
    }

    /* wait for spe_ls_Ptrs[] dma'ed in */
    wait_for_tag(TAG_CB);

    // ------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------

    /* do the station and baseline processing for all of the raw data from the PPU-side */
    for (timeperiod = 0; timeperiod < cb.ffts_total; )
    {

        // ---- STATION PROCESSING 0 START
        // ---- STATION PROCESSING 0 START
        // ---- STATION PROCESSING 0 START

        // --- UNPACK data into two 1024-pt complex bufs
        more_raws_to_complex(stationdata0);
        more_raws_to_complex(stationdata0+2*N_FFT_RVECS);

        // --- SINCOS phase rotation for 2 x (1024-point complexnum blocks)
        sincos_contiguous(/* startangle */ 0.1, /* angleinc */ 0.23,
            /* in1real */ stationdata0,     /* in1imag */ stationdata0+N_FFT_RVECS,
            /* in2real */ stationdata0+512, /* in2imag */ stationdata0+768);

#if 0
        // --- FFT : 1024-point FFT needs 2048 floats (re&im), is 8kB, so two FFTs per 16kB block
        // TODO: inplace FFT, needs [re im re im] vectors!!
        // _fft_1d_r2(vector float *out, vector float *in, vector float *W, int log2_size)
        _fft_1d_r2(stationdata0,               stationdata0,               fftTwiddleF, 9 /* 2^10=1024 */);
        _fft_1d_r2(stationdata0+2*N_FFT_RVECS, stationdata0+2*N_FFT_RVECS, fftTwiddleF, 9 /* 2^10=1024 */);
#endif


        // --- SINCOS fast and rather accurate
        sincos_fast(fmultipliers, stationdata0+N_FFT_RVECS /* *sin() */, stationdata0     /* *cos() */);
        sincos_fast(fmultipliers, stationdata0+512         /* *sin() */, stationdata0+768 /* *cos() */);

        // --- DONE : dump own stationprocessed data to all other 5 SPEs
        spread_out_stationdata(stationdata0, TAG_LOCKS);

        // --- Autocorrelation
        #if !USE_COMBINED_3BASELINE_CALCS
        complex_multiply_accumulate_autocorrelation(stationdata0,               autocorr);
        complex_multiply_accumulate_autocorrelation(stationdata0+2*N_FFT_RVECS, autocorr);
        #endif


        // ---- STATION PROCESSING 1 START
        // ---- STATION PROCESSING 1 START
        // ---- STATION PROCESSING 1 START

        // --- UNPACK more data into two 1024-pt complex bufs
        more_raws_to_complex(stationdata1);
        more_raws_to_complex(stationdata1+2*N_FFT_RVECS);

        // --- SINCOS phase rotation for 2 x (1024-point complexnum blocks)
        sincos_contiguous(/* startangle */ 0.1, /* angleinc */ 0.23,
            /* in1real */ stationdata1    , /* in1imag */ stationdata1+N_FFT_RVECS,
            /* in2real */ stationdata1+512, /* in2imag */ stationdata1+768);

#if 0
        // --- FFT : 1024-point FFT needs 2048 floats (re&im), is 8kB, so two FFTs per 16kB block
        // TODO: inplace FFT, needs [re im re im] vectors!!
        // _fft_1d_r2(vector float *out, vector float *in, vector float *W, int log2_size)
        _fft_1d_r2(stationdata1,               stationdata1,               fftTwiddleF, 10 /* 2^10=1024 */);
        _fft_1d_r2(stationdata1+2*N_FFT_RVECS, stationdata1+2*N_FFT_RVECS, fftTwiddleF, 10 /* 2^10=1024 */);
#endif

        // --- SINCOS fast and rather accurate
        sincos_fast(fmultipliers, stationdata1+256 /* *sin() */, stationdata1     /* *cos() */);
        sincos_fast(fmultipliers, stationdata1+512 /* *sin() */, stationdata1+768 /* *cos() */);

        // --- Autocorrelation
        #if !USE_COMBINED_3BASELINE_CALCS
        complex_multiply_accumulate_autocorrelation(stationdata1,               autocorr);
        complex_multiply_accumulate_autocorrelation(stationdata1+2*N_FFT_RVECS, autocorr);
        #endif

        // --- DONE : dump own stationprocessed data to all other 5 SPEs
        // after baseline 0 done


        // ---- BASELINE 0 START
        // ---- BASELINE 0 START
        // ---- BASELINE 0 START

        // wait for others to have sent their stationdata0
        join_otherSPUs_stationdata(TAG_LOCKS);

        #if USE_COMBINED_3BASELINE_CALCS
        three_baseline_naive_MACauto(stationdata0, bl0Pair, bl1Pair, bl2Pair);
        three_baseline_naive_MACauto(stationdata0 + 2*N_FFT_RVECS, bl0Pair + 2*N_FFT_RVECS,
                                 bl1Pair + 2*N_FFT_RVECS, bl2Pair + 2*N_FFT_RVECS);
        #else
        // xcorrelate first 1024-point complex data sets
        complex_multiply_accumulate(stationdata0, stationdata0+N_FFT_RVECS,
                                    bl0Pair,      bl0Pair+N_FFT_RVECS,
                                    baselines[0], baselines[0]+N_FFT_RVECS);
        complex_multiply_accumulate(stationdata0, stationdata0+N_FFT_RVECS,
                                    bl1Pair,      bl1Pair+N_FFT_RVECS,
                                    baselines[1], baselines[1]+N_FFT_RVECS);
        complex_multiply_accumulate(stationdata0, stationdata0+N_FFT_RVECS,
                                    bl2Pair,      bl2Pair+N_FFT_RVECS,
                                    baselines[2], baselines[2]+N_FFT_RVECS);
        // xcorrelate second 1024-point complex data sets
        complex_multiply_accumulate(stationdata0+2*N_FFT_RVECS, stationdata0+3*N_FFT_RVECS,
                                    bl0Pair+2*N_FFT_RVECS,      bl0Pair+3*N_FFT_RVECS,
                                    baselines[0],               baselines[0]+N_FFT_RVECS);
        complex_multiply_accumulate(stationdata0+2*N_FFT_RVECS, stationdata0+3*N_FFT_RVECS,
                                    bl1Pair+2*N_FFT_RVECS,      bl1Pair+3*N_FFT_RVECS,
                                    baselines[1],               baselines[1]+N_FFT_RVECS);
        complex_multiply_accumulate(stationdata0+2*N_FFT_RVECS, stationdata0+3*N_FFT_RVECS,
                                    bl2Pair+2*N_FFT_RVECS,      bl2Pair+3*N_FFT_RVECS,
                                    baselines[2],               baselines[2]+N_FFT_RVECS);
        #endif

        // check if one integration period done
        timeperiod += 2;
        check_baseline_integration_period(timeperiod);

        // dispatch previos station data
        spread_out_stationdata(stationdata1, TAG_LOCKS);


        // ---- BASELINE 1 START
        // ---- BASELINE 1 START
        // ---- BASELINE 1 START

        // wait for others to have sent their stationdata1
        join_otherSPUs_stationdata(TAG_LOCKS);

        #if USE_COMBINED_3BASELINE_CALCS
        three_baseline_naive_MACauto(stationdata1, bl0Pair, bl1Pair, bl2Pair);
        three_baseline_naive_MACauto(stationdata1 + 2*N_FFT_RVECS, bl0Pair + 2*N_FFT_RVECS,
                                 bl1Pair + 2*N_FFT_RVECS, bl2Pair + 2*N_FFT_RVECS);
        #else
        // xcorrelate first 1024-point complex data sets
        complex_multiply_accumulate(stationdata1, stationdata1+N_FFT_RVECS,
                                    bl0Pair,      bl0Pair+N_FFT_RVECS,
                                    baselines[0], baselines[0]+N_FFT_RVECS);
        complex_multiply_accumulate(stationdata1, stationdata1+N_FFT_RVECS,
                                    bl1Pair,      bl1Pair+N_FFT_RVECS,
                                    baselines[1], baselines[1]+N_FFT_RVECS);
        complex_multiply_accumulate(stationdata1, stationdata1+N_FFT_RVECS,
                                    bl2Pair,      bl2Pair+N_FFT_RVECS,
                                    baselines[2], baselines[2]+N_FFT_RVECS);
        // xcorrelate second 1024-point complex data sets
        complex_multiply_accumulate(stationdata1+2*N_FFT_RVECS, stationdata1+3*N_FFT_RVECS,
                                    bl0Pair+2*N_FFT_RVECS,      bl0Pair+3*N_FFT_RVECS,
                                    baselines[0],               baselines[0]+N_FFT_RVECS);
        complex_multiply_accumulate(stationdata1+2*N_FFT_RVECS, stationdata1+3*N_FFT_RVECS,
                                    bl1Pair+2*N_FFT_RVECS,      bl1Pair+3*N_FFT_RVECS,
                                    baselines[1],               baselines[1]+N_FFT_RVECS);
        complex_multiply_accumulate(stationdata1+2*N_FFT_RVECS, stationdata1+3*N_FFT_RVECS,
                                    bl2Pair+2*N_FFT_RVECS,      bl2Pair+3*N_FFT_RVECS,
                                    baselines[2],               baselines[2]+N_FFT_RVECS);
        #endif


        // send data back to RAM after integration period 
        timeperiod += 2;
        check_baseline_integration_period(timeperiod);

    } // timeperiod = 0..cb.ffts_total

    // ------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------

    /* wait until all DMAs finished */
    wait_for_mask(tag2mask(TAG_LOCKS));

#endif
    /* stop the SPE-side measurement of the execution time */
    stop_timer(&t_spu);

    spu_write_out_mbox(t_spu);
    spu_write_out_mbox(rawblocknr);

    return 0;
}


// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------


// ------------------------------------------------------------------------------------------------
// Continguous phase rotation using quadrature oscillator, for two subsequent 1024-pt complex
// Runs as either coupled oscillator (slow, drifts) or direct form IIR filter (fast, drifts less)
//
//
// in: phi, delta_phi    inout: 1024pt 4kB/256V re1,im1, re2,im2
//                              total inout 2 x 1024pt complex float 16kB
//
inline void sincos_contiguous(float startangle, float angleinc, vector float * in1real, vector float * in1imag,
vector float * in2real, vector float * in2imag) {

    int i;
    #if USE_COUPLED_OSCILLATOR
    vector float sinvec, cosvec;
    vector float sin4Y, cos4Y, argvec, oldsin, oldcos, anginc4;

    argvec  = spu_madd(spu_splats(angleinc), ((vector float) { 0.0, 1.0, 2.0, 3.0 }),
                            /* + */ spu_splats(startangle) );
    anginc4 = spu_mul(spu_splats(angleinc), ((vector float) { 4.0, 4.0, 4.0, 4.0 }) );
    sin4Y   = sin14_v(anginc4);
    cos4Y   = cos14_v(anginc4);
    sinvec  = sin14_v(argvec); // use cont. phase for both blocks below
    cosvec  = cos14_v(argvec);
    #else
    vector float oldsin[2], oldcos[2], twoCosOmega, nsin, ncos;
    vector float phase1inc = spu_splats(angleinc);
    vector float phase     = spu_madd(((vector float){ 0.0, 1.0, 2.0, 3.0 }), phase1inc, spu_splats(startangle));
    vector float phaseinc  = spu_mul(VEC_4f(4.0), phase1inc);
    oldsin[1] = sin18_v(phase);
    oldcos[1] = cos18_v(phase);
    oldsin[0] = sin18_v(spu_add(phase, phaseinc));
    oldcos[0] = cos18_v(spu_add(phase, phaseinc));
    twoCosOmega = spu_mul(VEC_4f(2.0), cos18_v(phaseinc));
    #endif

    #if USE_COUPLED_OSCILLATOR
    // first 1024-point complex data
    for (i=0; i<(N_FFT_RVECS/8); i++) {
        // 2 complex nums in each 128-bit vector, process 4 nums each iteration
        UNROLL_BY_8( \
            *in1real = spu_mul(*in1real, sinvec); \
            *in1imag = spu_mul(*in1imag, cosvec); \
            oldsin = sinvec; \
            sinvec = spu_madd(oldsin, cos4Y, spu_mul(cosvec, sin4Y)); \
            oldcos = cosvec; \
            cosvec = spu_nmsub(oldsin, sin4Y, spu_mul(oldcos, cos4Y)); \
            in1real++; in1imag++; \
        );
    }
    // second 1024-point complex data
    for (i=0; i<(N_FFT_RVECS/8); i++) {
        // 2 complex nums in each 128-bit vector, process 4 nums each iteration
        UNROLL_BY_8( \
            *in2real = spu_mul(*in2real, sinvec); \
            *in2imag = spu_mul(*in2imag, cosvec); \
            oldsin = sinvec; \
            sinvec = spu_madd(oldsin, cos4Y, spu_mul(cosvec, sin4Y)); \
            oldcos = cosvec; \
            cosvec = spu_nmsub(oldsin, sin4Y, spu_mul(oldcos, cos4Y)); \
            in2real++; in2imag++; \
        );
    }
    #else
    // first 1024-point complex data
    for (i=0; i<(N_FFT_RVECS/4); i++) {
        UNROLL_BY_4( \
            *in1real  = spu_mul(*in1real, oldsin[1]); \
            *in1imag  = spu_mul(*in1imag, oldcos[1]); \
            nsin      = spu_msub(twoCosOmega, oldsin[0], oldsin[1]); \
            ncos      = spu_msub(twoCosOmega, oldcos[0], oldcos[1]); \
            oldcos[1] = oldcos[0]; \
            oldsin[1] = oldsin[0]; \
            oldcos[0] = ncos; \
            oldsin[0] = nsin; \
            in1real++; in1imag++; \
        );
    }

    // second 1024-point complex data
    for (i=0; i<(N_FFT_RVECS/4); i++) {
        UNROLL_BY_4( \
            *in2real  = spu_mul(*in2real, oldsin[1]); \
            *in2imag  = spu_mul(*in2imag, oldcos[1]); \
            nsin      = spu_msub(twoCosOmega, oldsin[0], oldsin[1]); \
            ncos      = spu_msub(twoCosOmega, oldcos[0], oldcos[1]); \
            oldcos[1] = oldcos[0]; \
            oldsin[1] = oldsin[0]; \
            oldcos[0] = ncos; \
            oldsin[0] = nsin; \
            in1real++; in1imag++; \
        );
    }
    #endif
    return;
}


// ------------------------------------------------------------------------------------------------
// Multiply the value in *mulwithsin with sin(*phase) and *mulwithcos with cos(*phase)
// in: 1024pt float phase 4kB/256V    inout: 1024pt re float 4kB/256V, 1024pt im float 4kB/256V
//                                           total inout 1024pt complex float 8kB
// phase must be -pi .. +pi, NORMALIZED into -1.0 .. +1.0 !!
//
inline void sincos_fast(vector float* phase, vector float* mulwithsin, vector float* mulwithcos) {

    // return; // delta: 0.188763 sec went to 0.143851
    vector float stmp, ctmp, phi2, phi;
    // vector signed int iphase;
    int i;
    for (i=0; i<(N_FFT_RVECS/8); i++) { // SPE_FFT_RE_VECS/16 = vecs/unrollfactor
        /* wrap argument into -1.0..+1.0 using integer truncation (hmm this isn't quite proper.. +1->0 etc... , */
        //    phi = *(phase++);
        //    iphase = spu_convts(phi, 0);
        //    phi = spu_sub(phi, spu_convtf(iphase, 0));
        //    phi2 = spu_mul(phi, phi);
        /* Matlab coeffs:
            phi=linspace(-pi, pi, 64000);
            phinorm = phi ./ pi;
            trueSine = sin(phi);
            trueCosine = cos(phi);
            coeffsS = polyfit(phinorm, trueSine, 12);
            coeffsC = polyfit(phinorm, trueCosine, 12);
            figure(1), plot(phinorm, polyval(coeffsS,phinorm)-trueSine, 'g-');
            figure(2), plot(phinorm, polyval(coeffsC,phinorm)-trueCosine, 'g-');
        */
        UNROLL_BY_8( \
            phi = *(phase++); \
            /* iphase = spu_convts(phi, 0); */ \
            /* phi = spu_sub(phi, spu_convtf(iphase, 0)); */ \
            phi2 = spu_mul(phi, phi); \
            stmp = spu_madd(phi2, VEC_4f(0.080605269719834), VEC_4f(-0.006041231531886)); \
            stmp = spu_madd(stmp, phi2, VEC_4f(-0.598397824003969)); \
            stmp = spu_madd(stmp, phi2, VEC_4f(2.549926789216792));  \
            stmp = spu_madd(stmp, phi2, VEC_4f(-5.167685041675656)); \
            stmp = spu_madd(stmp, phi2, VEC_4f(3.141591733073902));  \
            stmp = spu_mul(stmp, phi); \
            *mulwithsin = spu_mul(stmp, *mulwithsin); /* imag*sin(phase) */ \
            mulwithsin++; \
            ctmp = spu_madd(phi2, VEC_4f(-0.025391123907667), VEC_4f(0.001605367647616)); \
            ctmp = spu_madd(ctmp, phi2, VEC_4f(0.235063383537898));  \
            ctmp = spu_madd(ctmp, phi2, VEC_4f(-1.335174458310010)); \
            ctmp = spu_madd(ctmp, phi2, VEC_4f(4.058698263116778));  \
            ctmp = spu_madd(ctmp, phi2, VEC_4f(-4.934801388410942)); \
            ctmp = spu_madd(ctmp, phi2, VEC_4f(0.999999992289180));  \
            *mulwithcos = spu_mul(ctmp, *mulwithcos); /* real*cos(phase) */  \
            mulwithcos++; \
        );
    }
    return;
}


// ------------------------------------------------------------------------------------------------
// Unpack raw data to vector floats 
// in: raw serial 2-bit data 1kB/64V   out: 4096 floats 16kB/1kV
//
void more_raws_to_complex(vector float* re_im) {
    static vector unsigned char *rawsrc = raw_in;
    vector unsigned char toComplexA = (vector unsigned char) { 0,1,2,3,     0x10,0x10,0x10,0x10,
                                                               4,5,6,7,     0x10,0x10,0x10,0x10, };
    vector unsigned char toComplexB = (vector unsigned char) { 8,9,10,11,   0x10,0x10,0x10,0x10,
                                                               12,13,14,15, 0x10,0x10,0x10,0x10, };
    vector float *re, *im;
    vector float tmp;
    int i;

    /* do unpacking into N_FFT_RVECS real vectors [re re re re], zero out or duplicate to the imag part */
    #define RAWBCONV_R(x) \
        tmp=byte_to_vecfloat[spu_extract(*rawsrc,x)];  *(re++)=tmp; *(im++)=tmp; /*=(VEC_4f(0.0))*/
    re = re_im;
    im = re_im + N_FFT_RVECS;
    wait_for_tag(TAG_RAW);
    for (i=0; i<(N_FFT_RVECS/16); i++) { // 16 new vecs per iter, 1 byte into 4 floats = 1 vector
        RAWBCONV_R(0);  RAWBCONV_R(1);  RAWBCONV_R(2);  RAWBCONV_R(3);
        RAWBCONV_R(4);  RAWBCONV_R(5);  RAWBCONV_R(6);  RAWBCONV_R(7);
        RAWBCONV_R(8);  RAWBCONV_R(9);  RAWBCONV_R(10); RAWBCONV_R(11);
        RAWBCONV_R(12); RAWBCONV_R(13); RAWBCONV_R(14); RAWBCONV_R(15);
        rawsrc++;
    }
    /* prefetch more data when required */
    if (rawsrc >= &raw_in[SPE_RAW_BYTES/16]) {
        rawsrc = raw_in;
        load_data(cb.rawdata_src, (float*)raw_in, rawblocknr, SPE_RAW_BYTES, TAG_RAW);
        rawblocknr++;
    }
    return;
}


// ------------------------------------------------------------------------------------------------
// Multiply the complex data sets (A1re,A11im), (B1re,B1im), then accumulate into (accre,accim)
// in: Re are 4kB/1024F/256V, Im are same size     inout: same size Re,Im as input
//     total 2 x 8kB                                      total 8kB
//
inline void complex_multiply_accumulate(vector float* ReA, vector float* ImA,
                                        vector float* ReB, vector float* ImB,
                                        vector float* ReAcc, vector float* ImAcc) {

    //return; // delta: 0.188763 sec went to 0.111209
    int i;
    for(i=0; i<(N_FFT_RVECS/8); i++) {
        // (ReA*ReB - ImA*ImB) + i(ImA*ReB + ImB*ReA)
        #if 0 // direct version
        vector float acctmp1, acctmp2;
        UNROLL_BY_8( \
            acctmp1 = spu_msub(*ReA, *ReB, spu_mul(*ImA, *ImB)); \
            *ReAcc  = spu_add(*ReAcc, acctmp1); \
            acctmp2 = spu_madd(*ImA, *ReB, spu_mul(*ImB, *ReA)); \
            *ImAcc  = spu_add(*ImAcc, acctmp2); \
            ReAcc++; ImAcc++; \
            ReA++; ReB++; ImA++; ImB++; \
        );
        #else // IBM version, possibly better interleaving
        // (a*c - b*d) + i(a*d + b*c) = (ReA*ReB - ImA*ImB) + i(ReA*ImB + ImA*ReB)
        vector float neg_bd, bc, adPbc, acSbd;
        UNROLL_BY_8( \
            neg_bd = spu_nmsub(*ImA, *ImB, VEC_4f(0.0)); \
            bc     = spu_madd (*ImA, *ReB, VEC_4f(0.0)); \
            adPbc  = spu_madd(*ReA, *ImB, bc); \
            acSbd  = spu_madd(*ReA, *ReB, neg_bd); \
            *ReAcc = spu_add(*ReAcc, acSbd); \
            *ImAcc = spu_add(*ImAcc, adPbc); \
            ReAcc++; ImAcc++; \
            ReA++; ReB++; ImA++; ImB++; \
        );
        #endif
    }
    return;
}


// ------------------------------------------------------------------------------------------------
// Multiply one base and three other complex datasets and accumulate the products,
// also calculate the autocorrelation simultaneously.
//
// The accumulated results are placed into the 'baselines[]' arrays and into 'autocorr'.
//
// Input and output data sets are vector float* ptr -> [N_FFT_RVECSx[Re Re Re Re] N_FFT_RVECSx[Im Im Im Im]]
//
void three_baseline_naive_MACauto(vector float* center, vector float* A, vector float* B, vector float* C) {
    vector float *ReBase = center, *ImBase = center + N_FFT_RVECS;
    vector float *ReA = A, *ImA = A + N_FFT_RVECS;
    vector float *ReB = B, *ImB = B + N_FFT_RVECS;
    vector float *ReC = C, *ImC = C + N_FFT_RVECS;
    vector float *ReB0 = baselines[0], *ImB0 = baselines[0] + N_FFT_RVECS;
    vector float *ReB1 = baselines[1], *ImB1 = baselines[1] + N_FFT_RVECS;
    vector float *ReB2 = baselines[2], *ImB2 = baselines[2] + N_FFT_RVECS;
    vector float *ReAuto = autocorr,   *ImAuto = autocorr + N_FFT_RVECS;
    int i;
    for(i=0; i<(N_FFT_RVECS/2); i++) {
        vector float neg_bd[3], bc[3], adPbc[3], acSbd[3];
        UNROLL_BY_2( \
            /* autocorr */ \
            *ReAuto   = spu_madd(*ReBase, *ReBase, spu_madd(*ImBase, *ImBase, *ReAuto)); \
            *ImAuto   = (VEC_4f(0.0)); \
            /* correlate base with 3 others */ \
            neg_bd[0] = spu_nmsub(*ImBase, *ImA, VEC_4f(0.0)); \
            neg_bd[1] = spu_nmsub(*ImBase, *ImB, VEC_4f(0.0)); \
            neg_bd[2] = spu_nmsub(*ImBase, *ImC, VEC_4f(0.0)); \
            bc[0]     = spu_madd (*ImBase, *ReA, VEC_4f(0.0)); \
            bc[1]     = spu_madd (*ImBase, *ReB, VEC_4f(0.0)); \
            bc[2]     = spu_madd (*ImBase, *ReC, VEC_4f(0.0)); \
            adPbc[0]  = spu_madd (*ReBase, *ImA, bc[0]); \
            adPbc[1]  = spu_madd (*ReBase, *ImA, bc[1]); \
            adPbc[2]  = spu_madd (*ReBase, *ImA, bc[2]); \
            acSbd[0]  = spu_madd (*ReBase, *ReA, neg_bd[0]); \
            acSbd[1]  = spu_madd (*ReBase, *ReB, neg_bd[1]); \
            acSbd[2]  = spu_madd (*ReBase, *ReC, neg_bd[2]); \
            *ReB0 = spu_add(*ReB0, acSbd[0]); \
            *ReB1 = spu_add(*ReB1, acSbd[1]); \
            *ReB2 = spu_add(*ReB2, acSbd[2]); \
            *ImB0 = spu_add(*ImB0, adPbc[0]); \
            *ImB1 = spu_add(*ImB1, adPbc[1]); \
            *ImB2 = spu_add(*ImB2, adPbc[2]); \
            ReB0++; ImB0++; ReB1++; ImB1++; ReB2++; ImB2++; \
            ReA++; ImA++; ReB++; ImB++; ReC++; ImC++; \
            ReAuto++; ImAuto++; \
            ReBase++; ImBase++; \
        );
    }
}


// ------------------------------------------------------------------------------------------------
// Multiply the complex data set (A1re,A1im) with its conjugate, then accumulate into (accre,accim)
//
// in: Re are 4kB/1024F/256V, Im are same size     inout: same size Re,Im as input
//     total 8kB                                          total 8kB
//
inline void complex_multiply_accumulate_autocorrelation(vector float* A, vector float* Acc) {
    vector float *ReA = A, *ImA = A + N_FFT_RVECS;
    vector float *ReAcc = Acc, *ImAcc = Acc + N_FFT_RVECS;

    // return; // delta: 0.188763 sec went to 0.160620
    int i;
    for(i=0; i<(N_FFT_RVECS/16); i++) {
        // (Re+iIm)*(Re-iIm) = (Re^2 + Im^2) + i0
        UNROLL_BY_16( \
            *ReAcc = spu_madd(*ReA, *ReA, spu_madd(*ImA, *ImA, *ReAcc)); \
            *ImAcc = (VEC_4f(0.0)); \
            ReAcc++; ImAcc++; \
            ReA++; ImA++; \
        );
    }
    return;
}


// ------------------------------------------------------------------------------------------------
// Set baseline data back to 0, reset autocorrelation
//
void zero_baselines_and_autocorr() {

    // return; // delta: 0.188763 sec went to 0.188773
    int i;
    vector float *bl0 = baselines[0], *bl1 = baselines[1], *bl2 = baselines[2];
    vector float *a0  = autocorr;
    for (i=0; i<SPU_CORRELBUF_VECS; i++) {
        *bl0 = VEC_4f(0.0); *bl1 = VEC_4f(0.0); *bl1 = VEC_4f(0.0); *a0  = VEC_4f(0.0);
        bl0++; bl1++; bl2++; a0++;
    }
    return;
}


// ------------------------------------------------------------------------------------------------
// Check whether an integration period has ended.
// If so, DMA the results back to RAM and reset everything to zero.
//
inline void check_baseline_integration_period(unsigned int timeperiod) {
    int b0, b1, b2;
    int block;

    if ((timeperiod % cb.ffts_to_integrate) == 0) {

        // determine where to write baselines
        b0 = g_rank;                        //   0..5
        b1 = (g_rank + NUM_SPE_THREADS);    //  6..11
        b2 = (g_rank + 2*NUM_SPE_THREADS);  // 12..17
        block = (timeperiod/cb.ffts_to_integrate) - 1;

        // send the autocorrelation
        store_data(cb.autocorrelation_out, (float*)autocorr, block, PPU_CORRELBUF_BYTES, TAG_PPUOUT);

        // send the baseline data
        store_data(cb.baselines_out[b0], (float*)baselines[0], block, PPU_CORRELBUF_BYTES, TAG_PPUOUT);
        store_data(cb.baselines_out[b1], (float*)baselines[1], block, PPU_CORRELBUF_BYTES, TAG_PPUOUT);

        // there'll be 6*3-15=3 leftover SPEs that shouldn't send, but
        // TODO improve/fix because some non-redundant baselines may be left out
        if (b2 < NUM_BASELINES) {
            store_data(cb.baselines_out[b2], (float*)baselines[2], block, PPU_CORRELBUF_BYTES, TAG_PPUOUT);
        }

        wait_for_tag(TAG_PPUOUT);

        // reset to zero
        zero_baselines_and_autocorr();
    }
    return;
}


// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
