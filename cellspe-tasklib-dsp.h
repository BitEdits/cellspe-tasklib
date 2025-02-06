
/*************************************************************************** 
 *  Copyright (C) 2007 by Jan Wagner                                       *
 *                                                                         *
 *  Digital Signal Processing Primitives Library                           *
 *  for IBM Cell Broadband Engine                                          *
 *                                                                         *
 *  Cell SPU processor Version                                             *
 *  Using IBM Cell SDK Version 2.0                                         *
 *                                                                         *
 *  License: GNU LGPL                                                      *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
 ***************************************************************************/

#ifndef __CELLSPE_TASKLIB_DSP_H
#define __CELLSPE_TASKLIB_DSP_H

#include "cellspe-tasklib.h"

//---------------------------------------------------------------------------------
// Common Declares
//---------------------------------------------------------------------------------

enum CellDspCommandEnum {
  CELLDSP_TASK_NOP = 0,         // nop operation
  CELLDSP_TASK_MAC,             // multiply accumulate : in1[], in2[], out[] = out[] + (in1[] * in2[])
  CELLDSP_TASK_MULC,            // complex multiply    : in1[], in2[], out[] = in1[] * in2[]
  CELLDSP_TASK_ADDC,            // complex add         : in1[], in2[], out[] = in1[] + in2[]
  CELLDSP_TASK_IOTEST,          // DMA in and out check
  CELLDSP_TASK_IOBENCH,         // DMA in and out larger benchmark
  CELLDSP_TASK_FFT_INIT,        // prepare continuable FFT (keeps SPU reserved)
  CELLDSP_TASK_FFT_FREE,        // release continuable FFT
  CELLDSP_TASK_FFT_EXEC,        // execute one FFT of continuables
  CELLDSP_TASK_FFT_EXEC_DYN,    // FFT without init required, runs the init if fft size changes
  CELLDSP_TASK_SINCOS,          // calc interleaved cosine+sine table from input vector
  CELLDSP_TASK_SINCOS_NORM,     // -"- from normalized input phases vector (1.0 ~= 2pi)
  CELLDSP_TASK_SINCOS_GENERTR,  // sine&cosine oscillator
  CELLDSP_TASK_UNPACK,          // split 2-bit input data into (N x NumChannels) float sample table
  CELLDSP_TASK_RUNALL,          // run all of the calculation funcs but use dummy data
  CELLDSP_TASK_DMASIMULATION,   // test DMA with DMA-list, needs one 16kB input and eight 8kB output arrays on PPU
  CELLDSP_TASK_DIFXPROCESS,     // -- this is the DiFX Mode::process() to write, step by step, ...
  CELLDSP_TASK_SIGTERM,         // terminate (always keep as last entry here)
};

//---------------------------------------------------------------------------------
// Global SPU funcs
//---------------------------------------------------------------------------------
#ifdef __SPU__

int do_in_out_bench();
int do_complex_mac();
int do_runall_test();
int do_io_simulation();
int cornerturn(vector unsigned short* indata, unsigned int ** separatedvectors, int count);
int init_FFT();
int free_FFT();
int exec_FFT();
int exec_FFT_dyn();
int exec_sincos();
int exec_sincosnormalized();
int exec_sincos_generator();
int exec_unpack();
int exec_unpack_improved();
int exec_vecrepack_versio2();
int exec_difxprocess();
int do_nop();
int do_complex_mul();
int do_complex_add();
int do_in_out_test();
#endif

//---------------------------------------------------------------------------------
// Global PPU funcs
//---------------------------------------------------------------------------------
#ifndef __SPU__
int speNOP();
int speRunAllCheck();
int speIOCheck();
int speIOSimulation();
int speIOBench(void* buffer, unsigned int bufbytes, int spes);
int speIOBenchStreamed(void* buffer, unsigned int bufbytes, int spes, int stream_length);
int speAddProduct_32fc(fc32 * src, fc32 * src2, fc32 * accumulator, int length);
int speFFT_init(int fft_len, int* speID);
int speFFT_exec(void* in, void* out, int fft_len, int* speID);
int speFFT_destroy(int* speID);
int speFFT_singleexec(void* in, void* out, int fft_len);
int speDIFX_process(void* sincosrotdone, void* complexfracmult, void* fftoutJ, void* conjfftoutJ, void* autocorr,
                    int nfft, int twicenumchannels, int useUpperOrLower);
int speSinCosGenerator(float startphase, float increment, fc32* cos_sin_out, int numiter);
int speSinCosAccurate(float * argument, cf32 * cos_sin_out, int argumentcount);
int speSinCosAccurateNorm(float * normalizedargument, cf32 * cos_sin_out, int argumentcount);
#endif

#endif // __CELLSPE_TASKLIB_DSP_H
