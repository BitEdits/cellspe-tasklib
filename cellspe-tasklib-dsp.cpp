
/*************************************************************************** 
 *  Copyright (C) 2007 by Jan Wagner                                       *
 *                                                                         *
 *  Digital Signal Processing Primitives Library                           *
 *  for IBM Cell Broadband Engine                                          *
 *                                                                         *
 *  Cell SPU Processor Version                                             *
 *  Contains: DSP Code Implementations                                     *
 *                                                                         *
 *  License: GNU LGPL                                                      *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
 ***************************************************************************/

#include "cellspe-tasklib.h"
#include "cellspe-tasklib-dsp.h"

//---------------------------------------------------------------------------------
// Global PPU funcs instantiated in this .cpp
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
int speSinCosGenerator(float startphase, float increment, void * cos_sin_out);
int speSinCosAccurate(void * argument, void * cos_sin_out);
int speSinCosAccurateNorm(void * normalizedargument, void * cos_sin_out);
int speDIFX_process(void* sincosrotdone, void* complexfracmult, void* fftoutJ, void* conjfftoutJ, void* autocorr,
                    int nfft, int twicenumchannels, int useUpperOrLower);
#else
//---------------------------------------------------------------------------------
// Global SPU funcs instantiated in this .cpp
//---------------------------------------------------------------------------------
int do_in_out_test();
int do_in_out_bench();
int do_complex_mac();
int do_runall_test();
int do_io_simulation();
int do_nop();
int cornerturn(vector unsigned short* indata, unsigned int ** separatedvectors, int count);
int init_FFT();
int free_FFT();
int exec_FFT();
int exec_sincos();
int exec_sincosnormalized();
int exec_sincos_generator();
int exec_unpack();
int exec_unpack_improved();
int exec_vecrepack_versio2();
#endif



//===========================================================================
//=== P P U   C O D E                                                     ===
//===========================================================================
#ifndef __SPU__

using namespace std;

//---------------------------------------------------------------------------------
// Local helpers
//---------------------------------------------------------------------------------
inline void PRINT_SPU_TASKSTATS(CellDSPTask_t * finishedtask) 
{
   unsigned int ticks = finishedtask->tick_count;
   cout << "SPU decremeter: " << ticks << " ticks ("
        << SPU_TICKS_TO_USEC(ticks) << " usec, "
        << "or " << SPU_TICKS_TO_CYCLES(ticks) << " cycles)"
        << endl << flush;
   return;
}

inline void PRINT_VECF(vector float *v) { 
   float * pv = (float*)v; 
   printf("vec = %f %f %f %f\n", pv[0], pv[1], pv[2], pv[3]); 
}

//---------------------------------------------------------------------------------
// PPU IMPLEMENTATIONS OF THE DSP FUNCTION CALLS
//---------------------------------------------------------------------------------

int speNOP() 
{
   // initialize task context
   CellDSPTaskClass nop(CELLDSP_TASK_NOP);
   // for additional DMA check:
   static unsigned char tmp_buf[16384] _QUAD_ALIGN;
   nop.addInVec((void*)&tmp_buf[0], 16384);
   // execute with blocking wait
   return celldsp_executeTask(nop.getTask());
}

int speRunAllCheck()
{
   CellDSPTaskClass all(CELLDSP_TASK_RUNALL);
   return celldsp_executeTask(all.getTask());
}

int speIOCheck() 
{
   int status;
   static float output[8] _QUAD_ALIGN;
   // initialize task context
   CellDSPTaskClass iotest(CELLDSP_TASK_IOTEST);
   iotest.addOutVec((void*)&output[0], 8 * sizeof(float));
   iotest >> cout; cout << endl << flush; // debug dump of context contents
   // execute with blocking wait
   status = celldsp_executeTask(iotest.getTask());
   cout << "speIOCheck wrote back : vect={" << output[0] << "," << output[1] 
        << ","<< output[2] << "," << output[3] << "}   expected would be {0.1,2.3,4.5,6.7}" << endl;
   return status;
}

int speIOSimulation() 
{
   unsigned char inbuf[16384] _QUAD_ALIGN;   // DMA-in to SPE
   unsigned char outbuf[8*8192] _QUAD_ALIGN; // DMA-out to PPU
   // initialize task context
   CellDSPTaskClass iosimu(CELLDSP_TASK_DMASIMULATION);
   iosimu.addInVec((void*)&inbuf[0], 16384);
   for(int k=0; k<8; k++) {
      iosimu.addOutVec((void*)&outbuf[k*8192], 8192);
   }
   iosimu >> cout; cout << endl << flush; // debug dump of context contents
   // execute with blocking wait
   return celldsp_executeTask(iosimu.getTask());
}

int speIOBench(void* buffer, unsigned int bufbytes, int spes)
{
   int status = 0;
   int * myspes = new int[spes];

   // initialize task context
   CellDSPTask_t ** IOBench = new CellDSPTask_t*[spes];
   for (int i=0; i<spes; i++) {
      IOBench[i] = new CellDSPTask_t;
      IOBench[i]->command                 = CELLDSP_TASK_IOBENCH;
      IOBench[i]->num_input_vectors       = 0;
      IOBench[i]->num_output_vectors      = 1;
      IOBench[i]->outputvector_ptrs[0].p  = buffer;
      IOBench[i]->outputvector_lengths[0] = bufbytes;
      IOBench[i]->next_task.p = NULL;
   }

   // reserve free SPE(s) and start the tasks on them
   for (int i=0; i<spes; i++) {
     myspes[i] = celldsp_getFreeSPE();
     if (myspes[i] == -1) {
        for (int j=0; j<i; j++) {
           celldsp_setFreeSPE(myspes[j]);
        }
        cout << "speIOBench failed: could not get enough free SPEs";
        return 1;
     }
     celldsp_startTask(myspes[i], IOBench[i]);
   }
   // wait for all to complete
   for (int i=0; i<spes; i++) {
     status = celldsp_waitTask(myspes[i]);
   }
   // release
   for (int i=0; i<spes; i++) {
     celldsp_setFreeSPE(myspes[i]);
     delete IOBench[i];
   }
   delete[] myspes;
   delete[] IOBench;
   return status;
}

int speIOBenchStreamed(void* buffer, unsigned int bufbytes, int spes, int stream_length)
{
   int status = 0;
   // initialize task context
   cout << " speIOBenchStreamed IOBench[] prepare ... " << flush;
   CellDSPTask_t * IOBench[stream_length] _QUAD_ALIGN;
   for (int i=0; i<stream_length; i++) {
      IOBench[i] = new CellDSPTask_t;
      IOBench[i]->command = CELLDSP_TASK_IOBENCH;
      IOBench[i]->num_input_vectors       = 0;
      IOBench[i]->num_output_vectors      = 1;
      IOBench[i]->outputvector_ptrs[0].p  = buffer;
      IOBench[i]->outputvector_lengths[0] = bufbytes;
      IOBench[i]->next_task.p = NULL; 
      if (i > 0) {
         IOBench[i-1]->next_task.p = (void*)IOBench[i];
      }
   }
   cout << " done." << endl << flush;
   // start the tasks on different SPE's, use 
   // different task descriptor offset for each SPE
   int * myspes = new int[spes];
   for (int i=0; i<spes; i++) {
     myspes[i] = celldsp_getFreeSPE();
     if (myspes[i] == -1) {
        for (int j=0; j<i; j++) {
           celldsp_setFreeSPE(myspes[j]);
        }
        cout << "speIOBenchStreamed failed: could not get enough free SPEs";
        return 1;
     }
     celldsp_startTask(myspes[i], IOBench[i % stream_length]);
   }
   cout << " speIOBenchStreamed wait for tasks to finish ... " << flush;
   for (int i=0; i<spes; i++) {
     status = celldsp_waitTask(myspes[i]);
   }
   cout << " done. " << endl << flush;
   for (int i=0; i<stream_length; i++) {
      delete IOBench[i];
   }
   // release
   for (int i=0; i<spes; i++) {
     celldsp_setFreeSPE(myspes[i]);
   }
   delete myspes;
   return status;
}

int speAddProduct_32fc(fc32 * src, fc32 * src2, fc32 * accumulator, int length) 
{
   // initialize task context
   CellDSPTaskClass mac(CELLDSP_TASK_MAC);
   mac.addInVec((void*)src, length * sizeof(fc32));
   mac.addInVec((void*)src2, length * sizeof(fc32));
   mac.addOutVec((void*)accumulator, length * sizeof(fc32));
   // execute with blocking wait
   return celldsp_executeTask(mac.getTask());
}

// Get a free SPE and write the ID into *speID, then initialize FFT on the SPE
int speFFT_init(int fft_len, int* speID)
{
   // initialize task context
   CellDSPTaskClass fftinit(CELLDSP_TASK_FFT_INIT);
   fftinit.setFFTparams(fft_len);
   // initialize the FFT
   int status = celldsp_startContinuableTask(fftinit.getTask(), speID);
   //PRINT_SPU_TASKSTATS(fftinit.getTask());
   return status;
}

// Execute an FFT on an earlier initialized SPE with ID in *speID
int speFFT_exec(void* in, void* out, int fft_len, int* speID)
{
   // initialize task context
   CellDSPTaskClass fft(CELLDSP_TASK_FFT_EXEC);
   fft.addInVec(in, fft_len * sizeof(fc32));
   fft.addOutVec(out, fft_len * sizeof(fc32));
   fft.setFFTparams(fft_len);
   // perform the FFT
   return celldsp_continueTask(fft.getTask(), speID);
}

// End an FFT task block on an earlier initialized SPE with ID in *speID
int speFFT_destroy(int* speID)
{
   // initialize task context
   CellDSPTaskClass fftfree(CELLDSP_TASK_FFT_FREE);
   int status = celldsp_continueTask(fftfree.getTask(), speID);
   // terminate sequence
   celldsp_stopContinuableTask(speID);
   return status;
}

// Execute an FFT on some of the SPEs, doesn't have to be continuable task
int speFFT_singleexec(void* in, void* out, int fft_len)
{
   // initialize task context
   CellDSPTaskClass fft(CELLDSP_TASK_FFT_EXEC_DYN);
   fft.addInVec(in, fft_len * sizeof(fc32));
   fft.addOutVec(out, fft_len * sizeof(fc32));
   fft.setFFTparams(fft_len);
   // perform the FFT
   return celldsp_executeTask(fft.getTask());
}

// Execute sine cosine generator (oscillator), amplitude is not normalized!!
int speSinCosGenerator(float startphase, float increment, fc32* cos_sin_out, int numiter)
{
   // init
   CellDSPTaskClass scgen(CELLDSP_TASK_SINCOS_GENERTR);
   scgen.addOutVec((void*)cos_sin_out, numiter*sizeof(cf32));
   CellDSPParams_Phase_tt * aux = &(scgen.getSubcontext()->Phase);
   aux->argcount = numiter;
   aux->phase_start = startphase;
   aux->phase_delta = increment;
   // perform
   return celldsp_executeTask(scgen.getTask());
}

// Execute sine cosine calculation from phase argument
int speSinCosAccurate(float * argument, cf32 * cos_sin_out, int argumentcount)
{
   // init
   CellDSPTaskClass sc(CELLDSP_TASK_SINCOS);
   sc.addInVec((void*)argument, argumentcount*sizeof(float));
   sc.addOutVec((void*)cos_sin_out, argumentcount*sizeof(cf32));
   CellDSPParams_Phase_tt * aux = &(sc.getSubcontext()->Phase);
   aux->argcount = argumentcount;
   // perform
   return celldsp_executeTask(sc.getTask());
}

// Execute sine cosine calculation from normalized phase argument
int speSinCosAccurateNorm(float * normalizedargument, cf32 * cos_sin_out, int argumentcount)
{
   // init
   CellDSPTaskClass sc(CELLDSP_TASK_SINCOS_NORM);
   sc.addInVec((void*)normalizedargument, argumentcount*sizeof(float));
   sc.addOutVec((void*)cos_sin_out, argumentcount*sizeof(cf32));
   CellDSPParams_Phase_tt * aux = &(sc.getSubcontext()->Phase);
   aux->argcount = argumentcount;
   // perform
   return celldsp_executeTask(sc.getTask());

}


// Execute some things moved to SPE side from DiFX Mode::process()
int speDIFX_process(void* sincosrotdone, void* complexfracmult, void* fftoutJ, void* conjfftoutJ, void* autocorr,
                    int nfft, int twicenumchannels, int useUpperOrLower)
{
   int dual_array_len = twicenumchannels*sizeof(cf32);
   int normal_array_len = (twicenumchannels/2 + 1 + 3)*sizeof(cf32);
//   printf("speDIFX_process(in0=%p, in1=%p, in2=%p, in3=%p, in4=%p, nfft=%d, 2nch=%d, usblsb=%d)\n", sincosrotdone, complexfracmult, fftoutJ, conjfftoutJ, autocorr, 
//          nfft, twicenumchannels, useUpperOrLower);
   // initialize task context
   CellDSPTaskClass difx(CELLDSP_TASK_DIFXPROCESS);
   difx.addInVec(sincosrotdone, dual_array_len);
   difx.addInVec(complexfracmult, normal_array_len);
   difx.addOutVec(fftoutJ, normal_array_len);
   difx.addOutVec(conjfftoutJ, normal_array_len);
   difx.addOutVec(autocorr, normal_array_len);

   // pass some parameters in the subcontext
   CellDSPSubcontext_t * subctx = difx.getSubcontext();
   CellDSPParams_Difx_t * sctx = &subctx->Difx;
   sctx->nfft = nfft;
   sctx->useUpperOrLower = useUpperOrLower;
   sctx->twicenumchannels = twicenumchannels;
   sctx->numchannels = twicenumchannels/2;

   // perform, blocking
   return celldsp_executeTask(difx.getTask());
}

/**************************************************************************/




//===========================================================================
//=== S P U   C O D E                                                     ===
//===========================================================================
#else

// function specific globals
unsigned int fft_N;         // FFT length N
unsigned int fft_N_log2;    // log2 of -"-
fc32* fft_W;                // FFT twiddle factors
unsigned char noptask_dummy[16384] _QUAD_ALIGN;

//---------------------------------------------------------------------------------
// SPU IMPLEMENTATIONS OF TESTED AND WORKING FUNCS
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// do_nop() : just for checking
//---------------------------------------------------------------------------------
int do_nop() {
    return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// do_complex_mul_inplace() : complex float product, in2out = in*in2out
//   current data layout in vector: [Re0 Im0 Re1 Im1]
//   # of complex values must be a multiple of 4
//---------------------------------------------------------------------------------
inline int do_complex_mul_inplace(vector float* in, vector float* in2out, int numcomplexes) {
    vector float A1, A2, B1, B2, I1, I2, Q1, Q2;  
    vector float *pD1; 
    vector float *pD2;
    vector float v_zero = (vector float){0,0,0,0};

    for (int i=0; i<numcomplexes/(2*2); i++) { // 2 cplx per vector, process 2 vectors at once
      A1 = *in++; A2 = *in++;
      pD1 = in2out; B1 = *in2out++; 
      pD2 = in2out; B2 = *in2out++;

			/* in-phase (real), quadrature (imag), temp, and output vectors*/
      vector unsigned char I_Perm_Vector = (vector unsigned char){0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27};
      vector unsigned char Q_Perm_Vector = (vector unsigned char){4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31};
      vector unsigned char vcvmrgh = (vector unsigned char) {0,1,2,3,16,17,18,19,4,5,6,7,20,21,22,23};
      vector unsigned char vcvmrgl = (vector unsigned char) {8,9,10,11,24,25,26,27,12,13,14,15,28,29,30,31};
  
      /* input vectors are in interleaved form in A1,A2 and B1,B2 with each input vector
         representing 2 complex numbers and thus this loop would repeat for N/4 iterations
       */    
      I1 = spu_shuffle(A1, A2, I_Perm_Vector); /* pulls out 1st and 3rd 4-byte element from vectors A1 and A2 */
      I2 = spu_shuffle(B1, B2, I_Perm_Vector); /* pulls out 1st and 3rd 4-byte element from vectors B1 and B2 */
      Q1 = spu_shuffle(A1, A2, Q_Perm_Vector); /* pulls out 2nd and 4th 4-byte element from vectors A1 and A2 */
      Q2 = spu_shuffle(B1, B2, Q_Perm_Vector); /* pulls out 3rd and 4th 4-byte element from vectors B1 and B2 */
      A1 = spu_nmsub(Q1, Q2, v_zero);          /* calculates -(bd + 0) for all four elements */
      A2 = spu_madd(Q1, I2, v_zero);           /* calculates (bc + 0) for all four elements */
      Q1 = spu_madd(I1, Q2, A2);               /* calculates ad + bc for all four elements */
      I1 = spu_madd(I1, I2, A1);               /* calculates ac - bd for all four elements */ 
      *pD1 = spu_shuffle(I1, Q1, vcvmrgh);       /* spreads the results back into interleaved format */
      *pD2 = spu_shuffle(I1, Q1, vcvmrgl);       /* spreads the results back into interleaved format */
    }
    return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// do_complex_autocorrelation() : complex autocorrelation with accumulation
//                                in autocorrelation the imag part ends up as zero
//---------------------------------------------------------------------------------
inline int do_complex_autocorrelation(vector float* in, vector float* accu, int numcomplexes) {
   vector unsigned char Im_into_ReLocation = (vector unsigned char) { 0x04,0x05,0x06,0x07,  0x80,0x80,0x80,0x80, 12,13,14,15, 0x80,0x80,0x80,0x80 };
   vector unsigned char Re_into_ReLocation = (vector unsigned char) { 0,1,2,3,  0x80,0x80,0x80,0x80, 8,9,10,11, 0x80,0x80,0x80,0x80 };
   for (int i=0; i<numcomplexes/2; i++) { // 2 cplx per vector
   /* version 1: 
      vector float tmpI = spu_shuffle(*in, *in, Im_into_ReLocation);
      vector float tmpOut = spu_madd(tmpI, tmpI, *accu);
      *accu = spu_madd(*in, *in, tmpOut);
      in++; accu++;
   */  
   /* version 2: better pipeline(?) */
      vector float tmpR = spu_shuffle(*in, *in, Re_into_ReLocation);
      *accu = spu_madd(tmpR, tmpR, *accu);
      vector float tmpI = spu_shuffle(*in, *in, Im_into_ReLocation);
      *accu = spu_madd(tmpI, tmpI, *accu);
      accu++; in++;
   }
   return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// do_complex_mac() : complex multiply-and-accumulate
//---------------------------------------------------------------------------------
int do_complex_mac() {

    vector float * in1 = (vector float*) args[0];
    vector float * in2 = (vector float*) args[1];
    vector float * acc = (vector float*) args[2];
    int elements = myTask.inputvector_lengths[0]/(2*sizeof(float)); // # of complex nums
  
    // what this func should do, is:
    // for (int i=0; i<elements/2; i++) {  
    //     re_tmp[0|2] = in1[0|2]*in2[0|2] - in1[1|3]*in2[1|3]; // re*re - im*im
    //     im_tmp[1|3] = in1[0|2]*in2[1|3] + in1[1|3]*in2[0|2]; // re*im + im*re
    //     acc[0|2] = acc[0|2] + re_tmp[0|2]; 
    //     acc[1|3] = acc[1|3] + re_tmp[0|2];
    // }
/* Jouko:
		 multiplyrealimmediate b,{1,-1,1,-1}
		 x1 = mask imaginary from a
		 x2 = rotate quadword a left by 4 bytes, mask imaginary
		 x3 = rotate quadword b left by 4 bytes, mask imaginary
		 x4 = rotate quadword a right by 4 bytes, mask real
		 x5 = rotate quadword b right by 4 bytes, mask real
		 multiplyandadd x1,b,result
		 multiplyandsubstract x2,x3,result
		 multiplyandadd x4,b,result
		 multiplyandadd x5,a,result
*/	

    /* credits go to IBM(?) for this one, copied it from one of their Powerpoint presentations */

    vector float A1, A2, B1, B2, I1, I2, Q1, Q2, D1, D2;  
    /* in-phase (real), quadrature (imag), temp, and output vectors*/
    vector float v_zero = (vector float){0,0,0,0};

    vector unsigned char I_Perm_Vector = (vector unsigned char){0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27};
    vector unsigned char Q_Perm_Vector = (vector unsigned char){4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31};
    vector unsigned char vcvmrgh = (vector unsigned char) {0,1,2,3,16,17,18,19,4,5,6,7,20,21,22,23};
    vector unsigned char vcvmrgl = (vector unsigned char) {8,9,10,11,24,25,26,27,12,13,14,15,28,29,30,31};

    for (int i=0; i<elements/4; i++) { // 2 complex nums in each 128-bit vector, process 4 nums each iteration
       A1 = *(in1++); A2 = *(in1++);
       B1 = *(in2++); B2 = *(in2++);
       /* input vectors are in interleaved form in A1,A2 and B1,B2 with each input vector
        * representing 2 complex numbers and thus this loop would repeat for N/4 iterations
        */    
       I1 = spu_shuffle(A1, A2, I_Perm_Vector); /* pulls out 1st and 3rd 4-byte element from vectors A1 and A2 */
       I2 = spu_shuffle(B1, B2, I_Perm_Vector); /* pulls out 1st and 3rd 4-byte element from vectors B1 and B2 */
       Q1 = spu_shuffle(A1, A2, Q_Perm_Vector); /* pulls out 2nd and 4th 4-byte element from vectors A1 and A2 */
       Q2 = spu_shuffle(B1, B2, Q_Perm_Vector); /* pulls out 3rd and 4th 4-byte element from vectors B1 and B2 */
       A1 = spu_nmsub(Q1, Q2, v_zero);          /* calculates  - (bd-0) for all four elements */
       A2 = spu_madd(Q1, I2, v_zero);           /* calculates (bc + 0) for all four elements */
       Q1 = spu_madd(I1, Q2, A2);               /* calculates ad + bc for all four elements */
       I1 = spu_madd(I1, I2, A1);               /* calculates ac - bd for all four elements */ 
       D1 = spu_shuffle(I1, Q1, vcvmrgh);       /* spreads the results back into interleaved format */
       D2 = spu_shuffle(I1, Q1, vcvmrgl);       /* spreads the results back into interleaved format */
       *acc = spu_add(*acc, D1); acc++;
       *acc = spu_add(*acc, D2); acc++;
    }

    write_results();
    return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// init_FFT() : initialize the FFT coefficients table W, and calculate log2(N)
//---------------------------------------------------------------------------------
int init_FFT() {
   unsigned int nfft;

   if (fft_W != NULL) {
      printf("SPE init_FFT(): fft_W!=NULL, will free up first\n");
      free_FFT();
   }
   nfft = myTask.subcontext.FFT.nfft;
   fft_N = nfft;
   fft_W = (fc32*)malloc_align((nfft/4 + 1) * 8, MALLOC_QUADALIGN); // 2 floats per N
   #if VERBOSE_MALLOC
   if (fft_W==NULL) {
     printf("init_FFT: fft_W %d byte malloc failed\n", (nfft/4 + 1) * 8);
     return CELLDSP_TASK_GENERICERROR;
   }
   #endif

   fft_W->re = 1.0;
   fft_W->im = 0.0;
   for (unsigned int i=1; i<nfft/4; i++) {
      (fft_W+i)->re = cos(i * 2*M_PI/nfft);
      (fft_W + nfft/4 - i)->im = -(fft_W + i)->re;
   }
   while (nfft > 1) {
     fft_N_log2++; 
     nfft = nfft>>1;
   }
   return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// free_FFT() : free up the coefficients table W allocation 
//---------------------------------------------------------------------------------
int free_FFT() {
   if (fft_W != NULL) {
      free_align((void*)fft_W);
      fft_W = NULL;
      fft_N = 0;
      fft_N_log2 = 0;
   }
   return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// exec_FFT() : perform the FFT
//---------------------------------------------------------------------------------
int exec_FFT()
{
   vector float* in  = (vector float*) args[0];
   vector float* out = (vector float*) args[1];

   // perform FFT
   _fft_1d_r2(out, in, (vector float*)fft_W, fft_N_log2);
   
   // scale by N (depends on how you want it... on FFT or IFFT)
//   vector float scale;
//   scale = spu_splats(1.0f/fft_N); 
//   for (unsigned int i=0; i<fft_N/2; i++) {
//      // needs only N/2 mul ops since there are two complex nums in each 128-bit vector
//      *(out+i) = spu_mul(*(out+i), scale);
//   }

   write_results();

   return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// exec_FFT_dyn() : perform the FFT, with dynamic W twiddlefactor re-calc
//---------------------------------------------------------------------------------
int exec_FFT_dyn()
{
   vector float* in  = (vector float*) args[0];
   vector float* out = (vector float*) args[1];
   if(__builtin_expect((myTask.subcontext.FFT.nfft != fft_N), 0)) {
      init_FFT();
   }

   // perform FFT
   _fft_1d_r2(out, in, (vector float*)fft_W, fft_N_log2);
   write_results();

   return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// exec_difxprocess() : perform some parts of the DiFX Mode::process()
//                      unfortunately this is a bit low-level, to take advantage
//                      of DMA & arithmetics interleaving...
//---------------------------------------------------------------------------------
int exec_difxprocess()
{
   int bytelen_dual = myTask.inputvector_lengths[0];
   int bytelen_single = myTask.outputvector_lengths[0]; 

   //  CellDSPTaskClass difx(CELLDSP_TASK_DIFXPROCESS);
   //  difx.addInVec(sincosrotdone, dual_array_len);  // full FFT on this
   //  difx.addInVec(complexfracmult, normal_array_len);
   //  difx.addOutVec(fftoutJ, normal_array_len);
   //  difx.addOutVec(conjfftoutJ, normal_array_len);
   //  difx.addOutVec(autocorr, normal_array_len);
   //  sctx->nfft = nfft;
   //  sctx->useUpperOrLower = useUpperOrLower
   //
   // sincos:
   // simdmath   -  sincosf4: Sine and Cosine of Float
   //      (void) sincosf4 (vector float x, vector float *sx, vector float *cx);

   // helpers for mapping of args[x] to something more descriptive:
   vector float* complexfracmult;  
   vector float* autocorr;
   vector float* fftd;

   // helpers for mapping output outputvector_ptrs[x] to something more descriptive:
   Cell_addr64*  ppuptr_fftout = (Cell_addr64*)&myTask.outputvector_ptrs[0]; // discarding 'volatile' with cast
   Cell_addr64*  ppuptr_conjfftout = (Cell_addr64*)&myTask.outputvector_ptrs[1];
   Cell_addr64*  ppuptr_autocorr = (Cell_addr64*)&myTask.outputvector_ptrs[2];
   // same for inputvector_ptrs[x]:
   Cell_addr64*  ppuptr_sincosrotatedoutput = (Cell_addr64*)&myTask.inputvector_ptrs[0];
   Cell_addr64*  ppuptr_complexfracmult = (Cell_addr64*)&myTask.inputvector_ptrs[1];

   /* start fetch of FFT input data */
   args[0] = malloc_align(bytelen_dual, MALLOC_QUADALIGN);
   //printf("exec_difxprocess - FFT(nfft=%d) data fetch from %p to %p len=%d\n",
   //     myTask.subcontext.Difx.nfft, 
   //    (void*)myTask.inputvector_ptrs[0].ui[0], args[0], bytelen_dual);
   fftd = (vector float*) args[0];
   spu_dma_nowait(args[0], ppuptr_sincosrotatedoutput, bytelen_dual, DMA_TAG_READABLES, MFC_GET_CMD);

   /* start fetch of fractional correction multipliers */
   args[1] = malloc_align(bytelen_single, MALLOC_QUADALIGN);   
   complexfracmult = (vector float*)args[1];
   spu_dma_nowait(args[1], ppuptr_complexfracmult, bytelen_single, DMA_TAG_READABLES_DB, MFC_GET_CMD);

   /* init FFT W matrix if necessary, wait for FFT input data */
   if(__builtin_expect((myTask.subcontext.Difx.nfft != fft_N), 0)) {
      init_FFT();
   }
   mfc_write_tag_mask(DMA_TAG_MASK_READABLES);
   mfc_read_tag_status_all();

   /* do the in-place FFT */
   _fft_1d_r2(fftd, fftd, (vector float*)fft_W, fft_N_log2);
   // select which sideband to use
   void * fftband;
   if(myTask.subcontext.Difx.useUpperOrLower) {
      fftband = (unsigned char*)args[0] + myTask.outputvector_lengths[0]; ((myTask.subcontext.Difx.twicenumchannels<<1) * sizeof(cf32));
   } else {
      fftband = args[0];
   }

   /* start fetch of previous autocorrelation data for MAC */
   args[2] = malloc_align(bytelen_single, MALLOC_QUADALIGN);
   autocorr = (vector float*)args[2];
   spu_dma_nowait(args[2], ppuptr_autocorr, bytelen_single, DMA_TAG_READABLES, MFC_GET_CMD);

   // wait for the complexfracmult data and do the complex multiply
   mfc_write_tag_mask(DMA_TAG_MASK_READABLES_DB);
   mfc_read_tag_status_all();
   do_complex_mul_inplace(complexfracmult, (vector float*)fftband, myTask.subcontext.Difx.numchannels);

   /* do complex conjugation, MAC accumulation into autocorrelation array */
   // todo Conjugation
   // wait for autocorrelation data to arrive
   mfc_write_tag_mask(DMA_TAG_MASK_READABLES);
   mfc_read_tag_status_all();
   do_complex_autocorrelation((vector float*)fftband, autocorr, myTask.subcontext.Difx.numchannels);
   spu_dma_nowait(autocorr, ppuptr_autocorr, myTask.outputvector_lengths[2], DMA_TAG_WRITEABLES, MFC_PUTF_CMD);

   /* write back the result(s) */
   // printf("exec_difxprocess - FFT(nfft=%d)  data return from %p to %p, len=%d\n", fft_N, fftoutJ, (void*)myTask.outputvector_ptrs[0].ui[0], myTask.inputvector_lengths[0]);
   spu_dma_nowait(fftband, ppuptr_fftout, myTask.outputvector_lengths[0], DMA_TAG_WRITEABLES, MFC_PUTF_CMD);
   mfc_write_tag_mask(DMA_TAG_MASK_WRITEABLES);
   mfc_read_tag_status_all();

   free_align(args[0]);
   free_align(args[1]);
   free_align(args[2]);
   return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// exec_sincos() : calculates sine and cosine of argument vector, and interleaves
//                 the results to [cos0 sin0 cos1 sin1]... formatted output vectors
//---------------------------------------------------------------------------------
int exec_sincos() {
  //sc.addInVec((void*)argument, argumentcount*sizeof(float));
  //sc.addOutVec((void*)cos_sin_out, argumentcount*sizeof(cf32));
  vector float * phase = (vector float*)args[0];
  vector float * cossin = (vector float*)args[1];
  vector float sx, cx;
  for(int i=0; i<myTask.subcontext.Phase.argcount / 4; i++) {
     // sincos4f(*phase++, &sx, &cx); // SPU simdmath SDK2.0 does not have sincos4f()!
     sx = sin14_v(*phase); // sinf4()
     cx = cos14_v(*phase++); // cosf4()
     *cossin++ = spu_shuffle(sx, cx, (vector unsigned char){ 0x10,0x11,0x12,0x13, 0x00,0x01,0x02,0x03, 0x14,0x15,0x16,0x17, 0x04,0x05,0x06,0x07 });
     *cossin++ = spu_shuffle(sx, cx, (vector unsigned char){ 0x18,0x19,0x1A,0x1B, 0x08,0x09,0x0A,0x0B, 0x1C,0x1D,0x1E,0x1F, 0x0C,0x0D,0x0E,0x0F });
  }
  write_results();
  return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// exec_sincos_generator() : calculates N points of a sine and cosine waveform
//                           given the start phase and phase increment
//---------------------------------------------------------------------------------
int exec_sincos_generator()
{
   float angle_start = myTask.subcontext.Phase.phase_start;
   float angle_delta = myTask.subcontext.Phase.phase_delta;
   vector float * cossinN = (vector float*)args[0];
   unsigned int num_pairs = myTask.subcontext.Phase.argcount / 2; // two cf32's per result vector
   unsigned int pair;
   
   //
   // General case:   Method I
   //
   //  [c_N+1] = [ +cos(phi)  -sin(phi) ] * [c_N]
   //  [s_N+1]   [ +sin(phi)  +cos(phi) ]   [s_N]
   //
   // If cos(phi)=1-sin(phi):   Method II (less accurate)
   //
   //  [c_N+1]  = [ cos(phi)        cos(phi) - 1  ]  * [c_N]
   //  [s_N+1]    [ cos(phi) + 1    cos(phi)      ]    [s_N]
   //  <=> c_N+1 = cos(phi)*c_N + (cos(phi)-1)*s_N = cos(phi)*c_N + cos(phi)*s_N - s_N
   //      s_N+1 = (cos(phi)+1)*c_N + cos(phi)*s_N = cos(phi)*c_N + cos(phi)*s_N + c_N
   //

   // seed
   *cossinN   = (vector float){
                  /* N=0 */  float(cosf(angle_start)), float(sinf(angle_start)), 
                  /* N=1 */  float(cosf(angle_start+angle_delta)), float(sinf(angle_start+angle_delta)) 
                 };
/* Method I */
#if 1

   vector float M_row1 = (vector float) { // for calculating next two cosines
            cosf(2*angle_delta), -sinf(2*angle_delta), 
            cosf(2*angle_delta), -sinf(2*angle_delta) };
   vector float M_row2 = (vector float) { // for calculating next two sines
            sinf(2*angle_delta), cosf(2*angle_delta), 
            sinf(2*angle_delta), cosf(2*angle_delta) };
   vector unsigned char concat_odds = (vector unsigned char) { // combine above cos and sin
          0x00,0x01,0x02,0x03, 0x10,0x11,0x12,0x13, 
          0x08,0x09,0x0A,0x0B, 0x18,0x19,0x1A,0x1B };

   vector float prev = *cossinN;
   for (pair = 1; pair < num_pairs; pair++) {
       vector float A = spu_mul(prev, M_row1);
       vector float Ashift = spu_rlqwbyte(A, 4);
       vector float B = spu_mul(prev, M_row2);
       vector float Bshift = spu_rlqwbyte(B, 4);
       prev = spu_shuffle(spu_add(A, Ashift), spu_add(B, Bshift), concat_odds);
       *(++cossinN) = prev;
   }

/* Method II */
#else
   // increment for next two value pairs, sine amplitude correction
   vector float cos2phi   = spu_splats(float(cosf(2*angle_delta)));
   vector float amplitudecorrect = (vector float){ 1.0, tanf(angle_delta), 1.0, tanf(angle_delta) };

   // iterator temps
   vector float sincosN, sincosNm, sincosSum, tmpout, prev;

   // sine cosine vector value flip and sign toggle for sine
   vector unsigned char cossin2sincos = 
           (vector unsigned char) { 4,5,6,7,  0,1,2,3,  12,13,14,15,  8,9,10,11 };
   vector float signflip = reinterpret_cast<vector float>(
           (vector unsigned char) { 0x80,0,0,0 , 0,0,0,0, 0x80,0,0,0, 0,0,0,0 }); // XOR=>'-+','+-','-+','+-'

   prev = *cossinN;
   for (pair = 1; pair < num_pairs; pair++) {
      sincosN = spu_shuffle(prev, prev, cossin2sincos);  // [Cn0 Sn0 Cn1 Sn1] => [Sn0 Cn0 Sn1 Cn1]
      sincosNm = spu_xor(sincosN, signflip);             // => [-Sn0 +Cn0 -Sn1 +Cn1]
      sincosSum = spu_add(prev, sincosN);                // = [Cn0+Sn0 Sn0+Cn0 Cn1+Sn1 Sn1+Cn1]
      prev = spu_madd(cos2phi, sincosSum, sincosNm);     // = [cos(2phi)*(Cn0+Sn0)-Sn0 cos(2phi)*(Cn0+Sn0)+Cn0 cos(2phi)*(Cn1+Sn1)-Sn1 cos(2phi)*(Cn1+Sn1)+Cn1]
      *cossinN = spu_mul(*cossinN, amplitudecorrect);    // amplitude-correct previous pair
      *(++cossinN) = prev; // store new pair
   }
   *cossinN = spu_mul(*cossinN, amplitudecorrect); // ampl-corr final pair
#endif

   write_results();
   return CELLDSP_TASK_OK;
}



//---------------------------------------------------------------------------------
// SPU IMPLEMENTATIONS OF SOME TOTALLY UNTESTED FUNCTIONS & OTHER INCOMPLETES
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// do_complex_mul() : complex product
//---------------------------------------------------------------------------------
int do_complex_mul() {
    // TODO ...
    return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// do_complex_add() : complex addition
//---------------------------------------------------------------------------------
int do_complex_add() {
    // TODO ... : simply spu_add()
    return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// cornerturn() : do corner turn / data transpose
//---------------------------------------------------------------------------------
int cornerturn(vector unsigned short* indata, unsigned int ** separatedvectors, int count) {
    vector unsigned short apuli1;
    vector unsigned short apuli2;
    vector unsigned short apuli3;
    vector unsigned short apuli4;
    vector unsigned short apuli6;

    vector unsigned short mask_unsigned3;
    vector signed short shifttab = {0,14,12,10,8,6,4,2};
    mask_unsigned3 = spu_splats((unsigned short)3);

    for (int i=0; i<count; i++) {
       for (int j=0; j<8; j++) {
          apuli1 = spu_sl(indata[i], 2*j);
          apuli2 = spu_and(apuli1, mask_unsigned3);
          apuli3 = spu_rl(apuli2, shifttab);

          apuli4 = spu_rlqwbyte(apuli3, 8);
          apuli6 = spu_or(apuli4, apuli3);
          apuli4 = spu_rlqwbyte(apuli4, 4);
          apuli6 = spu_or(apuli4, apuli6);
          apuli4 = spu_rlqwbyte(apuli4, 2);
          apuli6 = spu_or(apuli4, apuli6);

          separatedvectors [j][i] = spu_extract(apuli6, 0);
      }
   }
   return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// exec_unpack() : perform packed 2-bit samples into float array unpacking
//---------------------------------------------------------------------------------
int exec_unpack() {
   return CELLDSP_TASK_OK;
}

//---------------------------------------------------------------------------------
// exec_sincosnormalized() : same as exec_sincos(), but arguments expected in
//                           normalized format (1.0 ~= 2pi)
//---------------------------------------------------------------------------------
//#include <cos_sin18_v.h>

  static vector float cos18_A_vDSP[2] = {
     {(vector float){0.00012354,0.00036588,0.00059416,0.00079961}},
     {(vector float){0.00097433,0.00111160,0.00120616,0.00125436}} };
  static vector float cos18_B_vDSP[2] = {
     {(vector float){-0.01933826,-0.01896574,-0.01786437,-0.01607648}},
     {(vector float){-0.01367078,-0.01073973,-0.00739595,-0.00376795}} };
  static vector float cos18_C_vDSP[2] = {
     {(vector float){-0.00000000,-0.03830589,-0.07513972,-0.10908596}},
     {(vector float){-0.13884009,-0.16325867,-0.18140332,-0.19257674}} };
  static vector float cos18_D_vDSP[2] = {
     {(vector float){1.00000000,0.98078525,0.92387950,0.83146960}},
     {(vector float){0.70710677,0.55557024,0.38268343,0.19509032}} };

vector unsigned int glob_negflag = {1,1,1,1};
vector float glob_a = {1,1,1,1}, glob_b = {1,1,1,1}, glob_c = {1,1,1,1}, glob_d = {1,1,1,1};

static __inline vector float _cos_sin18_vrep_MOD(vector float angle) {
  // ***** This part for every four samples
  // ***** For best performance calculate four vectors in one loop
  // ***** Or still better, sin and cos for two samples...

  vector float ctmpl, ctmph;
  vector float result;
  vector float ang2;

  ang2 = spu_mul(angle, angle);
  ctmpl = spu_madd(angle, glob_c, glob_d);
  ctmph = spu_madd(angle, glob_a, glob_b);
  result = spu_madd(ang2, ctmph, ctmpl);
  result = (vector float)spu_or((vector unsigned int)result, glob_negflag);

  // ***** End of sample-specific part
  return result;
}

static __inline vector float _cos_sin18_vDSP(vector float angle)
{
  vector unsigned int idx;
  vector unsigned int quadrant23;
  vector float ctmpl, ctmph;
  vector float ang, result;
  vector float fang;
  vector float ang2;
  vector unsigned int iang;
  vector unsigned int mirror_iang;
  vector float mirror_angf;

  ang = angle;
  ang = (vector float)spu_rlmask(spu_sl((vector unsigned int)ang, 1), -1);

  iang  = spu_convtu(ang, 0);
  fang  = spu_convtf(iang, 0);
  ang   = spu_sub(ang, fang);

  /* Handle mirroring of the function */
  quadrant23 = spu_and(iang, 8);
  quadrant23 = spu_cmpeq(quadrant23, 8);

  mirror_iang = spu_xor(iang, 0x17);
  mirror_angf = spu_sub(VEC_SPLAT_F32(1.0), ang);

  iang = spu_sel(iang, mirror_iang, quadrant23);
  ang  = spu_sel(ang, mirror_angf, quadrant23);

  /* Determine correct resultant sign */
  glob_negflag  = spu_and(iang, 0x10);
  glob_negflag  = spu_sl(glob_negflag, 27);

  idx      = spu_and(iang, 7);
  idx      = spu_sl(idx, 2);
  idx      = spu_shuffle(idx, idx, VEC_LITERAL(vector unsigned char,
                                               0x03, 0x03, 0x03, 0x03, 0x07, 0x07, 0x07, 0x07,
                                               0x0b, 0x0b, 0x0b, 0x0b, 0x0f, 0x0f, 0x0f, 0x0f));
  idx      = spu_sel(VEC_SPLAT_U32(0x00010203), idx, VEC_SPLAT_U32(0x1c1c1c1c));

  /* Fetch coeeficients */
  glob_a = spu_shuffle(cos18_A_vDSP[0], cos18_A_vDSP[1], (vector unsigned char)(idx));
  glob_b = spu_shuffle(cos18_B_vDSP[0], cos18_B_vDSP[1], (vector unsigned char)(idx));
  glob_c = spu_shuffle(cos18_C_vDSP[0], cos18_C_vDSP[1], (vector unsigned char)(idx));
  glob_d = spu_shuffle(cos18_D_vDSP[0], cos18_D_vDSP[1], (vector unsigned char)(idx));

  ang2 = spu_mul(ang, ang);
  ctmpl = spu_madd(ang, glob_c, glob_d);
  ctmph = spu_madd(ang, glob_a, glob_b);
  result = spu_madd(ang2, ctmph, ctmpl);
  result = (vector float)spu_or((vector unsigned int)result, glob_negflag);
  return (result);
}

int exec_sincosnormalized() {
  //sc.addInVec((void*)argument, argumentcount*sizeof(float));
  //sc.addOutVec((void*)cos_sin_out, argumentcount*sizeof(cf32));
  vector float * phase = (vector float*)args[0];
  vector float * cossin = (vector float*)args[1];
  vector float sx;
  for(int i=0; i<myTask.subcontext.Phase.argcount / 4; i++) {
       sx = _cos_sin18_vrep_MOD(*phase);
//     *cossin++ = _cos_sin18_vDSP(*phase);
//     *phase = spu_rlqw(*phase, 16);
//     *cossin++ = _cos_sin18_vDSP(*phase++);
  }
  write_results();
  return CELLDSP_TASK_OK;
}



#define D_fftsize_vec 1024   // for code debug
#define D_numchannels_vec 8

/*
int exec_vecrepack() 
{
   vector unsigned char regvector; // to speed up calculation
   vector unsigned char indexpointers[D_fftsize_vec / 4];

   vector unsigned char lowernibbles, uppernibbles, evenindex, oddindex;
   vector unsigned char umask, lmask;

   umask = (vector unsigned char) { 0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0, 0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0,0xf0 };
   lmask = (vector unsigned char) { 0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f, 0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f,0x0f };
   //umask = spu_splats(0xf0);
   //lmask = spu_splats(0x0f);

   const vector unsigned char odd2evenFlip = (vector unsigned char) { 1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14 } ;
   const vector unsigned char reinterleave = (vector unsigned char) { 0, 17, 2, 19, 4, 21, 6, 23, 8, 25, 10, 27, 12, 29, 14, 31 };

   // alternate, possibly faster program
   vector unsigned char inter, rotint;

   for (unsigned short i=0; i < D_fftsize_vec/4; i++){

      regvector = indexpointers[i];
      lowernibbles = spu_and(regvector, lmask);
      uppernibbles = spu_and(regvector, umask);
      inter = spu_or ( lowernibbles, spu_shuffle(uppernibbles, uppernibbles, odd2evenFlip) );
      // rotint = spu_rl(inter, 4);  spu_rl not defined for 
//      rotint = spu_rl(inter, 4);
      indexpointers[i] = spu_shuffle(inter, rotint, reinterleave);
   }
   return 0;
}
*/

int exec_vecrepack_versio2() 
{
/*
   vector unsigned short regvector; // to speed up calculation
   vector unsigned short indexpointers[D_fftsize_vec / 4];

   vector unsigned short lowernibbles, uppernibbles, evenindex, oddindex;

   vector unsigned short umask, lmask;
   umask = (vector unsigned short) { 0xf0f0, 0xf0f0, 0xf0f0, 0xf0f0, 0xf0f0, 0xf0f0, 0xf0f0, 0xf0f0 };
   lmask = (vector unsigned short) { 0x0f0f, 0x0f0f, 0x0f0f, 0x0f0f, 0x0f0f, 0x0f0f, 0x0f0f, 0x0f0f };

   // alternate, possibly faster program
   vector unsigned char inter, rotint;
 
   for (unsigned long kk=0; kk < 3125000; kk++)
   for (unsigned short i=0; i < D_fftsize_vec/64; i++){

      regvector = indexpointers[i];

      lowernibbles = spu_and(regvector, lmask); // blank out lower
      uppernibbles = spu_and(regvector, umask); // blank out higher

      vector unsigned short odd_valid  = spu_or ( lowernibbles, spu_sl(lowernibbles, 4) );
      vector unsigned short even_valid = spu_or ( uppernibbles, spu_sl(uppernibbles, 8) );
 
      const vector unsigned char reinterleave = (vector unsigned char) 
         { 0, 16, 2, 18, 4, 20, 6, 22, 8, 24, 10, 26, 12, 28, 14, 30 };

      indexpointers[i] = reinterpret_cast<vector unsigned short> (
                         spu_shuffle(reinterpret_cast<vector unsigned char>(even_valid), 
                                     reinterpret_cast<vector unsigned char>(odd_valid), 
                                     reinterleave) );
   }
*/
   return CELLDSP_TASK_OK;
}

int exec_unpack_improved()
{
/*
   unsigned int packeddata[D_fftsize_vec];
   vector unsigned int * packedvectors;
   vector unsigned int regvector; // to speed up calculation

   unsigned int intermediate[D_numchannels_vec][D_fftsize_vec / 8];

   vector unsigned int trak[16];
   // packedvectors = (vector unsigned int *)packeddata;

   for (unsigned long kk=0; kk < 3125000; kk++)
   for (unsigned short i=0; i < D_fftsize_vec/64; i++){

     // real data:
     // regvector = packedvectors[i];
     // test data:
     regvector = (vector unsigned int) { 0x1A1A, 0x3CC3, 0x4567, 0xFFFF };

     #define TRAK_GATHER(x) trak[x]=spu_gather(spu_rl(regvector,x))    
     trak[0] = spu_gather(regvector);
     TRAK_GATHER(1); // trak1 = spu_gather(spu_rl(regvector,1));
     TRAK_GATHER(2); // ...
     TRAK_GATHER(3);
     TRAK_GATHER(4);
     TRAK_GATHER(5);
     TRAK_GATHER(6);
     TRAK_GATHER(7);
     TRAK_GATHER(8);
     TRAK_GATHER(9);
     TRAK_GATHER(10);
     TRAK_GATHER(11);
     TRAK_GATHER(12);
     TRAK_GATHER(13);
     TRAK_GATHER(14);
     TRAK_GATHER(15);

     #define GET_INTERMED(x) intermediate[x][i] = spu_extract(\
           spu_shuffle(reinterpret_cast<vector unsigned char>(trak[2*x]), reinterpret_cast<vector unsigned char>(trak[2*x+1]),\
                      (vector unsigned char){17,1,17,1,17,1,17,1,17,1,17,1,17,1,17,1}), 0 )

     //intermediate[0][i] = spu_extract((spu_shuffle(trak0, trak1,
     //                                           (){17,1,17,1,17,1,17,1,
     //                                            17,1,17,1,17,1,17,1}) ,0));
     GET_INTERMED(0);
     GET_INTERMED(1);
     GET_INTERMED(2);
     GET_INTERMED(3);
     GET_INTERMED(4);
     GET_INTERMED(5);
     GET_INTERMED(6);
     GET_INTERMED(7);
   }
*/
   return CELLDSP_TASK_OK;
}




//---------------------------------------------------------------------------------
// SPU IMPLEMENTATIONS OF SOME "GARBAGE" DEBUGGING AND TESTING FUNCS
//---------------------------------------------------------------------------------

int do_in_out_test() 
{
   static vector float out_data = { 0.1, 2.3, 4.5, 6.7 };
   // reset these just in case
   myTask.num_input_vectors = 0;
   myTask.num_output_vectors = 1;
   *((vector float*)args[0]) = out_data;
   // write out our struct
   write_results();
   return CELLDSP_TASK_OK;
}

int do_in_out_bench() 
{
/*
   unsigned int length = myTask.outputvector_lengths[0];
   #if XPU_64BIT
   unsigned long long outptr = myTask.outputvector_ptrs[0].ull;
   #else
   unsigned long long outptr = myTask.outputvector_ptrs[0].ui[0];
   #endif
   void * lsaddr = args[0];
   // get both buffers started
   mfc_put(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
   mfc_put(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
   // continue double-buffered
   for (register unsigned short i=0; i<65535; i++) {
      mfc_write_tag_mask(DMA_TAG_WRITEABLES);
      mfc_read_tag_status_all();
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0); // 1
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0); // 5
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0); // 9
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0); // 13
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES, 0, 0); // 16
      mfc_write_tag_mask(DMA_TAG_WRITEABLES_DB);
      mfc_read_tag_status_all();
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0); // 1
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0); // 5
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0); // 9
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0); // 13
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0);
      mfc_get(lsaddr, outptr, length, DMA_TAG_WRITEABLES_DB, 0, 0); // 16
   }
   mfc_write_tag_mask(DMA_TAG_MASK_WRITEABLES_DB + DMA_TAG_WRITEABLES);
   mfc_read_tag_status_all();
*/
   return CELLDSP_TASK_OK;
}


//---------------------------------------------------------------------------------
// do_io_simulation() : run all arithmetic funcs on dummy data
//---------------------------------------------------------------------------------

int do_io_simulation() 
{
  struct dma_list_elem {
     union {
        unsigned int all32;
        struct {
           unsigned nbytes: 31;
           unsigned stall: 1;
        } bits;
     } size;
     unsigned int ea_low;
  };

  struct dma_list_elem list[8] __attribute__ ((aligned (8)));

  /* actual big output buf, and small input buf */
  void * big_buffer = malloc_align(8*8192, MALLOC_QUADALIGN);
  if (big_buffer == NULL) {
     printf("Yeah right, ran out of memory...\n");
     return 0;
  }
  printf("do_io_simulation: 8*8kB buffer allocated\n");

  /* prepare the output DMA-List */
  unsigned int listsize  = 0;
  for (unsigned short kk=0; kk<8 && kk<myTask.num_output_vectors; kk++) {
     listsize += sizeof(dma_list_elem);
     list[kk].size.all32 = myTask.outputvector_lengths[kk]; // how much data
     list[kk].size.bits.stall = 0;
     #if XPU_64BIT
     list[kk].ea_low = myTask.outputvector_ptrs[kk].ull & 0xFFFFFFFF; // where to
     #else
     list[kk].ea_low = myTask.outputvector_ptrs[kk].ui[0]; // where to
     #endif
     printf("  DMA list  item %d ea_low=0x%08X len=0x%X\n", kk, list[kk].ea_low, list[kk].size.all32);
  }
  printf("  Final listsize = %d, sizeof(dma_list_elem) = %ld\n", listsize, sizeof(dma_list_elem));

  /* do the DMA's: one in, 8 times a DMA-List of 8 eight elements out */
  for (unsigned short blocks=0; blocks<5120; blocks++) {
      printf("iteration %d - in start\n" , blocks);

      // dma in
      #if XPU_64BIT
      mfc_getf(args[0], myTask.inputvector_ptrs[0].ull, myTask.inputvector_lengths[0], DMA_TAG_READABLES, 0, 0);
      #else
      mfc_getf(args[0], myTask.inputvector_ptrs[0].ui[0], myTask.inputvector_lengths[0], DMA_TAG_READABLES, 0, 0);
      #endif
      mfc_write_tag_mask(DMA_TAG_MASK_READABLES);
      mfc_read_tag_status_all();
      printf("iteration %d - in done, out start \n" , blocks);

      // dma gathered scatter out 8 x 8kB linked list
      for (unsigned short jj=0; jj<8; jj++) {
         spu_mfcdma32(big_buffer, (unsigned int) &list[0], 8, DMA_TAG_WRITEABLES, MFC_GETL_CMD);
         //spu_getl(bug_buffer, myTask.outputvector_ptrs[kk].ui[0], list, 
         printf("   list iter %d\n", jj);
         mfc_write_tag_mask(DMA_TAG_MASK_WRITEABLES);
         mfc_read_tag_status_all();
      }
   }

   free_align(big_buffer);
   return CELLDSP_TASK_OK;
}


//---------------------------------------------------------------------------------
// do_runall_test() : run all arithmetic funcs on dummy data
//---------------------------------------------------------------------------------

int do_runall_test() 
{
  init_FFT();

  #define ALL_NUMFREQS 2
  #define ALL_NUMITER  2500000
  #define ALL_NUMCHANNELS 64

  printf("SPE: do_runall_test: numfreqs=%d numchannels=%d iterations=%d\n", ALL_NUMFREQS, ALL_NUMCHANNELS, ALL_NUMITER);

  for (unsigned int iter=0; iter<ALL_NUMITER; iter++) {
    //if((iter%10000)==0) { printf("iteration %ld\n", iter); }

  // delay interpolation
  exec_sincos(); // close enough to 2nd order interpolate

  for(unsigned short numfreqs=0; numfreqs<ALL_NUMFREQS; numfreqs++) {
  //printf("freq %d\n" , numfreqs);
     // fringe rotation prepare (delays)
     do_complex_mac(); // close enough to freq multiply for currentchannelfreqptr
     do_complex_mac(); // close enough to freq multiply for xval and fringedelayarray
     exec_sincos(); // rotateargument

     // exec_sincos(); // fracmult, done post-FFT now

     for(unsigned short numchannels=0; numchannels<ALL_NUMCHANNELS; numchannels++) {

        // fringe rotation
        do_complex_mac(); // pre-f rotate

        // transform
        exec_FFT();

        // frac sample correct        
        do_complex_mac(); // complexfracmult, almost like multiply

        // autocorrelation
        do_complex_mac();
     }
     // crosspolar "auto"correlation
     do_complex_mac();
     do_complex_mac();
  }
  }
 
  free_FFT();
  return CELLDSP_TASK_OK;
}

#endif // __SPU__
