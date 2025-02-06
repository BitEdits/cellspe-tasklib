
/***************************************************************************
 *   Copyright (C) 2005 by Adam Deller                                     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/** \file architecture.h
 *  \brief File contains mapping for vector functions to specific architectures
 */
 
#ifndef ARCHITECTURE_H
#define ARCHITECTURE_H

#define INTEL   1
#define AMD     2
#define CELL    3
#define GENERIC 4

//define the MPI tags
#define CR_TERMINATE      0
#define CR_VALIDVIS       1
#define CR_RECEIVETIME    2
#define CR_PROCESSDATA    3
#define CR_PROCESSCONTROL 4
#define DS_TERMINATE      5
#define DS_PROCESS        6

//define the architecture to be compiled for here
//#define ARCH            INTEL
//JanW: now defined in Makefile.am / Makefile

//if no architecture is selected, exit with error
#ifndef ARCH
   #error "No architecture defined in GCC paramenters! Use for example '-D ARCH=GENERIC'."
#endif


// -------------------------------------- INTEL --------------------------------------
// -------------------------------------- INTEL --------------------------------------
// -------------------------------------- INTEL --------------------------------------

//set up the function mapping for the intel architecture
#if(ARCH == INTEL)
#include <ipps.h>
#include <ippvm.h>
#include <ippcore.h>

//start with the data types
#define u8                       Ipp8u
#define u16                      Ipp16u
#define s16                      Ipp16s
#define cs16                     Ipp16sc
#define s32                      Ipp32s
#define f32                      Ipp32f
#define cf32                     Ipp32fc
#define f64                      Ipp64f

//and the constant values
#define vecNoErr                 ippStsNoErr
#define vecFFTSpecR_f32          IppsFFTSpec_R_32f
#define vecFFTSpecC_f32          IppsFFTSpec_C_32f
#define vecFFTSpec_s16           IppsFFTSpec_R_16s
#define vecHintAlg               IppHintAlgorithm
#define vecRndZero               ippRndZero
#define vecRndNear               ippRndNear
#define vecHamming               ippWinHamming
#define vecTrue                  ippTrue
#define MAX_S32                  IPP_MAX_32S
#define MAX_S16                  IPP_MAX_16S
#define MAX_U16                  IPP_MAX_16U
#define vecFFT_NoReNorm          IPP_FFT_NODIV_BY_ANY
#define vecAlgHintFast           ippAlgHintFast
#define TWO_PI                   IPP_2PI

//now the vector functions themselves
#define vectorAlloc_u8(length)   ippsMalloc_8u(length)
#define vectorAlloc_s16(length)  ippsMalloc_16s(length)
#define vectorAlloc_cs16(length) ippsMalloc_16sc(length)
#define vectorAlloc_s32(length)  ippsMalloc_32s(length)
#define vectorAlloc_f32(length)  ippsMalloc_32f(length)
#define vectorAlloc_cf32(length) ippsMalloc_32fc(length)
#define vectorAlloc_f64(length)  ippsMalloc_64f(length)

#define vectorFree(memptr)       ippsFree(memptr)

#define vectorAdd_f32_I(src, srcdest, length)                               ippsAdd_32f_I(src, srcdest, length)
#define vectorAdd_s16_I(src, srcdest, length)                               ippsAdd_16s_I(src, srcdest, length)
#define vectorAdd_s32_I(src, srcdest, length)                               ippsAdd_32s_ISfs(src, srcdest, length, 0)
#define vectorAdd_cf32_I(src, srcdest, length)                              ippsAdd_32fc_I(src, srcdest, length)
#define vectorAddC_f64(src, val, dest, length)                              ippsAddC_64f(src, val, dest, length)
#define vectorAddC_f32(src, val, dest, length)                              ippsAddC_32f(src, val, dest, length)
#define vectorAddC_f32_I(val, srcdest, length)                              ippsAddC_32f_I(val, srcdest, length)
#define vectorAddC_s16_I(val, srcdest, length)                              ippsAddC_16s_I(val, srcdest, length)
#define vectorAddC_f64_I(val, srcdest, length)                              ippsAddC_64f_I(val, srcdest, length)

#define vectorAddProduct_cf32(src1, src2, accumulator, length)              ippsAddProduct_32fc(src1, src2, accumulator, length)

#define vectorConj_cf32(src, dest, length)                                  ippsConj_32fc(src, dest, length)
#define vectorConjFlip_cf32(src, dest, length)                              ippsConjFlip_32fc(src, dest, length)

#define vectorCopy_u8(src, dest, length)                                    ippsCopy_8u(src, dest, length)
#define vectorCopy_s16(src, dest, length)                                   ippsCopy_16s(src, dest, length)
#define vectorCopy_s32(src, dest, length)                                   ippsCopy_32f((f32*)src, (f32*)dest, length)
#define vectorCopy_f32(src, dest, length)                                   ippsCopy_32f(src, dest, length)
#define vectorCopy_cf32(src, dest, length)                                  ippsCopy_32fc(src, dest, length)
#define vectorCopy_f64(src, dest, length)                                   ippsCopy_64f(src, dest, length)

#define vectorCos_f32(src, dest, length)                                    ippsCos_32f_A11(src, dest, length)

#define vectorConvertScaled_s16f32(src, dest, length, scalefactor)          ippsConvert_16s32f_Sfs(src, dest, length, scalefactor)
#define vectorConvertScaled_f32s16(src, dest, length, rndmode, scalefactor) ippsConvert_32f16s_Sfs(src, dest, length, rndmode, scalefactor)
#define vectorConvertScaled_f32u8(src, dest, length, rndmode, scalefactor)  ippsConvert_32f8u_Sfs(src, dest, length, rndmode, scalefactor)
#define vectorConvert_f32s32(src, dest, length, rndmode)                    ippsConvert_32f32s_Sfs(src, dest, length, rndmode, 0)
#define vectorConvert_s16f32(src, dest, length)                             ippsConvert_16s32f(src, dest, length)
#define vectorConvert_s32f32(src, dest, length)                             ippsConvert_32s32f(src, dest, length)
#define vectorConvert_f64f32(src, dest, length)                             ippsConvert_64f32f(src, dest, length)

#define vectorDotProduct_f64(src1, src2, length, output)                    ippsDotProd_64f(src1, src2, length, output);

#define vectorInitFFTR_f32(fftspec, order, flag, hint)                      ippsFFTInitAlloc_R_32f(fftspec, order, flag, hint)
#define vectorInitFFTC_f32(fftspec, order, flag, hint)                      ippsFFTInitAlloc_C_32f(fftspec, order, flag, hint)
//#define vectorInitFFT_f32(fftspec, order, flag, hint)                       ippsDFTInitAlloc_R_32f(fftspec, order, flag, hint)
#define vectorInitFFT_s16(fftspec, order, flag, hint)                       ippsFFTInitAlloc_R_16s(fftspec, order, flag, hint)
#define vectorGetFFTBufSizeR_f32(fftspec, buffersize)                       ippsFFTGetBufSize_R_32f(fftspec, buffersize)
#define vectorGetFFTBufSizeC_f32(fftspec, buffersize)                       ippsFFTGetBufSize_C_32f(fftspec, buffersize)
//#define vectorGetFFTBufSize_f32(fftspec, buffersize)                        ippsDFTGetBufSize_R_32f(fftspec, buffersize)
#define vectorGetFFTBufSize_s16(fftspec, buffersize)                        ippsFFTGetBufSize_R_16s(fftspec, buffersize)
#define vectorFreeFFTR_f32(fftspec)                                         ippsFFTFree_R_32f(fftspec)
#define vectorFreeFFTC_f32(fftspec)                                         ippsFFTFree_C_32f(fftspec)
//#define vectorFreeFFT_f32(fftspec)                                          ippsDFTFree_R_32f(pFFTSpec)
#define vectorFFT_RtoC_f32(src, dest, fftspec, fftbuffer)                   ippsFFTFwd_RToCCS_32f(src, dest, fftspec, fftbuffer)
#define vectorFFT_CtoC_f32(srcre, srcim, destre, destim, fftspec, fftbuff)  ippsFFTFwd_CToC_32f(srcre, srcim, destre, destim, fftspec, fftbuff);
    // { fprintf(stderr, "  genericFFTexec_C Complex, len %d, in ", spec->order); for(int ijk=0;ijk<8; ijk++) { printf("re=%f ", srcre[ijk]); } printf("\n");  }

//#define vectorFFT_RtoC_f32(src, dest, fftspec, fftbuffer)                   ippsDFTFwd_RToCCS_32f(src, dest, fftspec, fftbuffer)
#define vectorScaledFFT_RtoC_s16(src, dest, fftspec, fftbuffer, scale)      ippsFFTFwd_RToCCS_16s_Sfs(src, dest, fftspec, fftbuffer, scale)

#define vectorFlip_f64_I(srcdest, length)                                   ippsFlip_64f_I(srcdest, length)

#define vectorMagnitude_cf32(src, dest, length)                             ippsMagnitude_32fc(src, dest, length)

#define vectorMean_cf32(src, length, mean, hint)                            ippsMean_32fc(src, length, mean, hint)

#define vectorMul_f32(src1, src2, dest, length)                             ippsMul_32f(src1, src2, dest, length)
#define vectorMul_f32_I(src, srcdest, length)                               ippsMul_32f_I(src, srcdest, length)
#define vectorMul_cf32_I(src, srcdest, length)                              ippsMul_32fc_I(src, srcdest, length)
#define vectorMul_cf32(src1, src2, dest, length)                            ippsMul_32fc(src1, src2, dest, length)
#define vectorMulC_f32(src, val, dest, length)                              ippsMulC_32f(src, val, dest, length)
#define vectorMulC_cs16_I(val, srcdest, length)                             ippsMulC_16sc_ISfs(val, srcdest, length, 0)
#define vectorMulC_f32_I(val, srcdest, length)                              ippsMulC_32f_I(val, srcdest, length)
#define vectorMulC_cf32_I(val, srcdest, length)                             ippsMulC_32fc_I(val, srcdest, length)
#define vectorMulC_f64_I(val, srcdest, length)                              ippsMulC_64f_I(val, srcdest, length)
#define vectorMulC_f64(src, val, dest, length)                              ippsMulC_64f(src, val, dest, length)

#define vectorPhase_cf32(src, dest, length)                                 ippsPhase_32fc(src, dest, length)

#define vectorRealToComplex_f32(real, imag, complex, length)                ippsRealToCplx_32f(real, imag, complex, length)

#define vectorSet_f32(val, dest, length)                                    ippsSet_32f(val, dest, length)

#define vectorSin_f32(src, dest, length)                                    ippsSin_32f_A11(src, dest, length)

#define vectorSinCos_f32(src, sin, cos, length)                             ippsSinCos_32f_A11(src, sin, cos, length)

#define vectorSplitScaled_s16f32(src, dest, numchannels, chanlen)           ippsSplitScaled_16s32f_D2L(src, dest, numchannels, chanlen)

#define vectorSquare_f64_I(srcdest, length)                                 ippsSqr_64f_I(srcdest, length)

#define vectorSub_f32_I(src, srcdest, length)                               ippsSub_32f_I(src, srcdest, length)
#define vectorSub_s32(src1, src2, dest, length)                             ippsSub_32s_Sfs(src1, src2, dest, length, 0)
#define vectorSub_cf32_I(src, srcdest, length)                              ippsSub_32fc_I(src, srcdest, length)

#define vectorGenerateFIRLowpass_f64(freq, taps, length, window, normalise) ippsFIRGenLowpass_64f(freq, taps, length, window, normalise)

#define vectorZero_cf32(dest, length)                                       ippsZero_32fc(dest, length)
#define vectorZero_f32(dest, length)                                        ippsZero_32f(dest, length)
#define vectorZero_s16(dest, length)                                        ippsZero_16s(dest, length)
#define vectorZero_s32(dest, length)                                        genericZero_32s(dest, length)

inline int genericZero_32s(s32 * dest, int length)
{ for(int i=0;i<length;i++) dest[i] = 0; return vecNoErr; }

#endif /*Architecture == Intel */





// ----------------------- GENERIC CPU - CAN OVERRIDE FURTHER BELOW --------------------
// ----------------------- GENERIC CPU - CAN OVERRIDE FURTHER BELOW --------------------
// ----------------------- GENERIC CPU - CAN OVERRIDE FURTHER BELOW --------------------

#if (ARCH == GENERIC) || (ARCH == CELL)

#if (ARCH==GENERIC)
  #include <cmath>  // IBM Cell SDK2.0 bug: include cmath after ppu_intrinsics.h
#else
  #include <ppu_intrinsics.h>
  #include <cmath>
#endif

#include <fftw3.h>       // use ffwtf_*** for single precision to match f32, fc32 datatypes
                         // NOTE: FFTW has to be './configure --enable-float'ed before make and install

#include <malloc.h>         // memalign() and free()
#include "types_generic.h"  // data types for Generic and Cell

//vector allocation and deletion routines
template<class dt> inline dt* vectorAlloc_XX(size_t length) {
  return (dt*)memalign(128, length*sizeof(dt));
}
#define vectorAlloc_u8(length)   vectorAlloc_XX<u8>(length)
#define vectorAlloc_s16(length)  vectorAlloc_XX<s16>(length)
#define vectorAlloc_cs16(length) vectorAlloc_XX<cs16>(length)
#define vectorAlloc_s32(length)  vectorAlloc_XX<s32>(length)
#define vectorAlloc_f32(length)  vectorAlloc_XX<f32>(length)
#define vectorAlloc_cf32(length) vectorAlloc_XX<cf32>(length)
#define vectorAlloc_f64(length)  vectorAlloc_XX<f64>(length)
#define vectorFree(memptr)       free(memptr)

//then the generic functions
template<class dt> inline int genericCopy_XX(dt * src, dt * dest, int length)
{
   size_t dtsize = sizeof(dt);
   // memcpy((void*)dest, (void*)src, dtsize * length); // C++ finds this not nice
   for(size_t i=0; i<length; i++) *(dest+i) = *(src+i);
   return vecNoErr;
}
#define vectorCopy_u8(src, dest, length)       genericCopy_XX<u8>(src, dest, length)
#define vectorCopy_s16(src, dest, length)      genericCopy_XX<s16>(src, dest, length)
#define vectorCopy_f32(src, dest, length)      genericCopy_XX<f32>(src, dest, length)
#define vectorCopy_f64(src, dest, length)      genericCopy_XX<f64>(src, dest, length)

#define vectorZero_s32(dest, length)           genericZero_32s(dest, length)
inline int genericZero_32s(s32 * dest, int length)
   { for(int i=0;i<length;i++) dest[i] = 0; return vecNoErr; }

#define vectorZero_cf32(dest, length)          genericZero_32cf(dest, length)
inline int genericZero_32cf(fc32 * dest, int length)    
   { for(int i=0;i<length;i++) {dest[i].re = 0; dest[i].im=0;} return vecNoErr; }



// vector by vector multiply, real
inline int genericMul_f32(f32 * src1, f32 * src2, f32 * dest, int length) 
   {  for(int i=0;i<length;i++) dest[i]=src1[i]*src2[i]; return vecNoErr; }
// #define vectorMul_f32(src1, src2, dest, length)                ippsMul_32f(src1, src2, dest, length)
#define vectorMul_f32(src1, src2, dest, length)                 genericMul_f32(src1, src2, dest, length)
// #define vectorMul_f32_I(src, srcdest, length)                  ippsMul_32f_I(src, srcdest, length)
#define vectorMul_f32_I(src, srcdest, length)                   genericMul_f32(src, srcdest, srcdest, length)


// #define vectorMul_cf32(src1, src2, dest, length)               ippsMul_32fc(src1, src2, dest, length)
#define vectorMul_cf32(src1, src2, dest, length)                genericMul_32fc(src1, src2, dest, length)
// vector by vector multiply, complex
inline int genericMul_32fc(fc32 * src, fc32 * src2, fc32 * dest, int length) 
{
  for(int i=0;i<length;i++) {
     dest[i].re = (src[i].re * src2[i].re) - (src[i].im * src2[i].im);
     dest[i].im = (src[i].im * src2[i].re) + (src[i].re * src2[i].im);
  } 
  return vecNoErr; 
}


#define vectorMul_cf32_conj(src1, src2, dest, length)                genericMul_32fc_conj(src1, src2, dest, length)
// vector by vector multiply, complex with second argument to be conjugated
inline int genericMul_32fc_conj(fc32 * src, fc32 * src2, fc32 * dest, int length) 
{
  for(int i=0;i<length;i++) {
     dest[i].re = (src[i].re * src2[i].re) + (src[i].im * src2[i].im);
     dest[i].im = (src[i].im * src2[i].re) - (src[i].re * src2[i].im);
  } 
  return vecNoErr; 
}


// #define vectorMul_cf32_I(src, srcdest, length)                 ippsMul_32fc_I(src, srcdest, length)
#define vectorMul_cf32_I(src1, srcdest, length)                 genericMul_32fc_I(src1, srcdest, length)
inline int genericMul_32fc_I(fc32 * src, fc32 * src2, int length) 
{
  // inplace version, src2=dest
  fc32 tmp;
  for(int i=0;i<length;i++) {
     tmp.re = (src[i].re * src2[i].re) - (src[i].im * src2[i].im);
     tmp.im = (src[i].im * src2[i].re) + (src[i].re * src2[i].im);
     src2[i].re = tmp.re;
     src2[i].im = tmp.im;
  } 
  return vecNoErr; 
}


// vector scale
template<class dt> inline int genericMulC_XXf(dt * src, dt val, dt * dest, int length)
   { for(int i=0;i<length;i++) dest[i] = val * src[i]; return vecNoErr; }
// #define vectorMulC_f32(src, val, dest, length)                 ippsMulC_32f(src, val, dest, length)
#define vectorMulC_f32(src, val, dest, length)                  genericMulC_XXf<f32>(src, val, dest, length)
// #define vectorMulC_f32_I(val, srcdest, length)                 ippsMulC_32f_I(val, srcdest, length)
#define vectorMulC_f32_I(val, srcdest, length)                  genericMulC_XXf<f32>(srcdest, val, srcdest, length)
// #define vectorMulC_f64(src, val, dest, length)                 ippsMulC_64f(src, val, dest, length)
#define vectorMulC_f64(src, val, dest, length)                  genericMulC_XXf<f64>(src, val, dest, length)
// #define vectorMulC_f64_I(val, srcdest, length)                 ippsMulC_64f_I(val, srcdest, length)
#define vectorMulC_f64_I(val, srcdest, length)                  genericMulC_XXf<f64>(srcdest, val, srcdest, length)

// vector scale by complex
//#define vectorMulC_cs16_I(val, srcdest, length)                  ippsMulC_16sc_ISfs(val, srcdest, length, 0)
 
// vector add vector
template<class dt> inline int genericAdd_XX_I(dt * src, dt * srcdest, int length)
   { for(int i=0;i<length;i++) srcdest[i] += src[i]; return vecNoErr; }
// #define vectorAdd_s32_I(src, srcdest, length)                  ippsAdd_32s_ISfs(src, srcdest, length, 0)
#define vectorAdd_s32_I(src, srcdest, length)                   genericAdd_XX_I<s32>(src, srcdest, length)
// #define vectorAdd_f32_I(src, srcdest, length)                  ippsAdd_32f_I(src, srcdest, length)
#define vectorAdd_f32_I(src, srcdest, length)                   genericAdd_XX_I<f32>(src, srcdest, length)
// #define vectorAdd_s16_I(src, srcdest, length)                  ippsAdd_16s_I(src, srcdest, length)
#define vectorAdd_s16_I(src, srcdest, length)                   genericAdd_XX_I<s16>(src, srcdest, length)

// complex vector add complex vector
inline int genericAdd_32fc_I(fc32 * src, fc32 * srcdest, int length)           
   { for(int i=0;i<length;i++) { srcdest[i].re += src[i].re; srcdest[i].im += src[i].im;} return vecNoErr; }       
// #define vectorAdd_cf32_I(src, srcdest, length)                 ippsAdd_32fc_I(src, srcdest, length) 
#define vectorAdd_cf32_I(src, srcdest, length)                  genericAdd_32fc_I(src, srcdest, length)

// vector add constant
template<class dt> inline int genericAddC_XX(dt * src, dt val, dt * dest, int length)
   { for(int i=0;i<length;i++) { dest[i] = src[i] + val; } return vecNoErr; }
// #define vectorAddC_f64(src, val, dest, length)                 ippsAddC_64f(src, val, dest, length)
#define vectorAddC_f64(src, val, dest, length)                  genericAddC_XX<f64>(src, val, dest, length)
// #define vectorAddC_f32(src, val, dest, length)                 ippsAddC_32f(src, val, dest, length)
#define vectorAddC_f32(src, val, dest, length)                  genericAddC_XX<f32>(src, val, dest, length)
// #define vectorAddC_f32_I(val, srcdest, length)                 ippsAddC_32f_I(val, srcdest, length)
#define vectorAddC_f32_I(val, srcdest, length)                  genericAddC_XX<f32>(srcdest, val, srcdest, length)
//#define vectorAddC_s16_I(val, srcdest, length)                  ippsAddC_16s_I(val, srcdest, length)
#define vectorAddC_s16_I(val, srcdest, length)                  genericAddC_XX<s16>(srcdest, val, srcdest, length)
//#define vectorAddC_f64_I(val, srcdest, length)                  ippsAddC_64f_I(val, srcdest, length)
#define vectorAddC_f64_I(val, srcdest, length)                  genericAddC_XX<f64>(srcdest, val, srcdest, length)


#define vectorConj_cf32(src, dest, length)        genericConj_32fc(src, dest, length)

// vector conjugate
inline int genericConj_32fc(fc32 * src, fc32 * dest, int length) 
{
  for(int i=0;i<length;i++) {
     dest[i].re = src[i].re;
     dest[i].im = src[i].im * (-1);
  } 
  return vecNoErr; 
}

#define vectorSquare_f64_I(srcdest, length)            genericVectorSqr_64f_I(srcdest, length)

inline int genericVectorSqr_64f_I(f64 * srcdest, int length)
{ for(int i=0;i<length;i++) srcdest[i] *= srcdest[i]; return vecNoErr; }

#define vectorRealToComplex_f32(real, imag, complex, length)    genericRealToCplx_32f(real, imag, complex, length)

// Merge real and imaginary parta to one
inline int genericRealToCplx_32f(f32 * real, f32 * imag, fc32 * complex, int length) 
{
  for(int i=0;i<length;i++) {
     complex[i].re = real[i];
     complex[i].im = imag[i];
  } 
  return vecNoErr; 
}

#define vectorConjFlip_cf32(src, dest, length)    genericConjFlip_32fc(src, dest, length)

// vector conjugate and flip
inline int genericConjFlip_32fc(fc32 * src, fc32 * dest, int length) 
{
  for(int i=0;i<length;i++) {
     dest[length - i - 1].re = src[i].re;
     dest[length - i - 1].im = src[i].im * (-1);
  } 
  return vecNoErr; 
}

#define vectorDotProduct_f64(src1, src2, length, output)    genericDotProduct_64f(src1, src2, length, output);
inline int genericDotProduct_64f(f64 * src, f64 * src2, int length, double * output)
{
  *output = 0;
  for(int i=0;i<length;i++) *output += src[i] * src2[i];
  return vecNoErr;
}

#define vectorAddProduct_cf32(src, src2, accumulator, length) \
        genericAddProduct_32fc(src, src2, accumulator, length)
inline int genericAddProduct_32fc(fc32 * src, fc32 * src2, fc32 * accumulator, int length)
{ 
  for(int i=0;i<length;i++) {
     /*
     // Vector processor optimized way, one multiplication exchanged for an addition 
     // and a substraction, idea stolen from some IBM's Cell presentation
     // (a+ib)*(c+id) = (ac-bd) + i(ad+bc) = (ac-bd) + i( (b+c)*(a+d) - ac - bd)
     ac = src[i].re * src2[i].re;
     bd = src[i].im * src2[i].im;
     accumulator[i].re += ac - bd;
     accumulator[i].im += (src[i].im + src2[i].re) * (src[i].re + src2[i].im) - ac - bd;
     */
     // the more natural way, more efficient on generic CPU:  
     accumulator[i].re += (src[i].re * src2[i].re) - (src[i].im * src2[i].im);
     accumulator[i].im += (src[i].im * src2[i].re) + (src[i].re * src2[i].im);
  }
  return vecNoErr;
}

#define vectorAddProduct_cf32_conj(src, src2, accumulator, length) \
        genericAddProduct_32fc_conj(src, src2, accumulator, length)
// Conjugate src2, then do MAC calculation acc += src*(conj src2) += (ac+bd) + i(bc+ad)
inline int genericAddProduct_32fc_conj(fc32 * src, fc32 * src2, fc32 * accumulator, int length)
{
  for(int i=0;i<length;i++) {
     accumulator[i].re += (src[i].re * src2[i].re) + (src[i].im * src2[i].im);
     accumulator[i].im += (src[i].im * src2[i].re) - (src[i].re * src2[i].im);
  }
  return vecNoErr;
}

#define vectorAddProduct_cf32_conjMul(src, src2, mult, accumulator, length) \
        genericAddProduct_cf32_conjMul(src, src2, mult, accumulator, length)
// Multiply src, src2 by mult and store results back. Then do MAC into accumulator
// with included conjugation of src2. (rescaled cross-correlation)
inline int genericAddProduct_32fc_conjMul(fc32 * src, fc32 * src2, fc32 * mult, fc32 * accumulator, int length)
{
  fc32 tmp;
  for(int i=0;i<length;i++) {
     tmp.re = (src[i].re * mult[i].re) - (src[i].im * mult[i].im);
     tmp.im = (src[i].im * mult[i].re) + (src[i].re * mult[i].im);
     src[i].re = tmp.re;
     src[i].im = tmp.im;  
     tmp.re = (src2[i].re * mult[i].re) - (src2[i].im * mult[i].im);
     tmp.im = (src2[i].im * mult[i].re) + (src2[i].re * mult[i].im);
     src2[i].re = tmp.re;
     src2[i].im = tmp.im;  
     accumulator[i].re += (src[i].re * src2[i].re) + (src[i].im * src2[i].im);
     accumulator[i].im += (src[i].im * src2[i].re) - (src[i].re * src2[i].im);
  }
  return vecNoErr;
}

#define vectorAddProduct_cf32_conjMul_self(src, mult, accumulator, length) \
        genericAddProduct_cf32_conjMul_self(src, src2, mult, accumulator, length)
// Multiply src by mult and store result back. Then do MAC into accumulator
// with included conjugation of src. (rescaled auto-correlation)
inline int genericAddProduct_cf32_conjMul_self(fc32 * src, fc32 * mult, fc32 * accumulator, int length)
{
  fc32 tmp;
  for(int i=0;i<length;i++) {
     tmp.re = (src[i].re * mult[i].re) - (src[i].im * mult[i].im);
     tmp.im = (src[i].im * mult[i].re) + (src[i].re * mult[i].im);
     src[i].re = tmp.re;
     src[i].im = tmp.im;  
     accumulator[i].re += (src[i].re * src[i].re) + (src[i].im * src[i].im);
     // accumulator[i].im += (src[i].im * src[i].re) - (src[i].re * src[i].im); // = 0
  }
  return vecNoErr;
}

#define vectorMean_cf32(src, length, mean, hint)     genericMean_32fc(src, length, mean, hint)
inline int genericMean_32fc(fc32 * src, int length, fc32 * mean, IppHintAlgorithm hint)  
{
   mean->re = 0; 
   mean->im = 0;
   if (length<=0) { return !vecNoErr; }
   for(int i=0;i<length;i++) { 
      mean->re += src[i].re; 
      mean->im += src[i].im; 
   }
   mean->re /= length; 
   mean->im /= length;
   return vecNoErr;
}

#define vectorSinCos_f32(src, dstsin, dstcos, length)                             genericSinCos_32f(src, dstsin, dstcos, length)
inline int genericSinCos_32f(f32 * src, f32 * dstsin, f32 * dstcos, int length) 
{ 
  // IPP docs: "function flavor ippsSinCos_32f_A11 guarantees 11 correctly rounded bits of significand, or at least 3 exact decimal digits;"
  for(int i=0;i<length;i++) {
     dstcos[2*i] = f32(cos(src[i]));
     dstcos[2*i+1] = f32(sin(src[i]));
  } 
  return vecNoErr; 
}

#define vectorSplitScaled_s16f32(src, dest, numchannels, chanlen)  genericSplitScaled_16s32f_D2L(src, dest, numchannels, chanlen)
inline int genericSplitScaled_16s32f_D2L(s16 * pSrc, f32 ** ppDst, int nChannels, int chanLen) 
{
  int linidx = 0;
  for (int j=0;j<chanLen; j++) {
     for (int i=0;i<nChannels;i++) {
        ppDst[i][j] = f32(pSrc[linidx++]) / f32(32768); // ipps.h: Region of the dst data is [-1.0,1.0].
     }
  }
  return vecNoErr;
}



// Fast Fourier functions

#define vectorInitFFTR_f32(fftspec, order, flag, hint)     genericFFTInitAlloc_X_32f<vecFFTSpecR_f32>(fftspec, order, flag, hint)
#define vectorInitFFTC_f32(fftspec, order, flag, hint)     genericFFTInitAlloc_X_32f<vecFFTSpecC_f32>(fftspec, order, flag, hint)
template <class dt>
inline int genericFFTInitAlloc_X_32f(dt ** fftspec, int order, int flag, IppHintAlgorithm hint)
{
   // need delayed init since difx wants to use alloc + plan create exactly the wrong way for FFTW...
   dt * spec = new dt;
   spec->in    = NULL;
   spec->out   = NULL;
   spec->order = 1 << order;
   *fftspec = spec;
   return vecNoErr;
}


// ... not sure what 'fftbuffer' in mode.cpp is used for, but it is allocated to GetFFTBufSize() ...
inline int vectorGetFFTBufSizeR_f32(vecFFTSpecR_f32 * fftspec, int * buffersize)
{
   *buffersize = sizeof(fftwf_real) * (fftspec->order);
   return vecNoErr;
}
inline int vectorGetFFTBufSizeC_f32(vecFFTSpecC_f32 * fftspec, int * buffersize)
{
   *buffersize = sizeof(fftwf_complex) * (fftspec->order);
   return vecNoErr;
}

// free up FFT plan on PPU - later on SPE, this should map to kill_thread() or similar
#define vectorFreeFFTR_f32(fftspec)                        genericFreeFFT_f32<vecFFTSpecR_f32>(fftspec)
#define vectorFreeFFTC_f32(fftspec)                        genericFreeFFT_f32<vecFFTSpecC_f32>(fftspec)
template <class dt>
inline int genericFreeFFT_f32(dt * fftspec) {
   //delete fftspec->in;  -- DiFX does these separately..
   //delete fftspec->out;
   fftspec->in = NULL;
   fftspec->out = NULL;
   fftspec->order = 0;
   fftwf_destroy_plan(fftspec->plan); 
   delete fftspec;
   return vecNoErr;
}


#define vectorFFT_RtoC_f32(src, dest, fftspec, fftbuffer) genericFFTexec_R(src, dest, fftspec, fftbuffer) 
//inline int genericFFTexec_R(fftwf_real * src, fftwf_complex * dest, vecFFTSpecR_f32 ** fftRspec, int fftbufsize)
inline int genericFFTexec_R(fftwf_real * src, fftwf_real * dest, vecFFTSpecR_f32 * fftRspec, u8 * fftbufsize)
{
   vecFFTSpecR_f32 * spec = fftRspec;
   if (spec->in == NULL) {
      // delayed init/plan creation
      spec->in   = reinterpret_cast<fftwf_real *>(src);
      spec->out  = reinterpret_cast<fftwf_complex *>(dest);
      spec->plan = fftwf_plan_dft_r2c_1d(spec->order, spec->in, spec->out, FFTW_MEASURE); // FFTW_ESTIMATE
   }
   spec->in  = reinterpret_cast<fftwf_real *>(src);
   spec->out = reinterpret_cast<fftwf_complex *>(dest);
   // using FFTW Guru interface, as src/dest may be different from those in plan
   fftwf_execute_dft_r2c(spec->plan, spec->in, spec->out);   
   *fftbufsize = sizeof(fftwf_complex) * spec->order;
   return vecNoErr;
}

#define vectorFFT_CtoC_f32(srcre, srcim, destre, destim, fftspec, fftbuff) genericFFTexec_C(srcre, destre, fftspec, fftbuff)
inline int genericFFTexec_C(fftwf_real * src, fftwf_real * dest, vecFFTSpecC_f32 * fftCspec, u8 * fftbufsize)
{
   vecFFTSpecC_f32 * spec = fftCspec;
   if (spec->in == NULL) {
      // delayed init/plan creation   
      spec->in   = reinterpret_cast<fftwf_complex *>(src);
      spec->out  = reinterpret_cast<fftwf_complex *>(dest);
      spec->plan = fftwf_plan_dft_1d(spec->order, spec->in, spec->out, FFTW_FORWARD, FFTW_MEASURE); // FFTW_ESTIMATE
   }
   spec->in   = reinterpret_cast<fftwf_complex *>(src);
   spec->out  = reinterpret_cast<fftwf_complex *>(dest);
   //printf(" spec->in : re={%e %e %e %e %e %e)\n", src[2*0], src[2*1], src[2*2], src[2*3], src[2*4], src[2*5]);   
   // using FFTW Guru interface, as src/dest may be different from those in plan
   fftwf_execute_dft(spec->plan, spec->in, spec->out);      
   //printf(" spec->out : re={%e %e %e %e %e %e)\n", dest[2*0], dest[2*1], dest[2*2], dest[2*3], dest[2*4], dest[2*5]);
   *fftbufsize = sizeof(fftwf_complex) * spec->order;
   return vecNoErr;
}


#define flipWord(w) ( u16((w)<<8) | u16((w)>>8) )
// #define flipWord(w) ( (w)<<8 | (w)>>8 )

inline float floatSwap(char *value) {
  char buffer[4];
  buffer[0] = value[3];
  buffer[1] = value[2];
  buffer[2] = value[1];
  buffer[3] = value[0];
  return *( (float *) &buffer );
}

#endif /* Generic Architecture */



// --------------------------------------- CELL ---------------------------------------
// --------------------------------------- CELL ---------------------------------------
// --------------------------------------- CELL ---------------------------------------

#if (ARCH == CELL)

#if USE_CELL_SPE
#include "cellspe-tasklib.h"

#undef vectorInitFFTC_f32
#define vectorInitFFTC_f32(fftspec, order, flag, hint)     speFFTInitAlloc_X_32f<vecFFTSpecC_f32>(fftspec, order, flag, hint)
template <class dt>
inline int speFFTInitAlloc_X_32f(dt ** fftspec, int order, int flag, IppHintAlgorithm hint)
{
   if(*fftspec==NULL) { // allocate new SPE only once
     dt * spec = new dt;
     spec->in    = NULL;
     spec->out   = NULL;
     spec->order = 1 << order;
     speFFT_init(spec->order, &spec->aux); // store SPE id into aux
     *fftspec = spec;
   } else {
     (*fftspec)->order = 1 << order;
   }
   return vecNoErr;
}

#undef vectorFreeFFTC_f32
#define vectorFreeFFTC_f32(fftspec)                         speFreeFFT_f32<vecFFTSpecC_f32>(fftspec)
template <class dt>
inline int speFreeFFT_f32(dt * fftspec) {
   //delete fftspec->in;  -- DiFX does these separately..
   //delete fftspec->out;
   if (fftspec != NULL) { // free up only the first time
     fftspec->in = NULL;
     fftspec->out = NULL;
     fftspec->order = 0;
     speFFT_destroy(&fftspec->aux);
     delete fftspec;
   }
   return vecNoErr;
}

#undef vectorFFT_CtoC_f32
#define vectorFFT_CtoC_f32(srcre, srcim, destre, destim, fftspec, fftbuff) genericFFTexec_C(srcre, destre, fftspec, fftbuff)
inline int speFFTexec_C(fftwf_real * src, fftwf_real * dest, vecFFTSpecC_f32 * fftCspec, u8 * fftbufsize)
{
   speFFT_exec((void*)src, (void*)dest, fftCspec->order, &fftCspec->aux);
   *fftbufsize = sizeof(fftwf_complex) * fftCspec->order;
   return vecNoErr;
}

#endif //USE_CELL_SPE


#endif /* Cell architecture */


// ------------------------------- POWERPC/ALTIVEC ------------------------------------
// ------------------------------- POWERPC/ALTIVEC ------------------------------------
// ------------------------------- POWERPC/ALTIVEC ------------------------------------

#ifdef HAVE_ALTIVEC_H
#include "altivec-macros.h"
#include "altivec-libfreevec.h"

#undef vectorCopy_u8
#define vectorCopy_u8(src, dest, length)       altivec_memcpy(reinterpret_cast<void*>(dest), reinterpret_cast<void*>(src), length)
#undef vectorCopy_s16
#define vectorCopy_s16(src, dest, length)       altivec_memcpy(reinterpret_cast<void*>(dest), reinterpret_cast<void*>(src), 2*length)
#undef vectorCopy_f32
#define vectorCopy_f32(src, dest, length)       altivec_memcpy(reinterpret_cast<void*>(dest), reinterpret_cast<void*>(src), 4*length)
#undef vectorCopy_f64
#define vectorCopy_f64(src, dest, length)       altivec_memcpy(reinterpret_cast<void*>(dest), reinterpret_cast<void*>(src), 8*length)

#endif


#endif /* Defined architecture header */
