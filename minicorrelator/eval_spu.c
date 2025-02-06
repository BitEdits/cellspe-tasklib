/********************************************************************************
 * Fast Matrix Multiplication for Cell BE Processors
 * Copyright (C) 2007  Daniel Hackenberg, ZIH, TU-Dresden
 * A comprehensive description is available at
 * http://tu-dresden.de/die_tu_dresden/zentrale_einrichtungen/zih/forschung/
 * architektur_und_leistungsanalyse_von_hochleistungsrechnern/cell
 * Please send your feedback to daniel.hackenberg@zih.tu-dresden.de
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

#include "minicorrelator.h"
#include <spu_mfcio.h>
#include <simdmath.h>

#include <stdio.h>

#define UNROLL_BY_2(x)  { x }{ x }
#define UNROLL_BY_4(x)  { x }{ x }{ x }{ x }
#define UNROLL_BY_8(x)  UNROLL_BY_4(x) UNROLL_BY_4(x)
#define UNROLL_BY_16(x) UNROLL_BY_8(x) UNROLL_BY_8(x)
#define VEC_4f(x)       (vector float){x, x, x, x}

#define __dma__(VAR1, VAR2) { \
  vlist[VAR1] = v0;           \
  vlist[VAR2] = v1;           \
  v0 = spu_add(v0, add);      \
  v1 = spu_add(v1, add);      \
}

control_block cb              __attribute__ ((aligned (128)));
addr64        spe_ls_Ptrs[16] __attribute__ ((aligned (128)));
int           g_rank;
int           g_targetspu;

vector unsigned char fake_signal;

/* Multibuffering is used to do the calculations efficiently.
 * There are two buffers for matrix A, two buffers for matrix B
 * and four buffers for matrix C.
 * Capital letters (A, B, C) refer to the matrix itself.
 * Lower case letters (a, b, c, d) refer to the buffers in the local store.
 *
 *          B
 *    ---------
 *   |   | a b |
 *   |---------|
 * A | a | a b | C
 *   | b | c d |
 *    ---------
 */


float rawdata_inA[SPE_FIXEDFLOATSIZE] __attribute__ ((aligned (128)));
float rawdata_inB[SPE_FIXEDFLOATSIZE] __attribute__ ((aligned (128)));
float float_out  [SPE_FIXEDFLOATSIZE] __attribute__ ((aligned (128)));
float xcorr_out  [SPE_FIXEDFLOATSIZE] __attribute__ ((aligned (128)));

vector float byte_to_vecfloat[256]; // 2-bit samples to floats lookup (1 byte to 4 float samples), uninitialized for now...

#define TAG_1 1
#define TAG_2 2
#define TAG_3 3
#define TAG_4 4
#define TAG_5 5

// mini....h #define __USE_BDL__ 1 // DMA pilkkominen listaksi jarkevaa vasta ei-2^N && yli 16kB siirroilla
   unsigned int dma_list1in[128];
   unsigned int dma_list2in[128];
   unsigned int dma_list3in[128];
   unsigned int dma_list1fout[128];
   unsigned int dma_list2fout[128];
   unsigned int dma_list3fout[128];
   unsigned int dma_list1xcout[128];
   unsigned int dma_list2xcout[128];
   unsigned int dma_list3xcout[128];

addr64 PPE_rawdata_A_Ptr, PPE_rawdata_B_Ptr, PPE_rawdata_C_Ptr, PPE_rawdata_D_Ptr;
addr64 PPE_speresult_Ptr, SPE_result_Ptr;

unsigned int size, target_spu;
int m_blocks, n_blocks, p_blocks;

/* This vector is used by the assembly-coded matmul function. */
vector unsigned int SIMDadd =  (vector unsigned int){0, 64, 128, 192};


inline addr64 get_output_buffer(void* localbuf) {
   addr64 out;
   #if 1
   g_targetspu = (g_targetspu + 1) % cb.num_spes;
   if (g_rank == g_targetspu)
       g_targetspu = (g_targetspu + 1) % cb.num_spes;
   out = spe_ls_Ptrs[target_spu];
   out.ull += (int)localbuf;
   #else
   out = PPE_speresult_Ptr;
   #endif
   return out;
}

inline void wait_for_mbox(int _tag_id) {
  mfc_write_tag_mask(1 << _tag_id);
  mfc_read_tag_status_all();
  return;
}

inline void load_data(addr64 _PPE_Ptr, float *_feld, unsigned int *_dma_list, int _blocknr, int _tag) {
#ifdef __USE_BDL__
  mfc_get(_feld, _PPE_Ptr.ull + (unsigned long long)(_blocknr * SPE_FIXEDBUFSIZE), SPE_FIXEDBUFSIZE, _tag, 0, 0);
#else
  fill_dma_list(_PPE_Ptr.ui[1], _dma_list, _offset_m, _offset_n);
  mfc_getl(_feld, _PPE_Ptr.ull, _dma_list, DMA_LIST_SIZE, _tag, 0, 0);
#endif
}

inline void store_data(addr64 _PPE_Ptr, float *_feld, unsigned int *_dma_list, int _blocknr, int _tag) {
#ifdef __USE_BDL__
  mfc_put(_feld, _PPE_Ptr.ull + (unsigned long long)(_blocknr * SPE_FIXEDBUFSIZE), SPE_FIXEDBUFSIZE, _tag, 0, 0);
#else
  fill_dma_list(_PPE_Ptr.ui[1], _dma_list, _offset_m, _offset_n);
  mfc_putl(_feld, _PPE_Ptr.ull, _dma_list, DMA_LIST_SIZE, _tag, 0, 0);
#endif
}

void print_floatvecs(vector float* in, unsigned int numvecs)
{
    float f1, f2, f3, f4;
    unsigned int i;
    for (i=0; i<numvecs; i++) {
        f1 = spu_extract(*in, 0);
        f2 = spu_extract(*in, 1);
        f3 = spu_extract(*in, 2);
        f4 = spu_extract(*in, 3);
        printf("%04d [%+f \t %+f \t  %+f \t %+f]\n", i, f1, f2, f3, f4);
        in++;
    }
    return;
}
void print_floatvecs_matlab(char* name, vector float* in, unsigned int numvecs)
{
    float f1, f2, f3, f4;
    unsigned int i;
    printf(" %s = [ ", name);
    for (i=0; i<numvecs; i++) {
        f1 = spu_extract(*in, 0);
        f2 = spu_extract(*in, 1);
        f3 = spu_extract(*in, 2);
        f4 = spu_extract(*in, 3);
        printf("%f %f %f %f ", f1, f2, f3, f4);
        in++;
    }
    printf("];\n");
    return;
}

int main(unsigned long long speid, addr64 argp, addr64 envp) {

  int m, n, p = 0;
  unsigned int t_start = 0, t_spu, count = 0, totaldmascount = 0;
  unsigned char * raw;
  vector float * out;
  int i, j;

  /* start the SPE-side measurement of the execution time */
  spu_write_decrementer(0x7fffffff);
  t_start = spu_read_decrementer();

  fake_signal = (vector unsigned char) { 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 };
  m = SPE_FIXEDBUFSIZE;
  n = SPE_FIXEDBUFSIZE;
  g_rank = 0;

  // ---------------- INSERT CODE HERE
  // ---------------- INSERT CODE HERE
  #if 1
  {
          // WARN: angles must be -pi .. +pi, normalized into -1.0 .. +1.0 !!
          /*
             phi=linspace(-pi, pi, 64000);
             phinorm = phi ./ pi;
             trueSine = sin(phi);
             trueCosine = cos(phi);
             coeffsS = polyfit(phinorm, trueSine, 12);
             coeffsC = polyfit(phinorm, trueCosine, 12);
             figure(1), plot(phinorm, polyval(coeffsS,phinorm)-trueSine, 'g-');
             figure(2), plot(phinorm, polyval(coeffsC,phinorm)-trueCosine, 'g-');
          */
          vector float *phase  = (vector float*) rawdata_inA; 
          vector float *sinout = (vector float*) float_out;
          vector float *cosout = (vector float*) xcorr_out;
          vector float stmp, ctmp, phi2, phi;
          vector signed int iphase;

          {
             vector float inc = VEC_4f(2.0 / SPE_FIXEDVECSIZE), sum;
             sum = (vector float){-1.0, -1.0, -1.0, -1.0};
             for (i=0; i<SPE_FIXEDVECSIZE; i++) { // complex nums/2 = vectors
                 sum = spu_add(sum, inc);
                 *(phase++) = sum; 
             }
             phase  = (vector float*) rawdata_inA;
          }
          for (i=0; i<(SPE_FIXEDVECSIZE/16); i++) {
              #if 1
              /* wrap argument into -1.0..+1.0 using integer truncation (hmm this isn't quite proper.. +1->0 etc... , */
              //    phi = *(phase++);
              //    iphase = spu_convts(phi, 0);
              //    phi = spu_sub(phi, spu_convtf(iphase, 0));
              //    phi2 = spu_mul(phi, phi);
              UNROLL_BY_16( \
                  phi = *(phase++); \
                  phi2 = spu_mul(phi, phi); \
                  stmp = spu_madd(phi2, VEC_4f(0.080605269719834), VEC_4f(-0.006041231531886)); \
                  stmp = spu_madd(stmp, phi2, VEC_4f(-0.598397824003969)); \
                  stmp = spu_madd(stmp, phi2, VEC_4f(2.549926789216792)); \
                  stmp = spu_madd(stmp, phi2, VEC_4f(-5.167685041675656)); \
                  stmp = spu_madd(stmp, phi2, VEC_4f(3.141591733073902)); \
                  *(sinout++) = spu_mul(stmp, phi); \
                  ctmp = spu_madd(phi2, VEC_4f(-0.025391123907667), VEC_4f(0.001605367647616)); \
                  ctmp = spu_madd(ctmp, phi2, VEC_4f(0.235063383537898)); \
                  ctmp = spu_madd(ctmp, phi2, VEC_4f(-1.335174458310010)); \
                  ctmp = spu_madd(ctmp, phi2, VEC_4f(4.058698263116778)); \
                  ctmp = spu_madd(ctmp, phi2, VEC_4f(-4.934801388410942)); \
                  *(cosout++) = spu_madd(ctmp, phi2, VEC_4f(0.999999992289180)); \
              );
              #else
              UNROLL_BY_16( \
                  *(sinout++) = sin8_v(*phase); \
                  *(cosout++) = cos8_v(*phase); \
                  phase++; \
              );
              #endif
          }
        printf("------------------ PHASE ------------------\n");
        print_floatvecs((vector float*)rawdata_inA, SPE_FIXEDVECSIZE);
        printf("------------------ SINE ------------------\n"); // 
        print_floatvecs((vector float*)float_out, SPE_FIXEDVECSIZE);
        printf("------------------ COSINE ------------------\n");
        print_floatvecs((vector float*)xcorr_out, SPE_FIXEDVECSIZE);
        printf("------------------ MATLAB ------------------\n");
        print_floatvecs_matlab("phi", (vector float*)rawdata_inA, SPE_FIXEDVECSIZE);
        print_floatvecs_matlab("vsin", (vector float*)float_out, SPE_FIXEDVECSIZE);
        print_floatvecs_matlab("vcos", (vector float*)xcorr_out, SPE_FIXEDVECSIZE);
  }
  #endif

  /* stop the SPE-side measurement of the execution time */
  t_spu = t_start - spu_read_decrementer();

  printf("elements = %d\n", SPE_FIXEDVECSIZE);
  printf("t_spu = %d\n", t_spu);

  return 0;
}
