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

// hard-coded size of data blocks that are processed
#define SPE_FIXEDBUFSIZE     16384                 // 16kB
#define SPE_FIXEDFLOATSIZE   (SPE_FIXEDBUFSIZE/4)  // same as # of floats
#define SPE_FIXEDVECSIZE     (SPE_FIXEDBUFSIZE/16) // same as # of vectors

// some unrolling helpers
#define UNROLL_BY_2(x)  { x }{ x }
#define UNROLL_BY_4(x)  { x }{ x }{ x }{ x }
#define UNROLL_BY_8(x)  UNROLL_BY_4(x) UNROLL_BY_4(x)
#define UNROLL_BY_16(x) UNROLL_BY_8(x) UNROLL_BY_8(x)

// constant float scalar to vector splat
#define VEC_4f(x)       (vector float){x, x, x, x}


// ------------------------------------------------------------------------
// Fast dual sine/cosine calculation
// Argument must be -PI..+PI but normalized to -1.0..+1.0
//
// Single SPE performance with sin8_v(), cos8_v():
//   3.96 Gbit/sec, exectime 4.037809 sec, dmacount 131072 phase arg blocks * 16kB
//
// Single SPE performance with code below:
//  53.43 Gbit/sec, exectime 0.299470 sec, dmacount 131072
//
// ------------------------------------------------------------------------

         /*  Coefficients were generated with Matlab:

             phi=linspace(-pi, pi, 64000);
             phinorm = phi ./ pi;
             trueSine = sin(phi);
             trueCosine = cos(phi);
             coeffsS = polyfit(phinorm, trueSine, 12);
             coeffsC = polyfit(phinorm, trueCosine, 12);
             figure(1), plot(phinorm, polyval(coeffsS,phinorm)-trueSine, 'g-');
             figure(2), plot(phinorm, polyval(coeffsC,phinorm)-trueCosine, 'g-');

          */
   // for(#blocks) {
          vector float *phase  = (vector float*) rawdata_inA;
          vector float *sinout = (vector float*) float_out;
          vector float *cosout = (vector float*) xcorr_out;
          vector float stmp, ctmp, phi2, phi;
          vector signed int iphase;
          for (i=0; i<(SPE_FIXEDVECSIZE/16); i++) {
              #if 1
              /* wrap argument into -1.0..+1.0 using integer truncation (hmm this isn't quite proper.. +1->0 etc... , */
              //    phi = *(phase++);
              //    iphase = spu_convts(phi, 0);
              //    phi = spu_sub(phi, spu_convtf(iphase, 0));
              //    phi2 = spu_mul(phi, phi);
              UNROLL_BY_16( \
                  phi = *(phase++); \
                  iphase = spu_convts(phi, 0); \
                  phi = spu_sub(phi, spu_convtf(iphase, 0)); \
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
              // for speed comparison only
              UNROLL_BY_16( \
                  *(sinout++) = sin8_v(*phase); \
                  *(cosout++) = cos8_v(*phase); \
                  phase++; \
              );
              #endif
          }

          totaldmascount ++;
     // }



// ------------------------------------------------------------------------
// Fast dual sine/cosine calculation for a continuous phase increment
//
// Single SPE performance with code below:
//  59.15 Gbit/sec, exectime 8.453110 sec, dmacount 4096000 
//  (and actually a bit faster than that, removed two multiply-adds from 
//   the loop...)
//
// ------------------------------------------------------------------------

   // for(#blocks) {
          // Start phase x, constant increment y, iterate:
          //    sin(x+y) = sin x cos y + cos x sin y
          //    cos(x+y) = cos x cos y - sin x sin y

          float angleinc = 0.23, startangle = 0.1;

          vector float argvec, oldsin, oldcos, anginc4;

          argvec  = spu_madd( spu_splats(angleinc), 
                              ((vector float) { 0.0, 1.0, 2.0, 3.0 }),
                              /* + */ spu_splats(startangle) );
          anginc4 = spu_mul(spu_splats(angleinc), ((vector float) { 4.0, 4.0, 4.0, 4.0 }) );

          vector float * sin_out = (vector float*) rawdata_inA;
          vector float * cos_out = (vector float*) rawdata_inB;

          vector float sin4Y     = sin14_v(anginc4);
          vector float cos4Y     = cos14_v(anginc4);
          vector float sinvec    = sin14_v(argvec);
          vector float cosvec    = cos14_v(argvec);

          for (i=0; i<SPE_FIXEDVECSIZE/8; i++) {
            UNROLL_BY_8( \
              oldsin = sinvec; \
              *(sin_out++) = oldsin; \
              sinvec = spu_madd(oldsin, cos4Y, spu_mul(cosvec, sin4Y)); \
              oldcos = cosvec; \
              *(cos_out++) = oldcos; \
              cosvec = spu_nmsub(oldsin, sin4Y, spu_mul(oldcos, cos4Y)); \
            );
          }

          totaldmascount ++;
     // }



// ------------------------------------------------------------------------
// Fast complex multiply-accumulate (MAC, FMA)
//
// Single SPE performance:
//  3 * 95.41 Gbit/sec, exectime 0.041925 sec, dmacount 32768 (*3 bufs)
//
// ------------------------------------------------------------------------

   // for(#blocks) {
          // --- MAC
          // 16kB data block: 4k floats: 2k complex: 2 1024-point FFTs: 2 x 512 vectors
          // two 16kB blocks A,B for 2x1024 samples, one 16kB for 2x1024 accu
          {
              vector float *in1re = (vector float*) rawdata_inA;
              vector float *in1im = in1re + 256;
              vector float *in2re = (vector float*) rawdata_inB;
              vector float *in2im = in2re + 256;
              vector float *accre = (vector float*) xcorr_out;
              vector float *accim = accre + 256;
              vector float acctmp1, acctmp2;
              for(i=0; i<512/8; i++) {
                 UNROLL_BY_8( \
                    acctmp1 = spu_msub(*in1im, *in2im, spu_mul(*in1re, *in2re)); \
                    *accre  = spu_add(*accre, acctmp1); \
                    acctmp2 = spu_madd(*in1im, *in2re, spu_mul(*in1re, *in2im)); \
                    *accim  = spu_add(*accim, acctmp2); \
                    accre++; accim++; in1im++; in2im++; in1re++; in2re++;
                 );
              }
              in1re += 256; in2re += 256; in1im += 256; in2im += 256;
              accre += 256; accim += 256;
              for(i=0; i<512/8; i++) {
                 UNROLL_BY_8( \
                    acctmp1 = spu_msub(*in1im, *in2im, spu_mul(*in1re, *in2re)); \
                    *accre  = spu_add(*accre, acctmp1); \
                    acctmp2 = spu_madd(*in1im, *in2re, spu_mul(*in1re, *in2im)); \
                    *accim  = spu_add(*accim, acctmp2); \
                    accre++; accim++; in1im++; in2im++; in1re++; in2re++;
                 );
              }
          }

          totaldmascount ++;
     // }



// ------------------------------------------------------------------------
// Fast unpacking of 2-bit raw data samples into float vectors
//
// Output is a 16kB array of vector floats, input 1 byte results in
// one vector float (16 byte) thus the input array size is 1kB
//
// Single SPE performance:
//  83.67 Gbit/sec, exectime 0.812698 sec, dmacount 557056
//
// ------------------------------------------------------------------------

          // 4 x 2-bit samples (1 byte) to 4 floats lookup, initialize!
          vector float byte_to_vecfloat[256]; 

          unsigned char c;
          vector unsigned char * v_raw = (vector unsigned char*) rawdata_inA;
          vector float *         out   = (vector float *)float_out;

          #define RCONV(x)  c=spu_extract(*v_raw,x); *(out++)=byte_to_vecfloat[c];
          for (i=0; i<SPE_FIXEDVECSIZE/16; i++) {
             RCONV(0);  RCONV(1);  RCONV(2);  RCONV(3);
             RCONV(4);  RCONV(5);  RCONV(6);  RCONV(7);
             RCONV(8);  RCONV(9);  RCONV(10); RCONV(11);
             RCONV(12); RCONV(13); RCONV(14); RCONV(15);
             v_raw ++;
          }
          totaldmascount ++;
