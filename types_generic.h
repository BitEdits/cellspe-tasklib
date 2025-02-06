
/***************************************************************************
 *   Copyright (C) 2007 by Jan Wagner                                      *
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

#ifndef _TYPESGENERIC_H
#define _TYPESGENERIC_H

#include <fftw3.h>       // use ffwtf_*** for single precision to match f32, fc32 datatypes
                         // NOTE: FFTW has to be './configure --enable-float'ed before make and install

enum IppHintAlgorithm {
  ippAlgHintNone = 0,
  ippAlgHintFast = 1,
  ippAlgHintAccurate = 2,
};

//start with the types
typedef unsigned char 	u8;
typedef unsigned short 	u16;
typedef short 		    s16;
typedef int 		    s32;
typedef float 		    f32;
typedef double 		    f64;

typedef struct { short re; short im; } sc16;
typedef sc16 		cs16;

typedef struct { f32 re; f32 im; } fc32;
typedef fc32 		cf32;

typedef struct { f64 re; f64 im; } fc64;
typedef fc64		cf64;


typedef f32 fftwf_real; // "FFTW 2 had an fftw_real typedef that was an alias for double (in double precision). In FFTW 3 you should just use double (or whatever precision you are employing)."

typedef struct { 
  fftwf_plan      plan; 
  fftwf_complex * in;
  fftwf_complex * out;
  int             order;
} vecFFTSpecC_f32;

typedef struct {
  fftwf_plan      plan; 
  fftwf_real    * in;
  fftwf_complex * out;
  int             order;
} vecFFTSpecR_f32;

// dummy definitions
#define vecFFT_NoReNorm          0

//then the constant values
#define vecNoErr                 0
#define vecTrue                  1
#define TWO_PI                   6.28318530716
#define MAX_S16                  ((1<<15) - 1)
#define MAX_U16                  ((1<<16) - 1)
#define vecHintAlg               IppHintAlgorithm 
#define vecAlgHintFast           ippAlgHintFast

#endif
