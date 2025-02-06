#ifndef __SPUNOPTEST_H
#define __SPUNOPTEST_H

#define SIM_64BIT 0  /* 1 to simulate 64-bit address transfer in mbox */

#ifdef __SPU__

  /* SPU headers */

  #include <spu_intrinsics.h>
  #include <spu_mfcio.h>
  #include <libmisc.h>
  #include <stdio.h>

#else

  /* PPU headers */

  #include <libspe.h>
  #include <iostream>
  #include <time.h>
  #include <sys/times.h>
  #include <unistd.h>
  #include <cbe_mfc.h>
  using namespace std;

#endif

/* Common */

#define _QUAD_ALIGN  __attribute__ ((aligned(128)))
typedef struct CellDSPTask_tt {
   // command to execute
   unsigned int  command       _QUAD_ALIGN;
   unsigned int  filler[31];
} CellDSPTask_t                _QUAD_ALIGN;

#endif
