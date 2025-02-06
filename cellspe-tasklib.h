
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

#ifndef __CELLSPE_TASKLIB_H
#define __CELLSPE_TASKLIB_H

//---------------------------------------------------------------------------
// Common Defines
//---------------------------------------------------------------------------

/* SPE hogging limit */
#define CELLDSP_MAX_SPE_NR    4     // PS3 max 6, Cell max 8, QS20 max 24(?)

#define MALLOC_QUADALIGN     7       // as log2()
#define ALIGNMENT_QUADALIGN  128     // same as decimal

#define _QUAD_ALIGN  __attribute__ ((aligned(128))) // data alignment, compiler extension

//---------------------------------------------------------------------------
// Common Task Context Descriptor
//---------------------------------------------------------------------------

/* Flags and return codes */
#define CELLDSP_TASK_INCOMPLETE      0xFFFF  // still running or not started yet
#define CELLDSP_TASK_OK              0x0000
#define CELLDSP_TASK_GENERICERROR    0x0001
#define CELLDSP_TASK_UNKOWNCOMMAND   0x0002

/* Address struct for use between 32/64-bit PPU and 32-bit SPU */
typedef union
{
  unsigned long long ull;
  unsigned int ui[2];
  void *p;
} Cell_addr64;

/* Task sub-contexts with task-specific extra parameters */
typedef struct CellDSPParams_FFT_tt   { unsigned int nfft, log2_nfft; } CellDSPParams_FFT_t;
typedef struct CellDSPParams_Phase_tt { unsigned int argcount; float phase_start, phase_delta; } CellDSPParams_Phase_t;
typedef struct CellDSPParams_Difx_tt  { unsigned int nfft, log2_nfft, useUpperOrLower, twicenumchannels, numchannels; } CellDSPParams_Difx_t;
typedef union CellDSKSubcontext_tt {
   CellDSPParams_FFT_t   FFT;
   CellDSPParams_Phase_t Phase;
   CellDSPParams_Difx_t  Difx;
} CellDSPSubcontext_t;

/* Main task context */
#define CELLDSP_MAX_NUM_IOVECTORS    8
typedef struct CellDSPTask_tt {

   // command to execute
   unsigned int    command              _QUAD_ALIGN;
   // completion flag writte by SPE, -1 means not completed, other is the return code
   unsigned short  completed_flag;
   // for benchmarks, 32-bit SPU timer tick counter (for timebase: cat /proc/cpuinfo !)
   unsigned int    tick_count;

   // crude "variable arguments" list: input & output vectors and their lengths
   unsigned int    num_input_vectors;
   Cell_addr64     inputvector_ptrs     [CELLDSP_MAX_NUM_IOVECTORS];
   unsigned int    inputvector_lengths  [CELLDSP_MAX_NUM_IOVECTORS];  // length in bytes and not actual array elements
   unsigned int    num_output_vectors;
   Cell_addr64     outputvector_ptrs    [CELLDSP_MAX_NUM_IOVECTORS];
   unsigned int    outputvector_lengths [CELLDSP_MAX_NUM_IOVECTORS];

   // task sub-context (additional command specific params)
   CellDSPSubcontext_t subcontext;

   // for queueing up multiple tasks ("streaming"):
   Cell_addr64     next_task;           // .p=NULL if none, self-ref to NOP if "wait" / tail of stream

} CellDSPTask_t _QUAD_ALIGN;



//---------------------------------------------------------------------------
// PPU Specific Header
//---------------------------------------------------------------------------
#ifndef __SPU__
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ppu_intrinsics.h>
#include <cmath>
#include <cbe_mfc.h>
#include <libspe.h>
#include <libsync.h> // for SDK2.0, need to do 'make' and 'make install' in the 'locate libsync.h' SDK dir first

#include <pthread.h> // mutexes
#include <malloc.h>  // memalign() and free()
#include <unistd.h>  // usleep() for blocking getFreeSPE

#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>

#include <simdmath.h> // for some supposedly "fast" sin(),cos() etc

#include "types_generic.h" // data types for Generic and Cell

/* Profiling related */
#define SPU_TIMEBASE 79800000  // from 'cat /proc/cpuinfo | grep timebase' : 79800000 on PS3, 14318000 on QS20
#define SPU_CLOCK_HZ 3192e6    // this too from -"-
#define SPU_TICKS_TO_USEC(x)   (x*1e6/SPU_TIMEBASE)
#define SPU_TICKS_TO_CYCLES(x) (x*SPU_CLOCK_HZ/SPU_TIMEBASE)

/* SPE sharing and task control funcs */
int celldsp_start();
int celldsp_stop();
int celldsp_getFreeSPE();
int celldsp_setFreeSPE(int targetSPE);
int celldsp_startTask(int targetSPE, CellDSPTask_t * task);
int celldsp_waitTask(int targetSPE);
int celldsp_executeTask(CellDSPTask_t * task);
int celldsp_startContinuableTask(CellDSPTask_t * task, int* speID);
int celldsp_continueTask(CellDSPTask_t * task, int* speID);
int celldsp_stopContinuableTask(int* speID);
void celldsp_setBehaviour(unsigned int flags);
unsigned int celldsp_getBehaviour();

/* SPE bookkeeping struct for SPE sharing between threads and apps */
typedef struct CellDSP_tt {

  // Mapping access
  volatile int         mutexinit;        // 1 after shared mutex created
  volatile int         num_spes;         // how many SPEs handled by this shared memory struct
  pthread_mutex_t      tablemutex;       // to access the above tables
  pthread_mutexattr_t  tablemutexattr;
  // Free SPE signaling used in blocking behaviour
  pthread_cond_t       newfreespecond;   // used to signal first free'd SPE if there are blocking waiters
  volatile int         free_spe_count;
  volatile int         free_spe_waiters; // 1 if someone in blocking-wait for free SPE, 0 if not
  /// Thread barriers (can't use pthread_barrier since number of clients now known in advance)
  volatile int          clientcount;      // number of clients that called celldsp_bootThemAll()

  // Mapping
  volatile speid_t                 Allocations[CELLDSP_MAX_SPE_NR];     // maps SPE Core ID to thread ID
  volatile spe_gid_t               GroupMembership[CELLDSP_MAX_SPE_NR]; // group ID of each SPE thread
  volatile CellDSPTask_t*          LiveContexts[CELLDSP_MAX_SPE_NR];    // (can later be used for non-blocking DSP func implementations)
  volatile spe_spu_control_area_t* ControlAreas[CELLDSP_MAX_SPE_NR];    // direct problem state areas
  volatile unsigned int            ExtraFlags[CELLDSP_MAX_SPE_NR];      // some additional flags like CONTINUABLE etc

  // Behaviour
  volatile unsigned int            behaviour;

} CellDSP_t _QUAD_ALIGN;

/* CellDSP_t flag definitions */
// use for CellDSP_t.ExtraFlags[]:
#define CDSP_EF_RESERVED     1
#define CDSP_EF_CONTINUABLE  2
#define CDSP_EF_DEFAULT      0 // startup state flags
// use for CellDSP_t.behaviour:
#define CDSP_BEHAV_BLOCKING  1
#define CDSP_BEHAV_DEFAULT   0 // startup state flags

/******************************************************************************/
/* Helper Class: Memory aligned task context creation, modification, deletion */
/******************************************************************************/

class CellDSPTaskClass {
  protected:
    CellDSPTask_t task _QUAD_ALIGN;
  public:
     CellDSPTaskClass(unsigned int dspcommand) {
        task.command = dspcommand;
        task.num_input_vectors = 0; task.num_output_vectors = 0;
        task.next_task.p = NULL;
     }
     ~CellDSPTaskClass() {
     }
  public:
     void addInVec(void *src, size_t bytelength) {
        task.inputvector_ptrs[task.num_input_vectors].p = src;
        task.inputvector_lengths[task.num_input_vectors++] = bytelength;
     }
     void addOutVec(void *dest, size_t bytelength) {
        task.outputvector_ptrs[task.num_output_vectors].p = dest;
        task.outputvector_lengths[task.num_output_vectors++] = bytelength;
     }
     void setNext(CellDSPTask_t * nexttask) { 
        task.next_task.p = (void*)nexttask;
     }
     void setFFTparams(unsigned int fftlen) {
         task.subcontext.FFT.nfft = fftlen;
     }
     CellDSPTask_t * getTask() { 
        return &task;
     }
     CellDSPSubcontext_t * getSubcontext() {
        return &task.subcontext;
     }
     std::ostream& operator>>(std::ostream& out) {
        CellDSPTask_t * t = &task;
        out << "cmd=" << t->command 
            << " invecs=" << t->num_input_vectors << " outvecs=" << t->num_output_vectors;
        for(unsigned int i=0; i<t->num_input_vectors; i++) {
           out << " in[" << i << "]={ptr " << t->inputvector_ptrs[i].p << ", len " 
               << t->inputvector_lengths[i] << "} ";
        }
        for(unsigned int i=0; i<t->num_output_vectors; i++) {
           out << " out[" << i << "]={ptr " << t->outputvector_ptrs[i].p << ", len " 
               << t->outputvector_lengths[i] << "} ";
        }
        return out;
     }
     // TODO: add better support for queueable tasks
}; // class CellDSPTaskClass

#include "cellspe-tasklib-dsp.h"

#endif // !__SPU__


//---------------------------------------------------------------------------
// SPU Specific Header
//---------------------------------------------------------------------------
#ifdef __SPU__
#include <stdio.h>
#include <spu_intrinsics.h>
#include <spu_mfcio.h>
#include <vec_literal.h>
#include <spu_intrinsics.h>
#include <math.h>
#include <libmisc.h> // malloc_align(size_t length, int log2_of_align)
                     // for SDK2.0, need to do 'make' and 'make install' in the 'locate libmisc.h' SDK dir first
#include <stdlib.h>  // sbrk()
#include <profile.h>

#include <simdmath.h>
#include <libgmath.h>

#include <fft_1d.h>
#include <fft_1d_r2.h>

#define MAX_DMA_SIZE 16384
#define DMA_TAG_READABLES           1
#define DMA_TAG_MASK_READABLES      (1<<DMA_TAG_READABLES)
#define DMA_TAG_READABLES_DB        2
#define DMA_TAG_MASK_READABLES_DB   (1<<DMA_TAG_READABLES)
#define DMA_TAG_WRITEABLES          3
#define DMA_TAG_MASK_WRITEABLES     (1<<DMA_TAG_WRITEABLES)
#define DMA_TAG_WRITEABLES_DB       4
#define DMA_TAG_MASK_WRITEABLES_DB  (1<<DMA_TAG_WRITEABLES_DB)
#define DMA_TAG_TASKCONTEXTS        5
#define DMA_TAG_MASK_TASKCONTEXTS   (1<<DMA_TAG_TASKCONTEXTS)
#define DMA_TAG_MASK_ALL            -1

typedef struct { 
   float re; 
   float im; 
} fc32;
typedef fc32 cf32;

typedef struct { 
   double re; 
   double im; 
} fc64;
typedef fc64 cf64;

// currently processed task context
extern volatile CellDSPTask_t myTask _QUAD_ALIGN;

// function argument data DMA'ed onto SPE:
extern void *args[2*CELLDSP_MAX_NUM_IOVECTORS+1];  // data
extern unsigned int total_arg_arrays;              // how many args malloc()'ed

// function specific globals
extern unsigned int fft_N;         // FFT length N
extern unsigned int fft_N_log2;    // log2 of -"-
extern fc32* fft_W;                // FFT twiddle factors
extern unsigned char noptask_dummy[16384] _QUAD_ALIGN;

// some global helper funcs
extern void fetch_in_args();       // allocate -"- and load data from RAM into them
extern void fetch_in_out_args();   // allocate -"- and load data from RAM into them,
                                   // also fetch output data array from RAM
extern void free_args();           // deallocate -"-
extern void write_results();       // copy back results into RAM

inline void spu_dma_nowait(void *spu_addr_arg, volatile Cell_addr64 *ppu_addr_arg, const size_t sz_arg, int tag, unsigned int cmd);

#endif // __SPU__

#endif //  __CELLSPE_TASKLIB_H
