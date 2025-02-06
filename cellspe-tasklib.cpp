
/*************************************************************************** 
 *  Copyright (C) 2007 by Jan Wagner                                       *
 *                                                                         *
 *  Digital Signal Processing Primitives Library                           *
 *  for IBM Cell Broadband Engine                                          *
 *                                                                         *
 *  Cell SPU Processor Version                                             *
 *                                                                         *
 *  License: GNU LGPL                                                      *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
 ***************************************************************************/


#define DO_CHEAT        0   // 0: normal mode, 1: internal benchmarking mode
#define XPU_64BIT       0   // 0: PPU code is 32-bit, 1: PPU code is 64-bit
#define VERBOSE_MALLOC  1   // 0: don't check for NULL, 1: print error if NULL
#define TASK_POLLING    1   // 0: use SPE outbound mailbox to check return value (fast, but no SPE task queue from threads possible)
                            // 1: SPE changes integer in task context which PPU polls (slow, but threaded SPE task queue usage works)
#define RET_SPU_CYCLES  1   // 0: don't write SPU cycle count of task back to task context
                            // 1: return number of SPU cycles of executed task inside task context

#define ALIGN_CHECKING  0   // 1: for debug, check some DMA alignments and complain if faulty...

#include "cellspe-tasklib.h"
#ifdef __SPU__
  #include "cellspe-tasklib-dsp.cpp" // yuck...
#else
  #include "cellspe-tasklib-dsp.h"
#endif



//===========================================================================
//=== C O M M O N                                                         ===
//===========================================================================

inline void CHECK_POINTER_ALIGN(void * src, void * dst) {
   #if ALIGN_CHECKING
   if(  ( ((unsigned int)src) & (ALIGNMENT_QUADALIGN-1) )
      ||( ((unsigned int)dst) & (ALIGNMENT_QUADALIGN-1) )
      )
  { printf("unaligned: CHECK_POINTER_ALIGN(src=%p dst=%p)\n", src, dst); }
   #endif
}



//===========================================================================
//=== P P U   C O D E                                                     ===
//===========================================================================
#ifndef __SPU__

// TODO: this gets really messy... if multiple instances of
//       host application, how to communicate SPE allocations
//       between the app instances?

using namespace std;

//---------------------------------------------------------------------------------
// Declares
//---------------------------------------------------------------------------------

// Binary-embed the DSP Primitives Library:
extern spe_program_handle_t cellspe_tasklib;

// SPE mapping and properties
CellDSP_t *CellDSP = NULL;
CellDSP_t mainCellDSP; // (currently not used..)

//---------------------------------------------------------------------------------
// Global funcs instantiated in this .cpp
//---------------------------------------------------------------------------------
int           celldsp_start();
int           celldsp_stop();
void          celldsp_setBehaviour(unsigned int flags);
unsigned int  celldsp_getBehaviour();

int           celldsp_getFreeSPE();
int           celldsp_setFreeSPE(int targetSPE);

int           celldsp_startTask(int targetSPE, CellDSPTask_t * task);
int           celldsp_waitTask(int targetSPE);

int           celldsp_executeTask(CellDSPTask_t * task);

int           celldsp_startContinuableTask(CellDSPTask_t * task, int* speID);
int           celldsp_continueTask(CellDSPTask_t * task, int* speID);
int           celldsp_stopContinuableTask(int* speID);


//---------------------------------------------------------------------------------
// Local helpers
//---------------------------------------------------------------------------------
inline void PRINT_SPU_TASKSTATS(CellDSPTask_t * finishedtask) {
   unsigned int ticks = finishedtask->tick_count;
   cout << "SPU decremeter: " << ticks << " ticks ("
        << SPU_TICKS_TO_USEC(ticks) << " usec, "
        << "or " << SPU_TICKS_TO_CYCLES(ticks) << " cycles)"
        << endl << flush;
   return;
}

//---------------------------------------------------------------------------------
// celldsp_start() : initialize all DSP threads so they get uploaded to SPE's
//---------------------------------------------------------------------------------
int celldsp_start() {

   // already running: re-use mutexes, just register as new client
   if(CellDSP != NULL) {
      pthread_mutex_lock(&CellDSP->tablemutex);
      CellDSP->clientcount++;
      pthread_mutex_unlock(&CellDSP->tablemutex);
      cerr << "celldsp_start: you're client " << CellDSP->clientcount << ", will re-use SPEs " << endl;
      return 0;
   }

   CellDSP = new CellDSP_t;

   // need to create mutexes
   pthread_mutexattr_init(&CellDSP->tablemutexattr);
   pthread_mutexattr_settype(&CellDSP->tablemutexattr, PTHREAD_MUTEX_DEFAULT);
   pthread_mutexattr_setpshared(&CellDSP->tablemutexattr, PTHREAD_PROCESS_SHARED);
   pthread_mutex_init(&CellDSP->tablemutex, &CellDSP->tablemutexattr);
   if(pthread_mutex_lock(&CellDSP->tablemutex) != 0) {
      cerr << "celldsp_start: could not lock newly created mutex!!" << endl;
   }
   CellDSP->mutexinit = 1;
   // also create wait condition
   pthread_cond_init(&CellDSP->newfreespecond, NULL);

   // start SPE's
   CellDSP->num_spes = CELLDSP_MAX_SPE_NR;
   CellDSP->free_spe_count = CELLDSP_MAX_SPE_NR;
   CellDSP->free_spe_waiters = 0;
   for (int i=0; i<CELLDSP_MAX_SPE_NR; i++) {
      speid_t speid = spe_create_thread(0, &cellspe_tasklib, NULL, NULL, -1, SPE_MAP_PS);
      CellDSP->Allocations[i]     = speid;
      CellDSP->GroupMembership[i] = spe_get_group(speid);
      CellDSP->ControlAreas[i]    = (spe_spu_control_area_t *)spe_get_ps_area(speid, SPE_CONTROL_AREA);
      CellDSP->LiveContexts[i]    = NULL;
      CellDSP->ExtraFlags[i]      = CDSP_EF_DEFAULT;
   }

   // set some flags
   CellDSP->behaviour = CDSP_BEHAV_DEFAULT;
   CellDSP->clientcount = 1;

   if(pthread_mutex_unlock(&CellDSP->tablemutex) != 0) {
      cerr << "celldsp_start: mutex unlock failed" << endl;
   }
   return 0;
}

//---------------------------------------------------------------------------------
// celldsp_stop() : kill all SPE threads, shut down task library
//---------------------------------------------------------------------------------
int celldsp_stop() {
   int status = 0;
   CellDSPTaskClass grimreaper(CELLDSP_TASK_SIGTERM);

   if(pthread_mutex_lock(&CellDSP->tablemutex) != 0) {
      cerr << "celldsp_stop: could not lock mutex" << endl;
      return -1;
  }

   // only allow the last caller to actually kill the threads
   CellDSP->clientcount--;
   if(CellDSP->clientcount>0) {
      cerr << "celldsp_stop: other clients still using this tasklib instance, deferring SPE thread kill" << endl;
      pthread_mutex_unlock(&CellDSP->tablemutex);
      return 0;
   }
   // send a terminate task to all SPE's
   for (int i=0; i<CELLDSP_MAX_SPE_NR; i++) {
      CellDSPTask_t* harvest = grimreaper.getTask();
      CellDSP->LiveContexts[i] = harvest;
      #if XPU_64BIT
      unsigned long long addr64 = (unsigned long long)harvest;
      _spe_in_mbox_write(CellDSP->ControlAreas[i], (unsigned int)(addr64 >> 32));
      _spe_in_mbox_write(CellDSP->ControlAreas[i], (unsigned int)(addr64 & 0xFFFFFFFF));
      #else
      _spe_in_mbox_write(CellDSP->ControlAreas[i], (unsigned int)harvest);
      #endif
   }
   for (int i=0; i<CELLDSP_MAX_SPE_NR; i++) {
      // wait for reply in messagebox, then wait until thread exits
      status = _spe_out_mbox_read(CellDSP->ControlAreas[i]);
      spe_wait(CellDSP->Allocations[i], &status, 0);
   }
   pthread_mutex_unlock(&CellDSP->tablemutex);
   pthread_mutex_destroy(&CellDSP->tablemutex);
   pthread_mutexattr_destroy(&CellDSP->tablemutexattr);
   pthread_cond_destroy(&CellDSP->newfreespecond);
   // release the shared mem segment
   delete CellDSP;
   CellDSP = NULL;
   return 0;
}

//---------------------------------------------------------------------------------
// celldsp_setBehaviour() : update task lib behaviour in different situations
//---------------------------------------------------------------------------------
void celldsp_setBehaviour(unsigned int flags) {
   if(CellDSP==NULL) return;
   pthread_mutex_lock(&CellDSP->tablemutex);
   CellDSP->behaviour = flags;
   pthread_mutex_unlock(&CellDSP->tablemutex);
}

//---------------------------------------------------------------------------------
// celldsp_getBehaviour() : returns the current task lib behaviour flags
//---------------------------------------------------------------------------------
unsigned int celldsp_getBehaviour() {
   if(CellDSP==NULL) return 0;
   return CellDSP->behaviour;
}

//---------------------------------------------------------------------------------
// celldsp_getFreeSPE() : returns the number of the first free SPE in the list,
//                        on error returns -1
//---------------------------------------------------------------------------------
int celldsp_getFreeSPE() {
   int ret = -1;
   if(pthread_mutex_lock(&CellDSP->tablemutex) != 0) {
      cerr << "celldsp_getFreeSPE: could not lock mutex" << endl;
      return -1;
   }
getfreespe_retry:
   // Locate a free SPE
   for (int i=0; i<CELLDSP_MAX_SPE_NR; i++) {
      if (!(CellDSP->ExtraFlags[i] & CDSP_EF_RESERVED)) {
         CellDSP->ExtraFlags[i] |= CDSP_EF_RESERVED;
         CellDSP->free_spe_count--;
         ret = i;
         break;
      } 
   }
   // If none are free, and blocking behaviour is enabled, wait until an SPE is free
   if((CellDSP->behaviour & CDSP_BEHAV_BLOCKING) && (ret==-1)) {
      CellDSP->free_spe_waiters = 1;
      pthread_cond_wait(&CellDSP->newfreespecond, &CellDSP->tablemutex);
      CellDSP->free_spe_waiters = 0;
      goto getfreespe_retry; // need at least one of these per C++ program :-)
   }
   if(pthread_mutex_unlock(&CellDSP->tablemutex) != 0) {
      cerr << "celldsp_getFreeSPE: mutex unlock failed" << endl;
   }
   return ret;
}

//---------------------------------------------------------------------------------
// celldsp_setFreeSPE() : frees up the given SPE (be careful from where you call 
//                        this..)
//---------------------------------------------------------------------------------
int celldsp_setFreeSPE(int targetSPE) {
   if(pthread_mutex_lock(&CellDSP->tablemutex) != 0) {
      cerr << "celldsp_setFreeSPE: could not lock mutex" << endl;
      return -1;
   }
   // Release the reservation (don't do 'double-free')
   if(CellDSP->ExtraFlags[targetSPE] & CDSP_EF_RESERVED) {
      CellDSP->ExtraFlags[targetSPE] = CDSP_EF_DEFAULT;
      CellDSP->LiveContexts[targetSPE] = NULL;
      CellDSP->free_spe_count++;
      // Notify possible blocking-waiters of first free'd SPE
      if((CellDSP->behaviour & CDSP_BEHAV_BLOCKING) 
         && (1==CellDSP->free_spe_waiters)
         && (1==CellDSP->free_spe_count)) {
         pthread_cond_broadcast(&CellDSP->newfreespecond);
      }
   } else {
      cerr << "warning: celldsp_setFreeSPE(SPE " << targetSPE << ") was already free" << endl;
   }
   if(pthread_mutex_unlock(&CellDSP->tablemutex) != 0) {
      cerr << "celldsp_setFreeSPE: mutex unlock failed" << endl;
   }
   return 0;
}

//---------------------------------------------------------------------------------
// celldsp_startTask() : sends a task to the specified SPE for execution
//                       returns 0 on success
//---------------------------------------------------------------------------------
int celldsp_startTask(int targetSPE, CellDSPTask_t * task) {
   int ret = -1;
   if(pthread_mutex_lock(&CellDSP->tablemutex) != 0) {
      cerr << "celldsp_startTask: could not lock mutex" << endl;
      return -1;
   }
   // check that SPE was pre-reserved, send task
   if (CellDSP->ExtraFlags[targetSPE] & CDSP_EF_RESERVED) {
      CellDSP->LiveContexts[targetSPE] = task;
      task->completed_flag = CELLDSP_TASK_INCOMPLETE;
      #if XPU_64BIT
      unsigned long long addr64 = (unsigned long long) task;
      _spe_in_mbox_write(CellDSP->ControlAreas[targetSPE], (unsigned int)(addr64 >> 32));
      _spe_in_mbox_write(CellDSP->ControlAreas[targetSPE], (unsigned int)(addr64 & 0xFFFFFFFF));
      #else
      _spe_in_mbox_write(CellDSP->ControlAreas[targetSPE], (unsigned int)task);
      #endif
      ret = 0;
   } else {
      cerr << "celldsp_startTask: CellDSP->LiveContexts["<<targetSPE<<"]=" << CellDSP->LiveContexts[targetSPE] << " != 1" << endl;
   }
   if(pthread_mutex_unlock(&CellDSP->tablemutex) != 0) {
      cerr << "celldsp_startTask: mutex unlock failed" << endl;
   }
   return ret;
}

//---------------------------------------------------------------------------------
// celldsp_waitTask() : waits until the specified SPE has finished the task,
//                      then passes back the SPE task return code, doesn't free
//                      up the SPE
//---------------------------------------------------------------------------------

int celldsp_waitTask(int targetSPE) {
   // Busy wait for completion and read result 
   int result = 0;
   #if TASK_POLLING
   while(CellDSP->LiveContexts[targetSPE]->completed_flag==CELLDSP_TASK_INCOMPLETE);
   result = CellDSP->LiveContexts[targetSPE]->completed_flag;
   #else
   result = _spe_out_mbox_read(CellDSP->ControlAreas[targetSPE]);
   #endif
   return result;
}

//---------------------------------------------------------------------------------
// celldsp_executeTask() : reserves SPE and tries to execute the task, waits until 
//                         completion, frees up the SPE and passes back the SPE task 
//                         return code
//---------------------------------------------------------------------------------
int celldsp_executeTask(CellDSPTask_t * task) {
   // get an SPE
   int freespe = celldsp_getFreeSPE();
   if (freespe == -1) {
      cerr << "celldsp_executeTask(cmd=" << task->command << ") : could not get any free SPE" << endl;
      return !vecNoErr;
   }
   // execute, then blocking-wait until finished
   if (celldsp_startTask(freespe, task) == -1) {
      cerr << "celldsp_executeTask(cmd=" << task->command << ") : task start failed" << endl;
      return !vecNoErr;
   }
   int result = celldsp_waitTask(freespe);
   if (result != CELLDSP_TASK_OK) {
      cerr << "celldsp_executeTask : task returned error code " << result << endl;
      celldsp_setFreeSPE(freespe);
      return !vecNoErr;
   }
   celldsp_setFreeSPE(freespe);
   return vecNoErr;
}

//---------------------------------------------------------------------------------
// celldsp_startContinueableTask() : gets and SPE ID into *speID and then tries to 
//                                   execute the first task, returns SPE return code
//---------------------------------------------------------------------------------
int celldsp_startContinuableTask(CellDSPTask_t * task, int* speID) {
   // get an SPE
   *speID = celldsp_getFreeSPE();
   if (*speID == -1) {
      cerr << "celldsp_startContinuableTask(cmd=" << task->command << ") : could not get any free SPE" << endl;
      return !vecNoErr;
   }
   // execute, then blocking-wait until finished
   if (celldsp_startTask(*speID, task) == -1) {
      cerr << "celldsp_startContinuableTask(cmd=" << task->command << ") : task start failed" << endl;
      return !vecNoErr;
   }
   int result = celldsp_waitTask(*speID);
   if (result != CELLDSP_TASK_OK) {
      cerr << "celldsp_startContinuableTask : task returned error code " << result << endl;
      celldsp_setFreeSPE(*speID);
      return !vecNoErr;
   }
   return vecNoErr;
}


//---------------------------------------------------------------------------------
// celldsp_continueTask() : tries to execute the given task on an earlier used SPE, 
//                          waits until completion then passes back the SPE task 
//                          return code
//---------------------------------------------------------------------------------
int celldsp_continueTask(CellDSPTask_t * task, int* speID) {
   // execute, then blocking-wait until finished
   if (celldsp_startTask(*speID, task) == -1) {
      cerr << "celldsp_continueTask(cmd=" << task->command << ") : task start failed" << endl;
      return !vecNoErr;
   }
   int result = celldsp_waitTask(*speID);
   if (result != CELLDSP_TASK_OK) {
      cerr << "celldsp_continueTask : task returned error code " << result << endl;
      return !vecNoErr;
   }
   return vecNoErr;
}

//---------------------------------------------------------------------------------
// celldsp_stopContinuableTask() : frees up the SPE
//---------------------------------------------------------------------------------
int celldsp_stopContinuableTask(int* speID) {
   celldsp_setFreeSPE(*speID);
   *speID = -1;
   return vecNoErr;
}


/**************************************************************************/




//===========================================================================
//=== S P U   C O D E                                                     ===
//===========================================================================
#else

// extern globals instantiated here
volatile CellDSPTask_t myTask _QUAD_ALIGN;   // currently fetched task
void *args[2*CELLDSP_MAX_NUM_IOVECTORS+1];   // command input data arrays
unsigned int total_arg_arrays;               // size of *[] -"-

// some things copied from Cell FFTW 3.2alpha (GNU GPL):
inline void spu_dma_nowait(void *spu_addr_arg, volatile Cell_addr64 *ppu_addr_arg, const size_t sz_arg, int tag, unsigned int cmd);

//---------------------------------------------------------------------------------
// main() : main task handler, waits for new tasks in the SPE mailbox,
//          then tries to execute them (can be streamed in "DMA scatter-gather" style), 
//          and afterwards sends back a return code into the PPU-side mailbox
//---------------------------------------------------------------------------------


int main(unsigned long long spu_id, unsigned long long parm, unsigned long long envp) {
   unsigned char keeprunning     = 1;  // a CELLDSP_TASK_SIGTERM task sets this to 0
   unsigned long long  mbox_data = 0;  // 64-bit ptr into PPU memory
   unsigned int  mbox_returncode = 0;  // return code
   unsigned int  i;

   /* Init */
   fft_W = NULL;
   fft_N = 0; fft_N_log2 = 0;
   parm = 0; keeprunning = 1;
   // sbrk(40*1024); // reserve 40kB heap
  
   printf("CellDSPTasklib now running on SPU with id %lld\n", spu_id);

   /* While the world is turning */
   while(keeprunning) {
       
      /* Wait for new task context address in the mailbox */
      mbox_data = spu_read_in_mbox();
      #if XPU_64BIT
         // get upper 32-bit and assemble into full 64-bit address
         mbox_data = (mbox_data << 32) | spu_read_in_mbox();
      #endif

      do {

        /* Fetch task context */
        mfc_getf((void *)(&myTask), mbox_data, sizeof(CellDSPTask_t), DMA_TAG_TASKCONTEXTS, 0, 0);
        mfc_write_tag_mask(DMA_TAG_MASK_TASKCONTEXTS);
        mfc_read_tag_status_all();
 
        /* Initialize task timing if desired for benchmarks */
        #if RET_SPU_CYCLES
        spu_write_decrementer(0xffffffff);
        #endif

        /* Handle the task */
        switch(myTask.command) {
            case CELLDSP_TASK_DIFXPROCESS:
                mbox_returncode = exec_difxprocess();
                break;
            case CELLDSP_TASK_FFT_EXEC_DYN:
                fetch_in_args();
                mbox_returncode = exec_FFT_dyn();
                free_args();
                break;
            case CELLDSP_TASK_FFT_EXEC:
                fetch_in_args();
                mbox_returncode = exec_FFT();
                free_args();
                break;
            case CELLDSP_TASK_MAC:
                fetch_in_out_args();
                mbox_returncode = do_complex_mac();
                free_args();
                break;
            case CELLDSP_TASK_FFT_INIT:
                mbox_returncode = init_FFT();
                break;
            case CELLDSP_TASK_FFT_FREE:
                mbox_returncode = free_FFT();
                break;
            case CELLDSP_TASK_MULC:
            case CELLDSP_TASK_ADDC:
                mbox_returncode = CELLDSP_TASK_OK;
                break;                        
            case CELLDSP_TASK_SINCOS_GENERTR:
                fetch_in_args();
                mbox_returncode = exec_sincos_generator();
                free_args();
                break;
            case CELLDSP_TASK_SINCOS:
                fetch_in_args();
                // TODO: on first and every N:th call update global coefficient tables
                mbox_returncode = exec_sincos();
                free_args();
                break;
            case CELLDSP_TASK_SINCOS_NORM:
                fetch_in_args();
                // TODO: on first and every N:th call update global coefficient tables
                mbox_returncode = exec_sincosnormalized();
                free_args();
                break;
            case CELLDSP_TASK_UNPACK:
                fetch_in_args();
                mbox_returncode = exec_unpack();
                free_args();
                break;            
            case CELLDSP_TASK_IOTEST:
                mbox_returncode = do_in_out_test();
                break;
            case CELLDSP_TASK_IOBENCH:
                fetch_in_out_args();
                mbox_returncode = do_in_out_bench();
                free_args();
                break;
            case CELLDSP_TASK_NOP:
                #if 0
                fetch_in_args();
                free_args();
                #endif
                #if 0
                   #if XPU_64BIT
                   mfc_getf((void*)&noptask_dummy[0], myTask.inputvector_ptrs[0].ull, 16384, 30, 0, 0);
                   #else
                   mfc_getf((void*)&noptask_dummy[0], myTask.inputvector_ptrs[0].ui[0], 16384, 30, 0, 0);
                   #endif
                   mfc_write_tag_mask(1<<30);
                   mfc_read_tag_status_all();
                #endif                
                mbox_returncode = CELLDSP_TASK_OK;
                break;
            case CELLDSP_TASK_RUNALL:
//                myTask.num_input_vectors = 4;
//                myTask.num_output_vectors = 0;
                // -- just to get the code into .s file for timing analysis
                exec_vecrepack_versio2();
//                exec_unpack_improved();
//                // --
//                for (int ii=0; ii<4; ii++) {
//                   args[ii] = malloc_align(16384, MALLOC_QUADALIGN);
//                   myTask.inputvector_lengths[ii] = 512;
//                   myTask.outputvector_lengths[ii] = 512;
//                }
//                mbox_returncode = do_runall_test();
//                for (int ii=0; ii<4; ii++) {
//                   free_align(args[ii]);
//                }
                break;
            case CELLDSP_TASK_DMASIMULATION:
                mbox_returncode = do_io_simulation();
                break;
            case CELLDSP_TASK_SIGTERM:
                keeprunning = 0;
                spu_write_out_mbox(CELLDSP_TASK_OK); 
                return 0;
                break;                            
            default:
                mbox_returncode = CELLDSP_TASK_UNKOWNCOMMAND;
                break;
        }

        /* Write some return values back to PPU side struct, if needed */
        #if RET_SPU_CYCLES
        myTask.tick_count = -spu_read_decrementer();
        #endif
        #if TASK_POLLING
        myTask.completed_flag = mbox_returncode;
        #endif
        #if (TASK_POLLING || RET_SPU_CYCLES)
        CHECK_POINTER_ALIGN((void *)(&myTask), (void*)mbox_data);
        mfc_putf((void *)(&myTask), mbox_data, 16, DMA_TAG_TASKCONTEXTS, 0, 0); // note: 16 bytes should fit flag and cyclecount in struct
        #endif

        /* Check for the next streamed task */
        #if XPU_64BIT
        mbox_data = myTask.next_task.ull;
        #else
        mbox_data = myTask.next_task.ui[0];
        #endif
        // mbox_data = 0; // set to disable streaming feature

      } while(mbox_data != 0);

      
      /* Return the result code in the outbound mailbox */
      #if !TASK_POLLING
      spu_write_out_mbox(mbox_returncode); 
      // TODO: right now only one PPU thread per SPE task at a time is possible,
      // maybe somehow make this so that different PPU threads can fill up
      // SPE inbox with their tasks, and result outbox gets mapped back to
      // correct PPU initiator thread..?
      #endif

   } //while(keeprunning)
   
   return 0;
}



//---------------------------------------------------------------------------------
// spu_dma_nowait() : read or write data using DMA, splits it into 16kB chunks,
//                    start addresses must have proper 128-byte alignment
// [function copied and improved from Cell FFTW 3.2alpha GNU GPL]
//---------------------------------------------------------------------------------
inline void spu_dma_nowait(void *spu_addr_arg,
                           volatile Cell_addr64 *ppu_addr_arg,
                           const size_t sz_arg, int tag, unsigned int cmd)
{
   char *spu_addr              = (char *)spu_addr_arg;
   #if XPU_64BIT
   unsigned long long ppu_addr = ppu_addr_arg->ull;
   #else
   unsigned long long ppu_addr = ppu_addr_arg->ui[0];
   #endif
   size_t sz                   = sz_arg;
   size_t chunk                = MAX_DMA_SIZE;
   CHECK_POINTER_ALIGN(spu_addr_arg, (void*)ppu_addr);
   while (sz > 0) {
     if (sz < MAX_DMA_SIZE) {
       chunk = sz;
     }
     spu_mfcdma64(spu_addr, mfc_ea2h(ppu_addr), mfc_ea2l(ppu_addr),
                  chunk, tag, cmd);
     sz -= chunk; ppu_addr += chunk; spu_addr += chunk;
   }
}


//---------------------------------------------------------------------------------
// fetch_in_args(): load all input data from main memory to Localstore
//---------------------------------------------------------------------------------
void fetch_in_args() {
   int index = 0;
   int length;

   /* How many entries into pointer array */
   total_arg_arrays = myTask.num_input_vectors + myTask.num_output_vectors;
   /* Queue up input data DMA's, fence DMA */
   for (unsigned int i=0; i<myTask.num_input_vectors; i++) {
      length = myTask.inputvector_lengths[i];
      args[index] = malloc_align(length, MALLOC_QUADALIGN);
      #if VERBOSE_MALLOC
      if (args[index]==NULL) {
        printf("fetch_in_args: input vector %d length %d malloc failed\n", i, length);
      }
      #endif
      spu_dma_nowait(args[index], &myTask.inputvector_ptrs[i], length, DMA_TAG_READABLES, MFC_GETF_CMD);
      index++;
   }
   /* Output data arrays, just allocate enough memory */
   for (unsigned int i=0; i<myTask.num_output_vectors; i++) {
      length = myTask.outputvector_lengths[i];
      args[index] = malloc_align(length, MALLOC_QUADALIGN);
      #if VERBOSE_MALLOC
      if (args[index]==NULL) {
        printf("fetch_in_args: output vector %d length %d malloc failed\n", i, length);
      }
      #endif
      index++;
   }

   /* Stall until DMA completion */
   mfc_write_tag_mask(DMA_TAG_MASK_READABLES);
   mfc_read_tag_status_all();
   return;
}

//---------------------------------------------------------------------------------
// fetch_in_out_args(): load all input and output data from main memory to Localstore
//---------------------------------------------------------------------------------
void fetch_in_out_args() {
   int index = 0;
   int length;

   /* How many entries into pointer array */
   total_arg_arrays = myTask.num_input_vectors + myTask.num_output_vectors;

   /* Queue up input data DMA's, fence DMA */
   for (unsigned int i=0; i<myTask.num_input_vectors; i++) {
      length = myTask.inputvector_lengths[i];
      args[index] = malloc_align(length, MALLOC_QUADALIGN);
      #if VERBOSE_MALLOC
      if (args[index]==NULL) {
        printf("fetch_in_out_args: input vector %d length %d malloc failed\n", i, length);
      }
      #endif
      spu_dma_nowait(args[index], &myTask.inputvector_ptrs[i], length, DMA_TAG_READABLES, MFC_GETF_CMD);
      index++;
   }
   /* "Output" data required as input for e.g. MAC accumulation, so get these too! */
   for (unsigned int i=0; i<myTask.num_output_vectors; i++) {
      length = myTask.outputvector_lengths[i];
      args[index] = malloc_align(length, MALLOC_QUADALIGN);
      #if VERBOSE_MALLOC
      if (args[index]==NULL) {
        printf("fetch_in_out_args: output vector %d length %d malloc failed\n", i, length);
      }
      #endif
      spu_dma_nowait(args[index], &myTask.outputvector_ptrs[i], length, DMA_TAG_WRITEABLES, MFC_GETF_CMD);
      index++;
   }
   
   /* Stall until DMA completion */
   mfc_write_tag_mask(DMA_TAG_MASK_READABLES | DMA_TAG_MASK_WRITEABLES);
   mfc_read_tag_status_all();
   return;
}

//---------------------------------------------------------------------------------
// write_results() : write back all calculation results into main memory
//---------------------------------------------------------------------------------
void write_results() {
   int index;

   /* "Output" data required as input for e.g. MAC accumulation, so get these too! */
   index = myTask.num_input_vectors; // output data arrays start after input arrays
   for (unsigned int i=0; i<myTask.num_output_vectors; i++) {
      unsigned int length = myTask.outputvector_lengths[i];
      spu_dma_nowait(args[index], &myTask.outputvector_ptrs[i], length, DMA_TAG_WRITEABLES, MFC_PUTF_CMD);
      index++;
   }
   
   /* Stall until DMA completion */
   mfc_write_tag_mask(DMA_TAG_MASK_WRITEABLES);
   mfc_read_tag_status_all();
   return;
}

//---------------------------------------------------------------------------------
// free_args() : deallocate the arrays in LocalStore
//---------------------------------------------------------------------------------
void free_args() {
   
   /* Free up earlier allocations */
   for (unsigned int i=0; i<total_arg_arrays; i++) {
      if(args[i]!=NULL) free_align(args[i]);
      args[i] = NULL;
   }

   total_arg_arrays = 0;
   return;
}

#endif //  __CELLSPE_TASKBLIB_H
