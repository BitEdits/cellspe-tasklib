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

#define DEBUG 0

#define min(C,A,B) if (A<B) C=A; else C=B;

extern spe_program_handle_t minicorrelator_spu;

control_block cb[NUM_SPE_THREADS] _CLINE_ALIGN;
addr64        spe_ls_Ptrs[16]     _CLINE_ALIGN;
addr64        spe_sig_Ptrs[16]    _CLINE_ALIGN;

// -- HARD-CODED 6 STATION EXPERIMENT : data arrays for 128-byte aligned data
char*  raw_sources[NUM_SPE_THREADS];
float* autocorrelations[NUM_SPE_THREADS];
float* crosscorrelations[NUM_BASELINES];
float  fmultipliers[SPE_FIXEDBUFSIZE] _CLINE_ALIGN;

// -- Shared DMA-done tags, syncs that all other SPEs have DMA'ed data in. Sized to 128 byte.
unsigned short stationdata_locks[64]  _CLINE_ALIGN;

// ----------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------

void print_usage(char* text) {
  printf("USAGE: %s [options]\n", text);
  printf("       Valid Options include:\n");
  printf("         -h     : Output this usage help screen.\n");
  exit(1);
}

void print_error(char* text) {
  fprintf(stderr, "%s", text);
  exit(1);
}

double my_gettimeofday(void) {
  struct timeval time;
  gettimeofday( &time, (void *) 0);
  return (double)(time.tv_sec) + (double)time.tv_usec * 1.0e-6;
}

void evaluate_args(int argc, char *argv[]) {
  int i;
  for (i = 1; i < argc; i++) {
    if (*argv[i] == '-') {
      switch (*(argv[i]+1)) {
      default:
        print_usage(argv[0]);
        break;
      }
    } else print_usage(argv[0]);
  }
  return;
}

inline void print4f(float* f) {
    printf("%+#7.4e %+#7.4e %+#7.4e %+#7.4e ", f[0], f[1], f[2], f[3]);
}

// ----------------------------------------------------------------------------------------------------------------

int main(int argc, char *argv[]) {

  spe_gid_t gid;
  speid_t speid[NUM_SPE_THREADS];

  int i, j, status, flags;
  unsigned int spe_time[NUM_SPE_THREADS], spe_count[NUM_SPE_THREADS];
  double performance_all = 0.0, performance, t_start = 0.0, t_ppu = 0.0, data_amount = 0.0;

  /* parse command line arguments */
  evaluate_args(argc, argv);

  printf("\nUsing %u SPEs/stations.\nDataset: SPE buffers are %u kB. PPU has %u kB raw data for one SPE. Total %u kB for all stations.\n",
        NUM_SPE_THREADS, SPE_FIXEDBUFSIZE/1024, PPU_RAW_1SPEBYTES/1024, NUM_SPE_THREADS*PPU_RAW_1SPEBYTES/1024);

  printf("Dataset: computing %u complex FFTs of %u points, integrating over %u FFTs, %u baselines\n",
        PPU_RAW_TOTAL_FFTS, N_FFT, INTEGRATE_NUM_FFTS, NUM_BASELINES);
  printf("Resultset: %u kB in each integrated cross- or autocorrelation result chunk\n",
        PPU_CORRELBUF_BYTES/1024);

  /* malloc all buffers aligned to 128 Byte -- HARD-CODED 6 STATION EXPERIMENT */
  for (i = 0; i < NUM_SPE_THREADS; i++) {
      raw_sources[i] = _malloc_align(PPU_RAW_1SPEBYTES, 7);
      if (NULL == raw_sources[i]) { printf("Malloc of raw data input buf %d failed\n", i); return -1; }
      autocorrelations[i] = _malloc_align(PPU_CORRELBUF_BYTES, 7);
      if (NULL == autocorrelations[i]) { printf("Malloc of autocorrelation result buf %d failed\n", i); return -1; }
  }
  for (i = 0; i < NUM_BASELINES; i++) {
      crosscorrelations[i] = _malloc_align(PPU_CORRELBUF_BYTES, 7);
      if (NULL == crosscorrelations[i]) { printf("Malloc of cross-correlation result buf %d failed\n", i); return -1; }
  }
  for (i = 0; i < SPE_FIXEDBUFSIZE; i++) {
      fmultipliers[i] = 1.0 * ((i%1024)/1024.0);
  }
  for (i = 0; i< 64; i++) {
      stationdata_locks[i] = 0;
  }

  /* initialize raw inputs with test data */
  printf("Initializing raw data arrays... ");
  fflush(stdout);
  for (i = 0; i < NUM_SPE_THREADS; i++)
      for (j = 0; j < PPU_RAW_1SPEBYTES; j++)
          raw_sources[i][j] = (char)((17*i+3*j) % 256);

  /* put control block information together */
  for (i = 0; i < NUM_SPE_THREADS; i++) {
    cb[i].spe_num              = (int) i;
    cb[i].num_spes             = (int) NUM_SPE_THREADS;
    cb[i].ffts_total           = PPU_RAW_TOTAL_FFTS;
    cb[i].ffts_to_integrate    = INTEGRATE_NUM_FFTS;
    cb[i].spe_ls_listOnPPU     = (addr64) (unsigned long long)&spe_ls_Ptrs[0];
    cb[i].spe_sig_listOnPPU    = (addr64) (unsigned long long)&spe_sig_Ptrs[0];
    cb[i].rawdata_src          = (addr64) (unsigned long long)raw_sources[i];
    cb[i].autocorrelation_out  = (addr64) (unsigned long long)autocorrelations[i];
    cb[i].fmultipliers_src     = (addr64) (unsigned long long)fmultipliers;
    cb[i].syncline             = (addr64) (unsigned long long)stationdata_locks;
    for (j = 0; j < NUM_BASELINES; j++)
        cb[i].baselines_out[j] = (addr64) (unsigned long long)crosscorrelations[j];
  }

  if(DEBUG) printf("\n --debug checks: sizeof(cb)=%ld sizeof(spe_ls_Ptrs[])=%ld\n\n", sizeof(struct _control_block), sizeof(spe_ls_Ptrs));

  /* Create an SPE group */
  gid = spe_create_group (SCHED_OTHER, 0, 1);
  if ((gid == NULL) || (spe_group_max (gid) < 1)) exit(1);

  /* allocate the SPE tasks */
  for (i = 0; i < NUM_SPE_THREADS; i++) {

    flags = SPE_MAP_PS;
    if (i==0) { flags |= SPE_CFG_SIGNOTIFY1_OR; }
    speid[i] = spe_create_thread (gid, &minicorrelator_spu, (unsigned long long *) &cb[i], NULL, -1, flags);
    if (speid[i]== NULL) print_error("FAILED: spe_create_thread");

    spe_sig_Ptrs[i] = (addr64) (unsigned long long)spe_get_ps_area(speid[i], SPE_SIG_NOTIFY_1_AREA);
    spe_ls_Ptrs[i]  = (addr64) (unsigned long long)spe_get_ls(speid[i]);
    if(DEBUG) printf(" --debug: LS spe_ls_Ptrs[%d]=%p  SIG1 spe_sig_Ptrs[%d]=%p\n",
        i, (void*)spe_ls_Ptrs[i].ull, i, (void*)spe_sig_Ptrs[i].ull
     );
  }

  /* wait for a synchronisation signal of each SPE */
  for (i = 0; i < NUM_SPE_THREADS; i++) {
    while (!spe_stat_out_mbox(speid[i]));
    spe_read_out_mbox(speid[i]);
  }

  printf("Starting SPE calculations... ");
  fflush(stdout);

  /* start PPE-side measurement of the execution time */
  t_start = my_gettimeofday();

  /* send a start signal to each SPE */
  for (i = 0; i < NUM_SPE_THREADS; i++) spe_write_in_mbox(speid[i], 0);

  /* get the performance data of each SPE */
  for (i = 0; i < NUM_SPE_THREADS; i++) {
    while (!spe_stat_out_mbox(speid[i]));
    spe_time[i] = spe_read_out_mbox(speid[i]);
    while (!spe_stat_out_mbox(speid[i]));
    spe_count[i] = spe_read_out_mbox(speid[i]);
  }

  /* stop PPE-side measurement of execution time */
  t_ppu = my_gettimeofday() - t_start;

  /* wait until all SPE threads have finished */
  for (i = 0; i < NUM_SPE_THREADS; i++)
     spe_wait(speid[i], &status, 0);

  if (!WIFEXITED(status)) 
    print_error("FAILED: SPE abnormally terminated\n");
  else 
    printf("done!\n");

  /* print performance results */
  printf("\nPerformance results assuming a clock frequency of %d MHz and a timebase of %d:\n", __freq__, __timebase__);
  performance_all = 0;
  // per SPE, single DMA:
  data_amount = 8.0 * ((double)SPE_RAW_BYTES)/(1024*1024);
  for(i = 0; i < NUM_SPE_THREADS; i++) {
    performance = spe_count[i] * data_amount / ((double)spe_time[i]/__timebase__); 
    performance_all += performance/1024.0;
    printf("SPE %2d: %.2lf Gbit/sec, exectime %f sec, dmacount %d\n",
            i,
            performance/1024.0,
            (double)spe_time[i]/__timebase__, spe_count[i]);
  }

  printf("PPU-side gettimeofday() time delta: %.4lf sec\n", t_ppu);
  printf("PPU-side corresponding throughput : %.2lf Gbit/s (sum of SPE reported's)\n", performance_all );
  printf("PPU-side corresponding throughput : %.2lf Gbit/s (#SPEs*PPUdataPerSPU/t_ppu)\n",
             (8.0*NUM_SPE_THREADS*PPU_RAW_1SPEBYTES/(1024.0*1024.0*1024.0*t_ppu)) );
  printf("\n");

  printf("Aggregated performance for all %d SPEs: %.2lf Gbit/s (of %.2lf Gbit/s theoretical max).\n",
         NUM_SPE_THREADS,
         performance_all,
         8.0 * (NUM_SPE_THREADS * 51.2 + 26) // SPE's and 26 GB/s RAM
  );

  printf("Total data moved around: %.2lf MByte\n", data_amount / 8.0);
  printf("Final state lock vector: %u %u %u %u %u %u | %u %u | %u \n",
        stationdata_locks[0], stationdata_locks[1], stationdata_locks[2],
        stationdata_locks[3], stationdata_locks[4], stationdata_locks[5],
        stationdata_locks[6], stationdata_locks[7], stationdata_locks[8] );
  printf("\n");
  fflush(stdout);

  /* print some final baseline results */
  for (i = 0; i < NUM_SPE_THREADS; i++) {
     printf("Auto for SPE %2d = re: ", i);
     print4f(autocorrelations[i]+0); printf("...\t\tim: "); print4f(autocorrelations[i]+N_FFT);
     printf("...\n");
  }
  for (i = 0; i < NUM_BASELINES; i++) {
     printf("Cross for BL %2d = re: ", i);
     print4f(crosscorrelations[i]+0); printf("...\t\tim: "); print4f(crosscorrelations[i]+N_FFT);
     printf("...\n"); 
  }

  /* do rather unnecessary memory free()s: */
  for (i = 0; i < NUM_SPE_THREADS; i++) {
      _free_align(raw_sources[i]);
      _free_align(autocorrelations[i]);
  }
  for (i = 0; i < NUM_BASELINES; i++) {
      _free_align(crosscorrelations[i]);
  }
  return 0;
}
