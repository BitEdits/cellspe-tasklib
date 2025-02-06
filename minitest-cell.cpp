
#define ARCH 3 // Cell

#include "cellspe-tasklib.h" /* IBM SDK 2.0 bug: need to include ppu_intrinsics.h before math.h or cmath */

#include "architecture.h"

#include <iostream>
#include <cmath>
using namespace std;

#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#define RUN_NOP           0
#define RUN_ARITHCHECK    0
#define RUN_IOTEST        0
#define RUN_IOBENCH       0
#define RUN_IOSTREAMBENCH 0
#define RUN_IODMASIMUL    0
#define RUN_MAC           0
#define RUN_MAC_BENCH     0
#define RUN_FFT_BASIC     0
#define RUN_FFT_BENCH     0
#define RUN_FFT_BTHREADED 0
#define RUN_STARTSTOP     0
#define RUN_SINCOS        1
#define RUN_DIFXPROCESS   0

#define NUM_PPU_THREADS   (CELLDSP_MAX_SPE_NR-2) // for threading tests and 'overbooking'
void* fft_test_thread(void * threadID);



int main(void)
{

  /***************** spe vs generics validation **************/
  //const int numinputbands = 7;
  //const int unpacksamples = 13;
  const int twicenumchannels = 1024;
  int ret = 0;
  // timing:
  clock_t ctimes[8];
  struct tms eclock_tms_dmy;
  long clockspersec = sysconf(_SC_CLK_TCK);
  times(&eclock_tms_dmy);
  // for FFT on SPE validation:
  vecFFTSpecC_f32 * c2cspec;
  int spe_id = 0;
  u8 fftbufsize = 0;
  int fftlen = 0;
  int fft_order = 0;
  int status = 0;
  // for threads:
  pthread_t mythreads[NUM_PPU_THREADS];

  // vectorAdd_cf32_I: performance and output SPE vs Generic
  cout << "Calling vectorAlloc_cf32 for 3 arrays" << endl;
  cf32 * testset1 = vectorAlloc_cf32(twicenumchannels*16);
  cf32 * testset2 = vectorAlloc_cf32(twicenumchannels*16);
  cf32 * resultset = vectorAlloc_cf32(twicenumchannels*16);
  cout << " got pointers: &s1=" << (void*)testset1 << " &s2=" << (void*)testset2 
       << " &results=" << (void*)resultset << endl; 
  if( ((unsigned long)testset1)%128 != 0 || ((unsigned long)testset2)%128 != 0 
       || ((unsigned long)resultset)%128 != 0) {
     cout << "   pointers not aligned to 128 bytes" << endl;
  } else {
     cout << "   pointers aligned to 128 bytes, OK" << endl;
  }

#if RUN_STARTSTOP
  cout << "Starting and killing SPE threads several times" << endl << flush;
  for (int i=0; i<10; i++) {
     cout << " start " << i;
     celldsp_start();
     cout << " stop " << i;
     celldsp_stop();
     cout << " | ";
  }
  cout << "done!" << endl << endl;
#endif


  cout << "Cell thread create" << endl;
  celldsp_start();


// ---------------------- NOP --------------------
// ---------------------- NOP --------------------
// ---------------------- NOP --------------------

#if RUN_NOP
  #define NOP_REP_NUM 2500000
  cout << "Cell NOP task repeated " << NOP_REP_NUM << " times" << endl;
  ctimes[0] = times(&eclock_tms_dmy);
  for (unsigned long i=0; i<NOP_REP_NUM; i++) {
     int nopres = 0;
     nopres = speNOP();
     // cout << "i=" << i << " nopres = " << speNOP() << endl;
  }
  ctimes[1] = times(&eclock_tms_dmy);
  cout << "   time delta = " << (ctimes[1]-ctimes[0]) << 
          ", clockspersec = " << clockspersec << endl;
  cout << endl;
#endif

// ---------------------- RUNALL/ARITHCHECK --------------------
// ---------------------- RUNALL/ARITHCHECK --------------------
// ---------------------- RUNALL/ARITHCHECK --------------------

#if RUN_ARITHCHECK
  #define RUNALL_REP_NUM 1
  cout << "Cell RUNALL task (some \"selected\" tasks ;-) repeated " << RUNALL_REP_NUM << " times" << endl;
  ctimes[0] = times(&eclock_tms_dmy);
  for (unsigned long i=0; i<RUNALL_REP_NUM; i++) {
     speRunAllCheck();
  }
  ctimes[1] = times(&eclock_tms_dmy);
  cout << "   time delta = " << (ctimes[1]-ctimes[0]) << 
          ", clockspersec = " << clockspersec << endl;
  cout << endl;
#endif

// ---------------------- MAC VALIDATION --------------------
// ---------------------- MAC VALIDATION --------------------
// ---------------------- MAC VALIDATION --------------------

  // MAC testing
#if RUN_MAC
  for (int turn=0; turn<4; turn++) {
     for (int i=0; i<twicenumchannels; i++) {
        testset1[i].re = 0.5*i;    
        testset1[i].im = -0.5*i;
        testset2[i].re = 0.15*i*i; 
        testset2[i].im = -0.25;
        resultset[i].re = 1.0; 
        resultset[i].im = 1.5;    
     }
     if (turn%2==0) {
        cout << "genericAddProduct_32fc()" << endl;
        genericAddProduct_32fc(testset1, testset2, resultset, twicenumchannels);
     } else {
        cout << "speAddProduct_32fc()" << endl;
        speAddProduct_32fc(testset1, testset2, resultset, twicenumchannels);
     }
     for (int i=0; i<8 && i<twicenumchannels; i++) {
        cout << "i=" << i << " : in1 [" << testset1[i].re << " " << testset1[i].im << "i ]";
        cout << " in2 [" << testset2[i].re << " i" << testset2[i].im << "]";
        cout << " out [" << resultset[i].re << " i" << resultset[i].im << "]" << endl;
     }
     cout << "..." << endl;
     for (int i=twicenumchannels-8; i<twicenumchannels; i++) {
        cout << "i=" << i << " : in1 [" << testset1[i].re << " " << testset1[i].im << "i ]";
        cout << " in2 [" << testset2[i].re << " i" << testset2[i].im << "]";
        cout << " out [" << resultset[i].re << " i" << resultset[i].im << "]" << endl;
     }
  }
  cout << endl;
#endif


// ---------------------- MAC THROUGHPUT --------------------
// ---------------------- MAC THROUGHPUT --------------------
// ---------------------- MAC THROUGHPUT --------------------

  // MAC benchmark testing
#if RUN_MAC_BENCH

  unsigned int mac_len = 1024; // has to be <=twicenumchannels
  unsigned long mac_num_iter = 1;// 1000000; // (1000000L*1024L/mac_len);
  cout << "MAC benchmark with complex MAC len=" << mac_len
       << ", num_iter=" << mac_num_iter << endl;
  for (int turn=0; turn<2; turn++) {

     for (int i=0; i<mac_len; i++) {
        testset1[i].re = 0.5*i;    
        testset1[i].im = -0.5*i;
        testset2[i].re = 0.15*i*i; 
        testset2[i].im = -0.25;
        resultset[i].re = 1.0; 
        resultset[i].im = 1.5;    
     }

     if (turn%2==0) {
        cout << "genericAddProduct_32fc()" << endl;
        ctimes[0] = times(&eclock_tms_dmy);
        for (unsigned long k=0; k<mac_num_iter; k++) {
//           genericAddProduct_32fc(testset1, testset2, resultset, mac_len);
        }
        ctimes[1] = times(&eclock_tms_dmy);
     } else {
        cout << "speAddProduct_32fc()" << endl;
        ctimes[0] = times(&eclock_tms_dmy);
        for (unsigned long k=0; k<mac_num_iter; k++) {
           speAddProduct_32fc(testset1, testset2, resultset, mac_len);
        }
        ctimes[1] = times(&eclock_tms_dmy);
     }
     cout << "  => CPU tick delta: " << (ctimes[1] - ctimes[0]) << ", clockspersec=" << double(clockspersec) << endl;
  }
  cout << endl;
#endif

  // Throughput testing 2
#if RUN_MAC_BENCH
  const unsigned int num_iter = 100000;
  cout << endl << "MAC throughput testing 2, " << num_iter << " iterations... " << endl;
  ctimes[0] = clock();
  for (unsigned int j=0; j<num_iter; j++) {
     genericAddProduct_32fc(testset1, testset2, resultset, twicenumchannels);
  }
  ctimes[1] = clock();
  for (unsigned int j=0; j<num_iter; j++) {
     speAddProduct_32fc(testset1, testset2, resultset, twicenumchannels);
  }
  ctimes[2] = clock();
  cout << "Generic cpu tick diff: " << (ctimes[1] - ctimes[0]) << " " << endl;
  cout << "SPE     cpu tick diff: " << (ctimes[2] - ctimes[1]) << " " << endl;
  cout << "CLOCKS_PER_SEC = " << CLOCKS_PER_SEC << endl;            
#endif


// ---------------------- SINE/COSINE VALIDATION --------------------
// ---------------------- SINE/COSINE VALIDATION --------------------
// ---------------------- SINE/COSINE VALIDATION --------------------

#if RUN_SINCOS
  cout << endl << "SINCOS validation test, argument count " << twicenumchannels << endl;
  float * phase_args = (float*)testset1;
  cout << " phase @ " << (void*)phase_args << " = {";
  for (int i=0; i<twicenumchannels; i++) {
     phase_args[i] =  (3.1415927/2) * i * 1.0/twicenumchannels;
     resultset[i].re = 0;
     resultset[i].im = 0;
     if(i<32) cout << " " << phase_args[i];
  }
  cout << " }" << endl;
  speSinCosAccurate(phase_args, resultset, twicenumchannels);
  cout << "speSinCosAccurate() => resultset(cos,sin) @ " << (void*)resultset << " = { ";
  for (int i=0; i<twicenumchannels && i<32; i++) {
        cout << " (" << resultset[i].re << "," << resultset[i].im << ") ";
  }
  cout << " } " << endl;

  float a = 0.1, b = 4 * 2*M_PI/twicenumchannels;
  cout << endl << "speSinCosGenerator(start=" << a << ", delta=" << b << ", iter=" << twicenumchannels/4 << ") ... " << endl;
  speSinCosGenerator(a, b, resultset, twicenumchannels/4);
  cout << "cosine = ";
  for (int i=0; i<(twicenumchannels/4); i++) {
        cout << resultset[i].re << " ";
  }
  cout << endl << "sine = ";
  for (int i=0; i<(twicenumchannels/4); i++) {
        cout << resultset[i].im << " ";
  }
  cout << endl << endl;
#endif

// ---------------------- C2C FFT VALIDATION --------------------
// ---------------------- C2C FFT VALIDATION --------------------
// ---------------------- C2C FFT VALIDATION --------------------

#if RUN_FFT_BASIC
  spe_id = -1;
  fftlen = 32;
  fft_order = 0; while (fftlen > (1<<fft_order)) { fft_order++; };
  cout << "FFT compare, order=" << fftlen << " or, log2(order)=" << fft_order << endl;
  status = speFFT_init(fftlen, &spe_id);
  cout << "speFFT_init(): " << status << " got spe_id=" << spe_id << endl;
  for (int turn=0; turn<2; turn++) {
     for (int i=0; i<fftlen; i++) {
        testset1[i].re = 0.5*sin(i*i/M_PI);    
        testset1[i].im = -0.5*i;
        resultset[i].re = 0.0;
        resultset[i].im = 0.0;
     }
     cout << "calling speFFT(), iteration " << turn << endl;
     speFFT_exec((void*)testset1, (void*)resultset, fftlen, &spe_id);
     cout << "testset = [ ";
     for (int i=0; i<fftlen && i<32; i++) {
        cout << "{" << testset1[i].re << "," << testset1[i].im << "i}, ";
     }
     cout << " ] " << endl << "result  = [ ";
     for (int i=0; i<fftlen && i<32; i++) {
        cout << "{" << resultset[i].re << "," << resultset[i].im << "i}, ";
     }
     cout << " ] " << endl;
  }
  cout << "speFFT_destroy()" << endl;
  speFFT_destroy(&spe_id);
  cout << endl;
#endif

// ---------------------- C2C FFT THROUGHPUT vs LENGTH, THREADED  --------------------
// ---------------------- C2C FFT THROUGHPUT vs LENGTH, THREADED  --------------------
// ---------------------- C2C FFT THROUGHPUT vs LENGTH, THREADED  --------------------

#if RUN_FFT_BTHREADED
  celldsp_setBehaviour(CDSP_BEHAV_BLOCKING); // enable new task wait if all SPEs busy
  for (int i=0; i<NUM_PPU_THREADS; i++) {
     cout << "FFT bench threaded, main(): launching thread " << i << endl << flush;
     pthread_create(&mythreads[i], NULL, fft_test_thread, (void*)i);
  }
  for (int i=0; i<NUM_PPU_THREADS; i++) {
     cout << "FFT bench threaded, main(): joining with thread " << i << endl << flush;
     int discard;
     pthread_join(mythreads[i], (void **)&discard);
  }
  cout << "FFT bench threaded, main(): done" << endl << flush;
#endif

// ---------------------- C2C FFT THROUGHPUT vs LENGTH, NON-THREADED  --------------------
// ---------------------- C2C FFT THROUGHPUT vs LENGTH, NON-THREADED  --------------------
// ---------------------- C2C FFT THROUGHPUT vs LENGTH, NON-THREADED  --------------------

#if RUN_FFT_BENCH
  spe_id = -1;
  for (fftlen=64; fftlen<=2*2048; fftlen*=2) {
  fft_order = 0; while (fftlen > (1<<fft_order)) { fft_order++; };
  const unsigned long fft_iterations = ((1250000L * 256L) / fftlen);

  cout << "FFT compare bench, order=" << fftlen << " or, log2(order)=" << fft_order << ", call iterations " << fft_iterations << endl;
  for (int i=0; i<fftlen; i++) {
     testset1[i].re = 0.5*sin(i*i/M_PI);    
     testset1[i].im = -0.5*i;
     resultset[i].re = 0.0;
     resultset[i].im = 0.0;
  }

  // -- CellDSP / SPE --
  status = speFFT_init(fftlen, &spe_id);
  cout << "speFFT_init(): " << status << "  got spe_id=" << spe_id << endl;
  cout << "speFFT() - " << (1250000 * 256 / fftlen) << " runs with " << fftlen << "-point FFT" << endl;

  ctimes[0] = times(&eclock_tms_dmy);
  for (unsigned long kk=0; kk<fft_iterations; kk++) {
    speFFT_exec((void*)testset1, (void*)resultset, fftlen, &spe_id);
  }
  ctimes[1] = times(&eclock_tms_dmy);

  cout << "CPU tick delta: " << (ctimes[1] - ctimes[0]) << ", clockspersec=" << double(clockspersec) << endl;
  cout << "  time per FFT [us]: "  << 1000000.0 * ((ctimes[1] - ctimes[0])/double(clockspersec))/double(fft_iterations) << endl;
  cout << "speFFT_destroy()" << endl;
  speFFT_destroy(&spe_id);
  cout << endl;

  cout << "speFFT_singleexec() - " << (1250000 * 256 / fftlen) << " runs with " << fftlen << "-point FFT" << endl;
  ctimes[0] = times(&eclock_tms_dmy);
  for (unsigned long kk=0; kk<fft_iterations; kk++) {
    speFFT_singleexec((void*)testset1, (void*)resultset, fftlen);
  }
  ctimes[1] = times(&eclock_tms_dmy);
  cout << "CPU tick delta: " << (ctimes[1] - ctimes[0]) << ", clockspersec=" << double(clockspersec) << endl;
  cout << "  time per FFT [us]: "  << 1000000.0 * ((ctimes[1] - ctimes[0])/double(clockspersec))/double(fft_iterations) << endl;
  cout << endl;

  }
  
  // -- FFTW -- 
  cout << "FFTW comparison" << endl;
  fftwf_plan plan;
  fftwf_complex * fftwIn = reinterpret_cast<fftwf_complex *>(testset1);
  fftwf_complex * fftwOut = reinterpret_cast<fftwf_complex *>(resultset);
  //fftwf_cell_set_nspe(3); // Cell FFTW 3.2alpha
  plan = fftwf_plan_dft_1d(fftlen, fftwIn, fftwOut, FFTW_FORWARD, FFTW_MEASURE);

  ctimes[0] = times(&eclock_tms_dmy);
  // for (unsigned long kk=0; kk<fft_iterations; kk++) {
  for (unsigned long kk=0; kk<10; kk++) {
     fftwf_execute_dft(plan, fftwIn, fftwOut);
  }
  ctimes[1] = times(&eclock_tms_dmy);

  cout << "CPU tick delta: " << (ctimes[1] - ctimes[0]) << ", clockspersec=" << double(clockspersec) << endl;
  fftwf_destroy_plan(plan);
#endif

// ---------------------- DIFX PROCESS() CHECKING  --------------------
// ---------------------- DIFX PROCESS() CHECKING  --------------------
// ---------------------- DIFX PROCESS() CHECKING  --------------------	

#if RUN_DIFXPROCESS
  cout << "Testing speDIFX_process() (remember to update minitest-cell.cpp after cellspe-tasklib-dsp.cpp mods!)" << endl << flush;
  void * sincosrotatedoutput = (void*)testset1;
  void * complexfracmult = NULL;
  void * fftd = (void*)resultset;
  void * fftoutputs = NULL;
  void * conjfftoutputs = NULL;
  void * autocorrelations = NULL;
  fftlen = 128;
  for (int i=0; i<64; i++) {
     int result = speDIFX_process((void*)sincosrotatedoutput,
                       (void*)complexfracmult, (void*)fftd /*(void*)fftoutputs[j]*/, (void*)conjfftoutputs,
                        (void*)autocorrelations,
                         fftlen, fftlen, 0);
     cout << "run " << i << " returned: " << result << endl << flush;
     for (int i=0; i<fftlen && i<32; i++) {
        cout << "{" << resultset[i].re << "," << resultset[i].im << "i}, ";
     }
     cout << endl;
  }
#endif


// ---------------------- IO DMA CHECK --------------------
// ---------------------- IO DMA CHECK --------------------
// ---------------------- IO DMA CHECK --------------------

#if RUN_IOTEST
  cout << "Cell IO task" << endl;
  speIOCheck();
  cout << endl;
#endif

// ---------------------- IO DMA THROUGHPUT 1 --------------------
// ---------------------- IO DMA THROUGHPUT 1 --------------------
// ---------------------- IO DMA THROUGHPUT 1 --------------------
	 
  // Throughput testing 1 -- get all SPEs busy with DMAing into same RAM location (SPU->PPU DMA only)
  //                         Inside the SPE code this does several (>1000) writes back to the RAM buffer per one exec.
#if RUN_IOBENCH
  u8 * biggerarray = vectorAlloc_u8(16384);
  #define MYSPES_NR    CELLDSP_MAX_SPE_NR
  cout << "Throughput test, writing a batch of 2*16 in-line DMAs 65535 times into main memory, double buffered" << endl;
  cout << "DMA is SPE->PPU, using " << MYSPES_NR << " SPE's" << endl;
  for (int i=1; i<=64; i*=2) {
     cout << "speIOBench with " << (i*16384/64) << " byte chunks" << endl;
     ctimes[0] = times(&eclock_tms_dmy);
     speIOBench((void*)(biggerarray),  (i*16384/64), MYSPES_NR);
     ctimes[1] = times(&eclock_tms_dmy);
     cout << "CPU tick delta: " << (ctimes[1] - ctimes[0]) << " or " 
          << double(clockspersec) * (MYSPES_NR * 2.0*16.0*65535.0) * (double(i)*16384.0/64.0)/((double)ctimes[1] - (double)ctimes[0]) << " bytes per sec" <<  endl;
  }
  cout << endl;
#endif


  // Throughput testing 3 -- get all SPEs busy with dumping sequential data blocks into same 
  //                         RAM location (SPU->PPU DMA only)
  //                         Inside the SPE code this does several (>1000) writes back to the RAM buffer per one exec.
#if RUN_IOSTREAMBENCH
  u8 * biggerarray = vectorAlloc_u8(16384);
  const int stream_spe_count = 2;
  #define NUM_STREAM_TASKS  16
  cout << "Streamed throughput test, writing a batch of 16 DMAs 1000 times into main memory, " 
       << "using " << NUM_STREAM_TASKS << " streamed tasks per call " << endl;
  cout << "DMA is SPE->PPU, using " << stream_spe_count << " SPEs" << endl;
  for (int i=2; i<=64; i*=2) {
     cout << "speIOBenchStreamed with " << (i*16384/64) << " byte chunks" << endl;
     ctimes[0] = times(&eclock_tms_dmy);
     speIOBenchStreamed((void*)(biggerarray),  (i*16384/64), stream_spe_count, NUM_STREAM_TASKS);
     ctimes[1] = times(&eclock_tms_dmy);
     cout << "CPU tick delta: " << (ctimes[1] - ctimes[0]) << " or " 
          << NUM_STREAM_TASKS * stream_spe_count * double(clockspersec) * (16*1000000) * (i*16384/64)/((double)ctimes[1] - (double)ctimes[0]) << " bytes per sec" <<  endl;
  }
  cout << endl;
#endif

  // Throughput testing 4 -- one SPE runs single DMA in ("bit-packed data") and dumps out
  //                         "result data" using a DMA list.
#if RUN_IODMASIMUL
  cout << "Running DMA simulation with 'bit-packed' input data and 'processed result data' out."  << endl;
  ctimes[0] = times(&eclock_tms_dmy);
  speIOSimulation();
  ctimes[1] = times(&eclock_tms_dmy);
  cout << "CPU tick delta: " << (ctimes[1] - ctimes[0]) << endl;
#endif

  cout << "celldsp_stop() ..." << flush;
  celldsp_stop();
  cout << "done." << endl;

  pthread_exit(NULL);
  return ret;
}



void* fft_test_thread(void* threadID)
{
  vecFFTSpecC_f32 * c2cspec;
  int fftlen = 0, fft_order = 0, status = 0;
  const int twicenumchannels = 1024;
  clock_t ctimes[8];
  struct tms eclock_tms_dmy;
  long clockspersec = sysconf(_SC_CLK_TCK);
  times(&eclock_tms_dmy);

  cf32 * testset1 = vectorAlloc_cf32(twicenumchannels*16);
  cf32 * testset2 = vectorAlloc_cf32(twicenumchannels*16);
  cf32 * resultset = vectorAlloc_cf32(twicenumchannels*16);

  int spe_id = -1;
  int thread_id = (int)threadID;

  for (fftlen=64; fftlen<=2*2048; fftlen*=2) {
    fft_order = 0; while (fftlen > (1<<fft_order)) { fft_order++; };
    const unsigned long fft_iterations = ((1250000L * 256L) / fftlen);

    cout << "pputhread" << thread_id << " - FFT compare bench, order=" << fftlen << " or, log2(order)=" << fft_order 
         << ", call iterations " << fft_iterations << endl << flush;
    for (int i=0; i<fftlen; i++) {
       testset1[i].re = 0.5*sin(i*i/M_PI);    
       testset1[i].im = -0.5*i;
       resultset[i].re = 0.0;
       resultset[i].im = 0.0;
    }

    status = speFFT_init(fftlen, &spe_id);
    cout << "pputhread" << thread_id << " - speFFT_init(): " << status << "  got spe_id=" << spe_id << endl
         << "pputhread" << thread_id << " - speFFT() - " << (1250000 * 256 / fftlen) << " runs with " << fftlen << "-point FFT" << endl << flush;

    ctimes[0] = times(&eclock_tms_dmy);
    for (unsigned long kk=0; kk<fft_iterations; kk++) {
      speFFT_exec((void*)testset1, (void*)resultset, fftlen, &spe_id);
    }
    ctimes[1] = times(&eclock_tms_dmy);

    cout << "pputhread" << thread_id << " - CPU tick delta: " << (ctimes[1] - ctimes[0]) << ", clockspersec=" << double(clockspersec) 
         << " ~ time per FFT [us]: "  << 1000000.0 * ((ctimes[1] - ctimes[0])/double(clockspersec))/double(fft_iterations) << endl << flush;
    speFFT_destroy(&spe_id);
    cout << "pputhread" << thread_id << " - speFFT_destroy()" << endl << flush;

    cout << "pputhread" << thread_id << " - speFFT_singleexec() - " << (1250000 * 256 / fftlen) << " runs with " << fftlen << "-point FFT" << endl;
    ctimes[0] = times(&eclock_tms_dmy);
    for (unsigned long kk=0; kk<fft_iterations; kk++) {
      speFFT_singleexec((void*)testset1, (void*)resultset, fftlen);
    }
    ctimes[1] = times(&eclock_tms_dmy);
    cout << "pputhread" << thread_id << " - CPU tick delta: " << (ctimes[1] - ctimes[0]) << ", clockspersec=" << double(clockspersec) 
         << " ~ time per FFT [us]: "  << 1000000.0 * ((ctimes[1] - ctimes[0])/double(clockspersec))/double(fft_iterations) << endl << flush;

  }
 
  delete testset1; delete testset2; delete resultset;
  cout << "pputhread" << thread_id << " - done, exiting." << endl << endl << flush;
  pthread_exit(NULL);
}
