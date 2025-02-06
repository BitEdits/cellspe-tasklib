//===========================================================================
//=== P P U   C O D E                                                     ===
//===========================================================================

#include "spu/spunoptest.h"

#define NOP_REP_NUM 2500000             /* how many PPU->SPU->PPU roundtrips to measure */
extern spe_program_handle_t spunoptest; /* binary embed of SPU code */

int main(void) 
{
  int result;

  CellDSPTask_t NOPTask _QUAD_ALIGN;
  spe_spu_control_area_t *ps_area;
  speid_t speid;

  clock_t ctimes[8];
  struct tms eclock_tms_dmy;
  long clockspersec = sysconf(_SC_CLK_TCK);

  times(&eclock_tms_dmy);

  /* Initialize SPE thread */
  NOPTask.command = 0x00;
  speid = spe_create_thread(0, &spunoptest, 0, 0, -1, SPE_MAP_PS);
  ps_area = (spe_spu_control_area_t *)spe_get_ps_area(speid, SPE_CONTROL_AREA);

  /* Method I : busy polling */
  cout << "Cell NOP task repeated " << NOP_REP_NUM << " times, busy polling " << endl;
  ctimes[0] = times(&eclock_tms_dmy);
  for (unsigned long i=0; i<NOP_REP_NUM; i++) {
    spe_write_in_mbox(speid, (unsigned int)&NOPTask); // 32-bit lo
    #if SIM_64BIT
    spe_write_in_mbox(speid, (unsigned int)&NOPTask); // 32-bit hi (fake)
    #endif
    do {
       result = spe_read_out_mbox(speid);
    } while (result == -1);
  }
  ctimes[1] = times(&eclock_tms_dmy);
  cout << "   time delta = " << (ctimes[1]-ctimes[0]) <<
          ", clockspersec = " << clockspersec << endl;
  cout << endl;

  /* Method II : direct problem state busy polling */
  cout << "Cell NOP task repeated, direct problem state polling " << NOP_REP_NUM << " times" << endl;
  ctimes[0] = times(&eclock_tms_dmy);
  for (unsigned long i=0; i<NOP_REP_NUM; i++) {
    _spe_in_mbox_write(ps_area, (unsigned int)&NOPTask); // 32-bit lo
    #if SIM_64BIT 
    _spe_in_mbox_write(ps_area, (unsigned int)&NOPTask); // 32-bit hi (fake)
    #endif
    result = _spe_out_mbox_read(ps_area);
  }
  ctimes[1] = times(&eclock_tms_dmy);
  cout << "   time delta = " << (ctimes[1]-ctimes[0]) <<
          ", clockspersec = " << clockspersec << endl;
  cout << endl;

  return 0;
}
