//===========================================================================
//=== S P U   C O D E                                                     ===
//===========================================================================

#include "spunoptest.h"

CellDSPTask_t myTask _QUAD_ALIGN;   // currently fetched task

int main(unsigned long long spu_id, unsigned long long parm, unsigned long long envp) {
   unsigned long mbox_data;

   printf("SPU started.\n");

   while(1) {
      
      // wait for new command (64-bit address to task context inside message box)
      mbox_data = spu_read_in_mbox(); // 32-bit lo
      #if SIM_64BIT
      mbox_data = spu_read_in_mbox(); // 32-bit hi
      #endif

      // get task context
      mfc_getf((void *)(&myTask), mbox_data, sizeof(CellDSPTask_t), 31, 0, 0);
      mfc_write_tag_mask(1<<31);
      mfc_read_tag_status_all();

      // do nothing, just return "OK"
      spu_write_out_mbox(0x0000); 

   } //while(keeprunning)
   
   return 0;
}
