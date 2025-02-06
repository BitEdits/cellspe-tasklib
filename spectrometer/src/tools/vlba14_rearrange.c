//--------------------------------------------------------------------------------------------------------------
//
// Rearranges VLBA 1:4 fan-out 2-bit 4-channel data to [Ch1_sign,Ch1_mag][s,m][s,m][s,m]
//
// Compiling:
//   gcc -g -O3 -Wall -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 vlba14_rearrange.c -o vlba14_rearrange
//
//--------------------------------------------------------------------------------------------------------------

//#define DEBUG 1 // define to get bit pattern printouts

#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64
#define _GNU_SOURCE 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <malloc.h>

#include <sys/types.h>

#define BUF_LENGTH 8192

typedef u_int32_t vsi_t;

#define bit0 1
#define bit1 2
#define bit2 4
#define bit3 8
#define bit4 16
#define bit5 32
#define bit6 64
#define bit7 128

void bin_prnt_byte(int x);

int main(int argc, char** argv) {

   FILE      *fd_in;
   FILE      *fd_out;

   vsi_t     raw_word;
   vsi_t    *raw_buf;

   size_t     raws_in_buf;
   size_t     raws_written;
   off64_t    bytecount;
   int        i, j, iter;

   if (argc < 3) {
       fprintf(stderr, "\n"
                       "Rearrange VLBA 1:4 fan-out 2-bit 4-channel data to [Ch1_sign,Ch1_mag][s,m][s,m][s,m]\n\n"
                       "vlba14_rearrange <input file> <output file>\n\n");
       return -1;
   }

   /* alloc */
   raw_buf = (vsi_t*)memalign(128, sizeof(vsi_t) * BUF_LENGTH);
   if (NULL == raw_buf) {
       fprintf(stderr, "Could not allocate conversion buffer!\n");
       return -1;
   }

   /* open the files */
   fprintf(stderr, "Reading '%s' and converting into '%s'...\n", argv[1], argv[2]);
   fd_in  = fopen64(argv[1], "rb");
   fd_out = fopen64(argv[2], "wb");
   if (NULL == fd_in || NULL == fd_out) {
       fprintf(stderr, "Couldn't open file(s)!\n");
       return -2;
   }

   /* process */
   raws_in_buf = BUF_LENGTH;
   bytecount   = 0;
   iter        = 0;
   while (raws_in_buf > 0) {

       /* read more data */
       raws_in_buf  = fread((void*)raw_buf, sizeof(vsi_t), BUF_LENGTH, fd_in);
       bytecount   += (off64_t)(sizeof(vsi_t)) * (off64_t)raws_in_buf;

       for (i=0; i<raws_in_buf; i++) {
           // Notation: <s/m>_ch1..4_time1..4
           // VLBA 1:4 1-bit 4 channel
           //                  time 1          time 2          time 3           time 4        VSI
           //   byte 3      s_1_1   s_3_1   s_1_2   s_3_2   s_1_3    s_3_3   s_1_4    s_3_4 : bit 0..7
           //   byte 2      m_1_1   m_3_1   m_1_2   m_3_2   m_1_3    m_3_3   m_1_4    m_3_4 : bit 8..15
           //   byte 1      s_2_1   s_4_1   s_2_2   s_4_2   s_2_3    s_4_3   s_2_4    s_4_4 : bit 16..23
           //   byte 0      m_2_1   m_4_1   m_2_2   m_4_2   m_2_3    m_4_3   m_2_4    m_4_4 : bit 24..31
           //           bit   0       1       2       3       4        5       6        7
           //
           // The output is essentially a bit transpose. Output should be:
           //   {s_1_1 m_1_1 s_2_1 m_2_1 s_3_1 m_3_1 s_4_1 m_4_1} {*_2} {*_3} {*_4}  ; {time1},{time2},{time3},{time4}
           u_char  timeX_4chs[4];
           u_char* raw_in = (u_char*)(&raw_buf[i]);

           #ifdef DEBUG
           printf("----------------\n");
           printf("in    : ");   bin_prnt_byte(raw_in[3]);
           printf("\n        "); bin_prnt_byte(raw_in[2]);
           printf("\n        "); bin_prnt_byte(raw_in[1]);
           printf("\n        "); bin_prnt_byte(raw_in[0]); printf("\n");
           #endif

           for (j=0; j<4; j++) {
               /* extract 4 channels at one timepoint */
               timeX_4chs[j] =   ((raw_in[3] & bit0) << (7-0)) | ((raw_in[2] & bit0) << (6-0)) /* ch1 s,m */
                               | ((raw_in[1] & bit0) << (5-0)) | ((raw_in[0] & bit0) << (4-0)) /* ch2 s,m */
                               | ((raw_in[3] & bit1) << (3-1)) | ((raw_in[2] & bit1) << (2-1)) /* ch3 s,m */
                               | ((raw_in[1] & bit1) << (1-1)) | ((raw_in[0] & bit1) >> (1  )) /* ch4 s,m */
                               ;

               #ifdef DEBUG
               printf("time %d: ", j); bin_prnt_byte(timeX_4chs[j]); printf("\n");
               #endif


               /* shift input data to next timepoint */
               raw_in[3] >>= 2; raw_in[2] >>= 2; raw_in[1] >>= 2; raw_in[0] >>= 2;
           }

           /* arrange back into input buffer */
           raw_in[0] = timeX_4chs[0];
           raw_in[1] = timeX_4chs[1];
           raw_in[2] = timeX_4chs[2];
           raw_in[3] = timeX_4chs[3];
       }

       /* write converted data */
       raws_written = fwrite((void*)raw_buf, sizeof(raw_word), raws_in_buf, fd_out);
       if (raws_written != raws_in_buf) {
           fprintf(stderr, "short write!\n");
       }

       /* progress */
       if (((++iter) % 512) == 0) {
           fprintf(stderr, "MByte converted: %6.2f\n", bytecount / (1024.0*1024.0));
       }
   }

   /* done */
   free(raw_buf);
   fclose(fd_in);
   fclose(fd_out);
   return 0;
}

// debug function
void bin_prnt_byte(int x) {
   int n;
   char bitstr[10];
   for(n=0; n<8; n++) {
      if((x & 0x80) !=0) { bitstr[n] = '1'; } else { bitstr[n] = '0'; }
      x = x<<1;
   }
   bitstr[8] = ' ';
   bitstr[9] = '\0';
   fwrite(bitstr, 1, 10, stdout);
}
