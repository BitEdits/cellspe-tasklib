/***************************************************************************
 *   Copyright (C) 2009 by Jan Wagner                                      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mark5access.h>

#define MAX_LIST_LEN (512*1024)

int main(int argc, char* argv[])
{
    struct mark5_stream *ms;
    float **unpacked;
    int rate, channels, fanout, bitreso, format;
    char *longlist, *longlist_alloc;
    int combinations_tested = 0;
    const int rates[5] = { 64, 128, 256, 512, 1024 }; // only the most common ones...

    if (argc <= 1) {
        printf("\n  Usage: determine_tape_format <filename>\n\n"
               "  Tries to open the file using MkIV and VLBA formats with a\n"
               "  zillion of possible fan-out, data rate, channel count and\n"
               "  bit resolution combinations. The first successful combination\n"
               "  is reported.\n\n");
        return -1;
    }

    longlist_alloc = malloc(MAX_LIST_LEN);
    longlist = longlist_alloc;
    longlist[0] = '\0';

    unpacked = (float**)malloc(64 * sizeof(float*));
    for (channels = 0; channels < 64; channels++) {
        unpacked[channels] = (float*)malloc(1024*sizeof(float));
    }

    for (format = 0; format <= 1; format++) {
        for (fanout = 1; fanout <= 8; fanout *= 2) {
            for (rate = 0; rate < sizeof(rates)/sizeof(int); rate++) {
                for (channels = 1; channels <= 64; channels *= 2) {
                    for (bitreso = 1; bitreso <= 4; bitreso *= 2) {
                         char formatstring[256]; // <FORMAT>-<Mbps>-<nchan>-<nbit>
                         switch (format) {
                             case 0:
                                 snprintf(formatstring, 256, "MKIV1_%d-%d-%d-%d", fanout, rates[rate], channels, bitreso);
                                 break;
                             case 1:
                                 snprintf(formatstring, 256, "VLBA1_%d-%d-%d-%d", fanout, rates[rate], channels, bitreso);
                                 break;
                             default:
                                 snprintf(formatstring, 256, "unknown format=%d", format);
                         }
                         combinations_tested++;
                         printf("Trying '%s' : ", formatstring);
                         ms = new_mark5_stream(
                            new_mark5_stream_file(argv[1], 0),
                            new_mark5_format_generic_from_string(formatstring)
                         );
                         if (!ms) {
                             printf("FAIL\n");
                         } else {
                             int status = mark5_stream_decode(ms, 1024LL, unpacked);
	                         if (status < 0) {
                                 printf("DECODE FAIL\n");
                             } else {
                                 printf("success!! decode status = %d/1024\n", status);
                                 longlist = strcat(longlist, formatstring);
                                 longlist = strcat(longlist, "  ");
                                 mark5_stream_print(ms);
                             }
                         }
                     }                    
                }
            }
        }
   }

   printf("Tested %d combinations.\n", combinations_tested);
   if (strlen(longlist) > 0) {
       printf("Candidate list of modes that could open the file:\n%s\n\n", longlist);
       return 0;
   } else {
       printf("Tough luck, no formats and modes worked.\n\n");
   }
   return -1;
}
