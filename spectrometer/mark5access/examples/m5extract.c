/***************************************************************************
 *   Copyright (C) 2007 by Walter Brisken                                  *
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

#include <stdlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

#include "mark5access/mark5_stream.h"

const char program[] = "m5extract";
const char author[]  = "Walter Brisken / Jan Wagner";
const char version[] = "1.0";
const char verdate[] = "2007 Oct 6";

int usage(const char *pgm)
{
	printf("\n");

	printf("%s ver. %s   %s  %s\n\n", program, version, author, verdate);
	printf("A Mark5 decoder.  Can decode VLBA, Mark3/4, and Mark5B "
		"formats using the\nmark5access library. Extracts the "
                "specified channel and writes the samples\nas 8-bit data to an output file.\n\n");
	printf("Usage : %s <file> <dataformat> <n> <c> <outfile> [<offset>]\n\n", pgm);
	printf("  <file> is the name of the input file\n\n");
	printf("  <dataformat> should be of the form: "
		"<FORMAT>-<Mbps>-<nchan>-<nbit>, e.g.:\n");
	printf("    VLBA1_2-256-8-2\n");
	printf("    MKIV1_4-128-2-1\n");
	printf("    Mark5B-512-16-2\n\n");
	printf("  <n> is the number of samples per channel to decode (-1 to decode until EOF)\n\n");
	printf("  <c> is the channel (1..nchan) to extract\n\n");
	printf("  <offset> is number of bytes into file to start decoding\n\n");

	return 0;
}

int decode(const char *filename, const char *formatname, const char *f,
	long long offset, long long n, int ce, const char* outfilename)
{
	struct mark5_stream *ms;
	struct mark5_format_generic *m5form;
	struct mark5_stream_generic *m5strm;
	struct mark5_stream_generic *m5unpacker;
	struct stat fileStatus;
	float **data;
	char *rawdata;
	signed char* dataout;
	int i, j, k, status;
	long long chunk = 1024;
	size_t wrotebytes = 0;
	size_t readbytes = 0;
	ssize_t nread;
	int fpi, fpo;

    ms = new_mark5_stream(
            new_mark5_stream_file(filename, offset),
            new_mark5_format_generic_from_string(formatname) );

	fpi = open64(filename, O_RDONLY);
	if(!ms || fpi<0)
	{
		printf("problem opening %s\n", filename);
		return 0;
	}
	lseek64(fpi, offset, SEEK_SET);
	fstat(fpi, &fileStatus);

	fpo = open64(outfilename, O_CREAT|O_WRONLY, S_IWUSR|S_IRUSR|S_IRGRP);
	if (fpo<0)
	{
		printf("problem opening output file '%s': %s\n", outfilename, strerror(errno));
	}

	if (ce<1 || ce>ms->nchan)
	{
		ce = ms->nchan;
		printf("illegal channel nr specified, using #%d instead\n", ce);
	}

	data = (float **)malloc(ms->nchan*sizeof(float *));
	for(i = 0; i < ms->nchan; i++)
	{
		data[i] = (float *)malloc(chunk*sizeof(float ));
	}

	mark5_stream_print(ms);

	if(n % ms->samplegranularity > 0)
	{
		n -= (n % ms->samplegranularity);
		printf("Warning -- reducing read size to %Ld\n", n);
	}
	if(chunk % ms->samplegranularity > 0)
	{
		printf("Warning -- chunk %lld not granular to %d\n", chunk, ms->samplegranularity);
	}
	dataout = (signed char*)malloc(chunk);

    for(; n > 0; n -= chunk)
    {
        if(n < chunk) { chunk = n; }
        status = mark5_stream_decode(ms, chunk, data);

        if(status < 0)
        {
            printf("<EOF> status=%d\n", status);
            //break;
        }
        else
        {
            readbytes += chunk;
        }

        for(j = 0; j < chunk; j++)
        {
            #if 0
            for(k = 0; k < ms->nchan; k++)
            {
                printf(f, data[k][j]);
            }
            printf("\n");
            #endif
            dataout[j] = (signed char)(10.0 * data[ce][j]); // +-3.5 => +-35
        }

		nread = write(fpo, dataout, chunk*sizeof(char));
		if (nread<0)
		{
			printf("write() error '%s'\n", strerror(errno));
			break;
		}
		wrotebytes += nread;
	}

	printf("Read %lld raw bytes and wrote %lld 8-bit samples to '%s'\n", (long long)readbytes, (long long)wrotebytes, outfilename);

	for(i = 0; i < ms->nchan; i++)
	{
		free(data[i]);
	}
	free(data);

	delete_mark5_stream(ms);
	close(fpo);
	close(fpi);
	printf("\n");

	return 0;
}

int main(int argc, char **argv)
{
	long long offset = 0;
	long long n;
	int c;

	if(argc < 6)
	{
		return usage(program /*argv[0]*/);
	}

	n = atol(argv[3]);
	c = atoi(argv[4]);

	if(argc > 6)
	{
		offset=atoll(argv[6]);
	}

	decode(argv[1], argv[2], "%+2.0f ", offset, n, c, argv[5]);

	return 0;
}

