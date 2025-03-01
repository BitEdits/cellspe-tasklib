/***************************************************************************
 *   Copyright (C) 2009-2011 by Walter Brisken                             *
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
//===========================================================================
// SVN properties (DO NOT CHANGE)
//
// $Id: mark5_format_vlba_nomod.c,v 1.3 2012/02/21 13:35:44 jwagnerhki Exp $
// $HeadURL: https://svn.atnf.csiro.au/difx/libraries/mark5access/trunk/mark5access/mark5_format_vlba_nomod.c $
// $LastChangedRevision: 4106 $
// $Author: jwagnerhki $
// $LastChangedDate: 2011-12-09 21:48:51 +0100 (Fri, 09 Dec 2011) $
//
//============================================================================

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "mark5access/mark5_stream.h"

#define PAYLOADSIZE 20000

/* the high mag value for 2-bit reconstruction */
static const float HiMag = OPTIMAL_2BIT_HIGH;

static float lut1bit[256][8];  /* For all 1-bit modes */
static float lut2bit1[256][4]; /* fanout 1 @ 8/16t, fanout 4 @ 32/64t ! */
static float lut2bit2[256][4]; /* fanout 2 @ 8/16t, fanout 1 @ 32/64t   */
static float lut2bit3[256][4]; /* fanout 4 @ 8/16t, fanout 2 @ 32/64t   */
static float zeros[8];

/* ! 2bit/fanout4 use the following in decoding 32 and 64 track data: */

#ifdef WORDS_BIGENDIAN

#define reorder32(x) ((((x) & 0x55AA55AAUL)) | \
                      (((x) & 0xAA00AA00UL) >> 9) | \
		      (((x) & 0x00550055UL) << 9))
#define reorder64(x) ((((x) & 0x55AA55AA55AA55AAULL)) | \
                      (((x) & 0xAA00AA00AA00AA00ULL) >> 9) | \
		      (((x) & 0x0055005500550055ULL) << 9))

#else

#define reorder32(x) ((((x) & 0xAA55AA55UL)) | \
                      (((x) & 0x55005500UL) >> 7) | \
		      (((x) & 0x00AA00AAUL) << 7))
#define reorder64(x) ((((x) & 0xAA55AA55AA55AA55ULL)) | \
                      (((x) & 0x5500550055005500ULL) >> 7) | \
		      (((x) & 0x00AA00AA00AA00AAULL) << 7))

#endif

struct mark5_format_vlba_nomod
{
	int ntrack;
	int fanout;
	int kday;	/* kilo-mjd days.  51000, 52000, ... */
};

int countbits(unsigned char v);
int countbits32(unsigned int v);

static void initluts()
{
	int b, i, s, m, l;
	const float lut2level[2] = {1.0, -1.0};
	const float lut4level[4] = {-HiMag, 1.0, -1.0, HiMag};

	for(i = 0; i < 8; i++)
	{
		zeros[i] = 0.0;
	}

	for(b = 0; b < 256; b++)
	{
		/* lut1bit */
		for(i = 0; i < 8; i++)
		{
			l = (b>>i)&1;
			lut1bit[b][i] = lut2level[l];
		}

		/* lut2bit1 */
		for(i = 0; i < 4; i++)
		{
			s = i*2;	/* 0, 2, 4, 6 */
			m = s+1;	/* 1, 3, 5, 7 */
			l = ((b>>s)&1) + (((b>>m)&1)<<1);
			lut2bit1[b][i] = lut4level[l];
		}

		/* lut2bit2 */
		for(i = 0; i < 4; i++)
		{
			s = i+(i/2)*2;	/* 0, 1, 4, 5 */
			m = s+2;	/* 2, 3, 6, 7 */
			l = ((b>>s)&1) + (((b>>m)&1)<<1);
			lut2bit2[b][i] = lut4level[l];
		}

		/* lut2bit3 */
		for(i = 0; i < 4; i++)
		{
			s = i;		/* 0, 1, 2, 3 */
			m = s+4;	/* 4, 5, 6, 7 */
			l = ((b>>s)&1) + (((b>>m)&1)<<1);
			lut2bit3[b][i] = lut4level[l];
		}
	}
}

/* Look for the first occurrence (lowest offset >= 0) of the following pattern:
 *
 * 32*tracks bits set at offset bytes
 * 32*tracks bits set at offset+2520*tracks bytes
 * 2*tracks bits unset at offset+2518*tracks bytes
 *
 * With up to tol bits not set correctly
 *
 * return offset;
 */
static int findfirstframe(const unsigned char *data, int bytes, int tracks, int tol)
{
	int offset;
	int wrong = 0;
	int i, a, b;

	if(bytes < 2600*tracks)
	{
		return -1;
	}
	
	bytes -= 2600*tracks;

	b = tracks*2520;
	a = b - tracks/4;

	for(i = 0; i < 4*tracks; i++)
	{
		wrong += countbits(~data[i]);
		wrong += countbits(~data[i+b]);
	}
	for(i = 0; i < tracks/4; i++)
	{
		wrong += countbits(data[i+a]);
	}

	for(offset = 0; offset < bytes; offset++)
	{
		if(wrong < tol)
		{
			return offset;
		}
		wrong -= countbits(~data[offset]);
		wrong += countbits(~data[offset+4*tracks]);
		wrong -= countbits(~data[offset+b]);
		wrong += countbits(~data[offset+b+4*tracks]);
		wrong -= countbits(data[offset+a]);
		wrong += countbits(data[offset+a+tracks/4]);
	}

	return -1;
}

/* look at encoded nibbles.  Count bits in each track, assume set if
 * more than half are
 */
static void extractnibbles(const unsigned char *data, int ntracks, int numnibbles,
	char *nibbles)
{
	int i, j, b, n, c;

	n = ntracks/8;

	for(i = 0; i < numnibbles; i++)
	{
		nibbles[i] = 0;
		for(b = 0; b < 4; b++)
		{
			c = 0;
			for(j = 0; j < n; j++)
			{
				c += countbits(data[n*(4*i+3-b)+j]);
			}
			nibbles[i] += (c > n/2) ? (1 << b) : 0;
		}
	}
}

static int mark5_format_vlba_nomod_frame_time_int(const struct mark5_stream *ms,
	int *mjd, int *sec, int *ns)
{
	char nibs[12];
	struct mark5_format_vlba_nomod *v;
	int nRealTrack;

	if(!ms)
	{
		return -1;
	}
	v = (struct mark5_format_vlba_nomod *)(ms->formatdata);

	nRealTrack = v->ntrack;
	if(nRealTrack < 8)
	{
		nRealTrack = 8;
	}

	extractnibbles(ms->frame + 4*nRealTrack, nRealTrack, 12, nibs);

	if(mjd)
	{
		*mjd = v->kday + nibs[0]*100 + nibs[1]*10 + nibs[2];
	}
	if(sec)
	{
		*sec = nibs[3]*10000 + nibs[4]*1000 + nibs[5]*100
			+ nibs[6]*10 + nibs[7];
	}
	if(ns)
	{
		*ns = nibs[8]*100000000 + nibs[9]*10000000
			+ nibs[10]*1000000 + nibs[11]*100000;
	}

	return 0;
}

/* return in more general double value for ns */
static int mark5_format_vlba_nomod_frame_time(const struct mark5_stream *ms,
	int *mjd, int *sec, double *ns)
{
	int ins, v;

	v = mark5_format_vlba_nomod_frame_time_int(ms, mjd, sec, &ins);

	*ns = ins;

	return v;
}

static int mark5_format_vlba_nomod_validate(const struct mark5_stream *ms)
{
	struct mark5_format_vlba_nomod *v;
	int ntrack, t, e=0;
	unsigned int *data;
	int mjd_d, mjd_t, sec_d, sec_t;
	int ns_d;
	long long ns_t;

	if(!ms)
	{
		fprintf(m5stdout, "mark5_format_vlba_nomod_validate: ms=0\n");
		return 0;
	}

	v = (struct mark5_format_vlba_nomod *)(ms->formatdata);
	ntrack = v->ntrack;
	data = (unsigned int *)ms->frame;
	for(t = 0; t < ntrack; t++)
	{
		/* allow 1 of every 32 bits to be incorrect */
		if(countbits32(data[t]) < 31)
		{
//			fprintf(m5stdout, "<%s %d %d>", ms->streamname, t, data[t]);
			e++;
		}
	}

	if(e > 0)
	{
//		fprintf(m5stdout, "mark5_format_vlba_nomod_validate[%s]: e=%d\n", ms->streamname, e);
		return 0;
	}

	if(ms->mjd && ms->framenum % ms->framegranularity == 0)
	{
		mark5_format_vlba_nomod_frame_time_int(ms, &mjd_d, &sec_d, &ns_d);

		ns_t = (long long)(ms->framenum)*(long long)(ms->gframens/ms->framegranularity) + (long long)(ms->ns);
		sec_t = ns_t / 1000000000L;
		ns_t -= (long long)sec_t * 1000000000L;
		sec_t += ms->sec;
		mjd_t = sec_t / 86400;
		sec_t -= mjd_t * 86400;
		mjd_t += ms->mjd;

		if(mjd_t != mjd_d || sec_t != sec_d || ns_t != ns_d)
		{
			return 0;
		}
	}

	return 1;
}

static int mark5_format_vlba_nomod_fixmjd(struct mark5_stream *ms, int refmjd)
{
	struct mark5_format_vlba_nomod *v;
	int n;

	if(!ms)
	{
		return -1;
	}

	v = (struct mark5_format_vlba_nomod *)(ms->formatdata);
	if(v->kday == 0)
	{
		n = (refmjd - ms->mjd + 500) / 1000;
		ms->mjd += n*1000;
		v->kday = n*1000;

		return 1;
	}

	return 0;
}

/*********************** data unpack routines **********************/

/* NOTE: all decimation4 decoders work for decimation = 2^k for k >= 2 */

static int vlba_nomod_decode_1bit_1track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_1track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += 2;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_1track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_2track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[1];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_2track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += 2;

		data[0][o] = fp[0];
		data[1][o] = fp[1];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_2track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];
		data[1][o] = fp[1];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_2track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		o++;
		data[0][o] = fp[1];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_2track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_2track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_4track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[1];
		data[2][o] = fp[2];
		data[3][o] = fp[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_4track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += 2;

		data[0][o] = fp[0];
		data[1][o] = fp[1];
		data[2][o] = fp[2];
		data[3][o] = fp[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_4track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];
		data[1][o] = fp[1];
		data[2][o] = fp[2];
		data[3][o] = fp[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_4track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[2];
		o++;
		data[0][o] = fp[1];
		data[1][o] = fp[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_4track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[2];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_4track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];
		data[1][o] = fp[2];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_4track_fanout4_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		o++;
		data[0][o] = fp[1];
		o++;
		data[0][o] = fp[2];
		o++;
		data[0][o] = fp[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_4track_fanout4_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		o++;
		data[0][o] = fp[2];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_4track_fanout4_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_8track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[1];
		data[2][o] = fp[2];
		data[3][o] = fp[3];
		data[4][o] = fp[4];
		data[5][o] = fp[5];
		data[6][o] = fp[6];
		data[7][o] = fp[7];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_8track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += 2;

		data[0][o] = fp[0];
		data[1][o] = fp[1];
		data[2][o] = fp[2];
		data[3][o] = fp[3];
		data[4][o] = fp[4];
		data[5][o] = fp[5];
		data[6][o] = fp[6];
		data[7][o] = fp[7];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_8track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];
		data[1][o] = fp[1];
		data[2][o] = fp[2];
		data[3][o] = fp[3];
		data[4][o] = fp[4];
		data[5][o] = fp[5];
		data[6][o] = fp[6];
		data[7][o] = fp[7];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_8track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[2];
		data[2][o] = fp[4];
		data[3][o] = fp[6];
		o++;
		data[0][o] = fp[1];
		data[1][o] = fp[3];
		data[2][o] = fp[5];
		data[3][o] = fp[7];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_8track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[2];
		data[2][o] = fp[4];
		data[3][o] = fp[6];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_8track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];
		data[1][o] = fp[2];
		data[2][o] = fp[4];
		data[3][o] = fp[6];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_8track_fanout4_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[4];
		o++;
		data[0][o] = fp[1];
		data[1][o] = fp[5];
		o++;
		data[0][o] = fp[2];
		data[1][o] = fp[6];
		o++;
		data[0][o] = fp[3];
		data[1][o] = fp[7];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_8track_fanout4_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[4];
		o++;
		data[0][o] = fp[2];
		data[1][o] = fp[6];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_8track_fanout4_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut1bit[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];
		data[1][o] = fp[4];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_16track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 2;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[1];
		data[2][o]  = fp0[2];
		data[3][o]  = fp0[3];
		data[4][o]  = fp0[4];
		data[5][o]  = fp0[5];
		data[6][o]  = fp0[6];
		data[7][o]  = fp0[7];
		data[8][o]  = fp1[0];
		data[9][o]  = fp1[1];
		data[10][o] = fp1[2];
		data[11][o] = fp1[3];
		data[12][o] = fp1[4];
		data[13][o] = fp1[5];
		data[14][o] = fp1[6];
		data[15][o] = fp1[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_16track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 4;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i += 3;
		}
		m += 2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[1];
		data[2][o]  = fp0[2];
		data[3][o]  = fp0[3];
		data[4][o]  = fp0[4];
		data[5][o]  = fp0[5];
		data[6][o]  = fp0[6];
		data[7][o]  = fp0[7];
		data[8][o]  = fp1[0];
		data[9][o]  = fp1[1];
		data[10][o] = fp1[2];
		data[11][o] = fp1[3];
		data[12][o] = fp1[4];
		data[13][o] = fp1[5];
		data[14][o] = fp1[6];
		data[15][o] = fp1[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_16track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;
	df = ms->decimation*2 - 1;
	df2 = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i++;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[1];
		data[2][o]  = fp0[2];
		data[3][o]  = fp0[3];
		data[4][o]  = fp0[4];
		data[5][o]  = fp0[5];
		data[6][o]  = fp0[6];
		data[7][o]  = fp0[7];
		data[8][o]  = fp1[0];
		data[9][o]  = fp1[1];
		data[10][o] = fp1[2];
		data[11][o] = fp1[3];
		data[12][o] = fp1[4];
		data[13][o] = fp1[5];
		data[14][o] = fp1[6];
		data[15][o] = fp1[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_16track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 2;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
		}
		m++;
		
		data[0][o] = fp0[0];
		data[1][o] = fp0[2];
		data[2][o] = fp0[4];
		data[3][o] = fp0[6];
		data[4][o] = fp1[0];
		data[5][o] = fp1[2];
		data[6][o] = fp1[4];
		data[7][o] = fp1[6];
		o++;
		data[0][o] = fp0[1];
		data[1][o] = fp0[3];
		data[2][o] = fp0[5];
		data[3][o] = fp0[7];
		data[4][o] = fp1[1];
		data[5][o] = fp1[3];
		data[6][o] = fp1[5];
		data[7][o] = fp1[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_16track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 2;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
		}
		m++;
		
		data[0][o] = fp0[0];
		data[1][o] = fp0[2];
		data[2][o] = fp0[4];
		data[3][o] = fp0[6];
		data[4][o] = fp1[0];
		data[5][o] = fp1[2];
		data[6][o] = fp1[4];
		data[7][o] = fp1[6];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_16track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;
	df = ms->decimation - 1;
	df2 = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i++;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
		}
		i += df;
		m += df2;
		
		data[0][o] = fp0[0];
		data[1][o] = fp0[2];
		data[2][o] = fp0[4];
		data[3][o] = fp0[6];
		data[4][o] = fp1[0];
		data[5][o] = fp1[2];
		data[6][o] = fp1[4];
		data[7][o] = fp1[6];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_16track_fanout4_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 2;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o] = fp0[0];
		data[1][o] = fp0[4];
		data[2][o] = fp1[0];
		data[3][o] = fp1[4];
		o++;
		data[0][o] = fp0[1];
		data[1][o] = fp0[5];
		data[2][o] = fp1[1];
		data[3][o] = fp1[5];
		o++;
		data[0][o] = fp0[2];
		data[1][o] = fp0[6];
		data[2][o] = fp1[2];
		data[3][o] = fp1[6];
		o++;
		data[0][o] = fp0[3];
		data[1][o] = fp0[7];
		data[2][o] = fp1[3];
		data[3][o] = fp1[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_16track_fanout4_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 2;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o] = fp0[0];
		data[1][o] = fp0[4];
		data[2][o] = fp1[0];
		data[3][o] = fp1[4];
		o++;
		data[0][o] = fp0[2];
		data[1][o] = fp0[6];
		data[2][o] = fp1[2];
		data[3][o] = fp1[6];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_16track_fanout4_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;
	df = ms->decimation/2 - 1;
	df2 = ms->decimation/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i++;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o] = fp0[0];
		data[1][o] = fp0[4];
		data[2][o] = fp1[0];
		data[3][o] = fp1[4];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_32track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 4;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp0[4];
		data[3][o]  = fp0[6];
		data[4][o]  = fp1[0];
		data[5][o]  = fp1[2];
		data[6][o]  = fp1[4];
		data[7][o]  = fp1[6];
		data[8][o]  = fp2[0];
		data[9][o]  = fp2[2];
		data[10][o] = fp2[4];
		data[11][o] = fp2[6];
		data[12][o] = fp3[0];
		data[13][o] = fp3[2];
		data[14][o] = fp3[4];
		data[15][o] = fp3[6];
		data[16][o] = fp0[1];
		data[17][o] = fp0[3];
		data[18][o] = fp0[5];
		data[19][o] = fp0[7];
		data[20][o] = fp1[1];
		data[21][o] = fp1[3];
		data[22][o] = fp1[5];
		data[23][o] = fp1[7];
		data[24][o] = fp2[1];
		data[25][o] = fp2[3];
		data[26][o] = fp2[5];
		data[27][o] = fp2[7];
		data[28][o] = fp3[1];
		data[29][o] = fp3[3];
		data[30][o] = fp3[5];
		data[31][o] = fp3[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_32track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 8;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i += 5;
		}
		m += 2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp0[4];
		data[3][o]  = fp0[6];
		data[4][o]  = fp1[0];
		data[5][o]  = fp1[2];
		data[6][o]  = fp1[4];
		data[7][o]  = fp1[6];
		data[8][o]  = fp2[0];
		data[9][o]  = fp2[2];
		data[10][o] = fp2[4];
		data[11][o] = fp2[6];
		data[12][o] = fp3[0];
		data[13][o] = fp3[2];
		data[14][o] = fp3[4];
		data[15][o] = fp3[6];
		data[16][o] = fp0[1];
		data[17][o] = fp0[3];
		data[18][o] = fp0[5];
		data[19][o] = fp0[7];
		data[20][o] = fp1[1];
		data[21][o] = fp1[3];
		data[22][o] = fp1[5];
		data[23][o] = fp1[7];
		data[24][o] = fp2[1];
		data[25][o] = fp2[3];
		data[26][o] = fp2[5];
		data[27][o] = fp2[7];
		data[28][o] = fp3[1];
		data[29][o] = fp3[3];
		data[30][o] = fp3[5];
		data[31][o] = fp3[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_32track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;
	df = ms->decimation*4 - 3;
	df2 = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 3;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp0[4];
		data[3][o]  = fp0[6];
		data[4][o]  = fp1[0];
		data[5][o]  = fp1[2];
		data[6][o]  = fp1[4];
		data[7][o]  = fp1[6];
		data[8][o]  = fp2[0];
		data[9][o]  = fp2[2];
		data[10][o] = fp2[4];
		data[11][o] = fp2[6];
		data[12][o] = fp3[0];
		data[13][o] = fp3[2];
		data[14][o] = fp3[4];
		data[15][o] = fp3[6];
		data[16][o] = fp0[1];
		data[17][o] = fp0[3];
		data[18][o] = fp0[5];
		data[19][o] = fp0[7];
		data[20][o] = fp1[1];
		data[21][o] = fp1[3];
		data[22][o] = fp1[5];
		data[23][o] = fp1[7];
		data[24][o] = fp2[1];
		data[25][o] = fp2[3];
		data[26][o] = fp2[5];
		data[27][o] = fp2[7];
		data[28][o] = fp3[1];
		data[29][o] = fp3[3];
		data[30][o] = fp3[5];
		data[31][o] = fp3[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_32track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 4;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[4];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[4];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[4];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[4];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[5];
		data[10][o] = fp1[1];
		data[11][o] = fp1[5];
		data[12][o] = fp2[1];
		data[13][o] = fp2[5];
		data[14][o] = fp3[1];
		data[15][o] = fp3[5];
		o++;
		data[0][o]  = fp0[2];
		data[1][o]  = fp0[6];
		data[2][o]  = fp1[2];
		data[3][o]  = fp1[6];
		data[4][o]  = fp2[2];
		data[5][o]  = fp2[6];
		data[6][o]  = fp3[2];
		data[7][o]  = fp3[6];
		data[8][o]  = fp0[3];
		data[9][o]  = fp0[7];
		data[10][o] = fp1[3];
		data[11][o] = fp1[7];
		data[12][o] = fp2[3];
		data[13][o] = fp2[7];
		data[14][o] = fp3[3];
		data[15][o] = fp3[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_32track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 4;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[4];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[4];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[4];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[4];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[5];
		data[10][o] = fp1[1];
		data[11][o] = fp1[5];
		data[12][o] = fp2[1];
		data[13][o] = fp2[5];
		data[14][o] = fp3[1];
		data[15][o] = fp3[5];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_32track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;
	df = ms->decimation*2 - 3;
	df2 = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 3;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[4];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[4];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[4];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[4];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[5];
		data[10][o] = fp1[1];
		data[11][o] = fp1[5];
		data[12][o] = fp2[1];
		data[13][o] = fp2[5];
		data[14][o] = fp3[1];
		data[15][o] = fp3[5];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_32track_fanout4_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 4;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o] = fp0[0];
		data[1][o] = fp1[0];
		data[2][o] = fp2[0];
		data[3][o] = fp3[0];
		data[4][o] = fp0[1];
		data[5][o] = fp1[1];
		data[6][o] = fp2[1];
		data[7][o] = fp3[1];
		o++;
		data[0][o] = fp0[2];
		data[1][o] = fp1[2];
		data[2][o] = fp2[2];
		data[3][o] = fp3[2];
		data[4][o] = fp0[3];
		data[5][o] = fp1[3];
		data[6][o] = fp2[3];
		data[7][o] = fp3[3];
		o++;
		data[0][o] = fp0[4];
		data[1][o] = fp1[4];
		data[2][o] = fp2[4];
		data[3][o] = fp3[4];
		data[4][o] = fp0[5];
		data[5][o] = fp1[5];
		data[6][o] = fp2[5];
		data[7][o] = fp3[5];
		o++;
		data[0][o] = fp0[6];
		data[1][o] = fp1[6];
		data[2][o] = fp2[6];
		data[3][o] = fp3[6];
		data[4][o] = fp0[7];
		data[5][o] = fp1[7];
		data[6][o] = fp2[7];
		data[7][o] = fp3[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_32track_fanout4_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 4;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o] = fp0[0];
		data[1][o] = fp1[0];
		data[2][o] = fp2[0];
		data[3][o] = fp3[0];
		data[4][o] = fp0[1];
		data[5][o] = fp1[1];
		data[6][o] = fp2[1];
		data[7][o] = fp3[1];
		o++;
		data[0][o] = fp0[4];
		data[1][o] = fp1[4];
		data[2][o] = fp2[4];
		data[3][o] = fp3[4];
		data[4][o] = fp0[5];
		data[5][o] = fp1[5];
		data[6][o] = fp2[5];
		data[7][o] = fp3[5];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_32track_fanout4_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;
	df = ms->decimation - 3;
	df2 = ms->decimation/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 3;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o] = fp0[0];
		data[1][o] = fp1[0];
		data[2][o] = fp2[0];
		data[3][o] = fp3[0];
		data[4][o] = fp0[1];
		data[5][o] = fp1[1];
		data[6][o] = fp2[1];
		data[7][o] = fp3[1];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_64track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 8;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
			fp4 = lut1bit[buf[i]];
			i++;
			fp5 = lut1bit[buf[i]];
			i++;
			fp6 = lut1bit[buf[i]];
			i++;
			fp7 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp0[4];
		data[3][o]  = fp0[6];
		data[4][o]  = fp1[0];
		data[5][o]  = fp1[2];
		data[6][o]  = fp1[4];
		data[7][o]  = fp1[6];
		data[8][o]  = fp2[0];
		data[9][o]  = fp2[2];
		data[10][o] = fp2[4];
		data[11][o] = fp2[6];
		data[12][o] = fp3[0];
		data[13][o] = fp3[2];
		data[14][o] = fp3[4];
		data[15][o] = fp3[6];
		data[16][o] = fp0[1];
		data[17][o] = fp0[3];
		data[18][o] = fp0[5];
		data[19][o] = fp0[7];
		data[20][o] = fp1[1];
		data[21][o] = fp1[3];
		data[22][o] = fp1[5];
		data[23][o] = fp1[7];
		data[24][o] = fp2[1];
		data[25][o] = fp2[3];
		data[26][o] = fp2[5];
		data[27][o] = fp2[7];
		data[28][o] = fp3[1];
		data[29][o] = fp3[3];
		data[30][o] = fp3[5];
		data[31][o] = fp3[7];
		data[32][o] = fp4[0];
		data[33][o] = fp4[2];
		data[34][o] = fp4[4];
		data[35][o] = fp4[6];
		data[36][o] = fp5[0];
		data[37][o] = fp5[2];
		data[38][o] = fp5[4];
		data[39][o] = fp5[6];
		data[40][o] = fp6[0];
		data[41][o] = fp6[2];
		data[42][o] = fp6[4];
		data[43][o] = fp6[6];
		data[44][o] = fp7[0];
		data[45][o] = fp7[2];
		data[46][o] = fp7[4];
		data[47][o] = fp7[6];
		data[48][o] = fp4[1];
		data[49][o] = fp4[3];
		data[50][o] = fp4[5];
		data[51][o] = fp4[7];
		data[52][o] = fp5[1];
		data[53][o] = fp5[3];
		data[54][o] = fp5[5];
		data[55][o] = fp5[7];
		data[56][o] = fp6[1];
		data[57][o] = fp6[3];
		data[58][o] = fp6[5];
		data[59][o] = fp6[7];
		data[60][o] = fp7[1];
		data[61][o] = fp7[3];
		data[62][o] = fp7[5];
		data[63][o] = fp7[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_64track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 16;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
			fp4 = lut1bit[buf[i]];
			i++;
			fp5 = lut1bit[buf[i]];
			i++;
			fp6 = lut1bit[buf[i]];
			i++;
			fp7 = lut1bit[buf[i]];
			i += 9;
		}
		m += 2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp0[4];
		data[3][o]  = fp0[6];
		data[4][o]  = fp1[0];
		data[5][o]  = fp1[2];
		data[6][o]  = fp1[4];
		data[7][o]  = fp1[6];
		data[8][o]  = fp2[0];
		data[9][o]  = fp2[2];
		data[10][o] = fp2[4];
		data[11][o] = fp2[6];
		data[12][o] = fp3[0];
		data[13][o] = fp3[2];
		data[14][o] = fp3[4];
		data[15][o] = fp3[6];
		data[16][o] = fp0[1];
		data[17][o] = fp0[3];
		data[18][o] = fp0[5];
		data[19][o] = fp0[7];
		data[20][o] = fp1[1];
		data[21][o] = fp1[3];
		data[22][o] = fp1[5];
		data[23][o] = fp1[7];
		data[24][o] = fp2[1];
		data[25][o] = fp2[3];
		data[26][o] = fp2[5];
		data[27][o] = fp2[7];
		data[28][o] = fp3[1];
		data[29][o] = fp3[3];
		data[30][o] = fp3[5];
		data[31][o] = fp3[7];
		data[32][o] = fp4[0];
		data[33][o] = fp4[2];
		data[34][o] = fp4[4];
		data[35][o] = fp4[6];
		data[36][o] = fp5[0];
		data[37][o] = fp5[2];
		data[38][o] = fp5[4];
		data[39][o] = fp5[6];
		data[40][o] = fp6[0];
		data[41][o] = fp6[2];
		data[42][o] = fp6[4];
		data[43][o] = fp6[6];
		data[44][o] = fp7[0];
		data[45][o] = fp7[2];
		data[46][o] = fp7[4];
		data[47][o] = fp7[6];
		data[48][o] = fp4[1];
		data[49][o] = fp4[3];
		data[50][o] = fp4[5];
		data[51][o] = fp4[7];
		data[52][o] = fp5[1];
		data[53][o] = fp5[3];
		data[54][o] = fp5[5];
		data[55][o] = fp5[7];
		data[56][o] = fp6[1];
		data[57][o] = fp6[3];
		data[58][o] = fp6[5];
		data[59][o] = fp6[7];
		data[60][o] = fp7[1];
		data[61][o] = fp7[3];
		data[62][o] = fp7[5];
		data[63][o] = fp7[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_64track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;
	df = ms->decimation*8 - 7;
	df2 = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 7;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
			fp4 = lut1bit[buf[i]];
			i++;
			fp5 = lut1bit[buf[i]];
			i++;
			fp6 = lut1bit[buf[i]];
			i++;
			fp7 = lut1bit[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp0[4];
		data[3][o]  = fp0[6];
		data[4][o]  = fp1[0];
		data[5][o]  = fp1[2];
		data[6][o]  = fp1[4];
		data[7][o]  = fp1[6];
		data[8][o]  = fp2[0];
		data[9][o]  = fp2[2];
		data[10][o] = fp2[4];
		data[11][o] = fp2[6];
		data[12][o] = fp3[0];
		data[13][o] = fp3[2];
		data[14][o] = fp3[4];
		data[15][o] = fp3[6];
		data[16][o] = fp0[1];
		data[17][o] = fp0[3];
		data[18][o] = fp0[5];
		data[19][o] = fp0[7];
		data[20][o] = fp1[1];
		data[21][o] = fp1[3];
		data[22][o] = fp1[5];
		data[23][o] = fp1[7];
		data[24][o] = fp2[1];
		data[25][o] = fp2[3];
		data[26][o] = fp2[5];
		data[27][o] = fp2[7];
		data[28][o] = fp3[1];
		data[29][o] = fp3[3];
		data[30][o] = fp3[5];
		data[31][o] = fp3[7];
		data[32][o] = fp4[0];
		data[33][o] = fp4[2];
		data[34][o] = fp4[4];
		data[35][o] = fp4[6];
		data[36][o] = fp5[0];
		data[37][o] = fp5[2];
		data[38][o] = fp5[4];
		data[39][o] = fp5[6];
		data[40][o] = fp6[0];
		data[41][o] = fp6[2];
		data[42][o] = fp6[4];
		data[43][o] = fp6[6];
		data[44][o] = fp7[0];
		data[45][o] = fp7[2];
		data[46][o] = fp7[4];
		data[47][o] = fp7[6];
		data[48][o] = fp4[1];
		data[49][o] = fp4[3];
		data[50][o] = fp4[5];
		data[51][o] = fp4[7];
		data[52][o] = fp5[1];
		data[53][o] = fp5[3];
		data[54][o] = fp5[5];
		data[55][o] = fp5[7];
		data[56][o] = fp6[1];
		data[57][o] = fp6[3];
		data[58][o] = fp6[5];
		data[59][o] = fp6[7];
		data[60][o] = fp7[1];
		data[61][o] = fp7[3];
		data[62][o] = fp7[5];
		data[63][o] = fp7[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_1bit_64track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 8;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
			fp4 = lut1bit[buf[i]];
			i++;
			fp5 = lut1bit[buf[i]];
			i++;
			fp6 = lut1bit[buf[i]];
			i++;
			fp7 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[4];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[4];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[4];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[4];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[5];
		data[10][o] = fp1[1];
		data[11][o] = fp1[5];
		data[12][o] = fp2[1];
		data[13][o] = fp2[5];
		data[14][o] = fp3[1];
		data[15][o] = fp3[5];
		data[16][o] = fp4[0];
		data[17][o] = fp4[4];
		data[18][o] = fp5[0];
		data[19][o] = fp5[4];
		data[20][o] = fp6[0];
		data[21][o] = fp6[4];
		data[22][o] = fp7[0];
		data[23][o] = fp7[4];
		data[24][o] = fp4[1];
		data[25][o] = fp4[5];
		data[26][o] = fp5[1];
		data[27][o] = fp5[5];
		data[28][o] = fp6[1];
		data[29][o] = fp6[5];
		data[30][o] = fp7[1];
		data[31][o] = fp7[5];
		o++;
		data[0][o]  = fp0[2];
		data[1][o]  = fp0[6];
		data[2][o]  = fp1[2];
		data[3][o]  = fp1[6];
		data[4][o]  = fp2[2];
		data[5][o]  = fp2[6];
		data[6][o]  = fp3[2];
		data[7][o]  = fp3[6];
		data[8][o]  = fp0[3];
		data[9][o]  = fp0[7];
		data[10][o] = fp1[3];
		data[11][o] = fp1[7];
		data[12][o] = fp2[3];
		data[13][o] = fp2[7];
		data[14][o] = fp3[3];
		data[15][o] = fp3[7];
		data[16][o] = fp4[2];
		data[17][o] = fp4[6];
		data[18][o] = fp5[2];
		data[19][o] = fp5[6];
		data[20][o] = fp6[2];
		data[21][o] = fp6[6];
		data[22][o] = fp7[2];
		data[23][o] = fp7[6];
		data[24][o] = fp4[3];
		data[25][o] = fp4[7];
		data[26][o] = fp5[3];
		data[27][o] = fp5[7];
		data[28][o] = fp6[3];
		data[29][o] = fp6[7];
		data[30][o] = fp7[3];
		data[31][o] = fp7[7];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_64track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 8;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
			fp4 = lut1bit[buf[i]];
			i++;
			fp5 = lut1bit[buf[i]];
			i++;
			fp6 = lut1bit[buf[i]];
			i++;
			fp7 = lut1bit[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[4];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[4];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[4];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[4];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[5];
		data[10][o] = fp1[1];
		data[11][o] = fp1[5];
		data[12][o] = fp2[1];
		data[13][o] = fp2[5];
		data[14][o] = fp3[1];
		data[15][o] = fp3[5];
		data[16][o] = fp4[0];
		data[17][o] = fp4[4];
		data[18][o] = fp5[0];
		data[19][o] = fp5[4];
		data[20][o] = fp6[0];
		data[21][o] = fp6[4];
		data[22][o] = fp7[0];
		data[23][o] = fp7[4];
		data[24][o] = fp4[1];
		data[25][o] = fp4[5];
		data[26][o] = fp5[1];
		data[27][o] = fp5[5];
		data[28][o] = fp6[1];
		data[29][o] = fp6[5];
		data[30][o] = fp7[1];
		data[31][o] = fp7[5];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_64track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;
	df = ms->decimation*4 - 7;
	df2 = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 7;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
			fp4 = lut1bit[buf[i]];
			i++;
			fp5 = lut1bit[buf[i]];
			i++;
			fp6 = lut1bit[buf[i]];
			i++;
			fp7 = lut1bit[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[4];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[4];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[4];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[4];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[5];
		data[10][o] = fp1[1];
		data[11][o] = fp1[5];
		data[12][o] = fp2[1];
		data[13][o] = fp2[5];
		data[14][o] = fp3[1];
		data[15][o] = fp3[5];
		data[16][o] = fp4[0];
		data[17][o] = fp4[4];
		data[18][o] = fp5[0];
		data[19][o] = fp5[4];
		data[20][o] = fp6[0];
		data[21][o] = fp6[4];
		data[22][o] = fp7[0];
		data[23][o] = fp7[4];
		data[24][o] = fp4[1];
		data[25][o] = fp4[5];
		data[26][o] = fp5[1];
		data[27][o] = fp5[5];
		data[28][o] = fp6[1];
		data[29][o] = fp6[5];
		data[30][o] = fp7[1];
		data[31][o] = fp7[5];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_1bit_64track_fanout4_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 8;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
			fp4 = lut1bit[buf[i]];
			i++;
			fp5 = lut1bit[buf[i]];
			i++;
			fp6 = lut1bit[buf[i]];
			i++;
			fp7 = lut1bit[buf[i]];
			i++;
		}
		m++;
		
		data[0][o]  = fp0[0];
		data[1][o]  = fp1[0];
		data[2][o]  = fp2[0];
		data[3][o]  = fp3[0];
		data[4][o]  = fp0[1];
		data[5][o]  = fp1[1];
		data[6][o]  = fp2[1];
		data[7][o]  = fp3[1];
		data[8][o]  = fp4[0];
		data[9][o]  = fp5[0];
		data[10][o] = fp6[0];
		data[11][o] = fp7[0];
		data[12][o] = fp4[1];
		data[13][o] = fp5[1];
		data[14][o] = fp6[1];
		data[15][o] = fp7[1];
		o++;
		data[0][o]  = fp0[2];
		data[1][o]  = fp1[2];
		data[2][o]  = fp2[2];
		data[3][o]  = fp3[2];
		data[4][o]  = fp0[3];
		data[5][o]  = fp1[3];
		data[6][o]  = fp2[3];
		data[7][o]  = fp3[3];
		data[8][o]  = fp4[2];
		data[9][o]  = fp5[2];
		data[10][o] = fp6[2];
		data[11][o] = fp7[2];
		data[12][o] = fp4[3];
		data[13][o] = fp5[3];
		data[14][o] = fp6[3];
		data[15][o] = fp7[3];
		o++;
		data[0][o]  = fp0[4];
		data[1][o]  = fp1[4];
		data[2][o]  = fp2[4];
		data[3][o]  = fp3[4];
		data[4][o]  = fp0[5];
		data[5][o]  = fp1[5];
		data[6][o]  = fp2[5];
		data[7][o]  = fp3[5];
		data[8][o]  = fp4[4];
		data[9][o]  = fp5[4];
		data[10][o] = fp6[4];
		data[11][o] = fp7[4];
		data[12][o] = fp4[5];
		data[13][o] = fp5[5];
		data[14][o] = fp6[5];
		data[15][o] = fp7[5];
		o++;
		data[0][o]  = fp0[6];
		data[1][o]  = fp1[6];
		data[2][o]  = fp2[6];
		data[3][o]  = fp3[6];
		data[4][o]  = fp0[7];
		data[5][o]  = fp1[7];
		data[6][o]  = fp2[7];
		data[7][o]  = fp3[7];
		data[8][o]  = fp4[6];
		data[9][o]  = fp5[6];
		data[10][o] = fp6[6];
		data[11][o] = fp7[6];
		data[12][o] = fp4[7];
		data[13][o] = fp5[7];
		data[14][o] = fp6[7];
		data[15][o] = fp7[7];


		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_64track_fanout4_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 8;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
			fp4 = lut1bit[buf[i]];
			i++;
			fp5 = lut1bit[buf[i]];
			i++;
			fp6 = lut1bit[buf[i]];
			i++;
			fp7 = lut1bit[buf[i]];
			i++;
		}
		m++;
		
		data[0][o]  = fp0[0];
		data[1][o]  = fp1[0];
		data[2][o]  = fp2[0];
		data[3][o]  = fp3[0];
		data[4][o]  = fp0[1];
		data[5][o]  = fp1[1];
		data[6][o]  = fp2[1];
		data[7][o]  = fp3[1];
		data[8][o]  = fp4[0];
		data[9][o]  = fp5[0];
		data[10][o] = fp6[0];
		data[11][o] = fp7[0];
		data[12][o] = fp4[1];
		data[13][o] = fp5[1];
		data[14][o] = fp6[1];
		data[15][o] = fp7[1];
		o++;
		data[0][o]  = fp0[4];
		data[1][o]  = fp1[4];
		data[2][o]  = fp2[4];
		data[3][o]  = fp3[4];
		data[4][o]  = fp0[5];
		data[5][o]  = fp1[5];
		data[6][o]  = fp2[5];
		data[7][o]  = fp3[5];
		data[8][o]  = fp4[4];
		data[9][o]  = fp5[4];
		data[10][o] = fp6[4];
		data[11][o] = fp7[4];
		data[12][o] = fp4[5];
		data[13][o] = fp5[5];
		data[14][o] = fp6[5];
		data[15][o] = fp7[5];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_1bit_64track_fanout4_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;
	df = ms->decimation*2 - 7;
	df2 = ms->decimation/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 7;
			nblank++;
		}
		else
		{
			fp0 = lut1bit[buf[i]];
			i++;
			fp1 = lut1bit[buf[i]];
			i++;
			fp2 = lut1bit[buf[i]];
			i++;
			fp3 = lut1bit[buf[i]];
			i++;
			fp4 = lut1bit[buf[i]];
			i++;
			fp5 = lut1bit[buf[i]];
			i++;
			fp6 = lut1bit[buf[i]];
			i++;
			fp7 = lut1bit[buf[i]];
		}
		i += df;
		m += df2;
		
		data[0][o]  = fp0[0];
		data[1][o]  = fp1[0];
		data[2][o]  = fp2[0];
		data[3][o]  = fp3[0];
		data[4][o]  = fp0[1];
		data[5][o]  = fp1[1];
		data[6][o]  = fp2[1];
		data[7][o]  = fp3[1];
		data[8][o]  = fp4[0];
		data[9][o]  = fp5[0];
		data[10][o] = fp6[0];
		data[11][o] = fp7[0];
		data[12][o] = fp4[1];
		data[13][o] = fp5[1];
		data[14][o] = fp6[1];
		data[15][o] = fp7[1];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

/************************ 2-bit decoders *********************/

static int vlba_nomod_decode_2bit_2track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit1[buf[i]];
		}
		i++;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_2track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit1[buf[i]];
		}
		i += 2;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_2track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit1[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_4track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit1[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[1];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_4track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit1[buf[i]];
		}
		i += 2;

		data[0][o] = fp[0];
		data[1][o] = fp[1];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_4track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation;
	
	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit1[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];
		data[1][o] = fp[1];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_4track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit2[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		o++;
		data[0][o] = fp[1];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_4track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit2[buf[i]];
		}
		i++;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_4track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit2[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_8track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit1[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[1];
		data[2][o] = fp[2];
		data[3][o] = fp[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}


static int vlba_nomod_decode_2bit_8track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit1[buf[i]];
		}
		i += 2;

		data[0][o] = fp[0];
		data[1][o] = fp[1];
		data[2][o] = fp[2];
		data[3][o] = fp[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}


static int vlba_nomod_decode_2bit_8track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit1[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];
		data[1][o] = fp[1];
		data[2][o] = fp[2];
		data[3][o] = fp[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_8track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit2[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[2];
		o++;
		data[0][o] = fp[1];
		data[1][o] = fp[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_8track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit2[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		data[1][o] = fp[2];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_8track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit2[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];
		data[1][o] = fp[2];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_8track_fanout4_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit3[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		o++;
		data[0][o] = fp[1];
		o++;
		data[0][o] = fp[2];
		o++;
		data[0][o] = fp[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_8track_fanout4_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit3[buf[i]];
		}
		i++;

		data[0][o] = fp[0];
		o++;
		data[0][o] = fp[2];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_8track_fanout4_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp;
	int o, i, df;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	df = ms->decimation/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp = zeros;
			nblank++;
		}
		else
		{
			fp = lut2bit3[buf[i]];
		}
		i += df;

		data[0][o] = fp[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_16track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 2;
			nblank++;
		}
		else
		{
			fp0 = lut2bit1[buf[i]];
			i++;
			fp1 = lut2bit1[buf[i]];
			i++;
		}
		m++;
	
		data[0][o]  = fp0[0];
		data[1][o]  = fp0[1];
		data[2][o]  = fp0[2];
		data[3][o]  = fp0[3];
		data[4][o]  = fp1[0];
		data[5][o]  = fp1[1];
		data[6][o]  = fp1[2];
		data[7][o]  = fp1[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_16track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 4;
			nblank++;
		}
		else
		{
			fp0 = lut2bit1[buf[i]];
			i++;
			fp1 = lut2bit1[buf[i]];
			i += 3;
		}
		m += 2;
	
		data[0][o]  = fp0[0];
		data[1][o]  = fp0[1];
		data[2][o]  = fp0[2];
		data[3][o]  = fp0[3];
		data[4][o]  = fp1[0];
		data[5][o]  = fp1[1];
		data[6][o]  = fp1[2];
		data[7][o]  = fp1[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_16track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;
	df = ms->decimation*2 - 1;
	df2 = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i++;
			nblank++;
		}
		else
		{
			fp0 = lut2bit1[buf[i]];
			i++;
			fp1 = lut2bit1[buf[i]];
		}
		i += df;
		m += df2;
	
		data[0][o]  = fp0[0];
		data[1][o]  = fp0[1];
		data[2][o]  = fp0[2];
		data[3][o]  = fp0[3];
		data[4][o]  = fp1[0];
		data[5][o]  = fp1[1];
		data[6][o]  = fp1[2];
		data[7][o]  = fp1[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_16track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 2;
			nblank++;
		}
		else
		{
			fp0 = lut2bit2[buf[i]];
			i++;
			fp1 = lut2bit2[buf[i]];
			i++;
		}
		m++;

		data[0][o] = fp0[0];
		data[1][o] = fp0[2];
		data[2][o] = fp1[0];
		data[3][o] = fp1[2];
		o++;
		data[0][o] = fp0[1];
		data[1][o] = fp0[3];
		data[2][o] = fp1[1];
		data[3][o] = fp1[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_16track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 2;
			nblank++;
		}
		else
		{
			fp0 = lut2bit2[buf[i]];
			i++;
			fp1 = lut2bit2[buf[i]];
			i++;
		}
		m++;

		data[0][o] = fp0[0];
		data[1][o] = fp0[2];
		data[2][o] = fp1[0];
		data[3][o] = fp1[2];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_16track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;
	df = ms->decimation - 1;
	df2 = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i++;
			nblank++;
		}
		else
		{
			fp0 = lut2bit2[buf[i]];
			i++;
			fp1 = lut2bit2[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o] = fp0[0];
		data[1][o] = fp0[2];
		data[2][o] = fp1[0];
		data[3][o] = fp1[2];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_16track_fanout4_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 2;
			nblank++;
		}
		else
		{
			fp0 = lut2bit3[buf[i]];
			i++;
			fp1 = lut2bit3[buf[i]];
			i++;
		}
		m++;

		data[0][o] = fp0[0];
		data[1][o] = fp1[0];
		o++;
		data[0][o] = fp0[1];
		data[1][o] = fp1[1];
		o++;
		data[0][o] = fp0[2];
		data[1][o] = fp1[2];
		o++;
		data[0][o] = fp0[3];
		data[1][o] = fp1[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_16track_fanout4_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i += 2;
			nblank++;
		}
		else
		{
			fp0 = lut2bit3[buf[i]];
			i++;
			fp1 = lut2bit3[buf[i]];
			i++;
		}
		m++;

		data[0][o] = fp0[0];
		data[1][o] = fp1[0];
		o++;
		data[0][o] = fp0[2];
		data[1][o] = fp1[2];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_16track_fanout4_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/2;
	df = ms->decimation/2 - 1;
	df2 = ms->decimation/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = zeros;
			i++;
			nblank++;
		}
		else
		{
			fp0 = lut2bit3[buf[i]];
			i++;
			fp1 = lut2bit3[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o] = fp0[0];
		data[1][o] = fp1[0];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_32track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 4;
			nblank++;
		}
		else
		{
			fp0 = lut2bit2[buf[i]];
			i++;
			fp1 = lut2bit2[buf[i]];
			i++;
			fp2 = lut2bit2[buf[i]];
			i++;
			fp3 = lut2bit2[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[2];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[2];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[2];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[3];
		data[10][o] = fp1[1];
		data[11][o] = fp1[3];
		data[12][o] = fp2[1];
		data[13][o] = fp2[3];
		data[14][o] = fp3[1];
		data[15][o] = fp3[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_32track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 8;
			nblank++;
		}
		else
		{
			fp0 = lut2bit2[buf[i]];
			i++;
			fp1 = lut2bit2[buf[i]];
			i++;
			fp2 = lut2bit2[buf[i]];
			i++;
			fp3 = lut2bit2[buf[i]];
			i += 5;
		}
		m += 2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[2];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[2];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[2];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[3];
		data[10][o] = fp1[1];
		data[11][o] = fp1[3];
		data[12][o] = fp2[1];
		data[13][o] = fp2[3];
		data[14][o] = fp3[1];
		data[15][o] = fp3[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_32track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;
	df = ms->decimation*4 - 3;
	df2 = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 3;
			nblank++;
		}
		else
		{
			fp0 = lut2bit2[buf[i]];
			i++;
			fp1 = lut2bit2[buf[i]];
			i++;
			fp2 = lut2bit2[buf[i]];
			i++;
			fp3 = lut2bit2[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[2];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[2];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[2];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[3];
		data[10][o] = fp1[1];
		data[11][o] = fp1[3];
		data[12][o] = fp2[1];
		data[13][o] = fp2[3];
		data[14][o] = fp3[1];
		data[15][o] = fp3[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_32track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 4;
			nblank++;
		}
		else
		{
			fp0 = lut2bit3[buf[i]];
			i++;
			fp1 = lut2bit3[buf[i]];
			i++;
			fp2 = lut2bit3[buf[i]];
			i++;
			fp3 = lut2bit3[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp1[0];
		data[2][o]  = fp2[0];
		data[3][o]  = fp3[0];
		data[4][o]  = fp0[1];
		data[5][o]  = fp1[1];
		data[6][o]  = fp2[1];
		data[7][o]  = fp3[1];
		o++;
		data[0][o]  = fp0[2];
		data[1][o]  = fp1[2];
		data[2][o]  = fp2[2];
		data[3][o]  = fp3[2];
		data[4][o]  = fp0[3];
		data[5][o]  = fp1[3];
		data[6][o]  = fp2[3];
		data[7][o]  = fp3[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_32track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 4;
			nblank++;
		}
		else
		{
			fp0 = lut2bit3[buf[i]];
			i++;
			fp1 = lut2bit3[buf[i]];
			i++;
			fp2 = lut2bit3[buf[i]];
			i++;
			fp3 = lut2bit3[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp1[0];
		data[2][o]  = fp2[0];
		data[3][o]  = fp3[0];
		data[4][o]  = fp0[1];
		data[5][o]  = fp1[1];
		data[6][o]  = fp2[1];
		data[7][o]  = fp3[1];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_32track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/4;
	df = ms->decimation*2 - 3;
	df2 = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			i += 3;
			nblank++;
		}
		else
		{
			fp0 = lut2bit3[buf[i]];
			i++;
			fp1 = lut2bit3[buf[i]];
			i++;
			fp2 = lut2bit3[buf[i]];
			i++;
			fp3 = lut2bit3[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp1[0];
		data[2][o]  = fp2[0];
		data[3][o]  = fp3[0];
		data[4][o]  = fp0[1];
		data[5][o]  = fp1[1];
		data[6][o]  = fp2[1];
		data[7][o]  = fp3[1];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_32track_fanout4_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned int *buf, bits;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i;
	unsigned char *bytes;
	int zone, l2;
	int nblank = 0;

	buf = (unsigned int *)(ms->payload);
	i = ms->readposition >> 2;  /* note here that i counts 32-bit words */
	l2 = ms->log2blankzonesize - 2;

	bytes = (unsigned char *)(& bits);

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> l2;

		if(4*i <  ms->blankzonestartvalid[zone] ||
		   4*i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			nblank++;
		}
		else
		{
			bits = reorder32(buf[i]);
			fp0 = lut2bit1[bytes[0]];
			fp1 = lut2bit1[bytes[1]];
			fp2 = lut2bit1[bytes[2]];
			fp3 = lut2bit1[bytes[3]];
		}
		i++;

		data[0][o] = fp0[0];
		data[1][o] = fp2[0];
		data[2][o] = fp1[0];
		data[3][o] = fp3[0];
		o++;
		data[0][o] = fp0[1];
		data[1][o] = fp2[1];
		data[2][o] = fp1[1];
		data[3][o] = fp3[1];
		o++;
		data[0][o] = fp0[2];
		data[1][o] = fp2[2];
		data[2][o] = fp1[2];
		data[3][o] = fp3[2];
		o++;
		data[0][o] = fp0[3];
		data[1][o] = fp2[3];
		data[2][o] = fp1[3];
		data[3][o] = fp3[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = (unsigned int *)(ms->payload);
			i = 0;
		}
	}

	ms->readposition = 4*i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_32track_fanout4_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned int *buf, bits;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i;
	unsigned char *bytes;
	int zone, l2;
	int nblank = 0;

	buf = (unsigned int *)(ms->payload);
	i = ms->readposition >> 2;  /* note here that i counts 32-bit words */
	l2 = ms->log2blankzonesize - 2;

	bytes = (unsigned char *)(& bits);

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> l2;

		if(4*i <  ms->blankzonestartvalid[zone] ||
		   4*i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			nblank++;
		}
		else
		{
			bits = reorder32(buf[i]);
			fp0 = lut2bit1[bytes[0]];
			fp1 = lut2bit1[bytes[1]];
			fp2 = lut2bit1[bytes[2]];
			fp3 = lut2bit1[bytes[3]];
		}
		i++;

		data[0][o] = fp0[0];
		data[1][o] = fp2[0];
		data[2][o] = fp1[0];
		data[3][o] = fp3[0];
		o++;
		data[0][o] = fp0[2];
		data[1][o] = fp2[2];
		data[2][o] = fp1[2];
		data[3][o] = fp3[2];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = (unsigned int *)(ms->payload);
			i = 0;
		}
	}

	ms->readposition = 4*i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_32track_fanout4_decimation4(struct mark5_stream *ms, int nsamp,
	float **data)
{
	unsigned int *buf, bits;
	float *fp0, *fp1, *fp2, *fp3;
	int o, i, df;
	unsigned char *bytes;
	int zone, l2;
	int nblank = 0;

	buf = (unsigned int *)(ms->payload);
	i = ms->readposition >> 2;  /* note here that i counts 32-bit words */
	l2 = ms->log2blankzonesize - 2;
	df = ms->decimation/4;

	bytes = (unsigned char *)(& bits);

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> l2;

		if(4*i <  ms->blankzonestartvalid[zone] ||
		   4*i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = zeros;
			nblank++;
		}
		else
		{
			bits = reorder32(buf[i]);
			fp0 = lut2bit1[bytes[0]];
			fp1 = lut2bit1[bytes[1]];
			fp2 = lut2bit1[bytes[2]];
			fp3 = lut2bit1[bytes[3]];
		}
		i += df;

		data[0][o] = fp0[0];
		data[1][o] = fp2[0];
		data[2][o] = fp1[0];
		data[3][o] = fp3[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = (unsigned int *)(ms->payload);
			i = 0;
		}
	}

	ms->readposition = 4*i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_64track_fanout1_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 8;
			nblank++;
		}
		else
		{
			fp0 = lut2bit2[buf[i]];
			i++;
			fp1 = lut2bit2[buf[i]];
			i++;
			fp2 = lut2bit2[buf[i]];
			i++;
			fp3 = lut2bit2[buf[i]];
			i++;
			fp4 = lut2bit2[buf[i]];
			i++;
			fp5 = lut2bit2[buf[i]];
			i++;
			fp6 = lut2bit2[buf[i]];
			i++;
			fp7 = lut2bit2[buf[i]];
			i++;
		}
		m++;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[2];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[2];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[2];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[3];
		data[10][o] = fp1[1];
		data[11][o] = fp1[3];
		data[12][o] = fp2[1];
		data[13][o] = fp2[3];
		data[14][o] = fp3[1];
		data[15][o] = fp3[3];
		data[16][o] = fp4[0];
		data[17][o] = fp4[2];
		data[18][o] = fp5[0];
		data[19][o] = fp5[2];
		data[20][o] = fp6[0];
		data[21][o] = fp6[2];
		data[22][o] = fp7[0];
		data[23][o] = fp7[2];
		data[24][o] = fp4[1];
		data[25][o] = fp4[3];
		data[26][o] = fp5[1];
		data[27][o] = fp5[3];
		data[28][o] = fp6[1];
		data[29][o] = fp6[3];
		data[30][o] = fp7[1];
		data[31][o] = fp7[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_64track_fanout1_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 16;
			nblank++;
		}
		else
		{
			fp0 = lut2bit2[buf[i]];
			i++;
			fp1 = lut2bit2[buf[i]];
			i++;
			fp2 = lut2bit2[buf[i]];
			i++;
			fp3 = lut2bit2[buf[i]];
			i++;
			fp4 = lut2bit2[buf[i]];
			i++;
			fp5 = lut2bit2[buf[i]];
			i++;
			fp6 = lut2bit2[buf[i]];
			i++;
			fp7 = lut2bit2[buf[i]];
			i += 9;
		}
		m += 2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[2];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[2];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[2];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[3];
		data[10][o] = fp1[1];
		data[11][o] = fp1[3];
		data[12][o] = fp2[1];
		data[13][o] = fp2[3];
		data[14][o] = fp3[1];
		data[15][o] = fp3[3];
		data[16][o] = fp4[0];
		data[17][o] = fp4[2];
		data[18][o] = fp5[0];
		data[19][o] = fp5[2];
		data[20][o] = fp6[0];
		data[21][o] = fp6[2];
		data[22][o] = fp7[0];
		data[23][o] = fp7[2];
		data[24][o] = fp4[1];
		data[25][o] = fp4[3];
		data[26][o] = fp5[1];
		data[27][o] = fp5[3];
		data[28][o] = fp6[1];
		data[29][o] = fp6[3];
		data[30][o] = fp7[1];
		data[31][o] = fp7[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_64track_fanout1_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;
	df = ms->decimation*8 - 7;
	df2 = ms->decimation;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 7;
			nblank++;
		}
		else
		{
			fp0 = lut2bit2[buf[i]];
			i++;
			fp1 = lut2bit2[buf[i]];
			i++;
			fp2 = lut2bit2[buf[i]];
			i++;
			fp3 = lut2bit2[buf[i]];
			i++;
			fp4 = lut2bit2[buf[i]];
			i++;
			fp5 = lut2bit2[buf[i]];
			i++;
			fp6 = lut2bit2[buf[i]];
			i++;
			fp7 = lut2bit2[buf[i]];
		}
		i += df;
		m += df2;

		data[0][o]  = fp0[0];
		data[1][o]  = fp0[2];
		data[2][o]  = fp1[0];
		data[3][o]  = fp1[2];
		data[4][o]  = fp2[0];
		data[5][o]  = fp2[2];
		data[6][o]  = fp3[0];
		data[7][o]  = fp3[2];
		data[8][o]  = fp0[1];
		data[9][o]  = fp0[3];
		data[10][o] = fp1[1];
		data[11][o] = fp1[3];
		data[12][o] = fp2[1];
		data[13][o] = fp2[3];
		data[14][o] = fp3[1];
		data[15][o] = fp3[3];
		data[16][o] = fp4[0];
		data[17][o] = fp4[2];
		data[18][o] = fp5[0];
		data[19][o] = fp5[2];
		data[20][o] = fp6[0];
		data[21][o] = fp6[2];
		data[22][o] = fp7[0];
		data[23][o] = fp7[2];
		data[24][o] = fp4[1];
		data[25][o] = fp4[3];
		data[26][o] = fp5[1];
		data[27][o] = fp5[3];
		data[28][o] = fp6[1];
		data[29][o] = fp6[3];
		data[30][o] = fp7[1];
		data[31][o] = fp7[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - nblank;
}

static int vlba_nomod_decode_2bit_64track_fanout2_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 8;
			nblank++;
		}
		else
		{
			fp0 = lut2bit3[buf[i]];
			i++;
			fp1 = lut2bit3[buf[i]];
			i++;
			fp2 = lut2bit3[buf[i]];
			i++;
			fp3 = lut2bit3[buf[i]];
			i++;
			fp4 = lut2bit3[buf[i]];
			i++;
			fp5 = lut2bit3[buf[i]];
			i++;
			fp6 = lut2bit3[buf[i]];
			i++;
			fp7 = lut2bit3[buf[i]];
			i++;
		}
		m++;
		
		data[0][o]  = fp0[0];
		data[1][o]  = fp1[0];
		data[2][o]  = fp2[0];
		data[3][o]  = fp3[0];
		data[4][o]  = fp0[1];
		data[5][o]  = fp1[1];
		data[6][o]  = fp2[1];
		data[7][o]  = fp3[1];
		data[8][o]  = fp4[0];
		data[9][o]  = fp5[0];
		data[10][o] = fp6[0];
		data[11][o] = fp7[0];
		data[12][o] = fp4[1];
		data[13][o] = fp5[1];
		data[14][o] = fp6[1];
		data[15][o] = fp7[1];
		o++;
		data[0][o]  = fp0[2];
		data[1][o]  = fp1[2];
		data[2][o]  = fp2[2];
		data[3][o]  = fp3[2];
		data[4][o]  = fp0[3];
		data[5][o]  = fp1[3];
		data[6][o]  = fp2[3];
		data[7][o]  = fp3[3];
		data[8][o]  = fp4[2];
		data[9][o]  = fp5[2];
		data[10][o] = fp6[2];
		data[11][o] = fp7[2];
		data[12][o] = fp4[3];
		data[13][o] = fp5[3];
		data[14][o] = fp6[3];
		data[15][o] = fp7[3];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_64track_fanout2_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 8;
			nblank++;
		}
		else
		{
			fp0 = lut2bit3[buf[i]];
			i++;
			fp1 = lut2bit3[buf[i]];
			i++;
			fp2 = lut2bit3[buf[i]];
			i++;
			fp3 = lut2bit3[buf[i]];
			i++;
			fp4 = lut2bit3[buf[i]];
			i++;
			fp5 = lut2bit3[buf[i]];
			i++;
			fp6 = lut2bit3[buf[i]];
			i++;
			fp7 = lut2bit3[buf[i]];
			i++;
		}
		m++;
		
		data[0][o]  = fp0[0];
		data[1][o]  = fp1[0];
		data[2][o]  = fp2[0];
		data[3][o]  = fp3[0];
		data[4][o]  = fp0[1];
		data[5][o]  = fp1[1];
		data[6][o]  = fp2[1];
		data[7][o]  = fp3[1];
		data[8][o]  = fp4[0];
		data[9][o]  = fp5[0];
		data[10][o] = fp6[0];
		data[11][o] = fp7[0];
		data[12][o] = fp4[1];
		data[13][o] = fp5[1];
		data[14][o] = fp6[1];
		data[15][o] = fp7[1];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_64track_fanout2_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned char *buf;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, m, df, df2;
	int zone;
	int nblank = 0;

	buf = ms->payload;
	i = ms->readposition;
	m = i/8;
	df = ms->decimation*4 - 7;
	df2 = ms->decimation/2;

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> ms->log2blankzonesize;

		if(i <  ms->blankzonestartvalid[zone] ||
		   i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			i += 7;
			nblank++;
		}
		else
		{
			fp0 = lut2bit3[buf[i]];
			i++;
			fp1 = lut2bit3[buf[i]];
			i++;
			fp2 = lut2bit3[buf[i]];
			i++;
			fp3 = lut2bit3[buf[i]];
			i++;
			fp4 = lut2bit3[buf[i]];
			i++;
			fp5 = lut2bit3[buf[i]];
			i++;
			fp6 = lut2bit3[buf[i]];
			i++;
			fp7 = lut2bit3[buf[i]];
		}
		i += df;
		m += df2;
		
		data[0][o]  = fp0[0];
		data[1][o]  = fp1[0];
		data[2][o]  = fp2[0];
		data[3][o]  = fp3[0];
		data[4][o]  = fp0[1];
		data[5][o]  = fp1[1];
		data[6][o]  = fp2[1];
		data[7][o]  = fp3[1];
		data[8][o]  = fp4[0];
		data[9][o]  = fp5[0];
		data[10][o] = fp6[0];
		data[11][o] = fp7[0];
		data[12][o] = fp4[1];
		data[13][o] = fp5[1];
		data[14][o] = fp6[1];
		data[15][o] = fp7[1];

		if(m >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = ms->payload;
			i = 0;
			m = 0;
		}
	}

	ms->readposition = i;

	return nsamp - 2*nblank;
}

static int vlba_nomod_decode_2bit_64track_fanout4_decimation1(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned long long *buf, bits;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i;
	unsigned char *bytes;
	int zone, l2;
	int nblank = 0;

	buf = (unsigned long long *)(ms->payload);
	i = ms->readposition >> 3;  /* note here that i counts 64-bit words */
	l2 = ms->log2blankzonesize - 3;

	bytes = (unsigned char *)(& bits);

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> l2;

		if(8*i <  ms->blankzonestartvalid[zone] ||
		   8*i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			nblank++;
		}
		else
		{
			bits = reorder64(buf[i]);
			fp0 = lut2bit1[bytes[0]];
			fp1 = lut2bit1[bytes[1]];
			fp2 = lut2bit1[bytes[2]];
			fp3 = lut2bit1[bytes[3]];
			fp4 = lut2bit1[bytes[4]];
			fp5 = lut2bit1[bytes[5]];
			fp6 = lut2bit1[bytes[6]];
			fp7 = lut2bit1[bytes[7]];
		}
		i++;
		
		data[0][o] = fp0[0];
		data[1][o] = fp2[0];
		data[2][o] = fp1[0];
		data[3][o] = fp3[0];
		data[4][o] = fp4[0];
		data[5][o] = fp6[0];
		data[6][o] = fp5[0];
		data[7][o] = fp7[0];
		o++;
		data[0][o] = fp0[1];
		data[1][o] = fp2[1];
		data[2][o] = fp1[1];
		data[3][o] = fp3[1];
		data[4][o] = fp4[1];
		data[5][o] = fp6[1];
		data[6][o] = fp5[1];
		data[7][o] = fp7[1];
		o++;
		data[0][o] = fp0[2];
		data[1][o] = fp2[2];
		data[2][o] = fp1[2];
		data[3][o] = fp3[2];
		data[4][o] = fp4[2];
		data[5][o] = fp6[2];
		data[6][o] = fp5[2];
		data[7][o] = fp7[2];
		o++;
		data[0][o] = fp0[3];
		data[1][o] = fp2[3];
		data[2][o] = fp1[3];
		data[3][o] = fp3[3];
		data[4][o] = fp4[3];
		data[5][o] = fp6[3];
		data[6][o] = fp5[3];
		data[7][o] = fp7[3];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = (unsigned long long *)(ms->payload);
			i = 0;
		}
	}

	ms->readposition = 8*i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_64track_fanout4_decimation2(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned long long *buf, bits;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i;
	unsigned char *bytes;
	int zone, l2;
	int nblank = 0;

	buf = (unsigned long long *)(ms->payload);
	i = ms->readposition >> 3;  /* note here that i counts 64-bit words */
	l2 = ms->log2blankzonesize - 3;

	bytes = (unsigned char *)(& bits);

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> l2;

		if(8*i <  ms->blankzonestartvalid[zone] ||
		   8*i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			nblank++;
		}
		else
		{
			bits = reorder64(buf[i]);
			fp0 = lut2bit1[bytes[0]];
			fp1 = lut2bit1[bytes[1]];
			fp2 = lut2bit1[bytes[2]];
			fp3 = lut2bit1[bytes[3]];
			fp4 = lut2bit1[bytes[4]];
			fp5 = lut2bit1[bytes[5]];
			fp6 = lut2bit1[bytes[6]];
			fp7 = lut2bit1[bytes[7]];
		}
		i++;
		
		data[0][o] = fp0[0];
		data[1][o] = fp2[0];
		data[2][o] = fp1[0];
		data[3][o] = fp3[0];
		data[4][o] = fp4[0];
		data[5][o] = fp6[0];
		data[6][o] = fp5[0];
		data[7][o] = fp7[0];
		o++;
		data[0][o] = fp0[2];
		data[1][o] = fp2[2];
		data[2][o] = fp1[2];
		data[3][o] = fp3[2];
		data[4][o] = fp4[2];
		data[5][o] = fp6[2];
		data[6][o] = fp5[2];
		data[7][o] = fp7[2];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = (unsigned long long *)(ms->payload);
			i = 0;
		}
	}

	ms->readposition = 8*i;

	return nsamp - 4*nblank;
}

static int vlba_nomod_decode_2bit_64track_fanout4_decimation4(struct mark5_stream *ms,
	int nsamp, float **data)
{
	unsigned long long *buf, bits;
	float *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7;
	int o, i, df;
	unsigned char *bytes;
	int zone, l2;
	int nblank = 0;

	buf = (unsigned long long *)(ms->payload);
	i = ms->readposition >> 3;  /* note here that i counts 64-bit words */
	l2 = ms->log2blankzonesize - 3;
	df = ms->decimation/4;

	bytes = (unsigned char *)(& bits);

	for(o = 0; o < nsamp; o++)
	{
		zone = i >> l2;

		if(8*i <  ms->blankzonestartvalid[zone] ||
		   8*i >= ms->blankzoneendvalid[zone])
		{
			fp0 = fp1 = fp2 = fp3 = fp4 = fp5 = fp6 = fp7 = zeros;
			nblank++;
		}
		else
		{
			bits = reorder64(buf[i]);
			fp0 = lut2bit1[bytes[0]];
			fp1 = lut2bit1[bytes[1]];
			fp2 = lut2bit1[bytes[2]];
			fp3 = lut2bit1[bytes[3]];
			fp4 = lut2bit1[bytes[4]];
			fp5 = lut2bit1[bytes[5]];
			fp6 = lut2bit1[bytes[6]];
			fp7 = lut2bit1[bytes[7]];
		}
		i += df;
		
		data[0][o] = fp0[0];
		data[1][o] = fp2[0];
		data[2][o] = fp1[0];
		data[3][o] = fp3[0];
		data[4][o] = fp4[0];
		data[5][o] = fp6[0];
		data[6][o] = fp5[0];
		data[7][o] = fp7[0];

		if(i >= PAYLOADSIZE)
		{
			if(mark5_stream_next_frame(ms) < 0)
			{
				return -1;
			}
			buf = (unsigned long long *)(ms->payload);
			i = 0;
		}
	}

	ms->readposition = 8*i;

	return nsamp - 4*nblank;
}

/******************************************************************/

static int mark5_format_vlba_nomod_make_formatname(struct mark5_stream *ms)
{
	struct mark5_format_vlba_nomod *f;

	f = (struct mark5_format_vlba_nomod *)(ms->formatdata);
	
	snprintf(ms->formatname, MARK5_STREAM_ID_LENGTH, "VLBN1_%d-%d-%d-%d/%d",
		f->fanout, ms->Mbps, ms->nchan, ms->nbit, ms->decimation);

	return 0;
}

static int mark5_format_vlba_nomod_init(struct mark5_stream *ms)
{
	struct mark5_format_vlba_nomod *f;
	int bytes;
	int tol=5;
	int mjd1, sec1, ns1;
	double dns1, dns;
	int datarate;
	int nRealTrack;

	if(!ms)
	{
		fprintf(m5stderr, "mark5_format_vlba_nomod_init: ms = 0\n");
		return -1;
	}

	f = (struct mark5_format_vlba_nomod *)(ms->formatdata);

	nRealTrack = f->ntrack;
	if(nRealTrack < 8)
	{
		nRealTrack = 8;
	}
	ms->samplegranularity = f->fanout/ms->decimation;
	if(ms->samplegranularity <= 0)
	{
		ms->samplegranularity = 1;
	}
	ms->framebytes = 20160*nRealTrack/8;
	ms->databytes = 20000*nRealTrack/8;
	ms->payloadoffset = 12*nRealTrack;
	ms->framesamples = 20000*f->fanout/ms->decimation;
	ms->format = MK5_FORMAT_VLBN;
	ms->framegranularity = 1;

	ms->blanker = blanker_mark5;

	if(ms->datawindow)
	{
		if(ms->datawindowsize < ms->framebytes)
		{
			return -1;
		}

		/* look through entire data window, up to 1Mibytes */
		bytes = ms->datawindowsize < (1<<20) ?
			ms->datawindowsize : (1<<20);
		ms->frameoffset = findfirstframe(ms->datawindow, bytes,
			nRealTrack, tol);
		if(ms->frameoffset < 0)
		{
			return -1;
		}

		ms->frame = ms->datawindow + ms->frameoffset;
		ms->payload = ms->frame + ms->payloadoffset;

		ms->gettime(ms, &ms->mjd, &ms->sec, &dns);
		ms->ns = (int)(dns+0.5);
		ms->frame += ms->framebytes;
		ms->gettime(ms, &mjd1, &sec1, &dns1);
		ns1 = (int)(dns1 + 0.5);
		ms->frame -= ms->framebytes;

		/* assume frame time less than 1 second, integer number of
		 * frames per second
		 */
		if(ns1 != ms->ns)
		{
			ms->framens = ns1 - ms->ns;
			if(ms->framens <= 0)
			{
				ms->framens += 1000000000;
			}
			ms->samprate = ms->framesamples*
				(1000000000/ms->framens);
			datarate = ms->samprate*ms->nbit*ms->nchan/1000000;
			if(datarate != ms->Mbps)
			{
				if(ms->Mbps > 0)
				{
					fprintf(m5stderr, "Warning: data rate "
						"disagrees : %d != %d\n",
						datarate, ms->Mbps);
				}
				ms->Mbps = datarate;
			}
		}
		else
		{
			fprintf(m5stderr, "Warning: rate calc. suspect\n");
		}
	}

	ms->gframens = (int)(ms->framegranularity*ms->framens + 0.5);

	mark5_format_vlba_nomod_make_formatname(ms);

	return 0;
}

static int mark5_format_vlba_nomod_final(struct mark5_stream *ms)
{
	if(!ms)
	{
		return -1;
	}

	if(ms->formatdata)
	{
		free(ms->formatdata);
	}

	return 0;
}

struct mark5_format_generic *new_mark5_format_vlba_nomod(int Mbps, int nchan,
	int nbit, int fanout, int decimation)
{
	static int first = 1;
	struct mark5_format_generic *f;
	struct mark5_format_vlba_nomod *v;
	int decoderindex=0;
	int ntrack;

	ntrack = nchan*fanout*nbit;

	if(first)
	{
		first = 0;
		initluts();
	}

	if(decimation == 1)
	{
		decoderindex += 0;
	}
	else if(decimation == 2)
	{
		decoderindex += 42;
	}
	else if(decimation % 4 == 0) /* all multiples of 4 are handed here */
	{
		decoderindex += 84;
	}
	else
	{
		fprintf(m5stderr, "decimation must be 1, 2 or a mult of 4\n");
	}

	if(nbit == 1)
	{
		decoderindex += 0;
	}
	else if(nbit == 2)
	{
		decoderindex += 21;
	}
	else
	{
		fprintf(m5stderr, "nbit must be 1 or 2\n");
		return 0;
	}

	if(ntrack == 1)
	{
		decoderindex += 0;
	}
	else if(ntrack == 2)
	{
		decoderindex += 3;
	}
	else if(ntrack == 4)
	{
		decoderindex += 6;
	}
	else if(ntrack == 8)
	{
		decoderindex += 9;
	}
	else if(ntrack == 16)
	{
		decoderindex += 12;
	}
	else if(ntrack == 32)
	{
		decoderindex += 15;
	}
	else if(ntrack == 64)
	{
		decoderindex += 18;
	}
	else
	{
		fprintf(m5stderr, "ntrack must be 2^n : n = 0..6\n");
		return 0;
	}

	if(fanout == 1)
	{
		decoderindex += 0;
	}
	else if(fanout == 2)
	{
		decoderindex += 1;
	}
	else if(fanout == 4)
	{
		decoderindex += 2;
	}
	else
	{
		fprintf(m5stderr, "fanout must be 1, 2 or 4\n");
		return 0;
	}

	v = (struct mark5_format_vlba_nomod *)calloc(1,
		sizeof(struct mark5_format_vlba_nomod));
	f = (struct mark5_format_generic *)calloc(1,
		sizeof(struct mark5_format_generic));

	v->ntrack = ntrack;
	v->fanout = fanout;
	v->kday = 0;

	f->Mbps = Mbps;
	f->nchan = nchan;
	f->nbit = nbit;
	f->formatdata = v;
	f->formatdatasize = sizeof(struct mark5_format_vlba_nomod);
	f->gettime = mark5_format_vlba_nomod_frame_time;
	f->fixmjd = mark5_format_vlba_nomod_fixmjd;
	f->init_format = mark5_format_vlba_nomod_init;
	f->final_format = mark5_format_vlba_nomod_final;
	f->validate = mark5_format_vlba_nomod_validate;
	f->decimation = decimation;
	f->complex_decode = 0;
	f->count = 0;
	switch(decoderindex)
	{
		case 0  : f->decode = vlba_nomod_decode_1bit_1track_fanout1_decimation1; break;


		case 3  : f->decode = vlba_nomod_decode_1bit_2track_fanout1_decimation1; break;
		case 4  : f->decode = vlba_nomod_decode_1bit_2track_fanout2_decimation1; break;
		
		case 6  : f->decode = vlba_nomod_decode_1bit_4track_fanout1_decimation1; break;
		case 7  : f->decode = vlba_nomod_decode_1bit_4track_fanout2_decimation1; break;
		case 8  : f->decode = vlba_nomod_decode_1bit_4track_fanout4_decimation1; break;
		case 9  : f->decode = vlba_nomod_decode_1bit_8track_fanout1_decimation1; break;
		case 10 : f->decode = vlba_nomod_decode_1bit_8track_fanout2_decimation1; break;
		case 11 : f->decode = vlba_nomod_decode_1bit_8track_fanout4_decimation1; break;
		case 12 : f->decode = vlba_nomod_decode_1bit_16track_fanout1_decimation1; break;
		case 13 : f->decode = vlba_nomod_decode_1bit_16track_fanout2_decimation1; break;
		case 14 : f->decode = vlba_nomod_decode_1bit_16track_fanout4_decimation1; break;
		case 15 : f->decode = vlba_nomod_decode_1bit_32track_fanout1_decimation1; break;
		case 16 : f->decode = vlba_nomod_decode_1bit_32track_fanout2_decimation1; break;
		case 17 : f->decode = vlba_nomod_decode_1bit_32track_fanout4_decimation1; break;
		case 18 : f->decode = vlba_nomod_decode_1bit_64track_fanout1_decimation1; break;
		case 19 : f->decode = vlba_nomod_decode_1bit_64track_fanout2_decimation1; break;
		case 20 : f->decode = vlba_nomod_decode_1bit_64track_fanout4_decimation1; break;
		
		
		
		case 24 : f->decode = vlba_nomod_decode_2bit_2track_fanout1_decimation1; break;
		
		
		case 27 : f->decode = vlba_nomod_decode_2bit_4track_fanout1_decimation1; break;
		case 28 : f->decode = vlba_nomod_decode_2bit_4track_fanout2_decimation1; break;
		
		case 30 : f->decode = vlba_nomod_decode_2bit_8track_fanout1_decimation1; break;
		case 31 : f->decode = vlba_nomod_decode_2bit_8track_fanout2_decimation1; break;
		case 32 : f->decode = vlba_nomod_decode_2bit_8track_fanout4_decimation1; break;
		case 33 : f->decode = vlba_nomod_decode_2bit_16track_fanout1_decimation1; break;
		case 34 : f->decode = vlba_nomod_decode_2bit_16track_fanout2_decimation1; break;
		case 35 : f->decode = vlba_nomod_decode_2bit_16track_fanout4_decimation1; break;
		case 36 : f->decode = vlba_nomod_decode_2bit_32track_fanout1_decimation1; break;
		case 37 : f->decode = vlba_nomod_decode_2bit_32track_fanout2_decimation1; break;
		case 38 : f->decode = vlba_nomod_decode_2bit_32track_fanout4_decimation1; break;
		case 39 : f->decode = vlba_nomod_decode_2bit_64track_fanout1_decimation1; break;
		case 40 : f->decode = vlba_nomod_decode_2bit_64track_fanout2_decimation1; break;
		case 41 : f->decode = vlba_nomod_decode_2bit_64track_fanout4_decimation1; break;
		case 42 : f->decode = vlba_nomod_decode_1bit_1track_fanout1_decimation2; break;
		
		
		case 45 : f->decode = vlba_nomod_decode_1bit_2track_fanout1_decimation2; break;
		case 46 : f->decode = vlba_nomod_decode_1bit_2track_fanout2_decimation2; break;
		
		case 48 : f->decode = vlba_nomod_decode_1bit_4track_fanout1_decimation2; break;
		case 49 : f->decode = vlba_nomod_decode_1bit_4track_fanout2_decimation2; break;
		case 50 : f->decode = vlba_nomod_decode_1bit_4track_fanout4_decimation2; break;
		case 51 : f->decode = vlba_nomod_decode_1bit_8track_fanout1_decimation2; break;
		case 52 : f->decode = vlba_nomod_decode_1bit_8track_fanout2_decimation2; break;
		case 53 : f->decode = vlba_nomod_decode_1bit_8track_fanout4_decimation2; break;
		case 54 : f->decode = vlba_nomod_decode_1bit_16track_fanout1_decimation2; break;
		case 55 : f->decode = vlba_nomod_decode_1bit_16track_fanout2_decimation2; break;
		case 56 : f->decode = vlba_nomod_decode_1bit_16track_fanout4_decimation2; break;
		case 57 : f->decode = vlba_nomod_decode_1bit_32track_fanout1_decimation2; break;
		case 58 : f->decode = vlba_nomod_decode_1bit_32track_fanout2_decimation2; break;
		case 59 : f->decode = vlba_nomod_decode_1bit_32track_fanout4_decimation2; break;
		case 60 : f->decode = vlba_nomod_decode_1bit_64track_fanout1_decimation2; break;
		case 61 : f->decode = vlba_nomod_decode_1bit_64track_fanout2_decimation2; break;
		case 62 : f->decode = vlba_nomod_decode_1bit_64track_fanout4_decimation2; break;
		
		
		
		case 66 : f->decode = vlba_nomod_decode_2bit_2track_fanout1_decimation2; break;
		
		
		case 69 : f->decode = vlba_nomod_decode_2bit_4track_fanout1_decimation2; break;
		case 70 : f->decode = vlba_nomod_decode_2bit_4track_fanout2_decimation2; break;
		
		case 72 : f->decode = vlba_nomod_decode_2bit_8track_fanout1_decimation2; break;
		case 73 : f->decode = vlba_nomod_decode_2bit_8track_fanout2_decimation2; break;
		case 74 : f->decode = vlba_nomod_decode_2bit_8track_fanout4_decimation2; break;
		case 75 : f->decode = vlba_nomod_decode_2bit_16track_fanout1_decimation2; break;
		case 76 : f->decode = vlba_nomod_decode_2bit_16track_fanout2_decimation2; break;
		case 77 : f->decode = vlba_nomod_decode_2bit_16track_fanout4_decimation2; break;
		case 78 : f->decode = vlba_nomod_decode_2bit_32track_fanout1_decimation2; break;
		case 79 : f->decode = vlba_nomod_decode_2bit_32track_fanout2_decimation2; break;
		case 80 : f->decode = vlba_nomod_decode_2bit_32track_fanout4_decimation2; break;
		case 81 : f->decode = vlba_nomod_decode_2bit_64track_fanout1_decimation2; break;
		case 82 : f->decode = vlba_nomod_decode_2bit_64track_fanout2_decimation2; break;
		case 83 : f->decode = vlba_nomod_decode_2bit_64track_fanout4_decimation2; break;
		case 84 : f->decode = vlba_nomod_decode_1bit_1track_fanout1_decimation4; break;
		
		
		case 87 : f->decode = vlba_nomod_decode_1bit_2track_fanout1_decimation4; break;
		case 88 : f->decode = vlba_nomod_decode_1bit_2track_fanout2_decimation4; break;
		
		case 90 : f->decode = vlba_nomod_decode_1bit_4track_fanout1_decimation4; break;
		case 91 : f->decode = vlba_nomod_decode_1bit_4track_fanout2_decimation4; break;
		case 92 : f->decode = vlba_nomod_decode_1bit_4track_fanout4_decimation4; break;
		case 93 : f->decode = vlba_nomod_decode_1bit_8track_fanout1_decimation4; break;
		case 94 : f->decode = vlba_nomod_decode_1bit_8track_fanout2_decimation4; break;
		case 95 : f->decode = vlba_nomod_decode_1bit_8track_fanout4_decimation4; break;
		case 96 : f->decode = vlba_nomod_decode_1bit_16track_fanout1_decimation4; break;
		case 97 : f->decode = vlba_nomod_decode_1bit_16track_fanout2_decimation4; break;
		case 98 : f->decode = vlba_nomod_decode_1bit_16track_fanout4_decimation4; break;
		case 99 : f->decode = vlba_nomod_decode_1bit_32track_fanout1_decimation4; break;
		case 100: f->decode = vlba_nomod_decode_1bit_32track_fanout2_decimation4; break;
		case 101: f->decode = vlba_nomod_decode_1bit_32track_fanout4_decimation4; break;
		case 102: f->decode = vlba_nomod_decode_1bit_64track_fanout1_decimation4; break;
		case 103: f->decode = vlba_nomod_decode_1bit_64track_fanout2_decimation4; break;
		case 104: f->decode = vlba_nomod_decode_1bit_64track_fanout4_decimation4; break;
		
		
		
		case 108: f->decode = vlba_nomod_decode_2bit_2track_fanout1_decimation4; break;
		
		
		case 111: f->decode = vlba_nomod_decode_2bit_4track_fanout1_decimation4; break;
		case 112: f->decode = vlba_nomod_decode_2bit_4track_fanout2_decimation4; break;
		
		case 114: f->decode = vlba_nomod_decode_2bit_8track_fanout1_decimation4; break;
		case 115: f->decode = vlba_nomod_decode_2bit_8track_fanout2_decimation4; break;
		case 116: f->decode = vlba_nomod_decode_2bit_8track_fanout4_decimation4; break;
		case 117: f->decode = vlba_nomod_decode_2bit_16track_fanout1_decimation4; break;
		case 118: f->decode = vlba_nomod_decode_2bit_16track_fanout2_decimation4; break;
		case 119: f->decode = vlba_nomod_decode_2bit_16track_fanout4_decimation4; break;
		case 120: f->decode = vlba_nomod_decode_2bit_32track_fanout1_decimation4; break;
		case 121: f->decode = vlba_nomod_decode_2bit_32track_fanout2_decimation4; break;
		case 122: f->decode = vlba_nomod_decode_2bit_32track_fanout4_decimation4; break;
		case 123: f->decode = vlba_nomod_decode_2bit_64track_fanout1_decimation4; break;
		case 124: f->decode = vlba_nomod_decode_2bit_64track_fanout2_decimation4; break;
		case 125: f->decode = vlba_nomod_decode_2bit_64track_fanout4_decimation4; break;
		default:  f->decode = 0;
	}

	if(f->decode == 0)
	{
		fprintf(m5stderr, "Illegal combination of fanout, tracks and bits\n");
		free(f);
		free(v);
		return 0;
	}

	return f;
}
