/***************************************************************************
 *   Copyright (C) 2006-2011 by Walter Brisken                             *
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
// $Id: mark5_stream.h,v 1.4 2012/02/21 13:35:44 jwagnerhki Exp $
// $HeadURL: https://svn.atnf.csiro.au/difx/libraries/mark5access/trunk/mark5access/mark5_stream.h $
// $LastChangedRevision: 3356 $
// $Author: jwagnerhki $
// $LastChangedDate: 2011-06-02 18:21:45 +0200 (Thu, 02 Jun 2011) $
//
//============================================================================

#ifndef __MARK5_STREAM_H__
#define __MARK5_STREAM_H__

/* needed for open64 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

// Pinched from Linux features.h, as OSX does not have it
#if defined __GNUC__ && defined __GNUC_MINOR__
# define __GNUC_PREREQ(maj, min) \
        ((__GNUC__ << 16) + __GNUC_MINOR__ >= ((maj) << 16) + (min))
#else
# define __GNUC_PREREQ(maj, min) 0
#endif

#if !defined(__cplusplus) || !defined(__GNUC_PREREQ) || __GNUC_PREREQ(4,3)
#include <complex.h>
typedef double complex mark5_double_complex;
typedef float  complex mark5_float_complex;
#else
typedef struct { double re, im; } mark5_double_complex;
typedef struct { float  re, im; } mark5_float_complex;
#endif

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

enum Mark5Format
{
	MK5_FORMAT_UNKNOWN = -1,
	MK5_FORMAT_VLBA    =  0,
	MK5_FORMAT_MARK4   =  1,
	MK5_FORMAT_MARK5B  =  2,
	MK5_FORMAT_VDIF    =  3,
	MK5_FORMAT_VDIFL   =  4,		/* Legacy headers on VDIF */
	MK5_FORMAT_K5      =  5,		/* Not Yet Implemented */
	MK5_FORMAT_VLBN    =  6
};

#define MAXBLANKZONES		32
#define OPTIMAL_2BIT_HIGH	3.3359
#define MARK5_STREAM_ID_LENGTH	256

enum Mark5Blanker
{
	MK5_BLANKER_NONE  = 0,
	MK5_BLANKER_MARK5 = 1,
	MK5_BLANKER_VDIF  = 2
};

struct mark5_stream
{
	/* globally readable values: should not be changed */
	char streamname[MARK5_STREAM_ID_LENGTH]; /* name of stream */
	char formatname[MARK5_STREAM_ID_LENGTH]; /* name of format */
	enum Mark5Format format;/* format id */
	int Mbps;		/* total data rate */
	int nchan;		/* # of data channels; all will be decoded */
	int nbit;		/* quantization bits of data */
	int samplegranularity;	/* decoding and copying must be in mults of */
	int framegranularity;	/* min num of frames to have int. ns length */
	int mjd;		/* date of first found frame */
	int sec;		/* time of first found frame */
	int ns;			/* ns portion of time of first frame */
	int samprate;		/* (Hz) of de-fanned stream */
	int frameoffset;	/* bytes into stream of first frame */
	int framesamples;	/* number of samples per chan in a frame */
	double framens;		/* nanoseconds per frame */
	int gframens;		/* integer ns for framegranularity frames */
	int framebytes;		/* total number of bytes in a frame */
	int databytes;		/* bytes of data in a frame, incl. data */
				/*   replacement headers */
	long long framenum;	/* current complete frame, start at 0 */
	int decimation;		/* decimation factor */
	int nvalidatefail;	/* number of times frame validation failed */
	int nvalidatepass;	/* number of times frame validation passed */
	int consecutivefails;	/* number of validations failed in a row */

	/* internal state parameters: not to be used by users */
	unsigned char *frame;
	unsigned char *payload;
	int payloadoffset;	    /* == payload - frame */
	long long datawindowsize;     /* number of bytes resident at a time */
	unsigned char *datawindow;	    /* pointer to data window */
	int readposition;	    /* index into frame of current read */

	/* data blanking */
	int log2blankzonesize;
	int blankzonestartvalid[MAXBLANKZONES];
	int blankzoneendvalid[MAXBLANKZONES];
	int (*blanker)(struct mark5_stream *ms);

	/* stream commands and data pointer */
	int (*init_stream)(struct mark5_stream *ms);
	int (*final_stream)(struct mark5_stream *ms);
	int (*next)(struct mark5_stream *ms);
	int (*seek)(struct mark5_stream *ms, long long framenum);
	void *inputdata;

	/* format commands and data pointer */
	int (*init_format)(struct mark5_stream *ms);
	int (*final_format)(struct mark5_stream *ms);
	int (*decode)(struct mark5_stream *ms, int nsamp, float **data);
	int (*count)(struct mark5_stream *ms, int nsamp, unsigned int *highstates);
        int (*complex_decode)(struct mark5_stream *ms, int nsamp, mark5_float_complex **data);
	int (*validate)(const struct mark5_stream *ms);
	int (*gettime)(const struct mark5_stream *ms, int *mjd, 
		int *sec, double *ns);
	int (*fixmjd)(struct mark5_stream *ms, int refmjd);
	void *formatdata;
};

struct mark5_stream_generic
{
	int (*init_stream)(struct mark5_stream *ms);	/* required */
	int (*final_stream)(struct mark5_stream *ms);	/* required */
	int (*next)(struct mark5_stream *ms);		/* required */
	int (*seek)(struct mark5_stream *ms, long long framenum);
	void *inputdata;
	int inputdatasize;
};

struct mark5_format_generic
{
	int (*init_format)(struct mark5_stream *ms);	/* required */
	int (*final_format)(struct mark5_stream *ms);	/* required */
	int (*decode)(struct mark5_stream *ms, 		/* required */
		int nsamp, float **data); 
	int (*count)(struct mark5_stream *ms,
		int nsamp, unsigned int *highstates); 
	int (*complex_decode)(struct mark5_stream *ms,
		int nsamp, mark5_float_complex **data); 
	int (*validate)(const struct mark5_stream *ms);
	int (*gettime)(const struct mark5_stream *ms, 	/* required */
		int *mjd, int *sec, double *ns);
	int (*fixmjd)(struct mark5_stream *ms, int refmjd);
	void *formatdata;
	int formatdatasize;
	int Mbps;
	int nchan;
	int nbit;
	int decimation;					/* decimationling factor */
};

void delete_mark5_stream_generic(struct mark5_stream_generic *s);

void delete_mark5_format_generic(struct mark5_format_generic *f);

struct mark5_stream *new_mark5_stream(const struct mark5_stream_generic *s,
	const struct mark5_format_generic *f);

/* same as new_mark5_stream() but deletes s and f upon creation */
struct mark5_stream *new_mark5_stream_absorb(struct mark5_stream_generic *s,
	struct mark5_format_generic *f);

void delete_mark5_stream(struct mark5_stream *ms);

int mark5_stream_print(const struct mark5_stream *ms);

int mark5_stream_get_frame_time(struct mark5_stream *ms, 
	int *mjd, int *sec, double *ns);

int mark5_stream_get_sample_time(struct mark5_stream *ms, 
	int *mjd, int *sec, double *ns);

int mark5_stream_fix_mjd(struct mark5_stream *ms, int refmjd);

int mark5_stream_seek(struct mark5_stream *ms, int mjd, int sec, double ns);

int mark5_stream_copy(struct mark5_stream *ms, int nbytes, char *data);

int mark5_stream_set_blanker(struct mark5_stream *ms,
	enum Mark5Blanker blanker);

int mark5_stream_decode(struct mark5_stream *ms, int nsamp, float **data);

int mark5_stream_decode_double(struct mark5_stream *ms, int nsamp, 
	double **data);

int mark5_stream_decode_complex(struct mark5_stream *ms, int nsamp, 
	mark5_float_complex **data);

int mark5_stream_decode_double_complex(struct mark5_stream *ms, int nsamp, 
	mark5_double_complex **data);

int mark5_stream_count_high_states(struct mark5_stream *ms, int nsamp,
	unsigned int *highstates);

/* SPECIFIC STREAM TYPES */

/*   Memory based stream */

struct mark5_stream_generic *new_mark5_stream_memory(void *data,
	unsigned int nbytes);

/*   File based stream */

struct mark5_stream_generic *new_mark5_stream_file(const char *filename,
	long long offset);

int mark5_stream_file_add_infile(struct mark5_stream *ms, 
	const char *filename);

/*   Just an unpacker: for repeated unpacking of a particular format from
 *	arbitrary memory locations 
 */

struct mark5_stream_generic *new_mark5_stream_unpacker(int noheaders);

int mark5_unpack(struct mark5_stream *ms, void *packed, float **unpacked,
	int nsamp);

int mark5_unpack_with_offset(struct mark5_stream *ms, void *packed,
	int offsetsamples, float **unpacked, int nsamp);

int mark5_unpack_complex(struct mark5_stream *ms, void *packed, 
			 mark5_float_complex **unpacked, int nsamp);

int mark5_unpack_complex_with_offset(struct mark5_stream *ms, void *packed,
	int offsetsamples, mark5_float_complex **unpacked, int nsamp);


/* SPECIFIC FORMAT TYPES */

/*   VLBA format */

struct mark5_format_generic *new_mark5_format_vlba(int Mbps, int nchan,
	int nbit, int fanout, int decimation);

struct mark5_format_generic *new_mark5_format_vlba_nomod(int Mbps, int nchan,
	int nbit, int fanout, int decimation);

/*   Mark4 format */

struct mark5_format_generic *new_mark5_format_mark4(int Mbps, int nchan,
	int nbit, int fanout, int decimation);

/*   Mark5B format */

struct mark5_format_generic *new_mark5_format_mark5b(int Mbps, 
	int nchan, int nbit, int decimation);

/*   VDIF format */

struct mark5_format_generic *new_mark5_format_vdif(int Mbps,
	int nchan, int nbit, int decimation, 
	int databytesperpacket, int frameheadersize, int usecomplex);

void mark5_format_vdif_set_leapsecs(struct mark5_stream *ms, int leapsecs);

/*   K5 format: not yet complete */

struct mark5_format_generic *new_mark5_format_k5(int Mbps, int nchan, int nbit,
	int submode);

/*   Generate format from a string description */

struct mark5_format_generic *new_mark5_format_generic_from_string(
	const char *formatname);


/* DATA BLANKING ALGORITHMS */

/* The null blanker */

int blanker_none(struct mark5_stream *ms);

/* Blanker for generic data stored on Mark5 modules */

int blanker_mark5(struct mark5_stream *ms);

/* Blanker for Mark4 data stored on Mark5 modules */

int blanker_mark4(struct mark5_stream *ms);

/* Blankser for VDIF data stored on Mark5 modules */

int blanker_vdif(struct mark5_stream *ms);

/* TO PARTIALLY DETERMINE DATA FORMAT FROM DATA OR DESCRIPTION */

/* contains information that can be determined by a glance at data or name */
struct mark5_format
{
	enum Mark5Format format;  /* format type */
	int Mbps, nchan, nbit;
	int frameoffset;	  /* bytes from stream start to 1st frame */
	int framebytes;		  /* bytes in a frame */
	int databytes;		  /* bytes of data in a frame */
	double framens;		  /* duration of a frame in nanosec */
	int mjd, sec, ns;	  /* date and time of first frame */
	int ntrack;		  /* for Mark4 and VLBA formats only */
	int fanout;		  /* for Mark4 and VLBA formats only */
	int decimation;
};

const char *mark5_stream_list_formats();

struct mark5_format *new_mark5_format_from_name(const char *formatname);

struct mark5_format *new_mark5_format_from_stream(
	struct mark5_stream_generic *s);

void delete_mark5_format(struct mark5_format *mf);

void print_mark5_format(const struct mark5_format *mf);


/* OTHER USEFUL FUNCTIONS */
double correct_2bit_power(double x);
double high_state_fraction_to_power(double x);


/* BELOW HERE USER BEWARE: If you are using any functionality defined
 * below this comment then your code may not function properly with
 * future versions of this library.
 */

void mark5_library_init(void);
int mark5_library_getoption(const int mk5option, void* result);
int mark5_library_setoption(const int mk5option, void* value);

#define M5A_OPT_STDOUTFD 1
#define M5A_OPT_STDERRFD 2
extern FILE* m5stderr;
extern FILE* m5stdout;

/* private functions, not intended for external use */

int mark5_stream_next_frame(struct mark5_stream *ms);


/* for compatibility */

struct mark5_stream *mark5_stream_open(const char *filename,
        int nbit, int fanout, long long offset);

#define PAYLOADSIZE 20000
#define FRAMESIZE   20160

#define VLBA_FRAMESIZE	20160
#define MARK4_FRAMESIZE	20000

#ifdef __cplusplus
}
#endif

#endif
