#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"

#include <signal.h>

// type def
/*typedef uint16_t u16_t;*/
/*typedef uint8_t u8_t;*/
typedef float real_t;
typedef float deal_t;
/*typedef double deal_t;*/

typedef enum { SEARCH, FOLD } fitsmode_t;


typedef struct {
	fitsfile *fits; 
	int status;

	double mjd;
	/* tsamp or tbin */
	double tsamp_tbin;
	double period;
	float dm;

	// RA and DEC in degrees
	float rad, decd;

	// 
	float bandwidth_mhz;
	unsigned fftint;

	unsigned nchan, nsblk, npol;
	unsigned nchan_npol;
	unsigned nsblk_nchan_npol;

	unsigned nbin;
	unsigned nbin_nchan_npol;
	unsigned bitdepth;

	float *center_freqs;

	// not more than that
	char source[16];
	char pol_type[16];

	char *filename;

	unsigned long nrow;

	// to manage IO
	float *weights, *scales, *offsets;
	char  *data;
	int16_t *bata;
	/*float *bata;*/

	// mode
	fitsmode_t obsmode;
} gmrtfits_t;

// XXX make it like FFTW
// returns pointer to gmrtfits_t
void gmrtfits_search_prepare ( gmrtfits_t *fits, const char *filename, double mjd, unsigned npol, unsigned nchan, float fedge, float bw, unsigned nsblk, unsigned fftint, unsigned bitdepth );
void gmrtfits_fold_prepare ( gmrtfits_t *fits, const char *filename, double mjd, unsigned npol, unsigned nchan, float fedge, float bw, unsigned nbin, double period );

/* IO methods */
void gmrtfits_open ( gmrtfits_t *fits );
void gmrtfits_data_table ( gmrtfits_t *fits );
void gmrtfits_close ( gmrtfits_t *fits );

/* *
 *
 * */
void gmrtfits_search_add ( gmrtfits_t *fits, real_t *data, unsigned int ngulp );
void gmrtfits_fold_add ( gmrtfits_t *fits, deal_t *data, unsigned nturns, double period );

// private methods

/*void normalize_data ( gmrtfits_t *fits, u16_t *data );*/
/*void normalize_data ( gmrtfits_t *fits, real_t *data );*/

