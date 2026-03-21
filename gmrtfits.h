#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <float.h>
#include <string.h>
#include "fitsio.h"

// type def
/*typedef uint16_t u16_t;*/
/*typedef uint8_t u8_t;*/
typedef float real_t;


typedef struct {
	fitsfile *fits; 
	int status;

	double mjd;
	double tsamp;

	// RA and DEC in degrees
	float rad, decd;

	// 
	float bandwidth_mhz;
	unsigned fftint;

	unsigned nchan, ngulp, npol;
	unsigned nchan_npol;
	unsigned ngulp_nchan_npol;
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

} gmrtfits_t;

void gmrtfits_prepare ( gmrtfits_t *fits, const char *filename, unsigned npol, unsigned nchan, float fedge, float bw, unsigned ngulp, unsigned fftint, unsigned bitdepth );

// public methods
void gmrtfits_create ( gmrtfits_t *fits );
void gmrtfits_close ( gmrtfits_t *fits );

void gmrtfits_subint_open ( gmrtfits_t *fits );

/*void add_subint_16bit ( gmrtfits_t *fits, u16_t *data );*/
/*void add_subint_8bit ( gmrtfits_t *fits, u8_t *data );*/
void gmrtfits_subint_real_single ( gmrtfits_t *fits, real_t *data );

// private methods

/*void normalize_data ( gmrtfits_t *fits, u16_t *data );*/
void normalize_data ( gmrtfits_t *fits, real_t *data );

