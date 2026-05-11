#include "gmrtfits.h"

#ifdef TIMING
#include <time.h>
#endif

#define FLUSH_TIMES 8
// flush once per FLUSH_TIMES rows

void gmrtfits_search_prepare ( gmrtfits_t *fits, const char *filename, double mjd, unsigned npol, unsigned nchan, float fedge, float bw, unsigned nsblk, unsigned fftint, unsigned bitdepth ) {
	// sign(bw) determines LSB or USB
	
	fits->obsmode  = SEARCH;
	
	fits->filename = strdup ( filename );
	
	fits->status  = 0;
	fits->nbin    = 1;
	fits->npol    = npol;
	fits->nchan   = nchan;
	fits->nsblk   = nsblk;

	if ( fits->npol == 4 ) {
		sprintf ( fits->pol_type, "AABBCRCI" );
	} else if ( fits->npol == 1 ) {
		sprintf ( fits->pol_type, "AA+BB" );
	}


	fits->nchan_npol       = nchan * npol;
	fits->nsblk_nchan_npol = nsblk * nchan * npol;

	fits->data      = (char*) malloc ( sizeof(char)  * nsblk * nchan * npol );
	fits->weights   = (float*)malloc ( sizeof(float) * nchan );
	fits->scales    = (float*)malloc ( sizeof(float) * nchan * npol );
	fits->offsets   = (float*)malloc ( sizeof(float) * nchan * npol );
	fits->bata      = NULL;

	fits->bitdepth  = bitdepth;
	fits->fftint    = fftint;


	if (bw < 0) fits->bandwidth_mhz = -1.0 * bw;
	fits->bandwidth_mhz  = bw;
	fits->center_freqs   = (float*) malloc ( sizeof(float) * nchan );

	int ichan = 0;
	float  cbw      = bw / nchan;
	for (;ichan < nchan; ichan++) {
		fits->center_freqs [ ichan ] = fedge + ( ichan + 0.5 )*cbw;
	}

	fits->tsamp_tbin = fftint * nchan / bw / 1E6;
	fits->mjd        = mjd;
}

void gmrtfits_fold_prepare ( gmrtfits_t *fits, const char *filename, double mjd, unsigned npol, unsigned nchan, float fedge, float bw, unsigned nbin, double period ) {
	// sign(bw) determines LSB or USB
	
	fits->period   = period;
	fits->obsmode  = FOLD;
	
	fits->filename = strdup ( filename );
	
	fits->status  = 0;
	fits->nbin    = nbin;
	fits->npol    = npol;
	fits->nchan   = nchan;
	fits->nsblk   = 1;

	if ( fits->npol == 4 ) {
		sprintf ( fits->pol_type, "AABBCRCI" );
	} else if ( fits->npol == 1 ) {
		sprintf ( fits->pol_type, "AA+BB" );
	}


	fits->nchan_npol       = nchan * npol;
	fits->nbin_nchan_npol  = nbin * nchan * npol;

	fits->bata      = (int16_t*) malloc ( sizeof(int16_t)  * nbin * nchan * npol );
	/*fits->bata      = (float*) malloc ( sizeof(float)  * nbin * nchan * npol );*/
	fits->weights   = (float*)malloc ( sizeof(float) * nchan );
	fits->scales    = (float*)malloc ( sizeof(float) * nchan * npol );
	fits->offsets   = (float*)malloc ( sizeof(float) * nchan * npol );
	fits->data      = NULL;

	/* default set to unity; data is saved as int16 */
	fits->bitdepth  = 1;

	if (bw < 0) fits->bandwidth_mhz = -1.0 * bw;
	else fits->bandwidth_mhz  = bw;

	int ichan            = 0;
	fits->center_freqs   = (float*) malloc ( sizeof(float) * nchan );
	float  cbw           = bw / nchan;
	for (;ichan < nchan; ichan++) {
		fits->center_freqs [ ichan ] = fedge + ( ichan + 0.5 )*cbw;
		/* default initialize wts to unity */
		fits->weights [ ichan ]      = 1.0;
		fits->scales[ ichan ]        = 1.0;
		fits->offsets[ ichan ]       = 0.0;
	}

	fits->tsamp_tbin = period / nbin;
	fits->mjd        = mjd;
}

void gmrtfits_open ( gmrtfits_t *fits ) {
	fits->status = 0;

	int bitpix   = 8;
	int naxis    = 0;
	long *naxes  = NULL;
	int ival     = 0;
	char   cval  = 0;
	double dval  = 0.0;
	float  fval  = 0.0;
	long   lval  = 0;

	fits_create_file(&fits->fits, fits->filename, &fits->status);
	fits_create_img(fits->fits, bitpix, naxis, naxes, &fits->status);

	ival = 1;
	fits_update_key(fits->fits, TSTRING, "HDRVER","5.4gmrt","Header version, modified for uGMRT", &fits->status);
	fits_update_key(fits->fits, TSTRING, "HDRVER","5.4gmrt","Header version, modified for uGMRT", &fits->status);
	fits_update_key(fits->fits,TSTRING,"FITSTYPE","GMRTFITS","PSRFITS definition for uGMRT beamform data",&fits->status);

	fits_update_key(fits->fits,TSTRING,"TELESCOP","GMRT","upgraded Giant Metrewave Radio Telescope",&fits->status);

	// taken from observatories.dat from tempo2
	dval = 1656342.30; // ANT_X
	fits_update_key(fits->fits,TDOUBLE,"ANT_X", &dval, "[m] Antenna ITRF X-coordinate",&fits->status);
	dval = 5797947.77; // ANT_Y
	fits_update_key(fits->fits,TDOUBLE,"ANT_Y", &dval,"[m] Antenna ITRF Y-coordinate",&fits->status);
	dval = 2073243.16; // ANT_Z
	fits_update_key(fits->fits,TDOUBLE,"ANT_Z", &dval,"[m] Antenna ITRF Z-coordinate",&fits->status);

	// polarization basis
	// always circular
	fits_update_key(fits->fits, TSTRING, "FD_POLN", "CIRC", "LIN or CIRC", &fits->status);
	ival = 1;
	fits_update_key(fits->fits, TINT,"FD_HAND",&ival,"+/- 1. +1 is LIN:A=X,B=Y, CIRC:A=L,B=R (I) ", &fits->status);
	fval = 0.0;
	/*
	 * i did not know why FD_SANG was zero in case of circular basis for the longest. 
	 * Answer:
	 * For equal signal in R and L (A and B), there should not be any circular polarization.
	 * As noise diode only injects linear polarization, we are good.
	 */
	fits_update_key(fits->fits, TFLOAT,"FD_SANG",&fval,"[deg] FA of E vect for equal sig in A&B (E) ", &fits->status);

	// obsmode
	if ( fits->obsmode == SEARCH )
		fits_update_key (fits->fits, TSTRING, "OBS_MODE", "SEARCH", "SEARCH, PSR or CAL", &fits->status );
	else if ( fits->obsmode == FOLD )
		fits_update_key (fits->fits, TSTRING, "OBS_MODE", "PSR", "SEARCH, PSR or CAL", &fits->status );
	else
		fprintf ( stderr, " fits->obsmode not understood \n" );

	// frequency
	dval = fits->center_freqs [ fits->nchan / 2 ];
	fits_update_key (fits->fits, TDOUBLE, "OBSFREQ", &dval, "[MHz] Centre frequency for observation", &fits->status);
	dval = fits->bandwidth_mhz;
	fits_update_key (fits->fits, TDOUBLE, "OBSBW", &dval, "[MHz] Bandwidth of observation", &fits->status);
	ival = fits->nchan;
	fits_update_key (fits->fits, TINT, "NCHAN", &ival, "Number of channels", &fits->status);
	ival = fits->nbin;
	fits_update_key(fits->fits, TINT, "NBIN", &ival, "Nr of bins (PSR/CAL mode; else 1) ", &fits->status);

	// coordinates
	// XXX As of now, uGMRT does not provide pointing information with the raw file
	// so maybe ignore this.
	/*
		fits_update_key ( fits->fits, TSTRING, "COORD_MD", "J2000.0", "ICRS J2000", &fits->status );
		fval = 2000.0;
		fits_update_key ( fits->fits, TFLOAT, "EQUINOX", &fval, "Equinox", &fits->status );
		fits_update_key ( fits->fits, TSTRING, "RA", fits->rastr, "Right Ascension (hh:mm:ss.sss)", &fits->status );
		fits_update_key ( fits->fits, TSTRING, "DEC", fits->decstr, "Declination ([+-]dd:mm:ss.sss)", &fits->status );
		fits_update_key ( fits->fits, TDOUBLE, "RADEG", fits->ra, "[deg] Right Ascension", &fits->status );
		fits_update_key ( fits->fits, TDOUBLE, "DECDEG", fits->dec, "[deg] Declination", &fits->status );
	*/
	fits_update_key(fits->fits, TSTRING, "SRC_NAME", fits->source, "Source name", &fits->status);
	
	// time
	dval = fits->mjd;
	fits_update_key ( fits->fits, TDOUBLE, "STT_MJD", &dval, "MJD (UTC)" ,&fits->status );
	/*fits_update_key ( fits->fits, TSTRING, "STT_OBS", fits->datetime, "Human readable (UTC)" ,&fits->status );*/

	long imjd           = floor ( fits->mjd );
	double frac_seconds = 86400 * ( fits->mjd - (double)imjd );
	long seconds        = floor ( frac_seconds );
	double off_seconds  = frac_seconds - (double)seconds;

	fits_update_key ( fits->fits, TLONG, "STT_IMJD", &imjd, "Start MJD day (UTC)" ,&fits->status );
	fits_update_key ( fits->fits, TLONG, "STT_SMJD", &seconds, "Start Seconds (past 00h, UTC)" ,&fits->status );
	fits_update_key ( fits->fits, TDOUBLE, "STT_OFFS", &off_seconds, "[s] Start time offset" ,&fits->status );

	fits_write_date(fits->fits, &fits->status);

	if ( fits->status ) {
		printf (" at gmrtfits_create ...\n");
		fits_report_error(stderr, fits->status);
	}
}

void gmrtfits_data_table ( gmrtfits_t *fits ) {
	/*char *ttypes, *tforms, *tunits;*/

	/* adding period column in both modes */
	char *ttypes[] = { "ISUBINT", "TSUBINT", "OFFS_SUB", "NPOL", "PERIOD", "DAT_FREQ", "DAT_WTS", "DAT_OFFS", "DAT_SCL", "DATA"};
	char *tunits[] = {"", "s", "s", "", "s", "MHz", "", "", "", ""};
	char *tforms[10];
	int i;
	for (i = 0; i < 10; i++) tforms[i] = malloc ( sizeof(char) * 16 );
	long naxes [3];

	int ival;
	float fval;

	sprintf (tforms[0], "1J");
	sprintf (tforms[1], "1D");
	sprintf (tforms[2], "1D");
	sprintf (tforms[3], "1J");
	sprintf (tforms[4], "1D");
	sprintf (tforms[5], "%dE", fits->nchan );
	sprintf (tforms[6], "%dE", fits->nchan_npol );
	sprintf (tforms[7], "%dE", fits->nchan_npol );
	sprintf (tforms[8], "%dE", fits->nchan_npol );
	if ( fits->obsmode == SEARCH ) sprintf (tforms[9], "%dB", fits->nsblk_nchan_npol );
	else if ( fits->obsmode == FOLD ) sprintf (tforms[9], "%dI", fits->nbin_nchan_npol );
	/*else if ( fits->obsmode == FOLD ) sprintf (tforms[9], "%dE", fits->nbin_nchan_npol );*/

	fits_create_tbl( fits->fits, BINARY_TBL, 0, 10, ttypes, tforms, tunits, "SUBINT", &fits->status);
	// free them
	for (i = 0; i < 10; i++) free ( tforms[i] );

	// cfitsio fits needs me to save in reverse order
	if ( fits->obsmode == SEARCH ) {
		/* documentation = (nchan, npol, nsblk) */
		/* data          = (nsblk, npol, nchan) */
		naxes[0] = fits->nsblk;
		naxes[1] = fits->npol;
		naxes[2] = fits->nchan;
	} else if ( fits->obsmode == FOLD ) {
		/* documentation = (nbin, nchan, npol) */
		/* data          = (nbin, nchan, npol) */
		naxes[0] = fits->nbin;
		naxes[1] = fits->nchan;
		naxes[2] = fits->npol;
	}

	fits_write_tdim( fits->fits, 10, 3, naxes, &fits->status );

	fits->nrow     = 1;

	// header
	fits_update_key(fits->fits,TSTRING,"INT_TYPE","TIME","search mode in time",&fits->status);
	fits_update_key(fits->fits,TSTRING,"INT_UNIT","SEC","Unit of time axis",&fits->status);
	fits_update_key(fits->fits,TSTRING,"SCALE","SAMPLES","Digital samples.",&fits->status);

	// npol, nchan
	ival  = fits->npol;
	fits_update_key(fits->fits,TINT,"NPOL",&ival,"Number of polarization.",&fits->status);
	ival  = fits->nchan;
	fits_update_key(fits->fits,TINT,"NCHAN",&ival,"Number of frequency channels.",&fits->status);
	ival  = fits->nsblk;
	fits_update_key(fits->fits,TINT,"NSBLK",&ival,"Number of samples per block.",&fits->status);
	ival  = fits->bitdepth;
	fits_update_key(fits->fits,TINT,"NBITS",&ival,"Number of bits per datum.",&fits->status);
	ival  = fits->nbin;
	fits_update_key(fits->fits, TINT, "NBIN", &ival, "Nr of bins (PSR/CAL mode; else 1) ", &fits->status);
	fval  = fits->tsamp_tbin;
	fits_update_key(fits->fits, TFLOAT, "TBIN", &fval, "[s] Time per bin or sample ", &fits->status);
	fval  = fits->dm;
	fits_update_key(fits->fits, TFLOAT, "DM", &fval, "[pc/cc] dispersion measure", &fits->status);

	fits_update_key(fits->fits,TSTRING,"POL_TYPE",fits->pol_type,"Polarization type.",&fits->status);

	//
	fits_flush_file (fits->fits, &fits->status);
}

void gmrtfits_close ( gmrtfits_t *fits ) {
	free ( fits->center_freqs );
	free ( fits->offsets );
	free ( fits->scales );
	free ( fits->weights );

	if (fits->obsmode == SEARCH) free ( fits->data );
	else if (fits->obsmode == FOLD) free ( fits->bata );

	free ( fits->filename );
	

	fits_close_file ( fits->fits, &fits->status );
	if ( fits->status ) {
		fits_report_error(stderr, fits->status);
	}
}

void gmrtfits_search_add ( gmrtfits_t *fits, real_t *data, unsigned int ngulp ) {
	// the data layout is confusing since CFITSIO behaves like FORTRAN with column major ordering
	// ngulp = k * nsblk for some k
	// will write multiple times in the same function call
	// need to write (start:(start+nsblk), npol, nchan)
	//
	
	int colnum  = 0;
	int    ival = 0;
	double dval = 0.0;

	// normalize data
	// scales and offsets are (nchan, npol)
	// DATA = DAT_WTS * ( BATA*SCALES + OFFSETS )
	
	// iterators
	unsigned isamp, osamp, ichan, ipol;
	/*unsigned stop  = start + fits->nsblk;*/
	unsigned int start, stop;

	// compute min and max (nchan, npol)
	// min is offset and max-min is scales
	real_t _offset, _scale;

#ifdef TIMING
	clock_t tstart, tstop;
	unsigned rowwrites  = 0;
	double time_process = 0., time_fitsio = 0., time_flush = 0.;
#endif

	for (start = 0; start < ngulp; start+=fits->nsblk) {
		stop    = start + fits->nsblk;
#ifdef TIMING
		tstart  = clock ();
#endif
		// data is (nsblk, npol, nchan)
		for (ipol = 0; ipol < fits->npol; ipol++) {
			for (ichan = 0; ichan < fits->nchan; ichan++) {

				real_t _min = FLT_MAX, _max = -FLT_MAX;

				// find min and max per (nchan,npol)
				for (isamp  = start; isamp < stop; isamp++) {
					real_t _d = data [ ichan + fits->nchan*ipol + fits->nchan_npol*isamp ];

					if ( _d < _min ) _min = _d;
					if ( _d > _max ) _max = _d;
				} // nsamp
					
				_offset    = _min;
				_scale     = (_max - _min) / 255.0;

				fits->offsets [ ichan + fits->nchan*ipol ] = _offset;
				fits->scales  [ ichan + fits->nchan*ipol ] = _scale;

				for (isamp = start,osamp = 0; isamp < stop; isamp++, osamp++) {
					real_t _d = data [ ichan + fits->nchan*ipol + fits->nchan_npol*isamp ];
					fits->data [ ichan + fits->nchan*ipol + fits->nchan_npol*osamp ] = (char) ( ( _d - _offset ) / _scale );
					/*fits->data [ ichan + fits->nchan*ipol + fits->nchan_npol*osamp ] = (char) ( ipol + (64*ichan/2048));*/
				} // nsamp
			} // freq
		} // pol
#ifdef TIMING
		tstop        = clock ();
		time_process += (tstop - tstart) / CLOCKS_PER_SEC;
		tstart       = clock ();
#endif

		// ISUBINT
		colnum  = 1;
		fits_write_col( fits->fits, TINT, colnum, fits->nrow, 1, 1, &fits->nrow, &fits->status );

		// TSUBINT
		colnum  = 2;
		dval    = fits->nsblk * fits->tsamp_tbin;
		fits_write_col( fits->fits, TDOUBLE, colnum, fits->nrow, 1, 1, &dval, &fits->status );
		
		// OFFS_SUB
		colnum  = 3;
		dval    =  (fits->nrow + 0.5) * fits->nsblk * fits->tsamp_tbin;
		fits_write_col( fits->fits, TDOUBLE, colnum, fits->nrow, 1, 1, &dval, &fits->status );
		// middle of the subint
		
		// NPOL
		colnum  = 4;
		ival    = fits->npol;
		fits_write_col( fits->fits, TINT, colnum, fits->nrow, 1, 1, &ival, &fits->status );

		// DAT_FREQ
		colnum  = 5;
		fits_write_col( fits->fits, TFLOAT, colnum, fits->nrow, 1, fits->nchan, fits->center_freqs, &fits->status );

		// DAT_WTS
		colnum  = 6;
		fits_write_col( fits->fits, TFLOAT, colnum, fits->nrow, 1, fits->nchan_npol, fits->weights, &fits->status );

		// DAT_OFFS
		colnum  = 7;
		fits_write_col( fits->fits, TFLOAT, colnum, fits->nrow, 1, fits->nchan_npol, fits->offsets, &fits->status );

		// DAT_SCL
		colnum  = 8;
		fits_write_col( fits->fits, TFLOAT, colnum, fits->nrow, 1, fits->nchan_npol, fits->scales, &fits->status );

		// DATA
		colnum  = 9;
		fits_write_col( fits->fits, TBYTE, colnum, fits->nrow, 1, fits->nsblk_nchan_npol, fits->data, &fits->status );


		if ( fits->status ) {
			printf ("\n at gmrtfits_subint_real_single row=%ld ...\n", fits->nrow);
			fits_report_error(stderr, fits->status);
			raise(SIGINT);
			/*exit (1);*/
		}

#ifdef TIMING
		tstop        = clock ();
		time_fitsio  += (tstop - tstart) / CLOCKS_PER_SEC;
		tstart       = clock ();
#endif
		if ( fits->nrow % FLUSH_TIMES == 0 ) {
			fits_flush_buffer( fits->fits, 0, &fits->status );
		}
#ifdef TIMING
		tstop        = clock ();
		time_flush   += (tstop - tstart) / CLOCKS_PER_SEC;
		rowwrites++;
#endif

		fits->nrow++;
	} // nsblk write loop
		
#ifdef TIMING
	printf ("[TIMING] written rows=%d process=%.5f fitsio=%.5f flush=%.5f\n", rowwrites, time_process, time_fitsio, time_flush);
#endif
} // search_add

void gmrtfits_fold_add ( gmrtfits_t *fits, deal_t *data, double tsubint, double period ) {
	/* data is (nbin, nchan, npol)  */
	
	int colnum  = 0;
	int    ival = 0;
	double dval = 0.0;

	// normalize data
	// scales and offsets are (nchan, npol)
	// DATA = DAT_WTS * ( BATA*SCALES + OFFSETS )
	
	// iterators
	unsigned ibin, osamp, ichan, ipol;
	unsigned int start, stop;

	// compute min and max (nchan, npol)
	// min is offset and max-min is scales
	real_t _offset, _scale;

#ifdef TIMING
	clock_t tstart, tstop;
	unsigned rowwrites  = 0;
	double time_process = 0., time_fitsio = 0., time_flush = 0.;
#endif

#ifdef TIMING
		tstart  = clock ();
#endif
	// data is (nbin, nchan, npol)
	for (ichan = 0; ichan < fits->nchan; ichan++) {
		for (ipol = 0; ipol < fits->npol; ipol++) {
			deal_t _min = FLT_MAX, _max = -FLT_MAX;

			// find min and max per (nchan,npol)
			for (ibin = 0; ibin < fits->nbin; ibin++) {
				deal_t _d = data [ ipol + fits->npol*ichan + fits->nchan_npol*ibin ];

				if ( _d < _min ) _min = _d;
				if ( _d > _max ) _max = _d;
			} // bin

			/**
			 * min -> INT16_MIN
			 * max -> INT16_MAX
			 *
			 * (b - bmin) / (d - min) = ( INT16_MAX - INT16_MIN ) / ( max - min )
			 * (b - bmin) = ( INT16_MAX - INT16_MIN ) / ( max - min ) * ( d - min )
			 * b  = (INT16_MAX-INT16_MIN)/(max-min) * d + ( bmin - min*(INT16_MAX-INT16_MIN)/(max-min) )
			 * ss = (INT16_MAX-INT16_MIN)/(max-min) 
			 * b  = ss * d + ( INT16_MIN - min*ss ) 
			 *
			 * bmin -> dmin
			 * bmax -> dmax
			 *
			 * slope = 
			 *
			 * data = scl*bata + offs
			 *
			 * scl  = (dmax - dmin) / (INT16_MAX - INT16_MIN)
			 *
			 * (d - dmin) / ( b - bmin ) = scl
			 * d - dmin = scl*b - scl*bmin
			 * d = scl*b + dmin - scl*bmin
			 *
			 *
			 **/
					
			_scale     = (_max - _min) / (INT16_MAX - INT16_MIN);
			_offset    = _min - _scale*INT16_MIN;

			fits->offsets [ ichan + fits->nchan*ipol ] = _offset;
			fits->scales  [ ichan + fits->nchan*ipol ] = _scale;

			for (ibin = 0; ibin < fits->nbin; ibin++) {
				size_t idx = ipol + fits->npol*ichan + fits->nchan_npol*ibin;
				size_t jdx = ibin + fits->nbin*ichan + fits->nbin*fits->nchan*ipol;
				deal_t _d  = data [ idx ];
				fits->bata [ jdx ] = (int16_t) ( ( _d - _offset ) / _scale );
			} // bin
		} // pol
	} // chan
#ifdef TIMING
		tstop        = clock ();
		time_process += (tstop - tstart) / CLOCKS_PER_SEC;
		tstart       = clock ();
#endif

		// ISUBINT
		colnum  = 1;
		fits_write_col( fits->fits, TINT, colnum, fits->nrow, 1, 1, &fits->nrow, &fits->status );

		// TSUBINT
		colnum  = 2;
		dval    = tsubint;
		fits_write_col( fits->fits, TDOUBLE, colnum, fits->nrow, 1, 1, &dval, &fits->status );
		
		// OFFS_SUB
		colnum  = 3;
		dval    =  (fits->nrow + 0.5) * tsubint;
		fits_write_col( fits->fits, TDOUBLE, colnum, fits->nrow, 1, 1, &dval, &fits->status );
		// middle of the subint
		
		// NPOL
		colnum  = 4;
		ival    = fits->npol;
		fits_write_col( fits->fits, TINT, colnum, fits->nrow, 1, 1, &ival, &fits->status );

		// PERIOD
		colnum  = 5;
		dval    =  period;
		fits_write_col( fits->fits, TDOUBLE, colnum, fits->nrow, 1, 1, &dval, &fits->status );

		// DAT_FREQ
		colnum  = 6;
		fits_write_col( fits->fits, TFLOAT, colnum, fits->nrow, 1, fits->nchan, fits->center_freqs, &fits->status );

		// DAT_WTS
		colnum  = 7;
		fits_write_col( fits->fits, TFLOAT, colnum, fits->nrow, 1, fits->nchan_npol, fits->weights, &fits->status );

		// DAT_OFFS
		colnum  = 8;
		fits_write_col( fits->fits, TFLOAT, colnum, fits->nrow, 1, fits->nchan_npol, fits->offsets, &fits->status );

		// DAT_SCL
		colnum  = 9;
		fits_write_col( fits->fits, TFLOAT, colnum, fits->nrow, 1, fits->nchan_npol, fits->scales, &fits->status );

		// DATA
		colnum  = 10;
		fits_write_col( fits->fits, TSHORT, colnum, fits->nrow, 1, fits->nbin_nchan_npol, fits->bata, &fits->status );
		/*fits_write_col( fits->fits, TFLOAT, colnum, fits->nrow, 1, fits->nbin_nchan_npol, fits->bata, &fits->status );*/


		if ( fits->status ) {
			printf ("\n at gmrtfits_subint_real_single row=%ld ...\n", fits->nrow);
			fits_report_error(stderr, fits->status);
			raise(SIGINT);
			/*exit (1);*/
		}

#ifdef TIMING
		tstop        = clock ();
		time_fitsio  += (tstop - tstart) / CLOCKS_PER_SEC;
		tstart       = clock ();
#endif
		if ( fits->nrow % FLUSH_TIMES == 0 ) {
			fits_flush_buffer( fits->fits, 0, &fits->status );
		}
#ifdef TIMING
		tstop        = clock ();
		time_flush   += (tstop - tstart) / CLOCKS_PER_SEC;
		rowwrites++;
#endif

		fits->nrow++;

#ifdef TIMING
	printf ("[TIMING] written rows=%d process=%.5f fitsio=%.5f flush=%.5f\n", rowwrites, time_process, time_fitsio, time_flush);
#endif
} // fold_add
