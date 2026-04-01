#include "gmrtfits.h"

#ifdef TIMING
#include <time.h>
#endif

void gmrtfits_prepare ( gmrtfits_t *fits, const char *filename, double mjd, unsigned npol, unsigned nchan, float fedge, float bw, unsigned nsblk, unsigned fftint, unsigned bitdepth ) {
	// sign(bw) determines LSB or USB
	
	fits->filename = strdup ( filename );
	
	fits->status  = 0;
	fits->npol    = npol;
	fits->nchan   = nchan;
	fits->nsblk   = nsblk;

	sprintf ( fits->pol_type, "AABBCRCI" );

	fits->nchan_npol       = nchan * npol;
	fits->nsblk_nchan_npol = nsblk * nchan * npol;

	fits->data      = (char*) malloc ( sizeof(char)  * nsblk * nchan * npol );
	fits->weights   = (float*)malloc ( sizeof(float) * nchan );
	fits->scales    = (float*)malloc ( sizeof(float) * nchan * npol );
	fits->offsets   = (float*)malloc ( sizeof(float) * nchan * npol );

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

	fits->tsamp     = fftint * nchan / bw / 1E6;
	fits->mjd       = mjd;
}

void gmrtfits_close ( gmrtfits_t *fits ) {
	free ( fits->center_freqs );
	free ( fits->offsets );
	free ( fits->scales );
	free ( fits->weights );
	free ( fits->data );

	free ( fits->filename );
	

	fits_close_file ( fits->fits, &fits->status );
	if ( fits->status ) {
		fits_report_error(stderr, fits->status);
	}
}

void gmrtfits_create ( gmrtfits_t *fits ) {
	fits->status = 0;

	int bitpix   = 8;
	int naxis    = 0;
	long *naxes  = NULL;
	int ival     = 0;
	double dval  = 0.0;
	float  fval  = 0.0;
	long   lval  = 0;

	fits_create_file(&fits->fits, fits->filename, &fits->status);
	fits_create_img(fits->fits, bitpix, naxis, naxes, &fits->status);

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
	fits_update_key (fits->fits, TSTRING, "OBS_MODE", "SEARCH", "SEARCH by default", &fits->status );

	// frequency
	dval = fits->center_freqs [ fits->nchan / 2 ];
	fits_update_key (fits->fits, TDOUBLE, "OBSFREQ", &dval, "[MHz] Centre frequency for observation", &fits->status);
	dval = fits->bandwidth_mhz;
	fits_update_key (fits->fits, TDOUBLE, "OBSBW", &dval, "[MHz] Bandwidth of observation", &fits->status);
	ival = fits->nchan;
	fits_update_key (fits->fits, TINT, "NCHAN", &ival, "Number of channels", &fits->status);
	ival = 1;
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

void gmrtfits_subint_open ( gmrtfits_t *fits ) {
	/*char *ttypes, *tforms, *tunits;*/

	char *ttypes[] = { "ISUBINT", "TSUBINT", "OFFS_SUB", "NPOL", "DAT_FREQ", "DAT_WTS", "DAT_OFFS", "DAT_SCL", "DATA"};
	char *tunits[] = {"", "s", "s", "", "MHz", "", "", "", ""};
	/*char tforms[9][16];*/
	char *tforms[9];
	int i;
	for (i = 0; i < 9; i++) tforms[i] = malloc ( sizeof(char) * 16 );
	long naxes [3];

	int ival;
	float fval;

	sprintf (tforms[0], "1J");
	sprintf (tforms[1], "1D");
	sprintf (tforms[2], "1D");
	sprintf (tforms[3], "1J");
	sprintf (tforms[4], "%dD", fits->nchan );
	sprintf (tforms[5], "%dD", fits->nchan_npol );
	sprintf (tforms[6], "%dD", fits->nchan_npol );
	sprintf (tforms[7], "%dD", fits->nchan_npol );
	sprintf (tforms[8], "%dB", fits->nsblk_nchan_npol );

	fits_create_tbl( fits->fits, BINARY_TBL, 0, 9, ttypes, tforms, tunits, "SUBINT", &fits->status);
	// free them
	for (i = 0; i < 9; i++) free ( tforms[i] );

	// astropy.io.fits needs me to save in reverse order
	naxes[0] = fits->nsblk;
	naxes[1] = fits->npol;
	naxes[2] = fits->nchan;

	fits_write_tdim( fits->fits, 9, 3, naxes, &fits->status );

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
	ival  = 1;
	fits_update_key(fits->fits, TINT, "NBIN", &ival, "Nr of bins (PSR/CAL mode; else 1) ", &fits->status);
	fval  = fits->tsamp;
	fits_update_key(fits->fits, TFLOAT, "TBIN", &fval, "[s] Time per bin or sample ", &fits->status);

	fits_update_key(fits->fits,TSTRING,"POL_TYPE",fits->pol_type,"Polarization type.",&fits->status);

	//
	fits_flush_file (fits->fits, &fits->status);

}

void gmrtfits_subint_real ( gmrtfits_t *fits, real_t *data, unsigned int start, unsigned int ngulp ) {
	// the data layout is confusing since CFITSIO behaves like FORTRAN with column major ordering
	// need to write (start:(start+nsblk), npol, nchan)
	
	int colnum  = 0;
	int    ival = 0;
	double dval = 0.0;

	// normalize data
	// scales and offsets are (nchan, npol)
	// DATA = DAT_WTS * ( BATA*SCALES + OFFSETS )
	
	// iterators
	unsigned isamp, osamp, ichan, ipol;
	unsigned stop  = start + fits->nsblk;

	// compute min and max (nchan, npol)
	// min is offset and max-min is scales
	real_t _offset, _scale;

#ifdef TIMING
	clock_t tstart, tstop;
	double time_process, time_fitsio, time_flush;
	tstart  = clock ();
#endif

	// data is (nsblk, npol, nchan)
	// i am sure of this
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
	time_process = (tstop - tstart) / CLOCKS_PER_SEC;
	tstart       = clock ();
#endif
	
	// ISUBINT
	colnum  = 1;
	fits_write_col( fits->fits, TINT, colnum, fits->nrow, 1, 1, &fits->nrow, &fits->status );

	// TSUBINT
	colnum  = 2;
	dval    = fits->nsblk * fits->tsamp;
	fits_write_col( fits->fits, TDOUBLE, colnum, fits->nrow, 1, 1, &dval, &fits->status );
	
	// OFFS_SUB
	colnum  = 3;
	dval    =  (fits->nrow + 0.5) * fits->nsblk * fits->tsamp;
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
	time_fitsio  = (tstop - tstart) / CLOCKS_PER_SEC;
	tstart       = clock ();
#endif

	fits_flush_buffer( fits->fits, 0, &fits->status );
	fits->nrow++;

#ifdef TIMING
	tstop        = clock ();
	time_flush   = (tstop - tstart) / CLOCKS_PER_SEC;

	printf ("[TIMING] process=%.5f fitsio=%.5f flush=%.5f\n", time_process, time_fitsio, time_flush);
#endif
}

void gmrtfits_subint_add ( gmrtfits_t *fits, real_t *data, unsigned int ngulp ) {
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
		dval    = fits->nsblk * fits->tsamp;
		fits_write_col( fits->fits, TDOUBLE, colnum, fits->nrow, 1, 1, &dval, &fits->status );
		
		// OFFS_SUB
		colnum  = 3;
		dval    =  (fits->nrow + 0.5) * fits->nsblk * fits->tsamp;
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

		fits_flush_buffer( fits->fits, 0, &fits->status );
		fits->nrow++;

#ifdef TIMING
		tstop        = clock ();
		time_flush   += (tstop - tstart) / CLOCKS_PER_SEC;
		rowwrites++;
#endif
	} // nsblk write loop
		
#ifdef TIMING
	printf ("[TIMING] written rows=%d process=%.5f fitsio=%.5f flush=%.5f\n", rowwrites, time_process, time_fitsio, time_flush);
#endif

} // function

