#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include "gmrtfits.h"

#define NPOL 4
#define NCHAN 2048
#define NSBLK 4096
#define NGULP 4096
// time integration
#define INTEGRATION 1
#define OGULP ( NGULP / INTEGRATION )
#define NSUBPERGULP ( OGULP / NSBLK )
// channelizer sizes
// These are per-pol
#define VREAD ( 2 * NCHAN * NGULP )
#define FFTCOMPLEXSIZE (1 + NCHAN*NGULP)

void print_help() {
	fprintf(stderr, "\n");
	fprintf(stderr, "  dualpol8bit - Transforms 8bit GMRT GWB baseband data into search mode PSRFITS\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "    -h Print this help\n");
	fprintf(stderr, "    -r <rraw> Right Circular Polarized baseband file\n");
	fprintf(stderr, "    -l <lraw> Left Circular Polarized baseband file\n");
	fprintf(stderr, "    -f <freq> Frequency edge of recording in MHz\n");
	fprintf(stderr, "    -b <fbw> Signed bandwidth in MHz where sign determines USB(positive) or LSB(negative)\n");
	fprintf(stderr, "    -t <mjd> MJD in maximum available precision\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " Notes:\n");
	fprintf(stderr, "    Time averaging is fixed at 16. Change in source and recompile if needed to change\n");
	fprintf(stderr, "    Bitdepth is 8. Changing it requires additional logic to pack/unpack data.\n");
}

int main() {
	print_help ();
	const char *pol1_infile_path = "/tmp/baseband/C02pol1200MHz8bitsCRAB550MHz_10s.raw";
	const char *pol2_infile_path = "/tmp/baseband/C02pol1200MHz8bitsCRAB550MHz_10s.raw";
	/*const char *toufile_path = "/tmp/baseband/testfb_float.raw";*/
	const char *oufile_path = "/tmp/baseband/testfs.fits";
	double mjd            = 60900.23741898148;

	/*FILE *infile, *outfile;*/
	FILE *pol1_infile, *pol2_infile;
	FILE *tou1, *tou2, *tou3;
	FILE *toutfile;

	char *pol1_volt_read  = (char*) malloc ( VREAD * sizeof (char) );
	char *pol2_volt_read  = (char*) malloc ( VREAD * sizeof (char) );
	double *volt_double   = (double*) malloc ( 2 * VREAD * sizeof(double) );
	float *outfb          = (float*)malloc ( NPOL * OGULP * NCHAN * sizeof(float) );
	//-------------------------//
	gmrtfits_t gf;
	gmrtfits_prepare ( &gf, oufile_path, mjd, NPOL, NCHAN, 550.0f, 200.0f, NSBLK, INTEGRATION, 8 );
	gmrtfits_create ( &gf );
	//-------------------------//


	fftw_plan pforward, pbackward1, pbackward2;
	fftw_complex *fdata;

	unsigned long read_bytes;
	unsigned long infilesize1, infilesize2;

	unsigned iread;
	unsigned nreads;
	unsigned ichan, isamp, iavg;

	fdata = (fftw_complex*) fftw_alloc_complex  ( 2 * FFTCOMPLEXSIZE );

	// files open
	pol1_infile  = fopen ( pol1_infile_path, "r" );
	pol2_infile  = fopen ( pol2_infile_path, "r" );

	/*tou1         = fopen ( "/tmp/baseband/testfb_data.raw" , "w" );*/
	/*tou2         = fopen ( "/tmp/baseband/testfb_scales.raw" , "w" );*/
	/*tou3         = fopen ( "/tmp/baseband/testfb_offsets.raw" , "w" );*/
	/*toutfile = fopen ( toufile_path, "w+" );*/

	// get file size
	fseek ( pol1_infile, 0L, SEEK_END );
	infilesize1 = ftell ( pol1_infile );
	fseek ( pol1_infile, 0L, SEEK_SET );
	fseek ( pol2_infile, 0L, SEEK_END );
	infilesize2 = ftell ( pol2_infile );
	fseek ( pol2_infile, 0L, SEEK_SET );
	//
	if ( infilesize1 != infilesize2 ) {
		fprintf(stderr, " filesizes pol1=%lu pol2=%lu\n", infilesize1, infilesize2);
		fprintf(stderr, " Unexpected case.");
	}

	/*nreads  = infilesize1 / VREAD;*/
	nreads  = 42;

	int ng = VREAD;
	// forward FFT
	pforward = fftw_plan_many_dft_r2c ( 1, &ng, 2, volt_double, NULL, 1, VREAD, fdata, NULL, 1, FFTCOMPLEXSIZE, FFTW_ESTIMATE  );
	// per pol :: real VREAD -> complex NCHAN*NGULP + 1
	// ignore the DC term
	ng = NGULP;
	/**
		CHECK 
	 * FFT R2C has the DC term as the first element.
	 * We ignore this DC term while coherent filterbanking.
	 * This messes up the array layout that forces us 
	 * to make two fftwplans for backward FFT
	 **/
	pbackward1 = fftw_plan_many_dft ( 1, &ng, NCHAN, fdata+1, NULL, 1, NGULP, fdata+1, NULL, 1, NGULP, FFTW_BACKWARD, FFTW_ESTIMATE );
	pbackward2 = fftw_plan_many_dft ( 1, &ng, NCHAN, fdata+FFTCOMPLEXSIZE+1, NULL, 1, NGULP, fdata+FFTCOMPLEXSIZE+1, NULL, 1, NGULP, FFTW_BACKWARD, FFTW_ESTIMATE );

	printf (" prepared plans .. starting ...\n");
	gmrtfits_subint_open ( &gf );

	for (iread = 0; iread < nreads; iread++) {
		printf( " iread=%d .. ", iread );
		// read
		read_bytes = fread ( pol1_volt_read, sizeof(char), VREAD, pol1_infile );
		read_bytes = fread ( pol2_volt_read, sizeof(char), VREAD, pol2_infile );
		// if need to unpack
		// unpack here
		for ( int i = 0; i < VREAD; i++ ) {
			volt_double[i]         = (double) pol1_volt_read[i];
			volt_double[VREAD + i] = (double) pol2_volt_read[i];
		}

		printf( " read .. " );

		// forward FFT
		fftw_execute ( pforward );

		printf( " forward FFT .. " );

		// backward FFT
		fftw_execute ( pbackward1 );
		fftw_execute ( pbackward2 );

		printf( " backward FFT .. " );

		// detection
		/**
		 *
		 * backward FFT
		 * skipping the DC element in both POLS
		 * POL1 || F1T1 F1T2 F1T3 .... F1TN
		 *      || F2T1 F2T2 F2T3 .... F2TN
		 *      || :::: :::: :::: .... ::::
		 *      || FNT1 FNT2 FNT3 .... FNTN
		 * --------------------------------
		 * POL2 || F1T1 F1T2 F1T3 .... F1TN
		 *      || F2T1 F2T2 F2T3 .... F2TN
		 *      || :::: :::: :::: .... ::::
		 *      || FNT1 FNT2 FNT3 .... FNTN
		 *
		 * Coherence products come from 
		 * same time-frequency complex pixel of the two pols
		 *
		 * (1) i understand now why PSRFITS searchmode has layout
		 * (nchan, npol, nsblk)
		 * (2) i do not. PSRFITS uses CFITSIO which is column ordered
		 * It should be sent to CFITSIO as (nsamp, npol, nchan)
		 * which will be stored as (nchan, npol, nsamp) in CFITSIO.
		 * (3) no layout seems to work, now going with 
		 **/

		// do direct read
		// zeroout the outfb
		memset ( outfb, 0, NPOL * OGULP * NCHAN * sizeof(float) );
		unsigned long int ii = 1, jj = 0, kk = 0;
		double p1real, p1imag, p2real, p2imag;
		float aa, bb, cr, ci;
		// (1) data is (nchan, npol, nsblk) 
		// (2) data is (nsblk, npol, nchan) 
		for (ichan    = 0; ichan < NCHAN; ichan++) {
			for (isamp  = 0; isamp < NGULP; isamp++) {
				// input index  - FT
				ii        = isamp + NGULP*ichan;
				// output index - FPT (1)
				// isamp + ipol*NGULP + NGULP*NPOL*ichan 
				/*jj        = isamp + NPOL*NGULP*ichan;*/
				// output index - TPF (2)
				// ichan + NCHAN*ipol + NCHAN*NPOL*isamp
				jj        = ichan +  NCHAN*NPOL*isamp;
				// time-averaged output index
				kk        = ichan + NCHAN*NPOL*isamp/INTEGRATION;
				/*kk        = isamp/INTEGRATION + NPOL*OGULP*ichan;*/
				// extract same time frequency pixels
				// do offbyone to ignore the DC term
				p1real    = fdata [ 1 + ii ][0];
				p1imag    = fdata [ 1 + ii ][1];
				p2real    = fdata [ 1 + FFTCOMPLEXSIZE + ii ][0];
				p2imag    = fdata [ 1 + FFTCOMPLEXSIZE + ii ][1];
				// compute coherences
				// (p1r,p1i) * (p2r,-p2i)
				aa        = p1real*p1real + p1imag*p1imag;
				bb        = p2real*p2real + p2imag*p2imag;
				cr        = p1real*p2real + p1imag*p2imag;
				ci        = p1imag*p2real - p1real*p2imag;
				// write into outfb
				outfb [ kk           ]  += aa;
				outfb [ kk +   NCHAN ]  += bb;
				outfb [ kk + 2*NCHAN ]  += cr;
				outfb [ kk + 3*NCHAN ]  += ci;
			} // samp
		} // channel

		printf( " detection .. " );

		// write to file
		/*fwrite ( outfb, sizeof(float), OGULP * NCHAN * NPOL, toutfile );*/
		/*goto exit;*/

		/*
		 * for the case when NGULP > NSBLK
		 * possibly bug here
		 * when OGULP does not nicely divide 
		 * outfb gets zeroed out and that data is lost
		 */
		for (unsigned int i = 0; i < OGULP; i+=NSBLK) {
			printf ("    writing subint i=%d\n", i);
			gmrtfits_subint_real ( &gf, outfb, i, OGULP );
		}

		printf( " written \n" );
		//
		/*fwrite ( gf.data, sizeof(char), OGULP * NCHAN * NPOL, tou1 );*/
		/*fwrite ( gf.scales, sizeof(float), NCHAN * NPOL, tou2 );*/
		/*fwrite ( gf.offsets, sizeof(float), NCHAN * NPOL, tou3 );*/
		//
		/*goto exit;*/
	}

exit:
	gmrtfits_close ( &gf );

	// destroy plans
	fftw_destroy_plan ( pforward );
	fftw_destroy_plan ( pbackward1 );
	fftw_destroy_plan ( pbackward2 );

	fftw_free ( fdata );
	free ( pol1_volt_read );
	free ( pol2_volt_read );
	free (volt_double);
	free (outfb);

	fclose ( pol1_infile );
	fclose ( pol2_infile );

	/*fclose ( tou1 );*/
	/*fclose ( tou2 );*/
	/*fclose ( tou3 );*/
	/*fclose ( toutfile );*/

	return 0;
}

