#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include "gmrtfits.h"

#define NPOL 4
#define NCHAN 2048
#define NSBLK 2048
#define NGULP 16384
// time integration
#define INTEGRATION 8
#define OGULP ( NGULP / INTEGRATION )
#define NSUBPERGULP ( OGULP / NSBLK )
// channelizer sizes
// These are per-pol
#define VREAD ( 2 * NCHAN * NGULP )
#define FFTCOMPLEXSIZE (1 + NCHAN*NGULP)

void print_help() {
	fprintf(stderr, "\n");
	fprintf(stderr, "dualpol8bit - Channelize and detect 8bit GMRT GWB baseband data\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "\t -h Print this help\n");
	fprintf(stderr, "\t -r <rraw> Right Circular Polarized baseband file\n");
	fprintf(stderr, "\t -l <lraw> Left Circular Polarized baseband file\n");
	fprintf(stderr, "\t -f <freq> Frequency edge of recording in MHz\n");
	fprintf(stderr, "\t -b <fbw> Signed bandwidth in MHz where sign determines USB(positive) or LSB(negative)\n");
	fprintf(stderr, "\t -t <mjd> MJD in maximum available precision\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Notes:\n");
	fprintf(stderr, "\t Time averaging is fixed at 16. Change in source and recompile if needed to change\n");
}

int main() {
	print_help ();
	const char *pol1_infile_path = "/tmp/baseband/C02pol1200MHz8bitsCRAB550MHz_10s.raw";
	const char *pol2_infile_path = "/tmp/baseband/C02pol1200MHz8bitsCRAB550MHz_10s.raw";
	/*const char *oufile_path = "/tmp/baseband/testfbfloat32.raw";*/
	const char *oufile_path = "/tmp/baseband/testfs.fits";
	double mjd            = 60900.23741898148;

	/*FILE *infile, *outfile;*/
	FILE *pol1_infile, *pol2_infile;

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
	/*outfile = fopen ( oufile_path, "w+" );*/

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
	nreads  = 16;

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
		 * i understand now why PSRFITS searchmode has layout
		 * (nchan, npol, nsblk)
		 **/

		// do direct read
		// zeroout the outfb
		memset ( outfb, 0, NPOL * OGULP * NCHAN * sizeof(float) );
		unsigned long int ii = 1, jj = 0, kk = 0;
		double p1real, p1imag, p2real, p2imag;
		float aa, bb, cr, ci;
		for (ichan = 0; ichan < NCHAN; ichan++) {
			for (isamp = 0; isamp < NGULP; isamp++) {
				// input index
				ii        = isamp + NGULP*ichan;
				// output index
				jj        = isamp + NGULP*ichan;
				// time-averaged output index
				kk        = isamp/INTEGRATION + OGULP*NPOL*ichan;
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
				outfb [ kk +   OGULP ]  += bb;
				outfb [ kk + 2*OGULP ]  += cr;
				outfb [ kk + 3*OGULP ]  += ci;
			} // samp
		} // channel

		printf( " detection .. " );

		// write to file
		/*fwrite ( outfb, sizeof(float), OGULP * NCHAN, outfile );*/
		/*
		 * (nchan, npol, ngulp) where ngulp > nsblk
		 *
		 * igulp + ngulp*ipol + ngulp*npol*ichan
		 */
		for (unsigned int i = 0; i < OGULP; i+=NSBLK) {
			gmrtfits_subint_real ( &gf, outfb, i, OGULP );
		}

		printf( " written \n" );
	}

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
	/*fclose ( outfile );*/

	return 0;
}

