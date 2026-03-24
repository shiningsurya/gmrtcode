

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include "gmrtfits.h"

#define NCHAN 2048
#define NGULP 4096
#define NTAVG 8
#define OGULP ( NGULP / NTAVG )
#define OSIZE ( OGULP * NCHAN )
#define VREAD ( 2 * NCHAN * NGULP )
#define FFTCOMPLEXSIZE (1 + NCHAN*NGULP)

/*
 * 2 * nchan * ngulp * ts / ngulp / ntavg
=> 2 * nchan * ts / ntavg

ts = 1 / (2*bw)
=> 2 * nchan / ntavg / 2 / bw
=> nchan / bw / ntavg

actual formula is
fftint * nchan / bandwidth
 *
 */


int main() {
	const char *infile_path = "/tmp/baseband/C02pol1200MHz8bitsCRAB550MHz_10s.raw";
	/*const char *oufile_path = "/tmp/baseband/testfbfloat32.raw";*/
	const char *oufile_path = "/tmp/baseband/test.fits";

	/*FILE *infile, *outfile;*/
	FILE *infile;

	char *volt_read = (char*) malloc ( VREAD * sizeof (char) );
	double *volt_double = (double*) malloc ( VREAD * sizeof(double) );
	float *outfb = (float*)malloc ( OSIZE * sizeof(float) );
	// fastest axis is freq
	// TF
	//
	//-------------------------//
	gmrtfits_t gf;
	gmrtfits_prepare ( &gf, oufile_path, 1, NCHAN, 550.0f, 200.0f, OGULP, NTAVG, 8 );
	gmrtfits_create ( &gf );
	//-------------------------//


	fftw_plan pforward, pbackward;
	fftw_complex *fdata, *fdata_off1;

	unsigned long read_bytes;
	unsigned long infilesize;

	unsigned iread;
	unsigned nreads;
	unsigned ichan, isamp, iavg;

	/*fdata = (fftw_complex*) fftw_malloc  ( (1 + NCHAN*NGULP ) * sizeof(fftw_complex) );*/
	/*size_t ac = FFTCOMPLEXSIZE;*/
	/*printf(" allocating size=%zu\n", ac);*/

	/*fdata = (fftw_complex*) fftw_alloc_complex  ( ac );*/
	fdata = (fftw_complex*) fftw_alloc_complex  ( FFTCOMPLEXSIZE );
	fdata_off1 = fdata + 1;

	// files open
	infile  = fopen ( infile_path, "r" );
	/*outfile = fopen ( oufile_path, "w+" );*/

	// get file size
	fseek ( infile, 0L, SEEK_END );
	infilesize = ftell ( infile );
	fseek ( infile, 0L, SEEK_SET );
	printf(" filesize=%lu\n", infilesize);

	nreads  = infilesize / VREAD;
	/*nreads  = 4;*/

	// forward FFT
	pforward = fftw_plan_dft_r2c_1d ( VREAD, volt_double, fdata, FFTW_ESTIMATE );
	// real VREAD -> complex NCHAN*NGULP + 1
	// ignore the DC term
	const int ng = NGULP;
	pbackward = fftw_plan_many_dft ( 1, &ng, NCHAN, fdata_off1, NULL, 1, NGULP, fdata_off1, NULL, 1, NGULP, FFTW_BACKWARD, FFTW_ESTIMATE );

	printf (" prepared plans .. starting ...\n");
	gmrtfits_subint_open ( &gf );

	for (iread = 0; iread < nreads; iread++) {
		printf( " iread=%d .. ", iread );
		// read
		read_bytes = fread ( volt_read, sizeof(char), VREAD, infile );
		for ( int i = 0; i < VREAD; i++ ) volt_double[i] = (double) volt_read[i];

		printf( " read .. " );

		// forward FFT
		fftw_execute ( pforward );

		printf( " forward FFT .. " );

		// backward FFT
		fftw_execute ( pbackward );

		printf( " backward FFT .. " );

		// detection
		unsigned long il = 0, jl = 0, kl = 0;
		float fil = 0.0f;
		for (ichan = 0; ichan < NCHAN; ichan++) {
			for (isamp = 0; isamp < OGULP; isamp++) {
				il = ichan*NGULP + isamp*NTAVG;
				kl = ichan + isamp*NCHAN;
				// time-averaging
				outfb[ kl ]  = 0.0;
				for ( iavg = 0; iavg < NTAVG; iavg++ ) {
					jl = il + iavg;
					fil = fdata_off1[jl][0]*fdata_off1[jl][0] + fdata_off1[jl][1]*fdata_off1[jl][1];
					//outfb [ ichan*NCHAN + isamp ] += fil;
					outfb [ kl ] += fil;
				}
			}
		}

		printf( " detection .. " );

		// write to file
		/*fwrite ( outfb, sizeof(float), OGULP * NCHAN, outfile );*/
		gmrtfits_subint_real_single ( &gf, outfb );

		printf( " written \n" );
	}

	gmrtfits_close ( &gf );

	// destroy plans
	fftw_destroy_plan ( pforward );
	fftw_destroy_plan ( pbackward );

	fftw_free ( fdata );
	free ( volt_read );
	free (volt_double);
	free (outfb);

	fclose ( infile );
	/*fclose ( outfile );*/

	return 0;
}
