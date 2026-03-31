#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include "gmrtfits.h"

// gpu
#include <cufft.h>

#ifdef TIMING
#include <time.h>
#endif

#define NPOL 4
#define NCHAN 2048
#define NSBLK 2048
#define NGULP 524288
// time integration
#define INTEGRATION 8
#define OGULP ( NGULP / INTEGRATION )
// ensure OGULP is integer multiple of NSBLK
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

__global__ 
void detect_integrate ( cufftComplex *in, cufftReal *out ) {
	/*
	 * in  is 2 x ( 1 + FFTCOMPLEXSIZE )
	 * out is NCHAN x OGULP x NPOL
	 * layout of out is (nsblk, npol, nchan)
	 *
	 * x axis is TIME
	 * y axis is FREQ
	 *
	 * dimensions output aligned
	 * <ogulp,nchan>
	 * <ogulp//16,nchans//16><16,16>
	 *
	 */

	// output aligned
	unsigned osamp  = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned ochan  = blockIdx.y * blockDim.y + threadIdx.y;
	
	unsigned long int ii = 0, kk = 0;

	ii    = NTAVG*osamp + NTAVG*OGULP*ichan;
	kk    = ochan + NCHAN*NPOL*osamp;

	__syncthreads();

	// definitions
	int    i = 0;
	float aa = 0.0, bb = 0.0, cr = 0.0, ci = 0.0;
	for ( i = 0; i < NTAVG; i++ ) {

		// load 2xNTAVG samples from pols
		// the pointer offset needed
		cufftComplex p1 = in [ 1 + ii + i ];
		cufftComplex p2 = in [ 1 + FFTCOMPLEXSIZE + ii + i ];

		aa += p1[0]*p1[0] + p1[1]*p1[1];
		bb += p2[0]*p2[0] + p2[1]*p2[1];
		cr += p1[0]*p2[0] + p1[1]*p2[1];
		ci += p1[1]*p2[0] - p1[0]*p2[1];
	}

	__syncthreads();

	// save in outfb
	out [ kk           ]  = aa;
	out [ kk +   NCHAN ]  = bb;
	out [ kk + 2*NCHAN ]  = cr;
	out [ kk + 3*NCHAN ]  = ci;
}

#ifdef TIMING
clock_t tstart, tstop;
double time_read, time_ffts, time_det, time_write;
double time_h2d, time_d2h;
#endif

int main() {
	print_help ();
	const char *pol1_infile_path = "/home/shining/C02pol1200MHz8bitsCRAB550MHz_10s.raw";
	const char *pol2_infile_path = "/home/shining/C02pol1200MHz8bitsCRAB550MHz_10s.raw";
	/*const char *toufile_path = "/tmp/baseband/testfb_float.raw";*/
	const char *oufile_path = "/tmp/baseband/testfs.fits";
	double mjd            = 60900.23741898148;

	// file pointers
	/*FILE *infile, *outfile;*/
	FILE *pol1_infile, *pol2_infile;
	FILE *tou1, *tou2, *tou3;
	/*FILE *toutfile;*/

	// variables
	unsigned long read_bytes;
	unsigned long infilesize1, infilesize2;

	unsigned int i;
	unsigned iread;
	unsigned nreads;
	unsigned ichan, isamp, osamp;

	// array mallocs
	float *volt_float;
	float *outfb;
	char *pol1_volt_read, *pol2_volt_read;

	//-------------------------//
	gmrtfits_t gf;
	gmrtfits_prepare ( &gf, oufile_path, mjd, NPOL, NCHAN, 550.0f, 200.0f, NSBLK, INTEGRATION, 8 );
	gmrtfits_create ( &gf );
	//-------------------------//

	// FFTW
	cufftResult cures;
	cudaError_t cuerr;
	/*cufftHandle pforward, pbackward1, pbackward2;*/
	cufftHandle pforward, pbackward;
	/*cudaStream_t stream   = NULL;*/
	cufftComplex *fdata_d;
	float *outfb_d;
	float *volt_d;

	dim3 diblock(16,16);
	dim3 digrid( OGULP/16, NCHAN/16 );

	// mallocs
	pol1_volt_read  = (char*) malloc ( VREAD * sizeof (char) );
	pol2_volt_read  = (char*) malloc ( VREAD * sizeof (char) );
	cudaHostMalloc ( &volt_float, 2 * VREAD * sizeof(float) );
	cudaMalloc ( &volt_d, 2 * VREAD * sizeof(float) );
	cudaMalloc ( &fdata_d, 2 * FFTCOMPLEXSIZE * sizeof(cufftComplex) );
	cudaHostMalloc ( &outfb, NPOL * NCHAN * OGULP * sizeof(float) );
	cudaMalloc ( &outfb_d, NPOL * OGULP * NCHAN * sizeof(float) );

	// files open
	pol1_infile  = fopen ( pol1_infile_path, "r" );
	pol2_infile  = fopen ( pol2_infile_path, "r" );

	/*tou1         = fopen ( "/tmp/baseband/testfb_data.raw" , "w" );*/
	/*tou2         = fopen ( "/tmp/baseband/testfb_scales.raw" , "w" );*/
	/*tou3         = fopen ( "/tmp/baseband/testfb_offsets.raw" , "w" );*/
	/*toutfile = fopen ( "/tmp/baseband/toutf32.raw", "w+" );*/

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
	nreads  = 2;

	int ng = VREAD;
	// forward FFT
	/*pforward = fftw_plan_many_dft_r2c ( 1, &ng, 2, volt_double, NULL, 1, VREAD, fdata, NULL, 1, FFTCOMPLEXSIZE, FFTW_ESTIMATE  );*/
	cures = cufftPlanMany ( &pforward, &ng, NULL, 1, VREAD, NULL, 1, FFTCOMPLEXSIZE, CUFFT_R2C, 2);


	// per pol :: real VREAD -> complex NCHAN*NGULP + 1
	// ignore the DC term
	/**
		CHECK 
	 * FFT R2C has the DC term as the first element.
	 * We ignore this DC term while coherent filterbanking.
	 * This messes up the array layout that forces us 
	 * to make two fftwplans for backward FFT
	 * ---
	 *  i am doing it twice with fftw, but one plan with cufft
	 **/
	/*pbackward1 = fftw_plan_many_dft ( 1, &ng, NCHAN, fdata+1, NULL, 1, NGULP, fdata+1, NULL, 1, NGULP, FFTW_BACKWARD, FFTW_ESTIMATE );*/
	/*pbackward2 = fftw_plan_many_dft ( 1, &ng, NCHAN, fdata+FFTCOMPLEXSIZE+1, NULL, 1, NGULP, fdata+FFTCOMPLEXSIZE+1, NULL, 1, NGULP, FFTW_BACKWARD, FFTW_ESTIMATE );*/
	ng    = NGULP;
	cures = cufftPlanMany (&pbackward, 1, &ng, NULL, 1, NGULP, NULL, 1, NGULP, CUFFT_C2C, NCHAN);

	printf (" prepared plans .. starting ... total nreads=%d\n", nreads);
	gmrtfits_subint_open ( &gf );

	for (iread = 0; iread < nreads; iread++) {
		printf( " iread=%d .. ", iread );

#ifdef TIMING
		tstart  = clock ();
#endif

		// read
		read_bytes = fread ( pol1_volt_read, sizeof(char), VREAD, pol1_infile );
		read_bytes = fread ( pol2_volt_read, sizeof(char), VREAD, pol2_infile );
		// if need to unpack
		// unpack here
		for (i = 0; i < VREAD; i++ ) {
			volt_float[i]         = (float) pol1_volt_read[i];
			volt_float[VREAD + i] = (float) pol2_volt_read[i];
		}

#ifdef TIMING
		tstop     = clock ();
		time_read = (tstop - tstart) / CLOCKS_PER_SEC;
#endif

		printf( " read .. " );

#ifdef TIMING
		tstart  = clock ();
#endif

		cures   = cudaMemcpy ( volt_d, volt_float, 2 * VREAD * sizeof(float), cudaMemcpyHostToDevice );

#ifdef TIMING
		tstop     = clock ();
		time_h2d  = (tstop - tstart) / CLOCKS_PER_SEC;
#endif

#ifdef TIMING
		tstart  = clock ();
#endif

		// forward FFT
		/*fftw_execute ( pforward );*/
		cures   = cufftExecR2C ( pforward, volt_d, fdata_d );
		if ( cures != CUFFT_SUCCESS ) {
			printf (" [!!] Forward FFT failed\n");
		}

		printf( " forward FFT .. " );

		// backward FFT
		// need to do two calls
		cures   = cufftExecC2C ( pbackward, fdata_d+1, fdata_d+1, CUFFT_INVERSE );
		if ( cures != CUFFT_SUCCESS ) {
			printf (" [!!] Backward FFT failed\n");
		}
		cures   = cufftExecC2C ( pbackward, fdata_d+FFTCOMPLEXSIZE+1, fdata_d+FFTCOMPLEXSIZE+1, CUFFT_INVERSE );
		if ( cures != CUFFT_SUCCESS ) {
			printf (" [!!] Backward FFT failed\n");
		}

		printf( " backward FFT .. " );

#ifdef TIMING
		tstop     = clock ();
		time_ffts = (tstop - tstart) / CLOCKS_PER_SEC;
#endif

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
		 * (ngulp, npol, nchan)
		 **/

#ifdef TIMING
		tstart  = clock ();
#endif
		detect_and_integrate <<<digrid,diblock>>> ( volt_d, outfb_d );
#ifdef TIMING
		tstop    = clock ();
		time_det = (tstop - tstart) / CLOCKS_PER_SEC;
#endif

		printf( " detection .. " );

#ifdef TIMING
		tstart  = clock ();
#endif

		cures   = cudaMemcpy ( outfb, outfb_d, OGULP * NCHAN * NPOL * sizeof(float), cudaMemcpyDeviceToHost );

#ifdef TIMING
		tstop     = clock ();
		time_d2h  = (tstop - tstart) / CLOCKS_PER_SEC;
#endif


		// write to file
		/*fwrite ( outfb, sizeof(float), NGULP * NCHAN * NPOL / INTEGRATION, toutfile );*/
		/*goto exit;*/

		/*
		 * for the case when OGULP > NSBLK
		 * possibly bug here
		 * when OGULP does not nicely divide 
		 * outfb gets zeroed out and that data is lost
		 * it isn't a bug when OGULP is an integer
		 * multiple of NSBLK
		 */
#ifdef TIMING
		tstart  = clock ();
#endif
		printf ("  writing subints=");
		for (i = 0; i < OGULP; i+=NSBLK) {
			printf ("%d ", i);
			gmrtfits_subint_real ( &gf, outfb, i, OGULP );
		}

#ifdef TIMING
		tstop      = clock ();
		time_write = (tstop - tstart) / CLOCKS_PER_SEC;
#endif

#ifdef TIMING
		printf ("[Timing] read=%.2f ffts=%.2f det=%.2f write=%.2f\n",time_read, time_ffts, time_det, time_write);
#endif

		//
		/*fwrite ( gf.data, sizeof(char), OGULP * NCHAN * NPOL, tou1 );*/
		/*fwrite ( gf.scales, sizeof(float), NCHAN * NPOL, tou2 );*/
		/*fwrite ( gf.offsets, sizeof(float), NCHAN * NPOL, tou3 );*/
		//
		/*goto exit;*/
	}

/*exit:*/
	gmrtfits_close ( &gf );

	// release memory
	cudaFreeHost ( volt_float );
	cudaFreeHost ( outfb );
	cudaFree ( volt_d );
	cudaFree ( fdata_d );
	cudaFree ( outfb_d );
	free ( pol1_volt_read );
	free ( pol2_volt_read );

	// destroy plans
	cures = cufftDestroy ( pforward );
	if ( cures != CUFFT_SUCCESS ) {
		printf (" [!!] Plan not destroyed\n");
	}
	cures = cufftDestroy ( pbackward );
	if ( cures != CUFFT_SUCCESS ) {
		printf (" [!!] Plan not destroyed\n");
	}

	// close files
	fclose ( pol1_infile );
	fclose ( pol2_infile );

	/*fclose ( tou1 );*/
	/*fclose ( tou2 );*/
	/*fclose ( tou3 );*/
	/*fclose ( toutfile );*/

	return 0;
}
