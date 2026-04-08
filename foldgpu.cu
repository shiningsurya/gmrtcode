/*
 * foldgpu - Use GPU to detect voltages into coherence products 
 * and directly fold using either the given period or polyco file.
 * If given, also calibrate the voltages using jones matrix.
 */
/* #include "foldgpu.hpp" */
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cmath>
#include <time.h>

/* gpu */
#include <cufft.h>

constexpr int NPOL                     = 4;
constexpr int NCHAN                    = 2048;
constexpr int NBIN                     = 4096;
constexpr unsigned NGULP               = 65536;
constexpr unsigned long VREAD          = 2 * NCHAN * NGULP;
constexpr unsigned long FFTCOMPLEXSIZE = 1 + NCHAN*NGULP;

constexpr double TAU        = 6.283185307179586;
constexpr double DMCONSTANT = 2.41E-10;


void print_help() {
	fprintf(stderr, "\n");
	fprintf(stderr, "  foldgpu - Transforms 8bit GMRT GWB baseband data into fold mode PSRFITS\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "    -h Print this help\n");
	fprintf(stderr, "    -o <sf> output PSRFITS file\n");
	fprintf(stderr, "    -r <raw> Right Circular Polarized baseband file\n");
	fprintf(stderr, "    -l <raw> Left Circular Polarized baseband file\n");
	fprintf(stderr, "    -f <freq> Frequency edge of recording in MHz\n");
	fprintf(stderr, "    -b <fbw> Signed bandwidth in MHz where sign determines USB(positive) or LSB(negative)\n");
	fprintf(stderr, "    -t <mjd> Start MJD in maximum available precision\n");
	fprintf(stderr, "    -d <dm> Dispersion measure to perform coherent de-dispersion\n");
	fprintf(stderr, "    -c <period> period to use while folding\n");
	fprintf(stderr, "    -p <polyco> polyco file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " Notes:\n");
	fprintf(stderr, " *  First polarization is RCP (-r) and second is LCP (-l)\n");
	fprintf(stderr, "    This causes V defined as RR - LL.\n");
	fprintf(stderr, "\n");
}

__global__ void dedisperse ( cufftComplex *in, cufftComplex *out, float fedge, float bw, float dm ) {
	/*
	 * in is [2xFFTCOMPLEXSIZE] of {pol1 ... pol2}
	 * ou is [2xFFTCOMPLEXSIZE] of {pol1 ... pol2}
	 *
	 * per pol, R2C is FFTCOMPLEXSIZE=1+NCHAN*NGULP
	 * ignore DC term
	 * remaining=NCHAN*NGULP
	 * 
	 * DC | fedge | -------- > | fedge + bw
	 * 0  |    1  | ...........| NCHAN*NGULP
	 *
	 * ichan should go from 
	 *    |    0  | ...........| NCHAN*NGULP - 1
	 *
	 *
	 *
	 */
	
	size_t ivolt = blockIdx.x * blockDim.x + threadIdx.x;

	int ichan    = (ivolt > FFTCOMPLEXSIZE) ? (ivolt - FFTCOMPLEXSIZE - 1) : ( ivolt - 1 ) ;

	//float fcen  = fedge + (0.5 * bw);
	double fbw  = bw / NCHAN / NGULP;
	float f     = ichan * fbw;
	float fecen = fedge + ( 0.5 * fbw );

	float phase  = -1.0 * f * f * dm * TAU * / DMCONSTANT / fecen / fecen / ( fecen + f );
	float cphase = cosf ( phase ); 
	float sphase = sinf ( phase );

	cufftComplex vin = in  [ ivolt ];
	cufftComplex vou = out [ ivolt ];

	vou.x    = ( vin.x * cphase ) - ( vin.y * sphase );
	vou.y    = ( vin.x * sphase ) + ( vin.y * cphase );
}

__global__ void detect_folder ( cufftComplex *in, double *out, double f0, double phase0, double tsamp ) {
	/*
	 * in  is 2 x ( 1 + FFTCOMPLEXSIZE )
	 * out is NBIN x NCHAN x NPOL
	 * layout of in is {1, NCHAN*NGULP, 1, NCHAN*NGULP}
	 * layout of out is (nbin, nchan, npol)
	 *
	 * dimensions of threadlayout: (freq, time)
	 * matches the in array
	 *
	 * x axis is TIME
	 * y axis is FREQ
	 *
	 * <ogulp,nchan>
	 * <ogulp//16,nchans//16><16,16>
	 *
	 * - replace (f0,phase0,tsamp) with binplan
	 *  which will be another device array
	 *
	 */

	/* thread dimensions */
	unsigned isamp  = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned ichan  = blockIdx.y * blockDim.y + threadIdx.y;

	/* fft indices */
	size_t iv1, iv2;
	iv1   = isamp + NGULP*ichan;
	iv2   = FFTCOMPLEXSIZE + ( isamp + NGULP*ichan; );

	/* fold indices */
	double phase  = modf (phase0 + ( f0 * isamp * tsamp ), 1.0);
	unsigned ibin = (unsigned) (phase * NBIN);
	size_t ob;
	/* (nbin, nchan, npol) */
	ob    = NPOL*ichan + NPOL*NCHAN*ibin;

	/* detection */
	cufftComplex p1 = in [ 1 + iv1 ];
	cufftComplex p2 = in [ 1 + iv2 ];
	float aa = 0., bb = 0., cr = 0., ci = 0.;

	aa = p1.x*p1.x + p1.y*p1.y;
	bb = p2.x*p2.x + p2.y*p2.y;
	cr = p1.x*p2.x + p1.y*p2.y;
	ci = p1.y*p2.x - p1.x*p2.y;

	/* folding */
	out [ ob ]     += aa;
	out [ ob + 1 ] += bb;
	out [ ob + 2 ] += cr;
	out [ ob + 3 ] += ci;
}

int main( int argc, char *argv[]) {
	// arguments
	int opt;
	int gid = 0;
	double mjd = -1.0;
	float fedge = -1.0, bw = 0.0;
	float dm = 0.0;
	float fold_period = 0.0;
	string pol1_infile_path, pol2_infile_path;
	string oufile_path;
	string polyco_file;
	//char *pol1_infile_path = nullptr, *pol2_infile_path = nullptr;
	//char *oufile_path = nullptr;
	//char *polyco_file = nullptr;

	while ( (opt = getopt ( argc, argv, "hr:l:f:b:t:o:c:p:d:g:" )) != -1 ) {
		switch (opt) {
			case 'h':
				print_help ();
				exit (EXIT_SUCCESS);
				break;
			case 't':
				mjd    = atof ( optarg );
				break;
			case 'g':
				gid    = atoi ( optarg );
				break;
			case 'f':
				fedge  = atof ( optarg );
				break;
			case 'd':
				dm  = atof ( optarg );
				break;
			case 'c':
				fold_period  = atof ( optarg );
				break;
			case 'b':
				bw     = atof ( optarg );
				break;
			case 'l':
				//pol2_infile_path = strdup ( optarg );
				pol2_infile_path = optarg;
				break;
			case 'r':
				//pol1_infile_path = strdup ( optarg );
				pol1_infile_path = optarg;
				break;
			case 'o':
				//oufile_path      = strdup ( optarg );
				oufile_path      =  optarg;
				break;
			case 'p':
				//oufile_path      = strdup ( optarg );
				polyco_file      =  optarg;
				break;
		} // switch
	} // arg loop
		
	// sanitizing arguments
	if ( mjd < 0.0 ) {
		fprintf(stderr, "MJD not set .. exiting .. ");
		exit (EXIT_SUCCESS);
	}
	if ( fedge < 0.0 ) {
		fprintf(stderr, "fedge not set .. exiting .. ");
		exit (EXIT_SUCCESS);
	}
	if ( bw == 0.0 ) {
		fprintf(stderr, "bw not set .. exiting .. ");
		exit (EXIT_SUCCESS);
	}
	if ( fold_period == 0.0 && polyco_file.length() == 0) {
		fprintf(stderr, "period cannot be computed .. exiting .. ");
		exit (EXIT_SUCCESS);
	}
	if ( pol1_infile_path.length() == 0 || pol2_infile_path.length() == 0 ) {
		fprintf(stderr, "raw input files missing .. exiting .. ");
		exit (EXIT_SUCCESS);
	}
	if ( oufile_path.length() == 0 ) {
		fprintf(stderr, "output file missing .. exiting .. ");
		exit (EXIT_SUCCESS);
	}


	/* file io */
	auto pol1_infile  = std::ifstream ( pol1_infile_path, std::ios::in | std::ios::binary );
	auto pol2_infile  = std::ifstream ( pol2_infile_path, std::ios::in | std::ios::binary );
	
	/* file size check */
	pol1_infile.seekg ( 0, std::ios::end );
	auto infilesize1 = pol1_infile.tellg ();
	pol1_infile.seekg ( 0, std::ios::beg );

	pol2_infile.seekg ( 0, std::ios::end );
	auto infilesize2 = pol2_infile.tellg ();
	pol2_infile.seekg ( 0, std::ios::beg );

	if ( infilesize1 != infilesize2 ) {
		fprintf(stderr, " filesizes pol1=%lu pol2=%lu\n", infilesize1, infilesize2);
		fprintf(stderr, " Unexpected case.");
		exit (EXIT_SUCCESS);
	}

	/* main group of variables */
	unsigned i;
	unsigned nreads;
	unsigned iread;

	/* main allocations */
	// std::array<char,VREAD> pol1_volt_read;
	// std::array<char,VREAD> pol2_volt_read;
	std::vector<char> pol1_volt_read ( VREAD );
	std::vector<char> pol2_volt_read ( VREAD );

	/* set cuda device */
	cudaSetDevice ( gid );

	cufftResult cures;
	cudaError_t cuerr;
	cufftHandle pforward, pbackward;

	/* cuda pointers */
	int ng     = VREAD;
	float *volt_float, *volt_d;
	cufftComplex *fdata_d;
	double *pdata, *pdata_d;
	volt_d     = nullptr;
	volt_float = nullptr;
	fdata_d    = nullptr;
	pdata      = nullptr;
	pdata_d    = nullptr;

	/* forward FFT */
	cures = cufftPlanMany ( &pforward, 1, &ng, NULL, 1, VREAD, NULL, 1, FFTCOMPLEXSIZE, CUFFT_R2C, 2);
	if ( cures != CUFFT_SUCCESS ) {
			printf (" [!!] Forward FFT plan failed code=%d\n", cures);
			goto exit;
	}
	/* backward FFT */
	ng    = NGULP;
	cures = cufftPlanMany (&pbackward, 1, &ng, NULL, 1, NGULP, NULL, 1, NGULP, CUFFT_C2C, NCHAN);
	if ( cures != CUFFT_SUCCESS ) {
			printf (" [!!] Backward FFT plan failed code=%d\n", cures);
			goto exit;
	}

	/* cuda mallocs */
	// asking thrust to manage this 
	// is ugly
	cudaMalloc ( &volt_d, 2 * VREAD * sizeof(float) );
	if ( volt_d == nullptr ) {
		fprintf(stderr, "cudaMalloc of volt_d failed\n");
		goto exit;
	}
	cudaHostAlloc ( &volt_float, 2 * VREAD * sizeof(float), cudaHostAllocDefault );
	if ( volt_float == nullptr ) {
		fprintf(stderr, "cudaHostAlloc of volt_float failed\n");
		goto exit;
	}

	cudaMalloc ( &fdata_d, 2 * FFTCOMPLEXSIZE * sizeof(cufftComplex) );
	if ( fdata_d == nullptr ) {
		fprintf(stderr, "cudaMalloc of fdata_d failed\n");
		goto exit;
	}

	/* fold-subint */
	size_t pdata_size = NBIN * NCHAN * NPOL * sizeof(double);
	cudaHostAlloc ( &pdata, pdata_size, cudaHostAllocDefault );
	if ( pdata == nullptr ) {
		fprintf(stderr, "cudaHostAlloc of pdata failed\n");
		goto exit;
	}
	cudaMalloc ( &pdata_d,  pdata_size );
	if ( pdata_d == nullptr ) {
		fprintf(stderr, "cudaMalloc of pdata_d failed\n");
		goto exit;
	}

	/* GPU thread layout */
	dim3 dd_block (256);
	dim3 dd_grid ( 2 * NCHAN * NGULP / 256);
	dim3 fold_block (16,16);
	dim3 fold_grid ( NGULP/16, NCHAN/16 );

	/*
	 * Main work loop
	 */
	for (iread = 0; iread < nreads; iread++) {
		printf( " iread=%d .. ", iread );

#ifdef TIMING
		tstart  = clock ();
#endif

		// read
		pol1_infile.read ( pol1_volt_read.data(), VREAD );
		pol2_infile.read ( pol2_volt_read.data(), VREAD );
		//read_bytes = fread ( pol1_volt_read, sizeof(char), VREAD, pol1_infile );
		//read_bytes = fread ( pol2_volt_read, sizeof(char), VREAD, pol2_infile );
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

		cuerr   = cudaMemcpy ( volt_d, volt_float, 2 * VREAD * sizeof(float), cudaMemcpyHostToDevice );

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
			printf (" [!!] Forward FFT failed code=%d\n", cures);
		}

		printf( " forward FFT .. " );

		// de-dispersion
		// doit with kernel
		dedisperse<<<dd_grid,dd_block>>> ( fdata_d, fdata_d, fedge, bw, dm );

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
		/*
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
		 */

		/* cuda memset */
		cudaMemset ( (void*) pdata_d, 0, pdata_size );

#ifdef TIMING
		tstart  = clock ();
#endif
		detect_folder<<<fold_grid,fold_block>>> ( fdata_d, pdata_d, f0, phase0, tsamp );
#ifdef TIMING
		tstop    = clock ();
		time_det = (tstop - tstart) / CLOCKS_PER_SEC;
#endif

		printf( " detection .. " );

#ifdef TIMING
		tstart  = clock ();
#endif

		cuerr   = cudaMemcpy ( pdata, pdata_d, pdata_size, cudaMemcpyDeviceToHost );

#ifdef TIMING
		tstop     = clock ();
		time_d2h  = (tstop - tstart) / CLOCKS_PER_SEC;
#endif


		/* write to file */
		/*fwrite ( outfb, sizeof(float), NGULP * NCHAN * NPOL / INTEGRATION, toutfile );*/
		/*goto exit;*/

#ifdef TIMING
		tstart  = clock ();
#endif
		//gmrtfits_subint_add ( &gf, outfb, OGULP );
		printf ("  writing subints\n");

#ifdef TIMING
		tstop      = clock ();
		time_write = (tstop - tstart) / CLOCKS_PER_SEC;
#endif

#ifdef TIMING
		printf ("[Timing] read=%.2f ffts=%.2f det=%.2f write=%.2f copies=%.2f\n",time_read, time_ffts, time_det, time_write, time_h2d + time_d2h);
#endif

		//
		/*fwrite ( gf.data, sizeof(char), OGULP * NCHAN * NPOL, tou1 );*/
		/*fwrite ( gf.scales, sizeof(float), NCHAN * NPOL, tou2 );*/
		/*fwrite ( gf.offsets, sizeof(float), NCHAN * NPOL, tou3 );*/
		//
		/*goto exit;*/
	}

exit:
	if ( pdata_d )    cudaFree ( pdata_d );
	if ( pdata )      cudaFreeHost ( pdata );
	if ( fdata_d )    cudaFree ( fdata_d );
	if ( volt_float ) cudaFreeHost ( volt_float );
	if ( volt_d )     cudaFree ( volt_d );

	/* destroy cufft plans */
	cures = cufftDestroy ( pforward );
	if ( cures != CUFFT_SUCCESS ) {
		printf (" [!!] Plan not destroyed\n");
	}
	cures = cufftDestroy ( pbackward );
	if ( cures != CUFFT_SUCCESS ) {
		printf (" [!!] Plan not destroyed\n");
	}

	return 0;
}
