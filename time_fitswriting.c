/*
 *
 * Profile gmrtfits interface
 *
 * preliminary GPU based testing suggests
 * (as expected) FFT computations are not the bottleneck
 * but the disk I/O is, specifically the subint writing
 *
 * This is confirmed with this test here
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#include "gmrtfits.h"

#define NPOL 4
#define NCHAN 2048
#define NSBLK 2048
#define NGULP 65536
#define INTEGRATION 8
#define OGULP ( NGULP / INTEGRATION )


int main () {
	gmrtfits_t gf;
	gmrtfits_prepare ( &gf, "/tmp/baseband/prof.fits", 61122.223, NPOL, NCHAN, 550.0f, 200.0f, NSBLK, INTEGRATION, 8 );
	gmrtfits_create ( &gf );
	gmrtfits_subint_open ( &gf );

	float *fb    = (float*) malloc ( sizeof(float) * OGULP * NCHAN * NPOL );
	srand ( (unsigned int) time ( NULL ) );
	unsigned long nsize = OGULP * NCHAN * NPOL;
	for ( unsigned long it = 0; it < nsize; it++ ) {
		fb [ it ] =  1000.0f * ( (float)rand()/(float)(RAND_MAX) );
	}
	///
	int ntimes = 10;

	while ( ntimes-- ) {
		printf ("  writing subints\n");
		for (int i = 0; i < OGULP; i+=NSBLK) {
			//printf ("%d ", i);
			gmrtfits_subint_real ( &gf, fb, i, OGULP );
		}
	}
	///
	free ( fb );
	return 0;
}
