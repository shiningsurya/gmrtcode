#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "gmrtfits.h"

#define NSBLK 2048

void print_help() {
	fprintf(stderr, "\n");
	fprintf(stderr, "  raw2sf - Convert GMRT GWB raw data into search mode PSRFITS\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "    -h Print this help\n");
	fprintf(stderr, "    -n <nchan> Number of channels in the raw file\n");
	fprintf(stderr, "    -f <freq> Frequency edge of recording in MHz\n");
	fprintf(stderr, "    -b <fbw> Signed bandwidth in MHz where sign determines USB(positive) or LSB(negative)\n");
	fprintf(stderr, "    -t <mjd> MJD in maximum available precision\n");
	fprintf(stderr, "    -o sf fits filename\n");
	fprintf(stderr, "\n");
	fprintf(stderr, " Notes:\n");
	fprintf(stderr, "    Currently set to work with 16bit filterbank raws.\n");
	fprintf(stderr, "    Extending it to 8bit is straightforward but not considered necessary now.\n");
	fprintf(stderr, "    Defaults to circular basis.\n");
}

int main (int argc, char *argv[]) {
	int opt;

	// args
	double mjd;
	float fedge, bw;
	unsigned nchan;

	// files
	char *infile;
	char *oufile;
	char *chunk;
	FILE *inraw;
	gmrtfits_t gf;

	while ( (opt = getopt ( argc, argv, "hn:f:b:t:o:" )) != -1 ) {
		switch (opt) {
			case 'h':
				print_help ();
				exit (EXIT_SUCCESS);
				break;
			case 't':
				mjd    = atof ( optarg );
				break;
			case 'f':
				fedge  = atof ( optarg );
				break;
			case 'b':
				bw     = atof ( optarg );
				break;
			case 'o':
				oufile = strdup ( optarg );
				break;
		} // switch
	} // loop
	
	// 


	free ( infile );
	free ( oufile );
	free ( chunk );

	return 0;
}
