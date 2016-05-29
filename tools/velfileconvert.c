////////////////////////////////////////////////////////////////////////////////
// velfileconvert.c - Atlee Brink
// converts from old velocity file format to new (smaller) file format
// warning: loads entire volume into memory; there is a better way but I'm lazy
////////////////////////////////////////////////////////////////////////////////

#include "boxfiler.h"

#include <stdio.h>


const char signature[4] = {'v', 'b', 'o', 'x'};


int
main (
    int argc,
    char *argv[]
)
{
    if( argc != 3 ) {
        printf( "velfileconvert: velocity file converter\n" );
        printf( "usage: %s <in:oldfile.txt> <out:newfile.vbox>\n", argv[0] );
        return 0;
    }

    char *infile = argv[1];
    char *outfile = argv[2];

	struct FLOATBOX box;

    printf( "reading old velocity model %s...\n", infile ); fflush( stdout );
	if (!boxfileloadtext( &box, infile ) ) {
		fprintf( stderr, "%s: error: boxfileloadtext(..) failed to load %s\n",
			argv[0], infile );
		return 1;
	}
    printf( "\tdone.\n" ); fflush( stdout );

    printf( "writing new velocity model %s...\n", outfile ); fflush( stdout );
    if( !boxfilestorebinary( outfile, signature, box ) ) {
		fprintf( stderr, "%s: error: boxfilestorebinary(..) failed to store %s\n",
			argv[0], outfile );	
		return 1;
	}
    printf( "\tdone.\n" ); fflush( stdout );

    boxfree( &box );

    return 0;
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
