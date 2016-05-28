// example_floatbox.c - Atlee Brink
// A working example of all the functions in ../include/floatbox.h

#include "floatbox.h"

#include <stdio.h> /* for printf and stdout, not necessary otherwise */

int main()
{
    struct FLOATBOX box;

    // allocate heap memory for box
    if( !boxalloc( &box, p3d( 241, 241, 51 ) ) ) { /* handle error */ }

    // print box metadata
    boxfprint( stdout, "example: ", "\t", box );

    // put the value 4.567f to coordinates [1][2][3]
    const struct POINT3D pos = {1, 2, 3};
    const float putval = 4.567f;
    boxput( box, pos, putval );

    // get a value from coordinates [pos.x][pos.y][pos.z] and print it
    printf( "value at (%d, %d, %d): %g (should be %g)\n",
        pos.x, pos.y, pos.z, boxget( box, pos ), putval );

    const float fillval = 9.876f;
    int bad = 0;
    printf( "calling boxsetall(..)...\n" ); fflush( stdout );
    boxsetall( box, fillval );
    printf( "validating values...\n" ); fflush( stdout );
    for( int x = 0; !bad && x < box.size.x; x++ ) {
        for( int y = 0; !bad && y < box.size.y; y++ ) {
            for( int z = 0; !bad && z < box.size.z; z++ ) {
                if( boxget( box, p3d( x, y, z ) ) != fillval ) {
                    printf( "error: boxsetall(..) or boxget(..) went bad at (%d, %d, %d)!\n",
                        x, y, z );
                    bad = 1;
                }
            }
        }
    }
    if( !bad ) printf( "boxsetall(..) passed\n" );

    // free heap memory
    boxfree( &box );
}
