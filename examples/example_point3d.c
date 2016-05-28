// example_point3d.c - Atlee Brink
// Demonstrates and tests point3d.h.

#include "point3d.h"

#include <stdio.h>

int main() {

    struct POINT3D a = {1, 2, 3};          // can initialize like this
    struct POINT3D b = p3d(10, 11, 12);    // or can initialize like this
    struct POINT3D c;

    int val = 5;

    printf( "a: (%d, %d, %d)\n", a.x, a.y, a.z );
    printf( "b: (%d, %d, %d)\n", b.x, b.y, b.z );
    
    // adding two points
    c = p3daddp3d( a, b );
    printf( "a + b: (%d, %d, %d)\n", c.x, c.y, c.z );

    // adding a value to all components of a point
    c = p3daddval( a, val );
    printf( "a + %d: (%d, %d, %d)\n", val, c.x, c.y, c.z );

    // subtracting one point from another 
    c = p3dsubp3d( a, b );
    printf( "a - b: (%d, %d, %d)\n", c.x, c.y, c.z );

    // boundary comparison 'isless'
    printf( "is a < b ? %s\n", p3disless( a, b ) ? "yes" : "no" );

    // boundary comparison 'ismore'
    printf( "is a > b ? %s\n", p3dismore( a, b ) ? "yes" : "no" );

    // compute volume from dimensions
    printf( "volume of box with size (%d, %d, %d): %zu\n",
        b.x, b.y, b.z, p3dcalcvolume( b ) );

    // compute size of region from a to b, inclusive
    c = p3dsizeofregion( a, b );
    printf( "size of (%d, %d, %d) to (%d, %d, %d): (%d, %d, %d)\n",
        a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z );

    printf( "all done\n" );
}
