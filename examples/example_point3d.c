#include "point3d.h"

#include <stdio.h>

int main() {

    struct POINT3D a = {1, 2, 3};
    struct POINT3D b = {10, 11, 12};
    struct POINT3D c;
    int val = 5;

    printf( "a: (%d, %d, %d)\n", a.x, a.y, a.z );
    printf( "b: (%d, %d, %d)\n", b.x, b.y, b.z );
    
    c = p3daddp3d( a, b );
    printf( "a + b: (%d, %d, %d)\n", c.x, c.y, c.z );

    c = p3daddval( a, val );
    printf( "a + %d: (%d, %d, %d)\n", val, c.x, c.y, c.z );

    c = p3daddval( b, val );
    printf( "b + %d: (%d, %d, %d)\n", val, c.x, c.y, c.z );

    c = p3dsubp3d( a, b );
    printf( "a - b: (%d, %d, %d)\n", c.x, c.y, c.z );

    printf( "all done\n" );

}
