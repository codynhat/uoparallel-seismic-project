////////////////////////////////////////////////////////////////////////////////
// example_floatbox.c - Atlee Brink
// A thorough test of all functions in floatbox.h
////////////////////////////////////////////////////////////////////////////////

#include "floatbox.h"
#include "point3d.h"

#include <stdio.h> // for printf and stdout, not necessary otherwise
#include <stdlib.h>

int main()
{
    struct FLOATBOX box;

    struct POINT3D omin = {1, -5, 2};
    struct POINT3D omax = {241, 241, 51};
    struct POINT3D imin = {121, 121, 5};
    struct POINT3D imax = {241, 241, 51};

    // prepare the box
    puts( "allocating box..." );
    if( !boxalloc( &box, omin, omax, imin, imax ) ) {
        puts( "\tboxalloc(..) failed!" );
    }
    puts( "\tsuccess!" );

    puts("");

    // show the box
    boxfprint( stdout, "box: ", "  ", box );

    puts("");

    // test boxsetall(..)
    puts( "testing boxsetall(..):" );
    {
        const float val = 1.7;
        puts( "\tcalling boxsetall(..)..." );
        boxsetall( box, val );
        puts( "\tchecking values..." );
        for( long i = 0; i <= box.offset.m; i++ ) {
            if( box.flat[i] != val ) {
                printf( "boxsetall: expected != stored at index %lu\n", i );
                exit(0);
            }
        }
    }
    puts( "\tsuccess!" );

    puts("");

    // give each element a unique value so we can check other functions
    puts( "directly setting values of box.flat[]..." );
    for( long i = 0; i <= box.offset.m; i++ ) {
        box.flat[i] = (float)i;
    }
    puts( "\tdone" );

    puts("");

    // test boxgetlocal(..), which goes from (0,0,0) to size-1
    puts( "testing boxgetlocal(..) and boxputlocal(..)..." );
    {
        long i = 0;
        for( int x = 0; x < box.size.x; x++ ) {
            for( int y = 0; y < box.size.y; y++ ) {
                for( int z = 0; z < box.size.z; z++ ) {
                    const float expected = (float)i;
                    const struct POINT3D pt = {x, y, z};
                    float stored;
                    // boxgetlocal
                    stored = boxgetlocal( box, pt );
                    if( expected != stored ) {
                        printf( "boxgetlocal: expected != stored at local (%d, %d, %d)\n",
                            x, y, z );
                        exit(0);
                    }
                    // boxputlocal
                    boxputlocal( box, pt, -expected );
                    stored = boxgetlocal( box, pt );
                    if( -expected != stored ) {
                        printf( "boxputlocal: expected != stored at local (%d, %d, %d)\n",
                            x, y, z );
                        exit(0);
                    }
                    // reset value
                    boxputlocal( box, pt, expected );
                    i++;
                }
            }
        }
    }
    puts( "\tsuccess!" );

    puts("");

    // test boxgetglobal(..), which goes from omin to omax
    puts( "testing boxgetglobal(..) and boxputglobal(..)..." );
    {
        long i = 0;
        for( int x = box.omin.x; x <= box.omax.x; x++ ) {
            for( int y = box.omin.y; y <= box.omax.y; y++ ) {
                for( int z = box.omin.z; z <= box.omax.z; z++ ) {
                    const float expected = (float)i;
                    const struct POINT3D pt = {x, y, z};
                    float stored;
                    // boxgetglobal
                    stored = boxgetglobal( box, pt );
                    if( expected != stored ) {
                        printf( "boxgetglobal: expected != stored at local (%d, %d, %d)\n",
                            x, y, z );
                        exit(0);
                    }
                    // boxputglobal
                    boxputglobal( box, pt, -expected );
                    stored = boxgetglobal( box, pt );
                    if( -expected != stored ) {
                        printf( "boxputglobal: expected != stored at local (%d, %d, %d)\n",
                            x, y, z );
                        exit(0);
                    }
                    // reset value
                    boxputglobal( box, pt, expected );
                    i++;
                }
            }
        }
    }
    puts( "\tsuccess!" );

    puts("");

    // test boxgetinner(..), which goes from (0,0,0) to (imax - imin + 1)
    // NOTE: it should still work for negative numbers
    puts( "testing boxgetinner(..) and boxputinner(..)..." );
    {
        long i = 0;
        for( int x = 0; x < box.size.x; x++ ) {
            for( int y = 0; y < box.size.y; y++ ) {
                for( int z = 0; z < box.size.z; z++ ) {
                    const float expected = (float)i;
                    const struct POINT3D pt = p3dsubp3d( p3d(x, y, z), box.imin );
                    float stored;
                    // boxgetinner
                    stored = boxgetinner( box, pt );
                    if( expected != stored ) {
                        printf( "boxgetinner: expected != stored at local (%d, %d, %d)\n",
                            x, y, z );
                        exit(0);
                    }
                    // boxputinner
                    boxputinner( box, pt, -expected );
                    stored = boxgetinner( box, pt );
                    if( -expected != stored ) {
                        printf( "boxputinner: expected != stored at local (%d, %d, %d)\n",
                            x, y, z );
                        exit(0);
                    }
                    // reset value
                    boxputinner( box, pt, expected );
                    i++;
                }
            }
        }
    }
    puts( "\tsuccess!" );

    puts("");

    puts( "freeing box..." );
    boxfree( &box );
    puts( "all done" );

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
