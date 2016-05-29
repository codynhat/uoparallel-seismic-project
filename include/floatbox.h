////////////////////////////////////////////////////////////////////////////////
// floatbox.h - 2016.05.06 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
//
// Requires:
//   point3d.h
//
//
// Data:
//   struct FLOATBOX
//
//
// Functions:
//   boxalloc
//   boxfree
//
//   boxindex
//   boxget
//   boxput
//   boxsetall
//
//   boxfprint
//
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "point3d.h"


#include <stdio.h>
#include <stdlib.h>


////////////////////////////////////////////////////////////////////////////////
// structs
////////////////////////////////////////////////////////////////////////////////

struct FLOATBOX {
    struct POINT3D size;               // total dimensions
    struct POINT3D omin, omax;         // outer minimum, maximum in box
    struct POINT3D imin, imax;         // inner minimum, maximum in box
    struct { long x, y, z; } stride;   // for quickly computing indices
    struct {
        long o;                        // hypothetical index of omin
        long i;                        // hypothetical index of imin
        long m;                        // actual index of size-1 (max index)
    } offset;
    float *flat;                       // [x][y][z] order
};


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

inline extern
long
boxindex (
    const struct FLOATBOX box,
    const struct POINT3D pt
)
// return an index into box.flat corresponding to the given coordinates
{
    return box.stride.x * pt.x + box.stride.y * pt.y + box.stride.z * pt.z;
}


int
boxalloc (
    struct FLOATBOX *box,
    const struct POINT3D omin,
    const struct POINT3D omax,
    const struct POINT3D imin,
    const struct POINT3D imax
)
// allocates memory for and prepares contents of 'box'
// return: 0 on error, non-0 on success
{
    struct POINT3D size = p3dsizeofregion( omin, omax );
    long numbytes = sizeof(*box->flat) * p3dcalcvolume( size );

    float *flat = malloc( numbytes );
    if( flat == NULL ) return 0;

    box->size = size;
    box->omin = omin;
    box->omax = omax;
    box->imin = imin;
    box->imax = imax;

    box->stride.x = (long)size.y * size.z;
    box->stride.y = size.z;
    box->stride.z = 1;

    box->offset.o = boxindex( *box, omin );
    box->offset.i = boxindex( *box, imin );
    box->offset.m = boxindex( *box, p3dsubp3d( omax, omin) );

    box->flat = flat;

    return 1;
}


void
boxfree (
    struct FLOATBOX *box
)
// releases heap memory associated with box
{
    if( box == NULL ) return;
    free( box->flat );
    box->flat = NULL;
}


inline extern
float
boxgetlocal (
    const struct FLOATBOX box,
    const struct POINT3D pt
)
// coordinates relative to total box (from 0,0,0 to size-1)
{
    return box.flat[ boxindex( box, pt ) ];
}


inline extern
void
boxputlocal (
    const struct FLOATBOX box,
    const struct POINT3D pt,
    const float val
)
// coordinates relative to total box (starts at 0,0,0)
{
    box.flat[ boxindex( box, pt ) ] = val;
}


inline extern
float
boxgetglobal (
    const struct FLOATBOX box,
    const struct POINT3D pt
)
// global coordinates (from box.omin to box.omax)
{
    return box.flat[ boxindex( box, pt ) - box.offset.o ];
}


inline extern
void
boxputglobal (
    const struct FLOATBOX box,
    const struct POINT3D pt,
    const float val
)
// global coordinates (from box.omin to box.omax)
{
    box.flat[ boxindex( box, pt ) - box.offset.o ] = val;
}


inline extern
float
boxgetinner (
    const struct FLOATBOX box,
    const struct POINT3D pt
)
// inner coordinates (from 0,0,0 to imax-imin+1)
{
    return box.flat[ boxindex( box, pt ) + box.offset.i ];
}


inline extern
float
boxputinner (
    const struct FLOATBOX box,
    const struct POINT3D pt,
    const float val
)
// inner coordinates (from 0,0,0 to imax-imin+1)
{
    box.flat[ boxindex( box, pt ) + box.offset.i ] = val;
}


void
boxsetall (
    const struct FLOATBOX box,
    const float val
)
// sets ALL values in the volume to the given value
{
    if( box.flat == NULL ) return;
    for( long i = p3dcalcvolume( box.size ); i-- > 0; box.flat[i] = val );
}


void
boxfprint (
    FILE *stream,
    const char *prefix,
    const char *indent,
    const struct FLOATBOX box
)
// metadata friendly-printing
{
    if( stream == NULL ) stream = stdout;
    if( prefix == NULL ) prefix = "";
    if( indent == NULL ) indent = "  ";

    fprintf( stream, "%sFLOATBOX {\n", prefix );
    fprintf( stream, "%s%ssize: (%d, %d, %d)\n",
        prefix, indent,
        box.size.x, box.size.y, box.size.z );
    fprintf( stream, "%s%somin: (%d, %d, %d), omax: (%d, %d, %d)\n",
        prefix, indent,
        box.omin.x, box.omin.y, box.omin.z,
        box.omax.x, box.omax.y, box.omax.z );
    fprintf( stream, "%s%simin: (%d, %d, %d), imax: (%d, %d, %d)\n",
        prefix, indent,
        box.imin.x, box.imin.y, box.imin.z,
        box.imax.x, box.imax.y, box.imax.z );
    fprintf( stream, "%s%sstride: (%zu, %zu, %zu)\n",
        prefix, indent,
        box.stride.x, box.stride.y, box.stride.z );
    fprintf( stream, "%s%soffset.o: %zu, offset.i: %zu, offset.m: %zu\n",
        prefix, indent,
        box.offset.o, box.offset.i, box.offset.m );
    fprintf( stream, "%s%sflat: %p\n",
        prefix, indent, (void*)box.flat );
    fprintf( stream, "%s}\n", prefix );
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
