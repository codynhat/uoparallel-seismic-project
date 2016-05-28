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
    size_t sx, sy, sz; // array strides in 'flat'
    struct POINT3D size; // (x,y,z) dimensions
    float *flat; // [x][y][z] order
};


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

int
boxalloc (
    struct FLOATBOX *box,
    const struct POINT3D size
)
// allocates box->flat and sets strides and dimensions appropriately;
// on error: returns 0 (failure to allocate memory)
// on success: returns non-zero
{
    size_t numbytes = sizeof(*box->flat) * p3dcalcvolume( size );
    float *flat = malloc( numbytes );

    if( flat == NULL ) return 0;

    box->sx = (size_t)size.y * size.z;
    box->sy = size.z;
    box->sz = 1;

    box->size = size;

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
    box->size = p3d(0, 0, 0);
    free( box->flat );
    box->flat = NULL;
}


inline extern
size_t
boxindex (
    const struct FLOATBOX box,
    const struct POINT3D pt
)
// return an index into box.flat corresponding to the given coordinates
{
    return box.sx * pt.x + box.sy * pt.y + box.sz * pt.z;
}


inline extern
float
boxget (
    const struct FLOATBOX box,
    const struct POINT3D pt
)
// returns a single value from the given coordinates
{
    return box.flat[ boxindex( box, pt ) ];
}


inline extern
void
boxput (
    const struct FLOATBOX box,
    const struct POINT3D pt,
    const float val
)
// stores a single value at the given coordinates
{
    box.flat[ boxindex( box, pt ) ] = val;
}


void
boxsetall (
    const struct FLOATBOX box,
    const float val
)
// sets ALL values in the volume to the given value
{
    if( box.flat == NULL ) return;
    for( size_t i = p3dcalcvolume( box.size ); i-- > 0; box.flat[i] = val );
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
    fprintf( stream, "%s%sstrides: (%zu, %zu, %zu)\n",
        prefix, indent, box.sx, box.sy, box.sz );
    fprintf( stream, "%s%ssize: (%d, %d, %d)\n",
        prefix, indent, box.size.x, box.size.y, box.size.z );
    fprintf( stream, "%s%sflat: %p\n",
        prefix, indent, (void*)box.flat );
    fprintf( stream, "%s}\n", prefix );
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
