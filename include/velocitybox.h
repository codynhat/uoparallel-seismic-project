////////////////////////////////////////////////////////////////////////////////
// velocitybox.h - 2016.05.26 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
// 
// A container for a regular-grid velocity volume using global coordinates.
//
//
// Requires:
//   floatbox.h
//   point3d.h
//
//
// Data:
//   struct VELOCITYBOX
//
//
// Functions:
//   vboxalloc
//   vboxfree
//
//   vboxget
//   vboxput
//
//   vboxfprint
//
// 
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "floatbox.h"
#include "point3d.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////
// structs 
////////////////////////////////////////////////////////////////////////////////

struct VELOCITYBOX {
    struct POINT3D min, max;
    struct FLOATBOX box; // contains dimensions and velocity data
};


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

int
vboxalloc (
    struct VELOCITYBOX *vbox,
    const struct POINT3D min,
    const struct POINT3D max
)
// allocates and initializes a new VELOCITYBOX
// on error: returns 0 (failure to allocate memory)
// on success: returns non-zero
// note: you should use vboxfree(&vbox) when you're done with it.
{
    if( !boxalloc( &vbox->box, p3dsizeofregion( min, max ) ) ) return 0;

    vbox->min = min;
    vbox->max = max;

    return 1;
}


void
vboxfree (
    struct VELOCITYBOX *vbox
)
// releases heap memory associated with vbox
{
    if( vbox == NULL ) return;
    boxfree( &vbox->box );
}


inline extern
float
vboxget (
    const struct VELOCITYBOX vbox,
    const struct POINT3D pt
)
// returns a single value from the given global coordinates
{
    return boxget( vbox.box, p3dsubp3d( pt, vbox.min ) );
}


inline extern
void
vboxput (
    const struct VELOCITYBOX vbox,
    const struct POINT3D pt,
    const float val
)
// stores a single value at the given global coordinates
{
    boxput( vbox.box, p3dsubp3d( pt, vbox.min ), val );
}


void
vboxfprint (
    FILE *stream,
    const char *prefix,
    const char *indent,
    const struct VELOCITYBOX vbox
)
// metadata friendly-printing
{
    if( stream == NULL ) stream = stdout;
    if( prefix == NULL ) prefix = "";
    if( indent == NULL ) indent = "  ";

    fprintf( stream, "%sVELOCITYBOX {\n", prefix );
    fprintf( stream, "%s%sminimum corner: (%d, %d, %d)\n",
        prefix, indent, vbox.min.x, vbox.min.y, vbox.min.z );
    fprintf( stream, "%s%smaximum corner: (%d, %d, %d)\n",
        prefix, indent, vbox.max.x, vbox.max.y, vbox.max.z );

    char *newprefix = malloc( strlen( prefix ) + strlen( indent ) + 1 );
    if( newprefix != NULL ) {
        strcpy( newprefix, prefix );
        strcpy( newprefix + strlen( prefix ), indent );
        boxfprint( stream, newprefix, indent, vbox.box );
        free( newprefix );
    } else {
        boxfprint( stream, prefix, indent, vbox.box );
    }

    fprintf( stream, "%s}\n", prefix );
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
