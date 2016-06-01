////////////////////////////////////////////////////////////////////////////////
// intersect.h - 2015.05.31 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
//
// Requires:
//   point3d.h
//
//
// Functions:
//   intersect1d
//   intersect3d
//
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "point3d.h"


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

int
intersect1d (
    int *min,
    int *max,
    int amin,
    int amax,
    int bmin,
    int bmax
)
{
    if( !min || !max ) return 0;
    if( amax < bmin || bmax < amin ) return 0;

    *min = amin <= bmin ? bmin : amin;
    *max = amax >= bmax ? bmax : amax;

    return 1;
}


int
intersect3d (
    struct POINT3D *min,  // minimum corner of intersection, if any
    struct POINT3D *max,  // maximum corner of intersection, if any
    struct POINT3D amin,
    struct POINT3D amax,
    struct POINT3D bmin,
    struct POINT3D bmax
)
// returns: 0 if non-intersecting: min, max will be unchanged
//          non-0 if intersecting: min, max set to intersection
{
    if( !min || !max ) return 0;

    if (
        intersect1d( &min->x, &max->x, amin.x, amax.x, bmin.x, bmax.x ) &&
        intersect1d( &min->y, &max->y, amin.y, amax.y, bmin.y, bmax.y ) &&
        intersect1d( &min->z, &max->z, amin.z, amax.z, bmin.z, bmax.z )
    ) {
        return 1;
    }

    return 0;
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
