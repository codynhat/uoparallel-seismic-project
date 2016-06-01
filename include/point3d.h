////////////////////////////////////////////////////////////////////////////////
// point3d.h - 2016.05.23 - Atlee Brink
////////////////////////////////////////////////////////////////////////////////
//
// Data:
//   struct POINT3D
//
//
// Functions:
//   p3d
//
//   p3daddp3d
//   p3daddval
//   p3dsubbp3d
//
//   p3disless
//   p3dismore
//   p3disnotequal
//
//   p3dcalcvolume
//   p3dsizeofregion
//
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// structs
////////////////////////////////////////////////////////////////////////////////

struct POINT3D {
    int x, y, z;
};


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

inline extern
struct POINT3D
p3d (
    const int x, const int y, const int z
)
// for convenience to prevent (struct POINT3D) everywhere
{
    return (struct POINT3D){x, y, z};
}


inline extern
struct POINT3D
p3daddp3d (
    const struct POINT3D pta,
    const struct POINT3D ptb
)
// pta + ptb
{
    return (struct POINT3D){pta.x + ptb.x, pta.y + ptb.y, pta.z + ptb.z};
}


inline extern
struct POINT3D
p3daddval (
    const struct POINT3D pt,
    const int val
)
// add integer to all components of point
// integer may be negative
{
    return p3daddp3d( pt, p3d(val, val, val) );
}


inline extern
int
p3disless (
    const struct POINT3D a,
    const struct POINT3D b
)
// weak comparison: a < b for any component
{
    return a.x < b.x || a.y < b.y || a.z < b.z;
}


inline extern
int
p3dismore (
    const struct POINT3D a,
    const struct POINT3D b
)
// weak comparison: a > b for any component
{
    return a.x > b.x || a.y > b.y || a.z > b.z;
}


inline extern
int
p3disnotequal (
    const struct POINT3D a,
    const struct POINT3D b
)
{
    return a.z != b.z || a.y != b.y || a.x != b.x;
}


inline extern
struct POINT3D
p3dsubp3d (
    const struct POINT3D pta,
    const struct POINT3D ptb
)
// pta - ptb
{
    return (struct POINT3D){pta.x - ptb.x, pta.y - ptb.y, pta.z - ptb.z};
}


inline extern
struct POINT3D
p3dsizeofregion (
    const struct POINT3D min,
    const struct POINT3D max
)
{
    return p3daddval( p3dsubp3d( max, min ), 1 );
}


inline extern
long
p3dcalcvolume (
    const struct POINT3D size
)
// correctly computes volume even if larger than int
{
    return (long)size.x * size.y * size.z;
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
