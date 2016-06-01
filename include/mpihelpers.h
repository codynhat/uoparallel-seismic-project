////////////////////////////////////////////////////////////////////////////////
// mpihelpers.h - 2016.05.25 - Cody Hatfield
////////////////////////////////////////////////////////////////////////////////
//
// Functions:
//   mpicalculatemycoordinates
//   mpifindneighborrank
//   mpigetsendcoordinates
//   mpigetreceivecoordinates
//
//
// Ghost cell ids:
//   Starting from SW, counter-clockwise
//   6 5 4
//   7   3
//   0 1 2
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "point3d.h"
#include "splitsquare.h"


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

void
mpicalculatemycoordinates (
    struct POINT3D *min,
    struct POINT3D *dims,
    int my_rank,
    int total_nodes,
    struct POINT3D src_min,
    struct POINT3D src_dims
)
// returns the global coordinates for a particular rank
// assumes z is depth
// assumes total_nodes is power of 2
{
  if( min == 0 || dims == 0 ) return;

  int num_x = splitsquare_numx(total_nodes);
  int num_y = (total_nodes/num_x);

  int dims_x = (src_dims.x/num_x);
  int dims_y = (src_dims.y/num_y);

  int col_rank = my_rank % num_x;
  int row_rank = my_rank / num_x;

  min->x = (dims_x * col_rank) + src_min.x;
  dims->x = col_rank < (num_x - 1) ? dims_x : src_dims.x - (dims_x * (num_x - 1));
  min->y = (dims_y * row_rank) + src_min.y;
  dims->y = row_rank < (num_y - 1) ? dims_y : src_dims.y - (dims_y * (num_y - 1));
  min->z = src_min.z;
  dims->z = src_dims.z;
}


void
mpifindregionfromrankcoords (
    struct POINT3D *min,
    struct POINT3D *max,
    const struct POINT3D gmin,
    const struct POINT3D gmax,
    const struct POINT3D rankdims,
    const struct POINT3D rankcoords
)
// to ensure no floating-point weirdness, block size is integral;
// excess goes to most positive rank coords
{
    struct POINT3D size = p3dsizeofregion( gmin, gmax );
    struct POINT3D block = { size.x / rankdims.x, size.y / rankdims.y, size.z / rankdims.z };

    *min = p3daddp3d (
        gmin,
        p3d( rankcoords.x * block.x, rankcoords.y * block.y, rankcoords.z * block.z )
    );
    *max = p3daddp3d( *min, p3daddval( block, -1 ) );

    if( rankcoords.x + 1 == rankdims.x ) max->x = gmax.x;
    if( rankcoords.y + 1 == rankdims.y ) max->y = gmax.y;
    if( rankcoords.z + 1 == rankdims.z ) max->z = gmax.z;
}


int
mpifindrankfromrankcoords (
    const struct POINT3D rankdims,   // number of ranks along each dimension
    const struct POINT3D rankcoords  // coordinates within rankdims
)
// returns -1 if no rank at rankcoords, else a valid rank number
{
    if (
        p3disless( rankcoords, p3d( 0, 0, 0 ) ) ||
        p3dismore( rankcoords, p3daddval( rankdims, -1 ) )
    ) {
        return -1;
    }

    int spanz = 1;
    int spany = rankdims.z * spanz;
    int spanx = rankdims.y * spany; 

    return rankcoords.x * spanx + rankcoords.y * spany + rankcoords.z * spanz;
}


struct POINT3D
mpifindrankcoordsfromrank (
    const struct POINT3D rankdims,
    int rank
)
{
    int spanz = 1;
    int spany = rankdims.z * spanz;
    int spanx = rankdims.y * spany; 

    int x = rank / spanx;
    rank -= x * spanx;

    int y = rank / spany;
    rank -= y * spany;

    return p3d( x, y, rank );
}


int
mpifindneighborrank_ (
    int my_rank,
    int ghost_id,   // see ^Ghost cell ids
    int width
)
// returns mpi rank for a ghost cell id
// < 0 is no neighbor
{
  int north     = (my_rank+width < width*width) ? my_rank+width : -1;
  int south     = (my_rank-width);
  int west      = (my_rank % width == 0) ? -1 : my_rank-1;
  int east      = (my_rank % width == width - 1) ? -1 : my_rank+1;

  switch (ghost_id) {
    case 0: // SW
      return (south >= 0 && west >= 0) ? (south-1) : -1;
      break;
    case 1: // S
      return south;
      break;
    case 2: // SE
      return (south >= 0 && east >= 0) ? (south+1) : -1;
      break;
    case 3: // E
      return east;
      break;
    case 4: // NE
      return (north >= 0 && east >= 0) ? (north+1) : -1;
      break;
    case 5: // N
      return north;
      break;
    case 6: // NW
      return (north >= 0 && west >= 0) ? (north-1) : -1;
      break;
    case 7: // W
      return west;
      break;
    default:
      return -1;
      break;
  }
}


void
mpigetsendcoordinates (
    struct POINT3D *min,
    struct POINT3D *dims,
    int ghost_id,       // see ^Ghost cell ids
    int ghost_size,      // probably want 7
    struct POINT3D imin,
    struct POINT3D imax
    //int min_x,
    //int min_y,
    //int max_x,
    //int max_y,
    //int min_z,
    //int max_z
    //int z
)
// returns send coordinates for a given ghost cell and dims
{
  if( min == 0 || dims == 0 ) return;

  switch (ghost_id) {
    case 0:
      min->x = imin.x - ghost_size;
      min->y = imin.y - ghost_size;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 1:
      min->x = imin.x;
      min->y = imin.y - ghost_size;
      dims->x = imax.x - imin.x;
      dims->y = ghost_size;
      break;
    case 2:
      min->x = imax.x;
      min->y = imin.y - ghost_size;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 3:
      min->x = imax.x;
      min->y = imin.y;
      dims->x = ghost_size;
      dims->y = imax.y - imin.y;
      break;
    case 4:
      min->x = imax.x;
      min->y = imax.y;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 5:
      min->x = imin.x;
      min->y = imax.y;
      dims->x = imax.x - imin.x;
      dims->y = ghost_size;
      break;
    case 6:
      min->x = imin.x - ghost_size;
      min->y = imax.y;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 7:
      min->x = imin.x - ghost_size;
      min->y = imin.y;
      dims->x = ghost_size;
      dims->y = imax.y - imin.y;
      break;
  }
  min->z = imin.z;
  dims->z = imax.z - imin.z + 1;
}


void
mpigetreceivecoordinates (
    struct POINT3D *min,     // out
    struct POINT3D *dims,    // out
    int ghost_id,            // see ^Ghost cell ids
    int ghost_size,          // probably want 7
    int min_x,
    int min_y,
    int max_x,
    int max_y,
    int z
)
// returns receive coordinates for a given ghost cell and dims
{
  if( min == 0 || dims == 0 ) return;

  switch (ghost_id) {
    case 0:
      min->x = min_x;
      min->y = min_y;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 1:
      min->x = min_x;
      min->y = min_y;
      dims->x = max_x-min_x;
      dims->y = ghost_size;
      break;
    case 2:
      min->x = max_x-ghost_size;
      min->y = min_y;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 3:
      min->x = max_x-ghost_size;
      min->y = min_y;
      dims->x = ghost_size;
      dims->y = max_y-min_y;
      break;
    case 4:
      min->x = max_x-ghost_size;
      min->y = max_y-ghost_size;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 5:
      min->x = min_x;
      min->y = max_y-ghost_size;
      dims->x = max_x-min_x;
      dims->y = ghost_size;
      break;
    case 6:
      min->x = min_x;
      min->y = max_y-ghost_size;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 7:
      min->x = min_x;
      min->y = min_y;
      dims->x = ghost_size;
      dims->y = max_y-min_y;
      break;
  }
  min->z = 0;
  dims->z = z;
}


////////////////////////////////////////////////////////////////////////////////
// vim: set tabstop=4 shiftwidth=4 softtabstop=4 expandtab:
// END
////////////////////////////////////////////////////////////////////////////////
