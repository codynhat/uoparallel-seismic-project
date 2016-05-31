////////////////////////////////////////////////////////////////////////////////
// sweepxyz.h - split from sweep-tt-multistart.c 2016.05.31
////////////////////////////////////////////////////////////////////////////////
//
// Requires:
//   forwardstar.h
//   point3d.h
//
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "forwardstar.h"


#include "point3d.h"


#include <stdio.h>


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

int
sweepxyz (
  const int FSMAX,
  const struct FORWARDSTAR fs[FSMAX],
  const struct MODEL model_new[130][130][51],
  const struct POINT3D min, // x,y,z of minimum corner (inclusive)
  const struct POINT3D max, // x,y,z of maximum corner (inclusive)
  int s, int starstart, int starstop
)
{
  int	oi, oj, ok;
  int	change = 0;
  float	delay = 0.0, tt = 0.0, tto = 0.0;

  printf (
    "sweepXYZ (\n"
    "  FSMAX: %d\n"
    "  fs[FSMAX]: (stuff)\n"
    "  min: (%d, %d, %d)\n"
    "  max: (%d, %d, %d)\n"
    "  s: %d, starstart: %d, starstop: %d\n"
    ")\n",
    FSMAX,
    min.x, min.y, min.z,
    max.x, max.y, max.z,
    s, starstart, starstop
  );

  #pragma omp parallel for private(pt, x, y, z, l, tt, tto, delay) \
    default(shared) reduction(+:change) schedule(dynamic) num_threads(16)
  //#pragma omp parallel for private(oi, oj, ok, x, y, z, l, tt, tto, delay) \
    default(shared) reduction(+:change) schedule(dynamic) num_threads(16)
  for( int x = min.x; x <= max.x; x++ ) {
    for( int y = min.y; y <= max.y; y++ ) {
      for( int z = min.z; z <= max.z; z++ ) {
  /* TODO: remove
  for (int x=startinew; x<size.x-stopinew; x++) {
    for (int y=startjnew; y<size.y-stopjnew; y++) {
      for (int z=0; z<size.z; z++) {
  */
        for (int l=starstart; l<starstop; l++) {
          /* find point in forward star based on offsets */

          // TODO: remove
          //oi = x + fs[l].pos.x;
          //oj = y + fs[l].pos.y;
          //ok = z + fs[l].pos.z;
          
          // TODO: remove
          //struct POINT3D pt = {oi, oj, ok};

          struct POINT3D pt = p3daddp3d( p3d( x, y, z ), fs[l].pos );


          /* if (oi,oj,ok) is outside the boundaries, then skip */
          if( p3disless( pt, min ) || p3dismore( pt, max ) ) continue;
          /* TODO: remove
          if ((oi < 0) || (oi > size.x-1)
              || (oj < 0) || (oj > size.y-1)
              || (ok < 0) || (ok > size.z-1)) {

            continue;
          }
          */
          /* compute delay from (x,y,z) to (pt.x,pt.y,pt.z) with end point average */
          delay = fs[l].distance * (model_new[x][y][z].v + model_new[pt.x][pt.y][pt.z].v) / 2.0;
          /* update travel times for all starting points */
          /* if (x,y,z) is starting point, then skip */

          if ((x == start[s].x) && (y == start[s].y) && (z == start[s].z)) {

            continue;
          }
          tt = model_new[x][y][z].tt[s];
          tto = model_new[pt.x][pt.y][pt.z].tt[s];
          /* if offset point has infinity travel time, then update */
          if ((tt == INFINITY) && (tto == INFINITY)) {

            continue;
          }
          if ((tt != INFINITY) && (tto == INFINITY)) {
            model_new[pt.x][pt.y][pt.z].tt[s] = delay + tt;
            change += 1;
            continue;
          }
          if ((tt == INFINITY) && (tto != INFINITY)) {
            model_new[x][y][z].tt[s] = delay + tto;
            change += 1;
            continue;
          }
          if ((tt != INFINITY) && (tto != INFINITY)) {
            /* if a shorter travel time through (pt.x,pt.y,pt.z), update (x,y,z) */
            if ((delay + tto) < tt) {
              model_new[x][y][z].tt[s] = delay + tto;
              change += 1;
            }
            /* if a shorter travel time through (x,y,z), update (pt.x,pt.y,pt.z) */
            else if ((delay + tt) < tto) {
              model_new[pt.x][pt.y][pt.z].tt[s] = delay + tt;
              change += 1;
            }
          }
        }
      }
    }
  }

  return(change);

}

////////////////////////////////////////////////////////////////////////////////
// vim: set tabstop=2 shiftwidth=2 softtabstop=2 expandtab:
// END
////////////////////////////////////////////////////////////////////////////////
