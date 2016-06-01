////////////////////////////////////////////////////////////////////////////////
// sweepxyz.h - split from sweep-tt-multistart.c 2016.05.31
////////////////////////////////////////////////////////////////////////////////
//
// Requires:
//   forwardstar.h
//
//   floatbox.h
//   point3d.h
//
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "forwardstar.h"


#include "floatbox.h"
#include "point3d.h"


#include <math.h>
#include <stdio.h>


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

long
sweepxyz (
  const int FSMAX,
  const struct FORWARDSTAR fs[FSMAX],
  const struct FLOATBOX vbox,  // velocity box: doesn't change
  const struct FLOATBOX ttbox, // travel times: ttbox.flat[..] may change
  const struct POINT3D startpoint,
  const int starstart, const int starstop
)
{
  long change = 0;

  #pragma omp parallel for private(here, there, l, vel_here, vel_there, tt_here, tt_there, delay)\
    default(shared) reduction(+:change) schedule(dynamic) num_threads(16)
  //#pragma omp parallel for private(oi, oj, ok, x, y, z, l, tt, tto, delay) \
    default(shared) reduction(+:change) schedule(dynamic) num_threads(16)
  for( struct POINT3D here = vbox.imin; here.x <= vbox.imax.x; here.x++ ) {
    for( here.y = vbox.imin.y; here.y <= vbox.imax.y; here.y++ ) {
      for( here.z = vbox.imin.z; here.z <= vbox.imax.z; here.z++ ) {
        for( int l = starstart; l < starstop; l++ ) {

          // find point in forward star based on offsets
          struct POINT3D there = p3daddp3d( here, fs[l].pos );

          // if 'there' is outside the boundaries, then skip
          if (
            p3disless( there, vbox.omin ) ||
            p3dismore( there, vbox.omax )
          ) {
            continue;
          }
          
          // compute delay from 'here' to 'there' with endpoint average
          float vel_here = boxgetglobal( vbox, here );
          float vel_there = boxgetglobal( vbox, there );
          float delay = fs[l].distance * (vel_here + vel_there) * 0.5f;
          
          // ignore the starting point
          if( p3disnotequal( here, startpoint ) ) {

            float tt_here = boxgetglobal( ttbox, here );
            float tt_there = boxgetglobal( ttbox, there );

            // if offset point has infinity travel time, then update
            if ((tt_here == INFINITY) && (tt_there == INFINITY)) {
              continue;
            }

            if ((tt_here != INFINITY) && (tt_there == INFINITY)) {
              boxputglobal( ttbox, there, delay + tt_here );
              change++;
              continue;
            }

            if ((tt_here == INFINITY) && (tt_there != INFINITY)) {
              boxputglobal( ttbox, here, delay + tt_there );
              change++;
              continue;
            }

            if ((tt_here != INFINITY) && (tt_there != INFINITY)) {
              // if a shorter travel time through 'there', update 'here'
              if ((delay + tt_there) < tt_here) {
                boxputglobal( ttbox, here, delay + tt_there );
                change++;
              }
              // if a shorter travel time through 'here', update 'there'
              else if ((delay + tt_here) < tt_there) {
                boxputglobal( ttbox, there, delay + tt_here );
                change++;
              }
            }
          }
        }
      }
    }
  }

  return change;
}

////////////////////////////////////////////////////////////////////////////////
// vim: set tabstop=2 shiftwidth=2 softtabstop=2 expandtab:
// END
////////////////////////////////////////////////////////////////////////////////
