////////////////////////////////////////////////////////////////////////////////
// sweep.h
////////////////////////////////////////////////////////////////////////////////
//
// Functions:
//   do_sweep( struct STATE *state )
//
//
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

long
do_sweep (
  struct STATE *state  
)
{
  // copy some state into local stack memory for fastness
  const struct FLOATBOX vbox = state->vbox;
  const struct FLOATBOX ttbox = state->ttbox;
  const struct POINT3D ttstart = state->ttstart;
  const struct FORWARDSTAR * const star = state->star;
  const int numinstar = state->numinstar;

  // count how many (if any) values we change
  long changes = 0;

  #pragma omp parallel for default(shared) \
    reduction(+:changes) schedule(dynamic) num_threads(16)
  for( int x = vbox.imin.x; x <= vbox.imax.x; x++ ) {
    printf( "%d: x = %d\n", state->myrank, x );
    for( int y = vbox.imin.y; y <= vbox.imax.y; y++ ) {
      for( int z = vbox.imin.z; z <= vbox.imax.z; z++ ) {

        const struct POINT3D here = p3d( x, y, z );

        const float vel_here = boxgetglobal( vbox, here );
        const float tt_here = boxgetglobal( ttbox, here );

        for( int l = 0; l < numinstar; l++ ) {

          // find point in forward star based on offsets
          const struct POINT3D there = p3daddp3d( here, star[l].pos );

          // if 'there' is outside the boundaries, then skip
          if (
            p3disless( there, vbox.omin ) ||
            p3dismore( there, vbox.omax )
          ) {
            continue;
          }
          
          // compute delay from 'here' to 'there' with endpoint average
          const float vel_there = boxgetglobal( vbox, there );
          const float delay = star[l].halfdistance * (vel_here + vel_there);
          
          // ignore the starting point
          if( p3disnotequal( here, ttstart ) ) {

            const float tt_there = boxgetglobal( ttbox, there );

            // if offset point has infinity travel time, then update
            if ((tt_here == INFINITY) && (tt_there == INFINITY)) {
              continue;
            }

            if ((tt_here != INFINITY) && (tt_there == INFINITY)) {
              boxputglobal( ttbox, there, delay + tt_here );
              changes++;
              continue;
            }

            if ((tt_here == INFINITY) && (tt_there != INFINITY)) {
              boxputglobal( ttbox, here, delay + tt_there );
              changes++;
              continue;
            }

            if ((tt_here != INFINITY) && (tt_there != INFINITY)) {
              // if a shorter travel time through 'there', update 'here'
              if ((delay + tt_there) < tt_here) {
                boxputglobal( ttbox, here, delay + tt_there );
                changes++;
              }
              // if a shorter travel time through 'here', update 'there'
              else if ((delay + tt_here) < tt_there) {
                boxputglobal( ttbox, there, delay + tt_here );
                changes++;
              }
            }
          }
        }
      }
    }
  }

  return changes;
}


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
