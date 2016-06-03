////////////////////////////////////////////////////////////////////////////////
// common.h
////////////////////////////////////////////////////////////////////////////////

#pragma once


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "boxfiler.h"
#include "floatbox.h"
#include "intersect.h"
#include "mpihelpers.h"
#include "point3d.h"
//#include "timing.h"


#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


////////////////////////////////////////////////////////////////////////////////
// structs
////////////////////////////////////////////////////////////////////////////////

struct FORWARDSTAR {
  struct POINT3D pos;
  float halfdistance;
};


struct NEIGHBOR {
  struct POINT3D relation; // (-1,-1,-1) to (1,1,1)
  int rank;                // MPI rank
  struct FLOATBOX send;    // send buffer
  struct FLOATBOX recv;    // receive buffer
};


struct STATE {
  struct ARGS args;              // parsed command-line arguments

  int numranks, myrank;          // MPI information
  struct POINT3D rankdims;       // number of ranks along each axis 
  struct POINT3D rankcoords;     // my rank coordinates

  struct POINT3D gmin;           // global minimum coordinate of any rank
  struct POINT3D gmax;           // global maximum coordinate of any rank
  struct FLOATBOX vbox;          // contains velocity data and region

  struct FLOATBOX ttbox;         // contains travel time data and region
  struct POINT3D ttstart;        // for now just one ttbox and one ttstart
  long numsweeps;                // counts how many sweeps this rank has done

  struct NEIGHBOR neighbors[26]; // could be 0 to 26 actual neighbors
  int numneighbors;              // 0 to 26

  struct FORWARDSTAR *star;      // remember to free later
  int numinstar;
};


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
