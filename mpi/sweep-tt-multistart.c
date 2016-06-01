////////////////////////////////////////////////////////////////////////////////
// sweep-tt-multistart.c
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// includes
////////////////////////////////////////////////////////////////////////////////

#include "parseargs.h"
#include "forwardstar.h"
#include "sweepxyz.h"


#include "boxfiler.h"
#include "floatbox.h"
#include "intersect.h"
#include "mpihelpers.h"
#include "point3d.h"


#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


////////////////////////////////////////////////////////////////////////////////
// constants
////////////////////////////////////////////////////////////////////////////////

#define	FSRADIUSMAX	7 // maximum radius forward star
#define FSMAX 818     //maximum # of points in a forward star */
//#define MODELMAX	250	/* maximum model dimension in X,Y,Z */
//#define STARTMAX 12 /* maximum starting points */
#define GHOSTDEPTH FSRADIUSMAX

// TODO: remove
//struct FORWARDSTAR {			/* forward start offset */
//  int		i, j, k;	/* point coordinates */
//  float		d;		/* distance to star center (0,0,0)*/
//};


// TODO: remove
//struct MODEL {			/* model point */
//  float		v;		/* velocity */
//  float		tt[STARTMAX];	/* travel time for starting points */
//};


// TODO: remove
//struct START {			/* starting point */
//  int		i, j, k;	/* point coordinates */
//};

//int		changed[STARTMAX];

//struct FORWARDSTAR fs[FSMAX];

//struct FLOATBOX model;
//struct FLOATBOX model_new;
//struct FLOATBOX model_send;
//struct FLOATBOX model_recv;

//struct MODEL model[MODELMAX][MODELMAX][MODELMAX]


//struct MODEL model_new[130][130][51];
//struct MODEL transferS[4][130][130][51];
//struct MODEL transferR[4][130][130][51];

// TODO: remove
//struct START start[STARTMAX], start_new;

//int startinew,startjnew,stopinew,stopjnew;


////////////////////////////////////////////////////////////////////////////////
// box file signatures
////////////////////////////////////////////////////////////////////////////////

const char vbox_sig[4] = {'v', 'b', 'o', 'x'};
const char ttbox_sig[4] = {'t', 't', 'b', 'x'};


////////////////////////////////////////////////////////////////////////////////
// structs
////////////////////////////////////////////////////////////////////////////////

struct STATE {
  struct ARGS args;          // parsed command-line arguments
  int numranks, myrank;      // MPI information
  struct POINT3D rankdims;   // number of ranks along each axis 
  struct POINT3D rankcoords; // my rank coordinates
  struct POINT3D gmin;       // global minimum coordinate of any rank
  struct POINT3D gmax;       // global maximum coordinate of any rank
  struct FLOATBOX vbox;      // contains velocity data and region
  struct {
    int rank;                // if < 0 then no neighbor here, else MPI rank
    struct FLOATBOX send;    // send buffer
    struct FLOATBOX recv;    // receive buffer
  } neighbors[3][3][3];      // for now assume only split along x,y plane
};


////////////////////////////////////////////////////////////////////////////////
// prototypes
// note: prefix do_* because of name conflicts with dynamic libraries
////////////////////////////////////////////////////////////////////////////////

void do_getargs( struct STATE *state, int argc, char *argv[] );
void do_initmpi( struct STATE *state, int argc, char *argv[] );
void do_initstate( struct STATE *state );
void do_loaddatafromfiles( struct STATE *state );
void do_preparebuffers( struct STATE *state );
void do_shutdown( struct STATE *state );


////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

int
main (
  int argc,
  char *argv[]
)
{
  struct STATE state;

  do_initstate( &state );
  do_initmpi( &state, argc, argv );
  do_getargs( &state, argc, argv );
  do_loaddatafromfiles( &state );
  do_preparebuffers( &state );

  do_shutdown( &state );

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
// other functions
////////////////////////////////////////////////////////////////////////////////

void
do_getargs (
  struct STATE *state,
  int argc,
  char *argv[]
)
{
  if( !parseargs( &state->args, argc, argv ) ) {
    if( state->myrank == 0 ) {
      printf(
        "usage: %s <in:velocity.vbox> <in:startpoints.txt> <out:traveltimes.ttbox>\n",
        argv[0]
      );
      fflush( stdout );
    }
    do_shutdown( state );
  }
}


void
do_initmpi (
  struct STATE *state,
  int argc,
  char *argv[]
)
{
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &state->numranks );
  MPI_Comm_rank( MPI_COMM_WORLD, &state->myrank );

  if( !state->myrank ) {
    printf( "MPI ranks: %d\n", state->numranks );
  }
  
  int bestx = splitsquare_numx( state->numranks );
  int besty = state->numranks / bestx;
  int bestz = 1; // TODO: allow other than 1

  state->rankdims = p3d( bestx, besty, bestz );
  state->rankcoords = mpifindrankcoordsfromrank( state->rankdims, state->myrank );
}


void
do_initstate (
  struct STATE *state
)
{
  boxinit( &state->vbox );
}


void
do_loaddatafromfiles (
  struct STATE *state
)
{
  // open file
  struct BOXOPENFILE vboxfile;
  if( !boxfileopenbinary( &vboxfile, state->args.velocityfilename, vbox_sig ) ) {
    fprintf (
      stderr,
      "%d: boxfileopenbinary(.., %s, ..) failed!\n",
      state->myrank, state->args.velocityfilename
    );
    fflush( stdout );
    do_shutdown( state );
  }

  // store file coordinates
  state->gmin = vboxfile.min;
  state->gmax = vboxfile.max;

  // compute inner coordinates for this rank
  struct POINT3D imin, imax;
  mpifindregionfromrankcoords (
    &imin, &imax,
    state->gmin, state->gmax,
    state->rankdims, state->rankcoords
  );

  // compute outer coordinates from this rank
  struct POINT3D omin, omax;
  p3dgrowinside( &omin, &omax, imin, imax, GHOSTDEPTH, state->gmin, state->gmax );

  // load volume of outer coordinates from velocity file
  if( !boxfileloadbinarysubset( &state->vbox, omin, omax, imin, imax, vboxfile ) ) {
    fprintf (
      stderr,
      "%d: boxfileloadbinarysubset(..) from %s failed!\n",
      state->myrank, state->args.velocityfilename
    );
    fflush( stdout );
    do_shutdown( state );
  }
  printf(
    "%d: loaded my region from velocity file: (%d, %d, %d) to (%d, %d, %d)\n",
    state->myrank, imin.x, imin.y, imin.z, imax.x, imax.y, imax.z
  );

  boxfileclosebinary( &vboxfile );
}


void
do_preparebuffers (
  struct STATE *state
)
// allocates send and receive buffers for adjacent neighbors
{
  for( int x = -1; x <= 1; x++ ) {
    for( int y = -1; y <= 1; y++ ) {
      for( int z = -1; z <= 1; z++ ) {
        if( x || y || z ) { // (0,0,0) is me
          struct POINT3D ncoords = p3daddp3d( state->rankcoords, p3d( x, y, z ) );
          int nrank = mpifindrankfromrankcoords( state->rankdims, ncoords );
          if( nrank >= 0 ) {

            // neighbor's inner region (no ghost)
            struct POINT3D nimin, nimax;
            mpifindregionfromrankcoords (
              &nimin, &nimax,
              state->gmin, state->gmax,
              state->rankdims, ncoords
            );

            // neighbor's outer region (including ghost)
            struct POINT3D nomin, nomax;
            p3dgrowinside (
              &nomin, &nomax, nimin, nimax,
              GHOSTDEPTH,
              state->gmin, state->gmax
            );

            // send intersection: neighbor's ghost overlapping my inner region
            struct POINT3D smin, smax;
            // recv intersection: my ghost with neighbor's inner region
            struct POINT3D rmin, rmax;
            // compute intersections
            if (
              !intersect3d (
                &smin, &smax, nomin, nomax, state->vbox.imin, state->vbox.imax
              ) ||
              !intersect3d (
                &rmin, &rmax, nimin, nimax, state->vbox.omin, state->vbox.omax
              )
            ) {
              fprintf (
                stderr,
                "%d: error: no intersection with neighbor %d, but there should be!\n",
                state->myrank, nrank
              );
            }
            else {
              printf (
                "%d: send buffer for neighbor %d: (%d, %d, %d) to (%d, %d, %d)\n",
                state->myrank, nrank,
                smin.x, smin.y, smin.z, smax.x, smax.y, smax.z
              );
              printf (
                "%d: recv buffer for neighbor %d: (%d, %d, %d) to (%d, %d, %d)\n",
                state->myrank, nrank,
                rmin.x, rmin.y, rmin.z, rmax.x, rmax.y, rmax.z
              );
            }
          }
        }
      }
    }
  }
}


void
do_shutdown (
  struct STATE *state
)
{
  boxfree( &state->vbox );
  printf( "%d: shutting down normally\n", state->myrank );
  fflush( stdout );
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize();
  exit(0);
}

////////////////////////////////////////////////////////////////////////////////
// cruft
////////////////////////////////////////////////////////////////////////////////

// TODO: remove
//int
//main_ (
//  int argc,
//  char* argv[]
//)
//{
//  //int		i, j, k, l, m, nx, ny, nz, oi, oj, ok, s;
//  int		numradius, starsize, anychange, numstart, numsweeps=0, numOfTasks=4;
//  int		fsindex[FSRADIUSMAX];
//  float		delta, delay;
//  FILE		*vfile, *fsfile, *ttfile, *startfile;
//  int           numtasks, taskid, len,tag;
//  taskid=0;
//  tag=1;
//  struct POINT3D start[STARTMAX], start_new;
//
//  struct timeval t1, t2,t3, sweep_begin, sweep_end;
//  double elapsedTime;
//  gettimeofday(&t3, NULL);
//
//  const char *velocityfilename = "../docs/velocity-241-241-51.txt";
//  const char *startfilename = "../docs/start-4.txt";
//
//  int sizeOfTasks[4][9]= {
//    {128, 128, 51, 0, 0, 0, 0,7,7},
//    {241, 128, 51, 114, 0, 7, 0, 0,7},
//    {128, 241, 52, 0, 114, 0, 7, 7,0},
//    {241, 241, 51, 114, 114, 7, 7,0,0}
//  };
//
//  int communicationMatrix[16][14]= {
//    {0,0,0,0,0,0,0,0,0,0,0,0,0,0},
//    {0,1,114,120,0,120,0,50,0,6,0,120,0,50},
//    {0,2,0,120,114,120,0,50,0,120,0,6,0,50},
//    {0,3,114,120,114,120,0,50,0,6,0,6,0,50},
//    {1,0,7,13,0,120,0,50,121,127,0,120,0,50},
//    {1,1,0,0,0,0,0,0,0,0,0,0,0,0},
//    {1,2,7,13,114,120,0,50,121,127,0,6,0,50},
//    {1,3,7,126,114,120,0,50,7,126,0,6,0,50},
//    {2,0,0,120,7,13,0,50,0,120,121,127,0,50},
//    {2,1,114,120,7,13,0,50,0,6,121,127,0,50},
//    {2,2,0,0,0,0,0,0,0,0,0,0,0,0},
//    {2,3,114,120,7,126,0,50,0,6,7,126,0,50},
//    {3,0,7,13,7,13,0,50,121,127,121,127,0,50},
//    {3,1,7,126,7,13,0,50,7,126,121,127,0,50},
//    {3,2,7,13,7,126,0,50,121,127,7,126,0,50},
//    {3,3,0,0,0,0,0,0,0,0,0,0,0,0}
//  };
//
//  MPI_Status stat;
//  MPI_Status stats[6];
//  MPI_Request reqs[6];
//  MPI_Init(&argc, &argv);
//  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
//  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
//  MPI_Datatype fstype, modeltype, starttype, oldtypes[2], oldtypes1[1];
//  MPI_Aint    offsets[2], extent, offsets1[1];
//  int          blockcounts[2], blockcounts1[1];
//
//  MPI_Aint intlb, intex;
//  MPI_Aint floatlb, floatex;
//  MPI_Type_get_extent( MPI_INT, &intlb, &intex );
//  MPI_Type_get_extent( MPI_INT, &floatlb, &floatex );
//
//  offsets[0] = 0;
//  oldtypes[0] = MPI_INT;
//  blockcounts[0] = 3;
//  MPI_Type_extent(MPI_INT, &extent);
//  offsets[1] = 3 * extent;
//  oldtypes[1] = MPI_FLOAT;
//  blockcounts[1] = 1;
//  MPI_Type_struct(2, blockcounts, offsets, oldtypes, &fstype);
//  MPI_Type_commit(&fstype);
//
//
//
//  //offset1 = 0;
//  //oldtype1 = MPI_FLOAT;
//  //blockcount1 = STARTMAX+1;
//  offsets1[0] = 0;
//  oldtypes1[0] = MPI_FLOAT;
//  blockcounts1[0] =STARTMAX+1;
//  MPI_Type_struct(1, blockcounts1, offsets1, oldtypes1, &modeltype);
//  MPI_Type_commit(&modeltype);
//
//  //offset1 = 0;
//  //oldtype1 = MPI_INT;
//  //blockcount1 = 3;
//
//  offsets1[0] = 0;
//  oldtypes1[0] = MPI_INT;
//  blockcounts1[0] = 3;
//
//  MPI_Type_struct(1, blockcounts1, offsets1, oldtypes1,  &starttype);
//  MPI_Type_commit(&starttype);
//
//  /* open velocity model file */
//  //printf("Cannot open velocity model file: %s\n", argv[1]);
//  if(taskid==0){
//    printf("\ntask number : %d \n", taskid);
//    //vfile = fopen(argv[1],"r");
//    vfile = fopen( velocityfilename, "r" );
//    if(vfile == NULL) {
//      printf("Cannot open velocity model file: %s\n", argv[1]);
//      exit(1);
//    }
//    printf("Velocity model file: %s\n", argv[1]);
//
//    /* open forward star offset file */
//    //fsfile = fopen(argv[2],"r");
//    fsfile = fopen("../docs/818-FS.txt","r");
//    if(fsfile == NULL) {
//      printf("Cannot open forward star offset file: %s\n", argv[2]);
//      exit(1);
//    }
//    printf("Forward star offset file: %s\n", argv[2]);
//
//    /* open file with starting points */
//    //startfile = fopen(argv[3],"r");
//    startfile = fopen( startfilename, "r" );
//    if(startfile == NULL) {
//      printf("Cannot open starting points file: %s\n", argv[4]);
//      exit(1);
//    }
//    printf("Starting points file: %s\n", argv[3]);
//
//    /* get delta */
//    delta = 10.0;
//    printf("Delta: %f\n", delta);
//
//    /* read velocity model (modified to read original input for Fortran) */
//    nx = 241; ny = 241; nz = 51;
//    printf("Velocity model dimensions: %i %i %i\n",nx, ny, nz);
//    int ir, jr, kr;   /* read indices */
//    for (i=0; i<nx; i++) {
//      for (j=0; j<ny; j++) {
//        for (k=0; k<nz; k++) {
//          fscanf(vfile, "%d,%d,%d,%f", &ir, &jr, &kr, &model[i][j][k].v);
//          if (ir != i+1 || jr != j+1 || kr != k+1) {
//            printf("ERROR: index error reading velocity model "
//                "(%d,%d,%d) (%d,%d,%d)\n", i,j,k, ir,jr,kr);
//            exit(1);
//          }
//        }
//      }
//    }
//    printf("Velocity data read\n");
//    //printf("..... Monil Velocity data read %f", model[241][241][51].v);
//
//    /* read forward star offsets */
//    starsize = 0;
//    fscanf(fsfile, "%i", &starsize);
//    printf("Forward star size: %d\n", starsize);
//
//    for (i=0; i<FSRADIUSMAX; i++) {
//      fsindex[i] = 0;
//    }
//    numradius = 0;
//    for (i=0; i<starsize; i++) {
//      fscanf(fsfile, "%i %i %i", &fs[i].pos.x, &fs[i].pos.y, &fs[i].pos.z);
//      fs[i].distance = sqrt (
//        fs[i].pos.x * fs[i].pos.x +
//        fs[i].pos.y * fs[i].pos.y +
//        fs[i].pos.z * fs[i].pos.z
//      );
//
//      if ((numradius+1) < fs[i].distance) {
//        fsindex[numradius] = i;
//        numradius++;
//      }
//      fs[i].distance = delta * fs[i].distance;
//    }
//    printf("Forward star offsets read\n");
//    //printf("...... Monil Forward star offsets %i %i %i", fs[starsize].i, &fs[starsize].j, &fs[starsize].k);
//    for (i=0; i<FSRADIUSMAX; i++) {
//      printf("numradius: %d, fsindex[%d]: %d\n", numradius, i, fsindex[i]);
//    }
//
//    /* read starting points */
//    fscanf(startfile, "%i", &numstart);
//    for (i=0; i<nx; i++) {
//      for (j=0; j<ny; j++) {
//        for (k=0; k<nz; k++) {
//          for (s=0; s<numstart; s++) {
//            model[i][j][k].tt[s] = INFINITY;
//          }
//        }
//      }
//    }
//    for (s=0; s<numstart; s++) {
//      fscanf(startfile, "%i %i %i", &start[s].x, &start[s].y, &start[s].z);
//      model[start[s].x][start[s].y][start[s].z].tt[s] = 0;
//      printf("starting point %d: %d %d %d\n", s, start[s].x, start[s].y, start[s].z);
//    }
//    printf("Starting points read\n");
//
//
//    ////////////////////////// Send the starting points //////////////////////////////
//
//
//    //start_new=start[0];
//    //if( numtasks>1) {
//    for (i=1; i<numtasks; i++)  {
//      MPI_Send(&start, STARTMAX, starttype, i, tag, MPI_COMM_WORLD);
//      //MPI_Send(&model, 250*250*250, modeltype, i, tag, MPI_COMM_WORLD);
//      //MPI_Send(&fs, FSMAX, fstype, i, tag, MPI_COMM_WORLD);
//    }
//    //}
//
//    //printf("Taskid : %d: starting point %d: %d %d %d\n", taskid, start[0].i, start[0].j, start[0].k);
//
//  } //if ends here
//  else
//  {//STARTMAX
//
//    ////////////////////////// Receive the starting points //////////////////////////////
//    //printf("\ntask number : %d \n", taskid);
//    MPI_Recv(&start, STARTMAX, starttype, 0, tag, MPI_COMM_WORLD, &stat);
//    //MPI_Recv(&model, 250*250*250, modeltype, 0, tag, MPI_COMM_WORLD, &stat);
//    //MPI_Recv(&fs, FSMAX, fstype, 0, tag, MPI_COMM_WORLD, &stat);
//    //printf (
//    //  "Taskid : %d: starting point %d: (%d, %d, %d)\n",
//    //  taskid, start[taskid].x, start[taskid].y, start[taskid].z
//    //);
//
//  }
//
//  ////////////////////////// Broadcast forward start and entire model //////////////////////////////
//
//  MPI_Bcast(fs, FSMAX, fstype, 0, MPI_COMM_WORLD);
//  MPI_Bcast(model, 250*250*250, modeltype, 0, MPI_COMM_WORLD);
//  MPI_Barrier(MPI_COMM_WORLD);
//
//  /*
//   * This section is for data copy to all tasks
//   */
//  int task=0;
//  int starti, startj;
//  numstart=1;
//  for(task=0;task<numOfTasks;task++){
//
//    if (task==taskid){
//      //nx = 127+1; ny = 127+1; nz = 51;
//      nx = sizeOfTasks[taskid][0];
//      ny = sizeOfTasks[taskid][1];
//      nz = sizeOfTasks[taskid][2];
//      //printf("Velocity model dimensions: %i %i %i\n",nx, ny, nz);
//      //int ir, jr, kr;   /* read indices */
//      starti=sizeOfTasks[taskid][3];
//      startj=sizeOfTasks[taskid][4];
//
//      ////////////////////////// Place just the data needed for my node in model_new //////////////////////////////
//
//      printf("for Task : %dVelocity model dimensions: %i %i %i %d %d\n",taskid, nx, ny, nz, starti,startj);
//      for (i=starti; i<nx; i++) {
//        for (j=startj; j<ny; j++) {
//          for (k=0; k<nz; k++) {
//            model_new[i-starti][j-startj][k].v=model[i][j][k].v;
//            for (s=0; s<numstart; s++) {
//              model_new[i-starti][j-startj][k].tt[s]=model[i][j][k].tt[s];
//            }
//          }
//        }
//      }
//
//      ////////////////////////// Set new start offsets //////////////////////////////
//      nx = sizeOfTasks[taskid][0]-sizeOfTasks[taskid][3];
//      ny = sizeOfTasks[taskid][1]-sizeOfTasks[taskid][4];
//      nz = 51;
//      
//      startinew=sizeOfTasks[taskid][5];  // to start from the desired position and ignore the ghost cells
//      startjnew=sizeOfTasks[taskid][6]; // to start from the desired position and ignore the ghost cells
//      stopinew=sizeOfTasks[taskid][7];
//      stopjnew=sizeOfTasks[taskid][8];
//    }
//  }
//
//  //exit(0);
//
//  //memory handling ends here
//
//  //sweep activities begin here
//
//  if (taskid==0)  gettimeofday(&t1, NULL); // time counting
//
//  anychange = 1;
//  s=0;
//  tag=0;
//  int send=0,receive=0,sender=0, reciever=0,reqnumber=0;
//  int new_changeS[4], new_changeR[4];
//  while (anychange) {
//    //for (m=0;m<5;m++){
//    gettimeofday(&sweep_begin, NULL);
//    numsweeps++;
//    anychange = 0;
//    //for (s=0; s<numstart; s++) {
//    for (task=0; task<numOfTasks; task++) {
//      if(taskid==task){
//        changed[s] = 0;
//
//        changed[s] = sweepxyz (
//          FSMAX,
//          fs,
//          model_new,
//          p3d( startinew, startjnew, 0 ),
//          p3d( nx - stopinew, ny - stopjnew, nz - 1 ),
//          s, 0, 818-1, noderegion
//        );
//
//        //printf(">>> start point %d: by Task: %d, changed == %d\n", s,taskid, changed[s]);
//        anychange = changed[s];
//        new_changeS[taskid]=anychange;
//        //if (anychange>0) MPI_Bcast(&anychange, 1, MPI_INT, taskid, MPI_COMM_WORLD);
//        //if (numsweeps<4) anychange=1;
//        gettimeofday(&sweep_end, NULL);
//        elapsedTime = (sweep_end.tv_sec - sweep_begin.tv_sec) + (sweep_end.tv_usec - sweep_begin.tv_usec) / 1000000.0;
//
//        printf("sweep %d finished by task: %d for >>> start point %d: anychange = %d\n",  numsweeps,taskid, s, changed[s]);
//        printf("Time: %f\n", elapsedTime);
//
//        //if (numsweeps<4) anychange=1;
//
//
//
//
//        ////////////////////////// SEND ghost cell data //////////////////////////////
//        ////////
//        ////////       for each send ghost cell region
//        ////////             Get the rank of the neighbor - mpifindneighborrank()
//        ////////             Get the coordinates for region - mpigetsendcoordinates()
//        ////////             Get the data at those coordinates
//        ////////             Send the data to the appropriate rank (non-blocking)
//        ////////
//        //////////////////////////////////////////////////////////////////////////////
//
//        reqnumber=0;
//        for (send=0; send<16; send++) {
//          sender= communicationMatrix[send][0];
//          reciever= communicationMatrix[send][1];
//          if (sender==taskid && sender!=reciever ){
//
//            //printf("sender: %d, reciever: %d, i %d i %d j %d j %d k %d k %d \n",sender,reciever, communicationMatrix[send][2], communicationMatrix[send][3]+1, communicationMatrix[send][4],communicationMatrix[send][5],communicationMatrix[send][6],communicationMatrix[send][7]+1 );
//
//            for (i=communicationMatrix[send][2]; i<communicationMatrix[send][3]+1; i++) {
//              for (j=communicationMatrix[send][4]; j<communicationMatrix[send][5]+1; j++) {
//                for (k=communicationMatrix[send][6]; k<communicationMatrix[send][7]+1; k++) {
//                  //for (s=0; s<numstart; s++) {
//
//                  ////////////////////////// transferS is send buffer //////////////////////////////
//
//                  transferS[reciever][i-communicationMatrix[send][2]][j-communicationMatrix[send][4]][k].tt[s]=model_new[i][j][k].tt[s];
//                  //printf("sender: %d, reciever: %d, i %d j %d k %d \n",sender,reciever, i, j, k);
//                  //}
//                }
//              }
//            }
//            MPI_Isend(&transferS[reciever], 130*130*51, modeltype, reciever, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
//            reqnumber++;
//
//
//
//          }//for comparing sender and reciever
//          // for task
//
//        } //sender for ends here
//
//
//
//        ////////////////////////// RECEIVE ghost cell data, transferR is receive buffer //////////////////////////////
//        ////////
//        ////////       for each receive ghost cell region
//        ////////             Get the rank of the neighbor - mpifindneighborrank()
//        ////////             Get the coordinates for region - mpigetreceivecoordinates()
//        ////////             Get the data at those coordinates
//        ////////             Receive the data from the appropriate rank
//        ////////                     Non-blocking - separate buffer needed for each region, wait for all to complete
//        ////////                     Blocking - single buffer needed, but has to wait in order
//        ////////
//        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//        for (receive=0; receive<16; receive++) {
//          sender= communicationMatrix[receive][0];
//          reciever= communicationMatrix[receive][1];
//          if (reciever==taskid && sender!=reciever){
//            MPI_Irecv(&transferR[sender], 130*130*51, modeltype, sender, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
//            reqnumber++;
//
//          }//for comparing sender and reciever
//
//
//        } //reciever for ends here
//
//
//        //now block until requests are complete
//        MPI_Waitall(reqnumber, reqs, stats);
//
//        reqnumber=0;
//
//        ////////////////////////// Send + receive anychange //////////////////////////////
//
//        for(i=0;i<4;i++){
//          if (i!=taskid) {
//            //MPI_Send(&anychange, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
//            MPI_Isend(&new_changeS[taskid], 1, MPI_INT, i, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
//            reqnumber++;
//          }
//        }
//        for(i=0;i<4;i++){
//          if (i!=taskid){
//            // MPI_Recv(&new_changeR[taskid], 1, MPI_INT, i, tag, MPI_COMM_WORLD);
//            MPI_Irecv(&new_changeR[i], 1, MPI_INT, i, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
//            reqnumber++;
//          }
//        }
//        //MPI_Recv(&start, STARTMAX, starttype, 0, tag, MPI_COMM_WORLD, &stat);
//
//        MPI_Waitall(reqnumber, reqs, stats);
//
//        for(i=0;i<4;i++){
//          if (i!=taskid && new_changeR[i]>0){ anychange=new_changeR[i];
//          }
//        }
//
//
//        int ii=0;
//
//
//        ////////////////////////// Go through receive buffers and place ghost data in model_new //////////////////////////////
//
//        for (receive=0; receive<16; receive++) {
//          sender= communicationMatrix[receive][0];
//          reciever= communicationMatrix[receive][1];
//          if (reciever==taskid && sender!=reciever){
//            //printf("sender: %d, reciever: %d, i %d i %d j %d j %d k %d k %d \n",sender,reciever, communicationMatrix[receive][8], communicationMatrix[receive][9]+1, communicationMatrix[receive][10],communicationMatrix[receive][11]+1,communicationMatrix[receive][12],communicationMatrix[receive][13]+1 );
//            ii=0;
//
//            for (i=communicationMatrix[receive][8]; i<communicationMatrix[receive][9]+1; i++) {
//              for (j=communicationMatrix[receive][10]; j<communicationMatrix[receive][11]+1; j++) {
//                for (k=communicationMatrix[receive][12]; k<communicationMatrix[receive][13]+1; k++) {
//
//                  //printf (" task : %d, vaule of tt %f \n", taskid,transferR[sender][i-communicationMatrix[receive][8]][j-communicationMatrix[receive][10]][k].tt[s]);
//                  //ii++;
//
//                  //printf (" task : %d, vaule of i:%d j:%d k:%d \n", taskid,i,j,k);
//                  model_new[i][j][k].tt[s]=transferR[sender][i-communicationMatrix[receive][8]][j-communicationMatrix[receive][10]][k].tt[s];
//                  //ii++;
//
//                }
//              }
//            }
//            //printf("task: %d sender:%d, reciever %d, ii: %d\n",taskid, sender, reciever,ii);
//
//          }//for comparing sender and reciever
//          //printf ( "\n value of i : %d",i);
//        } //reciever for ends here
//
//
//
//      } //task comparing if ends here
//    } //starting for ends here for ends here
//
//
//    //MPI_Barrier(MPI_COMM_WORLD);
//    //exit(0);
//    //anychange = 1;
//    //printf("Taskid: %d anychange: %d \n",taskid,anychange);
//  } //while ends here
//
//
//  //exit(0);
//
//  /**
//   * TODO: Rank 0 should gather results from other processes or it's OK to have
//   *       each process write its own output.  Make sure that the filenames
//   *       don't collide.  At least think about how rank 0 could gather output
//   *       results.
//   */
//
//  /* TODO: Remove exit statement so output can complete. */
//  MPI_Barrier(MPI_COMM_WORLD);
//
//  //time measurement for sweep ends here
//
//  if (taskid==0) {
//    gettimeofday(&t2, NULL);
//
//    // compute and print the elapsed time in millisec
//    elapsedTime = (t2.tv_sec - t1.tv_sec);      // sec to ms
//    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
//    printf( "\n elapsed Second for sweeps only: %f \n", elapsedTime);
//  }
//
//
//  if (taskid==0) {
//    /* print travel times */
//
//    ttfile = fopen("mpi/shortest-path/output.tt","w");
//    if(ttfile == NULL) {
//      printf(".......................Can not open travel time output file: %s\n", "output.tt");
//      exit(1);
//    }
//    fprintf(ttfile, "%d %d %d\n", nx, ny, nz);
//    for (s=0; s<numstart; s++) {
//      fprintf(ttfile, "starting point: %d\n", s);
//      for (i=0; i<nx; i++) {
//        for (j=0; j<ny; j++) {
//          for (k=0; k<nz; k++) {
//            /* use %g for doubles */
//            fprintf(ttfile, "travel time for (%d,%d,%d): %f %d %d %d\n",
//                i, j, k, model[i][j][k].tt[s], 0, 0, 0);
//          }
//        }
//      }
//    }
//  }
//
//  MPI_Barrier(MPI_COMM_WORLD);
//  MPI_Type_free(&modeltype);
//  MPI_Type_free(&fstype);
//  MPI_Type_free(&starttype);
//  MPI_Finalize();
//
//
//  gettimeofday(&t2, NULL);
//
//  // compute and print the elapsed time in millisec
//  elapsedTime = (t2.tv_sec - t3.tv_sec);      // sec to ms
//  elapsedTime += (t2.tv_usec - t3.tv_usec) / 1000000.0;   // us to ms
//  printf( "\n elapsed Second for whole program: %f \n", elapsedTime);
//}


////////////////////////////////////////////////////////////////////////////////
// vim: set tabstop=2 shiftwidth=2 softtabstop=2 expandtab:
// END
////////////////////////////////////////////////////////////////////////////////
