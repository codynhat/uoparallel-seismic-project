////////////////////////////////////////////////////////////////////////////////
// sweep.h
////////////////////////////////////////////////////////////////////////////////
//
// Functions:
//   do_sweep( struct STATE *state )
//
//
////////////////////////////////////////////////////////////////////////////////

#include "sweep-cuda.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "parseargs.h"
#include "common.h"

////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

__device__ long changes;


__global__
void
do_sweep_cuda(
  const struct FLOATBOX *vbox,
  const struct FLOATBOX *ttbox,
  const struct POINT3D *ttstart,
  const struct FORWARDSTAR * const star,
  const int *numinstar
)
{ 
  int x = blockIdx.x * blockDim.x + threadIdx.x + vbox->imin.x;
  int y = blockIdx.y * blockDim.y + threadIdx.y + vbox->imin.y;
  int z = blockIdx.z * blockDim.z + threadIdx.z + vbox->imin.z;
  
  const struct POINT3D here = p3d( x, y, z );
  
  const float vel_here = boxgetglobal( *vbox, here );
  const float tt_here = boxgetglobal( *ttbox, here );
  
   
  for( int l = 0; l < *numinstar; l++ ) {
    // find point in forward star based on offsets
    const struct POINT3D there = p3daddp3d( here, star[l].pos );
  
    // if 'there' is outside the boundaries, then skip
    if (
      p3disless( there, vbox->omin ) ||
      p3dismore( there, vbox->omax )
    ) {
      continue;
    }
    
    // compute delay from 'here' to 'there' with endpoint average
    const float vel_there = boxgetglobal( *vbox, there );
    const float delay = star[l].halfdistance * (vel_here + vel_there);
  
    // ignore the starting point
    if( p3disnotequal( here, *ttstart ) ) {
      const float tt_there = boxgetglobal( *ttbox, there );
      //printf("%f, %f\n", tt_here, tt_there); 
      
      // if offset point has infinity travel time, then update
      if ((tt_here == INFINITY) && (tt_there == INFINITY)) {
        continue;
      }
  
      if ((tt_here != INFINITY) && (tt_there == INFINITY)) {
        boxputglobal( *ttbox, there, delay + tt_here );
        changes++;
        continue;
      }
  
      if ((tt_here == INFINITY) && (tt_there != INFINITY)) {
        boxputglobal( *ttbox, here, delay + tt_there );
        changes++;
        continue;
      }
  
      if ((tt_here != INFINITY) && (tt_there != INFINITY)) {
        // if a shorter travel time through 'there', update 'here'
        if ((delay + tt_there) < tt_here) {
          boxputglobal( *ttbox, here, delay + tt_there );
  	  changes++;
        }
        // if a shorter travel time through 'here', update 'there'
        else if ((delay + tt_here) < tt_there) {
          boxputglobal( *ttbox, there, delay + tt_here );
  	  changes++;
        }
      }
    }
  }
  
}


void
boxcudacpy_host(
  struct FLOATBOX *box,
  struct FLOATBOX **d_box,
  float **d_flat
)
{
  long flat_size = sizeof(float)*p3dcalcvolume(box->size);
  cudaMalloc(d_box, sizeof(struct FLOATBOX));
  cudaMalloc(d_flat, flat_size);

  cudaMemcpy(*d_box, box, sizeof(struct FLOATBOX), cudaMemcpyHostToDevice);
  cudaMemcpy(*d_flat, box->flat, flat_size, cudaMemcpyHostToDevice);
  cudaMemcpy(&((*d_box)->flat), d_flat, sizeof(float*), cudaMemcpyHostToDevice);
}

void
boxcudacpy_dev(
  struct FLOATBOX **d_box,
  struct FLOATBOX *box,
  float **d_flat
)
{ 
  long flat_size = sizeof(float) *p3dcalcvolume(box->size);
  cudaMemcpy(box->flat, *d_flat, flat_size, cudaMemcpyDeviceToHost);
}


extern "C"
long
do_sweep (
  struct STATE *state
)
{ 
  struct FLOATBOX vbox = state->vbox;
  struct FLOATBOX ttbox = state->ttbox;
  struct POINT3D ttstart = state->ttstart;
  const struct FORWARDSTAR * const star = state->star;
  const int numinstar = state->numinstar;

  struct FLOATBOX *d_vbox;
  float *d_vflat;
  boxcudacpy_host(&vbox, &d_vbox, &d_vflat);
  
  struct FLOATBOX *d_ttbox;
  float *d_ttflat;
  boxcudacpy_host(&ttbox, &d_ttbox, &d_ttflat);
   
  struct POINT3D *d_ttstart;
  cudaMalloc(&d_ttstart, sizeof(struct POINT3D));
  cudaMemcpy(d_ttstart, &ttstart, sizeof(struct POINT3D), cudaMemcpyHostToDevice);
  
  struct FORWARDSTAR *d_star;
  cudaMalloc(&d_star, sizeof(struct FORWARDSTAR));
  cudaMemcpy(d_star, star, sizeof(struct FORWARDSTAR), cudaMemcpyHostToDevice);
  
  int *d_numinstar;
  cudaMalloc(&d_numinstar, sizeof(int));
  cudaMemcpy(d_numinstar, &numinstar, sizeof(int), cudaMemcpyHostToDevice);
  
  long anychange = 0;   
  
  cudaMemcpyToSymbol(changes, &anychange, sizeof(long));

  dim3 grid(61, 61, 3);
  dim3 block(4, 4, 17);

  do_sweep_cuda<<<1, block>>>(d_vbox, d_ttbox, d_ttstart, d_star, d_numinstar);
  
  boxcudacpy_dev(&d_ttbox, &ttbox, &d_ttflat);
  cudaMemcpyFromSymbol(&anychange, changes, sizeof(long));
 
  printf("CHANGES: %ld\n", anychange); 
 
  cudaFree(d_vbox);
  cudaFree(d_vflat);
  cudaFree(d_ttbox);
  cudaFree(d_ttflat);
  cudaFree(d_ttstart);
  cudaFree(d_star);
  cudaFree(d_numinstar);
  
  return anychange;
}

// long
// do_sweep (
//   struct STATE *state
// )
/*int 
main()
{
  // copy some state into local stack memory for fastness
  struct POINT3D omin = {1, -5, 2};
  struct POINT3D omax = {241, 241, 51};
  struct POINT3D imin = {121, 121, 5};
  struct POINT3D imax = {241, 241, 51};

  struct FLOATBOX vbox;
  boxalloc(&vbox, omin, omax, imin, imax);
  vbox.flat[10] = 20.0f;  
  
  struct FLOATBOX ttbox;
  boxalloc(&ttbox, omin, omax, imin, imax);
  ttbox.flat[10] = 5.0f;

  struct FLOATBOX *d_vbox;
  float *d_vflat;
  boxcudacpy_host(&vbox, &d_vbox, &d_vflat);
  
  struct FLOATBOX *d_ttbox;
  float *d_ttflat;
  boxcudacpy_host(&ttbox, &d_ttbox, &d_ttflat);

  printf("0: %f\n", vbox.flat[10]);  
  printf("1: %f\n", ttbox.flat[10]);
  // count how many (if any) values we change
  // long changes = 0;
  // long *d_changes;
  // 
  // cudaMalloc(&d_changes, sizeof(long));
  // cudaMemcpy(d_changes, &changes, sizeof(long), cudaMemcpyHostToDevice);
   
  do_sweep_cuda<<<1, 32>>>(d_vbox, d_ttbox);
  
  boxcudacpy_dev(&d_vbox, &vbox, &d_vflat);
  boxcudacpy_dev(&d_ttbox, &ttbox, &d_ttflat);

  printf("0: %f\n", vbox.flat[10]);
  printf("1: %f\n", ttbox.flat[10]); 
  
  cudaFree(d_vbox);
  cudaFree(d_vflat);
  boxfree(&vbox);
  cudaFree(d_ttbox);
  cudaFree(d_ttflat);
  boxfree(&ttbox);
  //cudaFree(d_changes); 
 // return changes;
}*/


////////////////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////////////////
