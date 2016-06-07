#include <stdio.h>

__global__
void
do_sweep_cuda(
  float *vbox
 // const struct FLOATBOX *ttbox,
 // const struct POINT3D *ttstart,
 // const struct FORWARDSTAR ** const star,
 // const int *numinstar_p,
 // long *changes
)
{
 vbox[threadIdx.x] *= 2;
 /* int numinstar = *numinstar_p;

  const struct POINT3D here = p3d( x, y, z );

  const float vel_here = boxgetglobal( *vbox, here );
  const float tt_here = boxgetglobal( *ttbox, here );  
  

  for( int l = 0; l < numinstar; l++ ) {

    // find point in forward star based on offsets
    const struct POINT3D there = p3daddp3d( here, *star[l].pos );

    // if 'there' is outside the boundaries, then skip
    if (
      p3disless( there, vbox->omin ) ||
      p3dismore( there, vbox->omax )
    ) {
      continue;
    }

    // compute delay from 'here' to 'there' with endpoint average
    const float vel_there = boxgetglobal( *vbox, there );
    const float delay = *star[l].halfdistance * (vel_here + vel_there);

    // ignore the starting point
    if( p3disnotequal( here, ttstart ) ) {

      const float tt_there = boxgetglobal( ttbox, there );

      // if offset point has infinity travel time, then update
      if ((tt_here == INFINITY) && (tt_there == INFINITY)) {
        continue;
      }

      if ((tt_here != INFINITY) && (tt_there == INFINITY)) {
        boxputglobal( *ttbox, there, delay + tt_here );
        *changes++;
        continue;
      }

      if ((tt_here == INFINITY) && (tt_there != INFINITY)) {
        boxputglobal( *ttbox, here, delay + tt_there );
        *changes++;
        continue;
      }

      if ((tt_here != INFINITY) && (tt_there != INFINITY)) {
        // if a shorter travel time through 'there', update 'here'
        if ((delay + tt_there) < tt_here) {
	  boxputglobal( *ttbox, here, delay + tt_there );
	  *changes++;
        }
        // if a shorter travel time through 'here', update 'there'
        else if ((delay + tt_here) < tt_there) {
	  boxputglobal( *ttbox, there, delay + tt_here );
	  *changes++;
        }
      }
    }
  }*/

}

long
do_sweep_cuda_init (
  const int N,
  float *vbox
)
{
  float *d_vbox;
  
  size_t size = N * sizeof(float);

  printf("1,1,1: %f", vbox[0]);
  
  
  cudaMalloc(&d_vbox, size); 
  cudaMemcpy(d_vbox, vbox, size, cudaMemcpyHostToDevice);

  do_sweep_cuda<<<1, 10>>>(d_vbox);
   
  cudaMemcpy(vbox, d_vbox, size, cudaMemcpyDeviceToHost);
  
  printf("1,1,1: %f", vbox[0]);
  return 0;
}


int
main (
  int argc,
  char *argv[]
)
{
  float vbox[10] = {1,2,3,4,5,6,7,8,9,10};
  do_sweep_cuda_init(10, vbox);  
}
