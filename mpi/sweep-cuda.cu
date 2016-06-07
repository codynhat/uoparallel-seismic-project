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

////////////////////////////////////////////////////////////////////////////////
// functions
////////////////////////////////////////////////////////////////////////////////

extern "C" long
do_sweep (
  struct STATE *state
)
{
  return 0;
}

/*void
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
*/

/*__global__
void
do_sweep_cuda(
  const struct FLOATBOX *vbox,
  const struct FLOATBOX *ttbox
  // const struct POINT3D *ttstart,
  // const struct FORWARDSTAR ** const star,
  // const int *numinstar,
  // long *changes
)
{
 //  *changes = 10; 
  ttbox->flat[10] *= 2; 
  vbox->flat[10] *= 3;

}*/

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
