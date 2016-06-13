/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(fft/grid,ComputeFFTGrid)

#else

#ifndef SPARTA_COMPUTE_FFT_GRID_H
#define SPARTA_COMPUTE_FFT_GRID_H

#include "compute.h"

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

namespace SPARTA_NS {

class ComputeFFTGrid : public Compute {
 public:
  ComputeFFTGrid(class SPARTA *, int, char **);
  ~ComputeFFTGrid();
  void init();
  void compute_per_grid();
  void reallocate();
  bigint memory_usage();

 private:
  int me,nprocs;
  int nvalues;
  int *which,*argindex,*value2index;
  char **ids;

  int nglocal;
  int dimension;
  int nx,ny,nz;
  int npx,npy,npz;
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;
  int nfft,nxfft,nyfft,nzfft;

  FFT_SCALAR *complexbuf;  // complex buffer to perform FFT in
  double *inbuf;        // buffer of grid values from compute,fix,variable
                        // may not be used, if point directly to compute/fix
  double *outbuf;       // buffer of grid values with result of FFT
  int maxgrid;          // max size of in/out buf
  double *fftbuf;       // buffer of grid values in FFT layout

  int *map1;            // mapping of received SPARTA grid values to FFT grid
                        // map1[i] = index into ordered FFT grid of 
                        //           Ith value in buffer received
                        //           from SPARTA grid via irregular comm
  int *map2;            // mapping of received FFT grid values to SPARTA grid
                        // map2[i] = index into SPARTA grid of Ith value
                        //           in buffer received from FFT grid via
                        //           irregular comm

  class FFT3D *fft3d;
  class FFT2D *fft2d;
  class Irregular *irregular1,*irregular2;

  void fft_create();
  void irregular_create();
  void procs2grid2d(int, int, int, int &, int &);
  int factorable(int);
  void debug(const char *, int, double *, int *, cellint *, int stride=1);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute grid mixture ID does not exist

Self-explanatory.

E: Number of groups in compute grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

E: Invalid call to ComputeGrid::post_process_grid()

This indicates a coding error.  Please report the issue to the SPARTA
developers.

*/
