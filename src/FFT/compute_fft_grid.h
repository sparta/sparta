/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
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
  int zeroflag,conjugate,sumflag,kx,ky,kz,kmag;
  double scalefactor;
  int ncol,startcol;

  int dimension;
  int nglocal;                         // # of owned grid decomp values
  int maxgrid;                         // max size of grid bufs
  int nx,ny,nz;                        // size of global regular grid
  int npx,npy,npz;                     // # of procs in each dim in FFT decomp
  int nxlo,nxhi,nylo,nyhi,nzlo,nzhi;   // bounds of FFT grid on this proc
  int nxfft,nyfft,nzfft;               // extent of FFT grid on this proc
  int nfft;                            // # of owned FFT decomp values

  double *ingrid;      // input grid values from compute,fix,variable
                       // may be NULL if ingridptr just points to c/f/v
  double *fftwork;     // work buf in FFT decomp, length = nfft
  double *gridwork;    // work buf in grid decomp, length = nglocal

  FFT_SCALAR *fft;     // complex buf for performing FFT, length = nfft
  FFT_SCALAR *gridworkcomplex;  // work buf in grid decomp, length = nglocal

  int *map1;            // mapping of received SPARTA grid values to FFT grid
                        // map1[i] = index into ordered FFT grid of
                        //           Ith value in buffer received
                        //           from SPARTA decomp via irregular comm
  int *map2;            // mapping of received FFT grid values to SPARTA grid
                        // map2[i] = index into SPARTA grid of Ith value
                        //           in buffer received from FFT decomp via
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
