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

#include "string.h"
#include "stdlib.h"
#include "compute_fft_grid.h"
#include "update.h"
#include "domain.h"
#include "grid.h"
#include "modify.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "irregular.h"
#include "fft3d_wrap.h"
#include "fft2d_wrap.h"
#include "memory.h"
#include "error.h"

#ifdef SPARTA_MAP
#include <map>
#elif defined SPARTA_UNORDERED_MAP
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

using namespace SPARTA_NS;

enum{COMPUTE,FIX,VARIABLE};
enum{FORWARD=1,BACWARD=-1};

#define INVOKED_PER_GRID 16

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

ComputeFFTGrid::ComputeFFTGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute fft/grid command");

  // errors and warnings

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  dimension = domain->dimension;

  int periodic[3];
  int pflag = domain->periodic(periodic);
  if (!pflag and me == 0)
    error->warning(FLERR,"Grid is not periodic for compute fft/grid");

  if (grid->maxlevel != 1)
    error->all(FLERR,"Compute fft/grid require uniform one-level grid");
  if (grid->nsplit)
    error->all(FLERR,"Compute fft/grid cannot use grid with split cells");
  if (grid->unx % 2 || grid->uny % 2)
    error->all(FLERR,"Compute fft/grid cannot use grid "
               "with odd cell count in a dimension");
  if (dimension == 3 && grid->unz % 2)
    error->all(FLERR,"Compute fft/grid cannot use grid "
               "with odd cell count in a dimension");

  // parse input values

  nvalues = narg - 2;
  which = new int[nvalues];
  argindex = new int[nvalues];
  ids = new char*[nvalues];
  value2index = new int[nvalues];

  nvalues = 0;
  int iarg = 2;
  while (iarg < narg) {
    if (strncmp(arg[iarg],"c_",2) == 0 ||
        strncmp(arg[iarg],"f_",2) == 0 ||
        strncmp(arg[iarg],"v_",2) == 0) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal compute fft/grid command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      delete [] suffix;

    } else break;

    nvalues++;
    iarg++;
  }

  // optional args

  zeroflag = 0;
  sumflag = 0;
  scalefactor = 1.0;
  conjugate = 0;
  kx = ky = kz = kmag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"zero") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute fft/grid command");
      if (strcmp(arg[iarg+1],"yes") == 0) zeroflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) zeroflag = 0;
      else error->all(FLERR,"Illegal compute fft/grid command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"sum") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute fft/grid command");
      if (strcmp(arg[iarg+1],"yes") == 0) sumflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) sumflag = 0;
      else error->all(FLERR,"Illegal compute fft/grid command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute fft/grid command");
      scalefactor = input->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"conjugate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute fft/grid command");
      if (strcmp(arg[iarg+1],"yes") == 0) conjugate = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) conjugate = 0;
      else error->all(FLERR,"Illegal compute fft/grid command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"kx") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute fft/grid command");
      if (strcmp(arg[iarg+1],"yes") == 0) kx = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) kx = 0;
      else error->all(FLERR,"Illegal compute fft/grid command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"ky") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute fft/grid command");
      if (strcmp(arg[iarg+1],"yes") == 0) ky = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) ky = 0;
      else error->all(FLERR,"Illegal compute fft/grid command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"kz") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute fft/grid command");
      if (strcmp(arg[iarg+1],"yes") == 0) kz = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) kz = 0;
      else error->all(FLERR,"Illegal compute fft/grid command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"kmag") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute fft/grid command");
      if (strcmp(arg[iarg+1],"yes") == 0) kmag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) kmag = 0;
      else error->all(FLERR,"Illegal compute fft/grid command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute fft/grid command");
  }

  // setup and error check

  if (kz && dimension == 2)
    error->all(FLERR,"Compute fft/grid cannot use kz for 2d simulation");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute fft/grid does not exist");
      if (!modify->compute[icompute]->per_grid_flag)
        error->all(FLERR,"Compute fft/grid requires a per-grid compute");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_per_grid_cols != 0)
        error->all(FLERR,"Compute fft/grid compute does not "
                   "calculate a per-grid vector");
      if (argindex[i] && modify->compute[icompute]->size_per_grid_cols == 0)
        error->all(FLERR,"Compute fft/grid compute does not "
                   "calculate a per-grid array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_per_grid_cols)
        error->all(FLERR,
                   "Compute fft/grid compute array is accessed out-of-range");
    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute fft/grid does not exist");
      if (!modify->fix[ifix]->per_grid_flag)
        error->all(FLERR,"Compute fft/grid requires a per-grid fix");
      if (argindex[i] == 0 &&
          modify->fix[ifix]->size_per_grid_cols != 0)
        error->all(FLERR,"Compute fft/grid fix does not "
                   "calculate a per-grid vector");
      if (argindex[i] && modify->fix[ifix]->size_per_grid_cols == 0)
        error->all(FLERR,"Compute fft/grid fix does not "
                   "calculate a per-grid array");
      if (argindex[i] &&
          argindex[i] > modify->fix[ifix]->size_per_grid_cols)
        error->all(FLERR,"Compute fft/rgid fix array is accessed out-of-range");
    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for compute fft/grid does not exist");
      if (!input->variable->grid_style(ivariable))
        error->all(FLERR,"Compute fft/grid requires a grid-style variable");
    }
  }

  // compute settings
  // ncol = # of columns of output

  per_grid_flag = 1;

  if (sumflag) ncol = 1;
  else ncol = nvalues;
  if (conjugate == 0) ncol *= 2;
  startcol = kx + ky + kz + kmag;
  ncol += startcol;

  if (ncol == 1) size_per_grid_cols = 0;
  else size_per_grid_cols = ncol;

  // partition for FFTs
  // allocate bufs for grid and FFT decomps
  // NOTE: could avoid allocating inbuf in some cases, depends on values

  fft_create();

  memory->create(fft,2*nfft,"fft/grid:fft");
  memory->create(fftwork,nfft,"fft/grid:fftwork");

  irregular1 = irregular2 = NULL;
  map1 = map2 = NULL;
  ingrid = gridwork = NULL;
  gridworkcomplex = NULL;
  vector_grid = NULL;
  array_grid = NULL;

  nglocal = 0;
  reallocate();
}

/* ---------------------------------------------------------------------- */

ComputeFFTGrid::~ComputeFFTGrid()
{
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] which;
  delete [] argindex;
  delete [] ids;
  delete [] value2index;

  if (dimension == 3) delete fft3d;
  else delete fft2d;

  memory->destroy(fft);
  memory->destroy(fftwork);

  memory->destroy(ingrid);
  memory->destroy(gridwork);
  memory->destroy(gridworkcomplex);

  memory->destroy(vector_grid);
  memory->destroy(array_grid);

  delete irregular1;
  delete irregular2;
  memory->destroy(map1);
  memory->destroy(map2);
}

/* ---------------------------------------------------------------------- */

void ComputeFFTGrid::init()
{
  // check that grid has not adapted
  // check that grid still has no split cells

  if (grid->maxlevel != 1)
    error->all(FLERR,"Compute fft/grid require uniform one-level grid");
  if (grid->nsplit)
    error->all(FLERR,"Compute fft/grid cannot use grid with split cells");

  // create two irregular comm patterns for moving data
  //   from/to SPARTA grid to/from FFT grid
  // do this at each init in case grid partitioning has changed
  // also reallocate grid-based memory if needed

  reallocate();

  // set indices of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute fft/grid does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute fft/grid does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for compute fft/grid does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeFFTGrid::compute_per_grid()
{
  int i,j,m,n,icol;
  double real,imag;
  double *ingridptr;

  invoked_per_grid = update->ntimestep;

  // check that grid has not adapted
  // NOTE: also need to check it has not been re-balanced?

  if (grid->maxlevel != 1)
    error->all(FLERR,"Compute fft/grid require uniform one-level grid");

  // if sumflag set, zero output vector/array, but not K-space indices
  // so can sum each value's result into it

  if (sumflag) {
    if (ncol == 1) {
      for (i = 0; i < nglocal; i++)
        vector_grid[i] = 0.0;
    } else {
      for (i = 0; i < nglocal; i++)
        for (m = startcol; m < ncol; m++)
          array_grid[i][m] = 0.0;
    }
  }

  // process values, one FFT per value

  for (m = 0; m < nvalues; m++) {
    int vidx = value2index[m];
    int aidx = argindex[m];

    // invoke compute if not previously invoked
    // for per-grid compute, invoke post_process_grid() if necessary

    if (which[m] == COMPUTE) {
      Compute *c = modify->compute[vidx];

      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        c->compute_per_grid();
        c->invoked_flag |= INVOKED_PER_GRID;
      }

      if (c->post_process_grid_flag)
        c->post_process_grid(aidx,1,NULL,NULL,NULL,1);

      if (aidx == 0 || c->post_process_grid_flag) {
        ingridptr = c->vector_grid;
      } else {
        double **carray = c->array_grid;
        int n = grid->nlocal;
        int aidxm1 = aidx - 1;
        for (i = 0; i < n; i++) ingrid[i] = carray[i][aidxm1];
        ingridptr = ingrid;
      }

    // access fix fields, check if fix frequency is a match

    } else if (which[m] == FIX) {
      Fix *fix = modify->fix[vidx];

      if (update->ntimestep % modify->fix[vidx]->per_grid_freq)
        error->all(FLERR,"Fix used in compute fft/grid not "
                   "computed at compatible time");
      if (aidx == 0) {
        ingridptr = fix->vector_grid;
      } else {
        double **farray = fix->array_grid;
        int n = grid->nlocal;
        int aidxm1 = aidx - 1;
        for (i = 0; i < n; i++) ingrid[i] = farray[i][aidxm1];
        ingridptr = ingrid;
      }

    // evaluate particle-style or grid-style variable

    } else if (which[m] == VARIABLE) {
      input->variable->compute_grid(vidx,ingrid,1,0);
      ingridptr = ingrid;
    }

    // ------------------------------
    // perform FFT on inbufptr values
    // ------------------------------

    // irregular comm to move grid values from SPARTA owners -> FFT owners

    irregular1->exchange_uniform((char *) ingridptr,sizeof(double),
                                 (char *) fftwork);

    // convert SPARTA grid value to FFT complex
    // use map1 to morph recvbuf to local FFT layout

    for (i = 0; i < nfft; i++) {
      n = 2*map1[i];
      fft[n] = fftwork[i];
      fft[n+1] = ZEROF;
    }

    // perform forward FFT

    if (dimension == 3) fft3d->compute(fft,fft,FORWARD);
    else fft2d->compute(fft,fft,FORWARD);

    // reverse irregular comm to move results from FFT grid -> SPARTA grid
    // if conjugate set:
    //   convert complex FFT datums back to real via c times c*
    //   comm single floating point value per grid cell
    // else: comm entire complex value
    // copy or sum received values into output vec or array via map2

    if (conjugate) {
      n = 0;
      for (i = 0; i < nfft; i++) {
        real = fft[n];
        imag = fft[n+1];
        fftwork[i] = real*real + imag*imag;
        n += 2;
      }

      irregular2->exchange_uniform((char *) fftwork,sizeof(double),
                                   (char *) gridwork);

      if (sumflag) {
        if (ncol == 1) {
          int n = grid->nlocal;
          for (i = 0; i < n; i++)
            vector_grid[map2[i]] += gridwork[i];
        } else {
          icol = startcol;
          int n = grid->nlocal;
          for (i = 0; i < n; i++)
            array_grid[map2[i]][icol] += gridwork[i];
        }
      } else {
        if (ncol == 1) {
          int n = grid->nlocal;
          for (i = 0; i < n; i++)
            vector_grid[map2[i]] = gridwork[i];
        } else {
          icol = m + startcol;
          int n = grid->nlocal;
          for (i = 0; i < n; i++)
            array_grid[map2[i]][icol] = gridwork[i];
        }
      }

    } else {
      irregular2->exchange_uniform((char *) fft,2*sizeof(FFT_SCALAR),
                                   (char *) gridworkcomplex);

      if (sumflag) {
        icol = startcol;
        int n = grid->nlocal;
        j = 0;
        for (i = 0; i < n; i++) {
          array_grid[map2[i]][icol] += gridworkcomplex[j];
          array_grid[map2[i]][icol+1] += gridworkcomplex[j+1];
          j += 2;
        }
      } else {
        icol = 2*m + startcol;
        int n = grid->nlocal;
        j = 0;
        for (i = 0; i < n; i++) {
          array_grid[map2[i]][icol] = gridworkcomplex[j];
          array_grid[map2[i]][icol+1] = gridworkcomplex[j+1];
          j += 2;
        }
      }
    }
  }

  // scale results if requested

  if (scalefactor == 1.0) return;

  if (ncol == 1) {
    for (i = 0; i < nglocal; i++)
      vector_grid[i] *= scalefactor;
  } else {
    for (i = 0; i < nglocal; i++)
      for (m = startcol; m < ncol; m++)
        array_grid[i][m] *= scalefactor;
  }
}

/* ----------------------------------------------------------------------
   reallocate irregular comm patterns if local grid storage changes
   called by init() and whenever grid is rebalanced
------------------------------------------------------------------------- */

void ComputeFFTGrid::reallocate()
{
  delete irregular1;
  delete irregular2;
  memory->destroy(map1);
  memory->destroy(map2);

  irregular_create();

  memory->destroy(ingrid);
  memory->destroy(gridwork);
  memory->destroy(gridworkcomplex);
  memory->destroy(vector_grid);
  memory->destroy(array_grid);

  nglocal = grid->nlocal;

  memory->create(ingrid,nglocal,"fft/grid:ingrid");
  gridwork = NULL;
  gridworkcomplex = NULL;

  if (startcol || conjugate)
    memory->create(gridwork,nglocal,"fft/grid:gridwork");
  if (!conjugate) memory->create(gridworkcomplex,2*nglocal,
                                 "fft/grid:gridworkcomplex");

  if (ncol == 1) memory->create(vector_grid,nglocal,"fft/grid:vector_grid");
  else memory->create(array_grid,nglocal,ncol,"fft/grid:array_grid");

  // one-time setup of vector of K-space vector magnitudes if requested
  // compute values in K-space, irregular comm to grid decomposition
  // kx,ky,kz = indices of FFT grid cell in K-space
  // convert to distance from (0,0,0) cell using PBC
  // klen = length of K-space vector

  if (!startcol) return;

  int i,j,k;
  double ikx,iky,ikz;
  double klen;

  int nxhalf = nx/2;
  int nyhalf = ny/2;
  int nzhalf = nz/2;
  if (dimension == 2) nzhalf = 1;

  int icol = 0;

  for (int m = 0; m < 4; m++) {
    if (m == 0 && !kx) continue;
    if (m == 1 && !ky) continue;
    if (m == 2 && !kz) continue;
    if (m == 3 && !kmag) continue;

    int n = 0;
    for (k = nzlo; k <= nzhi; k++) {
      if (k < nzhalf) ikz = k;
      else ikz = nz - k;
      for (j = nylo; j <= nyhi; j++) {
        if (j < nyhalf) iky = j;
        else iky = ny - j;
        for (i = nxlo; i <= nxhi; i++) {
          if (i < nxhalf) ikx = i;
          else ikx = nx - i;

          if (m == 0) klen = ikx;
          else if (m == 1) klen = iky;
          else if (m == 2) klen = ikz;
          else if (m == 3) klen = sqrt(ikx*ikx + iky*iky + ikz*ikz);

          fftwork[n++] = klen;
        }
      }
    }

    irregular2->exchange_uniform((char *) fftwork,sizeof(double),
                                 (char *) gridwork);

    for (i = 0; i < nglocal; i++)
      array_grid[map2[i]][icol] = gridwork[i];

    icol++;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

bigint ComputeFFTGrid::memory_usage()
{
  bigint bytes = 0;
  bytes += 2*nfft * sizeof(FFT_SCALAR);       // fft
  bytes += nfft * sizeof(double);             // fftwork
  bytes += nglocal * sizeof(double);          // ingrid
  if (conjugate) bytes += nglocal * sizeof(double);       // gridwork
  else bytes += 2*nglocal * sizeof(FFT_SCALAR);           // gridworkcomplex
  bytes += ncol*nglocal * sizeof(double);     // vector/array grid
  bytes += (nfft+nglocal) * sizeof(int);      // map1,map2
  return bytes;
}

/* ----------------------------------------------------------------------
   one-time creation and partitioning of the FFT grid
------------------------------------------------------------------------- */

void ComputeFFTGrid::fft_create()
{
  // create FFT grid and partition it
  // proc partitioning = Npx by Npy by Npz
  // this proc's index within partitioning = xme,yme,zme
  // sub-domain this proc owns (inclusive): (Nxlo,Nxhi) (Nylo,Nyhi) (Nzlo,Nzhi)
  // nfft = # of FFT points on this proc

  nx = grid->unx;
  ny = grid->uny;
  nz = grid->unz;

  // warn if any grid dimension is not factorable by 2,3,5

  int flag = 0;
  if (!factorable(nx)) flag = 1;
  if (!factorable(ny)) flag = 1;
  if (dimension == 3 && !factorable(nz)) flag = 1;
  if (flag && me == 0)
    error->warning(FLERR,"Compute fft/grid FFT grid is not "
                   "factorable by 2,3,5");

  if (dimension == 3) {
    npx = 1;
    if (nz >= nprocs) {
      npy = 1;
      npz = nprocs;
    } else procs2grid2d(nprocs,ny,nz,npy,npz);
  } else {
    npx = 1;
    npy = nprocs;
    npz = 1;
  }

  int xme = me % npx;
  int yme = (me/npx) % npy;
  int zme = me / (npx*npy);

  nxlo = xme*nx/npx;
  nxhi = (xme+1)*nx/npx - 1;
  nylo = yme*ny/npy;
  nyhi = (yme+1)*ny/npy - 1;
  nzlo = zme*nz/npz;
  nzhi = (zme+1)*nz/npz - 1;

  nxfft = nxhi - nxlo + 1;
  nyfft = nyhi - nylo + 1;
  nzfft = nzhi - nzlo + 1;

  nfft = nxfft * nyfft * nzfft;

  //printf("FFT %d: nxyz %d %d %d np xyz %d %d %d: "
  //       "x %d %d y %d %d z %d %d: %d\n",
  //       me,nx,ny,nz,npx,npy,npz,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,nfft);

  bigint nfft2 = nfft;
  nfft2 *= 2;
  if (nfft2 > MAXSMALLINT)
    error->all(FLERR,"Compute fft/grid FFT is too large per-processor");

  // create FFT plan

  int collective_flag;
#ifdef __bg__
  collective_flag = 1;
#else
  collective_flag = 0;
#endif

  int tmp;
  if (dimension == 3) {
    fft3d = new FFT3D(sparta,world,nx,ny,nz,
                      nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
                      nxlo,nxhi,nylo,nyhi,nzlo,nzhi,
                      0,0,&tmp,collective_flag);
  } else {
    fft2d = new FFT2D(sparta,world,nx,ny,
                      nxlo,nxhi,nylo,nyhi,nxlo,nxhi,nylo,nyhi,
                      0,0,&tmp,collective_flag);
  }
}

/* ----------------------------------------------------------------------
   create two irregular comm patterns
   first for comm of data from SPARTA grid -> FFT grid
   second for comm of data from FFT grid -> SPARTA grid
------------------------------------------------------------------------- */

void ComputeFFTGrid::irregular_create()
{
  int i,ix,iy,iz,ipy,ipz;
  cellint gid;
  int *proclist1,*proclist2,*proclist3;
  char *sbuf1,*rbuf1,*sbuf2,*rbuf2;

  // plan for moving data from SPARTA grid -> FFT grid
  // send my cell IDs in SPARTA grid->cells order
  //   to proc who owns each in FFT partition
  // create map1 = how to reorder them on receiving FFT proc

  irregular1 = new Irregular(sparta);

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  memory->create(proclist1,nglocal,"fft/grid:proclist1");

  // use cell ID to determine which proc owns it in FFT partitioning
  // while loops are for case where some procs own nothing in FFT partition,
  //   due to more procs that rows/columns of FFT partition
  //   in this case formula (iy/ny * npy) may undershoot correct ipy proc
  //   increment until find ipy consistent with
  //   ipy*ny/npy to (ipy+1)*ny/npy bounds of that proc's FFT partition

  for (i = 0; i < nglocal; i++) {
    gid = cells[i].id;
    iy = ((gid-1) / nx) % ny;
    iz = (gid-1) / (nx*ny);

    ipy = static_cast<int> (1.0*iy/ny * npy);
    while (1) {
      if (iy >= ipy*ny/npy && iy < (ipy+1)*ny/npy) break;
      ipy++;
    }
    ipz = static_cast<int> (1.0*iz/nz * npz);
    while (1) {
      if (iz >= ipz*nz/npz && iz < (ipz+1)*nz/npz) break;
      ipz++;
    }

    proclist1[i] = ipz*npy + ipy;
  }

  int nrecv = irregular1->create_data_uniform(nglocal,proclist1);

  if (nrecv != nfft)
    error->one(FLERR,"Compute fft/grid FFT mapping is inconsistent");

  memory->create(sbuf1,nglocal*sizeof(cellint),"fft/grid:sbuf1");
  memory->create(rbuf1,nfft*sizeof(cellint),"fft/grid:rbuf1");

  cellint *idsend = (cellint *) sbuf1;
  for (i = 0; i < nglocal; i++) idsend[i] = cells[i].id;

  irregular1->exchange_uniform(sbuf1,sizeof(cellint),rbuf1);

  memory->create(map1,nfft,"fft/grid:map1");

  cellint *idrecv = (cellint *) rbuf1;

  for (i = 0; i < nfft; i++) {
    gid = idrecv[i];
    ix = (gid-1) % nx;
    iy = ((gid-1) / nx) % ny;
    iz = (gid-1) / (nx*ny);
    map1[i] = (iz-nzlo)*nxfft*nyfft + (iy-nylo)*nxfft + (ix-nxlo);
  }

  // plan for moving data from FFT grid -> SPARTA grid
  // send my cell IDs in FFT grid order
  //   back to proc who owns each in SPARTA grid->cells
  // proclist3 generated from irregular1->reverse()
  //   must permute proclist3 -> proclist2 via map1
  //   same way that received grid data will be permuted by map1 to FFT data
  // create map2 = how to reorder received FFT data on SPARTA grid proc

  irregular2 = new Irregular(sparta);

  memory->create(proclist3,nfft,"fft/grid:proclist2");
  irregular1->reverse(nrecv,proclist3);

  memory->create(proclist2,nfft,"fft/grid:proclist2");
  for (i = 0; i < nfft; i++) proclist2[map1[i]] = proclist3[i];

  nrecv = irregular2->create_data_uniform(nfft,proclist2);
  if (nrecv != nglocal)
    error->one(FLERR,"Compute fft/grid FFT mapping is inconsistent");

  memory->create(sbuf2,nfft*sizeof(cellint),"fft/grid:sbuf2");
  memory->create(rbuf2,nglocal*sizeof(cellint),"fft/grid:rbuf2");

  idsend = (cellint *) sbuf2;
  for (i = 0; i < nfft; i++) idsend[map1[i]] = idrecv[i];

  irregular2->exchange_uniform(sbuf2,sizeof(cellint),rbuf2);

  // insure grid cell IDs are hashed, so can use them to build map2

  if (!grid->hashfilled) grid->rehash();

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  idrecv = (cellint *) rbuf2;

  memory->create(map2,nglocal,"fft/grid:map1");
  for (i = 0; i < nglocal; i++) {
    gid = idrecv[i];
    map2[i] = (*hash)[gid];
  }

  // clean up

  memory->destroy(proclist1);
  memory->destroy(proclist2);
  memory->destroy(proclist3);
  memory->destroy(sbuf1);
  memory->destroy(rbuf1);
  memory->destroy(sbuf2);
  memory->destroy(rbuf2);
}

/* ----------------------------------------------------------------------
   map nprocs to NX by NY grid as PX by PY procs - return optimal px,py
------------------------------------------------------------------------- */

void ComputeFFTGrid::procs2grid2d(int nprocs, int nx, int ny, int &px, int &py)
{
  // loop thru all possible factorizations of nprocs
  // surf = surface area of largest proc sub-domain
  // innermost if test minimizes surface area and surface/volume ratio

  int bestsurf = 2 * (nx + ny);
  int bestboxx = 0;
  int bestboxy = 0;

  int boxx,boxy,surf,ipx,ipy;

  ipx = 1;
  while (ipx <= nprocs) {
    if (nprocs % ipx == 0) {
      ipy = nprocs/ipx;
      boxx = nx/ipx;
      if (nx % ipx) boxx++;
      boxy = ny/ipy;
      if (ny % ipy) boxy++;
      surf = boxx + boxy;
      if (surf < bestsurf ||
          (surf == bestsurf && boxx*boxy > bestboxx*bestboxy)) {
        bestsurf = surf;
        bestboxx = boxx;
        bestboxy = boxy;
        px = ipx;
        py = ipy;
      }
    }
    ipx++;
  }
}

/* ----------------------------------------------------------------------
   check if all factors of n are in list of factors
   return 1 if yes, 0 if no
------------------------------------------------------------------------- */

int ComputeFFTGrid::factorable(int n)
{
  int nfactors = 3;
  int factors[3] = {2,3,5};

  int i;
  while (n > 1) {
    for (i = 0; i < nfactors; i++) {
      if (n % factors[i] == 0) {
        n /= factors[i];
        break;
      }
    }
    if (i == nfactors) return 0;
  }

  return 1;
}

/* ----------------------------------------------------------------------
   debug by printing out vectors involved in grid <-> FFT remapping
------------------------------------------------------------------------- */

void ComputeFFTGrid::debug(const char *str, int n,
                           double *dx, int *ix, cellint *cx, int stride)
{
  int i,j;
  char buf[1024],one[32];

  sprintf(buf,"%s: %d %d:",str,me,n);
  if (dx) {
    j = 0;
    for (i = 0; i < n; i++) {
      sprintf(one," %g",dx[j]);
      strcat(buf,one);
      j += stride;
    }
  }
  if (ix) {
    j = 0;
    for (i = 0; i < n; i++) {
      sprintf(one," %d",ix[j]);
      strcat(buf,one);
      j += stride;
    }
  }
  if (cx) {
    j = 0;
    for (i = 0; i < n; i++) {
      sprintf(one," %d",cx[j]);
      strcat(buf,one);
      j += stride;
    }
  }
  sprintf(one,"\n");
  strcat(buf,one);
  printf("%s",buf);
}
