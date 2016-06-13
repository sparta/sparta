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

  int periodic[3];
  int pflag = domain->periodic(periodic);
  if (!pflag and me == 0)
    error->warning(FLERR,"Grid is not periodic for compute fft/grid");

  if (grid->maxlevel != 1) 
    error->all(FLERR,"Compute fft/grid require uniform one-level grid");
  if (grid->nsplit) 
    error->all(FLERR,"Compute fft/grid cannot use grid with split cells");

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

    } else error->all(FLERR,"Illegal compute fft/grid command");

    nvalues++;
    iarg++;
  }

  // setup and error check

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

  per_grid_flag = 1;
  if (nvalues == 1) size_per_grid_cols = 0;
  else size_per_grid_cols = nvalues;

  // create FFT grid, partition it, and allocate grid/FFT memory
  // NOTE: could avoid allocating inbuf in some cases, depends on values

  fft_create();

  memory->create(complexbuf,2*nfft,"fft/grid:complexbuf");
  memory->create(fftbuf,nfft,"fft/grid:fftbuf");

  nglocal = grid->nlocal;
  memory->create(inbuf,nglocal,"fft/grid:inbuf");
  memory->create(outbuf,nglocal,"fft/grid:inbuf");

  vector_grid = NULL;
  array_grid = NULL;

  if (nvalues == 1) memory->create(vector_grid,nglocal,"fft/grid:vector_grid");
  else memory->create(array_grid,nglocal,nvalues,"fft/grid:array_grid");

  // initialize pointers that will be allocated in reallocate()

  irregular1 = irregular2 = NULL;
  map1 = map2 = NULL;
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

  memory->destroy(complexbuf);
  memory->destroy(fftbuf);

  memory->destroy(inbuf);
  memory->destroy(outbuf);
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
  int i,m,n;
  double *inbufptr;

  invoked_per_grid = update->ntimestep;

  // check that grid has not adapted

  if (grid->maxlevel != 1) 
    error->all(FLERR,"Compute fft/grid require uniform one-level grid");

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
        c->post_process_grid(aidx,-1,1,NULL,NULL,NULL,1);
      
      if (aidx == 0 || c->post_process_grid_flag) {
        inbufptr = c->vector_grid;
      } else {
        double **carray = c->array_grid;
        int n = grid->nlocal;
        int aidxm1 = aidx - 1;
        for (i = 0; i < n; i++) inbuf[i] = carray[i][aidxm1];
        inbufptr = inbuf;
      }

    // access fix fields, check if fix frequency is a match

    } else if (which[m] == FIX) {
      Fix *fix = modify->fix[vidx];

      if (update->ntimestep % modify->fix[vidx]->per_grid_freq)
        error->all(FLERR,"Fix used in compute fft/grid not "
                   "computed at compatible time");
      if (aidx == 0) {
        inbufptr = fix->vector_grid;
      } else {
        double **farray = fix->array_grid;
        int n = grid->nlocal;
        int aidxm1 = aidx - 1;
        for (i = 0; i < n; i++) inbuf[i] = farray[i][aidxm1];
        inbufptr = inbuf;
      }

    // evaluate particle-style or grid-style variable

    } else if (which[m] == VARIABLE) {
      input->variable->compute_grid(vidx,inbuf,1,0);
      inbufptr = inbuf;
    }

    // ------------------------------
    // perform FFT on inbufptr values
    // ------------------------------

    // irregular comm to move values from SPARTA grid -> FFT grid

    //debug("SEND INBUF",nglocal,inbufptr,NULL,NULL);

    irregular1->exchange_uniform((char *) inbufptr,sizeof(double),
                                 (char *) fftbuf);

    //debug("RECV FFTBUF",nfft,fftbuf,NULL,NULL);

    // convert scalar to FFT complex

    for (i = 0; i < nfft; i++) {
      n = 2*map1[i];
      complexbuf[n] = fftbuf[i];
      complexbuf[n+1] = ZEROF;
    }

    //debug("COMPLEX BUF",nfft,complexbuf,NULL,NULL,2);

    // perform FFT

    if (dimension == 3) fft3d->compute(complexbuf,complexbuf,1);
    else fft2d->compute(complexbuf,complexbuf,1);

    // NOTE: reverse FFT is just for debugging

    if (dimension == 3) fft3d->compute(complexbuf,complexbuf,-1);
    else fft2d->compute(complexbuf,complexbuf,-1);

    // convert FFT complex to scalar

    n = 0;
    for (i = 0; i < nfft; i++) {
      fftbuf[i] = complexbuf[n];
      n += 2;
    }

    //debug("SEND FFTBUF",nfft,fftbuf,NULL,NULL);

    // reverse irregular comm to move results from FFT grid -> SPARTA grid

    irregular2->exchange_uniform((char *) fftbuf,sizeof(double),
                                 (char *) outbuf);

    //debug("RECV OUTBUF",nglocal,outbuf,NULL,NULL);

    // copy output to vector/array via map2

    if (nvalues == 1) {
      int n = grid->nlocal;
      for (i = 0; i < n; i++) 
        vector_grid[map2[i]] = outbuf[i];
    } else {
      int n = grid->nlocal;
      for (i = 0; i < n; i++) 
        array_grid[map2[i]][m] = outbuf[i];
    }

    // may need to use values from array_grid
    //debug("FINAL DUMPBUF",nglocal,vector_grid,NULL,NULL);
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

  if (grid->nlocal == nglocal) return;

  memory->destroy(inbuf);
  memory->destroy(outbuf);
  memory->destroy(vector_grid);
  memory->destroy(array_grid);

  nglocal = grid->nlocal;

  memory->create(inbuf,grid->nlocal,"fft/grid:inbuf");
  memory->create(outbuf,grid->nlocal,"fft/grid:inbuf");

  if (nvalues == 1) memory->create(vector_grid,nglocal,"fft/grid:vector_grid");
  else memory->create(array_grid,nglocal,nvalues,"fft/grid:array_grid");
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based data
------------------------------------------------------------------------- */

bigint ComputeFFTGrid::memory_usage()
{
  bigint bytes = 0;
  bytes += 3*nfft * sizeof(double);           // complexbuf, fftbuf
  bytes += 2*nglocal * sizeof(double);        // inbuf, outbuf
  bytes += nvalues*nglocal * sizeof(double);  // vector/array grid
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

  dimension = domain->dimension;

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
  
  //printf("FFT %d: nxyz %d %d %d np xyz %d %d %d: x %d %d y %d %d z %d %d: %d\n",
  //         me,nx,ny,nz,npx,npy,npz,nxlo,nxhi,nylo,nyhi,nzlo,nzhi,nfft);

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

  // debug

  //cellint *ids;
  //memory->create(ids,grid->nlocal,"fft/grid:proclist1");
  //for (i = 0; i < grid->nlocal; i++) ids[i] = grid->cells[i].id;
  //debug("IDS",grid->nlocal,NULL,NULL,ids);

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

  //debug("PROCLIST1",nglocal,NULL,proclist1,NULL);

  int nrecv = irregular1->create_data_uniform(nglocal,proclist1);

  if (nrecv != nfft) 
    error->one(FLERR,"Compute fft/grid FFT mapping is inconsistent");

  memory->create(sbuf1,nglocal*sizeof(cellint),"fft/grid:sbuf1");
  memory->create(rbuf1,nfft*sizeof(cellint),"fft/grid:rbuf1");

  cellint *idsend = (cellint *) sbuf1;
  for (i = 0; i < nglocal; i++) idsend[i] = cells[i].id;

  //debug("SBUF1",nglocal,NULL,NULL,idsend);

  irregular1->exchange_uniform(sbuf1,sizeof(cellint),rbuf1);

  memory->create(map1,nfft,"fft/grid:map1");

  cellint *idrecv = (cellint *) rbuf1;

  //debug("RBUF1",nfft,NULL,NULL,idrecv);

  for (i = 0; i < nfft; i++) {
    gid = idrecv[i];
    ix = (gid-1) % nx;
    iy = ((gid-1) / nx) % ny;
    iz = (gid-1) / (nx*ny);
    map1[i] = (iz-nzlo)*nxfft*nyfft + (iy-nylo)*nxfft + (ix-nxlo);
  }

  //debug("MAP1",nfft,NULL,map1,NULL);

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

  //debug("PROCLIST2",nfft,NULL,proclist2,NULL);

  nrecv = irregular2->create_data_uniform(nfft,proclist2);
  if (nrecv != nglocal) 
    error->one(FLERR,"Compute fft/grid FFT mapping is inconsistent");

  memory->create(sbuf2,nfft*sizeof(cellint),"fft/grid:sbuf2");
  memory->create(rbuf2,nglocal*sizeof(cellint),"fft/grid:rbuf2");

  idsend = (cellint *) sbuf2;
  for (i = 0; i < nfft; i++) idsend[map1[i]] = idrecv[i];

  //debug("SBUF2",nfft,NULL,NULL,idsend);

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

  //debug("RBUF2",nglocal,NULL,NULL,idrecv);

  memory->create(map2,nglocal,"fft/grid:map1");
  for (i = 0; i < nglocal; i++) {
    gid = idrecv[i];
    map2[i] = (*hash)[gid] - 1;
  }

  //debug("MAP2",nglocal,NULL,map2,NULL);

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
