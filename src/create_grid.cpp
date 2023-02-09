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

#include "stdlib.h"
#include "string.h"
#include "create_grid.h"
#include "grid.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NOSTYLE,BLOCK,CLUMP,RANDOM,STRIDE};
enum{NOLEVEL,SUBSET,REGION};
enum{XYZ,XZY,YXZ,YZX,ZXY,ZYX};
enum{ANY,ALL};

/* ---------------------------------------------------------------------- */

CreateGrid::CreateGrid(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void CreateGrid::command(int narg, char **arg)
{
  if (!domain->box_exist)
    error->all(FLERR,"Cannot create grid before simulation box is defined");
  if (grid->exist)
    error->all(FLERR,"Cannot create grid when grid is already defined");

  me = comm->me;
  nprocs = comm->nprocs;
  grid->exist = 1;
  dimension = domain->dimension;

  if (narg < 3) error->all(FLERR,"Illegal create_grid command");

  nx = atoi(arg[0]);
  ny = atoi(arg[1]);
  nz = atoi(arg[2]);

  // pstyle and args

  pstyle = NOSTYLE;
  nlevels = 1;
  levels = NULL;
  stack = NULL;
  inside = ANY;

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"block") == 0) {
      if (pstyle != NOSTYLE) error->all(FLERR,"Illegal create_grid command");
      if (iarg+4 > narg) error->all(FLERR,"Illegal create_grid command");
      pstyle = BLOCK;
      if (strcmp(arg[iarg+1],"*") == 0) px = 0;
      else px = atoi(arg[iarg+1]);
      if (strcmp(arg[iarg+2],"*") == 0) py = 0;
      else py = atoi(arg[iarg+2]);
      if (strcmp(arg[iarg+3],"*") == 0) pz = 0;
      else pz = atoi(arg[iarg+3]);
      iarg += 4;

    } else if (strcmp(arg[iarg],"clump") == 0) {
      if (pstyle != NOSTYLE) error->all(FLERR,"Illegal create_grid command");
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_grid command");
      pstyle = CLUMP;
      if (strcmp(arg[iarg+1],"xyz") == 0) order = XYZ;
      else if (strcmp(arg[iarg+1],"xzy") == 0) order = XZY;
      else if (strcmp(arg[iarg+1],"yxz") == 0) order = YXZ;
      else if (strcmp(arg[iarg+1],"yzx") == 0) order = YZX;
      else if (strcmp(arg[iarg+1],"zxy") == 0) order = ZXY;
      else if (strcmp(arg[iarg+1],"zyx") == 0) order = ZYX;
      else error->all(FLERR,"Illegal create_grid command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"stride") == 0) {
      if (pstyle != NOSTYLE) error->all(FLERR,"Illegal create_grid command");
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_grid command");
      pstyle = STRIDE;
      if (strcmp(arg[iarg+1],"xyz") == 0) order = XYZ;
      else if (strcmp(arg[iarg+1],"xzy") == 0) order = XZY;
      else if (strcmp(arg[iarg+1],"yxz") == 0) order = YXZ;
      else if (strcmp(arg[iarg+1],"yzx") == 0) order = YZX;
      else if (strcmp(arg[iarg+1],"zxy") == 0) order = ZXY;
      else if (strcmp(arg[iarg+1],"zyx") == 0) order = ZYX;
      else error->all(FLERR,"Illegal create_grid command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"random") == 0) {
      if (pstyle != NOSTYLE) error->all(FLERR,"Illegal create_grid command");
      if (iarg+1 > narg) error->all(FLERR,"Illegal create_grid command");
      pstyle = RANDOM;
      iarg += 1;

    } else if (strcmp(arg[iarg],"levels") == 0) {
      if (nlevels > 1) error->all(FLERR,"Illegal create_grid command");
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_grid command");
      nlevels = atoi(arg[iarg+1]);
      if (nlevels < 2) error->all(FLERR,"Create_grid nlevels must be > 1");
      if (nlevels > grid->plevel_limit)
        error->all(FLERR,"Create_grid nlevels exceeds MAXLEVEL");
      stack = new Stack[nlevels];
      levels = new Level[nlevels];
      for (int i = 0; i < nlevels; i++) levels[i].setflag = 0;
      levels[0].setflag = 1;
      levels[0].style = NOLEVEL;
      levels[0].cx = nx;
      levels[0].cy = ny;
      levels[0].cz = nz;
      iarg += 2;

    } else if (strcmp(arg[iarg],"subset") == 0) {
      if (nlevels == 1) error->all(FLERR,"Illegal create_grid command");
      if (iarg+8 > narg) error->all(FLERR,"Illegal create_grid command");
      int nlo,nhi;
      if (strchr(arg[iarg+1],'*')) {
        bounds(arg[iarg+1],nlevels,nlo,nhi);
        nlo = MAX(nlo,2);
      } else {
        nlo = nhi = atoi(arg[iarg+1]);
        if (nlo < 2) error->all(FLERR,"Create grid subset level < 2");
      }
      for (int i = nlo-1; i <= nhi-1; i++) {
        if (levels[i].setflag)
          error->all(FLERR,"Create_grid subset is resetting a level");
        levels[i].setflag = 1;
        levels[i].style = SUBSET;
        bounds(arg[iarg+2],levels[i-1].cx,levels[i].ixlo,levels[i].ixhi);
        bounds(arg[iarg+3],levels[i-1].cy,levels[i].iylo,levels[i].iyhi);
        bounds(arg[iarg+4],levels[i-1].cz,levels[i].izlo,levels[i].izhi);
        levels[i].cx = atoi(arg[iarg+5]);
        levels[i].cy = atoi(arg[iarg+6]);
        levels[i].cz = atoi(arg[iarg+7]);
      }
      iarg += 8;

    } else if (strcmp(arg[iarg],"region") == 0) {
      if (nlevels == 1) error->all(FLERR,"Illegal create_grid command");
      if (iarg+6 > narg) error->all(FLERR,"Illegal create_grid command");
      int nlo,nhi;
      if (strchr(arg[iarg+1],'*')) {
        bounds(arg[iarg+1],nlevels,nlo,nhi);
        nlo = MAX(nlo,2);
      } else {
        nlo = nhi = atoi(arg[iarg+1]);
        if (nlo < 2) error->all(FLERR,"Create grid region level < 2");
      }
      for (int i = nlo-1; i <= nhi-1; i++) {
        if (levels[i].setflag)
          error->all(FLERR,"Create_grid region is resetting a level");
        levels[i].setflag = 1;
        levels[i].style = REGION;
        int iregion = domain->find_region(arg[iarg+2]);
        if (iregion < 0) error->all(FLERR,"Create_grid region ID does not exist");
        levels[i].region = domain->regions[iregion];
        levels[i].cx = atoi(arg[iarg+3]);
        levels[i].cy = atoi(arg[iarg+4]);
        levels[i].cz = atoi(arg[iarg+5]);
      }
      iarg += 6;

    } else if (strcmp(arg[iarg],"inside") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_grid command");
      if (strcmp(arg[iarg+1],"any") == 0) inside = ANY;
      else if (strcmp(arg[iarg+1],"all") == 0) inside = ALL;
      else error->all(FLERR,"Illegal create_grid command");
      iarg += 2;

    } else error->all(FLERR,"Illegal create_grid command");
  }

  // default pstyle = BLOCK

  if (pstyle == NOSTYLE) {
    //pstyle = BLOCK;
    //px = py = pz = 0;
    pstyle = STRIDE;
    order = XYZ;
  }

  // check that grid levels are valid

  if (nx < 1 || ny < 1 || nz < 1)
    error->all(FLERR,"Illegal create_grid command");
  if (dimension == 2 && nz != 1)
    error->all(FLERR,"Create_grid nz value must be 1 for a 2d simulation");

  for (int i = 1; i < nlevels; i++) {
    if (!levels[i].setflag) error->all(FLERR,"Create_grid level was not set");
    if (dimension == 2 && levels[i].cz != 1)
      error->all(FLERR,"Create_grid cz value must be 1 for a 2d simulation");
    if (levels[i].cx < 1 || levels[i].cy < 1 || levels[i].cz < 1)
      error->all(FLERR,"Create_grid cx,cy,cz cannot be < 1");
    if (levels[i].cx == 1 && levels[i].cy == 1 && levels[i].cz == 1)
      error->all(FLERR,"Create_grid cx,cy,cz cannot all be one");
  }

  // transfer level info into Grid data structs

  Grid::ParentLevel *plevels = grid->plevels;
  grid->maxlevel = nlevels;

  plevels[0].nx = nx;
  plevels[0].ny = ny;
  plevels[0].nz = nz;
  plevels[0].nxyz = (bigint) nx * ny * nz;

  for (int i = 1; i < nlevels; i++) {
    plevels[i].nx = levels[i].cx;
    plevels[i].ny = levels[i].cy;
    plevels[i].nz = levels[i].cz;
    plevels[i].nxyz = (bigint) levels[i].cx * levels[i].cy * levels[i].cz;
  }

  for (int i = 0; i < nlevels; i++) {
    if (i == 0) plevels[i].nbits = 0;
    else plevels[i].nbits = plevels[i-1].nbits + plevels[i-1].newbits;
    plevels[i].newbits = grid->id_bits(plevels[i].nx,plevels[i].ny,plevels[i].nz);
  }

  // error check on too many bits for cell IDs

  int nbits = plevels[nlevels-1].nbits + plevels[nlevels-1].newbits;
  if (nbits > sizeof(cellint)*8) {
    char str[128];
    sprintf(str,"Hierarchical grid induces cell IDs that exceed %d bits",
            (int) sizeof(cellint)*8);
    error->all(FLERR,str);
  }

  // create grid with specified partitioning style

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  if (pstyle == BLOCK) create_block();
  else if (pstyle == CLUMP) create_clump();
  else if (pstyle == STRIDE) create_stride();
  else if (pstyle == RANDOM) create_random();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  // invoke grid methods to complete grid setup

  if (nprocs == 1 || pstyle == CLUMP || pstyle == BLOCK) grid->clumped = 1;
  else grid->clumped = 0;

  grid->set_maxlevel();
  grid->setup_owned();
  grid->acquire_ghosts();
  grid->find_neighbors();
  grid->check_uniform();
  comm->reset_neighbors();

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // clean up

  delete [] levels;
  delete [] stack;

  // stats

  double time_total = time3-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Created " BIGINT_FORMAT " child grid cells\n",
              grid->ncell);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  create/ghost percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
    }

    if (logfile) {
      fprintf(logfile,"Created " BIGINT_FORMAT " child grid cells\n",
              grid->ncell);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  create/ghost percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
    }
  }
}

/* ---------------------------------------------------------------------- */

void CreateGrid::create_block()
{
  int ix,iy,iz,level;
  cellint childID;
  double lo[3],hi[3];
  Stack *s;

  // partition single-level grid across procs

  procs2grid(nx,ny,nz,px,py,pz);
  if (px*py*pz != comm->nprocs)
    error->all(FLERR,"Bad grid of processors for create_grid block");

  // my subset of block

  int ipx = me % px;
  int ipy = (me / px) % py;
  int ipz = me / (px*py);

  int ixlo = static_cast<int> (1.0*ipx/px*nx);
  int ixhi = static_cast<int> (1.0*(ipx+1)/px*nx) - 1;
  int iylo = static_cast<int> (1.0*ipy/py*ny);
  int iyhi = static_cast<int> (1.0*(ipy+1)/py*ny) - 1;
  int izlo = static_cast<int> (1.0*ipz/pz*nz);
  int izhi = static_cast<int> (1.0*(ipz+1)/pz*nz) - 1;

  // create my subset of cells

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  for (iz = izlo; iz <= izhi; iz++) {
    for (iy = iylo; iy <= iyhi; iy++) {
      for (ix = ixlo; ix <= ixhi; ix++) {
        childID = (cellint) iz*nx*ny + (cellint) iy*nx + ix + 1;
        grid->id_child_lohi(0,boxlo,boxhi,childID,lo,hi);

        if (nlevels == 1) {
          grid->add_child_cell(childID,1,lo,hi);
        } else {
          s = &stack[0];
          s->id = childID;
          s->level = 1;
          s->lo[0] = lo[0]; s->lo[1] = lo[1]; s->lo[2] = lo[2];
          s->hi[0] = hi[0]; s->hi[1] = hi[1]; s->hi[2] = hi[2];
          recurse_levels(0);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void CreateGrid::create_clump()
{
  int ix,iy,iz,level;
  int i1,i2,i3,n1,n2,n3;
  cellint childID;
  double lo[3],hi[3];
  Stack *s;

  if (order == XYZ) {
    n1 = nx; n2 = ny; n3 = nz;
  } else if (order == XZY) {
    n1 = nx; n2 = nz; n3 = ny;
  } else if (order == YXZ) {
    n1 = ny; n2 = nx; n3 = nz;
  } else if (order == YZX) {
    n1 = ny; n2 = nz; n3 = nx;
  } else if (order == ZXY) {
    n1 = nz; n2 = nx; n3 = ny;
  } else if (order == ZYX) {
    n1 = nz; n2 = ny; n3 = nx;
  }

  // loop over my clump of ordered cells in requested order

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  cellint ntotal = (cellint) nx * ny * nz;
  cellint clumplo = static_cast<int> (1.0*me/nprocs * ntotal);
  cellint clumphi = static_cast<int> (1.0*(me+1)/nprocs * ntotal) - 1;

  for (cellint m = clumplo; m <= clumphi; m++) {
    i1 = m % n1;
    i2 = (m / n1) % n2;
    i3 = m / ((cellint) n1*n2);

    if (order == XYZ) {
      ix = i1; iy = i2; iz = i3;
    } else if (order == XZY) {
      ix = i1; iy = i3; iz = i2;
    } else if (order == YXZ) {
      ix = i2; iy = i1; iz = i3;
    } else if (order == YZX) {
      ix = i3; iy = i1; iz = i2;
    } else if (order == ZXY) {
      ix = i2; iy = i3; iz = i1;
    } else if (order == ZYX) {
      ix = i3; iy = i2; iz = i1;
    }

    childID = (cellint) iz*ny*nx + (cellint) iy*nx + ix + 1;
    grid->id_child_lohi(0,boxlo,boxhi,childID,lo,hi);

    if (nlevels == 1) {
      grid->add_child_cell(childID,1,lo,hi);
    } else {
      s = &stack[0];
      s->id = childID;
      s->level = 1;
      s->lo[0] = lo[0]; s->lo[1] = lo[1]; s->lo[2] = lo[2];
      s->hi[0] = hi[0]; s->hi[1] = hi[1]; s->hi[2] = hi[2];
      recurse_levels(0);
    }
  }
}

/* ---------------------------------------------------------------------- */

void CreateGrid::create_stride()
{
  int ix,iy,iz,level;
  int i1,i2,i3,n1,n2,n3;
  cellint childID;
  double lo[3],hi[3];
  Stack *s;

  if (order == XYZ) {
    n1 = nx; n2 = ny; n3 = nz;
  } else if (order == XZY) {
    n1 = nx; n2 = nz; n3 = ny;
  } else if (order == YXZ) {
    n1 = ny; n2 = nx; n3 = nz;
  } else if (order == YZX) {
    n1 = ny; n2 = nz; n3 = nx;
  } else if (order == ZXY) {
    n1 = nz; n2 = nx; n3 = ny;
  } else if (order == ZYX) {
    n1 = nz; n2 = ny; n3 = nx;
  }

  // stride over all ntotal cells in requested order

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  cellint ntotal = (cellint) nx * ny * nz;

  for (cellint m = me; m < ntotal; m += nprocs) {
    i1 = m % n1;
    i2 = (m / n1) % n2;
    i3 = m / ((cellint) n1*n2);

    if (order == XYZ) {
      ix = i1; iy = i2; iz = i3;
    } else if (order == XZY) {
      ix = i1; iy = i3; iz = i2;
    } else if (order == YXZ) {
      ix = i2; iy = i1; iz = i3;
    } else if (order == YZX) {
      ix = i3; iy = i1; iz = i2;
    } else if (order == ZXY) {
      ix = i2; iy = i3; iz = i1;
    } else if (order == ZYX) {
      ix = i3; iy = i2; iz = i1;
    }

    childID = (cellint) iz*ny*nx + (cellint) iy*nx + ix + 1;
    grid->id_child_lohi(0,boxlo,boxhi,childID,lo,hi);

    //printf("MMM m %d ixyz %d %d %d childID %d\n",m,ix,iy,iz,childID);

    if (nlevels == 1) {
      grid->add_child_cell(childID,1,lo,hi);
    } else {
      s = &stack[0];
      s->id = childID;
      s->level = 1;
      s->lo[0] = lo[0]; s->lo[1] = lo[1]; s->lo[2] = lo[2];
      s->hi[0] = hi[0]; s->hi[1] = hi[1]; s->hi[2] = hi[2];
      recurse_levels(0);
    }
  }
}

/* ---------------------------------------------------------------------- */

void CreateGrid::create_random()
{
  int proc,level;
  cellint childID;
  double lo[3],hi[3];
  Stack *s;

  RanKnuth *random = new RanKnuth(update->ranmaster->uniform());

  // loop over all ntotal cells
  // only add cell if this proc is the random owner

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  cellint ntotal = (cellint) nx * ny * nz;

  for (cellint m = 0; m < ntotal; m++) {
    proc = static_cast<int> (nprocs*random->uniform());
    if (proc != me) continue;

    childID = m+1;
    grid->id_child_lohi(0,boxlo,boxhi,childID,lo,hi);

    if (nlevels == 1) {
      grid->add_child_cell(childID,1,lo,hi);
    } else {
      s = &stack[0];
      s->id = childID;
      s->level = 1;
      s->lo[0] = lo[0]; s->lo[1] = lo[1]; s->lo[2] = lo[2];
      s->hi[0] = hi[0]; s->hi[1] = hi[1]; s->hi[2] = hi[2];
      recurse_levels(0);
    }
  }

  delete random;
}

/* ----------------------------------------------------------------------
   recurse through grid strucure to find and add child cells
   istack = index to current top of stack
------------------------------------------------------------------------- */

void CreateGrid::recurse_levels(int istack)
{
  int ix,iy,iz;
  double lo[3],hi[3];

  Grid::ParentLevel *plevels = grid->plevels;
  Stack *s = &stack[istack];

  int level = s->level;
  int nbits = plevels[level-1].nbits;
  cellint ichild = s->id >> nbits;

  int parentflag;
  if (level == nlevels) {
    parentflag = 0;
  } else if (levels[level].style == SUBSET) {
    parentflag = 1;
    ichild--;
    int nx = plevels[level-1].nx;
    int ny = plevels[level-1].ny;
    ix = ichild % nx;
    iy = (ichild / nx) % ny;
    iz = ichild / ((cellint) nx*ny);
    if (ix+1 < levels[level].ixlo) parentflag = 0;
    if (ix+1 > levels[level].ixhi) parentflag = 0;
    if (iy+1 < levels[level].iylo) parentflag = 0;
    if (iy+1 > levels[level].iyhi) parentflag = 0;
    if (iz+1 < levels[level].izlo) parentflag = 0;
    if (iz+1 > levels[level].izhi) parentflag = 0;
  } else if (levels[level].style == REGION) {
    parentflag = cell_in_region(s->lo,s->hi,levels[level].region);
  }

  if (!parentflag) {
    grid->add_child_cell(s->id,level,s->lo,s->hi);
    return;
  }

  int nx = plevels[level].nx;
  int ny = plevels[level].ny;
  int nz = plevels[level].nz;
  nbits = plevels[level].nbits;
  Stack *sc = &stack[istack+1];

  cellint m = 0;
  for (iz = 0; iz < nz; iz++) {
    for (iy = 0; iy < ny; iy++) {
      for (ix = 0; ix < nx; ix++) {
        m++;
        sc->id = (m << nbits) | s->id;
        sc->level = level+1;
        grid->id_child_lohi(level,s->lo,s->hi,m,sc->lo,sc->hi);
        recurse_levels(istack+1);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   1 = lower bound, nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = 1 to nmax,
     (3) i* = i to nmax, (4) *j = 1 to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void CreateGrid::bounds(char *str, int nmax, int &nlo, int &nhi)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = nhi = atoi(str);
  } else if (strlen(str) == 1) {
    nlo = 1;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = 1;
    nhi = atoi(ptr+1);
  } else if (strlen(ptr+1) == 0) {
    nlo = atoi(str);
    nhi = nmax;
  } else {
    nlo = atoi(str);
    nhi = atoi(ptr+1);
  }

  if (nlo < 1 || nhi > nmax || nlo > nhi)
    error->all(FLERR,"Create_grid subset index is out of bounds");
}

/* ----------------------------------------------------------------------
   test if grid cell defined by lo/hi is inside region
   one corner point inside is enough for inside = ANY
   all corner points must be inside for inside = ALL
   return 1 if inside, 0 if not
------------------------------------------------------------------------- */

int CreateGrid::cell_in_region(double *lo, double *hi, Region *region)
{
  double x[3];

  int n = 0;

  x[0] = lo[0]; x[1] = lo[1]; x[2] = lo[2];
  n += region->match(x);
  x[0] = hi[0]; x[1] = lo[1]; x[2] = lo[2];
  n += region->match(x);
  x[0] = lo[0]; x[1] = hi[1]; x[2] = lo[2];
  n += region->match(x);
  x[0] = hi[0]; x[1] = hi[1]; x[2] = lo[2];
  n += region->match(x);

  if (dimension == 3) {
    x[0] = lo[0]; x[1] = lo[1]; x[2] = hi[2];
    n += region->match(x);
    x[0] = hi[0]; x[1] = lo[1]; x[2] = hi[2];
    n += region->match(x);
    x[0] = lo[0]; x[1] = hi[1]; x[2] = hi[2];
    n += region->match(x);
    x[0] = hi[0]; x[1] = hi[1]; x[2] = hi[2];
    n += region->match(x);
  }

  if (inside == ANY) {
    if (n) return 1;
  } else if (inside == ALL) {
    if (dimension == 2 && n == 4) return 1;
    if (dimension == 3 && n == 8) return 1;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d grid so as to minimize surface area
   area = surface area of each of 3 faces of simulation box
------------------------------------------------------------------------- */

void CreateGrid::procs2grid(int nx, int ny, int nz, int &px, int &py, int &pz)
{
  int upx = px;
  int upy = py;
  int upz = pz;

  int nprocs = comm->nprocs;

  // all 3 proc counts are specified

  if (px && py && pz) return;

  // 2 out of 3 proc counts are specified

  if (py > 0 && pz > 0) {
    px = nprocs/(py*pz);
    return;
  } else if (px > 0 && pz > 0) {
    py = nprocs/(px*pz);
    return;
  } else if (px > 0 && py > 0) {
    pz = nprocs/(px*py);
    return;
  }

  // determine cross-sectional areas
  // area[0] = xy, area[1] = xz, area[2] = yz

  double area[3];
  area[0] = (bigint) nx*ny;
  area[1] = (bigint) nx*nz;
  area[2] = (bigint) ny*nz;

  double bestsurf = 2.0 * (area[0]+area[1]+area[2]);

  // loop thru all possible factorizations of nprocs
  // only consider valid cases that match procgrid settings
  // surf = surface area of a proc sub-domain

  int ipx,ipy,ipz,valid;
  double surf;

  ipx = 1;
  while (ipx <= nprocs) {
    valid = 1;
    if (upx && ipx != upx) valid = 0;
    if (nprocs % ipx) valid = 0;
    if (!valid) {
      ipx++;
      continue;
    }

    ipy = 1;
    while (ipy <= nprocs/ipx) {
      valid = 1;
      if (upy && ipy != upy) valid = 0;
      if ((nprocs/ipx) % ipy) valid = 0;
      if (!valid) {
        ipy++;
        continue;
      }

      ipz = nprocs/ipx/ipy;
      valid = 1;
      if (upz && ipz != upz) valid = 0;
      if (dimension == 2 && ipz != 1) valid = 0;
      if (!valid) {
        ipy++;
        continue;
      }

      surf = area[0]/ipx/ipy + area[1]/ipx/ipz + area[2]/ipy/ipz;
      if (surf < bestsurf) {
        bestsurf = surf;
        px = ipx;
        py = ipy;
        pz = ipz;
      }
      ipy++;
    }

    ipx++;
  }
}
