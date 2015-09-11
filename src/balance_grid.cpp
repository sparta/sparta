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

#include "stdlib.h"
#include "string.h"
#include "balance_grid.h"
#include "update.h"
#include "grid.h"
#include "domain.h"
#include "modify.h"
#include "comm.h"
#include "rcb.h"
#include "output.h"
#include "dump.h"
#include "random_mars.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

//#define RCB_DEBUG 1     // un-comment to include RCB proc boxes in image

enum{NONE,STRIDE,CLUMP,BLOCK,RANDOM,PROC,BISECTION};
enum{XYZ,XZY,YXZ,YZX,ZXY,ZYX};
enum{CELL,PARTICLE};

#define ZEROPARTICLE 0.1

/* ---------------------------------------------------------------------- */

BalanceGrid::BalanceGrid(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void BalanceGrid::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot balance grid before grid is defined");

  if (narg < 1) error->all(FLERR,"Illegal balance_grid command");

  int bstyle,order;
  int px,py,pz;
  int rcbwt,rcbflip;

  if (strcmp(arg[0],"none") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal balance_grid command");
    bstyle = NONE;

  } else if (strcmp(arg[0],"stride") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal balance_grid command");
    bstyle = STRIDE;
    if (strcmp(arg[1],"xyz") == 0) order = XYZ;
    else if (strcmp(arg[1],"xzy") == 0) order = XZY;
    else if (strcmp(arg[1],"yxz") == 0) order = YXZ;
    else if (strcmp(arg[1],"yzx") == 0) order = YZX;
    else if (strcmp(arg[1],"zxy") == 0) order = ZXY;
    else if (strcmp(arg[1],"zyx") == 0) order = ZYX;
    else error->all(FLERR,"Illegal balance_grid command");

  } else if (strcmp(arg[0],"clump") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal balance_grid command");
    bstyle = CLUMP;
    if (strcmp(arg[1],"xyz") == 0) order = XYZ;
    else if (strcmp(arg[1],"xzy") == 0) order = XZY;
    else if (strcmp(arg[1],"yxz") == 0) order = YXZ;
    else if (strcmp(arg[1],"yzx") == 0) order = YZX;
    else if (strcmp(arg[1],"zxy") == 0) order = ZXY;
    else if (strcmp(arg[1],"zyx") == 0) order = ZYX;
    else error->all(FLERR,"Illegal balance_grid command");

  } else if (strcmp(arg[0],"block") == 0) {
    if (narg != 4) error->all(FLERR,"Illegal balance_grid command");
    bstyle = BLOCK;
    if (strcmp(arg[1],"*") == 0) px = 0;
    else px = atoi(arg[1]);
    if (strcmp(arg[2],"*") == 0) py = 0;
    else py = atoi(arg[2]);
    if (strcmp(arg[3],"*") == 0) pz = 0;
    else pz = atoi(arg[3]);
    
  } else if (strcmp(arg[0],"random") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal balance_grid command");
    bstyle = RANDOM;

  } else if (strcmp(arg[0],"proc") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal balance_grid command");
    bstyle = PROC;

  } else if (strcmp(arg[0],"rcb") == 0) {
    if (narg != 2 && narg != 3) 
      error->all(FLERR,"Illegal balance_grid command");
    bstyle = BISECTION;
    if (strcmp(arg[1],"cell") == 0) rcbwt = CELL;
    else if (strcmp(arg[1],"part") == 0) rcbwt = PARTICLE;
    else error->all(FLERR,"Illegal balance_grid command");
    // undocumented optional 3rd arg
    // rcbflip = 3rd arg = 1 forces rcb->compute() to flip sign
    //           of all grid cell "dots" to force totally different
    //           assignment of grid cells to procs and induce
    //           complete rebalance data migration
    rcbflip = 0;
    if (narg == 3) rcbflip = atoi(arg[2]);

  } else error->all(FLERR,"Illegal balance_grid command");

  // error check on methods only allowed for a uniform grid

  if (bstyle == STRIDE || bstyle == CLUMP || bstyle == BLOCK)
    if (!grid->uniform) 
      error->all(FLERR,"Invalid balance_grid style for non-uniform grid");

  // re-assign each of my local child cells to a proc
  // only assign unsplit and split cells
  // do not assign sub-cells since they migrate with their split cell
  // set nmigrate = # of cells that will migrate to a new proc
  // reset proc field in cells for migrating cells
  // style NONE performs no re-assignment

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  int nprocs = comm->nprocs;
  int newproc;
  int nmigrate = 0;

  if (bstyle == STRIDE) {
    cellint idm1,ix,iy,iz,nth;

    cellint nx = grid->unx;
    cellint ny = grid->uny;
    cellint nz = grid->unz;

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      idm1 = cells[icell].id - 1;
      ix = idm1 % nx;
      iy = (idm1 / nx) % ny;
      iz = idm1 / (nx*ny);
    
      if (order == XYZ) nth = iz*nx*ny + iy*nx + ix;
      else if (order == XZY) nth = iy*nx*nz + iz*nx + ix;
      else if (order == YXZ) nth = iz*ny*nx + ix*ny + iy;
      else if (order == YZX) nth = ix*ny*nz + iz*ny + iy;
      else if (order == ZXY) nth = iy*nz*nx + ix*nz + iz;
      else if (order == ZYX) nth = ix*nz*ny + iy*nz + iz;

      newproc = nth % nprocs;

      if (newproc != cells[icell].proc) nmigrate++;
      cells[icell].proc = newproc;
    }

  } else if (bstyle == CLUMP) {
    cellint idm1,ix,iy,iz,nth;

    cellint nx = grid->unx;
    cellint ny = grid->uny;
    cellint nz = grid->unz;
    bigint ntotal = (bigint) nx * ny * nz;

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      idm1 = cells[icell].id - 1;
      ix = idm1 % nx;
      iy = (idm1 / nx) % ny;
      iz = idm1 / (nx*ny);
    
      if (order == XYZ) nth = iz*nx*ny + iy*nx + ix;
      else if (order == XZY) nth = iy*nx*nz + iz*nx + ix;
      else if (order == YXZ) nth = iz*ny*nx + ix*ny + iy;
      else if (order == YZX) nth = ix*ny*nz + iz*ny + iy;
      else if (order == ZXY) nth = iy*nz*nx + ix*nz + iz;
      else if (order == ZYX) nth = ix*nz*ny + iy*nz + iz;

      newproc = static_cast<int> (1.0*nth/ntotal * nprocs);

      if (newproc != cells[icell].proc) nmigrate++;
      cells[icell].proc = newproc;
    }

  } else if (bstyle == BLOCK) {
    cellint idm1,ix,iy,iz;
    int ipx,ipy,ipz;

    int nx = grid->unx;
    int ny = grid->uny;
    int nz = grid->unz;

    procs2grid(nx,ny,nz,px,py,pz);
    if (px*py*pz != nprocs)
      error->all(FLERR,"Bad grid of processors for balance_grid block");

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      idm1 = cells[icell].id - 1;
      ix = idm1 % nx;
      iy = (idm1 / nx) % ny;
      iz = idm1 / (nx*ny);
      ipx = ix*px / nx;
      ipy = iy*py / ny;
      ipz = iz*pz / nz;

      newproc = ipz*px*py + ipy*px + ipx;

      if (newproc != cells[icell].proc) nmigrate++;
      cells[icell].proc = newproc;
    }

  } else if (bstyle == RANDOM) {
    int newproc;
    RanPark *random = new RanPark(update->ranmaster->uniform());
    double seed = update->ranmaster->uniform();
    random->reset(seed,comm->me,100);

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      newproc = nprocs * random->uniform();
      if (newproc != cells[icell].proc) nmigrate++;
      cells[icell].proc = newproc;
    }

    delete random;

  } else if (bstyle == PROC) {
    int newproc;
    RanPark *random = new RanPark(update->ranmaster->uniform());
    newproc = nprocs * random->uniform();

    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      if (newproc != cells[icell].proc) nmigrate++;
      cells[icell].proc = newproc;
      newproc++;
      if (newproc == nprocs) newproc = 0;
    }

    delete random;

  } else if (bstyle == BISECTION) {
    RCB *rcb = new RCB(sparta);

    double **x;
    memory->create(x,nglocal,3,"balance_grid:x");

    double *lo,*hi;

    int nbalance = 0;
    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      lo = cells[icell].lo;
      hi = cells[icell].hi;
      x[nbalance][0] = 0.5*(lo[0]+hi[0]);
      x[nbalance][1] = 0.5*(lo[1]+hi[1]);
      x[nbalance][2] = 0.5*(lo[2]+hi[2]);
      nbalance++;
    }

    double *wt = NULL;
    if (rcbwt == PARTICLE) {
      particle->sort();
      int zero = 0;
      int n;
      memory->create(wt,nglocal,"balance_grid:wt");
      nbalance = 0;
      for (int icell = 0; icell < nglocal; icell++) {
        if (cells[icell].nsplit <= 0) continue;
        n = cinfo[icell].count;
        if (n) wt[nbalance++] = n;
        else {
          wt[nbalance++] = ZEROPARTICLE;
          zero++;
        }
      }
    }

    rcb->compute(nbalance,x,wt,rcbflip);

    // DEBUG info for dump image

#ifdef RCB_DEBUG

    update->rcblo[0] = rcb->lo[0];
    update->rcblo[1] = rcb->lo[1];
    update->rcblo[2] = rcb->lo[2];
    update->rcbhi[0] = rcb->hi[0];
    update->rcbhi[1] = rcb->hi[1];
    update->rcbhi[2] = rcb->hi[2];

#endif 

    rcb->invert();

    nbalance = 0;
    int *sendproc = rcb->sendproc;
    for (int icell = 0; icell < nglocal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      cells[icell].proc = sendproc[nbalance++];
    }
    nmigrate = nbalance - rcb->nkeep;

    delete rcb;
    memory->destroy(x);
    memory->destroy(wt);
  }

  // set clumped of not, depending on style
  // NONE style does not change clumping

  if (nprocs == 1 || bstyle == CLUMP || bstyle == BLOCK || bstyle == BISECTION) 
    grid->clumped = 1;
  else if (bstyle != NONE) grid->clumped = 0;

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // sort particles
  // NOTE: not needed again if rcbwt = PARTICLE for bstyle = BISECTION ??

  particle->sort();

  // DEBUG

  /*
  char file[32];
  sprintf(file,"tmp.bef.%d",comm->me);
  FILE *fp = fopen(file,"w");

  fprintf(fp,"Cells %d %d\n",grid->nlocal,grid->nghost);
  for (int i = 0; i < grid->nlocal+grid->nghost; i++) {
    fprintf(fp,"cell %d " CELLINT_FORMAT ": %d : %d %d %d %d %d %d\n",
           i,grid->cells[i].id,
           grid->cells[i].nmask,
           grid->cells[i].neigh[0],
           grid->cells[i].neigh[1],
           grid->cells[i].neigh[2],
           grid->cells[i].neigh[3],
           grid->cells[i].neigh[4],
           grid->cells[i].neigh[5]);
  }
  fclose(fp);
  */

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // invoke init() so all grid cell info, including collide & fixes,
  //   is ready to migrate
  // migrate grid cells and their particles to new owners
  // invoke grid methods to complete grid setup

  int ghost_previous = grid->exist_ghost;

  sparta->init();
  grid->unset_neighbors();
  grid->remove_ghosts();
  comm->migrate_cells(nmigrate);

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  grid->setup_owned();
  grid->acquire_ghosts();
  if (ghost_previous) grid->reset_neighbors();
  else grid->find_neighbors();
  comm->reset_neighbors();

  // reallocate per grid arrays in per grid dumps

  for (int i = 0; i < output->ndump; i++)
    output->dump[i]->reset_grid();

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  // DEBUG

  /*
  sprintf(file,"tmp.aft.%d",comm->me);
  fp = fopen(file,"w");

  fprintf(fp,"Cells %d %d\n",grid->nlocal,grid->nghost);
  for (int i = 0; i < grid->nlocal+grid->nghost; i++) {
    fprintf(fp,"cell %d " CELLINT_FORMAT ": %d : %d %d %d %d %d %d\n",
           i,grid->cells[i].id,
           grid->cells[i].nmask,
           grid->cells[i].neigh[0],
           grid->cells[i].neigh[1],
           grid->cells[i].neigh[2],
           grid->cells[i].neigh[3],
           grid->cells[i].neigh[4],
           grid->cells[i].neigh[5]);
  }
  fclose(fp);
  */

  // stats on balance operation

  bigint count = nmigrate;
  bigint nmigrate_all;
  MPI_Allreduce(&count,&nmigrate_all,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  double time_total = time5-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Balance grid migrated " BIGINT_FORMAT " cells\n",
              nmigrate_all);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  reassign/sort/migrate/ghost percent = %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total);
    }
    if (logfile) {
      fprintf(logfile,"Balance grid migrated " BIGINT_FORMAT " cells\n",
              nmigrate_all);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  reassign/sort/migrate/ghost percent = %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   assign nprocs to 3d grid so as to minimize surface area 
   area = surface area of each of 3 faces of simulation box
------------------------------------------------------------------------- */

void BalanceGrid::procs2grid(int nx, int ny, int nz, 
                             int &px, int &py, int &pz)
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
  area[0] = nx*ny;
  area[1] = nx*nz;
  area[2] = ny*nz;

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
      if (domain->dimension == 2 && ipz != 1) valid = 0;
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
