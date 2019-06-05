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
#include "create_grid.h"
#include "grid.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "random_mars.h"
#include "random_park.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NONE,LEVEL,STRIDE,CLUMP,BLOCK,RANDOM};
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

  grid->exist = 1;

  if (narg < 3) error->all(FLERR,"Illegal create_grid command");

  int nx = atoi(arg[0]);
  int ny = atoi(arg[1]);
  int nz = atoi(arg[2]);

  if (nx < 1 || ny < 1 || nz < 1)
    error->all(FLERR,"Illegal create_grid command");
  if (domain->dimension == 2 && nz != 1)
    error->all(FLERR,"Create_grid nz value must be 1 for a 2d simulation");

  // optional args

  dimension = domain->dimension;

  int nlevels = 1;
  int bstyle = NONE;
  int px = 0;
  int py = 0;
  int pz = 0;
  int order;
  int inside = ANY;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"level") == 0) {
      if (iarg+8 > narg) error->all(FLERR,"Illegal create_grid command");
      if (bstyle != NONE && bstyle != LEVEL)
        error->all(FLERR,"Illegal create_grid command");
      bstyle = LEVEL;
      if (atoi(arg[iarg+1]) != nlevels+1) 
        error->all(FLERR,"Illegal create_grid command");
      nlevels++;
      iarg += 8;

    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal create_grid command");
      if (bstyle != NONE && bstyle != LEVEL)
        error->all(FLERR,"Illegal create_grid command");
      bstyle = LEVEL;
      if (atoi(arg[iarg+1]) != nlevels+1) 
        error->all(FLERR,"Illegal create_grid command");
      if (domain->find_region(arg[iarg+2]) < 0)
        error->all(FLERR,"Create_grid region ID does not exist");
      nlevels++;
      iarg += 6;

    } else if (strcmp(arg[iarg],"stride") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_grid command");
      if (bstyle != NONE) error->all(FLERR,"Illegal create_grid command");
      bstyle = STRIDE;
      if (strcmp(arg[iarg+1],"xyz") == 0) order = XYZ;
      else if (strcmp(arg[iarg+1],"xzy") == 0) order = XZY;
      else if (strcmp(arg[iarg+1],"yxz") == 0) order = YXZ;
      else if (strcmp(arg[iarg+1],"yzx") == 0) order = YZX;
      else if (strcmp(arg[iarg+1],"zxy") == 0) order = ZXY;
      else if (strcmp(arg[iarg+1],"zyx") == 0) order = ZYX;
      else error->all(FLERR,"Illegal create_grid command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"clump") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_grid command");
      if (bstyle != NONE) error->all(FLERR,"Illegal create_grid command");
      bstyle = CLUMP;
      if (strcmp(arg[iarg+1],"xyz") == 0) order = XYZ;
      else if (strcmp(arg[iarg+1],"xzy") == 0) order = XZY;
      else if (strcmp(arg[iarg+1],"yxz") == 0) order = YXZ;
      else if (strcmp(arg[iarg+1],"yzx") == 0) order = YZX;
      else if (strcmp(arg[iarg+1],"zxy") == 0) order = ZXY;
      else if (strcmp(arg[iarg+1],"zyx") == 0) order = ZYX;
      else error->all(FLERR,"Illegal create_grid command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"block") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal create_grid command");
      if (bstyle != NONE) error->all(FLERR,"Illegal create_grid command");
      bstyle = BLOCK;
      if (strcmp(arg[iarg+1],"*") == 0) px = 0;
      else px = atoi(arg[iarg+1]);
      if (strcmp(arg[iarg+2],"*") == 0) py = 0;
      else py = atoi(arg[iarg+2]);
      if (strcmp(arg[iarg+3],"*") == 0) pz = 0;
      else pz = atoi(arg[iarg+3]);
      iarg += 4;

    } else if (strcmp(arg[iarg],"random") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal create_grid command");
      if (bstyle != NONE) error->all(FLERR,"Illegal create_grid command");
      bstyle = RANDOM;
      iarg += 1;

    } else if (strcmp(arg[iarg],"inside") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal create_grid command");
      if (strcmp(arg[iarg+1],"any") == 0) inside = ANY;
      else if (strcmp(arg[iarg+1],"all") == 0) inside = ALL;
      else error->all(FLERR,"Illegal create_grid command");
      iarg += 2;

    } else error->all(FLERR,"Illegal create_grid command");
  }
  
  if (bstyle == NONE) bstyle = LEVEL;

  // partition domain across procs, only for BLOCK style

  if (bstyle == BLOCK) {
    procs2grid(nx,ny,nz,px,py,pz);
    if (px*py*pz != comm->nprocs)
      error->all(FLERR,"Bad grid of processors for create_grid");
  }

  // create root parent cell
  // treat first 3 args as if specified as "level 1 * * * Nx Ny Nz"

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  int level = 1;
  int xlo,xhi,ylo,yhi,zlo,zhi;
  xlo = xhi = ylo = yhi = zlo = zhi = 1;
  iarg = 3;
  Region *region = NULL;

  // loop over levels
  // new level determines assignment of previous-level cells
  //   to be parent cells vs child cells
  // parent cells are those which are further partitioned by new level
  // child cells are those that are not
  // add all cells in previous level as either parent or child cells
  // if this is last level, also add all current level cells as child cells

  int me = comm->me;
  int nprocs = comm->nprocs;
  bigint count = 0;

  int pnx,pny,pnz,ix,iy,iz,nbits,pflag,proc;
  cellint m,nth,idgrandparent,idparent,idchild;
  double lo[3],hi[3];
  Grid::ParentCell *p;

  while (1) {

    // add previous level cells as parent or child cells
    // loop over all parent cells to find ones two levels up
    // use their info to generate parent or child cells at previous level
    // decision on parent vs child in previous level depends on 
    //   pxyz lo/hi bounds in this level
    //   or on whether parent is in/out of region

    if (level == 1) {
      grid->add_parent_cell(0,-1,nx,ny,nz,domain->boxlo,domain->boxhi);

    } else {
      int nparent = grid->nparent;
      int prevlevel = level-2;
      
      for (int igrandparent = 0; igrandparent < nparent; igrandparent++) {
        if (grid->pcells[igrandparent].level != prevlevel) continue;
        p = &grid->pcells[igrandparent];

        idgrandparent = p->id;
        nbits = p->nbits;
        pnx = p->nx;
        pny = p->ny;
        pnz = p->nz;

        m = 0;
        for (iz = 0; iz < pnz; iz++)
          for (iy = 0; iy < pny; iy++)
            for (ix = 0; ix < pnx; ix++) {
              m++;
              idparent = idgrandparent | (m << nbits);
              grid->id_child_lohi(igrandparent,m,lo,hi);
              if (region) pflag = cell_in_region(lo,hi,region,inside);
              else {
                pflag = 1;
                if (ix+1 < xlo || ix+1 > xhi) pflag = 0;
                if (iy+1 < ylo || iy+1 > yhi) pflag = 0;
                if (iz+1 < zlo || iz+1 > zhi) pflag = 0;
              }
              if (pflag)
                grid->add_parent_cell(idparent,igrandparent,nx,ny,nz,lo,hi);
              else {
                if (count % nprocs == me)
                  grid->add_child_cell(idparent,igrandparent,lo,hi);
                count++;
              }
            }
      }
    }

    // final level, add current level cells as child cells
    // loop over all parent cells to find ones at previous level
    // use their info to generate my child cells at this level
    // if BSTYLE is set, there is only 1 level, create proc's cells directly

    if (level == nlevels) {
      Grid::ParentCell *pcells = grid->pcells;
      int nparent = grid->nparent;
      int prevlevel = level-1;
      
      for (int iparent = 0; iparent < nparent; iparent++) {
        if (pcells[iparent].level != prevlevel) continue;
        p = &pcells[iparent];
        idparent = p->id;
        nbits = p->nbits;
        nx = p->nx;
        ny = p->ny;
        nz = p->nz;

        if (bstyle == LEVEL) {
          cellint ntotal = (cellint) nx * ny * nz;
          int firstproc = count % nprocs;
          cellint ifirst = me - firstproc + 1;
          if (ifirst <= 0) ifirst += nprocs;
          for (m = ifirst; m <= ntotal; m += nprocs) {
            idchild = idparent | (m << nbits);
            grid->id_child_lohi(iparent,m,lo,hi);
            grid->add_child_cell(idchild,iparent,lo,hi);
          }
          count += ntotal;

        // loop over all child cells
        // convert M to Nth based on order
        // assign each cell to proc based on Nth and STRIDE or CLUMP

        } else if (bstyle == STRIDE || bstyle == CLUMP) {
          cellint ntotal = (cellint) nx * ny * nz;
          for (m = 0; m < ntotal; m++) {
            ix = m % nx;
            iy = (m / nx) % ny;
            iz = m / ((cellint) nx*ny);
            if (order == XYZ) 
              nth = (cellint) iz*nx*ny + (cellint) iy*nx + ix;
            else if (order == XZY) 
              nth = (cellint) iy*nx*nz + (cellint) iz*nx + ix;
            else if (order == YXZ) 
              nth = (cellint) iz*ny*nx + (cellint) ix*ny + iy;
            else if (order == YZX) 
              nth = (cellint) ix*ny*nz + (cellint) iz*ny + iy;
            else if (order == ZXY) 
              nth = (cellint) iy*nz*nx + (cellint) ix*nz + iz;
            else if (order == ZYX) 
              nth = (cellint) ix*nz*ny + (cellint) iy*nz + iz;
            nth++;
            if (bstyle == STRIDE) proc = nth % nprocs;
            else proc = static_cast<int> (1.0*nth/ntotal * nprocs);
            if (proc != me) continue;
            idchild = idparent | (nth << nbits);
            grid->id_child_lohi(iparent,nth,lo,hi);
            grid->add_child_cell(idchild,iparent,lo,hi);
          }
          count += ntotal;

        // loop over subset of cells in my BLOCK

        } else if (bstyle == BLOCK) {
          int ipx = me % px;
          int ipy = (me / px) % py;
          int ipz = me / (px*py);

          int ixstart = static_cast<int> (1.0*ipx/px*nx);
          int ixstop = static_cast<int> (1.0*(ipx+1)/px*nx);
          int iystart = static_cast<int> (1.0*ipy/py*ny);
          int iystop = static_cast<int> (1.0*(ipy+1)/py*ny);
          int izstart = static_cast<int> (1.0*ipz/pz*nz);
          int izstop = static_cast<int> (1.0*(ipz+1)/pz*nz);

          for (iz = izstart; iz < izstop; iz++) {
            for (iy = iystart; iy < iystop; iy++) {
              for (ix = ixstart; ix < ixstop; ix++) {
                m = (cellint) iz*nx*ny + (cellint) iy*nx + ix;
                m++;
                idchild = idparent | (m << nbits);
                grid->id_child_lohi(iparent,m,lo,hi);
                grid->add_child_cell(idchild,iparent,lo,hi);
              }
            }
          }

        // loop over all child cells, assign randomly to a proc

        } else if (bstyle == RANDOM) {
          RanPark *random = new RanPark(update->ranmaster->uniform());
          cellint ntotal = (cellint) nx * ny * nz;
          for (m = 0; m < ntotal; m++) {
            proc = static_cast<int> (nprocs*random->uniform());
            if (proc != me) continue;
            idchild = idparent | ((m+1) << nbits);
            grid->id_child_lohi(iparent,m+1,lo,hi);
            grid->add_child_cell(idchild,iparent,lo,hi);
          }
          count += ntotal;
          delete random;
        }

      }
      break;
    }
    
    if (level == nlevels) break;

    // args for next level
    // levels must be in ascending order starting at beginning of arg list

    level++;

    if (strcmp(arg[iarg],"level") == 0) {
      bounds(arg[iarg+2],nx,xlo,xhi);
      bounds(arg[iarg+3],ny,ylo,yhi);
      bounds(arg[iarg+4],nz,zlo,zhi);
      nx = atoi(arg[iarg+5]);
      ny = atoi(arg[iarg+6]);
      nz = atoi(arg[iarg+7]);
      iarg += 8;
    } else if (strcmp(arg[iarg],"region") == 0) {
      int iregion = domain->find_region(arg[iarg+2]);
      region = domain->regions[iregion];
      nx = atoi(arg[iarg+3]);
      ny = atoi(arg[iarg+4]);
      nz = atoi(arg[iarg+5]);
      iarg += 6;
    } else error->all(FLERR,"Illegal create_grid command");

    if (nx < 1 || ny < 1 || nz < 1)
      error->all(FLERR,"Illegal create_grid command");
    if (dimension == 2) {
      if (zlo != 1 || zhi != 1 || nz != 1) 
        error->all(FLERR,"Illegal create_grid command");
    }
  }

  // set grandparent flag for all parent cells
  
  Grid::ParentCell *pcells = grid->pcells;
  int nparent = grid->nparent;

  for (int i = 1; i < nparent; i++)
    pcells[pcells[i].iparent].grandparent = 1;

  // invoke grid methods to complete grid setup

  if (nprocs == 1 || bstyle == CLUMP || bstyle == BLOCK) grid->clumped = 1;
  else grid->clumped = 0;

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->find_neighbors();
  grid->check_uniform();
  comm->reset_neighbors();

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // stats

  double time_total = time3-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Created " BIGINT_FORMAT " child grid cells\n",
              grid->ncell);
      fprintf(screen,"  parent cells = %d\n",grid->nparent);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  create/ghost percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
    }

    if (logfile) {
      fprintf(logfile,"Created " BIGINT_FORMAT " child grid cells\n",
              grid->ncell);
      fprintf(logfile,"  parent cells = %d\n",grid->nparent);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  create/ghost percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
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

  //printf("BOUNDS %s nmax %d nlo/hi %d %d\n",str,nmax,nlo,nhi);

  if (nlo < 1 || nhi > nmax || nlo > nhi) 
    error->all(FLERR,"Numeric index is out of bounds");
}

/* ----------------------------------------------------------------------
   test if grid cell defined by lo/hi is inside region
   one corner point inside is enough for inside = ANY
   all corner points must be inside for inside = ALL
   return 1 if inside, 0 if not
------------------------------------------------------------------------- */

int CreateGrid::cell_in_region(double *lo, double *hi, 
                               Region *region, int inside)
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

