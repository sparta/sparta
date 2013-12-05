/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "create_grid.h"
#include "grid.h"
#include "domain.h"
#include "comm.h"
#include "error.h"

using namespace SPARTA_NS;

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

  // create root parent cell

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // treat first 3 args as if they were specified as "1 * * * Nx Ny Nz"

  int iarg = 3;
  int level = 1;
  int xlo,xhi,ylo,yhi,zlo,zhi;
  xlo = xhi = ylo = yhi = zlo = zhi = 1;

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

  int pnx,pny,pnz,ix,iy,iz,nbits,pflag;
  cellint m,idgrandparent,idparent,idchild;
  double lo[3],hi[3];
  Grid::ParentCell *p;

  while (1) {

    // add previous level cells as parent or child cells
    // loop over all parent cells to find ones two levels up
    // use their info to generate parent or child cells at previous level
    // decision on parent vs child in previous level depends on 
    //   pxyz lo/hi bounds in this level

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
              pflag = 1;
              if (ix+1 < xlo || ix+1 > xhi) pflag = 0;
              if (iy+1 < ylo || iy+1 > yhi) pflag = 0;
              if (iz+1 < zlo || iz+1 > zhi) pflag = 0;
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

    if (iarg == narg) {
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

        /*
        m = 0;
        for (iz = 0; iz < nz; iz++)
          for (iy = 0; iy < ny; iy++)
            for (ix = 0; ix < nx; ix++) {
              m++;
              idchild = idparent | (m << nbits);
              if (count % nprocs == me) {
                grid->id_child_lohi(iparent,m,lo,hi);
                grid->add_child_cell(idchild,iparent,lo,hi);
              }
              count++;
            }
        */
      }
      break;
    }
    
    level++;

    // args for next level

    if (iarg+7 > narg) error->all(FLERR,"Illegal create_grid command");
    int newlevel = atoi(arg[iarg]);
    if (newlevel != level) 
      error->all(FLERR,"Create grid levels must be in ascending order");
    bounds(arg[iarg+1],nx,xlo,xhi);
    bounds(arg[iarg+2],ny,ylo,yhi);
    bounds(arg[iarg+3],nz,zlo,zhi);
    nx = atoi(arg[iarg+4]);
    ny = atoi(arg[iarg+5]);
    nz = atoi(arg[iarg+6]);
    if (nx < 1 || ny < 1 || nz < 1)
      error->all(FLERR,"Illegal create_grid command");
    if (domain->dimension == 2) {
      if (zlo != 1 || zhi != 1 || nz != 1) 
        error->all(FLERR,"Illegal create_grid command");
    }
    iarg += 7;
  }

  // invoke grid methods to complete grid setup

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

  if (nlo < 1 || nhi > nmax) error->all(FLERR,"Numeric index is out of bounds");
}
