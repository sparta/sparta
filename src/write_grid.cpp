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

#include "spatype.h"
#include "string.h"
#include "write_grid.h"
#include "grid.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#ifdef SPARTA_MAP
#include <map>
#else
#include <tr1/unordered_map>
#endif

using namespace SPARTA_NS;

#define MAXLINE 256

/* ---------------------------------------------------------------------- */

WriteGrid::WriteGrid(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void WriteGrid::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot write grid when grid is not defined");

  if (narg != 1) error->all(FLERR,"Illegal write_grid command");

  // write file, create parent cells and then child cells

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  int me = comm->me;
  if (me == 0) {
    if (screen) fprintf(screen,"Writing grid file ...\n");
    fp = fopen(arg[0],"w");
    if (!fp) {
      char str[128];
      sprintf(str,"Cannot open file %s",arg[0]);
      error->one(FLERR,str);
    }
  }

  if (me == 0) header();

  // write Parents section

  if (me == 0) write_parents();

  // close file

  if (me == 0) fclose(fp);

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // stats

  double time_total = time2-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  parent cells = %d\n",grid->nparent);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
    }

    if (logfile) {
      fprintf(logfile,"  parent cells = %d\n",grid->nparent);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   write header of grid file
   only called by proc 0
------------------------------------------------------------------------- */

void WriteGrid::header()
{
  fprintf(fp,"# grid file written by SPARTA\n\n");
  fprintf(fp,"%d nparents\n",grid->nparent);
}

/* ----------------------------------------------------------------------
   write Parents section of grid file
   only called by proc 0
------------------------------------------------------------------------- */

void WriteGrid::write_parents()
{
  char str[32];

  // fill hash with parent IDs if necessary

  if (!grid->hashfilled) {

#ifdef SPARTA_MAP
    std::map<cellint,int> *hash = grid->hash;
#else
    std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

    hash->clear();

    Grid::ParentCell *pcells = grid->pcells;
    int nparent = grid->nparent;

    for (int icell = 0; icell < nparent; icell++)
      (*hash)[pcells[icell].id] = -(icell+1);
  }

  fprintf(fp,"\nParents\n\n");

  Grid::ParentCell *pcells = grid->pcells;
  int nparent = grid->nparent;

  // one parent cell per line

  for (int i = 0; i < nparent; i++) {
    grid->id_num2str(pcells[i].id,str);
    fprintf(fp,"%d %s %d %d %d\n",i+1,str,
            pcells[i].nx,pcells[i].ny,pcells[i].nz);
  }

  // clear hash if filled it

  if (!grid->hashfilled) grid->hash->clear();
}
