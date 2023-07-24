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
#include "remove_surf.h"
#include "surf.h"
#include "grid.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

RemoveSurf::RemoveSurf(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

RemoveSurf::~RemoveSurf() {}

/* ---------------------------------------------------------------------- */

void RemoveSurf::command(int narg, char **arg)
{
  if (!surf->exist)
    error->all(FLERR,"Cannot remove_surf before surf elements are defined");
  if (surf->implicit)
    error->all(FLERR,"Cannot remove_surf for implicit surfs");

  if (narg < 1) error->all(FLERR,"Illegal remove_surf command");

  int igroup = surf->find_group(arg[0]);
  if (igroup < 0) error->all(FLERR,"Remove surf group ID does not exist");
  int groupbit = surf->bitmask[igroup];

  if (comm->me == 0)
    if (screen) fprintf(screen,"Removing surfs ...\n");

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // sort particles

  if (particle->exist) particle->sort();

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // remove all surfs in group

  nsurf_old = surf->nsurf;
  nremove = remove(groupbit);
  nsurf_new = nsurf_old - nremove;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  removed " BIGINT_FORMAT " surfs\n",nremove);
      fprintf(screen,"  " BIGINT_FORMAT " surfs remain\n",nsurf_new);
    }
    if (logfile) {
      fprintf(logfile,"  removed " BIGINT_FORMAT " surfs\n",nremove);
      fprintf(logfile,"  " BIGINT_FORMAT " surfsremain\n",nsurf_new);
    }
  }

  // replace surfs in Surf with reduced set of local surfs and custom values
  // ncustom = # of current custom vecs/arrays
  // index_custom = indices for each custom vec/array in Surf custom list
  
  int ncustom = 0;
  int *index_custom = NULL;

  if (surf->ncustom) {
    for (int i = 0; i < surf->ncustom; i++) {
      if (!surf->ename[i]) continue;
      index_custom[ncustom++] = i;
    }
  }
  
  if (nremove) surf->add_surfs(1,nsurf,lines,tris,
			       ncustom,index_custom,cvalues);

  memory->sfree(lines);
  memory->sfree(tris);
  memory->destroy(cvalues);
  delete [] index_custom;
  
  // check that remaining surfs are still watertight

  if (domain->dimension == 2) surf->check_watertight_2d();
  else surf->check_watertight_3d();
  
  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // reset grid due to changing surfs
  // assign surfs to grid cells
  
  surf->setup_owned();
  grid->unset_neighbors();
  grid->remove_ghosts();

  // reassign split cell particles to parent split cell

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->combine_split_cell_particles(icell,1);
  }

  grid->clear_surf();
  if (surf->exist) grid->surf2grid(1);

  // reassign particles in split cells to sub cell owner

  if (particle->exist && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->assign_split_cell_particles(icell);
  }

  MPI_Barrier(world);
  double time4 = MPI_Wtime();

  // re-setup owned and ghost cell info

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  grid->set_inout();
  grid->type_check();

  MPI_Barrier(world);
  double time5 = MPI_Wtime();

  double time_total = time3-time1;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  sort/remove/surf2grid/ghost percent = %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total);
    }
    if (logfile) {
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  sort/remove/surf2grid/ghost percent = %g %g %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total,
              100.0*(time4-time3)/time_total,100.0*(time5-time4)/time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   copy all surfs and custom data into local data structs
   remove all surfs in surf group from local data structs
------------------------------------------------------------------------- */

bigint RemoveSurf::remove(int groupbit)
{
  int i,m,n;

  int dim = domain->dimension;
  int distributed = surf->distributed;

  int nbytes;
  if (dim == 2) nbytes = sizeof(Surf::Line);
  else nbytes = sizeof(Surf::Tri);
  
  // copy Surf::lines/tris to local lines/tris
  // local data is distributed via striding

  lines = NULL;
  tris = NULL;
  nsurf = surf->nown;

  int me = comm->me;
  int nprocs = comm->nprocs;
  
  if (dim == 2) {
    lines = (Surf::Line *) memory->smalloc(nsurf*nbytes,"remove/surf:lines");
    if (distributed) memcpy(lines,surf->mylines,nsurf*nbytes);
    else {
      int nslocal = surf->nlocal;
      m = 0;
      for (int i = me; i < nslocal; i += nprocs)
	memcpy(&lines[m],&surf->lines[i],nbytes);
    }
  } else {
    tris = (Surf::Tri *) memory->smalloc(nsurf*nbytes,"remove/surf:tris");
    if (distributed) memcpy(tris,surf->mytris,nsurf*nbytes);
    else {
      int nslocal = surf->nlocal;
      m = 0;
      for (int i = me; i < nslocal; i += nprocs)
	memcpy(&tris[m],&surf->tris[i],nbytes);
    }
  }
  
  // copy Surf custom data to local custom data

  cvalues = NULL;
  int ncustom = surf->ncustom;
  int ncbytes = 0;

  if (ncustom) {
    int nvalues_custom = surf->extract_custom(cvalues);
    ncbytes = (1+nvalues_custom) * sizeof(double);
  }
  
  // remove surfs in group, both from lines/tris and cvalues

  n = 0;
  for (i = 0; i < nsurf; i++) {
    if (dim == 2) {
      if (!(lines[i].mask & groupbit)) continue;
      if (i != n) memcpy(&lines[n],&lines[i],nbytes);
    } else {
      if (!(tris[i].mask & groupbit)) continue;
      if (i != n) memcpy(&tris[n],&tris[i],nbytes);
    }
    if (ncustom) memcpy(&cvalues[n],&cvalues[i],ncbytes);
    n++;
  }

  bigint ndiscard_me = nsurf - n;
  nsurf = n;

  // ndiscard = total # of removed surfs

  bigint ndiscard;
  MPI_Allreduce(&ndiscard_me,&ndiscard,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // if removed any surfs, renumber surf IDs across all procs

  if (ndiscard) {
    bigint bnsurf = nsurf;
    bigint offset;
    MPI_Scan(&bnsurf,&offset,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    offset -= bnsurf;

    if (dim == 2)
      for (i = 0; i < nsurf; i++)
	lines[i].id = static_cast<surfint> (offset + i + 1);
    else
      for (i = 0; i < nsurf; i++)
	tris[i].id = static_cast<surfint> (offset + i + 1);

    surfint id;
    if (ncustom)
      for (i = 0; i < nsurf; i++) {
	id = static_cast<surfint> (offset + i + 1);
	cvalues[i][0] = ubuf(id).d;
      }
  }

  // clean up from call to surf->extract_custom()
  
  if (ncustom) memory->destroy(cvalues);
  
  return ndiscard;
}
