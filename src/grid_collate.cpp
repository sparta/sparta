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
#include "grid.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

// ----------------------------------------------------------------------
// collate operations for implicit surfs, done as per-grid cell value
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   vector version of collate for implicit surf tallies into grid cells
   input tallies should be only for unsplit/split cells, not sub-cells
     caller can copy output to sub-cells if needed
   n = # of tallies for my owned+ghost cells
   ids = grid cell IDs for each tally (actually implicit surf IDs)
   in = value for each tally
   return out = summed values for my owned cells, including ghost contributions
     will only be values for unsplit and split cells, not sub-cells
   communication of ghost tallies done via rendezvous with irregular option
   called from fix ave/grid
------------------------------------------------------------------------- */

void Grid::collate_vector_implicit(int n, cellint *ids,
                                   double *in, double *out)
{
  int i,icell;
  cellint cellID;

  // if grid cell hash is not current, create it

  if (!hashfilled) rehash();

  // zero output values, only for owned cells

  if (nlocal) memset(out,0,nlocal*sizeof(double));

  // if I own grid cell, sum in value to out values directly
  // else nsend = # of tallies to contribute to irregular

  int nsend = 0;
  for (i = 0; i < n; i++) {
    icell = (*hash)[ids[i]];
    if (icell >= nlocal) nsend++;
    else out[icell] += in[i];
  }

  // done if just one proc

  if (comm->nprocs == 1) return;

  // create rvous inputs
  // nsend = # of datums to send

  int *proclist;
  memory->create(proclist,nsend,"grid:proclist");
  double *in_rvous;
  memory->create(in_rvous,2*nsend,"grid:in_rvous");

  int m = 0;
  nsend = 0;
  for (i = 0; i < n; i++) {
    icell = (*hash)[ids[i]];
    if (icell >= nlocal) {
      proclist[nsend] = cells[icell].proc;
      in_rvous[m++] = ubuf(ids[i]).d;
      in_rvous[m++] = in[i];
    }
  }

  // perform rendezvous operation
  // which = 0 for irregular comm, 1 for all2all comm
  // use irregular comm with clumped grid decomposition
  //   ghost tallies should need sending to a handful of nearby procs,
  //   even if cutoff is infinite
  // use all2all comm with dispersed grid decomposition
  //   ghost tallies could need sending to nearly all other procs

  int which = 0;
  if (!clumped) which = 1;
  
  char *buf;
  int nout = comm->rendezvous(which,nsend,(char *) in_rvous,
                              2*sizeof(double),
                              0,proclist,NULL,
                              0,buf,2*sizeof(double),(void *) this);
  double *out_rvous = (double *) buf;

  memory->destroy(proclist);
  memory->destroy(in_rvous);

  // sum tallies returned for grid cells I own into out

  m = 0;
  for (i = 0; i < nout; i++) {
    cellID = (cellint) ubuf(out_rvous[m++]).u;
    icell = (*hash)[cellID];
    out[icell] += out_rvous[m++];
  }
  
  // clean-up
  
  memory->destroy(out_rvous);
}

/* ----------------------------------------------------------------------
   array version of collate for implicit surf tallies into grid cells
   input tallies should be only for unsplit/split cells, not sub-cells
     caller can copy output to sub-cells if needed
   nrow = # of tallies for my owned+ghost cells
   ncol = # of values per tally
   ids = grid cell IDs for each tally (actually implicit surf IDs)
   in = ncol values for each tally
   return out = summed values for my owned cells, including ghost contributions
     will only be values for unsplit and split cells, not sub-cells
   communication of ghost tallies done via rendezvous with irregular option
   called from compute isurf/grid, compute react/isurf/grid, 
     fix ave/grid, surf react/implicit
------------------------------------------------------------------------- */

void Grid::collate_array_implicit(int nrow, int ncol, cellint *ids, 
                                  double **in, double **out)
{
  int i,j,icell;
  cellint cellID;

  // if grid cell hash is not current, create it

  if (!hashfilled) rehash();

  // zero output values for owned cells

  if (nlocal) memset(&out[0][0],0,nlocal*ncol*sizeof(double));

  // if I own grid cell, sum in values to out values directly
  // else nsend = # of tallies to contribute to rendezvous

  int nsend = 0;
  for (i = 0; i < nrow; i++) {
    icell = (*hash)[ids[i]];
    if (icell >= nlocal) nsend++;
    else {
      for (j = 0; j < ncol; j++)
        out[icell][j] += in[i][j];
    }
  }

  // done if just one proc

  if (comm->nprocs == 1) return;

  // create rvous inputs
  // nsend = # of datums to send

  int *proclist;
  memory->create(proclist,nsend,"grid:proclist");
  double *in_rvous;
  memory->create(in_rvous,(1+ncol)*nsend,"grid:in_rvous");

  int m = 0;
  nsend = 0;
  for (i = 0; i < nrow; i++) {
    icell = (*hash)[ids[i]];
    if (icell >= nlocal) {
      proclist[nsend] = cells[icell].proc;
      in_rvous[m++] = ubuf(ids[i]).d;
      for (j = 0; j < ncol; j++)
        in_rvous[m++] = in[i][j];
      nsend++;
    }
  }

  // perform rendezvous operation
  // which = 0 for irregular comm, 1 for all2all comm
  // use irregular comm with clumped grid decomposition
  //   ghost tallies should need sending to a handful of nearby procs,
  //   even if cutoff is infinite
  // use all2all comm with dispersed grid decomposition
  //   ghost tallies could need sending to nearly all other procs

  int which = 0;
  if (!clumped) which = 1;
  
  char *buf;
  int nout = comm->rendezvous(which,nsend,(char *) in_rvous,
                              (ncol+1)*sizeof(double),
                              0,proclist,NULL,
                              0,buf,(ncol+1)*sizeof(double),(void *) this);
  double *out_rvous = (double *) buf;

  memory->destroy(proclist);
  memory->destroy(in_rvous);

  // sum tallies returned for grid cells I own into out

  m = 0;
  for (i = 0; i < nout; i++) {
    cellID = (cellint) ubuf(out_rvous[m++]).u;
    icell = (*hash)[cellID];
    for (j = 0; j < ncol; j++)
      out[icell][j] += out_rvous[m++];
  }
  
  // clean-up
  
  memory->destroy(out_rvous);
}

// ----------------------------------------------------------------------
// operation to communicate grid data from owned cells to ghost cells
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   array version of owned->ghost comm with ncol values per grid cell
   requested ghost cells receive copy of data in corresponding owned cell
   nrequest = # of my ghost cells which need values
   ncol = # of values per cell
   ghostIDs = cell IDs of nrequest ghost cells
   in = ncol values for all owned cells
   return out = ncol values for all ghost cells (only ghostIDs are filled in)
   called from surf react/implicit, only if nprocs > 1
------------------------------------------------------------------------- */

void Grid::owned_to_ghost_array(int nrequest, int ncol, cellint *ghostIDs,
                                double **in, double **out)
{
  int i,j,icell,itally;
  cellint cellID;

  int me = comm->me;
  
  // if grid cell hash is not current, create it

  if (!hashfilled) rehash();
  
  // zero output values for ghost cells

  if (out) memset(&out[nlocal][0],0,nghost*ncol*sizeof(double));
  
  // send ghost cell IDs to owning procs as request for data
  // datum = sending proc, ghost cell ID
  
  int *proclist;
  memory->create(proclist,nrequest,"grid:proclist");
  double *in_rvous;
  memory->create(in_rvous,2*nrequest,"grid:in_rvous");

  int m = 0;
  int nsend = 0;
  for (i = 0; i < nrequest; i++) {
    icell = (*hash)[ghostIDs[i]];
    proclist[nsend] = cells[icell].proc;
    in_rvous[m++] = ubuf(me).d;
    in_rvous[m++] = ubuf(ghostIDs[i]).d;
    nsend++;
  }

  // perform rendezvous operation
  // which = 0 for irregular comm, 1 for all2all comm
  // use irregular comm with grid cutoff (requires clumped decomp)
  //   my owned cells only sent to handful of nearby procs
  // use all2all comm with infinite cutoff (clumped or dispersed grid)
  //   my owned cells will be sent to all other procs

  int which = 0;
  if (cutoff < 0.0) which = 1;

  ncol_rvous = ncol;
  owned_data_rvous = in;
  
  char *buf;
  int nout = comm->rendezvous(which,nsend,(char *) in_rvous,
                              2*sizeof(double),
                              0,proclist,rendezvous_owned_to_ghost,
                              0,buf,(ncol+1)*sizeof(double),(void *) this);
  double *out_rvous = (double *) buf;

  memory->destroy(proclist);
  memory->destroy(in_rvous);

  // copy returned out_rvous data to out for each of my ghost cells

  m = 0;
  for (i = 0; i < nout; i++) {
    cellID = (cellint) ubuf(out_rvous[m++]).u;
    icell = (*hash)[cellID];
    for (j = 0; j < ncol; j++)
      out[icell][j] += out_rvous[m++];
  }
  
  // clean-up
  
  memory->destroy(out_rvous);
}

/* ----------------------------------------------------------------------
   callback from owned_to_ghost_array_neighbors rendezvous operation
   process requests for ghost cell data
   inbuf = list of N Inbuf datums
   outbuf = cellID + Ncol values back to requesting proc
------------------------------------------------------------------------- */

int Grid::rendezvous_owned_to_ghost(int n, char *inbuf,
                                    int &flag, int *&proclist, char *&outbuf,
                                    void *ptr)
{
  int i,j,k,m,proc;
  cellint cellID;
  
  Grid *gptr = (Grid *) ptr;
  Memory *memory = gptr->memory;
  MyHash *hash = gptr->hash;
  int ncol = gptr->ncol_rvous;
  double **owned_data = gptr->owned_data_rvous;
  
  // allocate proclist & outbuf based on size of inbuf

  memory->create(proclist,n,"grid:proclist");
  double *out;
  memory->create(out,n*(ncol+1),"grid:out");

  // loop over (proc,cellID) datums in inbuf
  // copy data from corresponding owned cell ID into out
  
  double *in = (double *) inbuf;
  
  m = k = 0;
  for (i = 0; i < n; i++) {
    proc = (int) ubuf(in[m++]).i;
    cellID = (cellint) ubuf(in[m++]).u;
    int icell = (*hash)[cellID];
    
    proclist[i] = proc;
    out[k++] = ubuf(cellID).d;
    for (j = 0; j < ncol; j++)
      out[k++] = owned_data[icell][j];
  }

  // flag = 2: new outbuf

  flag = 2;
  outbuf = (char *) out;
  return n;
}
