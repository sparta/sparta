/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
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

#include <map>

using namespace SPARTA_NS;

// ----------------------------------------------------------------------
// collate operations for implicit surfs, done as per-grid cell value
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   helpers for collate of implicit surf tallies
   the tally IDs are grid cell IDs, but the cell may be neither owned nor
     ghost by the tallying proc anymore if the grid changed after the
     tallies were made (e.g. fix balance migrated cells on the tally step,
     or between steps of a fix ave/grid window)
   such tallies are routed to the cell's current owner via a rendezvous
     directory of owned cell IDs
   tallies for cells that no longer exist anywhere are dropped
------------------------------------------------------------------------- */

namespace {
  union ubufc {
    double d;
    int64_t i;
    ubufc(double arg) : d(arg) {}
    ubufc(int64_t arg) : i(arg) {}
  };

  struct RvousCollate {
    Memory *memory;
    int ncol;
  };

  inline int rvous_proc(cellint id, int nprocs)
  {
    uint64_t h = (uint64_t) id;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    return (int) (h % (uint64_t) nprocs);
  }

  // input datum = cell ID, owner proc (-1 for a tally), ncol values
  // output datum = cell ID, ncol values, sent to the owner proc

  int rvous_route(int n, char *inbuf, int &flag, int *&proclist,
                  char *&outbuf, void *ptr)
  {
    RvousCollate *rc = (RvousCollate *) ptr;
    int ncol = rc->ncol;
    int insize = ncol + 2;
    double *in = (double *) inbuf;

    std::map<cellint,int> owner;
    int ntally = 0;
    for (int i = 0; i < n; i++) {
      double *one = &in[(bigint) i*insize];
      int proc = (int) ubufc(one[1]).i;
      if (proc >= 0) owner[(cellint) ubufc(one[0]).i] = proc;
      else ntally++;
    }

    rc->memory->create(proclist,ntally,"grid:proclist");
    double *out = (double *)
      rc->memory->smalloc((bigint) ntally*(ncol+1)*sizeof(double),
                          "grid:outbuf");

    int nout = 0;
    bigint m = 0;
    for (int i = 0; i < n; i++) {
      double *one = &in[(bigint) i*insize];
      if ((int) ubufc(one[1]).i >= 0) continue;
      std::map<cellint,int>::iterator it =
        owner.find((cellint) ubufc(one[0]).i);
      if (it == owner.end()) continue;
      proclist[nout++] = it->second;
      out[m++] = one[0];
      for (int j = 0; j < ncol; j++) out[m++] = one[2+j];
    }

    flag = 2;
    outbuf = (char *) out;
    return nout;
  }
}

/* ----------------------------------------------------------------------
   common collate for vector (ncol = 1) and array versions
   in[i] = pointer to ncol values of tally i
   out[icell] = pointer to ncol values of owned cell icell, zeroed here
   fast path (every tally ID is an owned or ghost cell of its proc):
     identical to the original algorithm, plus one MPI_Allreduce
------------------------------------------------------------------------- */

void Grid::collate_implicit(int nrow, int ncol, cellint *ids,
                            double **in, double **out)
{
  int i,j,icell;
  cellint cellID;

  // if grid cell hash is not current, create it

  if (!hashfilled) rehash();

  // zero output values for owned cells

  for (icell = 0; icell < nlocal; icell++)
    memset(out[icell],0,ncol*sizeof(double));

  // look up each tally ID, never insert into hash
  // icellrow[i] = owned/ghost index, -1 if neither (unknown)
  // if I own grid cell, sum in values to out values directly
  // else nsend = # of tallies to contribute to rendezvous

  int *icellrow;
  memory->create(icellrow,MAX(nrow,1),"grid:icellrow");

  int nsend = 0;
  int nunknown = 0;
  for (i = 0; i < nrow; i++) {
    MyHash::iterator it = hash->find(ids[i]);
    if (it == hash->end()) {
      icellrow[i] = -1;
      nunknown++;
      continue;
    }
    icell = icellrow[i] = it->second;
    if (icell >= nlocal) nsend++;
    else {
      for (j = 0; j < ncol; j++)
        out[icell][j] += in[i][j];
    }
  }

  // done if just one proc
  // unknown tallies on one proc are for cells that no longer exist

  if (comm->nprocs == 1) {
    memory->destroy(icellrow);
    return;
  }

  int anyunknown;
  MPI_Allreduce(&nunknown,&anyunknown,1,MPI_INT,MPI_MAX,world);

  // create rvous inputs
  // nsend = # of datums to send

  int *proclist;
  memory->create(proclist,nsend,"grid:proclist");
  double *in_rvous;
  if ((bigint) (1+ncol)*nsend > MAXSMALLINT)
    error->one(FLERR,"Grid collate buffer exceeds 2 GB");
  memory->create(in_rvous,(1+ncol)*nsend,"grid:in_rvous");

  bigint m = 0;
  nsend = 0;
  for (i = 0; i < nrow; i++) {
    icell = icellrow[i];
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

  memory->destroy(out_rvous);

  if (!anyunknown) {
    memory->destroy(icellrow);
    return;
  }

  // route unknown tallies to the current owners of their cells
  // rendezvous input: directory of my owned cells (not sub cells)
  //   + my unknown tallies, both hashed by cell ID to a rendezvous proc

  int nprocs = comm->nprocs;
  int me = comm->me;

  int ndir = 0;
  for (icell = 0; icell < nlocal; icell++)
    if (cells[icell].nsplit > 0) ndir++;

  int nsend2 = ndir + nunknown;
  int insize = ncol + 2;
  if ((bigint) insize*nsend2 > MAXSMALLINT)
    error->one(FLERR,"Grid collate buffer exceeds 2 GB");
  double *in2;
  memory->create(proclist,MAX(nsend2,1),"grid:proclist");
  memory->create(in2,MAX(nsend2,1)*insize,"grid:in_rvous");

  int n = 0;
  m = 0;
  for (icell = 0; icell < nlocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    proclist[n++] = rvous_proc(cells[icell].id,nprocs);
    in2[m++] = ubufc((int64_t) cells[icell].id).d;
    in2[m++] = ubufc((int64_t) me).d;
    for (j = 0; j < ncol; j++) in2[m++] = 0.0;
  }

  for (i = 0; i < nrow; i++) {
    if (icellrow[i] >= 0) continue;
    proclist[n++] = rvous_proc(ids[i],nprocs);
    in2[m++] = ubufc((int64_t) ids[i]).d;
    in2[m++] = ubufc((int64_t) -1).d;
    for (j = 0; j < ncol; j++) in2[m++] = in[i][j];
  }

  memory->destroy(icellrow);

  RvousCollate rc;
  rc.memory = memory;
  rc.ncol = ncol;

  nout = comm->rendezvous(1,nsend2,(char *) in2,insize*sizeof(double),
                          0,proclist,rvous_route,
                          0,buf,(ncol+1)*sizeof(double),(void *) &rc);
  out_rvous = (double *) buf;

  memory->destroy(proclist);
  memory->destroy(in2);

  // sum routed tallies into my owned cells

  m = 0;
  for (i = 0; i < nout; i++) {
    cellID = (cellint) ubufc(out_rvous[m++]).i;
    MyHash::iterator it = hash->find(cellID);
    if (it == hash->end() || it->second >= nlocal) {
      m += ncol;
      continue;
    }
    icell = it->second;
    for (j = 0; j < ncol; j++)
      out[icell][j] += out_rvous[m++];
  }

  memory->sfree(out_rvous);
}

/* ----------------------------------------------------------------------
   vector version of collate for implicit surf tallies into grid cells
   input tallies should be only for unsplit/split cells, not sub-cells
     caller can copy output to sub-cells if needed
   n = # of tallies for my owned+ghost cells (when they were made)
   ids = grid cell IDs for each tally (actually implicit surf IDs)
   in = value for each tally
   return out = summed values for my owned cells, including ghost contributions
     and tallies whose cell has migrated since they were made
     will only be values for unsplit and split cells, not sub-cells
   communication of ghost tallies done via rendezvous with irregular option
   called from fix ave/grid
------------------------------------------------------------------------- */

void Grid::collate_vector_implicit(int n, cellint *ids,
                                   double *in, double *out)
{
  double **inrow = (double **)
    memory->smalloc((bigint) MAX(n,1)*sizeof(double *),"grid:inrow");
  double **outrow = (double **)
    memory->smalloc((bigint) MAX(nlocal,1)*sizeof(double *),"grid:outrow");
  for (int i = 0; i < n; i++) inrow[i] = &in[i];
  for (int i = 0; i < nlocal; i++) outrow[i] = &out[i];

  collate_implicit(n,1,ids,inrow,outrow);

  memory->sfree(inrow);
  memory->sfree(outrow);
}

/* ----------------------------------------------------------------------
   array version of collate for implicit surf tallies into grid cells
   input tallies should be only for unsplit/split cells, not sub-cells
     caller can copy output to sub-cells if needed
   nrow = # of tallies for my owned+ghost cells (when they were made)
   ncol = # of values per tally
   ids = grid cell IDs for each tally (actually implicit surf IDs)
   in = ncol values for each tally
   return out = summed values for my owned cells, including ghost contributions
     and tallies whose cell has migrated since they were made
     will only be values for unsplit and split cells, not sub-cells
   communication of ghost tallies done via rendezvous with irregular option
   called from compute isurf/grid, compute react/isurf/grid,
     fix ave/grid, surf react/implicit
------------------------------------------------------------------------- */

void Grid::collate_array_implicit(int nrow, int ncol, cellint *ids,
                                  double **in, double **out)
{
  collate_implicit(nrow,ncol,ids,in,out);
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
  int i,j,icell;
  cellint cellID;

  int me = comm->me;

  // if grid cell hash is not current, create it

  if (!hashfilled) rehash();

  // zero output values for ghost cells

  if (out) memset(&out[nlocal][0],0,(bigint) nghost*ncol*sizeof(double));

  // send ghost cell IDs to owning procs as request for data
  // datum = sending proc, ghost cell ID

  int *proclist;
  memory->create(proclist,nrequest,"grid:proclist");
  double *in_rvous;
  memory->create(in_rvous,2*nrequest,"grid:in_rvous");

  bigint m = 0;
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
  int i,j,proc;
  bigint k,m;
  cellint cellID;

  Grid *gptr = (Grid *) ptr;
  Memory *memory = gptr->memory;
  MyHash *hash = gptr->hash;
  int ncol = gptr->ncol_rvous;
  double **owned_data = gptr->owned_data_rvous;

  // allocate proclist & outbuf based on size of inbuf

  memory->create(proclist,n,"grid:proclist");
  double *out;
  if ((bigint) n*(ncol+1) > MAXSMALLINT)
    gptr->error->one(FLERR,"Grid collate buffer exceeds 2 GB");
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
