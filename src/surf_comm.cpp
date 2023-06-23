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

#include "surf.h"
#include "domain.h"
#include "comm.h"
#include "grid.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{TALLYAUTO,TALLYREDUCE,TALLYRVOUS};         // same as Update
enum{INT,DOUBLE};                      // several files

// ----------------------------------------------------------------------
// compress operations after load-balancing or grid adaptation
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   compress owned explicit distributed surfs to account for deleted grid cells
     either due to load-balancing migration or grid adapt coarsening
   called from Comm::migrate_cells() and AdaptGrid::coarsen()
     AFTER grid cells are compressed
   discard nlocal surfs that are no longer referenced by owned grid cells
   use hash to store referenced surfs
   only called for explicit distributed surfs
------------------------------------------------------------------------- */

void Surf::compress_explicit()
{
  int i,m,ns;
  surfint *csurfs;

  int dim = domain->dimension;

  // keep = 1 if a local surf is referenced by a compressed local grid cell

  int *keep;
  memory->create(keep,nlocal,"surf:keep");
  for (i = 0; i < nlocal; i++) keep[i] = 0;

  // convert grid cell csurfs to surf IDs so can reset after surf compression
  // skip cells with no surfs or sub-cells

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  for (i = 0; i < nglocal; i++) {
    if (!cells[i].nsurf) continue;
    if (cells[i].nsplit <= 0) continue;
    csurfs = cells[i].csurfs;
    ns = cells[i].nsurf;
    if (dim == 2) {
      for (m = 0; m < ns; m++) {
        keep[csurfs[m]] = 1;
        csurfs[m] = lines[csurfs[m]].id;
      }
    } else {
      for (m = 0; m < ns; m++) {
        keep[csurfs[m]] = 1;
        csurfs[m] = tris[csurfs[m]].id;
      }
    }
  }

  // compress nlocal surfs based on keep flags

  m = 0;
  while (i < nlocal) {
    if (!keep[i]) {
      if (dim == 2) memcpy(&lines[i],&lines[nlocal-1],sizeof(Line));
      else memcpy(&tris[i],&tris[nlocal-1],sizeof(Tri));
      keep[i] = keep[nlocal-1];
      nlocal--;
    } else i++;
  }

  memory->destroy(keep);

  // reset grid cell csurfs IDs back to local surf indices
  // hash compressed surf list, then clear hash
  // skip cells with no surfs or sub-cells

  rehash();

  for (i = 0; i < nglocal; i++) {
    if (!cells[i].nsurf) continue;
    if (cells[i].nsplit <= 0) continue;
    csurfs = cells[i].csurfs;
    ns = cells[i].nsurf;
    for (m = 0; m < ns; m++) csurfs[m] = (*hash)[csurfs[m]];
  }

  hash->clear();
  hashfilled = 0;
}

/* ----------------------------------------------------------------------
   compress owned implicit surfs to account for migrating grid cells
   called from Comm::migrate_cells() BEFORE grid cells are compressed
   migrating grid cells are ones with proc != me
   reset csurfs indices for kept cells
   only called for implicit surfs
------------------------------------------------------------------------- */

void Surf::compress_implicit()
{
  int j,ns,icell;
  cellint cellID;
  surfint *csurfs;

  if (!grid->hashfilled) grid->rehash();

  Grid::ChildCell *cells = grid->cells;
  Grid::MyHash *ghash = grid->hash;
  int me = comm->me;
  int n = 0;

  if (domain->dimension == 2) {
    for (int i = 0; i < nlocal; i++) {
      icell = (*ghash)[lines[i].id];
      if (cells[icell].proc != me) continue;
      if (i != n) {
        // compress my surf list
        memcpy(&lines[n],&lines[i],sizeof(Line));
        // reset matching csurfs index in grid cell from i to n
        csurfs = cells[icell].csurfs;
        ns = cells[icell].nsurf;
        for (j = 0; j < ns; j++)
          if (csurfs[j] == i) {
            csurfs[j] = n;
            break;
          }
      }
      n++;
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      icell = (*ghash)[tris[i].id];
      if (cells[icell].proc != me) continue;
      if (i != n) {
        // compress my surf list
        memcpy(&tris[n],&tris[i],sizeof(Tri));
        // reset matching csurfs index in grid cell from i to n
        csurfs = cells[icell].csurfs;
        ns = cells[icell].nsurf;
        for (j = 0; j < ns; j++)
          if (csurfs[j] == i) {
            csurfs[j] = n;
            break;
          }
      }
      n++;
    }
  }

  nlocal = n;
}

// ----------------------------------------------------------------------
// collate operations for summing tallies to owned surfs
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   comm of tallies across all procs for a single value (vector)
   nrow = # of tally entries in input vector
   tally2surf = surf index of each entry in input vector
   in = input vector of tallies
   instride = stride between entries in input vector
   return out = summed tallies for explicit surfs I own
   only called for explicit surfs
------------------------------------------------------------------------- */

void Surf::collate_vector(int nrow, surfint *tally2surf,
                          double *in, int instride, double *out)
{
  // collate algorithm depends on tally_comm setting

  if (tally_comm == TALLYAUTO) {
    if (nprocs > nsurf)
      collate_vector_reduce(nrow,tally2surf,in,instride,out);
    else collate_vector_rendezvous(nrow,tally2surf,in,instride,out);
  } else if (tally_comm == TALLYREDUCE) {
    collate_vector_reduce(nrow,tally2surf,in,instride,out);
  } else if (tally_comm == TALLYRVOUS) {
    collate_vector_rendezvous(nrow,tally2surf,in,instride,out);
  }
}

/* ----------------------------------------------------------------------
   allreduce algorithm for collate_vector
------------------------------------------------------------------------- */

void Surf::collate_vector_reduce(int nrow, surfint *tally2surf,
                                 double *in, int instride, double *out)
{
  int i,j,m;

  if (nsurf > MAXSMALLINT)
    error->all(FLERR,"Two many surfs to tally reduce - "
               "use global surf/comm auto or rvous");

  int nglobal = nsurf;

  double *one,*all;
  memory->create(one,nglobal,"surf:one");
  memory->create(all,nglobal,"surf:all");

  // zero all values and add in values I accumulated

  memset(one,0,nglobal*sizeof(double));

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int dim = domain->dimension;
  surfint id;

  j = 0;
  for (i = 0; i < nrow; i++) {
    m = (int) tally2surf[i] - 1;
    one[m] = in[j];
    j += instride;
  }

  // global allreduce

  MPI_Allreduce(one,all,nglobal,MPI_DOUBLE,MPI_SUM,world);

  // out = only surfs I own

  m = 0;
  for (i = me; i < nglobal; i += nprocs)
    out[m++] = all[i];

  // NOTE: could persist these for multiple invocations

  memory->destroy(one);
  memory->destroy(all);
}

/* ----------------------------------------------------------------------
   rendezvous algorithm for collate_vector
------------------------------------------------------------------------- */

void Surf::collate_vector_rendezvous(int nrow, surfint *tally2surf,
                                     double *in, int instride, double *out)
{
  // allocate memory for rvous input

  int *proclist;
  memory->create(proclist,nrow,"surf:proclist");
  InRvousVec *in_rvous =
    (InRvousVec *) memory->smalloc((bigint) nrow*sizeof(InRvousVec),
                                   "surf:in_rvous");

  // create rvous inputs
  // proclist = owner of each surf
  // logic of (id-1) % nprocs sends
  //   surf IDs 1,11,21,etc on 10 procs to proc 0

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int dim = domain->dimension;

  surfint id;

  int m = 0;
  for (int i = 0; i < nrow; i++) {
    id = tally2surf[i];
    proclist[i] = (id-1) % nprocs;
    in_rvous[i].id = id;
    in_rvous[i].value = in[m];
    m += instride;
  }

  // perform rendezvous operation
  // each proc owns subset of surfs
  // receives all tally contributions to surfs it owns

  out_rvous = out;

  char *buf;
  int nout = comm->rendezvous(1,nrow,(char *) in_rvous,sizeof(InRvousVec),
                              0,proclist,rendezvous_vector,
                              0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->destroy(in_rvous);
}

/* ----------------------------------------------------------------------
   callback from collate_rendezvous_vector rendezvous operatiion
   process tallies for surfs assigned to me
   inbuf = list of N Inbuf datums
   no outbuf
------------------------------------------------------------------------- */

int Surf::rendezvous_vector(int n, char *inbuf, int &flag, int *&proclist,
                            char *&outbuf, void *ptr)
{
  Surf *sptr = (Surf *) ptr;
  Memory *memory = sptr->memory;
  int nown = sptr->nown;
  double *out = sptr->out_rvous;
  int nprocs = sptr->comm->nprocs;
  int me = sptr->comm->me;

  // zero my owned surf values

  for (int i = 0; i < nown; i++) out[i] = 0.0;

  // accumulate per-surf values from different procs to my owned surfs
  // logic of (id-1-me) / nprocs maps
  //   surf IDs [1,11,21,...] on 10 procs to [0,1,2,...] on proc 0

  Surf::InRvousVec *in_rvous = (Surf::InRvousVec *) inbuf;

  int m;
  for (int i = 0; i < n; i++) {
    m = (in_rvous[i].id-1-me) / nprocs;
    out[m] += in_rvous[i].value;
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   comm of tallies across all procs for multiple values (array)
   nrow,ncol = # of entries and columns in input array
   tally2surf = global surf index of each entry in input array
   in = input array of tallies
   return out = summed tallies for explicit surfs I own
   only called for explicit surfs
------------------------------------------------------------------------- */

void Surf::collate_array(int nrow, int ncol, surfint *tally2surf,
                         double **in, double **out)
{
  // collate algorithm depends on tally_comm setting

  if (tally_comm == TALLYAUTO) {
    if (nprocs > nsurf)
      collate_array_reduce(nrow,ncol,tally2surf,in,out);
    else collate_array_rendezvous(nrow,ncol,tally2surf,in,out);
  } else if (tally_comm == TALLYREDUCE) {
    collate_array_reduce(nrow,ncol,tally2surf,in,out);
  } else if (tally_comm == TALLYRVOUS) {
    collate_array_rendezvous(nrow,ncol,tally2surf,in,out);
  }
}

/* ----------------------------------------------------------------------
   allreduce algorithm for collate_array()
------------------------------------------------------------------------- */

void Surf::collate_array_reduce(int nrow, int ncol, surfint *tally2surf,
                                double **in, double **out)
{
  int i,j,m;

  bigint ntotal = (bigint) nsurf * ncol;

  if (ntotal > MAXSMALLINT)
    error->all(FLERR,"Two many surfs to tally reduce - "
               "use global surf/comm auto or rvous");

  int nglobal = nsurf;

  double **one,**all;
  memory->create(one,nglobal,ncol,"surf:one");
  memory->create(all,nglobal,ncol,"surf:all");

  // zero all values and set values I accumulated

  memset(&one[0][0],0,nglobal*ncol*sizeof(double));

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int dim = domain->dimension;

  for (i = 0; i < nrow; i++) {
    m = (int) tally2surf[i] - 1;
    for (j = 0; j < ncol; j++)
      one[m][j] = in[i][j];
  }

  // global allreduce

  MPI_Allreduce(&one[0][0],&all[0][0],ntotal,MPI_DOUBLE,MPI_SUM,world);

  // out = only surfs I own

  m = 0;
  for (i = me; i < nglobal; i += nprocs) {
    for (j = 0; j < ncol; j++) out[m][j] = all[i][j];
    m++;
  }

  // NOTE: could persist these for multiple invocations

  memory->destroy(one);
  memory->destroy(all);
}

/* ----------------------------------------------------------------------
   rendezvous algorithm for collate_array
------------------------------------------------------------------------- */

void Surf::collate_array_rendezvous(int nrow, int ncol, surfint *tally2surf,
                                    double **in, double **out)
{
  int i,j,m;

  // allocate memory for rvous input

  int *proclist;
  memory->create(proclist,nrow,"surf:proclist");
  double *in_rvous = (double *)     // worry about overflow
    memory->smalloc(nrow*(ncol+1)*sizeof(double*),"surf:in_rvous");

  // create rvous inputs
  // proclist = owner of each surf
  // logic of (id-1) % nprocs sends
  //   surf IDs 1,11,21,etc on 10 procs to proc 0

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int dim = domain->dimension;
  surfint id;

  m = 0;
  for (int i = 0; i < nrow; i++) {
    id = tally2surf[i];
    proclist[i] = (id-1) % nprocs;
    in_rvous[m++] = ubuf(id).d;
    for (j = 0; j < ncol; j++)
      in_rvous[m++] = in[i][j];
  }

  // perform rendezvous operation
  // each proc owns subset of surfs
  // receives all tally contributions to surfs it owns

  ncol_rvous = ncol;
  if (out == NULL) out_rvous = NULL;
  else out_rvous = &out[0][0];
  int size = (ncol+1) * sizeof(double);

  char *buf;
  int nout = comm->rendezvous(1,nrow,(char *) in_rvous,size,
                              0,proclist,rendezvous_array,
                              0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->sfree(in_rvous);
}

/* ----------------------------------------------------------------------
   callback from collate_rendezvous_array rendezvous operation
   process tallies for surfs assigned to me
   inbuf = list of N Inbuf datums
   no outbuf
------------------------------------------------------------------------- */

int Surf::rendezvous_array(int n, char *inbuf,
                           int &flag, int *&proclist, char *&outbuf,
                           void *ptr)
{
  int i,j,k,m;

  Surf *sptr = (Surf *) ptr;
  Memory *memory = sptr->memory;
  int nown = sptr->nown;
  int ncol = sptr->ncol_rvous;
  double *out = sptr->out_rvous;
  int nprocs = sptr->comm->nprocs;
  int me = sptr->comm->me;

  // zero my owned surf values
  // NOTE: is this needed if caller zeroes ?

  int ntotal = nown*ncol;
  for (m = 0; m < ntotal; m++) out[m] = 0.0;

  // accumulate per-surf values from different procs to my owned surfs
  // logic of (id-1-me) / nprocs maps
  //   surf IDs [1,11,21,...] on 10 procs to [0,1,2,...] on proc 0

  double *in_rvous = (double *) inbuf;
  surfint id;

  m = 0;
  for (int i = 0; i < n; i++) {
    id = (surfint) ubuf(in_rvous[m++]).i;
    k = (id-1-me) / nprocs * ncol;
    for (j = 0; j < ncol; j++)
      out[k++] += in_rvous[m++];
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   comm of tallies across all procs for a single value (vector)
   called from compute isurf/grid and fix ave/grid
     for implicit surf tallies by grid cell
   nrow = # of tallies
   tally2surf = surf ID for each tally (same as cell ID)
   in = vectir of tally values
   return out = summed tallies for grid cells I own
   always done via rendezvous algorithm (no allreduce algorithm)
   only called for implicit surfs
------------------------------------------------------------------------- */

void Surf::collate_vector_implicit(int nrow, surfint *tally2surf,
                                   double *in, double *out)
{
  int i,j,m,icell;
  cellint cellID;

  int me = comm->me;
  int nprocs = comm->nprocs;

  // create a grid cell hash for only my owned cells

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  MyCellHash hash;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    hash[cells[icell].id] = icell;
  }

  // for implicit surfs, tally2surf stores cellIDs

  cellint *tally2cell = (cellint *) tally2surf;

  // if I own tally grid cell, sum tallies to out directly
  // else nsend = # of tallies to contribute to rendezvous

  int nsend = 0;
  for (i = 0; i < nrow; i++) {
    if (hash.find(tally2cell[i]) == hash.end()) nsend++;
    else {
      icell = hash[tally2cell[i]];
      out[icell] += in[i];
    }
  }

  // done if just one proc

  if (nprocs == 1) return;

  // ncell = # of owned grid cells with implicit surfs, excluding sub cells
  // NOTE: could limit to cell group of caller

  int ncell = 0;
  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsurf <= 0) continue;
    if (cells[icell].nsplit <= 0) continue;
    ncell++;
  }

  // allocate memory for rvous input
  // ncount = ncell + nsend
  // 3 doubles for each input = proc, cellID, tally

  int ncount = ncell + nsend;

  int *proclist;
  double *in_rvous;
  memory->create(proclist,ncount,"surf:proclist");
  memory->create(in_rvous,3*ncount,"surf:in_rvous");

  // create rvous inputs
  // owning proc for each datum = random hash of cellID
  // flavor 1: one per ncell with proc and cellID, no tally
  // flavor 2: one per nsend with proc = -1, cellID, one tally

  ncount = m = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsurf <= 0) continue;
    if (cells[icell].nsplit <= 0) continue;
    proclist[ncount] = hashlittle(&cells[icell].id,sizeof(cellint),0) % nprocs;
    in_rvous[m++] = me;
    in_rvous[m++] = cells[icell].id;    // NOTE: should use ubuf
    in_rvous[m++] = 0.0;
    ncount++;
  }

  for (i = 0; i < nrow; i++) {
    if (hash.find(tally2cell[i]) == hash.end()) {
      proclist[ncount] = hashlittle(&tally2cell[i],sizeof(cellint),0) % nprocs;
      in_rvous[m++] = -1;
      in_rvous[m++] = tally2cell[i];    // NOTE: should use ubuf
      in_rvous[m++] = in[i];
      ncount++;
    }
  }

  // perform rendezvous operation

  ncol_rvous = 1;
  char *buf;
  int nout = comm->rendezvous(1,ncount,(char *) in_rvous,3*sizeof(double),
                              0,proclist,rendezvous_implicit,
                              0,buf,2*sizeof(double),(void *) this);
  double *out_rvous = (double *) buf;

  memory->destroy(proclist);
  memory->destroy(in_rvous);

  // sum tallies returned for grid cells I own into out

  m = 0;
  for (i = 0; i < nout; i++) {
    cellID = out_rvous[m++];      // NOTE: should use ubuf
    icell = hash[cellID];
    out[icell] += out_rvous[m++];
  }

  // clean-up

  memory->destroy(out_rvous);
}

/* ----------------------------------------------------------------------
   comm of tallies across all procs for multiple values (array)
   called from compute isurf/grid and fix ave/grid
     for implicit surf tallies by grid cell
   nrow = # of tallies
   ncol = # of values per tally
   tally2surf = surf ID for each tally (same as cell ID)
   in = array of tally values, nrow by ncol
   return out = summed tallies for grid cells I own, nlocal by ncol
   always done via rendezvous algorithm (no allreduce algorithm)
   only called for implicit surfs
------------------------------------------------------------------------- */

void Surf::collate_array_implicit(int nrow, int ncol, surfint *tally2surf,
                                  double **in, double **out)
{
  int i,j,m,icell;
  cellint cellID;

  int me = comm->me;
  int nprocs = comm->nprocs;

  // create a grid cell hash for only my owned cells

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  MyCellHash hash;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    hash[cells[icell].id] = icell;
  }

  // for implicit surfs, tally2surf stores cellIDs

  cellint *tally2cell = (cellint *) tally2surf;

  // if I own tally grid cell, sum tallies to out directly
  // else nsend = # of tallies to contribute to rendezvous

  int nsend = 0;
  for (i = 0; i < nrow; i++) {
    if (hash.find(tally2cell[i]) == hash.end()) nsend++;
    else {
      icell = hash[tally2cell[i]];
      for (j = 0; j < ncol; j++)
        out[icell][j] += in[i][j];
    }
  }

  // done if just one proc

  if (nprocs == 1) return;

  // ncell = # of owned grid cells with implicit surfs, excluding sub cells
  // NOTE: could limit to cell group of caller

  int ncell = 0;
  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsurf <= 0) continue;
    if (cells[icell].nsplit <= 0) continue;
    ncell++;
  }

  // allocate memory for rvous input
  // ncount = ncell + nsend
  // ncol+2 doubles for each input = proc, cellID, ncol values

  int ncount = ncell + nsend;

  int *proclist;
  double *in_rvous;
  memory->create(proclist,ncount,"surf:proclist");
  memory->create(in_rvous,ncount*(ncol+2),"surf:in_rvous");

  // create rvous inputs
  // owning proc for each datum = random hash of cellID
  // flavor 1: one per ncell with proc and cellID, no tallies
  // flavor 2: one per nsend with proc = -1, cellID, tallies

  ncount = m = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsurf <= 0) continue;
    if (cells[icell].nsplit <= 0) continue;
    proclist[ncount] = hashlittle(&cells[icell].id,sizeof(cellint),0) % nprocs;
    in_rvous[m++] = me;
    in_rvous[m++] = cells[icell].id;    // NOTE: should use ubuf
    for (j = 0; j < ncol; j++)
      in_rvous[m++] = 0.0;
    ncount++;
  }

  for (i = 0; i < nrow; i++) {
    if (hash.find(tally2cell[i]) == hash.end()) {
      proclist[ncount] = hashlittle(&tally2cell[i],sizeof(cellint),0) % nprocs;
      in_rvous[m++] = -1;
      in_rvous[m++] = tally2cell[i];    // NOTE: should use ubuf
      for (j = 0; j < ncol; j++)
        in_rvous[m++] = in[i][j];
      ncount++;
    }
  }

  // perform rendezvous operation

  ncol_rvous = ncol;
  char *buf;
  int nout = comm->rendezvous(1,ncount,(char *) in_rvous,
                              (ncol+2)*sizeof(double),
                              0,proclist,rendezvous_implicit,
                              0,buf,(ncol+1)*sizeof(double),(void *) this);
  double *out_rvous = (double *) buf;

  memory->destroy(proclist);
  memory->destroy(in_rvous);

  // sum tallies returned for grid cells I own into out

  m = 0;
  for (i = 0; i < nout; i++) {
    cellID = out_rvous[m++];      // NOTE: should use ubuf
    icell = hash[cellID] - 1;     // subtract one for child cell index
    for (j = 0; j < ncol; j++)
      out[icell][j] += out_rvous[m++];
  }

  // clean-up

  memory->destroy(out_rvous);
}

/* ----------------------------------------------------------------------
   callback from collate_vector/array_implicit rendezvous operations
   create summed tallies for each grid cell assigned to me
   inbuf = list of N input datums
   send cellID + Ncol values back to owning proc of each grid cell
------------------------------------------------------------------------- */

int Surf::rendezvous_implicit(int n, char *inbuf,
                              int &flag, int *&proclist, char *&outbuf, void *ptr)
{
  int i,j,k,m,proc,iout;
  cellint cellID;

  Surf *sptr = (Surf *) ptr;
  Memory *memory = sptr->memory;
  int ncol = sptr->ncol_rvous;

  // scan inbuf for (proc,cellID) entries
  // create phash so can lookup the proc for each cellID

  double *in_rvous = (double *) inbuf;
  MyCellHash phash;

  m = 0;
  for (i = 0; i < n; i++) {
    proc = static_cast<int> (in_rvous[m++]);
    cellID = static_cast<cellint> (in_rvous[m++]);
    if (proc >= 0 && phash.find(cellID) == phash.end()) phash[cellID] = proc;
    m += ncol;
  }

  // allocate proclist & outbuf, based on size of max-size of phash

  int nmax = phash.size();
  memory->create(proclist,nmax,"surf:proclist");
  double *out;
  memory->create(out,nmax*(ncol+1),"surf:out");

  // scan inbuf for (cellID,tallies) entries
  // create a 2nd hash so can lookup the outbuf entry for each cellID
  // create proclist and outbuf with summed tallies for every cellID

  MyCellHash ohash;

  int nout = 0;
  k = m = 0;

  for (i = 0; i < n; i++) {
    proc = static_cast<int> (in_rvous[m++]);
    cellID = static_cast<cellint> (in_rvous[m++]);
    if (proc >= 0) {
      m += ncol;                         // skip entries with novalues
      continue;
    }
    if (ohash.find(cellID) == phash.end()) {
      ohash[cellID] = nout;              // add a new set of out values
      proclist[nout] = phash[cellID];
      out[k++] = cellID;
      for (j = 0; j < ncol; j++)
        out[k++] = in_rvous[m++];
      nout++;
    } else {
      iout = ohash[cellID] * (ncol+1);   // offset into existing out values
      iout++;                            // skip cellID;
      for (j = 0; j < ncol; j++)
        out[iout++] += in_rvous[m++];    // sum to existing values
    }
  }

  // flag = 2: new outbuf

  flag = 2;
  outbuf = (char *) out;
  return nout;
}

// ----------------------------------------------------------------------
// redistribute operations for explicit distributed surfs after clipping
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   redistribute newly created distributed lines to owing procs
   nold = original nown value before new surfs were read in
   nown = current nown value that includes my new surfs to redistribute
   nnew = nown value after new surfs from all procs are assigned to me
   called by ReadSurf:clip() after proc creates new surfs via clipping
   only called for explicit distributed surfs
------------------------------------------------------------------------- */

void Surf::redistribute_lines_clip(int nold, int nnew)
{
  // allocate memory for rvous input

  int nsend = nown - nold;

  int *proclist;
  memory->create(proclist,nsend,"surf:proclist");
  Line *in_rvous = (Line *) memory->smalloc(nsend*sizeof(Line),"surf:in_rvous");

  // create rvous inputs
  // proclist = owner of each surf = (id-1) % nprocs

  surfint id;

  int i = nold;
  for (int m = 0; m < nsend; m++) {
    id = mylines[i].id;
    proclist[m] = (id-1) % nprocs;
    memcpy(&in_rvous[m],&mylines[i],sizeof(Line));
    i++;
  }

  // insure mylines is allocated sufficient for new lines
  // reset nown to new value after rendezvous

  if (nnew > maxown) {
    int old = maxown;
    maxown = nnew;
    grow_own(old);
  }
  nown = nnew;

  // perform rendezvous operation
  // each proc owns subset of new surfs
  // receives them from other procs

  char *buf;
  int nout = comm->rendezvous(1,nsend,(char *) in_rvous,sizeof(Line),
                              0,proclist,rendezvous_lines,
                              0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->sfree(in_rvous);
}

/* ----------------------------------------------------------------------
   redistribute newly created distributed lines to owing procs
   nnew = nown value after new surfs from all procs are assigned to me
   called by ReadSurf:read_multiple()
   only called for explicit distributed surfs
------------------------------------------------------------------------- */

void Surf::redistribute_lines_temporary(int nnew)
{
  // allocate memory for rvous input

  int nsend = ntmp;

  int *proclist;
  memory->create(proclist,nsend,"surf:proclist");
  Line *in_rvous = (Line *) memory->smalloc(nsend*sizeof(Line),"surf:in_rvous");

  // create rvous inputs
  // proclist = owner of each surf = (id-1) % nprocs

  surfint id;

  for (int i = 0; i < nsend; i++) {
    id = tmplines[i].id;
    proclist[i] = (id-1) % nprocs;
    memcpy(&in_rvous[i],&tmplines[i],sizeof(Line));
  }

  // insure mylines is allocated sufficient for new lines
  // reset nown to new value after rendezvous

  if (nnew > maxown) {
    int old = maxown;
    maxown = nnew;
    grow_own(old);
  }
  nown = nnew;

  // perform rendezvous operation
  // each proc owns subset of new surfs
  // receives them from other procs

  char *buf;
  int nout = comm->rendezvous(1,nsend,(char *) in_rvous,sizeof(Line),
                              0,proclist,rendezvous_lines,
                              0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->sfree(in_rvous);
}

/* ----------------------------------------------------------------------
   callback from redistributed_lines_temporary rendezvous operation
   store received surfs assigned to me in correct location in mylines
   inbuf = list of N Inbuf datums
   no outbuf
------------------------------------------------------------------------- */

int Surf::rendezvous_lines(int n, char *inbuf,
                           int &flag, int *&proclist, char *&outbuf,
                           void *ptr)
{
  int i,j,k,m;

  Surf *sptr = (Surf *) ptr;
  Line *lines = sptr->mylines;
  int nprocs = sptr->comm->nprocs;
  int me = sptr->comm->me;

  // zero my owned surf values

  Line *in_rvous = (Line *) inbuf;
  surfint id;

  for (int i = 0; i < n; i++) {
    id = in_rvous[i].id;
    m = (id-1-me) / nprocs;
    memcpy(&lines[m],&in_rvous[i],sizeof(Line));
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   redistribute newly created distributed tris to owing procs
   nold = original nown value before new surfs were read in
   nown = current nown value that includes my new surfs to redistribute
   nnew = nown value after new surfs from all procs are assigned to me
   old = starting index that skips previously distributed surfs
   called by ReadSurf:clip() after proc create new surfs via clipping
   only called for explicit distributed surfs
------------------------------------------------------------------------- */

void Surf::redistribute_tris_clip(int nold, int nnew)
{
  // allocate memory for rvous input

  int nsend = nown - nold;

  int *proclist;
  memory->create(proclist,nsend,"surf:proclist");
  Tri *in_rvous = (Tri *) memory->smalloc(nsend*sizeof(Tri),"surf:in_rvous");

  // create rvous inputs
  // proclist = owner of each surf = (id-1) % nprocs

  surfint id;

  int i = nold;
  for (int m = 0; m < nsend; m++) {
    id = mytris[i].id;
    proclist[m] = (id-1) % nprocs;
    memcpy(&in_rvous[m],&mytris[i],sizeof(Tri));
    i++;
  }

  // insure mytris is allocated sufficient for new tris
  // reset nown to new value after rendezvous

  if (nnew > maxown) {
    int old = maxown;
    maxown = nnew;
    grow_own(old);
  }
  nown = nnew;

  // perform rendezvous operation
  // each proc owns subset of new surfs
  // receives them from other procs

  char *buf;
  int nout = comm->rendezvous(1,nsend,(char *) in_rvous,sizeof(Tri),
                              0,proclist,rendezvous_tris,
                              0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->sfree(in_rvous);
}

/* ----------------------------------------------------------------------
   redistribute newly created distributed tris to owing procs
   nnew = nown value after new surfs from all procs are assigned to me
   called by ReadSurf:read_multiple()
   only called for explicit distributed surfs
------------------------------------------------------------------------- */

void Surf::redistribute_tris_temporary(int nnew)
{
  // allocate memory for rvous input

  int nsend = ntmp;

  int *proclist;
  memory->create(proclist,nsend,"surf:proclist");
  Tri *in_rvous = (Tri *) memory->smalloc(nsend*sizeof(Tri),"surf:in_rvous");

  // create rvous inputs
  // proclist = owner of each surf = (id-1) % nprocs

  surfint id;

  for (int i = 0; i < nsend; i++) {
    id = tmptris[i].id;
    proclist[i] = (id-1) % nprocs;
    memcpy(&in_rvous[i],&tmptris[i],sizeof(Tri));
  }

  // insure mytris is allocated sufficient for new tris
  // reset nown to new value after rendezvous

  if (nnew > maxown) {
    int old = maxown;
    maxown = nnew;
    grow_own(old);
  }
  nown = nnew;

  // perform rendezvous operation
  // each proc owns subset of new surfs
  // receives them from other procs

  char *buf;
  int nout = comm->rendezvous(1,nsend,(char *) in_rvous,sizeof(Tri),
                              0,proclist,rendezvous_tris,
                              0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->sfree(in_rvous);
}

/* ----------------------------------------------------------------------  
   callback from redistributed_tris_temporary rendezvous operation
   store received surfs assigned to me in correct location in mytris
   inbuf = list of N Inbuf datums
   no outbuf
------------------------------------------------------------------------- */

int Surf::rendezvous_tris(int n, char *inbuf,
                          int &flag, int *&proclist, char *&outbuf,
                          void *ptr)
{
  int i,j,k,m;

  Surf *sptr = (Surf *) ptr;
  Tri *tris = sptr->mytris;
  int nprocs = sptr->comm->nprocs;
  int me = sptr->comm->me;

  // zero my owned surf values

  Tri *in_rvous = (Tri *) inbuf;
  surfint id;

  for (int i = 0; i < n; i++) {
    id = in_rvous[i].id;
    m = (id-1-me) / nprocs;
    memcpy(&tris[m],&in_rvous[i],sizeof(Tri));
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   spread values for N datums per surf from in vector to out vector
   type = 0/1 = INT or DOUBLE datums
   for explicit non-distributed surfs, use allreduce algorithm
     in = owned surfs
     out = nlocal surfs (all procs own copy)
   for explicit distributed surfs, use rendezvous algorithms
     in = owned surfs
     out = nlocal+nghost surfs (for my owned + ghost cells)
   only called for explicit surfs, non-distributed or distributed
   called by spread_custom() and surf_collide.cpp
------------------------------------------------------------------------- */

void Surf::spread_own2local(int n, int type, void *in, void *out)
{
  if (!distributed) spread_own2local_reduce(n,type,in,out);
  else spread_own2local_rendezvous(n,type,in,out);
}

/* ----------------------------------------------------------------------
   allreduce algorithm for spread_own2local
   owned surf data combined for data for all nlocal surfs on all procs
   NOTE: check for overflow of n*nall or n*nall*sizeof()
---------------------------------------------------------------------- */

void Surf::spread_own2local_reduce(int n, int type, void *in, void *out)
{
  int i,j,m,ij,mj;
  
  if (type == INT) {
    int *ivec = (int *) in;
    int *ovec = (int *) out;

    int *myvec;
    memory->create(myvec,n*nlocal,"surf/spread:myvec");
    memset(myvec,0,n*nlocal*sizeof(int));

    if (n == 1) {
      for (i = 0; i < nown; i++) {
	m = me + i*nprocs;
	myvec[m] = ivec[i];
      }
    } else {
      for (i = 0; i < nown; i++) {
	ij = i * n;
	mj = (me + i*nprocs) * n;
	for (j = 0; j < n; j++)
	  myvec[mj++] = ivec[ij++];
      }
    }
    
    MPI_Allreduce(myvec,ovec,n*nlocal,MPI_INT,MPI_SUM,world);

    memory->destroy(myvec);
    
  } else {
    double *ivec = (double *) in;
    double *ovec = (double *) out;

    double *myvec;
    memory->create(myvec,n*nlocal,"surf/spread:myvec");
    memset(myvec,0,n*nlocal*sizeof(double));

    if (n == 1) {
      for (i = 0; i < nown; i++) {
	m = me + i*nprocs;
	myvec[m] = ivec[i];
      }
    } else {
      for (i = 0; i < nown; i++) {
	ij = i * n;
	mj = (me + i*nprocs) * n;
	for (j = 0; j < n; j++)
	  myvec[mj++] = ivec[ij++];
      }
    }
    
    MPI_Allreduce(myvec,ovec,n*nlocal,MPI_DOUBLE,MPI_SUM,world);

    memory->destroy(myvec);
  }
}

/* ----------------------------------------------------------------------
   rendezvous algorithm for spread_own2local
   owned surf data becomes data for nlocal+nghost surfs on each proc
---------------------------------------------------------------------- */

void Surf::spread_own2local_rendezvous(int n, int type, void *in, void *out)
{
  int i,j,k,m,index;

  // allocate memory for rvous input

  int nall = nlocal + nghost;

  int *proclist;
  memory->create(proclist,nall,"spread/own2local:proclist");
  int *inbuf;
  memory->create(inbuf,3*nall,"spread/own2local:inbuf");

  // create rvous inputs
  // proclist = owner of each surf
  // inbuf = 3 ints per request = me, my index, index on owning proc

  int dim = domain->dimension;
  surfint surfID;

  m = 0;
  for (i = 0; i < nall; i++) {
    if (dim == 2) surfID = lines[i].id;
    else surfID = tris[i].id;
    proclist[i] = (surfID-1) % nprocs;
    inbuf[m] = me;
    inbuf[m+1] = i;
    inbuf[m+2] = (surfID-1) / nprocs;
    m += 3;
  }

  // perform rendezvous operation
  // each proc owns subset of surfs
  // receives all surf requests to return per-surf values to each proc who needs it

  spread_type = type;
  spread_size = n;
  spread_data = (void *) in;
  
  int outbytes;
  if (type == INT) outbytes = (n+1) * sizeof(int);
  else outbytes = (n+1) * sizeof(double);
  char *buf;

  int nreturn = comm->rendezvous(1,nall,(char *) inbuf,3*sizeof(int),
				 0,proclist,rendezvous_spread,
				 0,buf,outbytes,(void *) this);
  
  memory->destroy(proclist);
  memory->destroy(inbuf);

  // loop over received datums for nlocal+nghost surfs
  // copy per-surf values into out

  int *ibuf,*ioutbuf;
  double *dbuf,*doutbuf;

  if (type == INT) ibuf = (int *) buf;
  else if (type == DOUBLE) dbuf = (double *) dbuf;

  m = 0;
  if (type == INT) {
    for (i = 0; i < nall; i++) {
      index = ibuf[m++];
      if (n == 1)
	ioutbuf[index] = ibuf[m++];
      else {
	k = index * n;
	for (j = 0; j < n; j++)
	  ioutbuf[k++] = ibuf[m++];
      }
    }
    
  } else if (type == DOUBLE) {
    for (i = 0; i < nall; i++) {
      index = (int) ubuf(dbuf[m++]).i;
      if (n == 1)
	doutbuf[index] = dbuf[m++];
      else {
	k = index * n;
	for (j = 0; j < n; j++)
	  doutbuf[k++] = dbuf[m++];
      }
    }
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   callback from spread_own2local_rendezvous rendezvous operation
   process requests for data from surf elements I own
------------------------------------------------------------------------- */

int Surf::rendezvous_spread(int n, char *inbuf,
			    int &flag, int *&proclist, char *&outbuf,
			    void *ptr)
{
  int i,j,k,m,id;
  int *idata,*iout;
  double *ddata,*dout;

  Surf *sptr = (Surf *) ptr;
  Memory *memory = sptr->memory;
  int type = sptr->spread_type;
  int size = sptr->spread_size;
  if (type == INT) idata = (int *) sptr->spread_data;
  else if (type == DOUBLE) ddata = (double *) sptr->spread_data;

  // allocate proclist & iout/dout, based on n, type, size

  memory->create(proclist,n,"spread:proclist");
  if (type == INT) memory->create(iout,(size+1)*n,"spread:iout");
  else if (type == DOUBLE) memory->create(dout,(size+1)*n,"spread:dout");

  // loop over received requests, pack data into iout/dout
  
  int *in_rvous = (int *) inbuf;
  int oproc,oindex,index;
  m = 0;
  k = 0;
  
  if (type == INT) {
    for (int i = 0; i < n; i++) {
      oproc = in_rvous[m];
      oindex = in_rvous[m+1];
      index = in_rvous[m+2];
      m += 3;
      
      proclist[i] = oproc;
      iout[k++] = oindex;
      if (size == 1)
	iout[k++] = idata[index];
      else {
	id = index * size;
	for (j = 0; j < size; j++)
	  iout[k++] = idata[id++];
      }
    }
    
  } else if (type == DOUBLE) {
    for (int i = 0; i < n; i++) {
      oproc = in_rvous[m];
      oindex = in_rvous[m+1];
      index = in_rvous[m+2];
      m += 3;
      
      proclist[i] = oproc;
      dout[k++] = oindex;
      if (size == 1)
	dout[k++] = ddata[index];
      else {
	id = index * size;
	for (j = 0; j < size; j++)
	  dout[k++] = ddata[id++];
      }
    }
  }

  // flag = 2: new outbuf

  flag = 2;
  if (type == INT) outbuf = (char *) iout;
  else if (type == DOUBLE) outbuf = (char *) dout;
  return n;
}
