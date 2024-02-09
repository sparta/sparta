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
#include "comm.h"
#include "grid.h"
#include "memory.h"
#include "error.h"
#include "string.h"

using namespace SPARTA_NS;

enum{TALLYAUTO,TALLYREDUCE,TALLYRVOUS};         // same as Update

// ----------------------------------------------------------------------
// collate operations for explicit surfs
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
  memory->create(one,nglobal,"collate/vec:one");
  memory->create(all,nglobal,"collate/vec:all");

  // zero all values and add in values I accumulated

  bigint bbytes = (bigint) nglobal * sizeof(double);
  memset(one,0,bbytes);

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
  int k,m;

  // allocate memory for rvous input

  int *proclist;
  memory->create(proclist,nrow,"collate/vec:proclist");
  double *dbuf;
  memory->create(dbuf,2*nrow,"collate/vec:dbuf");

  // create rvous inputs
  // proclist = owner of each surf
  // logic of (id-1) % nprocs sends
  //   surf IDs 1,11,21,etc on 10 procs to proc 0

  surfint id;

  k = 0;
  m = 0;
  for (int i = 0; i < nrow; i++) {
    id = tally2surf[i];
    proclist[i] = (id-1) % nprocs;
    dbuf[k++] = ubuf(id).d;
    dbuf[k++] = in[m];
    m += instride;
  }

  // perform rendezvous operation
  // each proc owns subset of surfs
  // receives all tally contributions to surfs it owns

  out_rvous = out;

  char *buf;
  int nout = comm->rendezvous(1,nrow,(char *) dbuf,2*sizeof(double),
                              0,proclist,rendezvous_vector,
                              0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->destroy(dbuf);
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
  int i,k,m;

  Surf *sptr = (Surf *) ptr;
  int nown = sptr->nown;
  double *out = sptr->out_rvous;
  int nprocs = sptr->comm->nprocs;

  // zero my owned surf values

  bigint bbytes = (bigint) nown * sizeof(double);
  memset(out,0,bbytes);

  // accumulate per-surf values from different procs to my owned surfs

  double *dbuf = (double *) inbuf;
  surfint id;

  k = 0;
  for (i = 0; i < n; i++) {
    id = (surfint) ubuf(dbuf[k++]).i;
    m = (id-1) / nprocs;
    out[m] += dbuf[k++];
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

  bigint bbytes = (bigint) nglobal * ncol * sizeof(double);
  memset(&one[0][0],0,bbytes);

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

  memory->destroy(one);
  memory->destroy(all);
}

/* ----------------------------------------------------------------------
   rendezvous algorithm for collate_array
------------------------------------------------------------------------- */

void Surf::collate_array_rendezvous(int nrow, int ncol, surfint *tally2surf,
                                    double **in, double **out)
{
  int j,m;

  // allocate memory for rvous input

  int *proclist;
  memory->create(proclist,nrow,"surf:proclist");
  double *in_rvous = (double *)
    memory->smalloc(nrow*(ncol+1)*sizeof(double*),"surf:in_rvous");

  // create rvous inputs
  // proclist = owner of each surf
  // logic of (id-1) % nprocs sends
  //   surf IDs 1,11,21,etc on 10 procs to proc 0

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
  int j,k,m;

  Surf *sptr = (Surf *) ptr;
  int nown = sptr->nown;
  int ncol = sptr->ncol_rvous;
  double *out = sptr->out_rvous;
  int nprocs = sptr->comm->nprocs;

  // zero my owned surf values

  bigint bbytes = (bigint) nown * ncol * sizeof(double);
  memset(out,0,bbytes);

  // accumulate per-surf values from different procs to my owned surfs

  double *in_rvous = (double *) inbuf;
  surfint id;

  m = 0;
  for (int i = 0; i < n; i++) {
    id = (surfint) ubuf(in_rvous[m++]).i;
    k = (id-1) / nprocs * ncol;
    for (j = 0; j < ncol; j++)
      out[k++] += in_rvous[m++];
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

// OLD methods for implicit surf collates, now done in Grid class
// these used All2All rendezvous
// Grid class uses irregular rendezvous when appropriate (faster)

/* ----------------------------------------------------------------------
   comm of tallies across all procs for a single value (vector)
   called from compute isurf/grid and fix ave/grid
     for implicit surf tallies by grid cell
   nrow = # of tallies
   tally2surf = surf ID for each tally (same as cell ID)
   in = vectir of tally values
   return out = summed tallies for grid cells I own
     caller MUST zero out before call, b/c this method does not know size
   always done via rendezvous algorithm (no allreduce algorithm)
   only called for implicit surfs
------------------------------------------------------------------------- */

/*
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
    in_rvous[m++] = ubuf(me).d;
    in_rvous[m++] = ubuf(cells[icell].id).d;
    in_rvous[m++] = 0.0;
    ncount++;
  }

  for (i = 0; i < nrow; i++) {
    if (hash.find(tally2cell[i]) == hash.end()) {
      proclist[ncount] = hashlittle(&tally2cell[i],sizeof(cellint),0) % nprocs;
      in_rvous[m++] = ubuf(-1).d;
      in_rvous[m++] = ubuf(tally2cell[i]).d;
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
    cellID = (cellint) ubuf(out_rvous[m++]).u;
    icell = hash[cellID];
    out[icell] += out_rvous[m++];
  }

  // clean-up

  memory->destroy(out_rvous);
}
*/

/* ----------------------------------------------------------------------
   comm of tallies across all procs for multiple values (array)
   called from compute isurf/grid and fix ave/grid
     for implicit surf tallies by grid cell
   nrow = # of tallies
   ncol = # of values per tally
   tally2surf = surf ID for each tally (same as cell ID)
   in = array of tally values, nrow by ncol
   return out = summed tallies for grid cells I own, nlocal by ncol
     caller MUST zero out before call, b/c this method does not know size
   always done via rendezvous algorithm (no allreduce algorithm)
   only called for implicit surfs
------------------------------------------------------------------------- */

/*
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
    in_rvous[m++] = ubuf(me).d;
    in_rvous[m++] = ubuf(cells[icell].id).d;
    for (j = 0; j < ncol; j++)
      in_rvous[m++] = 0.0;
    ncount++;
  }

  for (i = 0; i < nrow; i++) {
    if (hash.find(tally2cell[i]) == hash.end()) {
      proclist[ncount] = hashlittle(&tally2cell[i],sizeof(cellint),0) % nprocs;
      in_rvous[m++] = ubuf(-1).d;
      in_rvous[m++] = ubuf(tally2cell[i]).d;
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
    cellID = (cellint) ubuf(out_rvous[m++]).u;
    icell = hash[cellID];
    for (j = 0; j < ncol; j++)
      out[icell][j] += out_rvous[m++];
  }

  // clean-up

  memory->destroy(out_rvous);
}
*/

/* ----------------------------------------------------------------------
   callback from collate_vector/array_implicit rendezvous operations
   create summed tallies for each grid cell assigned to me
   inbuf = list of N input datums
   send cellID + Ncol values back to owning proc of each grid cell
------------------------------------------------------------------------- */

/*
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
    proc = (int) ubuf(in_rvous[m++]).i;
    cellID = (cellint) ubuf(in_rvous[m++]).u;
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
    proc = (int) ubuf(in_rvous[m++]).i;
    cellID = (cellint) ubuf(in_rvous[m++]).u;
    if (proc >= 0) {
      m += ncol;                         // skip entries with novalues
      continue;
    }
    if (ohash.find(cellID) == phash.end()) {
      ohash[cellID] = nout;              // add a new set of out values
      proclist[nout] = phash[cellID];
      out[k++] = ubuf(cellID).d;
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
*/
