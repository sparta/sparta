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

#include "stdlib.h"
#include "string.h"
#include "surf.h"
#include "domain.h"
#include "comm.h"
#include "grid.h"
#include "update.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};                      // several files

// ----------------------------------------------------------------------
// redistribution of new surfs operation from ReadSurf or RemoveSurf
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   redistribute caller per-proc surfs to owning procs in Surf data structs
     for explicit surfs only, all or distributed
     use rendezvous algorithm to send surfs to owning procs
   result:
     all: Allgatherv gives all procs a copy in their lines/tris
     distributed: each proc has its owned subset in mylines/mytris
   n = # of new surfs contributed by this proc
   newlines/newtris = list of new surfs, other ptr is NULL
   nc = # of custom attributes (vecs/arrays) for each surf
   index_custom = index for each custom vec or array in Surf custom lists
   cvalues = custom values for each surf in same order as newlines/newtris
     1st value is surf ID, remaining values are for nc vecs/arrays
   nsurf_new, nsurf_old = global #
   called from add_surfs() via ReadSurf or RemoveSurf or ReadRestart
------------------------------------------------------------------------- */

void Surf::redistribute_surfs(int n, Line *newlines, Tri *newtris,
			      int nc, int *index_custom, double **cvalues,
			      bigint nsurf_new, bigint nsurf_old)
{
  int dim = domain->dimension;

  int nbytes;
  if (dim == 2) nbytes = sizeof(Line);
  else nbytes = sizeof(Tri);

  // ---------------------------------------------------------
  // redistribute newlines or newtris via rendezvous operation
  // ---------------------------------------------------------

  // setup rendezvous memory
  // for all:
  //   target is local lines_contig or tris_contig for contiguous IDs
  //   ncontig = # of contiguous ID surfs this proc owns
  //   ncontig = nsurf_all / nprocs, last proc owns any extra surfs
  // for distributed:
  //   target is mylines/mytris owning strided IDs
  //   no need to allocate here

  int *proclist;
  memory->create(proclist,n,"readsurf:proclist");

  int surfperproc = nsurf_new / nprocs;
  int ncontig = surfperproc;
  if (me == nprocs-1) ncontig = nsurf_new - (nprocs-1)*surfperproc;

  Line *lines_contig = NULL;
  Tri *tris_contig = NULL;
  bigint bbytes = (bigint) ncontig * nbytes;

  if (!distributed) {
    if (dim == 2) {
      lines_contig = (Line *)
	memory->smalloc(bbytes,"surf:lines_contig");
      memset(lines_contig,0,bbytes);
    } else {
      tris_contig = (Tri *)
	memory->smalloc(bbytes,"surf:tris_contig");
      memset(tris_contig,0,bbytes);
    }
  }

  // create rvous inputs
  // proclist = owner of each surf based on surf ID
  // receivers own contiguous IDs for all, strided for distributed

  surfint id;
  int proc;

  for (int i = 0; i < n; i++) {
    if (dim == 2) id = newlines[i].id;
    else id = newtris[i].id;
    if (!distributed) {
      if (surfperproc) {
        proc = (id-nsurf_old-1) / surfperproc;
        if (proc >= nprocs) proc = nprocs-1;
      } else proc = nprocs - 1;
    } else if (distributed) {
      proc = (id-1) % nprocs;
    }
    proclist[i] = proc;
  }

  // perform rendezvous operation
  // each proc receives subset of new surfs from other procs

  redistribute_nsurf_old = nsurf_old;
  redistribute_surfperproc = surfperproc;
  redistribute_lines_contig = lines_contig;
  redistribute_tris_contig = tris_contig;

  char *in_rvous;
  if (dim == 2) in_rvous = (char *) newlines;
  else in_rvous = (char *) newtris;

  char *buf;
  int nout = comm->rendezvous(1,n,in_rvous,nbytes,
                              0,proclist,rendezvous_redistribute_surfs,
                              0,buf,0,(void *) this);

  memory->destroy(proclist);

  // if distribution is all:
  // perform Allgatherv with lines/tris_contig
  // target is to append to previous surfs in lines/tris

  if (!distributed) {
    int *recvcounts,*displs;
    memory->create(recvcounts,nprocs,"readsurf:recvcounts");
    memory->create(displs,nprocs,"readsurf:displs");

    int nbsize;
    if (dim == 2) nbsize = ncontig * sizeof(Line);
    else nbsize = ncontig * sizeof(Tri);

    MPI_Allgather(&nbsize,1,MPI_INT,recvcounts,1,MPI_INT,world);
    displs[0] = 0;
    for (int i = 1; i < nprocs; i++) displs[i] = displs[i-1] + recvcounts[i-1];

    if (dim == 2)
      MPI_Allgatherv(lines_contig,ncontig*sizeof(Line),MPI_CHAR,
		     &lines[nsurf_old],recvcounts,displs,MPI_CHAR,world);
    else
      MPI_Allgatherv(tris_contig,ncontig*sizeof(Tri),MPI_CHAR,
		     &tris[nsurf_old],recvcounts,displs,MPI_CHAR,world);

    // clean up

    memory->sfree(lines_contig);
    memory->sfree(tris_contig);
    memory->destroy(recvcounts);
    memory->destroy(displs);
  }

  // done if no new custom data

  if (nc == 0) return;

  // --------------------------------------------------------------
  // redistribute new custom per-surf data via rendezvous operation
  // --------------------------------------------------------------

  // proclist = owner of each surf based on surf ID
  // surf IDs in cvalues are already augmented by old surfs

  memory->create(proclist,n,"surf:proclist");

  for (int i = 0; i < n; i++) {
    id = (surfint) ubuf(cvalues[i][0]).i;
    proclist[i] = (id-1) % nprocs;
  }

  // nvalues_custom = # of new custom values per surf

  int nvalues_custom = 0;
  for (int ic = 0; ic < nc; ic++) {
    if (esize[ic] == 0) nvalues_custom++;
    else nvalues_custom += esize[ic];
  }

  // perform rendezvous operation
  // each proc receives subset of custom values from other procs

  redistribute_nvalues_custom = nvalues_custom;
  redistribute_index_custom = index_custom;

  in_rvous = NULL;
  if (n) in_rvous = (char *) &cvalues[0][0];
  nbytes = (1+nvalues_custom) * sizeof(double);

  nout = comm->rendezvous(1,n,in_rvous,nbytes,
			  0,proclist,rendezvous_redistribute_custom,
			  0,buf,0,(void *) this);

  memory->destroy(proclist);
}

/* ----------------------------------------------------------------------
   callback from redistribute_surfs rendezvous operatiion
   inbuf = list of N datums, each is a line or tri
   no outbuf
------------------------------------------------------------------------- */

int Surf::rendezvous_redistribute_surfs(int n, char *inbuf, int &flag,
					int *&proclist, char *&outbuf, void *ptr)
{
  Surf *sptr = (Surf *) ptr;

  // generic Surf class variables

  int me = sptr->me;
  int nprocs = sptr->nprocs;
  int distributed = sptr->distributed;
  Line *mylines = sptr->mylines;
  Tri *mytris = sptr->mytris;
  int dim = sptr->domain->dimension;

  // Surf class variables peculiar to this rendevous operation

  bigint nsurf_old = sptr->redistribute_nsurf_old;
  int surfperproc = sptr->redistribute_surfperproc;
  Line *lines_contig = sptr->redistribute_lines_contig;
  Tri *tris_contig = sptr->redistribute_tris_contig;

  int nbytes;
  if (dim == 2) nbytes = sizeof(Line);
  else nbytes = sizeof(Tri);

  // loop over received surfs
  // copy each into apprpriate location in Surf data struct
  // for all: target = lines/tris_contig
  // for distributed: target = mylines/mytris

  Line *in_lines,*out_lines;
  Tri *in_tris,*out_tris;

  if (dim == 2) {
    in_lines = (Line *) inbuf;
    if (!distributed) out_lines = lines_contig;
    else out_lines = mylines;
  } else {
    in_tris = (Tri *) inbuf;
    if (!distributed) out_tris = tris_contig;
    else out_tris = mytris;
  }

  surfint id;
  int m;

  for (int i = 0; i < n; i++) {
    if (dim == 2) id = in_lines[i].id;
    else id = in_tris[i].id;
    if (distributed) m = (id-1) / nprocs;
    else m = (id - nsurf_old - 1) - me*surfperproc;
    if (dim == 2) memcpy(&out_lines[m],&in_lines[i],nbytes);
    else memcpy(&out_tris[m],&in_tris[i],nbytes);
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

/* ----------------------------------------------------------------------
   callback from redistribute_surfs rendezvous operatiion
   inbuf = list of N datums, each is a set of custom values for one surf
   no outbuf
------------------------------------------------------------------------- */

int Surf::rendezvous_redistribute_custom(int n, char *inbuf, int &flag,
					 int *&proclist, char *&outbuf, void *ptr)
{
  Surf *sptr = (Surf *) ptr;

  // generic Surf class variables

  int nprocs = sptr->nprocs;
  int ncustom = sptr->ncustom;
  int *etype = sptr->etype;
  int *esize = sptr->esize;
  int *ewhich = sptr->ewhich;
  int **eivec = sptr->eivec;
  int ***eiarray = sptr->eiarray;
  double **edvec = sptr->edvec;
  double ***edarray = sptr->edarray;

  // Surf class variables peculiar to this rendevous operation

  int nvalues_custom = sptr->redistribute_nvalues_custom;
  int *index_custom = sptr->redistribute_index_custom;

  // loop over received sets of custom values
  // copy each value into appropriate locataion in custom vec or array

  double *in_custom = (double *) inbuf;

  int i,j,k,m;
  int index,type,size;
  surfint id;

  int skip = 1 + nvalues_custom;
  int offset = 1;

  for (int ic = 0; ic < ncustom; ic++) {
    index = index_custom[ic];
    type = etype[index];
    size = esize[index];

    if (type == 0) {
      if (size == 0) {
	int *ivector = eivec[ewhich[index]];

	m = 0;
	for (i = 0; i < n; i++) {
	  id = (surfint) ubuf(in_custom[m]).i;
	  j = (id-1) / nprocs;
	  ivector[j] = static_cast<int> (in_custom[m+offset]);
	  m += skip;
	}

      } else {
	int **iarray = eiarray[ewhich[index]];

	m = 0;
	for (i = 0; i < n; i++) {
	  id = (surfint) ubuf(in_custom[m]).i;
	  j = (id-1) / nprocs;
	  for (k = 0; k < size; k++)
	    iarray[j][k] = static_cast<int> (in_custom[m+offset+k]);
	  m += skip;
	}
      }

    } else {
      if (size == 0) {
	double *dvector = edvec[ewhich[index]];

	m = 0;
	for (i = 0; i < n; i++) {
	  id = (surfint) ubuf(in_custom[m]).i;
	  j = (id-1) / nprocs;
	  dvector[j] = in_custom[m+offset];
	  m += skip;
	}
	
      } else {
	double **darray = edarray[ewhich[index]];

	m = 0;
	for (i = 0; i < n; i++) {
	  id = (surfint) ubuf(in_custom[m]).i;
	  j = (id-1) / nprocs;
	  for (k = 0; k < size; k++)
	    darray[j][k] = in_custom[m+offset+k];
	  m += skip;
	}
      }
    }

    if (size == 0) offset++;
    else offset += size;
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}

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

  // keep = 1 if a local surf is referenced by a compressed local grid cell

  int *keep;
  memory->create(keep,nlocal,"surf:keep");
  for (i = 0; i < nlocal; i++) keep[i] = 0;

  // convert grid cell csurfs to surf IDs so can reset after surf compression
  // skip cells with no surfs or sub-cells

  int dim = domain->dimension;
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
// spread operations for surf-style variables or custom per-surf vecs/arrays
// comm data from owned strided to nlocal+nghost for all or distributed
// comm data from nlocal+nghost to owned strided for all or distributed
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   spread values for N datums per surf from in vector to out vector
   type = 0/1 = INT or DOUBLE datums
   for explicit all surfs, use allreduce algorithm
     in = owned surfs
     out = nlocal surfs (all procs own copy)
   for explicit distributed surfs, use rendezvous algorithm
     in = owned surfs
     out = nlocal+nghost surfs (for my owned + ghost cells)
   only called for explicit surfs, all or distributed
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
---------------------------------------------------------------------- */

void Surf::spread_own2local_reduce(int n, int type, void *in, void *out)
{
  int i,j,m,ij,mj;

  if (type == INT) {
    int *ivec = (int *) in;
    int *ovec = (int *) out;

    int *myvec;
    bigint bcount = (bigint) nlocal * n;
    if (bcount > MAXSMALLINT)
      error->all(FLERR,"Overflow in spread_own2local_reduce");
    bigint bbytes = (bigint) nlocal * n * sizeof(int);

    memory->create(myvec,nlocal*n,"surf/spread:myvec");
    memset(myvec,0,bbytes);

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
    bigint bcount = (bigint) nlocal * n;
    if (bcount > MAXSMALLINT)
      error->all(FLERR,"Overflow in spread_own2local_reduce");
    bigint bbytes = (bigint) nlocal * n * sizeof(double);

    memory->create(myvec,nlocal*n,"surf/spread:myvec");
    memset(myvec,0,bbytes);

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

  int dim = domain->dimension;
  int nall = nlocal + nghost;

  int *proclist;
  memory->create(proclist,nall,"spread/own2local:proclist");
  int *inbuf;
  memory->create(inbuf,3*nall,"spread/own2local:inbuf");

  // create rvous inputs
  // proclist = owner of each surf
  // inbuf = 3 ints per request = me, my index, index on owning proc

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
  spread_data = in;

  int outbytes;
  if (type == INT) outbytes = (n+1) * sizeof(int);
  else outbytes = (n+1) * sizeof(double);
  char *buf;

  int nreturn = comm->rendezvous(1,nall,(char *) inbuf,3*sizeof(int),
				 0,proclist,rendezvous_own2local,
				 0,buf,outbytes,(void *) this);

  memory->destroy(proclist);
  memory->destroy(inbuf);

  // loop over received datums for nlocal+nghost surfs
  // copy per-surf values into out

  int *ibuf,*ioutbuf;
  double *dbuf,*doutbuf;

  if (type == INT) {
    ibuf = (int *) buf;
    ioutbuf = (int *) out;
  } else if (type == DOUBLE) {
    dbuf = (double *) buf;
    doutbuf = (double *) out;
  }

  m = 0;
  if (type == INT) {
    for (i = 0; i < nreturn; i++) {
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
    for (i = 0; i < nreturn; i++) {
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

int Surf::rendezvous_own2local(int n, char *inbuf,
			       int &flag, int *&proclist, char *&outbuf,
			       void *ptr)
{
  int j,k,m,idata;
  int *ibuf,*iout;
  double *dbuf,*dout;

  Surf *sptr = (Surf *) ptr;
  Memory *memory = sptr->memory;
  int type = sptr->spread_type;
  int size = sptr->spread_size;
  if (type == INT) ibuf = (int *) sptr->spread_data;
  else if (type == DOUBLE) dbuf = (double *) sptr->spread_data;

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
	iout[k++] = ibuf[index];
      else {
	idata = index * size;
	for (j = 0; j < size; j++)
	  iout[k++] = ibuf[idata++];
      }
    }

  } else if (type == DOUBLE) {
    for (int i = 0; i < n; i++) {
      oproc = in_rvous[m];
      oindex = in_rvous[m+1];
      index = in_rvous[m+2];
      m += 3;

      proclist[i] = oproc;
      dout[k++] = ubuf(oindex).d;
      if (size == 1)
	dout[k++] = dbuf[index];
      else {
	idata = index * size;
	for (j = 0; j < size; j++)
	  dout[k++] = dbuf[idata++];
      }
    }
  }

  // flag = 2: new outbuf

  flag = 2;
  if (type == INT) outbuf = (char *) iout;
  else if (type == DOUBLE) outbuf = (char *) dout;
  return n;
}

/* ----------------------------------------------------------------------
   set nunique and unique vector on each proc
   for distributed surfs:
     each proc owns Nsurf/Nprocs fraction
     each proc stores Nlocal+Nghost surfs which overlap its owned+ghost cells
     this means a single owned surf can be in Nlocal list of many procs
     if a single owned surf is in Nlocal list of M procs
       assign it uniquely to one of those M procs with 1/M random probability
   nunique = # of my Nlocal surfs which are assigned to a unique owned surfs
   unique = list of indices in my local list of those nunique surfs
   only need to call this method when Nlocal+Nghost list change (e.g. LB, adapt)
   only need to call this method when performing spread_inverse()
     to map local+ghost custom data back to owned surfs
   only called for distributed surfs
------------------------------------------------------------------------- */

void Surf::assign_unique()
{
  int i,m;

  // RNG for unique assignment among multiple copies

  if (!urandom) {
    urandom = new RanKnuth(update->ranmaster->uniform());
    double seed = update->ranmaster->uniform();
    urandom->reset(seed,me,100);
  }

  // allocate memory for rvous input

  int *proclist;
  memory->create(proclist,nlocal,"assign/unique:proclist");
  int *inbuf;
  memory->create(inbuf,3*nlocal,"assigh/unique:inbuf");

  // create rvous inputs
  // proclist = owner of each surf
  // inbuf = 3 ints per request = me, my index, index on owning proc

  int dim = domain->dimension;
  surfint surfID;

  for (i = 0; i < nlocal; i++) {
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
  // receives datum from each proc who onws a local copy of an owned surf

  int outbytes = sizeof(int);
  char *buf;

  int nout = comm->rendezvous(1,nlocal,(char *) inbuf,3*sizeof(int),
			      0,proclist,rendezvous_unique,
			      0,buf,outbytes,(void *) this);

  memory->destroy(proclist);
  memory->destroy(inbuf);

  // copy rendezvous output into unique vector
  // set corresponding surf IDs in uniqueID vector

  int *out_rvous = (int *) buf;

  nunique = nout;
  memory->destroy(unique);
  memory->destroy(uniqueID);
  memory->create(unique,nunique,"assign/unique:unique");
  memory->create(uniqueID,nunique,"assign/unique:uniqueID");

  for (i = 0; i < nunique; i++) {
    unique[i] = out_rvous[i];
    if (dim == 2) uniqueID[i] = lines[unique[i]].id;
    else uniqueID[i] = tris[unique[i]].id;
  }

  // clean-up

  memory->destroy(out_rvous);
}

/* ----------------------------------------------------------------------
   callback from assign_unique rendezvous operation
   process copies of surf elements I own
------------------------------------------------------------------------- */

int Surf::rendezvous_unique(int n, char *inbuf,
			    int &flag, int *&proclist, char *&outbuf,
			    void *ptr)
{
  int i,k,m;

  Surf *sptr = (Surf *) ptr;
  Memory *memory = sptr->memory;
  RanKnuth *urandom = sptr->urandom;
  int nown = sptr->nown;

  // allocate proclist & ownflags & duplicates

  int *ownflags,*duplicates,*selection;

  memory->create(proclist,nown,"unique/assign:proclist");
  memory->create(ownflags,nown,"unique/assign:ownflags");
  memory->create(duplicates,nown,"unique/assign:duplicates");
  memory->create(selection,nown,"unique/assign:selection");

  // loop over received requests
  // count duplicates of each owned surf

  bigint bbytes = (bigint) nown * sizeof(int);
  memset(duplicates,0,bbytes);

  int *in_rvous = (int *) inbuf;
  int index;

  m = 0;
  for (i = 0; i < n; i++) {
    index = in_rvous[m+2];
    duplicates[index]++;
    m += 3;
  }

  // select one randomly among duplicates for each owned surf

  for (i = 0; i < nown; i++)
    selection[i] = static_cast<int> (urandom->uniform() * duplicates[i]);

  // loop over received requests
  // recount to find selected duplicate of each owned surf
  // send ownflag datum back to owner of that duplicate

  memset(duplicates,0,bbytes);

  k = 0;
  m = 0;
  for (i = 0; i < n; i++) {
    index = in_rvous[m+2];
    if (duplicates[index] == selection[index]) {
      proclist[k] = in_rvous[m];
      ownflags[k] = in_rvous[m+1];
      k++;
    }
    duplicates[index]++;
    m += 3;
  }

  // clean up

  memory->destroy(duplicates);
  memory->destroy(selection);

  // flag = 2: new outbuf

  flag = 2;
  outbuf = (char *) ownflags;
  return nown;
}

/* ----------------------------------------------------------------------
   spread values for N datums per surf from in vector to out vector
   type = 0/1 = INT or DOUBLE datums
   comm done via rendezvous algorithm
   in = nlocal surfs (for my owned cells), only unique subset is used
   out = mylines or mytris
   only called for explicit distributed surfs
   called by spread_custom() and surf_collide.cpp
------------------------------------------------------------------------- */

void Surf::spread_local2own(int n, int type, void *in, void *out)
{
  int i,j,k,idata;
  int *iinput,*ibuf;
  double *dinput,*dbuf;

  // allocate memory for rvous input

  int *proclist;
  memory->create(proclist,nunique,"spread/local2own:proclist");

  ibuf = NULL;
  dbuf = NULL;
  if (type == INT) {
    iinput = (int *) in;
    memory->create(ibuf,(n+1)*nunique,"spread/local2own:ibuf");
  } else if (type == DOUBLE) {
    dinput = (double *) in;
    memory->create(dbuf,(n+1)*nunique,"spread/local2own:dbuf");
  }

  int isurf,index;
  surfint surfID;

  k = 0;
  for (i = 0; i < nunique; i++) {
    isurf = unique[i];
    surfID = uniqueID[i];
    proclist[i] = (surfID-1) % nprocs;
    index = (surfID-1) / nprocs;

    if (type == INT) {
      ibuf[k++] = index;
      if (n == 1)
	ibuf[k++] = iinput[isurf];
      else {
	idata = index * n;
	for (j = 0; j < n; j++)
	  ibuf[k++] = iinput[idata++];
      }
    } else {
      dbuf[k++] = index;
      if (n == 1)
	dbuf[k++] = dinput[isurf];
      else {
	idata = index * n;
	for (j = 0; j < n; j++)
	  dbuf[k++] = dinput[idata++];
      }
    }
  }

  // perform rendezvous operation

  spread_type = type;
  spread_size = n;
  spread_data = out;

  char *in_rvous;
  int inbytes;
  if (type == INT) {
    in_rvous = (char *) ibuf;
    inbytes = (n+1) * sizeof(int);
  } else if (type == DOUBLE) {
    in_rvous = (char *) dbuf;
    inbytes = (n+1) * sizeof(double);
  }

  char *buf;
  int nreturn = comm->rendezvous(1,nunique,in_rvous,inbytes,
				 0,proclist,rendezvous_local2own,
				 0,buf,0,(void *) this);

  memory->destroy(proclist);
  memory->destroy(ibuf);
  memory->destroy(dbuf);
}

/* ----------------------------------------------------------------------
   callback from spread_local2own rendezvous operation
   receive values for surf elements I own
------------------------------------------------------------------------- */

int Surf::rendezvous_local2own(int n, char *inbuf,
			       int &flag, int *&proclist, char *&outbuf,
			       void *ptr)
{
  int j,m,idata;
  int *ibuf,*iout;
  double *dbuf,*dout;

  Surf *sptr = (Surf *) ptr;
  int type = sptr->spread_type;
  int size = sptr->spread_size;

  if (type == INT) {
    ibuf = (int *) inbuf;
    iout = (int *) sptr->spread_data;
  } else if (type == DOUBLE) {
    dbuf = (double *) inbuf;
    dout = (double *) sptr->spread_data;
  }

  // loop over received datums, copy values to out

  int index;

  m = 0;
  if (type == INT) {
    for (int i = 0; i < n; i++) {
      index = ibuf[m++];
      if (size == 1)
	iout[index] = ibuf[m++];
      else {
	idata = index * size;
	for (j = 0; j < size; j++)
	  iout[idata++] = ibuf[m++];
      }
    }
  } else if (type == DOUBLE) {
    for (int i = 0; i < n; i++) {
      index = static_cast<int> (dbuf[m++]);
      if (size == 1)
	dout[index] = dbuf[m++];
      else {
	idata = index * size;
	for (j = 0; j < size; j++)
	  dout[idata++] = dbuf[m++];
      }
    }
  }

  // flag = 0: no second comm needed in rendezvous

  flag = 0;
  return 0;
}
