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

#include "stdlib.h"
#include "string.h"
#include "surf.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};                      // several files

/* ----------------------------------------------------------------------
   find custom per-atom vector/array with name
   return index if found
   return -1 if not found
------------------------------------------------------------------------- */

int Surf::find_custom(char *name)
{
  for (int i = 0; i < ncustom; i++)
    if (ename[i] && strcmp(ename[i],name) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   add a custom attribute with name
   assumes name does not already exist, else error
   type = 0/1 for int/double
   size = 0 for vector, size > 0 for array with size columns
   return index of its location
------------------------------------------------------------------------- */

int Surf::add_custom(char *name, int type, int size)
{
  int index;

  // error if name already exists

  index = find_custom(name);
  if (index >= 0)
    error->all(FLERR,"Custom surf attribute name already exists");

  // ensure all existing custom data is correct current length

  reallocate_custom();

  // assign index to new custom data
  // use first available NULL entry or allocate a new one

  for (index = 0; index < ncustom; index++)
    if (ename[index] == NULL) break;

  if (index == ncustom) {
    ncustom++;
    ename = (char **) memory->srealloc(ename,ncustom*sizeof(char *),
                                       "surf:ename");
    memory->grow(etype,ncustom,"surf:etype");
    memory->grow(esize,ncustom,"surf:esize");
    memory->grow(estatus,ncustom,"surf:estatus");
    memory->grow(ewhich,ncustom,"surf:ewhich");
    memory->grow(size_custom_local,ncustom,"surf:size_custom_local");
  }

  int n = strlen(name) + 1;
  ename[index] = new char[n];
  strcpy(ename[index],name);
  etype[index] = type;
  esize[index] = size;
  estatus[index] = 0;

  if (type == INT) {
    if (size == 0) {
      ewhich[index] = ncustom_ivec++;
      eivec = (int **)
        memory->srealloc(eivec,ncustom_ivec*sizeof(int *),"surf:eivec");
      eivec_local = (int **)
        memory->srealloc(eivec_local,ncustom_ivec*sizeof(int *),
                         "surf:eivec_local");
      memory->grow(icustom_ivec,ncustom_ivec,"surf:icustom_ivec");
      icustom_ivec[ncustom_ivec-1] = index;
    } else {
      ewhich[index] = ncustom_iarray++;
      eiarray = (int ***)
        memory->srealloc(eiarray,ncustom_iarray*sizeof(int **),
                         "surf:eiarray");
      eiarray_local = (int ***)
        memory->srealloc(eiarray_local,ncustom_iarray*sizeof(int **),
                         "surf:eiarray_local");
      memory->grow(icustom_iarray,ncustom_iarray,"surf:icustom_iarray");
      icustom_iarray[ncustom_iarray-1] = index;
      memory->grow(eicol,ncustom_iarray,"surf:eicol");
      eicol[ncustom_iarray-1] = size;
    }

  } else if (type == DOUBLE) {
    if (size == 0) {
      ewhich[index] = ncustom_dvec++;
      edvec = (double **)
        memory->srealloc(edvec,ncustom_dvec*sizeof(double *),"surf:edvec");
      edvec_local = (double **)
        memory->srealloc(edvec_local,ncustom_dvec*sizeof(double *),
                         "surf:edvec_local");
      memory->grow(icustom_dvec,ncustom_dvec,"surf:icustom_dvec");
      icustom_dvec[ncustom_dvec-1] = index;
    } else {
      ewhich[index] = ncustom_darray++;
      edarray = (double ***)
        memory->srealloc(edarray,ncustom_darray*sizeof(double **),
                         "surf:edarray");
      edarray_local = (double ***)
        memory->srealloc(edarray_local,ncustom_darray*sizeof(double **),
                         "surf:edarray_local");
      memory->grow(icustom_darray,ncustom_darray,"surf:icustom_darray");
      icustom_darray[ncustom_darray-1] = index;
      memory->grow(edcol,ncustom_darray,"surf:edcol");
      edcol[ncustom_darray-1] = size;
    }
  }

  // allocate new custom vector or array

  allocate_custom(index);

  return index;
}

/* ----------------------------------------------------------------------
   allocate ONE custom per-surf vector/array associated with new index
   via memory->create() to current size nown
   set all values to 0 via memset()
------------------------------------------------------------------------- */

void Surf::allocate_custom(int index)
{
  int n = nown;

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      int *ivector = memory->create(eivec[ewhich[index]],n,"surf:eivec");
      if (ivector) memset(ivector,0,n*sizeof(int));
      eivec_local[ewhich[index]] = NULL;
    } else {
      int **iarray = memory->create(eiarray[ewhich[index]],
                                    n,eicol[ewhich[index]],"surf:eiarray");
      if (iarray) memset(&iarray[0][0],0,n*eicol[ewhich[index]]*sizeof(int));
      eiarray_local[ewhich[index]] = NULL;
    }

  } else {
    if (esize[index] == 0) {
      double *dvector = memory->create(edvec[ewhich[index]],n,"surf:edvec");
      if (dvector) memset(dvector,0,n*sizeof(double));
      edvec_local[ewhich[index]] = NULL;
    } else {
      double **darray = memory->create(edarray[ewhich[index]],
                                       n,edcol[ewhich[index]],"surf:eearray");
      if (darray) memset(&darray[0][0],0,n*edcol[ewhich[index]]*sizeof(double));
      edarray_local[ewhich[index]] = NULL;
    }
  }

  size_custom_local[index] = 0;
}

/* ----------------------------------------------------------------------
   reallocate ALL custom per-surf vectors/arrays to current nown size
   via memory->grow() to grow (or shrink) nown versus previous size_custom
   if adding values beyond old size, initialize to 0 via memset()
------------------------------------------------------------------------- */

void Surf::reallocate_custom()
{
  int nold = size_custom;
  int nnew = nown;
  if (nnew == nold) return;

  for (int index = 0; index < ncustom; index++) {
    if (ename[index] == NULL) continue;

    if (etype[index] == INT) {
      if (esize[index] == 0) {
        int *ivector = memory->grow(eivec[ewhich[index]],nnew,"surf:eivec");
        if (nnew > nold) memset(&ivector[nold],0,(nnew-nold)*sizeof(int));
      } else {
        int **iarray = memory->grow(eiarray[ewhich[index]],
                                    nnew,eicol[ewhich[index]],"surf:eiarray");
        if (nnew > nold)
        memset(iarray[nold],0,(nnew-nold)*eicol[ewhich[index]]*sizeof(int));
      }

    } else {
      if (esize[index] == 0) {
        double *dvector = memory->grow(edvec[ewhich[index]],nnew,"surf:edvec");
        if (nnew > nold) memset(&dvector[nold],0,(nnew-nold)*sizeof(double));
      } else {
        double **darray = memory->grow(edarray[ewhich[index]],
                                       nnew,edcol[ewhich[index]],"surf:eearray");
        if (nnew > nold)
          memset(darray[nold],0,(nnew-nold)*edcol[ewhich[index]]*sizeof(double));
      }
    }
  }

  size_custom = nown;
}

/* ----------------------------------------------------------------------
   remove a custom attribute at location index
   free memory for name, set name ptr to NULL to indicate removed
   free vector/array or per-surf data
   ncustom does not shrink, but ncustom i/d vec/array lists do shrink
   cross indices bewteen ewhich and i/d vec/array lists must be reset
------------------------------------------------------------------------- */

void Surf::remove_custom(int index)
{
  delete [] ename[index];
  ename[index] = NULL;

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      memory->destroy(eivec[ewhich[index]]);
      ncustom_ivec--;
      for (int i = ewhich[index]; i < ncustom_ivec; i++) {
        icustom_ivec[i] = icustom_ivec[i+1];
        ewhich[icustom_ivec[i]] = i;
        eivec[i] = eivec[i+1];
      }
    } else {
      memory->destroy(eiarray[ewhich[index]]);
      ncustom_iarray--;
      for (int i = ewhich[index]; i < ncustom_iarray; i++) {
        icustom_iarray[i] = icustom_iarray[i+1];
        ewhich[icustom_iarray[i]] = i;
        eiarray[i] = eiarray[i+1];
        eicol[i] = eicol[i+1];
      }
    }

  } else if (etype[index] == DOUBLE) {
    if (esize[index] == 0) {
      memory->destroy(edvec[ewhich[index]]);
      ncustom_dvec--;
      for (int i = ewhich[index]; i < ncustom_dvec; i++) {
        icustom_dvec[i] = icustom_dvec[i+1];
        ewhich[icustom_dvec[i]] = i;
        edvec[i] = edvec[i+1];
      }
    } else{
      memory->destroy(edarray[ewhich[index]]);
      ncustom_darray--;
      for (int i = ewhich[index]; i < ncustom_darray; i++) {
        icustom_darray[i] = icustom_darray[i+1];
        ewhich[icustom_darray[i]] = i;
        edarray[i] = edarray[i+1];
        edcol[i] = edcol[i+1];
      }
    }
  }

  // set ncustom = 0 if custom list is now entirely empty

  int empty = 1;
  for (int i = 0; i < ncustom; i++)
    if (ename[i]) empty = 0;
  if (empty) ncustom = 0;
}

/* ----------------------------------------------------------------------
   spread values for a custom attribute from owned to local+ghost vec/array
   reallocate custom_local vec/array if nlocal+nghost has changed
------------------------------------------------------------------------- */

void Surf::spread_custom(int index)
{
  if (etype[index] == INT) {
    if (esize[index] == 0) {
      if (nlocal+nghost != size_custom_local[index]) {
        memory->destroy(eivec_local[ewhich[index]]);
        size_custom_local[index] = nlocal + nghost;
        memory->create(eivec_local[ewhich[index]],size_custom_local[index],
                       "surf/spread:eivec_local_vec");
      }

      spread_own2local(1,INT,eivec[ewhich[index]],
                       eivec_local[ewhich[index]]);

    } else if (esize[index]) {
      if (nlocal+nghost != size_custom_local[index]) {
        memory->destroy(eiarray_local[ewhich[index]]);
        size_custom_local[index] = nlocal + nghost;
        if (size_custom_local[index])
          memory->create(eiarray_local[ewhich[index]],size_custom_local[index],
                         esize[index],"surf/spread:eiarray_local_array");
        else eiarray_local[ewhich[index]] = NULL;
      }

      int *in,*out;
      if (nown == 0) in = NULL;
      else in = &eiarray[ewhich[index]][0][0];
      if (size_custom_local[index] == 0) out = NULL;
      else out = &eiarray_local[ewhich[index]][0][0];
      spread_own2local(esize[index],INT,in,out);
    }

  } else if (etype[index] == DOUBLE) {
    if (esize[index] == 0) {
      if (nlocal+nghost != size_custom_local[index]) {
        memory->destroy(edvec_local[ewhich[index]]);
        size_custom_local[index] = nlocal + nghost;
        memory->create(edvec_local[ewhich[index]],size_custom_local[index],
                       "surf/spread:edvec_local_vec");
      }

      spread_own2local(1,DOUBLE,edvec[ewhich[index]],
                       edvec_local[ewhich[index]]);

    } else if (esize[index]) {
      if (nlocal+nghost != size_custom_local[index]) {
        memory->destroy(edarray_local[ewhich[index]]);
        size_custom_local[index] = nlocal + nghost;
        if (size_custom_local[index])
          memory->create(edarray_local[ewhich[index]],size_custom_local[index],
                         esize[index],"surf/spread:edvec_local_array");
        else edarray_local[ewhich[index]] = NULL;
      }

      double *in,*out;
      if (nown == 0) in = NULL;
      else in = &edarray[ewhich[index]][0][0];
      if (size_custom_local[index] == 0) out = NULL;
      else out = &edarray_local[ewhich[index]][0][0];
      spread_own2local(esize[index],DOUBLE,in,out);
    }
  }

  estatus[index] = 1;
}

/* ----------------------------------------------------------------------
   spread values for a custom attribute from local+ghost to owned vec/array
   only unique surfs within local list are used
------------------------------------------------------------------------- */

void Surf::spread_inverse_custom(int index)
{
  if (etype[index] == INT) {
    if (esize[index] == 0) {
      spread_local2own(1,INT,eivec_local[ewhich[index]],
                       eivec[ewhich[index]]);

    } else if (esize[index]) {
      int *in,*out;
      if (nown == 0) in = NULL;
      else in = &eiarray_local[ewhich[index]][0][0];
      if (size_custom_local[index] == 0) out = NULL;
      else out = &eiarray[ewhich[index]][0][0];
      spread_local2own(esize[index],INT,in,out);
    }

  } else if (etype[index] == DOUBLE) {
    if (esize[index] == 0) {
      spread_local2own(1,DOUBLE,edvec_local[ewhich[index]],
                       &edvec[ewhich[index]]);

    } else if (esize[index]) {
      double *in,*out;
      if (nown == 0) in = NULL;
      else in = &edarray_local[ewhich[index]][0][0];
      if (size_custom_local[index] == 0) out = NULL;
      else out = &edarray[ewhich[index]][0][0];
      spread_local2own(esize[index],DOUBLE,in,out);
    }
  }

  estatus[index] = 0;
}

/* ----------------------------------------------------------------------
   extract all custom per-surf data as an array of doubles
   prepend the surfID to values for each surf
   return count of custom values per surf
   called by RemoveSurf
------------------------------------------------------------------------- */

int Surf::extract_custom(double **&cvalues)
{
  int i,j;

  // nvalues_custom = # of values per surf

  int nvalues_custom = 0;
  for (int ic = 0; ic < ncustom; ic++) {
    if (ename[ic] == NULL) continue;
    if (esize[ic] == 0) nvalues_custom++;
    else nvalues_custom += esize[ic];
  }

  memory->create(cvalues,nown,1+nvalues_custom,"surf:cvalues");

  // fill cvalues with surf ID + per-surf values

  surfint id;

  for (i = 0; i < nown; i++) {
    id = (surfint) i*nprocs + me + 1;
    cvalues[i][0] = ubuf(id).d;
  }

  int m = 1;
  for (int ic = 0; ic < ncustom; ic++) {
    if (ename[ic] == NULL) continue;

    if (etype[ic] == INT) {
      if (esize[ic] == 0) {
        int *ivector = eivec[ewhich[ic]];
        for (i = 0; i < nown; i++)
          cvalues[i][m] = ubuf(ivector[i]).d;
        m++;
      } else {
        int **iarray = eiarray[ewhich[ic]];
        int n = esize[ic];
        for (i = 0; i < nown; i++)
          for (j = 0; j < n; j++)
            cvalues[i][m+j] = ubuf(iarray[i][j]).d;
        m += esize[ic];
      }

    } else if (etype[ic] == DOUBLE) {
      if (esize[ic] == 0) {
        double *dvector = edvec[ewhich[ic]];
        for (i = 0; i < nown; i++)
          cvalues[i][m] = dvector[i];
        m++;
      } else {
        double **darray = edarray[ewhich[ic]];
        int n = esize[ic];
        for (i = 0; i < nown; i++)
          for (j = 0; j < n; j++)
            cvalues[i][m+j] = darray[i][j];
        m += esize[ic];
      }
    }
  }

  return nvalues_custom;
}

/* ----------------------------------------------------------------------
   proc 0 writes custom attribute definition info to restart file
------------------------------------------------------------------------- */

void Surf::write_restart_custom(FILE *fp)
{
  int m,index;

  // nactive = # of ncustom that have active vectors/arrays

  int nactive = 0;
  for (int i = 0; i < ncustom; i++)
    if (ename[i]) nactive++;

  fwrite(&nactive,sizeof(int),1,fp);

  // must write custom info in same order
  //   the per-surf custom values will be written into file
  // not necessarily the same as ncustom list, due to deletions & additions

  for (m = 0; m < ncustom_ivec; m++) {
    index = icustom_ivec[m];
    int n = strlen(ename[index]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(ename[index],sizeof(char),n,fp);
    fwrite(&etype[index],sizeof(int),1,fp);
    fwrite(&esize[index],sizeof(int),1,fp);
  }
  for (m = 0; m < ncustom_iarray; m++) {
    index = icustom_iarray[m];
    int n = strlen(ename[index]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(ename[index],sizeof(char),n,fp);
    fwrite(&etype[index],sizeof(int),1,fp);
    fwrite(&esize[index],sizeof(int),1,fp);
  }
  for (m = 0; m < ncustom_dvec; m++) {
    index = icustom_dvec[m];
    int n = strlen(ename[index]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(ename[index],sizeof(char),n,fp);
    fwrite(&etype[index],sizeof(int),1,fp);
    fwrite(&esize[index],sizeof(int),1,fp);
  }
  for (m = 0; m < ncustom_darray; m++) {
    index = icustom_darray[m];
    int n = strlen(ename[index]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(ename[index],sizeof(char),n,fp);
    fwrite(&etype[index],sizeof(int),1,fp);
    fwrite(&esize[index],sizeof(int),1,fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads custom attribute definition info from restart file
   bcast to other procs and all procs instantiate series of custom properties
------------------------------------------------------------------------- */

void Surf::read_restart_custom(FILE *fp)
{
  int tmp;

  // ncustom is 0 at time restart file is read
  // will be incremented via add_custom() for each nactive

  int nactive;
  if (me == 0) tmp = fread(&nactive,sizeof(int),1,fp);
  MPI_Bcast(&nactive,1,MPI_INT,0,world);
  if (nactive == 0) return;

  // order that custom vectors/arrays are in restart file
  //   matches order the per-particle custom values will be read from file

  int n,type,size,ghostflag;
  char *name;

  for (int i = 0; i < nactive; i++) {
    if (me == 0) tmp = fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    name = new char[n];
    if (me == 0) tmp = fread(name,sizeof(char),n,fp);
    MPI_Bcast(name,n,MPI_CHAR,0,world);
    if (me == 0) tmp = fread(&type,sizeof(int),1,fp);
    MPI_Bcast(&type,1,MPI_INT,0,world);
    if (me == 0) tmp = fread(&size,sizeof(int),1,fp);
    MPI_Bcast(&size,1,MPI_INT,0,world);

    // create the custom attribute

    add_custom(name,type,size);
    delete [] name;
  }
}

/* ----------------------------------------------------------------------
   return size of all custom attributes in bytes for one surface element
   used by callers to allocate buffer memory for surfs
   assume integer attributes can be put at start of buffer
   only alignment needed is between integers and doubles
------------------------------------------------------------------------- */

int Surf::sizeof_custom()
{
  int n = 0;

  n += ncustom_ivec*sizeof(int);
  if (ncustom_iarray)
    for (int i = 0; i < ncustom_iarray; i++)
      n += eicol[i]*sizeof(int);

  n = IROUNDUP(n);

  n += ncustom_dvec*sizeof(double);
  if (ncustom_darray)
    for (int i = 0; i < ncustom_darray; i++)
      n += edcol[i]*sizeof(double);

  return n;
}

/* ----------------------------------------------------------------------
   pack custom attributes for a single surf element ISURF into buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

int Surf::pack_custom(int isurf, char *buf)
{
  int i;
  char *ptr = buf;

  if (ncustom_ivec) {
    for (i = 0; i < ncustom_ivec; i++) {
      memcpy(ptr,&eivec[i][isurf],sizeof(int));
      ptr += sizeof(int);
    }
  }
  if (ncustom_iarray) {
    for (i = 0; i < ncustom_iarray; i++) {
      memcpy(ptr,eiarray[i][isurf],eicol[i]*sizeof(int));
      ptr += eicol[i]*sizeof(int);
    }
  }

  ptr = ROUNDUP(ptr);

  if (ncustom_dvec) {
    for (i = 0; i < ncustom_dvec; i++) {
      memcpy(ptr,&edvec[i][isurf],sizeof(double));
      ptr += sizeof(double);
    }
  }
  if (ncustom_darray) {
    for (i = 0; i < ncustom_darray; i++) {
      memcpy(ptr,edarray[i][isurf],edcol[i]*sizeof(double));
      ptr += edcol[i]*sizeof(double);
    }
  }

  return ptr - buf;
}

/* ----------------------------------------------------------------------
   unpack custom attributes from a single surf buf into custom
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

int Surf::unpack_custom(char *buf, double *custom)
{
  int i,j;
  char *ptr = buf;

  int ic = 0;

  if (ncustom_ivec) {
    int *ibuf = (int *) ptr;
    for (i = 0; i < ncustom_ivec; i++) {
      custom[ic++] = ibuf[i];
      ptr += sizeof(int);
    }
  }
  if (ncustom_iarray) {
    int *ibuf = (int *) ptr;
    int m = 0;
    for (i = 0; i < ncustom_iarray; i++) {
      for (j = 0; j < eicol[i]; j++)
        custom[ic++] = ibuf[m++];
      ptr += eicol[i]*sizeof(int);
    }
  }

  ptr = ROUNDUP(ptr);

  if (ncustom_dvec) {
    double *dbuf = (double *) ptr;
    for (i = 0; i < ncustom_dvec; i++) {
      custom[ic++] = dbuf[i];
      ptr += sizeof(double);
    }
  }
  if (ncustom_darray) {
    double *dbuf = (double *) ptr;
    int m = 0;
    for (i = 0; i < ncustom_darray; i++) {
      for (j = 0; j < edcol[i]; j++)
        custom[ic++] = dbuf[m++];
      ptr += edcol[i]*sizeof(double);
    }
  }

  return ptr - buf;
}
