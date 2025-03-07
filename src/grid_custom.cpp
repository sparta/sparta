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
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};                               // several files

// per grid cell custom attributes

/* ----------------------------------------------------------------------
   find custom per-atom vector/array with name
   return index if found
   return -1 if not found
------------------------------------------------------------------------- */

int Grid::find_custom(char *name)
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
   return index of its location;
------------------------------------------------------------------------- */

int Grid::add_custom(char *name, int type, int size)
{
  int index;

  // error if name already exists

  index = find_custom(name);
  if (index >= 0)
    error->all(FLERR,"Custom grid attribute name already exists");

  // use first available NULL entry or allocate a new one

  for (index = 0; index < ncustom; index++)
    if (ename[index] == NULL) break;

  if (index == ncustom) {
    ncustom++;
    ename = (char **) memory->srealloc(ename,ncustom*sizeof(char *),
                                       "grid:ename");
    memory->grow(etype,ncustom,"grid:etype");
    memory->grow(esize,ncustom,"grid:etype");
    memory->grow(estatus,ncustom,"grid:estatus");
    memory->grow(ewhich,ncustom,"grid:etype");
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
        memory->srealloc(eivec,ncustom_ivec*sizeof(int *),"grid:eivec");
      memory->grow(icustom_ivec,ncustom_ivec,"grid:icustom_ivec");
      icustom_ivec[ncustom_ivec-1] = index;
    } else {
      ewhich[index] = ncustom_iarray++;
      eiarray = (int ***)
        memory->srealloc(eiarray,ncustom_iarray*sizeof(int **),
                         "grid:eiarray");
      memory->grow(icustom_iarray,ncustom_iarray,"grid:icustom_iarray");
      icustom_iarray[ncustom_iarray-1] = index;
      memory->grow(eicol,ncustom_iarray,"grid:eicol");
      eicol[ncustom_iarray-1] = size;
    }
  } else if (type == DOUBLE) {
    if (size == 0) {
      ewhich[index] = ncustom_dvec++;
      edvec = (double **)
        memory->srealloc(edvec,ncustom_dvec*sizeof(double *),"grid:edvec");
      memory->grow(icustom_dvec,ncustom_dvec,"grid:icustom_dvec");
      icustom_dvec[ncustom_dvec-1] = index;
    } else {
      ewhich[index] = ncustom_darray++;
      edarray = (double ***)
        memory->srealloc(edarray,ncustom_darray*sizeof(double **),
                         "grid:edarray");
      memory->grow(icustom_darray,ncustom_darray,"grid:icustom_darray");
      icustom_darray[ncustom_darray-1] = index;
      memory->grow(edcol,ncustom_darray,"grid:edcol");
      edcol[ncustom_darray-1] = size;
    }
  }

  allocate_custom(index);

  return index;
}

/* ----------------------------------------------------------------------
   allocate vector/array associated with custom attribute with index
   allocate to local + ghost in case any custom attributes set ghost
   set new values to 0 via memset()
------------------------------------------------------------------------- */

void Grid::allocate_custom(int index)
{
  int n = maxcell;

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      int *ivector = memory->create(eivec[ewhich[index]],n,"grid:eivec");
      if (ivector) memset(ivector,0,n*sizeof(int));
    } else {
      int **iarray = memory->create(eiarray[ewhich[index]],
                                    n,eicol[ewhich[index]],"grid:eiarray");
      if (iarray) memset(&iarray[0][0],0,n*eicol[ewhich[index]]*sizeof(int));
    }

  } else {
    if (esize[index] == 0) {
      double *dvector = memory->create(edvec[ewhich[index]],n,"grid:edvec");
      if (dvector) memset(dvector,0,n*sizeof(double));
    } else {
      double **darray = memory->create(edarray[ewhich[index]],
                                       n,edcol[ewhich[index]],"grid:eiarray");
      if (darray) memset(&darray[0][0],0,n*edcol[ewhich[index]]*sizeof(double));
    }
  }
}

/* ----------------------------------------------------------------------
   reallocate vector/array associated with custom attribute with index
   set new values to 0 via memset()
------------------------------------------------------------------------- */

void Grid::reallocate_custom(int nold, int nnew)
{
  for (int ic = 0; ic < ncustom; ic++) {
    if (etype[ic] == INT) {
      if (esize[ic] == 0) {
        int *ivector = memory->grow(eivec[ewhich[ic]],nnew,"grid:eivec");
        if (nnew > nold) memset(&ivector[nold],0,(nnew-nold)*sizeof(int));
      } else {
        int **iarray = memory->grow(eiarray[ewhich[ic]],
                                    nnew,eicol[ewhich[ic]],"grid:eiarray");
        if (nnew > nold)
          memset(&iarray[nold][0],0,
                 (nnew-nold)*eicol[ewhich[ic]]*sizeof(int));
      }

    } else {
      if (esize[ic] == 0) {
        double *dvector = memory->grow(edvec[ewhich[ic]],nnew,"grid:edvec");
        if (nnew > nold) memset(&dvector[nold],0,(nnew-nold)*sizeof(double));
      } else {
        double **darray = memory->grow(edarray[ewhich[ic]],
                                       nnew,edcol[ewhich[ic]],"grid:edarray");
        if (nnew - nold)
          memset(&darray[nold][0],0,
                 (nnew-nold)*edcol[ewhich[ic]]*sizeof(double));
      }
    }
  }
}

/* ----------------------------------------------------------------------
   remove a custom attribute at location index
   free memory for name and vector/array and set ptrs to NULL
   ncustom lists never shrink, but indices stored between
     the ncustom list and the dense vector/array lists must be reset
------------------------------------------------------------------------- */

void Grid::remove_custom(int index)
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
    } else {
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
   copy custom values from Icell to Jcell
   called whenever a grid cell is removed from this processor's list
   caller checks that Icell != Jcell
------------------------------------------------------------------------- */

void Grid::copy_custom(int icell, int jcell)
{
  int i,j,size;

  if (ncustom_ivec) {
    for (i = 0; i < ncustom_ivec; i++) {
      int *ivector = eivec[i];
      ivector[jcell] = ivector[icell];
    }
  }
  if (ncustom_iarray) {
    for (i = 0; i < ncustom_iarray; i++) {
      int **iarray = eiarray[i];
      size = esize[icustom_iarray[i]];
      for (j = 0; j < size; j++)
        iarray[jcell][j] = iarray[icell][j];
    }
  }

  if (ncustom_dvec) {
    for (i = 0; i < ncustom_dvec; i++) {
      double *dvector = edvec[i];
      dvector[jcell] = dvector[icell];
    }
  }
  if (ncustom_darray) {
    for (i = 0; i < ncustom_darray; i++) {
      double **darray = edarray[i];
      size = esize[icustom_darray[i]];
      for (j = 0; j < size; j++)
        darray[jcell][j] = darray[icell][j];
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes custom attribute definition info to restart file
------------------------------------------------------------------------- */

void Grid::write_restart_custom(FILE *fp)
{
  int m,index;

  // nactive = # of ncustom that have active vectors/arrays

  int nactive = 0;
  for (int i = 0; i < ncustom; i++)
    if (ename[i]) nactive++;

  fwrite(&nactive,sizeof(int),1,fp);

  // must write custom info in same order
  //   the per-grid custom values will be written into file
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

void Grid::read_restart_custom(FILE *fp)
{
  int tmp;

  // ncustom is 0 at time restart file is read
  // will be incremented as add_custom() for each nactive

  int nactive;
  if (me == 0) tmp = fread(&nactive,sizeof(int),1,fp);
  MPI_Bcast(&nactive,1,MPI_INT,0,world);
  if (nactive == 0) return;

  // order that custom vectors/arrays are in restart file
  //   matches order the per-particle custom values will be read from file

  int n,type,size;
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
   return size of all custom attributes in bytes for one grid cell
   used by callers to allocate buffer memory for grid cells
   assume integer attributes can be put at start of buffer
   only alignment needed is between integers and doubles
------------------------------------------------------------------------- */

int Grid::sizeof_custom()
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
   pack a custom attributes for a single grid cell ICELL into buf
   memflag = 0/1 = no/yes to actually pack into buf, 0 = just length
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

int Grid::pack_custom(int icell, char *buf, int memflag)
{
  int i;
  char *ptr = buf;

  if (ncustom_ivec) {
    for (i = 0; i < ncustom_ivec; i++) {
      if (memflag) memcpy(ptr,&eivec[i][icell],sizeof(int));
      ptr += sizeof(int);
    }
  }
  if (ncustom_iarray) {
    for (i = 0; i < ncustom_iarray; i++) {
      if (memflag) memcpy(ptr,eiarray[i][icell],eicol[i]*sizeof(int));
      ptr += eicol[i]*sizeof(int);
    }
  }

  ptr = ROUNDUP(ptr);

  if (ncustom_dvec) {
    for (i = 0; i < ncustom_dvec; i++) {
      if (memflag) memcpy(ptr,&edvec[i][icell],sizeof(double));
      ptr += sizeof(double);
    }
  }
  if (ncustom_darray) {
    for (i = 0; i < ncustom_darray; i++) {
      if (memflag) memcpy(ptr,edarray[i][icell],edcol[i]*sizeof(double));
      ptr += edcol[i]*sizeof(double);
    }
  }

  return ptr - buf;
}

/* ----------------------------------------------------------------------
   unpack custom attributes for a single grid cell ICELL from buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

int Grid::unpack_custom(char *buf, int icell)
{
  int i;
  char *ptr = buf;

  if (ncustom_ivec) {
    for (i = 0; i < ncustom_ivec; i++) {
      memcpy(&eivec[i][icell],ptr,sizeof(int));
      ptr += sizeof(int);
    }
  }
  if (ncustom_iarray) {
    for (i = 0; i < ncustom_iarray; i++) {
      memcpy(eiarray[i][icell],ptr,eicol[i]*sizeof(int));
      ptr += eicol[i]*sizeof(int);
    }
  }

  ptr = ROUNDUP(ptr);

  if (ncustom_dvec) {
    for (i = 0; i < ncustom_dvec; i++) {
      memcpy(&edvec[i][icell],ptr,sizeof(double));
      ptr += sizeof(double);
    }
  }
  if (ncustom_darray) {
    for (i = 0; i < ncustom_darray; i++) {
      memcpy(edarray[i][icell],ptr,edcol[i]*sizeof(double));
      ptr += edcol[i]*sizeof(double);
    }
  }

  return ptr - buf;
}
