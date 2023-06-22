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
   assumes name does not already exist, except in case of restart
   type = 0/1 for int/double
   size = 0 for vector, size > 0 for array with size columns
   return index of its location
------------------------------------------------------------------------- */

int Surf::add_custom(char *name, int type, int size)
{
  int index;

  // if name already exists
  // just return index if a restart script and re-defining the name
  // else error

  index = find_custom(name);
  if (index >= 0) {
    if (custom_restart_flag == NULL || custom_restart_flag[index] == 1)
      error->all(FLERR,"Custom surf attribute name already exists");
    custom_restart_flag[index] = 1;
    return index;
  }

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
    memory->grow(esize,ncustom,"surf:etype");
    memory->grow(ewhich,ncustom,"surf:etype");
  }

  int n = strlen(name) + 1;
  ename[index] = new char[n];
  strcpy(ename[index],name);
  etype[index] = type;
  esize[index] = size;

  if (type == INT) {
    if (size == 0) {
      ewhich[index] = ncustom_ivec++;
      eivec = (int **)
        memory->srealloc(eivec,ncustom_ivec*sizeof(int *),"surf:eivec");
      memory->grow(icustom_ivec,ncustom_ivec,"surf:icustom_ivec");
      icustom_ivec[ncustom_ivec-1] = index;
    } else {
      ewhich[index] = ncustom_iarray++;
      eiarray = (int ***)
        memory->srealloc(eiarray,ncustom_iarray*sizeof(int **),
                         "surf:eiarray");
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
      memory->grow(icustom_dvec,ncustom_dvec,"surf:icustom_dvec");
      icustom_dvec[ncustom_dvec-1] = index;
    } else {
      ewhich[index] = ncustom_darray++;
      edarray = (double ***)
        memory->srealloc(edarray,ncustom_darray*sizeof(double **),
                         "surf:edarray");
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
  int n = nlocal;
  
  if (etype[index] == INT) {
    if (esize[index] == 0) {
      int *ivector = memory->create(eivec[ewhich[index]],n,"surf:eivec");
      if (ivector) memset(ivector,0,n*sizeof(int));
    } else {
      int **iarray = memory->create(eiarray[ewhich[index]],
                                    n,eicol[ewhich[index]],"surf:eiarray");
      if (iarray) memset(&iarray[0][0],0,n*eicol[ewhich[index]]*sizeof(int));
    }

  } else {
    if (esize[index] == 0) {
      double *dvector = memory->create(edvec[ewhich[index]],n,"surf:edvec");
      if (dvector) memset(dvector,0,n*sizeof(double));
    } else {
      double **darray = memory->create(edarray[ewhich[index]],
                                       n,edcol[ewhich[index]],"surf:eearray");
      if (darray) memset(&darray[0][0],0,n*edcol[ewhich[index]]*sizeof(double));
    }
  }
}

/* ----------------------------------------------------------------------
   reallocate ALL custom per-surf vectors/arrays to current nown size
   via memory->grow() to grow or shrink nown versus previous size_custom
   if adding storage beyond size_custom, initialize to 0 via memset()
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
   copy custom per-surf data from location I to location J in vectors/arrays
   called when adding/removing lines/triangles
   reallocflag = 1 if new line/triangle just added and being copied to
------------------------------------------------------------------------- */

void Surf::copy_custom(int i, int j, int reallocflag)
{
  if (reallocflag) reallocate_custom();
  
  for (int ic = 0; ic < ncustom; ic++) {
    if (!ename[ic]) continue;
    
    if (etype[ic] == 0) {
      if (esize[ic] == 0) {
	int *ivector = eivec[ewhich[ic]];
	ivector[j] = ivector[i];
      } else {
	int **iarray = eiarray[ewhich[ic]];
	for (int m = 0; m < esize[ic]; m++)
	  iarray[j][m] = iarray[i][m];
      }
    } else {
      if (esize[ic] == 0) {
	double *dvector = edvec[ewhich[ic]];
	dvector[j] = dvector[i];
      } else {
	double **darray = edarray[ewhich[ic]];
	for (int m = 0; m < esize[ic]; m++)
	  darray[j][m] = darray[i][m];
      }
    }
  }
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
    } else{
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
   spread values for a custom attribute from owned to local+ghost
------------------------------------------------------------------------- */

// NOTE: need to (re)allocate custom data for nlocal+nghost
// NOTE: worry about NULL array ptrs, hence no [0][0] access

void Surf::spread_custom(int index)
{
  if (etype[index] == INT) {
    if (esize[index] == 0) {
      spread_own2local(1,INT,eivec[ewhich[index]],
		       eivec_local[ewhich[index]]);
    } else if (esize[index]) {
      spread_own2local(esize[index],INT,&eiarray[ewhich[index]][0][0],
		       &eiarray_local[ewhich[index]][0][0]);
    }
  } else if (etype[index] == DOUBLE) {
    if (esize[index] == 0) {
      spread_own2local(1,DOUBLE,edvec[ewhich[index]],
		       &edvec_local[ewhich[index]]);
    } else if (esize[index]) {
      spread_own2local(esize[index],DOUBLE,&edarray[ewhich[index]][0][0],
		       &edarray_local[ewhich[index]][0][0]);
    }
  }
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
  // will be incremented as add_custom() for each nactive

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
    MPI_Bcast(&type,n,MPI_CHAR,0,world);
    if (me == 0) tmp = fread(&size,sizeof(int),1,fp);
    MPI_Bcast(&size,n,MPI_CHAR,0,world);

    // create the custom attribute

    add_custom(name,type,size);
    delete [] name;
  }

  // set flag for each newly created custom attribute to 0
  // will be reset to 1 if restart script redefines attribute with same name

  custom_restart_flag = new int[ncustom];
  for (int i = 0; i < ncustom; i++) custom_restart_flag[i] = 0;
}
