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

#include "mpi.h"
#include "string.h"
#include "particle.h"
#include "collide.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NONE,DISCRETE,SMOOTH};            // several files
enum{INT,DOUBLE};                      // several files

// per particle custom attributes

/* ----------------------------------------------------------------------
   find custom per-atom vector/array with name
   return index if found
   return -1 if not found
------------------------------------------------------------------------- */

int Particle::find_custom(char *name)
{
  for (int i = 0; i < ncustom; i++)
    if (ename[i] && strcmp(ename[i],name) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   error checks on existence of custom vectors/arrays
------------------------------------------------------------------------- */

void Particle::error_custom()
{
  if (collide && collide->vibstyle == DISCRETE && maxvibmode > 1) {
    int index = find_custom((char *) "vibmode");
    if (index < 0)
      error->all(FLERR,"No custom particle vibmode array defined");
    if (esize[index] != maxvibmode)
      error->all(FLERR,"Custom particle vibmode array is wrong size");
  }
}

/* ----------------------------------------------------------------------
   add a custom attribute with name
   assumes name does not already exist, else error
   type = 0/1 for int/double
   size = 0 for vector, size > 0 for array with size columns
   allocate the vector or array to current maxlocal via grow_custom()
   return index of its location;
------------------------------------------------------------------------- */

int Particle::add_custom(char *name, int type, int size)
{
  int index;

  // error if name already exists

  index = find_custom(name);
  if (index >= 0)
    error->all(FLERR,"Custom particle attribute name already exists");

  // use first available NULL entry or allocate a new one

  for (index = 0; index < ncustom; index++)
    if (ename[index] == NULL) break;

  if (index == ncustom) {
    ncustom++;
    ename = (char **) memory->srealloc(ename,ncustom*sizeof(char *),
                                       "particle:ename");
    memory->grow(etype,ncustom,"particle:etype");
    memory->grow(esize,ncustom,"particle:esize");
    memory->grow(ewhich,ncustom,"particle:ewhich");
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
        memory->srealloc(eivec,ncustom_ivec*sizeof(int *),"particle:eivec");
      eivec[ncustom_ivec-1] = NULL;
      memory->grow(icustom_ivec,ncustom_ivec,"particle:icustom_ivec");
      icustom_ivec[ncustom_ivec-1] = index;
    } else {
      ewhich[index] = ncustom_iarray++;
      eiarray = (int ***)
        memory->srealloc(eiarray,ncustom_iarray*sizeof(int **),
                         "particle:eiarray");
      eiarray[ncustom_iarray-1] = NULL;
      memory->grow(icustom_iarray,ncustom_iarray,"particle:icustom_iarray");
      icustom_iarray[ncustom_iarray-1] = index;
      memory->grow(eicol,ncustom_iarray,"particle:eicol");
      eicol[ncustom_iarray-1] = size;
    }
  } else if (type == DOUBLE) {
    if (size == 0) {
      ewhich[index] = ncustom_dvec++;
      edvec = (double **)
        memory->srealloc(edvec,ncustom_dvec*sizeof(double *),"particle:edvec");
      edvec[ncustom_dvec-1] = NULL;
      memory->grow(icustom_dvec,ncustom_dvec,"particle:icustom_dvec");
      icustom_dvec[ncustom_dvec-1] = index;
    } else {
      ewhich[index] = ncustom_darray++;
      edarray = (double ***)
        memory->srealloc(edarray,ncustom_darray*sizeof(double **),
                         "particle:edarray");
      edarray[ncustom_darray-1] = NULL;
      memory->grow(icustom_darray,ncustom_darray,"particle:icustom_darray");
      icustom_darray[ncustom_darray-1] = index;
      memory->grow(edcol,ncustom_darray,"particle:edcol");
      edcol[ncustom_darray-1] = size;
    }
  }

  grow_custom(index,0,maxlocal);

  return index;
}

/* ----------------------------------------------------------------------
   grow the vector/array associated with custom attribute with index
   nold = old length, nnew = new length (typically maxlocal)
   set new values to 0 via memset()
------------------------------------------------------------------------- */

void Particle::grow_custom(int index, int nold, int nnew)
{
  if (etype[index] == INT) {
    if (esize[index] == 0) {
      int *ivector = eivec[ewhich[index]];
      memory->grow(ivector,nnew,"particle:eivec");
      if (ivector) memset(&ivector[nold],0,(nnew-nold)*sizeof(int));
      eivec[ewhich[index]] = ivector;
    } else {
      int **iarray = eiarray[ewhich[index]];
      memory->grow(iarray,nnew,esize[index],"particle:eiarray");
      if (iarray)
        memset(&iarray[nold][0],0,(nnew-nold)*esize[index]*sizeof(int));
      eiarray[ewhich[index]] = iarray;
    }

  } else {
    if (esize[index] == 0) {
      double *dvector = edvec[ewhich[index]];
      memory->grow(dvector,nnew,"particle:edvec");
      if (dvector) memset(&dvector[nold],0,(nnew-nold)*sizeof(double));
      edvec[ewhich[index]] = dvector;
    } else {
      double **darray = edarray[ewhich[index]];
      memory->grow(darray,nnew,esize[index],"particle:edarray");
      if (darray)
        memset(&darray[nold][0],0,(nnew-nold)*esize[index]*sizeof(double));
      edarray[ewhich[index]] = darray;
    }
  }
}

/* ----------------------------------------------------------------------
   remove a custom attribute at location index
   free memory for name and vector/array and set ptrs to NULL
   ncustom lists never shrink, but indices stored between
     the ncustom list and the dense vector/array lists must be reset
------------------------------------------------------------------------- */

void Particle::remove_custom(int index)
{
  if (!ename || !ename[index]) return;

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
   zero info for particle I in custom attribute vectors/arrays
   called from add_particle() when a new particle is added
------------------------------------------------------------------------- */

void Particle::zero_custom(int i)
{
  int j,m;

  // 4 flavors of vectors/arrays

  if (ncustom_ivec) {
    for (m = 0; m < ncustom_ivec; m++) eivec[m][i] = 0;
  }
  if (ncustom_iarray) {
    for (m = 0; m < ncustom_iarray; m++)
      for (j = 0; j < eicol[m]; j++) eiarray[m][i][j] = 0;
  }
  if (ncustom_dvec) {
    for (m = 0; m < ncustom_dvec; m++) edvec[m][i] = 0.0;
  }
  if (ncustom_darray) {
    for (m = 0; m < ncustom_darray; m++)
      for (j = 0; j < edcol[m]; j++) edarray[m][i][j] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   copy info for one particle in custom attribute vectors/arrays
   into location I from location J
------------------------------------------------------------------------- */

void Particle::copy_custom(int i, int j)
{
  int m;

  // caller does not always check this
  // shouldn't be a problem, but valgrind can complain if memcpy to self
  // oddly memcpy(&particles[i],&particles[j],sizeof(OnePart)) seems OK

  if (i == j) return;

  // 4 flavors of vectors/arrays

  if (ncustom_ivec) {
    for (m = 0; m < ncustom_ivec; m++) eivec[m][i] = eivec[m][j];
  }
  if (ncustom_iarray) {
    for (m = 0; m < ncustom_iarray; m++)
      memcpy(eiarray[m][i],eiarray[m][j],eicol[m]*sizeof(int));
  }
  if (ncustom_dvec) {
    for (m = 0; m < ncustom_dvec; m++) edvec[m][i] = edvec[m][j];
  }
  if (ncustom_darray) {
    for (m = 0; m < ncustom_darray; m++)
      memcpy(edarray[m][i],edarray[m][j],edcol[m]*sizeof(double));
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes custom attribute definition info to restart file
------------------------------------------------------------------------- */

void Particle::write_restart_custom(FILE *fp)
{
  int m,index;

  // nactive = # of ncustom that have active vectors/arrays

  int nactive = 0;
  for (int i = 0; i < ncustom; i++)
    if (ename[i]) nactive++;

  fwrite(&nactive,sizeof(int),1,fp);

  // must write custom info in same order
  //   the per-particle custom values will be written into file
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

void Particle::read_restart_custom(FILE *fp)
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
   return size of all custom attributes in bytes for one particle
   used by callers to allocate buffer memory for particles
   assume integer attributes can be put at start of buffer
   only alignment needed is between integers and doubles
------------------------------------------------------------------------- */

int Particle::sizeof_custom()
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
   pack a custom attributes for a single particle N into buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

void Particle::pack_custom(int n, char *buf)
{
  int i;
  char *ptr = buf;

  if (ncustom_ivec) {
    for (i = 0; i < ncustom_ivec; i++) {
      memcpy(ptr,&eivec[i][n],sizeof(int));
      ptr += sizeof(int);
    }
  }
  if (ncustom_iarray) {
    for (i = 0; i < ncustom_iarray; i++) {
      memcpy(ptr,eiarray[i][n],eicol[i]*sizeof(int));
      ptr += eicol[i]*sizeof(int);
    }
  }

  ptr = ROUNDUP(ptr);

  if (ncustom_dvec) {
    for (i = 0; i < ncustom_dvec; i++) {
      memcpy(ptr,&edvec[i][n],sizeof(double));
      ptr += sizeof(double);
    }
  }
  if (ncustom_darray) {
    for (i = 0; i < ncustom_darray; i++) {
      memcpy(ptr,edarray[i][n],edcol[i]*sizeof(double));
      ptr += edcol[i]*sizeof(double);
    }
  }
}

/* ----------------------------------------------------------------------
   unpack custom attributes for a single particle N from buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

void Particle::unpack_custom(char *buf, int n)
{
  int i;
  char *ptr = buf;

  if (ncustom_ivec) {
    for (i = 0; i < ncustom_ivec; i++) {
      memcpy(&eivec[i][n],ptr,sizeof(int));
      ptr += sizeof(int);
    }
  }
  if (ncustom_iarray) {
    for (i = 0; i < ncustom_iarray; i++) {
      memcpy(eiarray[i][n],ptr,eicol[i]*sizeof(int));
      ptr += eicol[i]*sizeof(int);
    }
  }

  ptr = ROUNDUP(ptr);

  if (ncustom_dvec) {
    for (i = 0; i < ncustom_dvec; i++) {
      memcpy(&edvec[i][n],ptr,sizeof(double));
      ptr += sizeof(double);
    }
  }
  if (ncustom_darray) {
    for (i = 0; i < ncustom_darray; i++) {
      memcpy(edarray[i][n],ptr,edcol[i]*sizeof(double));
      ptr += edcol[i]*sizeof(double);
    }
  }
}
