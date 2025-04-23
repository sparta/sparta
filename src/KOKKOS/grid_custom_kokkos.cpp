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
#include "grid_kokkos.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};                               // several files

// per grid cell custom attributes

/* ----------------------------------------------------------------------
   add a custom attribute with name
   assumes name does not already exist, else error
   type = 0/1 for int/double
   size = 0 for vector, size > 0 for array with size columns
   return index of its location;
------------------------------------------------------------------------- */

int GridKokkos::add_custom(char *name, int type, int size)
{
  // modifies eivec,eiarray,edvec,edarray on either host or device, probably device since host isn't modified. May just want to use host
  // modifies ewhich on host, sync to device here since it is never modified on the device

  // force resize on host

  k_eivec.modify_host();
  k_eiarray.modify_host();
  k_edvec.modify_host();
  k_edarray.modify_host();

  k_ewhich.modify_host();
  k_eicol.modify_host();
  k_edcol.modify_host();

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
    memory->grow(esize,ncustom,"grid:esize");
    memoryKK->grow_kokkos(k_ewhich,ewhich,ncustom,"grid:ewhich");
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
        memory->srealloc(eivec,ncustom_ivec*sizeof(int *),"grid:eivec");
      eivec[ncustom_ivec-1] = NULL;
      auto h_eivec = k_eivec.h_view;
      k_eivec.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_ivec);
      memory->grow(icustom_ivec,ncustom_ivec,"grid:icustom_ivec");
      icustom_ivec[ncustom_ivec-1] = index;
    } else {
      ewhich[index] = ncustom_iarray++;
      eiarray = (int ***)
        memory->srealloc(eiarray,ncustom_iarray*sizeof(int **),
                         "grid:eiarray");
      eiarray[ncustom_iarray-1] = NULL;
      auto h_eiarray = k_eiarray.h_view;
      k_eiarray.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_iarray);
      memory->grow(icustom_iarray,ncustom_iarray,"grid:icustom_iarray");
      icustom_iarray[ncustom_iarray-1] = index;
      memoryKK->grow_kokkos(k_eicol,eicol,ncustom_iarray,"grid:eicol");
      eicol[ncustom_iarray-1] = size;
    }
  } else if (type == DOUBLE) {
    if (size == 0) {
      ewhich[index] = ncustom_dvec++;
      edvec = (double **)
        memory->srealloc(edvec,ncustom_dvec*sizeof(double *),"grid:edvec");
      edvec[ncustom_dvec-1] = NULL;
      auto h_edvec = k_edvec.h_view;
      k_edvec.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_dvec);
      memory->grow(icustom_dvec,ncustom_dvec,"grid:icustom_dvec");
      icustom_dvec[ncustom_dvec-1] = index;
    } else {
      ewhich[index] = ncustom_darray++;
      edarray = (double ***)
        memory->srealloc(edarray,ncustom_darray*sizeof(double **),
                         "grid:edarray");
      edarray[ncustom_darray-1] = NULL;
      auto h_edarray = k_edarray.h_view;
      k_edarray.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_darray);
      memory->grow(icustom_darray,ncustom_darray,"grid:icustom_darray");
      icustom_darray[ncustom_darray-1] = index;
      memoryKK->grow_kokkos(k_edcol,edcol,ncustom_darray,"grid:edcol");
      edcol[ncustom_darray-1] = size;
    }
  }

  // ewhich,eicol,edcol never modified on the device, so sync here

  k_ewhich.modify_host();
  k_ewhich.sync_device();

  k_eicol.modify_host();
  k_eicol.sync_device();

  k_edcol.modify_host();
  k_edcol.sync_device();

  allocate_custom(index);

  return index;
}

/* ----------------------------------------------------------------------
   allocate ONE custom per-surf vector/array associated with new index
   via memory->create() to current size nown
   set all values to 0 via memset()
------------------------------------------------------------------------- */

void GridKokkos::allocate_custom(int index)
{
  // modifies the inner part of eivec,eiarray,edvec,edarray on whatever, and the outer view on the host

  k_eivec.sync_host();
  k_eiarray.sync_host();
  k_edvec.sync_host();
  k_edarray.sync_host();

  int n = maxcell;

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      int *ivector = eivec[ewhich[index]];
      auto k_ivector = k_eivec.h_view[ewhich[index]].k_view;
      k_ivector.modify_host(); // force resize on host
      memoryKK->grow_kokkos(k_ivector,ivector,n,"surf:eivec");
      k_eivec.h_view[ewhich[index]].k_view = k_ivector;
      eivec[ewhich[index]] = ivector;
    } else {
      int **iarray = eiarray[ewhich[index]];
      auto k_iarray = k_eiarray.h_view[ewhich[index]].k_view;
      k_iarray.modify_host(); // force resize on host
      memoryKK->grow_kokkos(k_iarray,iarray,n,esize[index],"surf:eiarray");
      k_eiarray.h_view[ewhich[index]].k_view = k_iarray;
      eiarray[ewhich[index]] = iarray;
    }

  } else {
    if (esize[index] == 0) {
      double *dvector = edvec[ewhich[index]];
      auto k_dvector = k_edvec.h_view[ewhich[index]].k_view;
      k_dvector.modify_host(); // force resize on host
      memoryKK->grow_kokkos(k_dvector,dvector,n,"surf:edvec");
      k_edvec.h_view[ewhich[index]].k_view = k_dvector;
      edvec[ewhich[index]] = dvector;
    } else {
      double **darray = edarray[ewhich[index]];
      auto k_darray = k_edarray.h_view[ewhich[index]].k_view;
      k_darray.modify_host(); // force resize on host
      memoryKK->grow_kokkos(k_darray,darray,n,esize[index],"surf:edarray");
      k_edarray.h_view[ewhich[index]].k_view = k_darray;
      edarray[ewhich[index]] = darray;
    }
  }

  k_eivec.modify_host();
  k_eiarray.modify_host();
  k_edvec.modify_host();
  k_edarray.modify_host();

  k_eivec.sync_device();
  k_eiarray.sync_device();
  k_edvec.sync_device();
  k_edarray.sync_device();
}

/* ----------------------------------------------------------------------
   reallocate vector/array associated with custom attribute with index
   set new values to 0 via memset()
------------------------------------------------------------------------- */

void GridKokkos::reallocate_custom(int /*nold*/, int nnew)
{
  // modifies the inner part of eivec,eiarray,edvec,edarray on whatever, and the outer view on the host

  k_eivec.sync_host();
  k_eiarray.sync_host();
  k_edvec.sync_host();
  k_edarray.sync_host();

  for (int ic = 0; ic < ncustom; ic++) {
    if (etype[ic] == INT) {
      if (esize[ic] == 0) {
        int *ivector = eivec[ewhich[ic]];
        auto k_ivector = k_eivec.h_view[ewhich[ic]].k_view;
        k_ivector.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_ivector,ivector,nnew,"surf:eivec");
        k_eivec.h_view[ewhich[ic]].k_view = k_ivector;
        eivec[ewhich[ic]] = ivector;
      } else {
        int **iarray = eiarray[ewhich[ic]];
        auto k_iarray = k_eiarray.h_view[ewhich[ic]].k_view;
        k_iarray.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_iarray,iarray,nnew,esize[ic],"surf:eiarray");
        k_eiarray.h_view[ewhich[ic]].k_view = k_iarray;
        eiarray[ewhich[ic]] = iarray;
      }

    } else {
      if (esize[ic] == 0) {
        double *dvector = edvec[ewhich[ic]];
        auto k_dvector = k_edvec.h_view[ewhich[ic]].k_view;
        k_dvector.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_dvector,dvector,nnew,"surf:edvec");
        k_edvec.h_view[ewhich[ic]].k_view = k_dvector;
        edvec[ewhich[ic]] = dvector;
      } else {
        double **darray = edarray[ewhich[ic]];
        auto k_darray = k_edarray.h_view[ewhich[ic]].k_view;
        k_darray.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_darray,darray,nnew,esize[ic],"surf:edarray");
        k_edarray.h_view[ewhich[ic]].k_view = k_darray;
        edarray[ewhich[ic]] = darray;
      }
    }
  }

  k_eivec.modify_host();
  k_eiarray.modify_host();
  k_edvec.modify_host();
  k_edarray.modify_host();

  k_eivec.sync_device();
  k_eiarray.sync_device();
  k_edvec.sync_device();
  k_edarray.sync_device();
}

/* ----------------------------------------------------------------------
   remove a custom attribute at location index
   free memory for name and vector/array and set ptrs to NULL
   ncustom lists never shrink, but indices stored between
     the ncustom list and the dense vector/array lists must be reset
------------------------------------------------------------------------- */

void GridKokkos::remove_custom(int index)
{
  // modifies the outer host view, deletes the inner dual view

  if (!ename || !ename[index]) return;

  delete [] ename[index];
  ename[index] = NULL;

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      memoryKK->destroy_kokkos(k_eivec.h_view[ewhich[index]].k_view,eivec[ewhich[index]]);
      ncustom_ivec--;
      for (int i = ewhich[index]; i < ncustom_ivec; i++) {
        icustom_ivec[i] = icustom_ivec[i+1];
        ewhich[icustom_ivec[i]] = i;
        eivec[i] = eivec[i+1];
        k_eivec.h_view[i] = k_eivec.h_view[i+1];
      }
    } else {
      memoryKK->destroy_kokkos(k_eiarray.h_view[ewhich[index]].k_view,eiarray[ewhich[index]]);
      ncustom_iarray--;
      for (int i = ewhich[index]; i < ncustom_iarray; i++) {
        icustom_iarray[i] = icustom_iarray[i+1];
        ewhich[icustom_iarray[i]] = i;
        eiarray[i] = eiarray[i+1];
        eicol[i] = eicol[i+1];
        k_eiarray.h_view[i] = k_eiarray.h_view[i+1];
      }
    }
  } else if (etype[index] == DOUBLE) {
    if (esize[index] == 0) {
      memoryKK->destroy_kokkos(k_edvec.h_view[ewhich[index]].k_view,edvec[ewhich[index]]);
      ncustom_dvec--;
      for (int i = ewhich[index]; i < ncustom_dvec; i++) {
        icustom_dvec[i] = icustom_dvec[i+1];
        ewhich[icustom_dvec[i]] = i;
        edvec[i] = edvec[i+1];
        k_edvec.h_view[i] = k_edvec.h_view[i+1];
      }
      k_edvec.modify_host();
    } else {
      memoryKK->destroy_kokkos(k_edarray.h_view[ewhich[index]].k_view,edarray[ewhich[index]]);
      ncustom_darray--;
      for (int i = ewhich[index]; i < ncustom_darray; i++) {
        icustom_darray[i] = icustom_darray[i+1];
        ewhich[icustom_darray[i]] = i;
        edarray[i] = edarray[i+1];
        edcol[i] = edcol[i+1];
        k_edarray.h_view[i] = k_edarray.h_view[i+1];
      }
      k_edarray.modify_host();
    }
  }

  // set ncustom = 0 if custom list is now entirely empty

  int empty = 1;
  for (int i = 0; i < ncustom; i++)
    if (ename[i]) empty = 0;
  if (empty) ncustom = 0;

  k_eivec.sync_device();
  k_eiarray.sync_device();
  k_edvec.sync_device();
  k_edarray.sync_device();
}

/* ----------------------------------------------------------------------
   copy custom values from Icell to Jcell
   called whenever a grid cell is removed from this processor's list
   caller checks that Icell != Jcell
------------------------------------------------------------------------- */

void GridKokkos::copy_custom(int icell, int jcell)
{
  this->sync(Host,CUSTOM_MASK);
  Grid::copy_custom(icell,jcell);
  this->modify(Host,CUSTOM_MASK);
}

/* ----------------------------------------------------------------------
   pack a custom attributes for a single grid cell ICELL into buf
   memflag = 0/1 = no/yes to actually pack into buf, 0 = just length
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

int GridKokkos::pack_custom(int icell, char *buf, int memflag)
{
  this->sync(Host,CUSTOM_MASK);
  return Grid::pack_custom(icell,buf,memflag);
}

/* ----------------------------------------------------------------------
   unpack custom attributes for a single grid cell ICELL from buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

int GridKokkos::unpack_custom(char *buf, int icell)
{
  int n = Grid::unpack_custom(buf,icell);
  this->modify(Host,CUSTOM_MASK);
  return n;
}
