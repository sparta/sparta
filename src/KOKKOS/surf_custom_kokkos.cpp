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
#include "kokkos.h"
#include "surf_kokkos.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{INT,DOUBLE};                      // several files

/* ----------------------------------------------------------------------
   add a custom attribute with name
   assumes name does not already exist, else error
   type = 0/1 for int/double
   size = 0 for vector, size > 0 for array with size columns
   return index of its location
------------------------------------------------------------------------- */

int SurfKokkos::add_custom(char *name, int type, int size)
{
  // modifies eivec,eiarray,edvec,edarray on either host or device, probably device since host isn't modified. May just want to use host
  // modifies ewhich on host, sync to device here since it is never modified on the device

  // force resize on host

  k_eivec.modify_host();
  k_eiarray.modify_host();
  k_edvec.modify_host();
  k_edarray.modify_host();

  k_eivec_local.modify_host();
  k_eiarray_local.modify_host();
  k_edvec_local.modify_host();
  k_edarray_local.modify_host();

  k_ewhich.modify_host();
  k_eicol.modify_host();
  k_edcol.modify_host();

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
    memoryKK->grow_kokkos(k_ewhich,ewhich,ncustom,"surf:ewhich");
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
      eivec[ncustom_ivec-1] = NULL;
      auto h_eivec = k_eivec.h_view;
      k_eivec.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_ivec);

      eivec_local = (int **)
        memory->srealloc(eivec_local,ncustom_ivec*sizeof(int *),"surf:eivec_local");
      eivec_local[ncustom_ivec-1] = NULL;
      auto h_eivec_local = k_eivec_local.h_view;
      k_eivec_local.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_ivec);

      memory->grow(icustom_ivec,ncustom_ivec,"surf:icustom_ivec");
      icustom_ivec[ncustom_ivec-1] = index;
    } else {
      ewhich[index] = ncustom_iarray++;

      eiarray = (int ***)
        memory->srealloc(eiarray,ncustom_iarray*sizeof(int **),
                         "surf:eiarray");
      eiarray[ncustom_iarray-1] = NULL;
      auto h_eiarray = k_eiarray.h_view;
      k_eiarray.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_iarray);

      eiarray_local = (int ***)
        memory->srealloc(eiarray_local,ncustom_iarray*sizeof(int **),
                         "surf:eiarray_local");
      eiarray_local[ncustom_iarray-1] = NULL;
      auto h_eiarray_local = k_eiarray_local.h_view;
      k_eiarray_local.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_iarray);

      memory->grow(icustom_iarray,ncustom_iarray,"surf:icustom_iarray");
      icustom_iarray[ncustom_iarray-1] = index;
      memoryKK->grow_kokkos(k_eicol,eicol,ncustom_iarray,"surf:eicol");
      eicol[ncustom_iarray-1] = size;
    }
  } else if (type == DOUBLE) {
    if (size == 0) {
      ewhich[index] = ncustom_dvec++;

      edvec = (double **)
        memory->srealloc(edvec,ncustom_dvec*sizeof(double *),"surf:edvec");
      edvec[ncustom_dvec-1] = NULL;
      auto h_edvec = k_edvec.h_view;
      k_edvec.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_dvec);

      edvec_local = (double **)
        memory->srealloc(edvec_local,ncustom_dvec*sizeof(double *),"surf:edvec_local");
      edvec_local[ncustom_dvec-1] = NULL;
      auto h_edvec_local = k_edvec_local.h_view;
      k_edvec_local.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_dvec);

      memory->grow(icustom_dvec,ncustom_dvec,"surf:icustom_dvec");
      icustom_dvec[ncustom_dvec-1] = index;
    } else {
      ewhich[index] = ncustom_darray++;

      edarray = (double ***)
        memory->srealloc(edarray,ncustom_darray*sizeof(double **),
                         "surf:edarray");
      edarray[ncustom_darray-1] = NULL;
      auto h_edarray = k_edarray.h_view;
      k_edarray.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_darray);

      edarray_local = (double ***)
        memory->srealloc(edarray_local,ncustom_darray*sizeof(double **),
                         "surf:edarray_local");
      edarray_local[ncustom_darray-1] = NULL;
      auto h_edarray_local = k_edarray_local.h_view;
      k_edarray_local.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_darray);

      memory->grow(icustom_darray,ncustom_darray,"surf:icustom_darray");
      icustom_darray[ncustom_darray-1] = index;
      memoryKK->grow_kokkos(k_edcol,edcol,ncustom_darray,"surf:edcol");
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

void SurfKokkos::allocate_custom(int index)
{
  // modifies the inner part of eivec,eiarray,edvec,edarray on whatever, and the outer view on the host

  k_eivec.sync_host();
  k_eiarray.sync_host();
  k_edvec.sync_host();
  k_edarray.sync_host();

  int n = nown;

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

  size_custom_local[index] = 0;

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
   reallocate ALL custom per-surf vectors/arrays to current nown size
   via memory->grow() to grow (or shrink) nown versus previous size_custom
   if adding values beyond old size, initialize to 0 via memset()
------------------------------------------------------------------------- */

void SurfKokkos::reallocate_custom()
{
  // modifies the inner part of eivec,eiarray,edvec,edarray on whatever, and the outer view on the host

  k_eivec.sync_host();
  k_eiarray.sync_host();
  k_edvec.sync_host();
  k_edarray.sync_host();

  int nold = size_custom;
  int nnew = nown;
  if (nnew == nold) return;

  for (int index = 0; index < ncustom; index++) {
    if (ename[index] == NULL) continue;

    if (etype[index] == INT) {
      if (esize[index] == 0) {
        int *ivector = eivec[ewhich[index]];
        auto k_ivector = k_eivec.h_view[ewhich[index]].k_view;
        k_ivector.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_ivector,ivector,nnew,"surf:eivec");
        k_eivec.h_view[ewhich[index]].k_view = k_ivector;
        eivec[ewhich[index]] = ivector;
      } else {
        int **iarray = eiarray[ewhich[index]];
        auto k_iarray = k_eiarray.h_view[ewhich[index]].k_view;
        k_iarray.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_iarray,iarray,nnew,esize[index],"surf:eiarray");
        k_eiarray.h_view[ewhich[index]].k_view = k_iarray;
        eiarray[ewhich[index]] = iarray;
      }

    } else {
      if (esize[index] == 0) {
        double *dvector = edvec[ewhich[index]];
        auto k_dvector = k_edvec.h_view[ewhich[index]].k_view;
        k_dvector.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_dvector,dvector,nnew,"surf:edvec");
        k_edvec.h_view[ewhich[index]].k_view = k_dvector;
        edvec[ewhich[index]] = dvector;
      } else {
        double **darray = edarray[ewhich[index]];
        auto k_darray = k_edarray.h_view[ewhich[index]].k_view;
        k_darray.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_darray,darray,nnew,esize[index],"surf:edarray");
        k_edarray.h_view[ewhich[index]].k_view = k_darray;
        edarray[ewhich[index]] = darray;
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
   free memory for name, set name ptr to NULL to indicate removed
   free vector/array or per-surf data
   ncustom does not shrink, but ncustom i/d vec/array lists do shrink
   cross indices bewteen ewhich and i/d vec/array lists must be reset
------------------------------------------------------------------------- */

void SurfKokkos::remove_custom(int index)
{
  // modifies the outer host view, deletes the inner dual view

  if (!ename || !ename[index]) return;

  delete [] ename[index];
  ename[index] = NULL;

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      memoryKK->destroy_kokkos(k_eivec.h_view[ewhich[index]].k_view,eivec[ewhich[index]]);
      memoryKK->destroy_kokkos(k_eivec_local.h_view[ewhich[index]].k_view,eivec_local[ewhich[index]]);
      ncustom_ivec--;
      for (int i = ewhich[index]; i < ncustom_ivec; i++) {
        icustom_ivec[i] = icustom_ivec[i+1];
        ewhich[icustom_ivec[i]] = i;
        eivec[i] = eivec[i+1];
        eivec_local[i] = eivec_local[i+1];
        k_eivec.h_view[i] = k_eivec.h_view[i+1];
        k_eivec_local.h_view[i] = k_eivec_local.h_view[i+1];
      }
      k_eivec.modify_host();
    } else {
      memoryKK->destroy_kokkos(k_eiarray.h_view[ewhich[index]].k_view,eiarray[ewhich[index]]);
      memoryKK->destroy_kokkos(k_eiarray_local.h_view[ewhich[index]].k_view,eiarray_local[ewhich[index]]);
      ncustom_iarray--;
      for (int i = ewhich[index]; i < ncustom_iarray; i++) {
        icustom_iarray[i] = icustom_iarray[i+1];
        ewhich[icustom_iarray[i]] = i;
        eicol[i] = eicol[i+1];
        eiarray[i] = eiarray[i+1];
        eiarray_local[i] = eiarray_local[i+1];
        k_eiarray.h_view[i] = k_eiarray.h_view[i+1];
        k_eiarray_local.h_view[i] = k_eiarray_local.h_view[i+1];
      }
      k_eiarray.modify_host();
    }
  } else if (etype[index] == DOUBLE) {
    if (esize[index] == 0) {
      memoryKK->destroy_kokkos(k_edvec.h_view[ewhich[index]].k_view,edvec[ewhich[index]]);
      ncustom_dvec--;
      for (int i = ewhich[index]; i < ncustom_dvec; i++) {
        icustom_dvec[i] = icustom_dvec[i+1];
        ewhich[icustom_dvec[i]] = i;
        edvec[i] = edvec[i+1];
        edvec_local[i] = edvec_local[i+1];
        k_edvec.h_view[i] = k_edvec.h_view[i+1];
        k_edvec_local.h_view[i] = k_edvec_local.h_view[i+1];
      }
      k_edvec.modify_host();
    } else {
      memoryKK->destroy_kokkos(k_edarray.h_view[ewhich[index]].k_view,edarray[ewhich[index]]);
      memoryKK->destroy_kokkos(k_edarray_local.h_view[ewhich[index]].k_view,edarray_local[ewhich[index]]);
      ncustom_darray--;
      for (int i = ewhich[index]; i < ncustom_darray; i++) {
        icustom_darray[i] = icustom_darray[i+1];
        ewhich[icustom_darray[i]] = i;
        edcol[i] = edcol[i+1];
        edarray[i] = edarray[i+1];
        edarray_local[i] = edarray_local[i+1];
        k_edarray.h_view[i] = k_edarray.h_view[i+1];
        k_edarray_local.h_view[i] = k_edarray_local.h_view[i+1];
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
   spread values for a custom attribute from owned to local+ghost vec/array
   reallocate custom_local vec/array if nlocal+nghost has changed
------------------------------------------------------------------------- */

void SurfKokkos::spread_custom(int index)
{
  // modifies the inner part of eivec,eiarray,edvec,edarray on whatever, and the outer view on the host


  if (etype[index] == INT) {
    if (esize[index] == 0) {
      k_eivec.sync_host();
      k_eivec_local.sync_host();

      if (nlocal+nghost != size_custom_local[index]) {

        size_custom_local[index] = nlocal + nghost;
        int *ivector_local = eivec_local[ewhich[index]];
        auto k_ivector_local = k_eivec_local.h_view[ewhich[index]].k_view;
        k_ivector_local.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_ivector_local,ivector_local,size_custom_local[index],
                              "surf/spread:eivec_local_vec");
        k_eivec_local.h_view[ewhich[index]].k_view = k_ivector_local;
        eivec_local[ewhich[index]] = ivector_local;

        k_eivec_local.modify_host();
      }

      k_eivec.h_view[ewhich[index]].k_view.sync_host();
      k_eivec_local.h_view[ewhich[index]].k_view.sync_host();

      spread_own2local(1,INT,eivec[ewhich[index]],
                       eivec_local[ewhich[index]]);

      k_eivec_local.h_view[ewhich[index]].k_view.modify_host();

    } else if (esize[index]) {
      k_eiarray.sync_host();
      k_eiarray_local.sync_host();

      if (nlocal+nghost != size_custom_local[index]) {

        size_custom_local[index] = nlocal + nghost;
        int **iarray_local = eiarray_local[ewhich[index]];
        auto k_iarray_local = k_eiarray_local.h_view[ewhich[index]].k_view;
        k_iarray_local.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_iarray_local,iarray_local,size_custom_local[index],
                              esize[index],"surf/spread:eiarray_local_array");
        k_eiarray_local.h_view[ewhich[index]].k_view = k_iarray_local;
        eiarray_local[ewhich[index]] = iarray_local;

        k_eiarray_local.modify_host();
      }

      k_eiarray.h_view[ewhich[index]].k_view.sync_host();
      k_eiarray_local.h_view[ewhich[index]].k_view.sync_host();

      int *in,*out;
      if (nown == 0) in = NULL;
      else in = &eiarray[ewhich[index]][0][0];
      if (size_custom_local[index] == 0) out = NULL;
      else out = &eiarray_local[ewhich[index]][0][0];
      spread_own2local(esize[index],INT,in,out);

      k_eiarray_local.h_view[ewhich[index]].k_view.modify_host();
    }

  } else if (etype[index] == DOUBLE) {
    if (esize[index] == 0) {
      k_edvec.sync_host();
      k_edvec_local.sync_host();

      if (nlocal+nghost != size_custom_local[index]) {

        size_custom_local[index] = nlocal + nghost;
        double *dvector_local = edvec_local[ewhich[index]];
        auto k_dvector_local = k_edvec_local.h_view[ewhich[index]].k_view;
        k_dvector_local.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_dvector_local,dvector_local,size_custom_local[index],
                              "surf/spread:edvec_local_vec");
        k_edvec_local.h_view[ewhich[index]].k_view = k_dvector_local;
        edvec_local[ewhich[index]] = dvector_local;

        k_edvec_local.modify_host();
      }

      k_edvec.h_view[ewhich[index]].k_view.sync_host();
      k_edvec_local.h_view[ewhich[index]].k_view.sync_host();

      spread_own2local(1,DOUBLE,edvec[ewhich[index]],
                       edvec_local[ewhich[index]]);

      k_edvec_local.h_view[ewhich[index]].k_view.modify_host();

    } else if (esize[index]) {
      k_edarray.sync_host();
      k_edarray_local.sync_host();

      if (nlocal+nghost != size_custom_local[index]) {

        size_custom_local[index] = nlocal + nghost;
        double **darray_local = edarray_local[ewhich[index]];
        auto k_darray_local = k_edarray_local.h_view[ewhich[index]].k_view;
        k_darray_local.modify_host(); // force resize on host
        memoryKK->grow_kokkos(k_darray_local,darray_local,size_custom_local[index],
                              esize[index],"surf/spread:edarray_local_array");
        k_edarray_local.h_view[ewhich[index]].k_view = k_darray_local;
        edarray_local[ewhich[index]] = darray_local;

        k_edarray_local.modify_host();
      }

      k_edarray.h_view[ewhich[index]].k_view.sync_host();
      k_edarray_local.h_view[ewhich[index]].k_view.sync_host();

      double *in,*out;
      if (nown == 0) in = NULL;
      else in = &edarray[ewhich[index]][0][0];
      if (size_custom_local[index] == 0) out = NULL;
      else out = &edarray_local[ewhich[index]][0][0];
      spread_own2local(esize[index],DOUBLE,in,out);

      k_edarray_local.h_view[ewhich[index]].k_view.modify_host();
    }
  }

  estatus[index] = 1;
}

/* ----------------------------------------------------------------------
   spread values for a custom attribute from local+ghost to owned vec/array
   only unique surfs within local list are used
------------------------------------------------------------------------- */

void SurfKokkos::spread_inverse_custom(int index)
{
  this->sync(Host,CUSTOM_MASK);
  Surf::spread_inverse_custom(index);
  this->modify(Host,CUSTOM_MASK);
}

/* ----------------------------------------------------------------------
   pack custom attributes for a single surf element ISURF into buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

int SurfKokkos::pack_custom(int isurf, char *buf)
{
  this->sync(Host,CUSTOM_MASK);
  return Surf::pack_custom(isurf,buf);
}

/* ----------------------------------------------------------------------
   unpack custom attributes from a single surf buf into custom
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

int SurfKokkos::unpack_custom(char *buf, double *custom)
{
  this->sync(Host,CUSTOM_MASK);
  int n = Surf::unpack_custom(buf,custom);
  this->modify(Host,CUSTOM_MASK);
  return n;
}
