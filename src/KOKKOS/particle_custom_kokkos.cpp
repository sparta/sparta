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
#include "particle_kokkos.h"
#include "collide.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

enum{NONE,DISCRETE,SMOOTH};            // several files
enum{INT,DOUBLE};                      // several files

// per particle custom attributes

/* ----------------------------------------------------------------------
   add a custom attribute with name
   assumes name does not already exist, except in case of restart
   type = 0/1 for int/double
   size = 0 for vector, size > 0 for array with size columns
   allocate the vector or array to current maxlocal via grow_custom()
   return index of its location;
------------------------------------------------------------------------- */

int ParticleKokkos::add_custom(char *name, int type, int size)
{
  // modifies eivec,eiarray,edvec,edarray on host
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

  // if name already exists
  // just return index if a restart script and re-defining the name
  // else error

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
    memoryKK->grow_kokkos(k_ewhich,ewhich,ncustom,"particle:ewhich");
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
      k_eivec.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_ivec);
      memory->grow(icustom_ivec,ncustom_ivec,"particle:icustom_ivec");
      icustom_ivec[ncustom_ivec-1] = index;
    } else {
      ewhich[index] = ncustom_iarray++;
      eiarray = (int ***)
        memory->srealloc(eiarray,ncustom_iarray*sizeof(int **),
                         "particle:eiarray");
      eiarray[ncustom_iarray-1] = NULL;
      k_eiarray.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_iarray);
      memory->grow(icustom_iarray,ncustom_iarray,"particle:icustom_iarray");
      icustom_iarray[ncustom_iarray-1] = index;
      memoryKK->grow_kokkos(k_eicol,eicol,ncustom_iarray,"particle:eicol");
      eicol[ncustom_iarray-1] = size;
    }
  } else if (type == DOUBLE) {
    if (size == 0) {
      ewhich[index] = ncustom_dvec++;
      edvec = (double **)
        memory->srealloc(edvec,ncustom_dvec*sizeof(double *),"particle:edvec");
      edvec[ncustom_dvec-1] = NULL;
      k_edvec.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_dvec);
      memory->grow(icustom_dvec,ncustom_dvec,"particle:icustom_dvec");
      icustom_dvec[ncustom_dvec-1] = index;
    } else {
      ewhich[index] = ncustom_darray++;
      edarray = (double ***)
        memory->srealloc(edarray,ncustom_darray*sizeof(double **),
                         "particle:edarray");
      edarray[ncustom_darray-1] = NULL;
      auto h_edarray = k_edarray.view_host();
      k_edarray.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit),ncustom_darray);
      memory->grow(icustom_darray,ncustom_darray,"particle:icustom_darray");
      icustom_darray[ncustom_darray-1] = index;
      memoryKK->grow_kokkos(k_edcol,edcol,ncustom_darray,"particle:edcol");
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

  grow_custom(index,0,maxlocal);

  return index;
}

/* ----------------------------------------------------------------------
   grow the vector/array associated with custom attribute with index
   nold = old length, nnew = new length (typically maxlocal)
   nold is unused: unlike memory->grow(), the resize() inside grow_kokkos()
     value initializes, so the new values are already 0
------------------------------------------------------------------------- */

void ParticleKokkos::grow_custom(int index, int /*nold*/, int nnew)
{
  // modifies the inner part of eivec,eiarray,edvec,edarray on host, and the outer view on device

  if (sparta->kokkos->prewrap) {
    sync(Host,CUSTOM_MASK);
    modify(Host,CUSTOM_MASK);
  } else
    sync(Device,CUSTOM_MASK);

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      int *ivector = eivec[ewhich[index]];
      auto k_ivector = k_eivec.view_host()[ewhich[index]].k_view;
      memoryKK->grow_kokkos(k_ivector,ivector,nnew,"particle:ivector");
      k_eivec.view_host()[ewhich[index]].k_view = k_ivector;
      eivec[ewhich[index]] = ivector;
    } else {
      int **iarray = eiarray[ewhich[index]];
      auto k_iarray = k_eiarray.view_host()[ewhich[index]].k_view;
      memoryKK->grow_kokkos(k_iarray,iarray,nnew,esize[index],"particle:iarray");
      k_eiarray.view_host()[ewhich[index]].k_view = k_iarray;
      eiarray[ewhich[index]] = iarray;
    }

  } else {
    if (esize[index] == 0) {
      double *dvector = edvec[ewhich[index]];
      auto k_dvector = k_edvec.view_host()[ewhich[index]].k_view;
      memoryKK->grow_kokkos(k_dvector,dvector,nnew,"particle:dvector");
      k_edvec.view_host()[ewhich[index]].k_view = k_dvector;
      edvec[ewhich[index]] = dvector;
    } else {
      double **darray = edarray[ewhich[index]];
      auto k_darray = k_edarray.view_host()[ewhich[index]].k_view;
      memoryKK->grow_kokkos(k_darray,darray,nnew,esize[index],"particle:darray");
      k_edarray.view_host()[ewhich[index]].k_view = k_darray;
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
   remove a custom attribute at location index
   free memory for name and vector/array and set ptrs to NULL
   ncustom lists never shrink, but indices stored between
     the ncustom list and the dense vector/array lists must be reset
------------------------------------------------------------------------- */

void ParticleKokkos::remove_custom(int index)
{
  // modifies the outer host view, deletes the inner dual view

  if (!ename || !ename[index]) return;

  delete [] ename[index];
  ename[index] = NULL;

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      memoryKK->destroy_kokkos(k_eivec.view_host()[ewhich[index]].k_view,eivec[ewhich[index]]);
      ncustom_ivec--;
      for (int i = ewhich[index]; i < ncustom_ivec; i++) {
        icustom_ivec[i] = icustom_ivec[i+1];
        ewhich[icustom_ivec[i]] = i;
        eivec[i] = eivec[i+1];
        k_eivec.view_host()[i] = k_eivec.view_host()[i+1];
      }
    } else {
      memoryKK->destroy_kokkos(k_eiarray.view_host()[ewhich[index]].k_view,eiarray[ewhich[index]]);
      ncustom_iarray--;
      for (int i = ewhich[index]; i < ncustom_iarray; i++) {
        icustom_iarray[i] = icustom_iarray[i+1];
        ewhich[icustom_iarray[i]] = i;
        eiarray[i] = eiarray[i+1];
        eicol[i] = eicol[i+1];
        k_eiarray.view_host()[i] = k_eiarray.view_host()[i+1];
      }
    }
  } else if (etype[index] == DOUBLE) {
    if (esize[index] == 0) {
      memoryKK->destroy_kokkos(k_edvec.view_host()[ewhich[index]].k_view,edvec[ewhich[index]]);
      ncustom_dvec--;
      for (int i = ewhich[index]; i < ncustom_dvec; i++) {
        icustom_dvec[i] = icustom_dvec[i+1];
        ewhich[icustom_dvec[i]] = i;
        edvec[i] = edvec[i+1];
        k_edvec.view_host()[i] = k_edvec.view_host()[i+1];
      }
      k_edvec.modify_host();
    } else {
      memoryKK->destroy_kokkos(k_edarray.view_host()[ewhich[index]].k_view,edarray[ewhich[index]]);
      ncustom_darray--;
      for (int i = ewhich[index]; i < ncustom_darray; i++) {
        icustom_darray[i] = icustom_darray[i+1];
        ewhich[icustom_darray[i]] = i;
        edarray[i] = edarray[i+1];
        edcol[i] = edcol[i+1];
        k_edarray.view_host()[i] = k_edarray.view_host()[i+1];
      }
      k_edarray.modify_host();
    }
  }

  // set ncustom = 0 if custom list is now entirely empty

  int empty = 1;
  for (int i = 0; i < ncustom; i++)
    if (ename[i]) empty = 0;
  if (empty) ncustom = 0;

  // all four outer views may have been compacted above
  // must flag them modified on host or the syncs below are no-ops
  //   and the device keeps the stale pre-removal ordering

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
   zero the custom attributes of particles LO through HI-1, on the device
   Particle::add_particle() zeroes them for every particle it creates, so a
     device path which adds particles has to do the same.  a slot at or
     above nlocal still holds whatever the last particle there left behind,
     and without this the new particle silently inherits it
------------------------------------------------------------------------- */

void ParticleKokkos::zero_custom_kokkos(int lo, int hi)
{
  if (!ncustom) return;
  const int n = hi - lo;
  if (n <= 0) return;

  this->sync(Device,CUSTOM_MASK);

  if (ncustom_ivec) {
    auto d_ivec = k_eivec.view_device();
    const int nvec = ncustom_ivec;
    Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int m) {
      const int i = lo + m;
      for (int k = 0; k < nvec; k++)
        d_ivec[k].k_view.view_device()[i] = 0;
    });
  }

  if (ncustom_iarray) {
    auto d_iarray = k_eiarray.view_device();
    auto d_icol = k_eicol.view_device();
    const int narray = ncustom_iarray;
    Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int m) {
      const int i = lo + m;
      for (int k = 0; k < narray; k++)
        for (int c = 0; c < d_icol[k]; c++)
          d_iarray[k].k_view.view_device()(i,c) = 0;
    });
  }

  if (ncustom_dvec) {
    auto d_dvec = k_edvec.view_device();
    const int nvec = ncustom_dvec;
    Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int m) {
      const int i = lo + m;
      for (int k = 0; k < nvec; k++)
        d_dvec[k].k_view.view_device()[i] = 0.0;
    });
  }

  if (ncustom_darray) {
    auto d_darray = k_edarray.view_device();
    auto d_dcol = k_edcol.view_device();
    const int narray = ncustom_darray;
    Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int m) {
      const int i = lo + m;
      for (int k = 0; k < narray; k++)
        for (int c = 0; c < d_dcol[k]; c++)
          d_darray[k].k_view.view_device()(i,c) = 0.0;
    });
  }

  this->modify(Device,CUSTOM_MASK);
}

/* ----------------------------------------------------------------------
   zero the custom attributes of the unused slots nlocal through maxlocal-1
   for a kernel which only creates particles, zeroing the new particles once
     the kernel is done is enough.  but a kernel may also set custom
     attributes of a particle it just created, e.g. SurfCollide calls
     FixAmbipolar::update_custom_kokkos() for the products of a surf
     reaction, and zeroing afterwards would wipe those values out
   so instead zero every slot such a kernel could fill before launching it,
     which is all of nlocal to maxlocal-1: add_particle_kokkos() hands out
     slots above nlocal and sets the retry flag once maxlocal is reached
------------------------------------------------------------------------- */

void ParticleKokkos::zero_custom_kokkos()
{
  zero_custom_kokkos(nlocal,maxlocal);
}

/* ----------------------------------------------------------------------
   zero the custom attributes of particle I
   Particle::add_particle() calls this for every particle a host caller
     creates, writing through the raw eivec/edvec/... pointers.  Register
     that write, or a later sync(Device,CUSTOM_MASK) is a no-op and the
     device keeps whatever the slot last held.  Not covered by the caller
     in every case: SurfReactAdsorbKokkos inserts PS-chemistry particles on
     the host and marks only PARTICLE_MASK
   this is the host-side counterpart of zero_custom_kokkos()
------------------------------------------------------------------------- */

void ParticleKokkos::zero_custom(int i)
{
  sync(Host,CUSTOM_MASK);
  Particle::zero_custom(i);
  modify(Host,CUSTOM_MASK);
}

/* ----------------------------------------------------------------------
   copy info for one particle in custom attribute vectors/arrays
   into location I from location J
------------------------------------------------------------------------- */

void ParticleKokkos::copy_custom(int i, int j)
{
  sync(Host,CUSTOM_MASK);
  Particle::copy_custom(i,j);
  modify(Host,CUSTOM_MASK);
}

/* ----------------------------------------------------------------------
   pack a custom attributes for a single particle N into buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

void ParticleKokkos::pack_custom(int n, char *buf)
{
  sync(Host,CUSTOM_MASK);
  Particle::pack_custom(n,buf);
}

/* ----------------------------------------------------------------------
   unpack custom attributes for a single particle N from buf
   this is done in order of 4 styles of vectors/arrays, not in ncustom order
------------------------------------------------------------------------- */

void ParticleKokkos::unpack_custom(char *buf, int n)
{
  sync(Host,CUSTOM_MASK);
  Particle::unpack_custom(buf,n);
  modify(Host,CUSTOM_MASK);
}

