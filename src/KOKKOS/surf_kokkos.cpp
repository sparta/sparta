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

#include "ctype.h"
#include "surf_kokkos.h"
#include "style_surf_collide.h"
#include "style_surf_react.h"
#include "surf_collide.h"
#include "surf_react.h"
#include "domain.h"
#include "region.h"
#include "comm.h"
#include "geometry.h"
#include "input.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{TALLYAUTO,TALLYREDUCE,TALLYLOCAL};         // same as Update
enum{REGION_ALL,REGION_ONE,REGION_CENTER};      // same as Grid
enum{TYPE,MOLECULE,ID};
enum{LT,LE,GT,GE,EQ,NEQ,BETWEEN};
enum{INT,DOUBLE};                      // several files

#define DELTA 4
#define EPSSQ 1.0e-12
#define BIG 1.0e20
#define MAXGROUP 32

/* ---------------------------------------------------------------------- */

SurfKokkos::SurfKokkos(SPARTA *sparta) : Surf(sparta)
{


}

/* ---------------------------------------------------------------------- */

SurfKokkos::~SurfKokkos()
{
  lines = NULL;
  tris = NULL;

  mylines = NULL;
  mytris = NULL;


  // deallocate views of views in serial to prevent race condition in profiling tools

  for (int i = 0; i < k_eivec.extent(0); i++)
    k_eivec.h_view(i).k_view = decltype(k_eivec.h_view(i).k_view)();

  for (int i = 0; i < k_eiarray.extent(0); i++)
    k_eiarray.h_view(i).k_view = decltype(k_eiarray.h_view(i).k_view)();

  for (int i = 0; i < k_edvec.extent(0); i++)
    k_edvec.h_view(i).k_view = decltype(k_edvec.h_view(i).k_view)();

  for (int i = 0; i < k_edarray.extent(0); i++)
    k_edarray.h_view(i).k_view = decltype(k_edarray.h_view(i).k_view)();

  ewhich = NULL;
  eicol = NULL;
  edcol = NULL;
}

/* ---------------------------------------------------------------------- */

void SurfKokkos::wrap_kokkos()
{
  if (domain->dimension == 2) {
    if (lines != NULL && nmax > 0) {
      if (lines != k_lines.h_view.data()) {
        memoryKK->wrap_kokkos(k_lines,lines,nmax,"surf:lines");
        k_lines.modify_host();
        k_lines.sync_device();
        memory->sfree(lines);
        lines = k_lines.h_view.data();
      }
    }
    if (mylines != NULL && nown > 0) {
      if (mylines != k_mylines.h_view.data()) {
        memoryKK->wrap_kokkos(k_mylines,mylines,nown,"surf:lines");
        k_mylines.modify_host();
        k_mylines.sync_device();
        memory->sfree(mylines);
        mylines = k_mylines.h_view.data();
      }
    }
  } else {
    if (tris != NULL && nmax > 0) {
      if (tris != k_tris.h_view.data()) {
        memoryKK->wrap_kokkos(k_tris,tris,nmax,"surf:tris");
        k_tris.modify_host();
        k_tris.sync_device();
        memory->sfree(tris);
        tris = k_tris.h_view.data();
      }
    }
    if (mytris != NULL && nown > 0) {
      if (mytris != k_mytris.h_view.data()) {
        memoryKK->wrap_kokkos(k_mytris,mytris,nown,"surf:mytris");
        k_mytris.modify_host();
        k_mytris.sync_device();
        memory->sfree(mytris);
        mytris = k_mytris.h_view.data();
      }
    }
  }
}

/* ----------------------------------------------------------------------
   add a custom attribute with name
   assumes name does not already exist, except in case of restart
   type = 0/1 for int/double
   size = 0 for vector, size > 0 for array with size columns
   allocate the vector or array to current maxlocal via grow_custom()
   return index of its location;
------------------------------------------------------------------------- */

int SurfKokkos::add_custom(char *name, int type, int size)
{
  ///modifies eivec,eiarray,edvec,edarray on either host or device, probably device since host isn't modified. May just want to use host
  ///modifies ewhich on host, sync to device here since it is never modified on the device
  //

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
  if (index >= 0) {
    if (custom_restart_flag == NULL || custom_restart_flag[index] == 1)
      error->all(FLERR,"Custom surf attribute name already exists");
    custom_restart_flag[index] = 1;
    return index;
  }

  // use first available NULL entry or allocate a new one

  for (index = 0; index < ncustom; index++)
    if (ename[index] == NULL) break;

  if (index == ncustom) {
    ncustom++;
    ename = (char **) memory->srealloc(ename,ncustom*sizeof(char *),
                                       "surf:ename");
    memory->grow(etype,ncustom,"surf:etype");
    memory->grow(esize,ncustom,"surf:esize");
    memoryKK->grow_kokkos(k_ewhich,ewhich,ncustom,"surf:ewhich");
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
      eivec[ncustom_ivec-1] = NULL;
      k_eivec.resize(ncustom_ivec);
      memory->grow(icustom_ivec,ncustom_ivec,"surf:icustom_ivec");
      icustom_ivec[ncustom_ivec-1] = index;
    } else {
      ewhich[index] = ncustom_iarray++;
      eiarray = (int ***)
        memory->srealloc(eiarray,ncustom_iarray*sizeof(int **),
                         "surf:eiarray");
      eiarray[ncustom_iarray-1] = NULL;
      k_eiarray.resize(ncustom_iarray);
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
      k_edvec.resize(ncustom_dvec);
      memory->grow(icustom_dvec,ncustom_dvec,"surf:icustom_dvec");
      icustom_dvec[ncustom_dvec-1] = index;
    } else {
      ewhich[index] = ncustom_darray++;
      edarray = (double ***)
        memory->srealloc(edarray,ncustom_darray*sizeof(double **),
                         "surf:edarray");
      edarray[ncustom_darray-1] = NULL;
      k_edarray.resize(ncustom_darray);
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

  allocate_custom(index,nlocal);

  return index;
}

/* ----------------------------------------------------------------------
   allocate vector/array associated with custom attribute with index
   set new values to 0 via memset()
------------------------------------------------------------------------- */

void SurfKokkos::allocate_custom(int index, int n)
{
  // modifies the inner part of eivec,eiarray,edvec,edarray on whatever, and the outer view on the host

  k_eivec.sync_host();
  k_eiarray.sync_host();
  k_edvec.sync_host();
  k_edarray.sync_host();

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
   remove a custom attribute at location index
   free memory for name and vector/array and set ptrs to NULL
   ncustom lists never shrink, but indices stored between
     the ncustom list and the dense vector/array lists must be reset
------------------------------------------------------------------------- */

void SurfKokkos::remove_custom(int index)
{
  // modifies the outer host view, deletes the inner dual view
  //
  delete [] ename[index];
  ename[index] = NULL;

  if (etype[index] == INT) {
    if (esize[index] == 0) {
      ncustom_ivec--;
      for (int i = ewhich[index]; i < ncustom_ivec; i++) {
        icustom_ivec[i] = icustom_ivec[i+1];
        ewhich[icustom_ivec[i]] = i;
        eivec[i] = eivec[i+1];
        k_eivec.h_view[i] = k_eivec.h_view[i+1];
      }
    } else {
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
      ncustom_dvec--;
      for (int i = ewhich[index]; i < ncustom_dvec; i++) {
        icustom_dvec[i] = icustom_dvec[i+1];
        ewhich[icustom_dvec[i]] = i;
        edvec[i] = edvec[i+1];
        k_edvec.h_view[i] = k_edvec.h_view[i+1];
      }
      k_edvec.modify_host();
    } else {
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
   grow surface data structures
------------------------------------------------------------------------- */
void SurfKokkos::grow(int old)
{
  if (sparta->kokkos->prewrap) {
    Surf::grow(old);
  } else {
    SurfKokkos* surf_kk = (SurfKokkos*) surf;

    if (domain->dimension == 2) {
      if (lines == NULL)
          surf_kk->k_lines = tdual_line_1d("surf:lines",nmax);
      else {
        surf_kk->sync(Host,LINE_MASK);
        surf_kk->modify(Host,LINE_MASK); // force resize on host
        surf_kk->k_lines.resize(nmax);
      }
      lines = surf_kk->k_lines.h_view.data();
    } else {
      if (tris == NULL)
          surf_kk->k_tris = tdual_tri_1d("surf:tris",nmax);
      else {
        surf_kk->sync(Host,TRI_MASK);
        surf_kk->modify(Host,TRI_MASK); // force resize on host
        surf_kk->k_tris.resize(nmax);
      }
      tris = surf_kk->k_tris.h_view.data();
    }
  }
}


/* ----------------------------------------------------------------------
   grow surface data structures
------------------------------------------------------------------------- */
void SurfKokkos::grow_own(int old)
{
  if (sparta->kokkos->prewrap) {
    Surf::grow_own(old);
  } else {
    SurfKokkos* surf_kk = (SurfKokkos*) surf;

    if (domain->dimension == 2) {
      if (mylines == NULL)
          surf_kk->k_mylines = tdual_line_1d("surf:mylines",nown);
      else {
        surf_kk->sync(Host,LINE_MASK);
        surf_kk->modify(Host,LINE_MASK); // force resize on host
        surf_kk->k_mylines.resize(nown);
      }
      mylines = surf_kk->k_mylines.h_view.data();
    } else {
      if (mytris == NULL)
          surf_kk->k_mytris = tdual_tri_1d("surf:mytris",nown);
      else {
        surf_kk->sync(Host,TRI_MASK);
        surf_kk->modify(Host,TRI_MASK); // force resize on host
        surf_kk->k_mytris.resize(nown);
      }
      mytris = surf_kk->k_mytris.h_view.data();
    }
  }
}

/* ---------------------------------------------------------------------- */

void SurfKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (sparta->kokkos->prewrap)
    if (space == Device)
      error->one(FLERR,"Sync Device before wrap");

  if (space == Device) {
    if (sparta->kokkos->auto_sync)
      modify(Host,mask);
    if (mask & LINE_MASK) k_lines.sync_device();
    if (mask & TRI_MASK) k_tris.sync_device();
    if (mask & SURF_CUSTOM_MASK) {
      if (ncustom) {
        if (ncustom_ivec)
          for (int i = 0; i < ncustom_ivec; i++)
            k_eivec.h_view[i].k_view.sync_device();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.h_view[i].k_view.sync_device();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.h_view[i].k_view.sync_device();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.h_view[i].k_view.sync_device();
      }
    }
  } else {
    if (mask & LINE_MASK) k_lines.sync_host();
    if (mask & TRI_MASK) k_tris.sync_host();
    if (mask & SURF_CUSTOM_MASK) {
      if (ncustom_ivec)
        for (int i = 0; i < ncustom_ivec; i++)
          k_eivec.h_view[i].k_view.sync_host();

      if (ncustom_iarray)
        for (int i = 0; i < ncustom_iarray; i++)
          k_eiarray.h_view[i].k_view.sync_host();

      if (ncustom_dvec)
        for (int i = 0; i < ncustom_dvec; i++)
          k_edvec.h_view[i].k_view.sync_host();

      if (ncustom_darray)
        for (int i = 0; i < ncustom_darray; i++)
          k_edarray.h_view[i].k_view.sync_host();
    }
  }
}

/* ---------------------------------------------------------------------- */

void SurfKokkos::modify(ExecutionSpace space, unsigned int mask)
{
  if (sparta->kokkos->prewrap)
    if (space == Device)
      error->one(FLERR,"Modify Device before wrap");

  if (space == Device) {
    if (mask & LINE_MASK) k_lines.modify_device();
    if (mask & TRI_MASK) k_tris.modify_device();
    if (mask & SURF_CUSTOM_MASK) {
      if (ncustom) {
        if (ncustom_ivec)
          for (int i = 0; i < ncustom_ivec; i++)
            k_eivec.h_view[i].k_view.modify_device();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.h_view[i].k_view.modify_device();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.h_view[i].k_view.modify_device();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.h_view[i].k_view.modify_device();
      }
    }

    if (sparta->kokkos->auto_sync)
      sync(Host,mask);
  } else {
    if (mask & LINE_MASK) k_lines.modify_host();
    if (mask & TRI_MASK) k_tris.modify_host();
    if (mask & SURF_CUSTOM_MASK) {
      if (ncustom) {
        if (ncustom_ivec)
          for (int i = 0; i < ncustom_ivec; i++)
            k_eivec.h_view[i].k_view.modify_host();

        if (ncustom_iarray)
          for (int i = 0; i < ncustom_iarray; i++)
            k_eiarray.h_view[i].k_view.modify_host();

        if (ncustom_dvec)
          for (int i = 0; i < ncustom_dvec; i++)
            k_edvec.h_view[i].k_view.modify_host();

        if (ncustom_darray)
          for (int i = 0; i < ncustom_darray; i++)
            k_edarray.h_view[i].k_view.modify_host();
      }
    }
  }
}
