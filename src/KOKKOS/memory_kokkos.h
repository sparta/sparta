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

#ifndef SPARTA_MEMORY_KOKKOS_H
#define SPARTA_MEMORY_KOKKOS_H

#include "memory.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

typedef MemoryKokkos MemKK;

class MemoryKokkos : public Memory {
 public:
  MemoryKokkos(class SPARTA *sparta) : Memory(sparta) {}

/* ----------------------------------------------------------------------
   Kokkos versions of create/grow/destroy multi-dimensional arrays
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   create a 1d array
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type *&array,
                   int n1, const char *name)
{
  data = TYPE(name,n1);
  array = data.h_view.data();
  return data;
}

template <typename TYPE>
TYPE wrap_kokkos(TYPE &data, const typename TYPE::value_type *array,
                   int n, const char *name)
{
  data = TYPE(name,n);
  for (int i=0; i<n; i++) {
    data.h_view(i) = array[i];
  }
  return data;
}

template <typename TYPE>
TYPE unwrap_kokkos(const TYPE &data, typename TYPE::value_type *&array,
                   int n)
{
  for (int i=0; i<n; i++)
    array[i] = data.h_view(i);
  return data;
}

template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data,
                     typename TYPE::value_type *&array, int n1,
                     const char *name)
{
  data = TYPE(std::string(name),n1);
  h_data = Kokkos::create_mirror_view(data);
  array = h_data.data();
  return data;
}


template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data,
                     int n1, const char *name)
{
  data = TYPE(std::string(name),n1);
  h_data = Kokkos::create_mirror_view(data);
  return data;
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 1d array
   last dim must stay the same
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE grow_kokkos(TYPE &data, typename TYPE::value_type *&array,
                 int n1, const char *name)
{
  if (array == NULL) return create_kokkos(data,array,n1,name);

  data.resize(n1);
  array = data.h_view.data();
  return data;
}

template <typename TYPE>
void destroy_kokkos(TYPE data, typename TYPE::value_type* &array)
{
  if (array == NULL) return;

  if (!data.d_view.data()) {
    destroy(array);
    return;
  }

  data = TYPE();
  array = NULL;
}

/* ----------------------------------------------------------------------
   create a 2d array
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE destroy_kokkos(TYPE &data)
{
  /*if(data.data()!=NULL)
    free(data.data());*/
  data = TYPE();
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2,n3);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3, int n4 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*n4*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2,n3,n4);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3, int n4, int n5 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*n4*n5*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name),n1,n2,n3,n4,n5);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, int n1, int n2, int n3, int n4, int n5 , int n6 ,const char *name)
{
  /*typename TYPE::non_const_value_type* ptr = (typename TYPE::non_const_value_type*)
    malloc(n1*n2*n3*n4*n5*n6*sizeof(typename TYPE::non_const_value_type)*4);*/
  data = TYPE(std::string(name) ,n1,n2,n3,n4,n5,n6);
  return data;
}



template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data, int n1, int n2,
                     const char *name)
{
  data = TYPE(std::string(name),n1,n2);
  h_data = Kokkos::create_mirror_view(data);
  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type **&array,
                   int n1, int n2, const char *name)
{
  data = TYPE(std::string(name),n1,n2);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);

  bigint n = 0;
  for (int i = 0; i < n1; i++) {
    if(n2==0)
      array[i] = NULL;
    else
      array[i] = &data.h_view(i,0);
    n += n2;
  }
  return data;
}

template <typename TYPE, typename HTYPE>
  TYPE create_kokkos(TYPE &data, HTYPE &h_data,
                     typename TYPE::value_type **&array, int n1, int n2,
                     const char *name)
{
  data = TYPE(std::string(name),n1,n2);
  h_data = Kokkos::create_mirror_view(data);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);

  bigint n = 0;
  for (int i = 0; i < n1; i++) {
    if(n2==0)
      array[i] = NULL;
    else
      array[i] = &h_data(i,0);
    n += n2;
  }
  return data;
}

/* ----------------------------------------------------------------------
   grow or shrink 1st dim of a 2d array
   last dim must stay the same
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE grow_kokkos(TYPE &data, typename TYPE::value_type **&array,
                 int n1, int n2, const char *name)
{
  if (array == NULL) return create_kokkos(data,array,n1,n2,name);
  data.resize(n1,n2);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type**) srealloc(array,nbytes,name);

  for (int i = 0; i < n1; i++)
    if(n2==0)
      array[i] = NULL;
    else
      array[i] = &data.h_view(i,0);

  return data;
}

template <typename TYPE>
TYPE create_kokkos(TYPE &data, typename TYPE::value_type **&array,
                   int n1, const char *name)
{
  data = TYPE(std::string(name),n1);
  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) smalloc(nbytes,name);

  for (int i = 0; i < n1; i++)
    if(data.h_view.extent(1)==0)
      array[i] = NULL;
    else
      array[i] = &data.h_view(i,0);

  return data;
}

template <typename TYPE>
TYPE grow_kokkos(TYPE &data, typename TYPE::value_type **&array,
                 int n1, const char *name)
{
  if (array == NULL) return create_kokkos(data,array,n1,name);

  data.resize(n1);

  bigint nbytes = ((bigint) sizeof(typename TYPE::value_type *)) * n1;
  array = (typename TYPE::value_type **) srealloc(array,nbytes,name);

  for (int i = 0; i < n1; i++)
    if(data.h_view.extent(1)==0)
      array[i] = NULL;
    else
      array[i] = &data.h_view(i,0);

  return data;
}

/* ----------------------------------------------------------------------
   destroy a 2d array
------------------------------------------------------------------------- */

template <typename TYPE>
void destroy_kokkos(TYPE data, typename TYPE::value_type** &array)
{
  if (array == NULL) return;

  if (!data.d_view.data()) {
    destroy(array);
    return;
  }

  data = TYPE();
  sfree(array);
  array = NULL;
}

/* ----------------------------------------------------------------------
   reallocate Kokkos views without initialization
   deallocate first to reduce memory use
------------------------------------------------------------------------- */

template <typename TYPE, typename... Indices>
static void realloc_kokkos(TYPE &data, const char *name, Indices... ns)
{
  data = TYPE();
  data = TYPE(Kokkos::NoInit(std::string(name)), ns...);
}

/* ----------------------------------------------------------------------
   get memory usage of Kokkos view in bytes
------------------------------------------------------------------------- */

template <typename TYPE>
static double memory_usage(TYPE &data)
{
  return data.span() * sizeof(typename TYPE::value_type);
}

};

}

#endif

