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

#include "spatype.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"
#include "error.h"

#if defined(__APPLE__)
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

Memory::Memory(SPARTA *sparta) : Pointers(sparta) {}

/* ----------------------------------------------------------------------
   safe malloc
------------------------------------------------------------------------- */

void *Memory::smalloc(bigint nbytes, const char *name, int align)
{
  if (nbytes == 0) return NULL;

  void *ptr;

  if (align) {
    int retval = posix_memalign(&ptr, align, nbytes);
    if (retval) ptr = NULL;
  } else {
    ptr = malloc(nbytes);
  }

  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to allocate " BIGINT_FORMAT " bytes for array %s",
            nbytes,name);
    error->one(FLERR,str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe realloc
------------------------------------------------------------------------- */

void *Memory::srealloc(void *ptr, bigint nbytes, const char *name, int align)
{
  if (nbytes == 0) {
    destroy(ptr);
    return NULL;
  }

  if (align) {
    ptr = realloc(ptr, nbytes);
    uintptr_t offset = ((uintptr_t)(const void *)(ptr)) % align;
    if (offset) {
      void *optr = ptr;
      ptr = smalloc(nbytes, name, align);
#if defined(__APPLE__)
      memcpy(ptr, optr, MIN(nbytes,malloc_size(optr)));
#else
      memcpy(ptr, optr, MIN(nbytes,malloc_usable_size(optr)));
#endif
      free(optr);
    }
  } else
    ptr = realloc(ptr,nbytes);

  if (ptr == NULL) {
    char str[128];
    sprintf(str,"Failed to reallocate " BIGINT_FORMAT " bytes for array %s",
            nbytes,name);
    error->one(FLERR,str);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   safe free
------------------------------------------------------------------------- */

void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------
   erroneous usage of templated create/grow functions
------------------------------------------------------------------------- */

void Memory::fail(const char *name)
{
  char str[128];
  sprintf(str,"Cannot create/grow a vector/array of pointers for %s",name);
  error->one(FLERR,str);
}
