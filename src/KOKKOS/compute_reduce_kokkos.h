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

#ifdef COMPUTE_CLASS

ComputeStyle(reduce/kk,ComputeReduceKokkos)

#else

#ifndef SPARTA_COMPUTE_REDUCE_KOKKOS_H
#define SPARTA_COMPUTE_REDUCE_KOKKOS_H

#include "compute_reduce.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace SPARTA_NS {

class ComputeReduceKokkos : public ComputeReduce, public KokkosBase {
 public:
  ComputeReduceKokkos(class SPARTA *, int, char **);
  ~ComputeReduceKokkos();
  void init();
  double compute_scalar();
  void compute_vector();

  // NOTE: the KokkosBase per-grid/per-particle views are intentionally left
  //   unused by this style.  compute reduce only produces a global scalar or
  //   global vector, so there is no device-resident output for a downstream
  //   Kokkos style to read.  KokkosBase is inherited so that this compute
  //   looks like every other /kk compute to a dynamic_cast<KokkosBase*>

  // ---------------------------------------------------------------------
  // DANGER, DO NOT "SIMPLIFY" THE KERNELS BELOW INTO TAGGED FUNCTORS
  //
  // ComputeReduce::~ComputeReduce() has NO "if (copymode) return;" guard
  //   (verified against src/compute_reduce.cpp: it unconditionally deletes
  //   which/argindex/flavor/ids/value2index/replace/vector/onevec/indices/
  //   owner/subsetID and destroys varparticle/vargrid/varsurf/smasks/
  //   areasurf).  Handing *this to Kokkos::parallel_for/parallel_reduce
  //   copies the object; when that copy is destroyed the base destructor
  //   runs again and double-frees all of the above.  Setting copymode=1
  //   does NOT help, because the base destructor never tests it.
  //
  // Therefore every kernel in the .cpp is a KOKKOS_LAMBDA over LOCAL copies
  //   of the views it touches.  No lambda may name a data member: KOKKOS_
  //   LAMBDA expands to [=] and naming a member would capture "this", which
  //   on a CUDA/HIP backend is a host pointer dereferenced in device code.
  //   Copy what you need into a local first.
  //
  // These helpers are public on purpose.  nvcc's extended-lambda rules
  //   forbid defining a __device__ lambda inside a member function that has
  //   private or protected access in its class, and each of them defines
  //   one.  Do not move them into the private section.
  // ---------------------------------------------------------------------

  double compute_one_kokkos(int, int);
  int setup_values(int);
  void build_include(int);
  void gather_float(DAT::t_float_1d_strided);
  void gather_int_vec(DAT::t_int_1d);
  void gather_int_array(DAT::t_int_2d_lr, int);
  double reduce_values();
  bigint count_included_kokkos();
  void grow_scratch(int);
  void sync_host_for_fallback(int);

 private:

  // gathered per-element values for the input currently being reduced, plus
  //   a 0/1 flag per element for membership in the subset (particle mixture
  //   or grid group).  both are indexed 0 <= i < nelements, where nelements
  //   is particle->nlocal for PARTICLE inputs and grid->nlocal for GRID ones

  DAT::t_float_1d d_values;
  DAT::t_int_1d d_include;
  int nelements,maxelements;

  // device copy of ComputeReduce::s2g, refreshed from the host pointer at
  //   every use so it can never disagree with, or be shorter than, the
  //   mixture map the non-Kokkos path would have read

  DAT::tdual_int_1d k_s2g;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Fix used in compute reduce not computed at compatible time

Fixes generate their values on specific timesteps.  Compute reduce is
requesting a value on a non-allowed timestep.

*/
