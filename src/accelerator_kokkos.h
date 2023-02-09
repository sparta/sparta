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

#ifndef SPARTA_ACCELERATOR_KOKKOS_H
#define SPARTA_ACCELERATOR_KOKKOS_H

#ifdef SPARTA_KOKKOS

// true interface to KOKKOS
// used when KOKKOS is installed

#include "kokkos.h"
#include "sparta_masks.h"
#include "update_kokkos.h"
#include "grid_kokkos.h"
#include "particle_kokkos.h"
#include "comm_kokkos.h"
#include "domain_kokkos.h"
#include "surf_kokkos.h"
#include "modify_kokkos.h"

#else

// dummy interface to KOKKOS
// needed for compiling when KOKKOS is not installed

#include "update.h"
#include "grid.h"
#include "particle.h"
#include "comm.h"
#include "domain.h"
#include "surf.h"
#include "modify.h"

namespace SPARTA_NS {

class KokkosSPARTA {
 public:
  int kokkos_exists;
  int num_threads;
  int numa;

  KokkosSPARTA(class SPARTA *, int, char **) {kokkos_exists = 0;}
  ~KokkosSPARTA() {}
  void accelerator(int, char **) {}
};

class Kokkos {
 public:
  static void finalize() {}
};

class UpdateKokkos : public Update {
 public:
  UpdateKokkos(class SPARTA *sparta) : Update(sparta) {}
  ~UpdateKokkos() {}
};

class ParticleKokkos : public Particle {
 public:
  ParticleKokkos(class SPARTA *sparta) : Particle(sparta) {}
  ~ParticleKokkos() {}
};

class CommKokkos : public Comm {
 public:
  CommKokkos(class SPARTA *sparta) : Comm(sparta) {}
  ~CommKokkos() {}
};

class DomainKokkos : public Domain {
 public:
  DomainKokkos(class SPARTA *sparta) : Domain(sparta) {}
  ~DomainKokkos() {}
};

class GridKokkos : public Grid {
 public:
  GridKokkos(class SPARTA *sparta) : Grid(sparta) {}
  ~GridKokkos() {}
};

class SurfKokkos : public Surf {
 public:
  SurfKokkos(class SPARTA *sparta) : Surf(sparta) {}
  ~SurfKokkos() {}
};

class ModifyKokkos : public Modify {
 public:
  ModifyKokkos(class SPARTA *sparta) : Modify(sparta) {}
  ~ModifyKokkos() {}
};

}

#endif
#endif
