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

#ifndef SPARTA_REACT_H
#define SPARTA_REACT_H

#include "pointers.h"
#include "particle.h"

namespace SPARTA_NS {

class React : protected Pointers {
 public:
  char *style;
  int nlist;                 // # of reactions read from file

  int recombflag;            // 1 if any recombination reactions defined
  int recombflag_user;       // 0 if user has turned off recomb reactions
  int recomb_species;        // species of 3rd particle in recomb reaction
  int computeChemRates;      // 1 if only computing a TCE rate without
                             // actually doing reaction

  int partialEnergy;         // 1 if using rDOF model, 0 if using all energy

  enum{VIB_DOF,VIB_MICRO};
  int vibEnergyMode;         // how discrete vibrational energy couples to
                             // the TCE reaction probability:
                             // VIB_DOF (default) = instantaneous
                             //   vibrational DOF heuristic
                             //   (2*i*ln(1+1/i), newtonTvib); deviates
                             //   from the input Arrhenius rate by ~10-30%
                             //   at 10-20 kK
                             // VIB_MICRO = SHO ladder folded into the
                             //   microcanonical TCE energy factor
                             //   (jointly with the electronic ladder when
                             //   elec_energy micro); keeps the equilibrium
                             //   rate on the input Arrhenius rate

  enum{ELEC_EXCLUDE,ELEC_INCLUDE,ELEC_MICRO};
  int elecEnergyMode;        // how electronic energy couples to the TCE
                             // reaction probability (energy conservation via
                             // pre/post etotal is unaffected by the mode):
                             // ELEC_EXCLUDE = excluded from the
                             //   reaction energy; equilibrium
                             //   consistency with the input Arrhenius rates
                             //   but excited states react at same rate as
                             //   ground state
                             // ELEC_INCLUDE = added to the reaction energy
                             //   with per-state DOF from the elecfile,
                             //   requires recalibration of Arrhenius rates
                             //   and per-state DOF with electronic mode
                             //   awareness. Arrhenius rates are interpreted
                             //   as rates for ground state, and using existing
                             //   rate models will overpredict overall chemical
                             //   rates.
                             // ELEC_MICRO (default) = added to the reaction
                             //   energy with the TCE energy factor replaced by
                             //   its microcanonical average over the pair's
                             //   electronic ladder: state-sensitive while
                             //   keeping the equilibrium rate on the input
                             //   Arrhenius rate; reduces to the standard
                             //   factor when no electronic states are present
  double recomb_density;     // num density of particles in collision grid cell
  double recomb_boost;       // rate boost param for recombination reactions
  double recomb_boost_inverse;   // inverse of boost parameter
  Particle::OnePart *recomb_part3;  // ptr to 3rd particle in recomb reaction

  int copy,copymode;  // prevent deallocation of
                      //  base class when child copy is destroyed

  React(class SPARTA *, int, char **);
  React(class SPARTA *sparta) : Pointers(sparta) // needed for Kokkos
    { style = NULL; random = NULL; }
  virtual ~React();
  virtual void init() {}
  virtual int recomb_exist(int, int) = 0;
  virtual void ambi_check() = 0;
  virtual int attempt(Particle::OnePart *, Particle::OnePart *,
                      double, double, double, double, double &, int &) = 0;
  virtual char *reactionID(int) = 0;
  virtual double extract_tally(int) = 0;

  void modify_params(int, char **);

 protected:
  class RanKnuth *random;
};

}

#endif
