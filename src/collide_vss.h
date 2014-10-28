/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COLLIDE_CLASS

CollideStyle(vss,CollideVSS)

#else

#ifndef SPARTA_COLLIDE_VSS_H
#define SPARTA_COLLIDE_VSS_H

#include "collide.h"
#include "particle.h"

namespace SPARTA_NS {

class CollideVSS : public Collide {
 public:
  CollideVSS(class SPARTA *, int, char **);
  ~CollideVSS();
  void init();

  double vremax_init(int, int);
  double attempt_collision(int, int, double);
  double attempt_collision(int, int, int, double);
  int test_collision(int, int, int, Particle::OnePart *, Particle::OnePart *);
  void setup_collision(Particle::OnePart *, Particle::OnePart *);
  Particle::OnePart *
    perform_collision(Particle::OnePart *&, Particle::OnePart *&);

  double extract(int, const char *);

 private:
  int eng_exchange;
  double vr_indice;
  double **prefactor; // static portion of collision attempt frequency

  struct State {      // two-particle state
    double vr2;
    double vr;
    double mr;
    double rotdof_i;
    double rotdof_j;
    double vibdof_i;
    double vibdof_j;
    double ave_rotdof;
    double ave_vibdof;
    double ave_dof;
    double etrans;
    double erot;
    double evib;
    double eexchange;
    double eint;
    double etotal;
    double mass_i;
    double mass_j;
  };
 
  struct State precoln;       // state before collision
  struct State postcoln;      // state after collision

  struct Params {             // VSS model parameters
    double diam;
    double omega;
    double tref;
    double alpha;
  };
  
  Params *params;             // VSS params for each species
  int nparams;                // # of per-species params read in

  void SCATTER_TwoBodyScattering(Particle::OnePart *, 
				 Particle::OnePart *);
  void EEXCHANGE_NonReactingEDisposal(Particle::OnePart *, 
				      Particle::OnePart *);
  void SCATTER_ThreeBodyScattering(Particle::OnePart *, 
                                   Particle::OnePart *,
                                   Particle::OnePart *);
  void EEXCHANGE_ReactingEDisposal(Particle::OnePart *, 
                                   Particle::OnePart *,
                                   Particle::OnePart *);

  double sample_bl(RanPark *, double, double);
  void read_param_file(char *);
  int wordcount(char *, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Species %s did not appear in VSS parameter file

Self-explanatory.

E: VSS parameters do not match current species

Species cannot be added after VSS colision file is read.

E: Cannot open VSS parameter file %s

Self-explantory.

E: Incorrect line format in VSS parameter file

Number of parameters in a line read from file is not valid.

E: Request for unknown parameter from collide

VSS model does not have the parameter being requested.

*/
