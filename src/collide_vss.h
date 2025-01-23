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
  virtual ~CollideVSS();
  virtual void init();

  double vremax_init(int, int);
  // Virgile - Modif Start - 20/10/23
  // ========================================================================
  // Add count_wi in the variable definition.
  // ========================================================================
  // Baseline code:
  // virtual double attempt_collision(int, int, double);
  // double attempt_collision(int, int, int, double);
  // Modified code:
  virtual double attempt_collision(int, int, double, double, double);
  double attempt_collision(int, int, int, double);
  // Virgile - Modif End - 20/10/23
  
  // Takato Morimoto - Modif Start - 20/10/23
  // Declarations for two speices Millikan-White VT mode are made
  struct Mwcoeff {             // Coefficient of MW-park model parameters
    double a;
    double b;
  };
  void read_param_file_tv_mw(char *);
  protected:
    int TV_MWflag; 
  Mwcoeff **mwcoeff;
  // Takato Morimoto - Modif End - 22/05/24
  
  // Virgile - Modif Start - 19/12/24
  // ========================================================================
  // Add maxwi to the test_collision function for SWSmax filter.
  // ========================================================================
  // Baseline code:
  // virtual int test_collision(int, int, int, Particle::OnePart *, Particle::OnePart *);
  // Modified code:
  virtual int test_collision(int, int, int, Particle::OnePart *, Particle::OnePart *, double);
  // Virgile - Modif End - 19/12/24
  
  virtual void setup_collision(Particle::OnePart *, Particle::OnePart *);
  //Takato Morimoto - Modif Start - 05/06/24
  //virtual int perform_collision(Particle::OnePart *&, Particle::OnePart *&,
  //                      Particle::OnePart *&);
  virtual int perform_collision(Particle::OnePart *&, Particle::OnePart *&,
                        Particle::OnePart *&, int &,int &,int &,int &);
  //Takato Morimoto - Modif End - 05/06/24
  double extract(int, int, const char *);

  struct State {      // two-particle state
    double vr2;
    double vr;
    double imass,jmass;
    double ave_rotdof;
    double ave_vibdof;
    double ave_dof;
    double etrans;
    double erot;
    double evib;
    double eexchange;
    double eint;
    double etotal;
    double ucmf;
    double vcmf;
    double wcmf;
  };

  struct Params {             // VSS model parameters
    double diam;
    double omega;
    double tref;
    double alpha;
    double rotc1;
    double rotc2;
    double rotc3;
    double vibc1;
    double vibc2;
    double mr;
  };

 protected:
  int relaxflag,eng_exchange;
  double vr_indice;
  double **prefactor; // static portion of collision attempt frequency

  struct State precoln;       // state before collision
  struct State postcoln;      // state after collision

  Params **params;             // VSS params for each species
  int nparams;                // # of per-species params read in

  // Virgile - Modif Start - 24/10/2024
  // Baseline code:
  // void SCATTER_TwoBodyScattering(Particle::OnePart *,
	// 			 Particle::OnePart *);
  // Modified code:
  void SCATTER_TwoBodyScattering(Particle::OnePart *,
				 Particle::OnePart *, int);
  // Virgile - Modif End - 24/10/2024
  void EEXCHANGE_NonReactingEDisposal(Particle::OnePart *,
                                      Particle::OnePart *);
  void SCATTER_ThreeBodyScattering(Particle::OnePart *,
                                   Particle::OnePart *,
                                   Particle::OnePart *);
  void EEXCHANGE_ReactingEDisposal(Particle::OnePart *,
                                   Particle::OnePart *,
                                   Particle::OnePart *);

  double sample_bl(RanKnuth *, double, double);
  double rotrel (int, double);
  // Takato Morimoto - Modif Start - 24/05/23
  //double vibrel (int, double); //baseline code
  double vibrel (int, int, double); //input is increased
  // Takato Morimoto - Modif End - 24/05/23

  void read_param_file(char *);
  int wordparse(int, char *, char **);
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

Self-explanatory.

E: Incorrect line format in VSS parameter file

Number of parameters in a line read from file is not valid.

E: Request for unknown parameter from collide

VSS model does not have the parameter being requested.

*/
