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

#ifdef REACT_CLASS

ReactStyle(tce,ReactTCE)

#else

#ifndef SPARTA_REACT_TCE_H
#define SPARTA_REACT_TCE_H

#include "react.h"
#include "particle.h"

namespace SPARTA_NS {

class ReactTCE : public React {
 public:
  ReactTCE(class SPARTA *, int, char **);
  ~ReactTCE();
  void init();
  int attempt(Particle::OnePart *, Particle::OnePart *, 
              double, double, double, double &, int &);

 private:
  FILE *fp;

  struct OneReaction {
    int active;                    // 1 if reaction is active
    int initflag;                  // 1 if reaction params have been init
    int type;                      // reaction type = DISSOCIATION, etc
    int style;                     // reaction style = ARRHENIUS, etc
    int ncoeff;                    // # of numerical coeffs
    int nreactant,nproduct;        // # of reactants and products
    char **id_reactants,**id_products;  // species IDs of reactants/products
    int *reactants,*products;      // species indices of reactants/products
    double *coeff;                 // numerical coeffs for reaction
  };

  OneReaction *rlist;              // list of all reactions read from file
  int nlist;                       // # of reactions read from file
  int maxlist;                     // max # of reactions in rlist

  // possible reactions a pair of reactant species is part of

  struct ReactionIJ {
    int *list;           // list of indices into rlist, ptr into indices
    int n;               // # of reactions in list
  };

  ReactionIJ **reactions;     // reactions for all IJ pairs of reactant species
  int *indices;               // master list of indices

  void readfile(char *);
  int readone(char *, char *, int &, int &);
};

}

#endif
#endif
