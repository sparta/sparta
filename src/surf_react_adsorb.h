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

#ifdef SURF_REACT_CLASS

SurfReactStyle(adsorb,SurfReactAdsorb)

#else

#ifndef SPARTA_SURF_REACT_ADSORB_H
#define SPARTA_SURF_REACT_ADSORB_H

#include "surf_react.h"

namespace SPARTA_NS {

class SurfReactAdsorb : public SurfReact {
 public:
  SurfReactAdsorb(class SPARTA *, int, char **);
  ~SurfReactAdsorb();
  void init();
  int react(Particle::OnePart *&, int, double *, Particle::OnePart *&);

  char *reactionID(int);
  int match_reactant(char *, int);
  int match_product(char *, int);

 private:
  int model,mode,nsp,nface;
  double twall,max_cover;
  int nspecies_surf;
  char **species_surf;       // list of surface species

  FILE *fp;                  // NOTE: dup of what is in parent class?
  class RanPark *random;     // RNG for reaction probabilities

  // surface state
  // either for box faces or per surface element, depending on mode

  int custom_owner;
  int nstick_index,nstick_direct;
  int nstick_total_index,nstick_total_direct;
  int area_index,area_direct;
  int weight_index,weight_direct;

  int **nstick_face;
  int nstick_total_face[6];
  double area_face[6];
  double weight_face[6];

  // GS react model
  
  struct OneReaction_GS {
    char *id;                      // reaction ID (formula)
    int active;                    // 1 if reaction is active
    int type;                      // reaction type = DISSOCIATION, etc
    int style;                     // reaction style = ARRHENIUS, etc
    int ncoeff;                    // # of numerical coeffs
    int nreactant,nproduct;        // # of reactants and products
    char **id_reactants,**id_products;  // species IDs of reactants/products
    char **state_reactants,**state_products; // state of reactants and products
    int *part_reactants,*part_products; // participation of reactants & products
                                 // 0 for catalyst 1 for participating reactions
    int *stoich_reactants, *stoich_products; // stoichiometric coefficients of
                                             // reactants and products 
    int *reactants,*products;      // species indices of reactants/products
    int *reactants_ad_index,*products_ad_index;
    double *coeff;                 // numerical coeffs for reaction
    double k_react;
    int cmodel;                    // style for post-reaction surf collisions
    int *cmodel_flags;             // integer flags to pass to SC class
    double *cmodel_coeffs;         // double coeffs to pass to SC class
  };

  OneReaction_GS *rlist_gs;           // list of all reactions read from file
  int nlist_gs;                       // # of reactions read from file
  int maxlist_gs;                     // max # of reactions in rlist

  struct ReactionI_GS {  // possible reactions a reactant species is part of
    int *list;           // list of indices into rlist, ptr into indices
    int n;               // # of reactions in list
  };

  ReactionI_GS *reactions_gs;    // reactions for all species
  int *indices_gs;               // master list of indices

 // PS react model
  
 struct OneReaction_PS {
    int index;                          // index of the reaction
    int active;                        // 1 if reaction is active
    int type;                          // reaction type = DISSOCIATION, etc
    int style;                         // reaction style = ARRHENIUS, etc
    int ncoeff;                        // # of numerical coeffs
    int nreactant,nproduct;            // # of reactants and products 
    char **id_reactants,**id_products; // species IDs of reactants & products
    char **state_reactants,**state_products;  // state of reactants & products
    int *part_reactants,*part_products; // participation of reactants & products
                                 // 0 for catalyst 1 for participating reactions
    int *stoich_reactants, *stoich_products; // stoichiometric coefficients of 
                                             // reactants and products 
    int *reactants,*products;          // species indices of reactants/products
    int *reactants_ad_index, *products_ad_index;   // adsorbed species index of
                                       // reactants and products
    double *coeff;                     // numerical coeffs for reaction
    double k_react;
    char *id;                          // reaction ID (formula)
  };

  OneReaction_PS *rlist_ps;           // list of all reactions read from file
  int nlist_ps;                       // # of reactions read from file
  int maxlist_ps;                     // max # of reactions in rlist
  int nactive_ps;
  int n_PS_react;

  // surface collision models, one per supported style
  // only if appears in reaction file

  class SurfCollide **cmodels;

  // GS methods

  void init_reactions_gs();
  void readfile_gs(char *);

  // PS methods
 
  //void init_reactions_ps();
  //void readfile_ps(char *);
  //void random_point(int, int, int, double*);
  //int find_cell(int, int, double*); 
  
  // methods common to both GS and PS

  void create_per_face_state();
  void create_per_surf_state();
  void energy_barrier_scatter(Particle::OnePart*, double *, 
                              double, double, double);
  void non_thermal_scatter(Particle::OnePart*, double *, 
                           double, double, double, double);
  void cll(Particle::OnePart *, double *, double, double, double);

  double stoich_pow(int, int);
  int find_surf_species(char *);
  void print_reaction(char *, char *, char *);
  int readone(char *, char *, char *, int &, int &, int &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Surf_adsorb ID must be alphanumeric or underscore characters

Self-explanatory.

*/
