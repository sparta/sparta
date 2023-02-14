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
  int react(Particle::OnePart *&, int, double *, Particle::OnePart *&, int &);

  char *reactionID(int);
  int match_reactant(char *, int);
  int match_product(char *, int);

  void tally_update();

 private:
  int me,nprocs;
  int gsflag,psflag;                // 0/1 if gas and/or surf chem enabled
  int mode;                         // FACE or SURF
  int nsync;                        // synchronize surf state
                                    // every this many steps
  double twall;                     // temperature of face or surf
  double max_cover;
  int this_index;                   // index of this surf reaction model
                                    // in Surf list of all reaction models

  class RanKnuth *random;     // RNG for reaction probabilities

  int nspecies_surf;         // number of surface species
  char **species_surf;       // list of surface species

  // mode = FACE for box faces
  // all this data is allocated here

  int nface;       // # of box faces, 4 (2d) or 6 (3d)

  int **face_species_state;     // 4 state quantities for up to 6 box faces
  int *face_total_state;
  double *face_area;
  double *face_weight;
  double **face_tau;

  int **face_species_delta;     // changes to state between syncs
  int **face_sum_delta;         // delta summed across all procs
  double **face_norm;           // norm of each face

  // mode = SURF for surface elements (lines or tris)

  int nstick_species_custom;    // indices to custom state in Surf
  int nstick_total_custom;
  int area_custom,weight_custom;
  int tau_custom;
  int first_owner;       // 1 if this instance of SRA allocates custom Surf data

  int **surf_species_state;     // ptrs to custom state vecs/arrays in Surf class
  int *surf_total_state;
  double *surf_area;
  double *surf_weight;
  double **surf_tau;

  int **surf_species_delta;     // changes to state between syncs

  int *mark;               // per-surf mark = 1 if reaction has occured, else 0
  surfint *tally2surf;     // global surf index for each entry in incollate
  int **intally,**outtally;      // used for Allreduce of state changes
  double **incollate,**outcollate;   // used to collate state changes across procs
  int maxtally;                  // allocated size of intally

  // ptrs to data for each box face or surface element
  // used in react() and react_PS() and sync operations

  int **species_delta;       // change in perspecies count since last sync
  int **species_state;       // perspecies count at last sync
  int *total_state;          // total count at last sync
  double *area;              // area of surf
  double *weight;            // weight of surf
  double **tau;              // PS time of surf

  // GS (gas/surface) reaction model

  struct OneReaction_GS {
    char *id;                      // reaction ID (formula)
    int active;                    // 1 if reaction is active
    int type;                      // reaction type = DISSOCIATION, etc
    int style;                     // reaction style = ARRHENIUS, etc
    int ncoeff;                    // # of numerical coeffs
    int nreactant,nproduct;        // # of reactants and products
    int nprod_g,nprod_g_tot;       // # of products which are gaseous species
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
    int kisliuk_flag, energy_flag;
    double kisliuk_coeff[3], energy_coeff[2];
    int cmodel_ip;                  // style for I's post-reaction surf collision
    int *cmodel_ip_flags;           // integer flags to pass to SC class
    double *cmodel_ip_coeffs;       // double coeffs to pass to SC class
    int cmodel_jp;                  // difto for J particle
    int *cmodel_jp_flags;
    double *cmodel_jp_coeffs;
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

 // PS (on-surf) reaction model

 struct OneReaction_PS {
    char *id;                          // reaction ID (formula)
    int index;                         // index of the reaction
    int active;                        // 1 if reaction is active
    int type;                          // reaction type = DISSOCIATION, etc
    int style;                         // reaction style = ARRHENIUS, etc
    int ncoeff;                        // # of numerical coeffs
    int nreactant,nproduct;            // # of reactants and products
    int nprod_g, nprod_g_tot;          // # of products which are gaseous species
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
    int cmodel_ip;                  // style for I's post-reaction surf collision
    int *cmodel_ip_flags;           // integer flags to pass to SC class
    double *cmodel_ip_coeffs;       // double coeffs to pass to SC class
    int cmodel_jp;                  // difto for J particle
    int *cmodel_jp_flags;
    double *cmodel_jp_coeffs;
  };

  OneReaction_PS *rlist_ps;           // list of all reactions read from file
  int nlist_ps;                       // # of reactions read from file
  int maxlist_ps;                     // max # of reactions in rlist

  int nactive_ps;
  int *reactions_ps_list;
  // SGK check
  int n_PS_react;

  // particles added to gas flow by PS surface reactions

  struct AddParticle {
    int id;                 // particle ID
    int ispecies;           // particle species index
    double x[3];            // particle position
    double v[3];            // particle velocity
    double erot;            // rotational energy
    double evib;            // vibrational energy
    double dtremain;        // fraction of timestep
  };

  AddParticle *mypart;      // particles this proc adds
  AddParticle *allpart;     // gathered list of all particles from all procs
  int *recvcounts,*displs;  // Nproc-length vectors for Allgatherv
  int npart;                // # of particles this proc adds
  int maxmypart;            // allocated size of mypart
  int maxallpart;           // allocated size of allpart

  // surface collision models, one per supported SC style
  // only non-NULL if the SC style appears in GS/PS reaction files

  class SurfCollide **cmodels;

  // pointers to Compute instances which tally reactions on per-surf basis
  // extracted from Update class

  int ncompute_tally;
  class Compute **clist_active;

  // GS methods

  void init_reactions_gs();
  void readfile_gs(char *);

  // PS methods

  void init_reactions_ps();
  void readfile_ps(char *);
  void PS_react(int, int, double *);
  void add_particle_mine(Particle::OnePart *);
  void PS_chemistry();
  void random_point(int, double*);

  // methods common to both GS and PS

  void create_per_face_state();
  void create_per_surf_state();

  void update_state_face();
  void update_state_surf();

  // NOTE: can remove these 3 at some point
  /*
  void energy_barrier_scatter(Particle::OnePart*, double *,
                              double, double, double);
  void non_thermal_scatter(Particle::OnePart*, double *,
                           double, double, double, double);
  void cll(Particle::OnePart *, double *, double, double, double);
  */

  double stoich_pow(int, int);
  int find_surf_species(char *);
  void print_reaction(char *, char *);
  int readone(char *, char *, int &, int &);
  int readextra(int, char *, char *, int &, int &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Surf_adsorb ID must be alphanumeric or underscore characters

Self-explanatory.

*/
