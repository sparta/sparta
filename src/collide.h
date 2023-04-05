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

#ifndef SPARTA_COLLIDE_H
#define SPARTA_COLLIDE_H

#include "pointers.h"
#include "memory.h"
#include "particle.h"

namespace SPARTA_NS {

#define DELTAPART 128

class Collide : protected Pointers {
 public:
  char *style;
  int rotstyle;       // none/smooth rotational modes
  int vibstyle;       // none/discrete/smooth vibrational modes
  int nearcp;         // 1 for near neighbor collisions
  int nearlimit;      // limit on neighbor serach for near neigh collisions

  int ncollide_one,nattempt_one,nreact_one;
  bigint ncollide_running,nattempt_running,nreact_running;

  Collide(class SPARTA *, int, char **);
  virtual ~Collide();
  virtual void init();
  void modify_params(int, char **);
  void reset_vremax();
  virtual void collisions();

  virtual double vremax_init(int, int) = 0;
  virtual double attempt_collision(int, int, double) = 0;
  virtual double attempt_collision(int, int, int, double) = 0;
  virtual int test_collision(int, int, int,
                             Particle::OnePart *, Particle::OnePart *) = 0;
  virtual void setup_collision(Particle::OnePart *, Particle::OnePart *) = 0;
  virtual int perform_collision(Particle::OnePart *&, Particle::OnePart *&,
                                Particle::OnePart *&) = 0;

  virtual double extract(int, int, const char *) {return 0.0;}

  virtual int pack_grid_one(int, char *, int);
  virtual int unpack_grid_one(int, char *);
  virtual void copy_grid_one(int, int);
  virtual void reset_grid_count(int);
  virtual void add_grid_one();
  virtual void adapt_grid();

  int ngroups;        // # of groups

 protected:
  int npmax;          // max # of particles in plist
  int *plist;         // list of particle indices for the entire cell

  int nglocal;        // current size of per-cell arrays
  int nglocalmax;     // max allocated size of per-cell arrays (vremax, remain)

  int *ngroup;        // # of particles in each group
  int *maxgroup;      // max # of particles allocated per group
  int **glist;        // indices into plist of particles in each group
  int **p2g;          // for each plist entry: 0 = igroup, 1 = index within glist

  int npair;          // # of group pairs to do collisions for
  int **gpair;        // Npairx3 list of group pairs to do collisions for
                      // 0 = igroup, 1 = jgroup, 2 = # of attempt collisions

  int max_nn;             // allocated size of nn_last_partner
  int *nn_last_partner;   // plist index+1 of last collision partner for each particle
                          // 0 = no collision yet (on this step)
  int *nn_last_partner_igroup;   // ditto for two groups of particles
  int *nn_last_partner_jgroup;

  int ndelete,maxdelete;      // # of particles removed by chemistry
  int *dellist;               // list of particle indices to delete

  char *mixID;               // ID of mixture to use for groups
  class Mixture *mixture;    // ptr to mixture
  class RanKnuth *random;     // RNG for collision generation

  int vre_first;      // 1 for first run after collision style is defined
  int vre_start;      // 1 if reset vre params at start of each run
  int vre_every;      // reset vre params every this many steps
  bigint vre_next;    // next timestep to reset vre params on
  int remainflag;     // 1 if remain defined, else use random fraction

  double ***vremax;   // max relative velocity, per cell, per group pair
  double ***remain;   // collision number remainder, per cell, per group pair
  double **vremax_initial;   // initial vremax value, per group pair

  // recombination reactions

  int recombflag;               // 1 if recomb reactions enabled, 0 if not
  double recomb_boost_inverse;  // recombination rate boost factor from React
  int **recomb_ijflag;          // 1 if species I,J have recomb reaction(s)

  // discrete vibrational energy data structs

  int index_vibmode;   // index to custom vibmode vector

  // ambipolar approximation data structs

  int ambiflag;       // 1 if ambipolar option is enabled
  int ambispecies;    // species for ambipolar electrons
  int index_ionambi;  // 2 custom ambipolar vectors
  int index_velambi;

  int maxelectron;              // max # elist can hold
  Particle::OnePart *elist;     // list of ambipolar electrons
                                // for one grid cell or pair of groups in cell
  // Kokkos data

  int oldgroups;         // pass from parent to child class
  int copymode;          // 1 if copy of class (prevents deallocation of
                         //   base class when child copy is destroyed)
  int kokkos_flag;        // 1 if collide method supports Kokkos

  // inline functions
  // add particle N to Igroup and set its g2p entry in plist to K
  // delete Ith entry in Igroup and reset g2p entries as well

  inline void addgroup(int igroup, int pindex)
  {
    if (ngroup[igroup] == maxgroup[igroup]) {
      maxgroup[igroup] += DELTAPART;
      memory->grow(glist[igroup],maxgroup[igroup],"collide:grouplist");
    }
    int ng = ngroup[igroup];
    glist[igroup][ng] = pindex;
    p2g[pindex][0] = igroup;
    p2g[pindex][1] = ng;
    ngroup[igroup]++;
  }

  inline void delgroup(int igroup, int i)
  {
    int ng = ngroup[igroup];
    if (i < ng-1) {
      glist[igroup][i] = glist[igroup][ng-1];
      int pindex = glist[igroup][i];
      p2g[pindex][0] = igroup;
      p2g[pindex][1] = i;
    }
    ngroup[igroup]--;
  }

  template < int > void collisions_one();
  template < int > void collisions_group();
  void collisions_one_ambipolar();
  void collisions_group_ambipolar();
  void ambi_reset(int, int, int, Particle::OnePart *, Particle::OnePart *,
                  Particle::OnePart *, int *);
  void ambi_check();
  void grow_percell(int);

  int find_nn(int, int);
  int find_nn_group(int, int *, int, int *, int *, int *, int *);
  void realloc_nn(int, int *&);
  void set_nn(int);
  void set_nn_group(int);
};

}

#endif

/* ERROR/WARNING messages:

E: Collision mixture does not exist

Self-explanatory.

E: Collision mixture does not contain all species

The specified mixture must contain all species in the simulation so
that they can be assigned to collision groups.

E: Must use Kokkos-supported collision style if Kokkos is enabled

Self-explanatory.

E: Cannot (yet) use KOKKOS package with 'collide_modify vibrate discrete'

This feature is not yet supported.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
