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

CollideStyle(vss/kk,CollideVSSKokkos)

#else

#ifndef SPARTA_COLLIDE_VSS_KOKKOS_H
#define SPARTA_COLLIDE_VSS_KOKKOS_H

#include "collide_vss.h"
#include "collide_vss_kokkos.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "react_tce_kokkos.h"
#include "react_qk_kokkos.h"
#include "react_tce_qk_kokkos.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "rand_pool_wrap.h"
#include "kokkos_copy.h"
#include "compute_gas_collision_grid_kokkos.h"
#include "compute_gas_reaction_grid_kokkos.h"
#include "compute_gas_collision_tally_kokkos.h"
#include "compute_gas_reaction_tally_kokkos.h"


namespace SPARTA_NS {

struct s_COLLIDE_REDUCE {
  // bigint since can exceed 2^31 in one step
  //   at large per-proc particle counts
  bigint nattempt_one,ncollide_one,nreact_one;
  KOKKOS_INLINE_FUNCTION
  s_COLLIDE_REDUCE() {
    nattempt_one = 0;
    ncollide_one = 0;
    nreact_one = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const s_COLLIDE_REDUCE &rhs) {
    nattempt_one += rhs.nattempt_one;
    ncollide_one += rhs.ncollide_one;
    nreact_one += rhs.nreact_one;
  }
};
typedef struct s_COLLIDE_REDUCE COLLIDE_REDUCE;

struct TagCollideResetVremax{};
struct TagCollideZeroNN{};

template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
struct TagCollideCollisionsOne{};

template < int DIM, int GASTALLY, int ATOMIC_REDUCTION >
struct TagCollideCollisionsOneSubcell{};

template < int GASTALLY, int ATOMIC_REDUCTION >
struct TagCollideCollisionsOneAmbipolar{};

template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
struct TagCollideCollisionsGroup{};

template < int GASTALLY, int ATOMIC_REDUCTION >
struct TagCollideCollisionsGroupAmbipolar{};

class CollideVSSKokkos : public CollideVSS {
 public:
  typedef COLLIDE_REDUCE value_type;

  CollideVSSKokkos(class SPARTA *, int, char **);
  ~CollideVSSKokkos();
  void init();
  void collisions();
  void sync(ExecutionSpace, unsigned int);
  void modified(ExecutionSpace, unsigned int);

#ifndef SPARTA_KOKKOS_EXACT
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

  //Kokkos::Random_XorShift1024_Pool<DeviceType> rand_pool;
  //typedef typename Kokkos::Random_XorShift1024_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif

  KOKKOS_INLINE_FUNCTION
  double attempt_collision_kokkos(int, int, double, rand_type &) const;
  KOKKOS_INLINE_FUNCTION
  double attempt_collision_kokkos(int, int, int, int, int, double, rand_type &) const;
  KOKKOS_INLINE_FUNCTION
  double poisson_kokkos(double, rand_type &) const;
  KOKKOS_INLINE_FUNCTION
  int test_collision_kokkos(int, int, int, Particle::OnePart *, Particle::OnePart *, struct State &, rand_type &) const;
  KOKKOS_INLINE_FUNCTION
  void setup_collision_kokkos(Particle::OnePart *, Particle::OnePart *, struct State &, struct State &) const;
  KOKKOS_INLINE_FUNCTION
  int perform_collision_kokkos(int, Particle::OnePart *&, Particle::OnePart *&,
                        Particle::OnePart *&, struct State &, struct State &, rand_type &,
                        Particle::OnePart *&, int &, double &,
                        int &) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideResetVremax, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideZeroNN, const int&) const;

  template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsOne< NEARCP, GASTALLY, ATOMIC_REDUCTION >, const int&) const;

  template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsOne< NEARCP, GASTALLY, ATOMIC_REDUCTION >, const int&, COLLIDE_REDUCE&) const;

  template < int DIM, int GASTALLY, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsOneSubcell< DIM, GASTALLY, ATOMIC_REDUCTION >, const int&) const;

  template < int DIM, int GASTALLY, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsOneSubcell< DIM, GASTALLY, ATOMIC_REDUCTION >, const int&, COLLIDE_REDUCE&) const;

  template < int GASTALLY, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsOneAmbipolar< GASTALLY, ATOMIC_REDUCTION >, const int&) const;

  template < int GASTALLY, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsOneAmbipolar< GASTALLY, ATOMIC_REDUCTION >, const int&, COLLIDE_REDUCE&) const;

  template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsGroup< NEARCP, GASTALLY, ATOMIC_REDUCTION >, const int&) const;

  template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsGroup< NEARCP, GASTALLY, ATOMIC_REDUCTION >, const int&, COLLIDE_REDUCE&) const;

  template < int GASTALLY, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsGroupAmbipolar< GASTALLY, ATOMIC_REDUCTION >, const int&) const;

  template < int GASTALLY, int ATOMIC_REDUCTION >
  KOKKOS_INLINE_FUNCTION
  void operator()(TagCollideCollisionsGroupAmbipolar< GASTALLY, ATOMIC_REDUCTION >, const int&, COLLIDE_REDUCE&) const;

  typedef Kokkos::
    DualView<Params**, Kokkos::LayoutRight, DeviceType> tdual_params_2d;
  typedef tdual_params_2d::t_dev t_params_2d;
  typedef tdual_params_2d::t_dev_const t_params_2d_const;
  t_params_2d_const d_params_const;

 private:
  KOKKOS_INLINE_FUNCTION
  void ambi_reset_kokkos(int, int, int, int,
                         Particle::OnePart *, Particle::OnePart *,
                         Particle::OnePart *, const DAT::t_int_1d &) const;
  void reset_vremax();
  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void copy_grid_one(int, int);
  void reset_grid_count(int);
  void add_grid_one();
  void adapt_grid();
  void grow_percell(int);

  KKCopy<GridKokkos> grid_kk_copy;
  KKCopy<ReactTCEKokkos> react_kk_copy;
  KKCopy<ReactQKKokkos> react_qk_kk_copy;
  KKCopy<ReactTCEQKKokkos> react_tceqk_kk_copy;
  int react_style;   // 0=TCE, 1=QK, 2=TCEQK (set in setup)

  // active gas/gas tally computes, partitioned by type.  Two representations,
  //   selected by SPARTA_KOKKOS_FIXED_LISTS (see kokkos_type.h); the kernel
  //   dispatch sites are written once against the CVK_* accessors below.

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  KKCopy<ComputeGasCollisionGridKokkos> glist_collision_copy[KOKKOS_MAX_GLIST];
  KKCopy<ComputeGasCollisionTallyKokkos> glist_coll_tally_copy[KOKKOS_MAX_GLIST];
  KKCopy<ComputeGasReactionTallyKokkos> glist_react_tally_copy[KOKKOS_MAX_GLIST];
  KKCopy<ComputeGasReactionGridKokkos> glist_reaction_copy[KOKKOS_MAX_GLIST];
  ComputeGasCollisionGridKokkos tmp_compute_gas_collision_kk;
  ComputeGasReactionGridKokkos tmp_compute_gas_reaction_kk;
  ComputeGasCollisionTallyKokkos tmp_compute_gas_coll_tally_kk;
  ComputeGasReactionTallyKokkos tmp_compute_gas_react_tally_kk;

#define CVK_GLIST_COLLISION(m)   glist_collision_copy[m].obj
#define CVK_GLIST_REACTION(m)    glist_reaction_copy[m].obj
#define CVK_GLIST_COLL_TALLY(m)  glist_coll_tally_copy[m].obj
#define CVK_GLIST_REACT_TALLY(m) glist_react_tally_copy[m].obj

#else
  DAT::tdual_char_1d k_glist_collision, k_glist_reaction,
                     k_glist_coll_tally, k_glist_react_tally;
  DAT::t_char_1d d_glist_collision, d_glist_reaction,
                 d_glist_coll_tally, d_glist_react_tally;

#define CVK_GLIST_COLLISION(m)   ((const ComputeGasCollisionGridKokkos *) d_glist_collision.data())[m]
#define CVK_GLIST_REACTION(m)    ((const ComputeGasReactionGridKokkos *) d_glist_reaction.data())[m]
#define CVK_GLIST_COLL_TALLY(m)  ((const ComputeGasCollisionTallyKokkos *) d_glist_coll_tally.data())[m]
#define CVK_GLIST_REACT_TALLY(m) ((const ComputeGasReactionTallyKokkos *) d_glist_react_tally.data())[m]
#endif

  int nglist_collision,nglist_reaction;
  int nglist_coll_tally,nglist_react_tally;
  DAT::t_int_scalar d_tally_overflow;
  HAT::t_int_scalar h_tally_overflow;
  void grow_gas_tally_computes();
  void rewind_gas_tally_computes(int);
  void setup_gas_tally();
  void finish_gas_tally();
  void clear_gas_tally();

  t_particle_1d d_particles;
  t_species_1d_const d_species;
  DAT::t_int_2d d_plist;

  DAT::t_int_1d d_nelecstates;
  t_elecstate_2d d_elecstates;
  DAT::t_float_2d d_elec_default_rels;
  DAT::t_float_3d d_elec_species_rels;
  DAT::t_int_2d d_enforce_spin_conservation;
  DAT::t_float_2d d_cumulative_probabilities;

  // group collision scratch (ngroups > 1)
  DAT::t_int_1d d_species2group;
  // reacting group collisions mutate group membership inside the kernel, so
  //   the per-group lists cannot be one group-contiguous row per cell as they
  //   were when only the non-reacting case was supported.  Each group gets its
  //   own region of capacity d_plist.extent(1) -- a group can never hold more
  //   than the cell does -- and d_p2g is the reverse map the host keeps in
  //   Collide::p2g, needed so a swap-remove can fix the moved entry's owner.

  Kokkos::View<int***,DeviceType> d_glist;   // (cell, group, k) -> plist index
  Kokkos::View<int***,DeviceType> d_p2g;     // (cell, plist index) -> group, k

  // per-group counters, formerly per-thread stack arrays dimensioned by a
  //   compile-time MAXGROUP.  Both group kernels are one work item per grid
  //   cell (RangePolicy over 0..nglocal indexed by icell), so icell is already
  //   a unique, deterministic, contention-free row index -- no UniqueToken
  //   needed.  Sized (nglocal, ngroups) alongside d_glist, which they cost
  //   1/d_plist.extent(1) as much as.
  //   both group kernels build their lists with addgroup_kk(), which keeps
  //   d_gcount and d_p2g in step, so no separate fill cursor is needed.

  Kokkos::View<int**,DeviceType> d_gcount;   // (cell, group) -> # in group

  // near-neighbor partner history for the two groups of the current pair;
  //   the host reallocates these per pair via set_nn_group()

  DAT::t_int_2d d_nn_igroup;
  DAT::t_int_2d d_nn_jgroup;

 public:

  // mirror Collide::addgroup / delgroup (collide.h:157-179) exactly, including
  //   the swap-with-last order: a reaction that rebins a particle changes which
  //   index a later random draw lands on, so any deviation diverges from the
  //   host rather than merely reordering

  // the group counts live in d_gcount(icell,*), so icell is all these need

  KOKKOS_INLINE_FUNCTION
  void addgroup_kk(const int icell, const int igroup, const int pindex) const
  {
    const int ng = d_gcount(icell,igroup);
    d_glist(icell,igroup,ng) = pindex;
    d_p2g(icell,pindex,0) = igroup;
    d_p2g(icell,pindex,1) = ng;
    d_gcount(icell,igroup)++;
  }

  KOKKOS_INLINE_FUNCTION
  void delgroup_kk(const int icell, const int igroup, const int i) const
  {
    const int ng = d_gcount(icell,igroup);
    if (i < ng-1) {
      d_glist(icell,igroup,i) = d_glist(icell,igroup,ng-1);
      const int pindex = d_glist(icell,igroup,i);
      d_p2g(icell,pindex,0) = igroup;
      d_p2g(icell,pindex,1) = i;
    }
    d_gcount(icell,igroup)--;
  }

 private:

  Kokkos::View<int***,DeviceType> d_nattempt_pair;

  DAT::t_int_1d d_ewhich;
  tdual_struct_tdual_int_1d_1d k_eivec;
  tdual_struct_tdual_int_2d_1d k_eiarray;
  tdual_struct_tdual_float_1d_1d k_edvec;
  tdual_struct_tdual_float_2d_1d k_edarray;
  DAT::t_int_1d d_ionambi;
  DAT::t_int_1d d_ions;
  DAT::t_float_2d_lr d_velambi;
  t_particle_2d d_elist;

  DAT::tdual_float_2d k_vremax_initial;
  DAT::t_float_2d d_vremax_initial;
  DAT::tdual_float_3d k_vremax;
  DAT::t_float_3d d_vremax;
  DAT::tdual_float_3d k_remain;
  DAT::t_float_3d d_remain;

  // int scalars = flags and view-size counters, must stay int
  // bigint scalars = per-step statistics counters, can exceed 2^31
  //   in one step at large per-proc particle counts

  typedef Kokkos::DualView<int[9], DeviceType::array_layout, DeviceType> tdual_int_8;
  typedef tdual_int_8::t_dev t_int_8;
  typedef tdual_int_8::t_host t_host_int_8;
  t_int_8 d_scalars;
  t_host_int_8 h_scalars;

  typedef Kokkos::DualView<bigint[3], DeviceType::array_layout, DeviceType> tdual_bigint_3;
  typedef tdual_bigint_3::t_dev t_bigint_3;
  typedef tdual_bigint_3::t_host t_host_bigint_3;
  t_bigint_3 d_scalars_big;
  t_host_bigint_3 h_scalars_big;

  DAT::t_bigint_scalar d_nattempt_one;
  HAT::t_bigint_scalar h_nattempt_one;

  DAT::t_bigint_scalar d_ncollide_one;
  HAT::t_bigint_scalar h_ncollide_one;

  DAT::t_bigint_scalar d_nreact_one;
  HAT::t_bigint_scalar h_nreact_one;

  DAT::t_int_scalar d_error_flag;
  HAT::t_int_scalar h_error_flag;

  DAT::t_int_scalar d_retry;
  HAT::t_int_scalar h_retry;

  DAT::t_int_scalar d_maxdelete;
  HAT::t_int_scalar h_maxdelete;

  DAT::t_int_scalar d_maxcellcount;
  HAT::t_int_scalar h_maxcellcount;

  DAT::t_int_scalar d_part_grow;
  HAT::t_int_scalar h_part_grow;

  DAT::t_int_scalar d_ndelete;
  HAT::t_int_scalar h_ndelete;

  DAT::t_int_scalar d_nlocal;
  HAT::t_int_scalar h_nlocal;

  DAT::t_int_scalar d_maxelectron;
  HAT::t_int_scalar h_maxelectron;

  DAT::tdual_int_1d k_dellist;
  DAT::t_int_1d d_dellist;

  DAT::t_float_2d d_recomb_ijflag;

  DAT::t_int_2d d_nn_last_partner;

  template < int NEARCP, int GASTALLY > void collisions_one(COLLIDE_REDUCE&);
  template < int DIM, int GASTALLY > void collisions_one_subcell(COLLIDE_REDUCE&);
  template < int GASTALLY > void collisions_one_ambipolar(COLLIDE_REDUCE&);
  template < int NEARCP, int GASTALLY > void collisions_group(COLLIDE_REDUCE&);
  template < int GASTALLY > void collisions_group_ambipolar(COLLIDE_REDUCE&);
  int egroup;        // mixture group containing the ambipolar electron species

  // transient subcell method, per-cell scratch indexed by (icell,index)
  // subcell_id/next indexed by particle; count/first/ring indexed by subcell

  DAT::t_int_2d d_subcell_id;
  DAT::t_int_2d d_subcell_count;
  DAT::t_int_2d d_subcell_first;
  DAT::t_int_2d d_subcell_next;
  DAT::t_int_2d d_subcell_ring;

  template < int DIM >
  KOKKOS_INLINE_FUNCTION
  void rebin_subcell(int, int, int, const double *, const double *) const;

  template < int DIM >
  KOKKOS_INLINE_FUNCTION
  void bin_one_subcell(int, int, int, const double *, const double *) const;

  KOKKOS_INLINE_FUNCTION
  void unbin_one_subcell(int, int, int) const;

  template < int DIM >
  KOKKOS_INLINE_FUNCTION
  int find_nn_subcell(rand_type &, int, int, int, int, int) const;

  void grow_subcell_views(int, int);

  // VSS specific

  DAT::tdual_float_2d k_prefactor;
  DAT::t_float_2d d_prefactor;

  tdual_params_2d k_params;
  t_params_2d d_params;

  double dt,fnum,boltz;
  int maxcellcount,react_defined;

  KOKKOS_INLINE_FUNCTION
  void SCATTER_TwoBodyScattering(Particle::OnePart *,
                                 Particle::OnePart *,
                                 struct State &, struct State &, rand_type &) const;
  KOKKOS_INLINE_FUNCTION
  void EEXCHANGE_NonReactingEDisposal(int icell, Particle::OnePart *,
                                      Particle::OnePart *,
                                      struct State &, struct State &, rand_type &) const;

  KOKKOS_INLINE_FUNCTION
  void SCATTER_ThreeBodyScattering(Particle::OnePart *,
                                   Particle::OnePart *,
                                   Particle::OnePart *,
                                   struct State &, struct State &, rand_type &) const;
  KOKKOS_INLINE_FUNCTION
  void EEXCHANGE_ReactingEDisposal(int, Particle::OnePart *,
                                   Particle::OnePart *,
                                   Particle::OnePart *,
                                   struct State &, struct State &, rand_type &) const;

  KOKKOS_INLINE_FUNCTION
  void relax_electronic_mode(int, Particle::OnePart *, Particle::OnePart *,
                             double&, double, rand_type &rand_gen, bool) const;
  KOKKOS_INLINE_FUNCTION
  void zero_elec(Particle::OnePart *) const;
  KOKKOS_INLINE_FUNCTION
  int elec_exchange(Particle::OnePart *, Particle::OnePart *) const;
  KOKKOS_INLINE_FUNCTION
  double get_elec_phi(int, int, int, double) const;
  KOKKOS_INLINE_FUNCTION
  int select_elec_state(int, Particle::OnePart *, Particle::OnePart *,
                        double, double, bool, rand_type &rand_gen, bool) const;

  KOKKOS_INLINE_FUNCTION
  double sample_bl(rand_type &, double, double) const;
  KOKKOS_INLINE_FUNCTION
  double eff_vib_dof(double, double) const;
  KOKKOS_INLINE_FUNCTION
  double vib_pool_temp(double, int, double *, double) const;
  KOKKOS_INLINE_FUNCTION
  double rotrel (int, double) const;
  KOKKOS_INLINE_FUNCTION
  double vibrel (int, double) const;

  KOKKOS_INLINE_FUNCTION
  int set_nn(int, int) const;
  KOKKOS_INLINE_FUNCTION
  int find_nn(rand_type &, int, int, int) const;

  void grow_group_lists();

  KOKKOS_INLINE_FUNCTION
  int find_nn_group(rand_type &, int, int, int, int, int, int) const;

  void backup();
  void restore();

  t_particle_1d d_particles_backup;
  DAT::t_int_2d d_plist_backup;
  DAT::t_float_3d d_vremax_backup;
  DAT::t_float_3d d_remain_backup;
  DAT::t_int_2d d_nn_last_partner_backup;
  DAT::t_int_1d d_ionambi_backup;
  DAT::t_float_2d_lr d_velambi_backup;
  DAT::t_int_2d d_vibmode_backup;
  DAT::t_float_1d d_eelec_backup;
  DAT::t_int_1d d_elecstate_backup;
  RanKnuth* random_backup;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
