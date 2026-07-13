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

#include "math.h"
#include "string.h"
#include "collide_vss_kokkos.h"
#include "collide_vss_kokkos_inline_impl.h"
#include "collide.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "update.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"
#include "sparta_masks.h"
#include "math_const.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum { SWPM_ENERGY=0, SWPM_HEAT=1, SWPM_STRESS=2 };

#define DELTADELETE 1024
#define DELTACELLCOUNT 2
#define SWPM_SMALL 1.0e-16

/* ---------------------------------------------------------------------- */

template < int NEARCP, int GASTALLY >
void CollideVSSKokkos::collisions_one_swpm(COLLIDE_REDUCE &reduce)
{
  this->sync(Device,ALL_MASK);

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK|CUSTOM_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species = particle_kk->k_species.view_device();

  // get index into edvec for the stochastic weight custom attribute
  index_sw_d = particle->ewhich[index_stochastic_weight];
  d_sw = particle_kk->k_edvec.view_host()[index_sw_d].k_view.view_device();

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_plist = grid_kk->d_plist;

  // copy SWPM params to device-accessible class members
  wtf_d = wtf;
  Ncmin_d = Ncmin;
  Ncmax_d = Ncmax;
  Ngmin_d = Ngmin;
  Ngmax_d = Ngmax;
  reduceflag_d = reduceflag;
  reduction_type_d = reduction_type;

  dt = update->dt;
  fnum = update->fnum;
  boltz = update->boltz;

  copymode = 1;

  if (NEARCP) {
    if (int(d_nn_last_partner.extent(0)) < nglocal ||
        int(d_nn_last_partner.extent(1)) < (int)d_plist.extent(1))
      MemKK::realloc_kokkos(d_nn_last_partner,"collide:nn_last_partner",
                            nglocal,d_plist.extent(1));
  }

  // pre-allocate with headroom for particles created by split
  int maxcellcount = particle_kk->get_maxcellcount();

  auto nlocal_extra = (bigint)(particle->nlocal) * 2 + 1024;
  if ((bigint)d_particles.extent(0) < nlocal_extra) {
    particle->grow((int)(nlocal_extra - particle->nlocal));
    d_particles = particle_kk->k_particles.view_device();
    d_sw = particle_kk->k_edvec.view_host()[index_sw_d].k_view.view_device();
  }

  if ((int)d_plist.extent(1) < maxcellcount + 4) {
    d_plist = {};
    Kokkos::resize(grid_kk->d_plist, nglocal, maxcellcount + 4);
    d_plist = grid_kk->d_plist;
    if (NEARCP)
      MemKK::realloc_kokkos(d_nn_last_partner,"collide:nn_last_partner",
                            nglocal, d_plist.extent(1));
  }

  auto maxdelete_extra = maxdelete + DELTADELETE;
  if ((int)d_dellist.extent(0) < maxdelete_extra) {
    memoryKK->destroy_kokkos(k_dellist,dellist);
    memoryKK->create_kokkos(k_dellist,dellist,maxdelete_extra,"collide:dellist");
    d_dellist = k_dellist.view_device();
  }

  h_retry() = 1;

  while (h_retry()) {

    // back up all mutable state so a buffer-overflow retry is idempotent.
    // unlike the non-SWPM path this is unconditional: split (particle
    // creation) and reduction (particle deletion + weight/velocity edits)
    // run without reactions, so retries are expected and must fully revert.
    backup();

    h_retry() = 0;
    h_maxdelete() = maxdelete;
    h_maxcellcount() = maxcellcount;
    h_part_grow() = 0;
    h_ndelete() = 0;
    h_nlocal() = particle->nlocal;

    if (sparta->kokkos->atomic_reduction) {
      h_nattempt_one() = 0;
      h_ncollide_one() = 0;
    }

    Kokkos::deep_copy(d_scalars, h_scalars);

    grid_kk_copy.copy(grid_kk);

    if (sparta->kokkos->atomic_reduction) {
      if (sparta->kokkos->need_atomics)
        Kokkos::parallel_for(
          Kokkos::RangePolicy<DeviceType, TagCollideSWPMOne<NEARCP,GASTALLY,1>>(0,nglocal),
          *this);
      else
        Kokkos::parallel_for(
          Kokkos::RangePolicy<DeviceType, TagCollideSWPMOne<NEARCP,GASTALLY,0>>(0,nglocal),
          *this);
    } else
      Kokkos::parallel_reduce(
        Kokkos::RangePolicy<DeviceType, TagCollideSWPMOne<NEARCP,GASTALLY,-1>>(0,nglocal),
        *this, reduce);

    Kokkos::deep_copy(h_scalars, d_scalars);

    if (h_retry()) {
      restore();
      reduce = COLLIDE_REDUCE();

      maxdelete = h_maxdelete();
      if ((int)d_dellist.extent(0) < maxdelete) {
        memoryKK->destroy_kokkos(k_dellist,dellist);
        memoryKK->grow_kokkos(k_dellist,dellist,maxdelete,"collide:dellist");
        d_dellist = k_dellist.view_device();
      }

      maxcellcount = h_maxcellcount();
      particle_kk->set_maxcellcount(maxcellcount);
      if ((int)d_plist.extent(1) < maxcellcount) {
        d_plist = {};
        Kokkos::resize(grid_kk->d_plist, nglocal, maxcellcount);
        d_plist = grid_kk->d_plist;
      }

      auto nlocal_new = h_nlocal();
      if ((int)d_particles.extent(0) < nlocal_new) {
        particle->grow(nlocal_new - particle->nlocal);
        d_particles = particle_kk->k_particles.view_device();
        d_sw = particle_kk->k_edvec.view_host()[index_sw_d].k_view.view_device();
      }
    }
  }

  ndelete = h_ndelete();
  particle->nlocal = h_nlocal();

  // counter accumulation and deletion of reduced particles (via dellist)
  // are performed by the shared CollideVSSKokkos::collisions() after this
  // returns; doing them here as well would double-compress the particle list

  copymode = 0;

  if (h_error_flag())
    error->one(FLERR,"Collision cell volume is zero");

  // mark custom (stochastic weight) and particle data modified on device
  particle_kk->k_edvec.view_host()[index_sw_d].k_view.modify_device();
  particle_kk->modify(Device, PARTICLE_MASK|CUSTOM_MASK);

  // SWPM creates/deletes particles without reactions, so the particle list is
  // no longer sorted; collisions() only clears the sort flag when react is set
  particle->sorted = 0;
  particle_kk->sorted_kk = 0;

  d_particles = t_particle_1d();
  d_plist = {};
}

/* ---------------------------------------------------------------------- */

template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideSWPMOne<NEARCP,GASTALLY,ATOMIC_REDUCTION>,
                                   const int &icell) const {
  COLLIDE_REDUCE reduce;
  this->template operator()<NEARCP,GASTALLY,ATOMIC_REDUCTION>(
    TagCollideSWPMOne<NEARCP,GASTALLY,ATOMIC_REDUCTION>(), icell, reduce);
}

template < int NEARCP, int GASTALLY, int ATOMIC_REDUCTION >
KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::operator()(TagCollideSWPMOne<NEARCP,GASTALLY,ATOMIC_REDUCTION>,
                                   const int &icell, COLLIDE_REDUCE &reduce) const
{
  if (d_retry()) return;

  int np = grid_kk_copy.obj.d_cellcount[icell];
  if (np <= 1) return;

  const double volume = grid_kk_copy.obj.k_cinfo.view_device()[icell].volume /
                        grid_kk_copy.obj.k_cinfo.view_device()[icell].weight;
  if (volume == 0.0) { d_error_flag() = 1; return; }

  rand_type rand_gen = rand_pool.get_state();

  // per-cell pre_wtf: active only when cell is below Ncmin threshold
  const double pre_wtf_cell = (Ncmin_d > 0 && np < Ncmin_d) ? 1.0 : 0.0;

  // find max stochastic weight in this cell
  double max_sw = 0.0;
  for (int k = 0; k < np; k++) {
    double sw = d_sw(d_plist(icell,k));
    if (sw > max_sw) max_sw = sw;
  }

  // NTC attempt count with SWPM-modified effective fnum
  const double fnum_eff = max_sw * fnum * (1.0 + pre_wtf_cell * wtf_d);
  double attempt_d;
  if (remainflag) {
    attempt_d = 0.5*np*(np-1)*d_vremax(icell,0,0)*dt*fnum_eff/volume
                + d_remain(icell,0,0);
    d_remain(icell,0,0) = attempt_d - static_cast<int>(attempt_d);
  } else {
    attempt_d = 0.5*np*(np-1)*d_vremax(icell,0,0)*dt*fnum_eff/volume
                + rand_gen.drand();
  }
  const int nattempt = static_cast<int>(attempt_d);
  if (!nattempt) { rand_pool.free_state(rand_gen); return; }

  if (ATOMIC_REDUCTION == 1)
    Kokkos::atomic_add(&d_nattempt_one(), nattempt);
  else if (ATOMIC_REDUCTION == 0)
    d_nattempt_one() += nattempt;
  else
    reduce.nattempt_one += nattempt;

  struct State precoln, postcoln;

  for (int iattempt = 0; iattempt < nattempt; iattempt++) {

    // SWPM uses random partner selection regardless of NEARCP: the
    // reference Collide::collisions_one_stochastic_weighting ignores the
    // NEARCP template parameter and always pairs randomly (nearest-neighbor
    // selection is not meaningful once split clones particles at identical
    // positions).  d_nn_last_partner is therefore never used here.
    int i = (int)(np * rand_gen.drand());
    int j = (int)(np * rand_gen.drand());
    while (i == j) j = (int)(np * rand_gen.drand());

    Particle::OnePart *ipart = &d_particles[d_plist(icell,i)];
    Particle::OnePart *jpart = &d_particles[d_plist(icell,j)];

    // SWPM rejection factor: the NTC attempt count above is inflated by
    // max_sw (fnum_eff), so each pair must be down-selected by
    // MAX(isw,jsw)/max_sw to recover the correct per-pair collision rate
    // (matches CollideVSS::test_collision with stochastic_weight_flag)
    const double isw = d_sw(d_plist(icell,i));
    const double jsw = d_sw(d_plist(icell,j));
    const double ijsw = (isw > jsw ? isw : jsw) / max_sw;

    if (!test_collision_kokkos(icell,0,0,ipart,jpart,precoln,rand_gen,ijsw)) continue;

    // split: adjust weights, create remnant particles
    split_kokkos(icell, np, ipart, jpart, pre_wtf_cell, rand_gen);
    if (d_retry()) { rand_pool.free_state(rand_gen); return; }

    // perform the actual VSS collision (reactions are disabled for SWPM)
    Particle::OnePart *kpart = NULL;
    Particle::OnePart *p3 = NULL;
    int index_kpart = -1;
    double recomb_density = 0.0;
    int recomb_species = -1;
    setup_collision_kokkos(ipart, jpart, precoln, postcoln);
    perform_collision_kokkos(ipart, jpart, kpart, precoln, postcoln, rand_gen,
                             p3, recomb_species, recomb_density, index_kpart);

    if (ATOMIC_REDUCTION == 1)
      Kokkos::atomic_inc(&d_ncollide_one());
    else if (ATOMIC_REDUCTION == 0)
      d_ncollide_one()++;
    else
      reduce.ncollide_one++;
  }

  // particle reduction
  if (reduceflag_d) {
    group_swpm_kokkos(icell, np, rand_gen);
    if (d_retry()) { rand_pool.free_state(rand_gen); return; }
  }

  rand_pool.free_state(rand_gen);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::split_kokkos(int icell, int &np,
    Particle::OnePart *&ip, Particle::OnePart *&jp,
    double pre_wtf_cell, rand_type &rand_gen) const
{
  const int ip_idx = (int)(ip - d_particles.data());
  const int jp_idx = (int)(jp - d_particles.data());
  const double isw = d_sw(ip_idx);
  const double jsw = d_sw(jp_idx);

  if (isw <= 0.0 || jsw <= 0.0) { d_error_flag() = 1; return; }

  double Gwtf, ksw, lsw;
  int ks, ls, kcell, lcell;
  double xk[3], vk[3], xl[3], vl[3];
  double erotk, erotl, evibk, evibl;

  const double divisor = 1.0 + pre_wtf_cell * wtf_d;

  if (isw >= jsw) {
    Gwtf = jsw / divisor;
    ksw  = isw - Gwtf;
    lsw  = jsw - Gwtf;
    ks = ip->ispecies;  ls = jp->ispecies;
    kcell = ip->icell;  lcell = jp->icell;
    for (int d = 0; d < 3; d++) { xk[d]=ip->x[d]; vk[d]=ip->v[d]; xl[d]=jp->x[d]; vl[d]=jp->v[d]; }
    erotk=ip->erot; erotl=jp->erot; evibk=ip->evib; evibl=jp->evib;
  } else {
    Gwtf = isw / divisor;
    ksw  = jsw - Gwtf;
    lsw  = isw - Gwtf;
    ks = jp->ispecies;  ls = ip->ispecies;
    kcell = jp->icell;  lcell = ip->icell;
    for (int d = 0; d < 3; d++) { xk[d]=jp->x[d]; vk[d]=jp->v[d]; xl[d]=ip->x[d]; vl[d]=ip->v[d]; }
    erotk=jp->erot; erotl=ip->erot; evibk=jp->evib; evibl=ip->evib;
  }

  // set both to equal weight Gwtf so the subsequent collision is valid
  d_sw(ip_idx) = Gwtf;
  d_sw(jp_idx) = Gwtf;

  if (Gwtf <= 0.0) { d_error_flag() = 1; return; }

  // create remnant k (copy of the heavier parent with remaining weight)
  if (ksw > 0.0) {
    const int id = (int)(MAXSMALLINT * rand_gen.drand());
    const int slot = Kokkos::atomic_fetch_add(&d_nlocal(), 1);
    const int reallocflag =
      ParticleKokkos::add_particle_kokkos(d_particles, slot, id, ks, kcell, xk, vk, erotk, evibk);
    if (reallocflag) { d_retry() = 1; d_part_grow() = 1; return; }
    d_sw(slot) = ksw;
    if (np < (int)d_plist.extent(1)) {
      d_plist(icell, np++) = slot;
    } else {
      d_retry() = 1;
      d_maxcellcount() += DELTACELLCOUNT;
      return;
    }
    // refresh pointers after plist update (d_particles base is stable within retry cycle)
    ip = &d_particles[ip_idx];
    jp = &d_particles[jp_idx];
  }

  // create remnant l (copy of the lighter parent with remaining weight)
  if (lsw > 0.0) {
    const int id = (int)(MAXSMALLINT * rand_gen.drand());
    const int slot = Kokkos::atomic_fetch_add(&d_nlocal(), 1);
    const int reallocflag =
      ParticleKokkos::add_particle_kokkos(d_particles, slot, id, ls, lcell, xl, vl, erotl, evibl);
    if (reallocflag) { d_retry() = 1; d_part_grow() = 1; return; }
    d_sw(slot) = lsw;
    if (np < (int)d_plist.extent(1)) {
      d_plist(icell, np++) = slot;
    } else {
      d_retry() = 1;
      d_maxcellcount() += DELTACELLCOUNT;
      return;
    }
    ip = &d_particles[ip_idx];
    jp = &d_particles[jp_idx];
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::group_swpm_kokkos(int icell, int &np, rand_type &rand_gen) const
{
  if (np <= Ncmax_d) return;

  // compact d_plist: keep only positive-weight particles
  int n = 0;
  for (int k = 0; k < np; k++) {
    const int pidx = d_plist(icell, k);
    if (d_sw(pidx) > 0.0) d_plist(icell, n++) = pidx;
  }
  np = n;

  int group_size_buffer = 0;

  while (n > Ncmax_d) {
    const int nold = n;

    group_bt_kokkos(icell, 0, n, group_size_buffer, rand_gen);
    if (d_retry()) return;

    // recompact after reduction (deleted particles have sw = -1)
    int nnew = 0;
    for (int k = 0; k < n; k++) {
      const int pidx = d_plist(icell, k);
      if (d_sw(pidx) > 0.0) d_plist(icell, nnew++) = pidx;
    }
    n = nnew;
    np = n;

    if (n == nold) group_size_buffer += 2;
    if (group_size_buffer > n) break;
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_FUNCTION
void CollideVSSKokkos::group_bt_kokkos(int icell, int istart, int iend,
                                        int gsb, rand_type &rand_gen) const
{
  const int np = iend - istart;
  if (np <= Ngmin_d) return;

  // --- first pass: mass-weighted accumulation ---
  double gsum = 0.0, msum = 0.0, Erot = 0.0, Evib = 0.0;
  double mV[3] = {0.0, 0.0, 0.0};

  for (int p = istart; p < iend; p++) {
    const int pidx = d_plist(icell, p);
    const Particle::OnePart *part = &d_particles[pidx];
    const double mass = d_species[part->ispecies].mass;
    const double psw  = d_sw(pidx);
    const double pmsw = psw * mass;
    gsum += psw;
    msum += pmsw;
    Erot += psw * part->erot;
    Evib += psw * part->evib;
    mV[0] += pmsw * part->v[0];
    mV[1] += pmsw * part->v[1];
    mV[2] += pmsw * part->v[2];
  }

  double V[3];
  V[0] = mV[0]/msum; V[1] = mV[1]/msum; V[2] = mV[2]/msum;

  const bool is_leaf = (np <= Ngmax_d + gsb);
  const bool need_q  = is_leaf && (reduction_type_d != SWPM_ENERGY);

  double pij[3][3];
  double q[3] = {0.0, 0.0, 0.0};
  for (int d = 0; d < 3; d++)
    for (int e = 0; e < 3; e++) pij[d][e] = 0.0;

  for (int p = istart; p < iend; p++) {
    const int pidx = d_plist(icell, p);
    const Particle::OnePart *part = &d_particles[pidx];
    const double mass = d_species[part->ispecies].mass;
    const double psw  = d_sw(pidx);
    const double pmsw = psw * mass;
    double dv[3];
    dv[0]=part->v[0]-V[0]; dv[1]=part->v[1]-V[1]; dv[2]=part->v[2]-V[2];
    if (need_q) {
      const double dvsq = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
      for (int d = 0; d < 3; d++) {
        q[d] += 0.5*pmsw*dv[d]*dvsq;
        for (int e = 0; e < 3; e++) pij[d][e] += pmsw*dv[d]*dv[e];
      }
    } else {
      for (int d = 0; d < 3; d++)
        for (int e = 0; e < 3; e++) pij[d][e] += pmsw*dv[d]*dv[e];
    }
  }

  if (is_leaf) {
    double T = (pij[0][0]+pij[1][1]+pij[2][2])/(3.0*gsum*boltz);
    const double mbar = msum/gsum;
    T *= boltz/mbar;
    for (int d = 0; d < 3; d++) {
      q[d] /= mbar;
      for (int e = 0; e < 3; e++) pij[d][e] /= mbar;
    }
    if (reduction_type_d == SWPM_ENERGY)
      reduce_energy_kokkos(icell, istart, iend, gsum, V, T, Erot, Evib, rand_gen);
    else if (reduction_type_d == SWPM_HEAT)
      reduce_heat_kokkos(icell, istart, iend, gsum, V, T, Erot, Evib, q, rand_gen);
    else
      reduce_stress_kokkos(icell, istart, iend, gsum, V, T, Erot, Evib, q, pij, rand_gen);
  } else {
    // compute covariance for eigenvector split
    double Rij[3][3];
    for (int d = 0; d < 3; d++)
      for (int e = 0; e < 3; e++) Rij[d][e] = pij[d][e]/gsum;

    double eval[3], evec[3][3];
    jacobi3_kokkos(Rij, eval, evec);

    // find dominant eigenvector (largest absolute eigenvalue)
    double maxeval = -1.0;
    double maxevec[3] = {0.0, 0.0, 1.0};
    for (int i = 0; i < 3; i++) {
      const double ae = eval[i] < 0.0 ? -eval[i] : eval[i];
      if (ae > maxeval) {
        maxeval = ae;
        maxevec[0]=evec[0][i]; maxevec[1]=evec[1][i]; maxevec[2]=evec[2][i];
      }
    }

    const double center = V[0]*maxevec[0] + V[1]*maxevec[1] + V[2]*maxevec[2];

    // insertion sort by signed projection of velocity onto dominant eigenvector
    for (int k = istart+1; k < iend; k++) {
      const int key = d_plist(icell, k);
      const Particle::OnePart *pk = &d_particles[key];
      const double dk = pk->v[0]*maxevec[0]+pk->v[1]*maxevec[1]+pk->v[2]*maxevec[2]-center;
      int m = k-1;
      while (m >= istart) {
        const Particle::OnePart *pm = &d_particles[d_plist(icell,m)];
        const double dm = pm->v[0]*maxevec[0]+pm->v[1]*maxevec[1]+pm->v[2]*maxevec[2]-center;
        if (dm <= dk) break;
        d_plist(icell, m+1) = d_plist(icell, m);
        m--;
      }
      d_plist(icell, m+1) = key;
    }

    // find split point at median (velocity projection changes sign)
    int cp = 0;
    for (int k = istart; k < iend; k++) {
      const Particle::OnePart *pk = &d_particles[d_plist(icell,k)];
      const double dist = pk->v[0]*maxevec[0]+pk->v[1]*maxevec[1]+pk->v[2]*maxevec[2]-center;
      if (dist < 0.0) cp++;
      else break;
    }
    if (cp == 0 || cp == np) cp = np/2;

    if (cp > Ngmin_d)
      group_bt_kokkos(icell, istart, istart+cp, gsb, rand_gen);
    if (d_retry()) return;
    if (np-cp > Ngmin_d)
      group_bt_kokkos(icell, istart+cp, iend, gsb, rand_gen);
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::jacobi3_kokkos(double m[3][3], double eval[3], double evec[3][3])
{
  double M[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) M[i][j] = m[i][j];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) evec[i][j] = (i==j) ? 1.0 : 0.0;

  // cyclic off-diagonal pairs
  const int pi[3] = {0, 0, 1};
  const int pj[3] = {1, 2, 2};

  for (int sweep = 0; sweep < 50; sweep++) {
    for (int p = 0; p < 3; p++) {
      const int ii = pi[p], jj = pj[p];
      const double Mij = M[ii][jj];
      if (Mij == 0.0) continue;
      const double Mii = M[ii][ii], Mjj = M[jj][jj];
      const double theta = 0.5*(Mjj - Mii) / Mij;
      const double t = (theta >= 0.0)
        ?  1.0/(theta + sqrt(1.0 + theta*theta))
        :  1.0/(theta - sqrt(1.0 + theta*theta));
      const double c = 1.0/sqrt(1.0 + t*t);
      const double s = t*c;
      const double tau = s/(1.0+c);

      M[ii][ii] -= t*Mij;
      M[jj][jj] += t*Mij;
      M[ii][jj] = 0.0;
      M[jj][ii] = 0.0;

      for (int k = 0; k < 3; k++) {
        if (k == ii || k == jj) continue;
        const double Mki = M[k][ii], Mkj = M[k][jj];
        M[k][ii] = M[ii][k] = Mki - s*(Mkj + tau*Mki);
        M[k][jj] = M[jj][k] = Mkj + s*(Mki - tau*Mkj);
      }
      for (int k = 0; k < 3; k++) {
        const double eki = evec[k][ii], ekj = evec[k][jj];
        evec[k][ii] = eki - s*(ekj + tau*eki);
        evec[k][jj] = ekj + s*(eki - tau*ekj);
      }
    }
  }
  eval[0] = M[0][0]; eval[1] = M[1][1]; eval[2] = M[2][2];
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::sample_unit_sphere_kokkos(rand_type &rand_gen, double *uvec)
{
  const double theta = 2.0*MY_PI*rand_gen.drand();
  const double phi   = acos(1.0 - 2.0*rand_gen.drand());
  uvec[0] = sin(phi)*cos(theta);
  uvec[1] = sin(phi)*sin(theta);
  uvec[2] = cos(phi);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::reduce_delete_kokkos(int particle_idx) const
{
  d_sw(particle_idx) = -1.0;
  const int slot = Kokkos::atomic_fetch_add(&d_ndelete(), 1);
  if (slot < (int)d_dellist.extent(0)) {
    d_dellist(slot) = particle_idx;
  } else {
    d_retry() = 1;
    d_maxdelete() += DELTADELETE;
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::reduce_energy_kokkos(int icell, int istart, int iend,
    double rho, double *V, double T, double Erot, double Evib,
    rand_type &rand_gen) const
{
  const int np = iend - istart;
  int ip = (int)(np*rand_gen.drand()) + istart;
  int jp = (int)(np*rand_gen.drand()) + istart;
  while (ip == jp) jp = (int)(np*rand_gen.drand()) + istart;

  Particle::OnePart *ipart = &d_particles[d_plist(icell,ip)];
  Particle::OnePart *jpart = &d_particles[d_plist(icell,jp)];

  double uvec[3];
  sample_unit_sphere_kokkos(rand_gen, uvec);

  const double sqT = sqrt(3.0*T);
  for (int d = 0; d < 3; d++) {
    ipart->v[d] = V[d] + sqT*uvec[d];
    jpart->v[d] = V[d] - sqT*uvec[d];
  }
  ipart->erot = Erot/rho; jpart->erot = Erot/rho;
  ipart->evib = Evib/rho; jpart->evib = Evib/rho;

  d_sw(d_plist(icell,ip)) = rho*0.5;
  d_sw(d_plist(icell,jp)) = rho*0.5;

  for (int i = istart; i < iend; i++) {
    if (i == ip || i == jp) continue;
    reduce_delete_kokkos(d_plist(icell, i));
    if (d_retry()) return;
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::reduce_heat_kokkos(int icell, int istart, int iend,
    double rho, double *V, double T, double Erot, double Evib,
    double *q, rand_type &rand_gen) const
{
  const int np = iend - istart;
  int ip = (int)(np*rand_gen.drand()) + istart;
  int jp = (int)(np*rand_gen.drand()) + istart;
  while (ip == jp) jp = (int)(np*rand_gen.drand()) + istart;

  Particle::OnePart *ipart = &d_particles[d_plist(icell,ip)];
  Particle::OnePart *jpart = &d_particles[d_plist(icell,jp)];

  const double sqT = sqrt(3.0*T);
  const double qmag = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
  double qge;
  if (sqT < SWPM_SMALL) qge = 0.0;
  else qge = qmag/(rho*(sqT*sqT*sqT));
  if (qge*qge > 1.0/SWPM_SMALL) qge = 0.0;
  const double itheta = qge + sqrt(1.0+qge*qge);
  const double alpha  = sqT*itheta;
  const double beta   = sqT/itheta;

  double uvec[3];
  if (qmag < SWPM_SMALL) sample_unit_sphere_kokkos(rand_gen, uvec);
  else { uvec[0]=q[0]/qmag; uvec[1]=q[1]/qmag; uvec[2]=q[2]/qmag; }

  for (int d = 0; d < 3; d++) {
    ipart->v[d] = V[d] + alpha*uvec[d];
    jpart->v[d] = V[d] - beta*uvec[d];
  }

  const double isw = rho/(1.0+itheta*itheta);
  const double jsw = rho - isw;

  ipart->erot = Erot/isw*0.5; jpart->erot = Erot/jsw*0.5;
  ipart->evib = Evib/isw*0.5; jpart->evib = Evib/jsw*0.5;

  d_sw(d_plist(icell,ip)) = isw;
  d_sw(d_plist(icell,jp)) = jsw;

  for (int i = istart; i < iend; i++) {
    if (i == ip || i == jp) continue;
    reduce_delete_kokkos(d_plist(icell, i));
    if (d_retry()) return;
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::reduce_stress_kokkos(int icell, int istart, int iend,
    double rho, double *V, double T, double Erot, double Evib,
    double *q, double pij[3][3], rand_type &rand_gen) const
{
  double pij_local[3][3];
  for (int d = 0; d < 3; d++)
    for (int e = 0; e < 3; e++) pij_local[d][e] = pij[d][e];

  double eval[3], evec[3][3];
  jacobi3_kokkos(pij_local, eval, evec);

  const double eval_thresh = pow(SWPM_SMALL, 2.0/3.0);
  int nK = 0;
  for (int i = 0; i < 3; i++) {
    if (eval[i] >= eval_thresh) {
      eval[nK] = eval[i];
      for (int d = 0; d < 3; d++) evec[d][nK] = evec[d][i];
      nK++;
    }
  }

  if (nK == 0) {
    reduce_energy_kokkos(icell, istart, iend, rho, V, T, Erot, Evib, rand_gen);
    return;
  }

  for (int iK = 0; iK < nK; iK++) {
    Particle::OnePart *ipart = &d_particles[d_plist(icell, 2*iK+istart)];
    Particle::OnePart *jpart = &d_particles[d_plist(icell, 2*iK+1+istart)];

    double qli = evec[0][iK]*q[0] + evec[1][iK]*q[1] + evec[2][iK]*q[2];
    if (qli < 0.0) {
      for (int d = 0; d < 3; d++) evec[d][iK] *= -1.0;
      qli = -qli;
    }

    const double sqE = sqrt(eval[iK]);
    double A = (qli < SWPM_SMALL) ? 0.0
               : sqrt(rho)*qli/(sqrt((double)nK)*sqE*eval[iK]);
    if (A*A > 1.0/SWPM_SMALL) A = 0.0;
    const double itheta = A + sqrt(1.0+A*A);
    const double scale  = sqrt((double)nK*eval[iK]/rho);

    for (int d = 0; d < 3; d++) {
      ipart->v[d] = V[d] + itheta*scale*evec[d][iK];
      jpart->v[d] = V[d] - (1.0/itheta)*scale*evec[d][iK];
    }

    const double isw = rho/((double)nK*(1.0+itheta*itheta));
    const double jsw = rho/(double)nK - isw;

    ipart->erot = Erot/isw*0.5/(double)nK;
    jpart->erot = Erot/jsw*0.5/(double)nK;
    ipart->evib = Evib/isw*0.5/(double)nK;
    jpart->evib = Evib/jsw*0.5/(double)nK;

    d_sw(d_plist(icell, 2*iK+istart)) = isw;
    d_sw(d_plist(icell, 2*iK+1+istart)) = jsw;
  }

  for (int i = istart+2*nK; i < iend; i++) {
    reduce_delete_kokkos(d_plist(icell, i));
    if (d_retry()) return;
  }
}

/* ---------------------------------------------------------------------- */

namespace SPARTA_NS {

template void CollideVSSKokkos::collisions_one_swpm<0,0>(COLLIDE_REDUCE&);
template void CollideVSSKokkos::collisions_one_swpm<0,1>(COLLIDE_REDUCE&);
template void CollideVSSKokkos::collisions_one_swpm<1,0>(COLLIDE_REDUCE&);
template void CollideVSSKokkos::collisions_one_swpm<1,1>(COLLIDE_REDUCE&);

}
