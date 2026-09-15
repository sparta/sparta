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

#include "stdlib.h"
#include "string.h"
#include "fix_emit_surf_kokkos.h"
#include "compute_isurf_grid_kokkos.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "grid.h"
#include "kokkos.h"
#include "surf_kokkos.h"
#include "particle.h"
#include "mixture.h"
#include "comm.h"
#include "modify.h"
#include "geometry.h"
#include "input.h"
#include "random_knuth.h"
#include "math_const.h"
#include "math_extra_kokkos.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos_type.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "fix_emit_kokkos.h"
#include "sparta_masks.h"
#include "variable.h"
#include "Kokkos_Random.hpp"

using namespace SPARTA_NS;

#ifndef SPARTA_KOKKOS_FIXED_LISTS
namespace {

  // the surf-tally computes this fix drives are held in one runtime-sized
  //   device byte buffer instead of a fixed KKCopy array, so their number is
  //   not capped.  the three steps mirror KKCopy exactly: size the buffer,
  //   blit each live object's bytes into its slot and mark the image a copy
  //   so its destructor frees nothing, then push the buffer to the device.
  //   this is the same code UpdateKokkos uses for its own tally lists
  //   (update_kokkos.cpp:94-108); it is duplicated rather than shared
  //   because both live in anonymous namespaces in their own translation
  //   units, and the whole point is that neither is a class member

  void slist_buf_grow(DAT::tdual_char_1d &k, int n, size_t bytes)
  {
    const size_t need = (size_t) n * bytes;
    if (k.view_host().extent(0) >= need) return;
    k = DAT::tdual_char_1d("fix_emit_surf:slist_surf",need);
  }

  template<class T>
  void slist_buf_blit(DAT::tdual_char_1d &k, int slot, T *obj)
  {
    char *dst = k.view_host().data() + (size_t) slot*sizeof(T);
    memcpy((void*) dst, (const void*) obj, sizeof(T));
    ((T *) dst)->copy = 1;
  }

  void slist_buf_sync(DAT::tdual_char_1d &k, DAT::t_char_1d &d)
  {
    if (k.view_device().extent(0) == 0) return;
    k.modify_host();
    k.sync_device();
    d = k.view_device();
  }
}
#endif
using namespace MathConst;

enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{NOSUBSONIC,PTBOTH,PONLY,MFLOW};   // same as fix_emit_surf.cpp
enum{FLOW,CONSTANT,VARIABLE};
enum{INT,DOUBLE};                                        // several files

#define DELTATASK 256
#define TEMPLIMIT 1.0e5

#define VAL_1(X) X
#define VAL_2(X) VAL_1(X), VAL_1(X)

/* ----------------------------------------------------------------------
   insert particles in grid cells with surfs touching inflow boundaries
------------------------------------------------------------------------- */

FixEmitSurfKokkos::FixEmitSurfKokkos(SPARTA *sparta, int narg, char **arg) :
  FixEmitSurf(sparta, narg, arg),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            ),
  particle_kk_copy(sparta)
#ifdef SPARTA_KOKKOS_FIXED_LISTS
  , slist_active_copy{VAL_2(KKCopy<ComputeSurfKokkos>(sparta))}
  , tmp_compute_surf_kk(sparta)
#endif
{
  kokkos_flag = 1;
  execution_space = Device;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  region_flag = 0;
  nregion_token = 0;

}

/* ---------------------------------------------------------------------- */

FixEmitSurfKokkos::~FixEmitSurfKokkos()
{
  if (copymode) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
#endif

  // the loop below reads and writes the host side of k_tasks, which the
  // emission kernels leave claimed for the device.  path and fracarea are host
  // allocations the device never writes, so the host values are the right ones
  // and no copy back is wanted here; say so rather than reading a side the
  // coherence state calls stale.

  k_tasks.clear_sync_state();

  for (int i = 0; i < ntaskmax; i++) {
    tasks[i].ntargetsp = NULL;
    tasks[i].vscale = NULL;
    delete [] tasks[i].path;
    delete [] tasks[i].fracarea;
  }

  tasks = NULL;
}

/* ---------------------------------------------------------------------- */

void FixEmitSurfKokkos::init()
{
  k_tasks.sync_host();
  if (perspecies) k_ntargetsp.sync_host();
  if (subsonic_style == PONLY || subsonic_style == MFLOW || temp_custom_flag) k_vscale.sync_host();

  k_path.sync_host();

  if (dimension != 2)
    k_fracarea.sync_host();

  FixEmitSurf::init();

  k_tasks.modify_host();
  if (perspecies) k_ntargetsp.modify_host();
  if (subsonic_style == PONLY || subsonic_style == MFLOW || temp_custom_flag) k_vscale.modify_host();

  k_path.modify_host();

  if (dimension != 2)
    k_fracarea.modify_host();

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif

  k_vscale_mix = DAT::tdual_float_1d("vscale_mix", nspecies);

  if (!(fractions_custom_flag && !perspecies))
    k_cummulative_mix = DAT::tdual_float_1d("cummulative", nspecies);

  k_mspecies = DAT::tdual_int_1d("species", nspecies);
  k_fraction = DAT::tdual_float_1d("fraction", nspecies);

  d_vscale_mix = k_vscale_mix .view_device();
  d_cummulative_mix = k_cummulative_mix.view_device();
  d_mspecies = k_mspecies.view_device();
  d_fraction = k_fraction.view_device();

  auto h_vscale_mix  = k_vscale_mix.view_host();
  auto h_cummulative_mix = k_cummulative_mix.view_host();
  auto h_species = k_mspecies.view_host();
  auto h_fraction = k_fraction.view_host();

  for (int isp = 0; isp < nspecies; ++isp) {
    h_vscale_mix(isp) = particle->mixture[imix]->vscale[isp];

    if (!(fractions_custom_flag && !perspecies))
      h_cummulative_mix(isp) = particle->mixture[imix]->cummulative[isp];

    h_species(isp) = particle->mixture[imix]->species[isp];
    h_fraction(isp) = particle->mixture[imix]->fraction[isp];
  }

  k_vscale_mix.modify_host();

  if (!(fractions_custom_flag && !perspecies))
    k_cummulative_mix.modify_host();

  k_mspecies.modify_host();
  k_fraction.modify_host();
  // the region is fixed for the run; flatten it once here rather than
  //   rebuilding and re-uploading the token stream in perform_task()

  flatten_region();

}

/* ----------------------------------------------------------------------
   grid changed operation
   invoke create_tasks() to rebuild entire task list
   invoked after per-processor list of grid cells has changed
   invoked after custom per-surf attributes have changed (fix custom)
------------------------------------------------------------------------- */

void FixEmitSurfKokkos::grid_changed()
{
  FixEmitSurf::grid_changed();

  // if custom fractions requested and perspecies = 0,
  // setup cummulative_custom array for nlocal surfs

  if (fractions_custom_flag && !perspecies) {
    if (k_cummulative_custom.extent(0) > max_cummulative)
      MemKK::realloc_kokkos(k_cummulative_custom,"fix/emit/surf:cummulative_custom",max_cummulative,nspecies);

    for (int isurf = 0; isurf < max_cummulative; isurf++) {
      for (int isp = 0; isp < nspecies; isp++) {
        k_cummulative_custom.view_host()(isurf,isp) = cummulative_custom[isurf][isp];
      }
    }

    k_cummulative_custom.modify_host();
  }
}

/* ----------------------------------------------------------------------
   create tasks for one grid cell
   add them to tasks list and increment ntasks
------------------------------------------------------------------------- */

void FixEmitSurfKokkos::create_tasks()
{
  k_tasks.sync_host();
  if (perspecies) k_ntargetsp.sync_host();
  if (subsonic_style == PONLY || subsonic_style == MFLOW || temp_custom_flag) k_vscale.sync_host();

  FixEmit::create_tasks();

  k_tasks.modify_host();
  if (perspecies) k_ntargetsp.modify_host();
  if (subsonic_style == PONLY || subsonic_style == MFLOW || temp_custom_flag) k_vscale.modify_host();

  if (ntaskmax > k_path.extent(0) || max_npoint > k_path.extent(1)/3) {
    MemKK::realloc_kokkos(k_path,"fix_emit_surf:path",ntaskmax,max_npoint*3);
    d_path = k_path.view_device();

    if (dimension != 2) {
      MemKK::realloc_kokkos(k_fracarea,"fix_emit_surf:fracarea",ntaskmax,max_npoint-2);
      d_fracarea = k_fracarea.view_device();
    }
  }

  for (int i = 0; i < ntask; i++) {
    auto task_i = tasks[i];
    int npoint = task_i.npoint;

    for (int n = 0; n < npoint*3; n++)
      k_path.view_host()(i,n) = task_i.path[n];

    if (dimension != 2)
      for (int n = 0; n < npoint-2; n++)
        k_fracarea.view_host()(i,n) = task_i.fracarea[n];
  }

  k_path.modify_host();

  if (dimension != 2)
    k_fracarea.modify_host();
}

/* ----------------------------------------------------------------------
   flatten the region to a device-resident postfix token stream, so the
   emit kernel needs no virtual dispatch and no typed copy per region style.
   the stream carries each sub-region's interior/exterior sense and the
   composite's own, so nothing else needs to be passed along.
   see region_prim_kokkos.h

   a region is fixed for the duration of a run -- Region exposes no move or
   rotate, and "region" is an input command -- so this runs from init()
   rather than from perform_task(), which would rebuild the stream on the
   host and re-upload it every step
------------------------------------------------------------------------- */

void FixEmitSurfKokkos::flatten_region()
{
  region_flag = 0;
  nregion_token = 0;
  if (region) {
    KokkosBase* region_kkbase = dynamic_cast<KokkosBase*>(region);
    if (!region->kokkos_flag || !region_kkbase)
      error->all(FLERR,"KOKKOS package does not (yet) support chosen region style");
    nregion_token = region_kkbase->flatten_region_kokkos(k_region_tokens);
    if (nregion_token <= 0)
      error->all(FLERR,"KOKKOS package does not (yet) support chosen region style");
    d_region_tokens = k_region_tokens.view_device();
    region_flag = 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixEmitSurfKokkos::perform_task()
{
  dt = update->dt;

  // if subsonic, re-compute particle inflow counts for each task
  // also computes current per-task temp_thermal and vstream

  if (subsonic) subsonic_inflow();

  // if npmode = VARIABLE, set npcurrent to variable evaluation

  if (npmode == VARIABLE) {
    npcurrent = input->variable->compute_equal(npvar);
    if (npcurrent <= 0.0) error->all(FLERR,"Fix emit/surf Np <= 0.0");
  }

  // insert particles for each task = cell/surf pair
  // ntarget/ninsert is either perspecies or for all species
  // for one particle:
  //   x = random position with overlap of surf with cell
  //   v = randomized thermal velocity + vstream
  //       if normalflag, mag of vstream is applied to surf normal dir
  //       first stage: normal dimension (normal)
  //       second stage: parallel dimensions (tan1,tan2)

  // double while loop until randomized particle velocity meets 2 criteria
  // inner do-while loop:
  //   v = vstream-component + vthermal is into simulation box
  //   see Bird 1994, p 425
  // outer do-while loop:
  //   shift Maxwellian distribution by stream velocity component
  //   see Bird 1994, p 259, eq 12.5

  // copy needed task data to device
  // the kernels below read d_tasks whether or not perspecies is set, so tasks
  //   is synced unconditionally and ntargetsp in addition, the same shape as
  //   the second copy further down this routine.  The if/else this replaces
  //   left tasks unsynced under perspecies and was only harmless because that
  //   later copy covered it

  k_tasks.sync_device();
  if (perspecies) k_ntargetsp.sync_device();

  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,ALL_MASK);
  d_lines = surf_kk->k_lines.view_device();
  d_tris = surf_kk->k_tris.view_device();
  nlocal_surf = surf->nlocal;

  nsurf_tally = update->nsurf_tally;
  Compute **slist_active = update->slist_active;

#ifdef SPARTA_KOKKOS_FIXED_LISTS
  if (nsurf_tally > KOKKOS_MAX_SLIST)
    error->all(FLERR,"Kokkos currently only supports two instances of compute "
               "surface with fix emit/surf; rebuild without "
               "-DSPARTA_KOKKOS_FIXED_LISTS to lift the limit");
#else
  slist_buf_grow(k_slist_surf,nsurf_tally,sizeof(ComputeSurfKokkos));
#endif

  // dispatch by dynamic_cast, not by style string: the Kokkos computes are
  //   registered under their "/kk" names too (isurf/grid/kk et al), so a
  //   style compare never matches a compute the user typed with the suffix
  //   and the caller then falls through to an unrelated error message.
  //   see the same note in UpdateKokkos::setup_surf_tally_copies()

  for (int i = 0; i < nsurf_tally; i++) {
    ComputeSurfKokkos* compute_surf_kk = dynamic_cast<ComputeSurfKokkos*>(slist_active[i]);
    if (!compute_surf_kk) {
      if (dynamic_cast<ComputeISurfGridKokkos*>(slist_active[i]))
        error->all(FLERR,"Kokkos does not (yet) support compute isurf/grid "
                   "with fix emit/surf");
      error->all(FLERR,"Kokkos does not (yet) support this surf tally compute "
                 "with fix emit/surf; use a Kokkos-enabled surf tally compute "
                 "(-sf kk)");
    }
    compute_surf_kk->pre_surf_tally();
#ifdef SPARTA_KOKKOS_FIXED_LISTS
    slist_active_copy[i].copy(compute_surf_kk);
#else
    slist_buf_blit(k_slist_surf,i,compute_surf_kk);
#endif
  }

#ifdef SPARTA_KOKKOS_FIXED_LISTS

  // every Kokkos functor captures the whole array by value, so the unused
  //  slots must not alias a compute that may be reallocated or deleted while
  //  they still reference count it: point them at a temporary that lives as
  //  long as this class

  for (int i = nsurf_tally; i < KOKKOS_MAX_SLIST; i++)
    slist_active_copy[i].copy(&tmp_compute_surf_kk);
#else
  slist_buf_sync(k_slist_surf,d_slist_surf);
#endif

  auto ninsert_dim1 = perspecies ? nspecies : 1;
  if (d_ninsert.extent(0) < ntask * ninsert_dim1)
    d_ninsert = DAT::t_int_1d("ninsert", ntask * ninsert_dim1);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitSurf_ninsert>(0,ntask),*this);
  copymode = 0;

  int ncands;
  d_task2cand = offset_scan(d_ninsert, ncands);

  if (ncands == 0) return;

  if (d_x.extent(0) < ncands || d_x.extent(1) < dimension)
    d_x = DAT::t_float_2d("x", ncands, dimension);

  if (d_v.extent(0) < ncands)
    d_v = DAT::t_float_2d("v", ncands, 3);

  if (d_task.extent(0) < ncands) {
    d_erot     = DAT::t_float_1d("erot", ncands);
    d_evib     = DAT::t_float_1d("evib", ncands);
    d_dtremain = DAT::t_float_1d("dtremain", ncands);
    d_id       = DAT::t_int_1d("id", ncands);
    d_isp      = DAT::t_int_1d("isp", ncands);
    d_task     = DAT::t_int_1d("task", ncands);
    d_keep     = DAT::t_int_1d("keep", ncands);
  }
  Kokkos::deep_copy(d_keep,0); // needs to be initialized with zeros

  // copy needed task data to device

  k_tasks.sync_device();
  if (perspecies) k_ntargetsp.sync_device();
  if (subsonic_style == PONLY || subsonic_style == MFLOW || temp_custom_flag) k_vscale.sync_device();

  k_path.sync_device();

  if (dimension != 2)
    k_fracarea.sync_device();

  // copy needed mixture data to device

  k_vscale_mix .sync_device();
  k_mspecies    .sync_device();

  if (fractions_custom_flag && !perspecies)
    k_cummulative_custom.sync_device();
  else
    k_cummulative_mix.sync_device();

  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);
  particle_kk->update_class_variables();
  particle_kk_copy.copy(particle_kk);

  // flatten the region to a device-resident postfix token stream, so the
  //   kernel below needs no virtual dispatch and no typed copy per region
  //   style.  the stream carries each sub-region's interior/exterior sense
  //   and the composite's own, so nothing else needs to be passed along.
  //   see region_prim_kokkos.h


  int nsingle_reduce = 0;
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEmitSurf_perform_task>(0,ntask),*this,nsingle_reduce);
  copymode = 0;
  nsingle += nsingle_reduce;

  int nnew;
  d_cands2new = offset_scan(d_keep, nnew);

  auto particleKK = dynamic_cast<ParticleKokkos*>(particle);
  nlocal_before = particleKK->nlocal;
  particleKK->grow(nnew);
  particleKK->sync(SPARTA_NS::Device, PARTICLE_MASK);
  d_particles = particleKK->k_particles.view_device();

  /* ATOMIC_REDUCTION: 1 = use atomics
                       0 = don't need atomics
  */

  copymode = 1;
  if (sparta->kokkos->need_atomics)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitSurf_insert_particles<1> >(0,ncands),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitSurf_insert_particles<0> >(0,ncands),*this);
  copymode = 0;

  particleKK->nlocal = nlocal_before + nnew;
  particleKK->zero_custom_kokkos(nlocal_before,particleKK->nlocal);
  particleKK->modify(SPARTA_NS::Device, PARTICLE_MASK);

  if (nsurf_tally) {
    for (int m = 0; m < nsurf_tally; m++) {
      ComputeSurfKokkos* compute_surf_kk = (ComputeSurfKokkos*)(slist_active[m]);
      compute_surf_kk->post_surf_tally();
    }
  }

  // if using per-surf custom attributes,
  // temps/vstream already set to custom attributes in create_task

  if (modify->n_update_custom) {
    auto h_keep = Kokkos::create_mirror_view(d_keep);
    auto h_task = Kokkos::create_mirror_view(d_task);
    Kokkos::deep_copy(h_keep, d_keep);
    Kokkos::deep_copy(h_task, d_task);

    // copy needed task data to host

    k_tasks.sync_host();

    auto h_cands2new = Kokkos::create_mirror_view(d_cands2new);
    Kokkos::deep_copy(h_cands2new, d_cands2new);
    for (int cand = 0; cand < ncands; ++cand) {
      if (!h_keep(cand)) continue;

      auto task = h_task(cand);

      auto temp_thermal = tasks[task].temp_thermal;
      auto temp_rot = tasks[task].temp_rot;
      auto temp_vib = tasks[task].temp_vib;
      auto vstream = tasks[task].vstream;

      auto inew = h_cands2new(cand);
      auto ilocal = nlocal_before + inew;

      modify->update_custom(ilocal,temp_thermal,
          temp_rot,temp_vib,vstream);
    }
  }
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixEmitSurfKokkos::operator()(TagFixEmitSurf_ninsert, const int &i) const
{
  int ninsert;

  rand_type rand_gen = rand_pool.get_state();

  if (perspecies) {
    for (int isp = 0; isp < nspecies; isp++) {
      auto ntarget = d_ntargetsp(i,isp) + rand_gen.drand();
      ninsert = static_cast<int> (ntarget);
      d_ninsert(i * nspecies + isp) = ninsert;
    }
  } else {
    // set ntarget for insertion mode FLOW, CONSTANT, or VARIABLE
    // for FLOW: ntarget is already set within task
    // for CONSTANT or VARIABLE: task narget is fraction of its surf's area
    //   scale fraction by np or npcurrent (variable evaluation)
    // ninsert = rounded-down (ntarget + random number)

    double ntarget;
    auto task_i = d_tasks(i);
    if (npmode == FLOW) ntarget = task_i.ntarget;
    else if (npmode == CONSTANT) ntarget = np * task_i.ntarget;
    else if (npmode == VARIABLE) ntarget = npcurrent * task_i.ntarget;
    if (ntarget >= (double) MAXSMALLINT)
      Kokkos::abort("Fix emit/surf insertion count per task exceeds 2^31");
    ninsert = static_cast<int> (ntarget + rand_gen.drand());
    d_ninsert(i) = ninsert;
  }

  rand_pool.free_state(rand_gen);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixEmitSurfKokkos::operator()(TagFixEmitSurf_perform_task, const int &i, int &nsingle) const
{
  double *vstream,*normal,*atan,*btan;

  rand_type rand_gen = rand_pool.get_state();

  auto task_i = d_tasks(i);

  const double temp_rot = task_i.temp_rot;
  const double temp_vib = task_i.temp_vib;
  const double magvstream = task_i.magvstream;
  vstream = task_i.vstream;

  const surfint isurf = task_i.isurf;
  if (isurf >= nlocal_surf) Kokkos::abort("BAD surf index");
  if (dimension == 2) normal = d_lines[isurf].norm;
  else normal = d_tris[isurf].norm;
  atan = task_i.tan1;
  btan = task_i.tan2;

  double indot;
  if (normalflag) indot = magvstream;
  else indot = vstream[0]*normal[0] + vstream[1]*normal[1] + vstream[2]*normal[2];

  // perspecies yes

  if (perspecies) {
    for (int isp = 0; isp < nspecies; isp++) {

      double vscale = (subsonic_style == PONLY || subsonic_style == MFLOW || temp_custom_flag) ?
        d_vscale(i, isp) : d_vscale_mix(isp);

      const int ispecies = d_mspecies[isp];
      const double ninsert = d_ninsert(i * nspecies + isp);
      const int start = d_task2cand(i * nspecies + isp);
      const double scosine = indot / vscale;

      // loop over ninsert for each species

      int nactual = 0;
      for (int m = 0; m < ninsert; m++) {
        auto cand = start + m;

        double x[3];
        if (dimension == 2) {
          double rn = rand_gen.drand();
          double* p1 = &d_path(i,0);
          double* p2 = &d_path(i,3);
          if (axisymmetric && p1[1] != p2[1])
            rn = (sqrt(p1[1]*p1[1] + rn*(p2[1]*p2[1]-p1[1]*p1[1])) - p1[1]) /
              (p2[1]-p1[1]);
          x[0] = p1[0] + rn * (p2[0]-p1[0]);
          x[1] = p1[1] + rn * (p2[1]-p1[1]);
          x[2] = 0.0;
        } else {
          const double rn = rand_gen.drand();
          int ntri = task_i.npoint - 2;
          int n;
          for (n = 0; n < ntri; n++)
            if (rn < d_fracarea(i,n)) break;
          double* p1 = &d_path(i,0);
          double* p2 = &d_path(i,3*(n+1));
          double* p3 = &d_path(i,3*(n+2));
          double e1[3],e2[3];
          MathExtraKokkos::sub3(p2,p1,e1);
          MathExtraKokkos::sub3(p3,p1,e2);
          double alpha = rand_gen.drand();
          double beta = rand_gen.drand();
          if (alpha+beta > 1.0) {
            alpha = 1.0 - alpha;
            beta = 1.0 - beta;
          }
          x[0] = p1[0] + alpha*e1[0] + beta*e2[0];
          x[1] = p1[1] + alpha*e1[1] + beta*e2[1];
          x[2] = p1[2] + alpha*e1[2] + beta*e2[2];
        }

        if (region_flag &&
            !region_match_kk(d_region_tokens,nregion_token,
                             x[0],x[1],x[2])) continue;

        nactual++;
        d_keep(cand) = 1;
        d_task(cand) = i;
        d_isp(cand) = isp;
        for (int d = 0; d < dimension; ++d) d_x(cand, d) = x[d];

        double beta_un, normalized_distbn_fn;
        do {
          do beta_un = (6.0*rand_gen.drand() - 3.0);
          while (beta_un + scosine < 0.0);
          normalized_distbn_fn = 2.0 * (beta_un + scosine) /
            (scosine + sqrt(scosine*scosine + 2.0)) *
            exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) -
                beta_un*beta_un);
        } while (normalized_distbn_fn < rand_gen.drand());

        double vnmag;
        if (normalflag) vnmag = beta_un*vscale + magvstream;
        else vnmag = beta_un*vscale + indot;

        const double theta = MY_2PI * rand_gen.drand();
        const double vr = vscale * sqrt(-log(rand_gen.drand()));

        double vamag,vbmag;
        if (normalflag) {
          vamag = vr * sin(theta);
          vbmag = vr * cos(theta);
        } else {
          vamag = vr * sin(theta) + MathExtraKokkos::dot3(vstream,atan);
          vbmag = vr * cos(theta) + MathExtraKokkos::dot3(vstream,btan);
        }

        d_v(cand, 0) = vnmag*normal[0] + vamag*atan[0] + vbmag*btan[0];
        d_v(cand, 1) = vnmag*normal[1] + vamag*atan[1] + vbmag*btan[1];
        d_v(cand, 2) = vnmag*normal[2] + vamag*atan[2] + vbmag*btan[2];

        d_erot(cand) = particle_kk_copy.obj.erot(ispecies,temp_rot,rand_gen);
        d_evib(cand) = particle_kk_copy.obj.evib(ispecies,temp_vib,rand_gen);
        d_id(cand) = MAXSMALLINT*rand_gen.drand();
        d_dtremain(cand) = dt * rand_gen.drand();
      }

      nsingle += nactual;
    }

  // perspecies no

  } else {
    int ninsert = d_ninsert(i);
    int start = d_task2cand(i);

    // loop over ninsert for all species
    // use cummulative fractions to assign species for each insertion
    // if requested, override cummulative from mixture with cummulative for isurf

    double* cummulative = &d_cummulative_mix[0];
    if (fractions_custom_flag) cummulative = &d_cummulative_custom(isurf,0);

    int nactual = 0;
    for (int m = 0; m < ninsert; m++) {
      const int cand = start + m;
      const double rn = rand_gen.drand();
      int isp = 0;
      while (cummulative[isp] < rn) isp++;

      double vscale = (subsonic_style == PONLY || subsonic_style == MFLOW || temp_custom_flag) ?
        d_vscale(i, isp) : d_vscale_mix(isp);

      const int ispecies = d_mspecies[isp];
      const double scosine = indot / vscale;

      double x[3];
      if (dimension == 2) {
        double rn = rand_gen.drand();
        double* p1 = &d_path(i,0);
        double* p2 = &d_path(i,3);
        if (axisymmetric && p1[1] != p2[1])
          rn = (sqrt(p1[1]*p1[1] + rn*(p2[1]*p2[1]-p1[1]*p1[1])) - p1[1]) /
            (p2[1]-p1[1]);
        x[0] = p1[0] + rn * (p2[0]-p1[0]);
        x[1] = p1[1] + rn * (p2[1]-p1[1]);
        x[2] = 0.0;
      } else {
        const double rn = rand_gen.drand();
        int ntri = task_i.npoint - 2;
        int n;
        for (n = 0; n < ntri; n++)
          if (rn < d_fracarea(i,n)) break;
        double* p1 = &d_path(i,0);
        double* p2 = &d_path(i,3*(n+1));
        double* p3 = &d_path(i,3*(n+2));
        double e1[3],e2[3];
        MathExtraKokkos::sub3(p2,p1,e1);
        MathExtraKokkos::sub3(p3,p1,e2);
        double alpha = rand_gen.drand();
        double beta = rand_gen.drand();
        if (alpha+beta > 1.0) {
          alpha = 1.0 - alpha;
          beta = 1.0 - beta;
        }
        x[0] = p1[0] + alpha*e1[0] + beta*e2[0];
        x[1] = p1[1] + alpha*e1[1] + beta*e2[1];
        x[2] = p1[2] + alpha*e1[2] + beta*e2[2];
      }

      if (region_flag &&
          !region_match_kk(d_region_tokens,nregion_token,
                           x[0],x[1],x[2])) continue;

      nactual++;
      d_keep(cand) = 1;
      d_task(cand) = i;
      d_isp(cand) = isp;
      for (int d = 0; d < dimension; ++d) d_x(cand, d) = x[d];

      double beta_un, normalized_distbn_fn;
      do {
        do {
          beta_un = (6.0*rand_gen.drand() - 3.0);
        } while (beta_un + scosine < 0.0);
        normalized_distbn_fn = 2.0 * (beta_un + scosine) /
          (scosine + sqrt(scosine*scosine + 2.0)) *
          exp(0.5 + (0.5*scosine)*(scosine-sqrt(scosine*scosine + 2.0)) -
              beta_un*beta_un);
      } while (normalized_distbn_fn < rand_gen.drand());

      double vnmag;
      if (normalflag) vnmag = beta_un*vscale + magvstream;
      else vnmag = beta_un*vscale + indot;

      const double theta = MY_2PI * rand_gen.drand();
      const double vr = vscale * sqrt(-log(rand_gen.drand()));

      double vamag,vbmag;
      if (normalflag) {
        vamag = vr * sin(theta);
        vbmag = vr * cos(theta);
      } else {
        vamag = vr * sin(theta) + MathExtraKokkos::dot3(vstream,atan);
        vbmag = vr * cos(theta) + MathExtraKokkos::dot3(vstream,btan);
      }

      d_v(cand, 0) = vnmag*normal[0] + vamag*atan[0] + vbmag*btan[0];
      d_v(cand, 1) = vnmag*normal[1] + vamag*atan[1] + vbmag*btan[1];
      d_v(cand, 2) = vnmag*normal[2] + vamag*atan[2] + vbmag*btan[2];

      d_erot(cand) = particle_kk_copy.obj.erot(ispecies,temp_rot,rand_gen);
      d_evib(cand) = particle_kk_copy.obj.evib(ispecies,temp_vib,rand_gen);
      d_id(cand) = MAXSMALLINT*rand_gen.drand();
      d_dtremain(cand) = dt * rand_gen.drand();
    }

    nsingle += nactual;
  }

  rand_pool.free_state(rand_gen);
}

/* ---------------------------------------------------------------------- */

template<int ATOMIC_REDUCTION>
KOKKOS_INLINE_FUNCTION
void FixEmitSurfKokkos::operator()(TagFixEmitSurf_insert_particles<ATOMIC_REDUCTION>, const int &cand) const
{
  if (!d_keep(cand)) return;

  auto i = d_task(cand);
  Task task_i = d_tasks(i);

  const int pcell = task_i.pcell;
  const surfint isurf = task_i.isurf;

  auto isp = d_isp(cand);
  auto ispecies = d_mspecies(isp);

  double x[3];
  for (int d = 0; d < dimension; ++d) x[d] = d_x(cand, d);
  for (int d = dimension; d < 3; ++d) x[d] = 0;

  auto erot = d_erot(cand);
  auto evib = d_evib(cand);
  auto id = d_id(cand);
  auto dtremain = d_dtremain(cand);

  double v[3];
  for (int d = 0; d < 3; ++d) v[d] = d_v(cand, d);

  auto inew = d_cands2new(cand);
  auto ilocal = nlocal_before + inew;

  ParticleKokkos::add_particle_kokkos(d_particles,ilocal,
      id,ispecies,pcell,x,v,erot,evib);

  auto p = &d_particles(ilocal);

  p->flag = PSURF + 1 + isurf;
  p->dtremain = dtremain;

  if (nsurf_tally)
    for (int k = 0; k < nsurf_tally; k++)
      FES_SLIST(k).
            surf_tally_kk<ATOMIC_REDUCTION>(dtremain,isurf,pcell,0,NULL,p,NULL);
}

/* ----------------------------------------------------------------------
   recalculate task properties based on subsonic BC
------------------------------------------------------------------------- */

void FixEmitSurfKokkos::subsonic_inflow()
{
  // for grid cells that are part of tasks:
  // calculate local nrho, vstream, and thermal temperature
  // if needed sort particles for grid cells with tasks

  subsonic_sort();
  if (subsonic_style == MFLOW) mflow_grid();
  else subsonic_grid();

  // recalculate particle insertion counts for each task
  // recompute mixture vscale, since depends on temp_thermal

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,SPECIES_MASK);
  d_species_all = particle_kk->k_species.view_device();

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.view_device();

  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,ALL_MASK);
  d_lines = surf_kk->k_lines.view_device();
  d_tris = surf_kk->k_tris.view_device();
  nlocal_surf = surf->nlocal;

  k_tasks.sync_device();
  if (perspecies) k_ntargetsp.sync_device();
  k_mspecies.sync_device();
  k_fraction.sync_device();

  boltz = update->boltz;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitSurf_subsonic_inflow>(0,ntask),*this);
  copymode = 0;

  k_tasks.modify_device();
  if (perspecies) k_ntargetsp.modify_device();
}

KOKKOS_INLINE_FUNCTION
void FixEmitSurfKokkos::operator()(TagFixEmitSurf_subsonic_inflow, const int &i) const
{
  double *vstream = d_tasks(i).vstream;

  // indot = vstream dotted into inward surf normal
  // depends on normalflag, same as FixEmitSurf::subsonic_inflow()

  double indot;
  if (normalflag) indot = magvstream;
  else {
    const surfint isurf = d_tasks(i).isurf;
    double *normal = (dimension == 2) ? d_lines[isurf].norm : d_tris[isurf].norm;
    indot = vstream[0]*normal[0] + vstream[1]*normal[1];
    if (dimension != 2) indot += vstream[2]*normal[2];
  }

  const double area = d_tasks(i).area;
  const double nrho = d_tasks(i).nrho;
  const double temp_thermal = d_tasks(i).temp_thermal;
  const int icell = d_tasks(i).icell;

  double ntarget = 0.0;
  for (int isp = 0; isp < nspecies; isp++) {
    const double mass = d_species_all[d_mspecies[isp]].mass;
    const double vscale = sqrt(2.0 * boltz * temp_thermal / mass);
    double ntargetsp = mol_inflow_kokkos(indot,vscale,d_fraction[isp]);
    ntargetsp *= nrho*area*dt / fnum;
    ntargetsp /= d_cinfo[icell].weight;
    ntarget += ntargetsp;
    if (perspecies) d_ntargetsp(i,isp) = ntargetsp;
  }
  d_tasks(i).ntarget = ntarget;
  if (ntarget >= MAXSMALLINT)
    Kokkos::abort("Fix emit/surf subsonic insertion count exceeds 32-bit int");
}

/* ----------------------------------------------------------------------
   sort particles into grid cells on device
   same compressed per-cell particle lists as built for collisions,
   used in lieu of the linked lists built by FixEmitSurf::subsonic_sort()
------------------------------------------------------------------------- */

void FixEmitSurfKokkos::subsonic_sort()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;

  // sorted_kk mirrors the host Particle::sorted flag the non-Kokkos path
  //   tests here.  Record it BEFORE sorting: like the non-Kokkos
  //   subsonic_sort(), a sort done on behalf of this fix builds a list that
  //   is walked in decreasing particle index, while an already-sorted list
  //   is walked in increasing index.

  plist_descending = !particle_kk->sorted_kk;
  if (!particle_kk->sorted_kk) particle_kk->sort_kokkos();
}

/* ----------------------------------------------------------------------
   compute number density, thermal temperature, stream velocity
   only for grid cells associated with a task
   first compute for grid cells, then adjust due to boundary conditions
------------------------------------------------------------------------- */

void FixEmitSurfKokkos::subsonic_grid()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_subsonic_particles = particle_kk->k_particles.view_device();
  d_species_all = particle_kk->k_species.view_device();

  // refresh particle_kk_copy since particle data structures may
  //   have changed since the last copy, e.g. by sort or grow
  // reset surf compute copies carried in *this to a valid empty copy
  //   for the subsonic kernel dispatch (real setup happens after)

  particle_kk->update_class_variables();
  particle_kk_copy.copy(particle_kk);
#ifdef SPARTA_KOKKOS_FIXED_LISTS
  for (int n = 0; n < KOKKOS_MAX_SLIST; n++)
    slist_active_copy[n].copy(&tmp_compute_surf_kk);
#endif

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.view_device();
  d_plist = grid_kk->d_plist;
  d_cellcount = grid_kk->d_cellcount;

  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,ALL_MASK);
  d_lines = surf_kk->k_lines.view_device();
  d_tris = surf_kk->k_tris.view_device();
  nlocal_surf = surf->nlocal;

  k_tasks.sync_device();
  if (subsonic_style == PONLY) {
    k_vscale.sync_device();
    k_mspecies.sync_device();
  }

  boltz = update->boltz;
  temp_thermal_mix = particle->mixture[imix]->temp_thermal;

  // acoef = 1.0 (window = 0) -> vcom = current cell COM velocity

  acoef = 1.0;
  if (subsonic_window > 0) acoef = 1.0 / (subsonic_window + 1.0);

  // only track max thermal temp until the one-time warning has fired
  // avoids a per-step device->host fence once subsonic_warning is set

  if (!subsonic_warning) {
    if (d_tempmax.data() == nullptr)
      d_tempmax = DAT::t_float_scalar("emit/surf:tempmax");
    Kokkos::deep_copy(d_tempmax,0.0);
  }

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitSurf_subsonic_grid>(0,ntask),*this);
  copymode = 0;

  k_tasks.modify_device();
  if (subsonic_style == PONLY) k_vscale.modify_device();

  // release references to reduce memory use

  d_subsonic_particles = t_particle_1d();
  d_species_all = t_species_1d();
  d_plist = {};
  d_cellcount = {};
  d_cinfo = {};

  // test if any task has invalid thermal temperature for first time

  if (!subsonic_warning) {
    double tempmax = 0.0;
    Kokkos::deep_copy(tempmax,d_tempmax);
    int temp_exceed_flag = 0;
    if (tempmax > TEMPLIMIT) temp_exceed_flag = 1;
    subsonic_warning = subsonic_temperature_check(temp_exceed_flag,tempmax);
  }
}

KOKKOS_INLINE_FUNCTION
void FixEmitSurfKokkos::operator()(TagFixEmitSurf_subsonic_grid, const int &i) const
{
  const int icell = d_tasks(i).pcell;
  const int np = d_cellcount(icell);

  // accumulate needed per-particle quantities
  // mv = mass*velocity terms, masstot = total mass
  // gamma = rotational/tranlational DOFs

  double mv[4];
  mv[0] = mv[1] = mv[2] = mv[3] = 0.0;
  double masstot = 0.0;
  double gamma = 0.0;

  // d_plist orders particles by increasing index.  The non-Kokkos path walks
  // whichever linked list is current: the one Particle::sort() builds (head =
  // lowest index, so INcreasing order) when the particles were already
  // sorted, else the one subsonic_sort() builds itself (head = highest index,
  // so DEcreasing order).  For SPARTA_KOKKOS_EXACT match that order so the
  // per-cell moment sums are bit-identical (serial, single thread, host).

#ifdef SPARTA_KOKKOS_EXACT
  const int nbeg = plist_descending ? np-1 : 0;
  const int nend = plist_descending ? -1 : np;
  const int ninc = plist_descending ? -1 : 1;
  for (int n = nbeg; n != nend; n += ninc) {
#else
  for (int n = 0; n < np; n++) {
#endif
    const int ip = d_plist(icell,n);
    const int ispecies = d_subsonic_particles[ip].ispecies;
    const double mass = d_species_all[ispecies].mass;
    const double *v = d_subsonic_particles[ip].v;
    mv[0] += mass*v[0];
    mv[1] += mass*v[1];
    mv[2] += mass*v[2];
    mv[3] += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    masstot += mass;
    gamma += 1.0 + 2.0 / (3.0 + d_species_all[ispecies].rotdof);
  }

  // compute/store nrho, 3 temps, vstream for task
  // also vscale for PONLY
  // if sound speed = 0.0 due to <= 1 particle in cell or
  //   all particles having COM velocity, set via mixture properties

  // vcom = cell COM velocity, optionally time-averaged (window > 0)
  // acoef = 1.0 for window = 0, which reproduces the instantaneous value
  // the average is kept in vcom, not vstream, so the PONLY pressure
  //   correction added below is applied once and not integrated over steps

  double vnew[3];
  if (np && masstot > 0.0) {
    vnew[0] = mv[0] / masstot;
    vnew[1] = mv[1] / masstot;
    vnew[2] = mv[2] / masstot;
  } else vnew[0] = vnew[1] = vnew[2] = 0.0;

  double *vcom = d_tasks(i).vcom;
  vcom[0] = acoef*vnew[0] + (1.0-acoef)*vcom[0];
  vcom[1] = acoef*vnew[1] + (1.0-acoef)*vcom[1];
  vcom[2] = acoef*vnew[2] + (1.0-acoef)*vcom[2];

  double *vstream = d_tasks(i).vstream;
  vstream[0] = vcom[0];
  vstream[1] = vcom[1];
  vstream[2] = vcom[2];

  double temp_thermal_cell;

  if (subsonic_style == PTBOTH) {
    d_tasks(i).nrho = nsubsonic;
    temp_thermal_cell = tsubsonic;

  } else {
    const double nrho_cell = np * fnum / d_cinfo[icell].volume;
    const double massrho_cell = masstot * fnum / d_cinfo[icell].volume;
    if (np > 1) {
      const double ke = mv[3]/np -
        (mv[0]*mv[0] + mv[1]*mv[1] + mv[2]*mv[2])/np/masstot;
      temp_thermal_cell = tprefactor * ke;
    } else temp_thermal_cell = temp_thermal_mix;

    const double press_cell = nrho_cell * boltz * temp_thermal_cell;
    double soundspeed_cell;
    if (np) {
      const double mass_cell = masstot / np;
      const double gamma_cell = gamma / np;
      soundspeed_cell = sqrt(gamma_cell*boltz*temp_thermal_cell / mass_cell);
    } else soundspeed_cell = soundspeed_mixture;

    d_tasks(i).nrho = nrho_cell +
      (psubsonic - press_cell) / (soundspeed_cell*soundspeed_cell);
    temp_thermal_cell = psubsonic / (boltz * d_tasks(i).nrho);
    if (!subsonic_warning && temp_thermal_cell > TEMPLIMIT)
      Kokkos::atomic_max(&d_tempmax(),temp_thermal_cell);

    // adjust COM vstream by difference between
    //   cell pressure and subsonic target pressure
    // normal = direction of difference, depends on normalflag

    const double *normal;
    if (normalflag) {
      const surfint isurf = d_tasks(i).isurf;
      normal = (dimension == 2) ? d_lines[isurf].norm : d_tris[isurf].norm;
    } else normal = norm_vstream;

    if (np) {
      const double vsmag = (psubsonic - press_cell) / (massrho_cell*soundspeed_cell);
      vstream[0] += vsmag*normal[0];
      vstream[1] += vsmag*normal[1];
      vstream[2] += vsmag*normal[2];
    }

    for (int m = 0; m < nspecies; m++) {
      const int ispecies = d_mspecies[m];
      d_vscale(i,m) = sqrt(2.0 * boltz * temp_thermal_cell /
                           d_species_all[ispecies].mass);
    }
  }

  d_tasks(i).temp_thermal = temp_thermal_cell;
  d_tasks(i).temp_rot = d_tasks(i).temp_vib = temp_thermal_cell;
}

/* ----------------------------------------------------------------------
   mass-flow inlet: choose a single global number density nrho_mflow each
   step so the expected total inserted mass rate equals the target mflow
   Kokkos port of FixEmitSurf::mflow_grid(): the per-task pass (cell COM
   velocity, fixed inlet temps and vscale, S accumulation) runs on the
   device as a parallel_reduce; the cross-proc S reduction stays on the host
------------------------------------------------------------------------- */

void FixEmitSurfKokkos::mflow_grid()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_subsonic_particles = particle_kk->k_particles.view_device();
  d_species_all = particle_kk->k_species.view_device();

  // refresh particle_kk_copy since particle data structures may
  //   have changed since the last copy, e.g. by sort or grow
  // reset surf compute copies carried in *this to a valid empty copy
  //   for the kernel dispatch (real setup happens after)

  particle_kk->update_class_variables();
  particle_kk_copy.copy(particle_kk);
#ifdef SPARTA_KOKKOS_FIXED_LISTS
  for (int n = 0; n < KOKKOS_MAX_SLIST; n++)
    slist_active_copy[n].copy(&tmp_compute_surf_kk);
#endif

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.view_device();
  d_plist = grid_kk->d_plist;
  d_cellcount = grid_kk->d_cellcount;

  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,ALL_MASK);
  d_lines = surf_kk->k_lines.view_device();
  d_tris = surf_kk->k_tris.view_device();
  nlocal_surf = surf->nlocal;

  k_tasks.sync_device();
  k_vscale.sync_device();
  k_mspecies.sync_device();
  k_fraction.sync_device();

  boltz = update->boltz;

  // acoef = 1.0 (window = 0) -> vcom = current cell COM velocity

  acoef = 1.0;
  if (subsonic_window > 0) acoef = 1.0 / (subsonic_window + 1.0);

  copymode = 1;
  double S_me = 0.0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEmitSurf_mflow_grid>(0,ntask),*this,Kokkos::Sum<double>(S_me));
  copymode = 0;

  double S;
  MPI_Allreduce(&S_me,&S,1,MPI_DOUBLE,MPI_SUM,world);

  // global number density that makes expected inserted mass rate = mflow
  // S == 0 (no inflow on any proc) -> insert nothing this step, self-corrects

  nrho_mflow = 0.0;
  if (S > 0.0) nrho_mflow = mflow / S;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitSurf_mflow_nrho>(0,ntask),*this);
  copymode = 0;

  k_tasks.modify_device();
  k_vscale.modify_device();

  // release references to reduce memory use

  d_subsonic_particles = t_particle_1d();
  d_species_all = t_species_1d();
  d_plist = {};
  d_cellcount = {};
  d_cinfo = {};
}

KOKKOS_INLINE_FUNCTION
void FixEmitSurfKokkos::operator()(TagFixEmitSurf_mflow_grid, const int &i, double &S_me) const
{
  const int icell = d_tasks(i).pcell;
  const int np = d_cellcount(icell);

  // accumulate needed per-particle quantities
  // mv = mass*velocity terms, masstot = total mass

  double mv[3];
  mv[0] = mv[1] = mv[2] = 0.0;
  double masstot = 0.0;

  // d_plist orders particles by increasing index.  The non-Kokkos path walks
  // whichever linked list is current: the one Particle::sort() builds (head =
  // lowest index, so INcreasing order) when the particles were already
  // sorted, else the one subsonic_sort() builds itself (head = highest index,
  // so DEcreasing order).  For SPARTA_KOKKOS_EXACT match that order so the
  // per-cell moment sums are bit-identical (serial, single thread, host).

#ifdef SPARTA_KOKKOS_EXACT
  const int nbeg = plist_descending ? np-1 : 0;
  const int nend = plist_descending ? -1 : np;
  const int ninc = plist_descending ? -1 : 1;
  for (int n = nbeg; n != nend; n += ninc) {
#else
  for (int n = 0; n < np; n++) {
#endif
    const int ip = d_plist(icell,n);
    const int ispecies = d_subsonic_particles[ip].ispecies;
    const double mass = d_species_all[ispecies].mass;
    const double *v = d_subsonic_particles[ip].v;
    mv[0] += mass*v[0];
    mv[1] += mass*v[1];
    mv[2] += mass*v[2];
    masstot += mass;
  }

  // vcom = cell COM velocity, optionally time-averaged (window > 0)
  // acoef = 1.0 for window = 0, which reproduces the instantaneous value

  double vnew[3];
  if (np && masstot > 0.0) {
    vnew[0] = mv[0] / masstot;
    vnew[1] = mv[1] / masstot;
    vnew[2] = mv[2] / masstot;
  } else vnew[0] = vnew[1] = vnew[2] = 0.0;

  double *vcom = d_tasks(i).vcom;
  vcom[0] = acoef*vnew[0] + (1.0-acoef)*vcom[0];
  vcom[1] = acoef*vnew[1] + (1.0-acoef)*vcom[1];
  vcom[2] = acoef*vnew[2] + (1.0-acoef)*vcom[2];

  double *vstream = d_tasks(i).vstream;
  vstream[0] = vcom[0];
  vstream[1] = vcom[1];
  vstream[2] = vcom[2];

  // fixed inlet temperatures and per-species vscale at tmflow
  // vscale depends only on the species, so the values are identical to the
  //   once-per-call precompute the non-Kokkos mflow_grid() performs

  d_tasks(i).temp_thermal = tmflow;
  d_tasks(i).temp_rot = d_tasks(i).temp_vib = tmflow;

  for (int isp = 0; isp < nspecies; isp++)
    d_vscale(i,isp) = sqrt(2.0 * boltz * tmflow /
                           d_species_all[d_mspecies[isp]].mass);

  // indot identical to the subsonic_inflow tail and perform_task

  double indot;
  if (normalflag) indot = magvstream;
  else {
    const surfint isurf = d_tasks(i).isurf;
    double *normal = (dimension == 2) ? d_lines[isurf].norm : d_tris[isurf].norm;
    indot = vstream[0]*normal[0] + vstream[1]*normal[1];
    if (dimension != 2) indot += vstream[2]*normal[2];
  }

  // mass weighting depends on how perform_task picks the species

  double fluxmass = 0.0;
  for (int isp = 0; isp < nspecies; isp++) {
    const double flux = mol_inflow_kokkos(indot,d_vscale(i,isp),d_fraction[isp]);
    if (perspecies) fluxmass += flux * d_species_all[d_mspecies[isp]].mass;
    else fluxmass += flux * avemass_mixture;
  }

  S_me += fluxmass * d_tasks(i).area;
}

KOKKOS_INLINE_FUNCTION
void FixEmitSurfKokkos::operator()(TagFixEmitSurf_mflow_nrho, const int &i) const
{
  d_tasks(i).nrho = nrho_mflow;
}

/* ----------------------------------------------------------------------
   grow task list
------------------------------------------------------------------------- */

void FixEmitSurfKokkos::grow_task()
{
  int oldmax = ntaskmax;
  ntaskmax += DELTATASK;

  k_tasks.sync_host();
  k_tasks.modify_host(); // force resize on host
  memoryKK->grow_kokkos(k_tasks,tasks,ntaskmax,"emit/surf:tasks");
  d_tasks = k_tasks.view_device();

  // allocate vectors in each new task or set to NULL

  if (perspecies) {
    k_ntargetsp.sync_host();
    k_ntargetsp.modify_host(); // force resize on host
    k_ntargetsp.resize(ntaskmax,nspecies);
    d_ntargetsp = k_ntargetsp.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].ntargetsp = &k_ntargetsp.view_host()(i,0);
  }

  if (subsonic_style == PONLY || subsonic_style == MFLOW || temp_custom_flag) {
    k_vscale.sync_host();
    k_vscale.modify_host(); // force resize on host
    k_vscale.resize(ntaskmax,nspecies);
    d_vscale = k_vscale.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].vscale = &k_vscale.view_host()(i,0);
  }

  for (int i = oldmax; i < ntaskmax; i++) {
    tasks[i].path = NULL;
    tasks[i].fracarea = NULL;
  }
}

/* ----------------------------------------------------------------------
   reallocate nspecies arrays
------------------------------------------------------------------------- */

void FixEmitSurfKokkos::realloc_nspecies()
{
  if (perspecies) {
    k_ntargetsp = DAT::tdual_float_2d_lr("emit/surf:ntargetsp",ntaskmax,nspecies);
    d_ntargetsp = k_ntargetsp.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].ntargetsp = &k_ntargetsp.view_host()(i,0);
  }
  if (subsonic_style == PONLY || subsonic_style == MFLOW || temp_custom_flag) {
    k_vscale = DAT::tdual_float_2d_lr("emit/surf:vscale",ntaskmax,nspecies);
    d_vscale = k_vscale.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].vscale = &k_vscale.view_host()(i,0);
  }
}
