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
#include "sparta_masks.h"
#include "variable.h"
#include "Kokkos_Random.hpp"

using namespace SPARTA_NS;
using namespace MathConst;

enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{NOSUBSONIC,PTBOTH,PONLY};
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
  particle_kk_copy(sparta),
  slist_active_copy{VAL_2(KKCopy<ComputeSurfKokkos>(sparta))},
  tmp_compute_surf_kk(sparta)
{
  kokkos_flag = 1;
  execution_space = Device;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

FixEmitSurfKokkos::~FixEmitSurfKokkos()
{
  if (copymode) return;

  particle_kk_copy.uncopy();

  for (int i=0; i<KOKKOS_MAX_SLIST; i++) {
    slist_active_copy[i].uncopy();
  }

  tmp_compute_surf_kk.uncopy = 1;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
#endif

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
  if (subsonic_style == PONLY) k_vscale.sync_host();

  k_path.sync_host();

  if (dimension != 2)
    k_fracarea.sync_host();

  FixEmitSurf::init();

  k_tasks.modify_host();
  if (perspecies) k_ntargetsp.modify_host();
  if (subsonic_style == PONLY) k_vscale.modify_host();

  k_path.modify_host();

  if (dimension != 2)
    k_fracarea.modify_host();

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif

  k_vscale_mix = DAT::tdual_float_1d("vscale_mix", nspecies);

  if (!(fractions_custom_flag && !perspecies))
    k_cummulative_mix = DAT::tdual_float_1d("cummulative", nspecies);

  k_species = DAT::tdual_int_1d("species", nspecies);

  d_vscale_mix = k_vscale_mix .view_device();
  d_cummulative_mix = k_cummulative_mix.view_device();
  d_species = k_species.view_device();

  auto h_vscale_mix  = k_vscale_mix.view_host();
  auto h_cummulative_mix = k_cummulative_mix.view_host();
  auto h_species = k_species.view_host();

  for (int isp = 0; isp < nspecies; ++isp) {
    h_vscale_mix(isp) = particle->mixture[imix]->vscale[isp];

    if (!(fractions_custom_flag && !perspecies))
      h_cummulative_mix(isp) = particle->mixture[imix]->cummulative[isp];

    h_species(isp) = particle->mixture[imix]->species[isp];
  }

  k_vscale_mix.modify_host();

  if (!(fractions_custom_flag && !perspecies))
    k_cummulative_mix.modify_host();

  k_species.modify_host();
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
  FixEmit::create_tasks();

  k_tasks.modify_host();
  k_ntargetsp.modify_host();
  k_vscale.modify_host();

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

/* ---------------------------------------------------------------------- */

void FixEmitSurfKokkos::perform_task()
{
  dt = update->dt;

  // if subsonic, re-compute particle inflow counts for each task
  // also computes current per-task temp_thermal and vstream

  if (subsonic)
    error->one(FLERR,"Cannot (yet) use subsonic emission in fix emit/surf/kk");

  //if (subsonic) subsonic_inflow(); ////////////////////////

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

  if (perspecies) k_ntargetsp.sync_device();
  else k_tasks.sync_device();

  SurfKokkos* surf_kk = (SurfKokkos*) surf;
  surf_kk->sync(Device,ALL_MASK);
  d_lines = surf_kk->k_lines.view_device();
  d_tris = surf_kk->k_tris.view_device();
  nlocal_surf = surf->nlocal;

  nsurf_tally = update->nsurf_tally;
  Compute **slist_active = update->slist_active;

  if (nsurf_tally) {
    for (int i = 0; i < nsurf_tally; i++) {
      if (strcmp(slist_active[i]->style,"isurf/grid") == 0)
        error->all(FLERR,"Kokkos doesn't yet support compute isurf/grid");
      ComputeSurfKokkos* compute_surf_kk = dynamic_cast<ComputeSurfKokkos*>(slist_active[i]);
      if (!compute_surf_kk)
        error->all(FLERR,"Kokkos does not (yet) support compute surf/collision/tally or compute surf/reaction/tally");
      compute_surf_kk->pre_surf_tally();
      slist_active_copy[i].copy(compute_surf_kk);
    }
  } else {
    for (int i = 0; i < KOKKOS_MAX_SLIST; i++) {

      // use temporary to avoid the copy getting stale leading to an issue
      //  with view reference counting

      slist_active_copy[i].copy(&tmp_compute_surf_kk);
    }
  }

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
  if (subsonic_style == PONLY) k_vscale.sync_device();

  k_path.sync_device();

  if (dimension != 2)
    k_fracarea.sync_device();

  // copy needed mixture data to device

  k_vscale_mix .sync_device();
  k_species    .sync_device();

  if (fractions_custom_flag && !perspecies)
    k_cummulative_custom.sync_device();
  else
    k_cummulative_mix.sync_device();

  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);
  particle_kk->update_class_variables();
  particle_kk_copy.copy(particle_kk);

  if (region)
    error->one(FLERR,"Cannot yet use fix emit/surf/kk with regions");

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

      double vscale = (subsonic_style == PONLY) ?
        d_vscale(i, isp) : d_vscale_mix(isp);

      const int ispecies = d_species[isp];
      const double ninsert = d_ninsert(i * nspecies + isp);
      const int start = d_task2cand(i * nspecies + isp);
      const double scosine = indot / vscale;

      // loop over ninsert for each species

      int nactual = 0;
      for (int m = 0; m < ninsert; m++) {
        auto cand = start + m;

        double x[3];
        if (dimension == 2) {
          const double rn = rand_gen.drand();
          double* p1 = &d_path(i,0);
          double* p2 = &d_path(i,3);
          x[0] = p1[0] + rn * (p2[0]-p1[0]);
          x[1] = p1[1] + rn * (p2[1]-p1[1]);
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

        //if (region && !region->match(x)) continue; ////////////////////////
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

      double vscale = (subsonic_style == PONLY) ?
        d_vscale(i, isp) : d_vscale_mix(isp);

      const int ispecies = d_species[isp];
      const double scosine = indot / vscale;

      double x[3];
      if (dimension == 2) {
        const double rn = rand_gen.drand();
        double* p1 = &d_path(i,0);
        double* p2 = &d_path(i,3);
        x[0] = p1[0] + rn * (p2[0]-p1[0]);
        x[1] = p1[1] + rn * (p2[1]-p1[1]);
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

      //if (region && !region->match(x)) continue; ////////////////////////
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
  auto ispecies = d_species(isp);

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
      slist_active_copy[k].obj.
            surf_tally_kk<ATOMIC_REDUCTION>(dtremain,isurf,pcell,0,NULL,p,NULL);
}

/* ----------------------------------------------------------------------
   grow task list
------------------------------------------------------------------------- */

void FixEmitSurfKokkos::grow_task()
{
  int oldmax = ntaskmax;
  ntaskmax += DELTATASK;

  if (tasks == NULL)
    k_tasks = tdual_task_1d("emit/surf:tasks",ntaskmax);
  else {
    k_tasks.sync_host();
    k_tasks.modify_host(); // force resize on host
    k_tasks.resize(ntaskmax);
  }
  d_tasks = k_tasks.view_device();
  tasks = k_tasks.view_host().data();

  // set all new task bytes to 0 so valgrind won't complain
  // if bytes between fields are uninitialized

  //memset(&tasks[oldmax],0,(ntaskmax-oldmax)*sizeof(Task));

  // allocate vectors in each new task or set to NULL

  if (perspecies) {
    k_ntargetsp.sync_host();
    k_ntargetsp.modify_host(); // force resize on host
    k_ntargetsp.resize(ntaskmax,nspecies);
    d_ntargetsp = k_ntargetsp.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].ntargetsp = k_ntargetsp.view_host().data() + i*k_ntargetsp.view_host().extent(1);
  }

  if (subsonic_style == PONLY) {
    k_vscale.modify_host(); // force resize on host
    k_vscale.sync_host();
    k_vscale.resize(ntaskmax,nspecies);
    d_vscale = k_vscale.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].vscale = k_vscale.view_host().data() + i*k_vscale.view_host().extent(1);
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
      tasks[i].ntargetsp = k_ntargetsp.view_host().data() + i*k_ntargetsp.view_host().extent(1);
  }
  if (subsonic_style == PONLY) {
    k_vscale = DAT::tdual_float_2d_lr("emit/surf:vscale",ntaskmax,nspecies);
    d_vscale = k_vscale.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].vscale = k_vscale.view_host().data() + i*k_vscale.view_host().extent(1);
  }
}
