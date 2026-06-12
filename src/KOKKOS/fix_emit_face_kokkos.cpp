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
#include "fix_emit_face_kokkos.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "grid.h"
#include "surf.h"
#include "particle.h"
#include "mixture.h"
#include "comm.h"
#include "modify.h"
#include "geometry.h"
#include "input.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos_type.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "sparta_masks.h"
#include "variable.h"
#include "Kokkos_Random.hpp"

using namespace SPARTA_NS;
using namespace MathConst;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // same as Grid
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Grid
enum{NOSUBSONIC,PTBOTH,PONLY};

#define DELTATASK 256
#define TEMPLIMIT 1.0e5

#define VAL_1(X) X
#define VAL_2(X) VAL_1(X), VAL_1(X)

/* ----------------------------------------------------------------------
   insert particles in grid cells with faces touching inflow boundaries
------------------------------------------------------------------------- */

FixEmitFaceKokkos::FixEmitFaceKokkos(SPARTA *sparta, int narg, char **arg) :
  FixEmitFace(sparta, narg, arg),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            ),
  particle_kk_copy(sparta),
  regblock_kk_copy(sparta),
  regcylinder_kk_copy(sparta),
  regplane_kk_copy(sparta),
  regsphere_kk_copy(sparta)
{
  kokkos_flag = 1;
  execution_space = Device;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  region_flag = 0;
}

/* ---------------------------------------------------------------------- */

FixEmitFaceKokkos::~FixEmitFaceKokkos()
{
  if (copymode) return;

  particle_kk_copy.uncopy();
  regblock_kk_copy.uncopy();
  regcylinder_kk_copy.uncopy();
  regplane_kk_copy.uncopy();
  regsphere_kk_copy.uncopy();

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
#endif

  for (int i = 0; i < ntaskmax; i++) {
    tasks[i].ntargetsp = NULL;
    tasks[i].vscale = NULL;
  }

  tasks = NULL;
}

/* ---------------------------------------------------------------------- */

void FixEmitFaceKokkos::init()
{
  k_tasks.sync_host();
  if (perspecies) k_ntargetsp.sync_host();
  if (subsonic_style == PONLY) k_vscale.sync_host();

  FixEmitFace::init();

  k_tasks.modify_host();
  if (perspecies) k_ntargetsp.modify_host();
  if (subsonic_style == PONLY) k_vscale.modify_host();

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif

  k_mix_vscale  = DAT::tdual_float_1d("mix_vscale", nspecies);
  k_cummulative = DAT::tdual_float_1d("cummulative", nspecies);
  k_mspecies    = DAT::tdual_int_1d("mspecies", nspecies);
  k_fraction    = DAT::tdual_float_1d("fraction", nspecies);

  d_mix_vscale  = k_mix_vscale .view_device();
  d_cummulative = k_cummulative.view_device();
  d_mspecies    = k_mspecies   .view_device();
  d_fraction    = k_fraction   .view_device();

  auto h_mix_vscale  = k_mix_vscale .view_host();
  auto h_cummulative = k_cummulative.view_host();
  auto h_mspecies    = k_mspecies   .view_host();
  auto h_fraction    = k_fraction   .view_host();

  for (int isp = 0; isp < nspecies; ++isp) {
    h_mix_vscale(isp) = particle->mixture[imix]->vscale[isp];
    h_cummulative(isp) = particle->mixture[imix]->cummulative[isp];
    h_mspecies(isp) = particle->mixture[imix]->species[isp];
    h_fraction(isp) = particle->mixture[imix]->fraction[isp];
  }

  k_mix_vscale .modify_host();
  k_cummulative.modify_host();
  k_mspecies   .modify_host();
  k_fraction   .modify_host();
}

/* ----------------------------------------------------------------------
   create tasks for one grid cell
   add them to tasks list and increment ntasks
------------------------------------------------------------------------- */

void FixEmitFaceKokkos::create_tasks()
{
  k_tasks.sync_host();
  if (perspecies) k_ntargetsp.sync_host();
  if (subsonic_style == PONLY) k_vscale.sync_host();

  FixEmitFace::create_tasks();

  k_tasks.modify_host();
  if (perspecies) k_ntargetsp.modify_host();
  if (subsonic_style == PONLY) k_vscale.modify_host();
}

/* ---------------------------------------------------------------------- */

void FixEmitFaceKokkos::perform_task()
{
  dt = update->dt;
  auto l_dimension = this->dimension;
  auto l_subsonic_style = this->subsonic_style;

  // if subsonic, re-compute particle inflow counts for each task
  // also computes current per-task temp_thermal and vstream

  if (subsonic) subsonic_inflow();

  // if modulate variable set, evaluate it as prefactor for this timestep

  prefactor = 1.0;
  if (modvar) {
    prefactor = input->variable->compute_equal(imodvar);
    if (prefactor < 0.0) error->all(FLERR,"Fix emit/face modulation < 0.0");
  }

  // insert particles for each task = cell/face pair
  // ntarget/ninsert is either perspecies or for all species

  // copy needed task data to device

  if (perspecies) k_ntargetsp.sync_device();
  else k_tasks.sync_device();

  auto ninsert_dim1 = perspecies ? nspecies : 1;
  if (d_ninsert.extent(0) < ntask * ninsert_dim1)
    d_ninsert = DAT::t_int_1d("ninsert", ntask * ninsert_dim1);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitFace_ninsert>(0,ntask),*this);
  copymode = 0;

  int ncands;
  d_task2cand = offset_scan(d_ninsert, ncands);

  if (ncands == 0) return;

  // for one particle:
  //   x = random position on face
  //   v = randomized thermal velocity + vstream
  //       first stage: normal dimension (ndim)
  //       second stage: parallel dimensions (pdim,qdim)

  // double while loop until randomized particle velocity meets 2 criteria
  // inner do-while loop:
  //   v = vstream-component + vthermal is into simulation box
  //   see Bird 1994, p 425
  // outer do-while loop:
  //   shift Maxwellian distribution by stream velocity component
  //   see Bird 1994, p 259, eq 12.5

  if (d_x.extent(0) < ncands || d_x.extent(1) < l_dimension)
    d_x = DAT::t_float_2d("x", ncands, l_dimension);

  if (d_task.extent(0) < ncands) {
    d_beta_un  = DAT::t_float_1d("beta_un", ncands);
    d_theta    = DAT::t_float_1d("theta", ncands);
    d_vr       = DAT::t_float_1d("vr", ncands);
    d_erot     = DAT::t_float_1d("erot", ncands);
    d_evib     = DAT::t_float_1d("evib", ncands);
    d_dtremain = DAT::t_float_1d("dtremain", ncands);
    d_id       = DAT::t_int_1d("id", ncands);
    d_isp      = DAT::t_int_1d("isp", ncands);
    d_task     = DAT::t_int_1d("task", ncands);
    d_keep     = DAT::t_int_1d("keep", ncands);
  }
  Kokkos::deep_copy(d_keep,0); // needs to be initialized with zeros

  auto ld_x        = d_x       ;
  auto ld_beta_un  = d_beta_un ;
  auto ld_theta    = d_theta   ;
  auto ld_vr       = d_vr      ;
  auto ld_erot     = d_erot    ;
  auto ld_evib     = d_evib    ;
  auto ld_dtremain = d_dtremain;
  auto ld_id       = d_id      ;
  auto ld_isp      = d_isp     ;
  auto ld_task     = d_task    ;
  auto ld_keep     = d_keep    ;

  // copy needed task data to device

  k_tasks.sync_device();
  if (perspecies) k_ntargetsp.sync_device();
  if (subsonic_style == PONLY) k_vscale.sync_device();
  auto ld_tasks = d_tasks;
  auto ld_vscale = d_vscale;

  // copy needed mixture data to device

  k_mix_vscale .sync_device();
  k_mspecies   .sync_device();
  k_cummulative.sync_device();

  auto ld_mix_vscale = d_mix_vscale;
  auto ld_mspecies   = d_mspecies  ;

  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);
  particle_kk->update_class_variables();
  particle_kk_copy.copy(particle_kk);

  if (region && !region->kokkos_flag)
    error->all(FLERR,"KOKKOS package does not (yet) support chosen region style");

  region_flag = 0;
  if (region) {
    if (strstr(region->style,"block") != NULL) {
      RegBlockKokkos* region_kk = ((RegBlockKokkos*)region);
      regblock_kk_copy.copy(region_kk);
      region_flag = 1;
    } else if (strstr(region->style,"cylinder") != NULL) {
      RegCylinderKokkos* region_kk = ((RegCylinderKokkos*)region);
      regcylinder_kk_copy.copy(region_kk);
      region_flag = 2;
    } else if (strstr(region->style,"plane") != NULL) {
      RegPlaneKokkos* region_kk = ((RegPlaneKokkos*)region);
      regplane_kk_copy.copy(region_kk);
      region_flag = 3;
    } else if (strstr(region->style,"sphere") != NULL) {
      RegSphereKokkos* region_kk = ((RegSphereKokkos*)region);
      regsphere_kk_copy.copy(region_kk);
      region_flag = 4;
    } else {
      error->all(FLERR,"KOKKOS package does not (yet) support chosen region style");
    }
  }

  int nsingle_reduce = 0;
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEmitFace_perform_task>(0,ntask),*this,nsingle_reduce);
  copymode = 0;
  nsingle += nsingle_reduce;

  int nnew;
  auto ld_cands2new = offset_scan(d_keep, nnew);

  auto particleKK = dynamic_cast<ParticleKokkos*>(particle);
  auto nlocal_before = particleKK->nlocal;
  particleKK->grow(nnew);
  particleKK->sync(SPARTA_NS::Device, PARTICLE_MASK);
  auto ld_particles = particleKK->k_particles.view_device();

  Kokkos::parallel_for(ncands, SPARTA_LAMBDA(int cand) {
    if (!ld_keep(cand)) return;

    double *normal,*vstream;

    auto i = ld_task(cand);
    Task task_i = ld_tasks(i);

    const int pcell = task_i.pcell;
    const int ndim = task_i.ndim;
    const int pdim = task_i.pdim;
    const int qdim = task_i.qdim;
    normal = task_i.normal;
    vstream = task_i.vstream;

    auto isp = ld_isp(cand);

    auto vscale_val = (l_subsonic_style == PONLY) ?
      ld_vscale(i, isp) : ld_mix_vscale(isp);

    auto ispecies = ld_mspecies(isp);

    double x[3];
    for (int d = 0; d < l_dimension; ++d) x[d] = ld_x(cand, d);
    for (int d = l_dimension; d < 3; ++d) x[d] = 0;

    auto beta_un = ld_beta_un(cand);

    auto theta = ld_theta(cand);
    auto vr = ld_vr(cand);
    auto erot = ld_erot(cand);
    auto evib = ld_evib(cand);
    auto id = ld_id(cand);
    auto dtremain = ld_dtremain(cand);

    double v[3];
    v[ndim] = beta_un*vscale_val*normal[ndim] + vstream[ndim];
    v[pdim] = vr * sin(theta) + vstream[pdim];
    v[qdim] = vr * cos(theta) + vstream[qdim];

    auto inew = ld_cands2new(cand);
    auto ilocal = nlocal_before + inew;

    ParticleKokkos::add_particle_kokkos(ld_particles,ilocal,
        id,ispecies,pcell,x,v,erot,evib);

    ld_particles(ilocal).flag = PINSERT;
    ld_particles(ilocal).dtremain = dtremain;
  });
  particleKK->nlocal = nlocal_before + nnew;
  particleKK->modify(SPARTA_NS::Device, PARTICLE_MASK);

  if (modify->n_update_custom) {
    auto h_keep = Kokkos::create_mirror_view(d_keep);
    auto h_task = Kokkos::create_mirror_view(d_task);
    Kokkos::deep_copy(h_keep, d_keep);
    Kokkos::deep_copy(h_task, d_task);

    // copy needed task data to host

    k_tasks.sync_host();

    auto h_cands2new = Kokkos::create_mirror_view(ld_cands2new);
    Kokkos::deep_copy(h_cands2new, ld_cands2new);
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

KOKKOS_INLINE_FUNCTION
void FixEmitFaceKokkos::operator()(TagFixEmitFace_ninsert, const int &i) const
{
  int ninsert;

  rand_type rand_gen = rand_pool.get_state();

  if (perspecies) {
    for (int isp = 0; isp < nspecies; isp++) {
      auto ntarget = prefactor*d_ntargetsp(i,isp) + rand_gen.drand();
      ninsert = static_cast<int> (ntarget);
      d_ninsert(i * nspecies + isp) = ninsert;
    }
  } else {
    if (np == 0) {
      auto ntarget = prefactor*d_tasks(i).ntarget + rand_gen.drand();
      ninsert = static_cast<int> (ntarget);
    } else {
      ninsert = npertask;
      if (i >= nthresh) ninsert++;
    }
    d_ninsert(i) = ninsert;
  }

  rand_pool.free_state(rand_gen);
}

KOKKOS_INLINE_FUNCTION
void FixEmitFaceKokkos::operator()(TagFixEmitFace_perform_task, const int &i, int &nsingle) const
{
  double *lo,*hi,*normal,*vstream;

  rand_type rand_gen = rand_pool.get_state();

  Task task_i = d_tasks(i);

  lo = task_i.lo;
  hi = task_i.hi;
  normal = task_i.normal;

  const double temp_rot = task_i.temp_rot;
  const double temp_vib = task_i.temp_vib;
  vstream = task_i.vstream;

  auto indot = vstream[0]*normal[0] + vstream[1]*normal[1] + vstream[2]*normal[2];

  if (perspecies) {
    for (int isp = 0; isp < nspecies; isp++) {
      auto vscale_val = (subsonic_style == PONLY) ?
        d_vscale(i, isp) : d_mix_vscale(isp);

      auto ispecies = d_mspecies[isp];
      auto ninsert = d_ninsert(i * nspecies + isp);
      auto start = d_task2cand(i * nspecies + isp);
      auto scosine = indot / vscale_val;

      int nactual = 0;
      for (int m = 0; m < ninsert; m++) {
        auto cand = start + m;
        double x[3];
        x[0] = lo[0] + rand_gen.drand() * (hi[0]-lo[0]);
        x[1] = lo[1] + rand_gen.drand() * (hi[1]-lo[1]);
        if (dimension == 3) x[2] = lo[2] + rand_gen.drand() * (hi[2]-lo[2]);
        else x[2] = 0.0;

        if (region_flag == 1) {
          if (!regblock_kk_copy.obj.match_kokkos(x[0], x[1], x[2])) continue;
        } else if (region_flag == 2) {
          if (!regcylinder_kk_copy.obj.match_kokkos(x[0], x[1], x[2])) continue;
        } else if (region_flag == 3) {
          if (!regplane_kk_copy.obj.match_kokkos(x[0], x[1], x[2])) continue;
        } else if (region_flag == 4) {
          if (!regsphere_kk_copy.obj.match_kokkos(x[0], x[1], x[2])) continue;
        }

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

        d_beta_un(cand) = beta_un;

        d_theta(cand) = MY_2PI * rand_gen.drand();
        d_vr(cand) = vscale_val * sqrt(-log(rand_gen.drand()));
        d_erot(cand) = particle_kk_copy.obj.erot(ispecies,temp_rot,rand_gen);
        d_evib(cand) = particle_kk_copy.obj.evib(ispecies,temp_vib,rand_gen);
        d_id(cand) = MAXSMALLINT*rand_gen.drand();
        d_dtremain(cand) = dt * rand_gen.drand();
      }
      nsingle += nactual;
    }
  } else {
    auto ninsert = d_ninsert(i);
    auto start = d_task2cand(i);

    int nactual = 0;
    for (int m = 0; m < ninsert; m++) {
      auto cand = start + m;
      auto rn = rand_gen.drand();
      int isp = 0;
      while (d_cummulative[isp] < rn) isp++;
      auto vscale_val = (subsonic_style == PONLY) ?
        d_vscale(i, isp) : d_mix_vscale(isp);
      auto ispecies = d_mspecies[isp];
      auto scosine = indot / vscale_val;

      double x[3];
      x[0] = lo[0] + rand_gen.drand() * (hi[0]-lo[0]);
      x[1] = lo[1] + rand_gen.drand() * (hi[1]-lo[1]);
      if (dimension == 3) x[2] = lo[2] + rand_gen.drand() * (hi[2]-lo[2]);
      else x[2] = 0.0;

      if (region_flag == 1) {
        if (!regblock_kk_copy.obj.match_kokkos(x[0], x[1], x[2])) continue;
      } else if (region_flag == 2) {
        if (!regcylinder_kk_copy.obj.match_kokkos(x[0], x[1], x[2])) continue;
      } else if (region_flag == 3) {
        if (!regplane_kk_copy.obj.match_kokkos(x[0], x[1], x[2])) continue;
      } else if (region_flag == 4) {
        if (!regsphere_kk_copy.obj.match_kokkos(x[0], x[1], x[2])) continue;
      }

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

      d_beta_un(cand) = beta_un;

      d_theta(cand) = MY_2PI * rand_gen.drand();
      d_vr(cand) = vscale_val * sqrt(-log(rand_gen.drand()));
      d_erot(cand) = particle_kk_copy.obj.erot(ispecies,temp_rot,rand_gen);
      d_evib(cand) = particle_kk_copy.obj.evib(ispecies,temp_vib,rand_gen);
      d_id(cand) = MAXSMALLINT*rand_gen.drand();
      d_dtremain(cand) = dt * rand_gen.drand();
    }

    nsingle += nactual;
  }

  rand_pool.free_state(rand_gen);
}

/* ----------------------------------------------------------------------
   recalculate task properties based on subsonic BC
------------------------------------------------------------------------- */

void FixEmitFaceKokkos::subsonic_inflow()
{
  // for grid cells that are part of tasks:
  // calculate local nrho, vstream, and thermal temperature
  // if needed sort particles for grid cells with tasks

  subsonic_sort();
  subsonic_grid();

  // recalculate particle insertion counts for each task
  // recompute mixture vscale, since depends on temp_thermal

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,SPECIES_MASK);
  d_species = particle_kk->k_species.view_device();

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.view_device();

  k_tasks.sync_device();
  if (perspecies) k_ntargetsp.sync_device();
  k_mspecies.sync_device();
  k_fraction.sync_device();

  boltz = update->boltz;

  int errflag = 0;
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEmitFace_subsonic_inflow>(0,ntask),*this,errflag);
  copymode = 0;

  k_tasks.modify_device();
  if (perspecies) k_ntargetsp.modify_device();

  if (errflag)
    error->one(FLERR,
               "Fix emit/face subsonic insertion count exceeds 32-bit int");
}

KOKKOS_INLINE_FUNCTION
void FixEmitFaceKokkos::operator()(TagFixEmitFace_subsonic_inflow, const int &i, int &errflag) const
{
  double *vstream = d_tasks(i).vstream;
  double *normal = d_tasks(i).normal;
  const double indot = vstream[0]*normal[0] + vstream[1]*normal[1] +
    vstream[2]*normal[2];

  const double area = d_tasks(i).area;
  const double nrho = d_tasks(i).nrho;
  const double temp_thermal = d_tasks(i).temp_thermal;
  const int icell = d_tasks(i).icell;

  double ntarget = 0.0;
  for (int isp = 0; isp < nspecies; isp++) {
    const double mass = d_species[d_mspecies[isp]].mass;
    const double vscale = sqrt(2.0 * boltz * temp_thermal / mass);
    double ntargetsp = mol_inflow_kokkos(indot,vscale,d_fraction[isp]);
    ntargetsp *= nrho*area*dt / fnum;
    ntargetsp /= d_cinfo[icell].weight;
    ntarget += ntargetsp;
    if (perspecies) d_ntargetsp(i,isp) = ntargetsp;
  }
  d_tasks(i).ntarget = ntarget;
  if (ntarget >= MAXSMALLINT) errflag++;
}

/* ----------------------------------------------------------------------
   sort particles into grid cells on device
   same compressed per-cell particle lists as built for collisions,
   used in lieu of the linked lists built by FixEmitFace::subsonic_sort()
------------------------------------------------------------------------- */

void FixEmitFaceKokkos::subsonic_sort()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  if (!particle_kk->sorted_kk) particle_kk->sort_kokkos();
}

/* ----------------------------------------------------------------------
   compute number density, thermal temperature, stream velocity
   only for grid cells associated with a task
   first compute for grid cells, then adjust due to boundary conditions
------------------------------------------------------------------------- */

void FixEmitFaceKokkos::subsonic_grid()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species = particle_kk->k_species.view_device();

  // refresh particle_kk_copy since particle data structures may
  //   have changed since the last copy, e.g. by sort or grow

  particle_kk->update_class_variables();
  particle_kk_copy.copy(particle_kk);

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.view_device();
  d_plist = grid_kk->d_plist;
  d_cellcount = grid_kk->d_cellcount;

  k_tasks.sync_device();
  if (subsonic_style == PONLY) {
    k_vscale.sync_device();
    k_mspecies.sync_device();
  }

  boltz = update->boltz;
  temp_thermal_mix = particle->mixture[imix]->temp_thermal;

  if (d_tempmax.data() == nullptr)
    d_tempmax = DAT::t_float_scalar("emit/face:tempmax");
  Kokkos::deep_copy(d_tempmax,0.0);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitFace_subsonic_grid>(0,ntask),*this);
  copymode = 0;

  k_tasks.modify_device();
  if (subsonic_style == PONLY) k_vscale.modify_device();

  d_particles = t_particle_1d(); // destroy reference to reduce memory use
  d_plist = {};

  // test if any task has invalid thermal temperature for first time

  double tempmax = 0.0;
  Kokkos::deep_copy(tempmax,d_tempmax);
  int temp_exceed_flag = 0;
  if (tempmax > TEMPLIMIT) temp_exceed_flag = 1;

  if (!subsonic_warning)
    subsonic_warning = subsonic_temperature_check(temp_exceed_flag,tempmax);
}

KOKKOS_INLINE_FUNCTION
void FixEmitFaceKokkos::operator()(TagFixEmitFace_subsonic_grid, const int &i) const
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

  for (int n = 0; n < np; n++) {
    const int ip = d_plist(icell,n);
    const int ispecies = d_particles[ip].ispecies;
    const double mass = d_species[ispecies].mass;
    const double *v = d_particles[ip].v;
    mv[0] += mass*v[0];
    mv[1] += mass*v[1];
    mv[2] += mass*v[2];
    mv[3] += mass * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    masstot += mass;
    gamma += 1.0 + 2.0 / (3.0 + d_species[ispecies].rotdof);
  }

  // compute/store nrho, 3 temps, vstream for task
  // also vscale for PONLY
  // if sound speed = 0.0 due to <= 1 particle in cell or
  //   all particles having COM velocity, set via mixture properties

  double *vstream = d_tasks(i).vstream;
  if (np) {
    vstream[0] = mv[0] / masstot;
    vstream[1] = mv[1] / masstot;
    vstream[2] = mv[2] / masstot;
  } else vstream[0] = vstream[1] = vstream[2] = 0.0;

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
    if (temp_thermal_cell > TEMPLIMIT)
      Kokkos::atomic_max(&d_tempmax(),temp_thermal_cell);

    if (np) {
      const int ndim = d_tasks(i).ndim;
      const double sign = d_tasks(i).normal[ndim];
      vstream[ndim] += sign *
        (psubsonic - press_cell) / (massrho_cell*soundspeed_cell);
    }

    for (int m = 0; m < nspecies; m++) {
      const int ispecies = d_mspecies[m];
      d_vscale(i,m) = sqrt(2.0 * boltz * temp_thermal_cell /
                           d_species[ispecies].mass);
    }
  }

  d_tasks(i).temp_thermal = temp_thermal_cell;
  d_tasks(i).temp_rot = d_tasks(i).temp_vib = temp_thermal_cell;
}

/* ----------------------------------------------------------------------
   device version of FixEmit::mol_inflow()
   calculate flux of particles of a species with vscale/fraction
     entering a grid cell
   see comments for FixEmit::mol_inflow()
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double FixEmitFaceKokkos::mol_inflow_kokkos(double indot, double vscale,
                                            double fraction) const
{
  const double scosine = indot / vscale;
  if (scosine < -3.0) return 0.0;
  const double inward_number_flux = vscale*fraction *
    (exp(-scosine*scosine) + MY_PIS*scosine*(1.0 + erf(scosine))) /
    (2*MY_PIS);
  return inward_number_flux;
}

/* ----------------------------------------------------------------------
   grow task list
------------------------------------------------------------------------- */

void FixEmitFaceKokkos::grow_task()
{
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

  if (subsonic_style == PONLY) {
    k_vscale.modify_host(); // force resize on host
    k_vscale.sync_host();
    k_vscale.resize(ntaskmax,nspecies);
    d_vscale = k_vscale.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].vscale = &k_vscale.view_host()(i,0);
  }
}

/* ----------------------------------------------------------------------
   reallocate nspecies arrays
------------------------------------------------------------------------- */

void FixEmitFaceKokkos::realloc_nspecies()
{
  if (perspecies) {
    k_ntargetsp = DAT::tdual_float_2d_lr("emit/face:ntargetsp",ntaskmax,nspecies);
    d_ntargetsp = k_ntargetsp.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].ntargetsp = &k_ntargetsp.view_host()(i,0);
  }
  if (subsonic_style == PONLY) {
    k_vscale = DAT::tdual_float_2d_lr("emit/face:vscale",ntaskmax,nspecies);
    d_vscale = k_vscale.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].vscale = &k_vscale.view_host()(i,0);
  }
}
