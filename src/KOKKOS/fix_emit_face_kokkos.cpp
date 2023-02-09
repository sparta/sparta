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
#include "sparta_masks.h"
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
  particle_kk_copy(sparta)
{
  kokkos_flag = 1;
  execution_space = Device;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

FixEmitFaceKokkos::~FixEmitFaceKokkos()
{
  if (copymode) return;

  particle_kk_copy.uncopy();

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
  k_species     = DAT::tdual_int_1d("species", nspecies);

  d_mix_vscale  = k_mix_vscale .d_view;
  d_cummulative = k_cummulative.d_view;
  d_species     = k_species    .d_view;

  auto h_mix_vscale  = k_mix_vscale .h_view;
  auto h_cummulative = k_cummulative.h_view;
  auto h_species     = k_species    .h_view;

  for (int isp = 0; isp < nspecies; ++isp) {
    h_mix_vscale(isp) = particle->mixture[imix]->vscale[isp];
    h_cummulative(isp) = particle->mixture[imix]->cummulative[isp];
    h_species(isp) = particle->mixture[imix]->species[isp];
  }

  k_mix_vscale .modify_host();
  k_cummulative.modify_host();
  k_species    .modify_host();
}

/* ----------------------------------------------------------------------
   create tasks for one grid cell
   add them to tasks list and increment ntasks
------------------------------------------------------------------------- */

void FixEmitFaceKokkos::create_task(int icell)
{
  FixEmitFace::create_task(icell);
  k_tasks.modify_host();
  k_ntargetsp.modify_host();
  k_vscale.modify_host();
}

/* ---------------------------------------------------------------------- */

void FixEmitFaceKokkos::perform_task()
{
  dt = update->dt;
  auto l_dimension = this->dimension;
  auto l_subsonic_style = this->subsonic_style;

  // if subsonic, re-compute particle inflow counts for each task
  // also computes current per-task temp_thermal and vstream

  if (subsonic)
    error->one(FLERR,"Cannot yet use fix emit/face/kk with subsonic emission");
  //if (subsonic) subsonic_inflow(); ////////////////////////

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
  k_species    .sync_device();
  k_cummulative.sync_device();

  auto ld_mix_vscale = d_mix_vscale;
  auto ld_species    = d_species   ;

  ParticleKokkos* particle_kk = ((ParticleKokkos*)particle);
  particle_kk->update_class_variables();
  particle_kk_copy.copy(particle_kk);

  if (region)
    error->one(FLERR,"Cannot yet use fix emit/face/kk with regions");

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
  auto ld_particles = particleKK->k_particles.d_view;

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

    auto ispecies = ld_species(isp);

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
      auto ntarget = d_ntargetsp(i,isp)+rand_gen.drand();
      ninsert = static_cast<int> (ntarget);
      d_ninsert(i * nspecies + isp) = ninsert;
    }
  } else {
    if (np == 0) {
      auto ntarget = d_tasks(i).ntarget+rand_gen.drand();
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

      auto ispecies = d_species[isp];
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
      auto ispecies = d_species[isp];
      auto scosine = indot / vscale_val;

      double x[3];
      x[0] = lo[0] + rand_gen.drand() * (hi[0]-lo[0]);
      x[1] = lo[1] + rand_gen.drand() * (hi[1]-lo[1]);
      if (dimension == 3) x[2] = lo[2] + rand_gen.drand() * (hi[2]-lo[2]);
      else x[2] = 0.0;

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
   grow task list
------------------------------------------------------------------------- */

void FixEmitFaceKokkos::grow_task()
{
  ntaskmax += DELTATASK;

  if (tasks == NULL)
    k_tasks = tdual_task_1d("emit/face:tasks",ntaskmax);
  else {
    k_tasks.sync_host();
    k_tasks.modify_host(); // force resize on host
    k_tasks.resize(ntaskmax);
  }
  d_tasks = k_tasks.d_view;
  tasks = k_tasks.h_view.data();

  // set all new task bytes to 0 so valgrind won't complain
  // if bytes between fields are uninitialized

  //memset(&tasks[oldmax],0,(ntaskmax-oldmax)*sizeof(Task));

  // allocate vectors in each new task or set to NULL

  if (perspecies) {
    k_ntargetsp.sync_host();
    k_ntargetsp.modify_host(); // force resize on host
    k_ntargetsp.resize(ntaskmax,nspecies);
    d_ntargetsp = k_ntargetsp.d_view;
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].ntargetsp = k_ntargetsp.h_view.data() + i*k_ntargetsp.h_view.extent(1);
  }

  if (subsonic_style == PONLY) {
    k_vscale.modify_host(); // force resize on host
    k_vscale.sync_host();
    k_vscale.resize(ntaskmax,nspecies);
    d_vscale = k_vscale.d_view;
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].vscale = k_vscale.h_view.data() + i*k_vscale.h_view.extent(1);
  }
}

/* ----------------------------------------------------------------------
   reallocate nspecies arrays
------------------------------------------------------------------------- */

void FixEmitFaceKokkos::realloc_nspecies()
{
  if (perspecies) {
    k_ntargetsp = DAT::tdual_float_2d_lr("emit/face:ntargetsp",ntaskmax,nspecies);
    d_ntargetsp = k_ntargetsp.d_view;
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].ntargetsp = k_ntargetsp.h_view.data() + i*k_ntargetsp.h_view.extent(1);
  }
  if (subsonic_style == PONLY) {
    k_vscale = DAT::tdual_float_2d_lr("emit/face:vscale",ntaskmax,nspecies);
    d_vscale = k_vscale.d_view;
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].vscale = k_vscale.h_view.data() + i*k_vscale.h_view.extent(1);
  }
}
