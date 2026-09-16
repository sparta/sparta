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
#include "fix_emit_face_file_kokkos.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "grid.h"
#include "surf.h"
#include "particle.h"
#include "mixture.h"
#include "comm.h"
#include "modify.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos_type.h"
#include "particle_kokkos.h"
#include "grid_kokkos.h"
#include "fix_emit_kokkos.h"
#include "sparta_masks.h"
#include "Kokkos_Random.hpp"

using namespace SPARTA_NS;
using namespace MathConst;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // same as Grid
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Grid
enum{NRHO,TEMP_THERMAL,TEMP_ROT,TEMP_VIB,VX,VY,VZ,PRESS,SPECIES};
enum{NOSUBSONIC,PTBOTH,PONLY};

#define DELTATASK 256
#define TEMPLIMIT 1.0e5

/* ----------------------------------------------------------------------
   insert particles on a boundary face, with per-face flow properties
     interpolated from a file, on the device
------------------------------------------------------------------------- */

FixEmitFaceFileKokkos::FixEmitFaceFileKokkos(SPARTA *sparta, int narg, char **arg) :
  FixEmitFaceFile(sparta, narg, arg),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            ),
  particle_kk_copy(sparta)
{
  kokkos_flag = 1;
  execution_space = Device;

  // this fix runs at START_OF_STEP.  ModifyKokkos::start_of_step() calls
  //   ParticleKokkos::sync(execution_space,datamask_read) before the fix and
  //   ::modify(execution_space,datamask_modify) after it.  Both masks are
  //   EMPTY here because this fix does all of its own syncing: it needs
  //   particles on the device only when subsonic, and it appends to the
  //   particle list itself (grow + modify(Device,PARTICLE_MASK)) inside
  //   perform_task().  Declaring PARTICLE_MASK here would force a
  //   host<->device round trip of the whole particle list on every step of
  //   every run that uses this fix.

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  region_flag = 0;
  nregion_token = 0;
  axisymmetric = 0;
  dt_step = 0.0;
  boltz = 0.0;
  plist_descending = 0;
}

/* ---------------------------------------------------------------------- */

FixEmitFaceFileKokkos::~FixEmitFaceFileKokkos()
{
  if (copymode) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
#endif

  // the per-task vectors point into Kokkos DualViews, not into new[] memory,
  //   so hide them from ~FixEmitFaceFile(), which would delete[] them

  if (tasks) {
    for (int i = 0; i < ntaskmax; i++) {
      tasks[i].ntargetsp = NULL;
      tasks[i].vscale = NULL;
      tasks[i].fraction = NULL;
      tasks[i].cummulative = NULL;
    }
  }

  // tasks itself is the host half of k_tasks, so it must not reach
  //   memory->sfree() in ~FixEmitFaceFile().  zero ntaskmax as well: unlike
  //   ~FixEmitFace(), ~FixEmitFaceFile() has no "if (tasks)" guard around its
  //   delete[] loop, so a NULL tasks with a nonzero ntaskmax would fault

  tasks = NULL;
  ntaskmax = 0;
}

/* ---------------------------------------------------------------------- */

void FixEmitFaceFileKokkos::init()
{
  // fix emit/face/file supports a one-pass and a two-pass insertion loop.
  //   only the two-pass one draws every task's insertion count before it
  //   generates any particle, which is the order the Kokkos kernel pair is
  //   forced into: the offset scan has to know all the counts before the
  //   candidate arrays can be sized.  So under SPARTA_KOKKOS_EXACT, where
  //   the whole point is a bit-for-bit match against the non-Kokkos run,
  //   refuse to run without it.  Outside SPARTA_KOKKOS_EXACT stay silent,
  //   the same as fix emit/face/kk: the random streams differ anyway.

#ifdef SPARTA_KOKKOS_EXACT
  if (!twopass)
    error->all(FLERR,"Fix emit/face/file/kk requires the twopass keyword "
               "under SPARTA_KOKKOS_EXACT: without it the Kokkos and "
               "non-Kokkos runs consume random numbers in a different order "
               "and will not produce the same particles");
#endif

  // pull anything the device wrote (subsonic updates the tasks and vscale)
  //   back to the host before host-side code reads or overwrites it, and
  //   before any of these DualViews is replaced below

  k_tasks.sync_host();
  if (perspecies) k_ntargetsp.sync_host();
  k_vscale.sync_host();
  k_cummulative.sync_host();
  k_fraction.sync_host();

  // FixEmitFaceFile::init() delete[]s and re-new[]s tasks[i].fraction,
  //   .cummulative, .vscale and .ntargetsp for the ntask tasks left over from
  //   the previous run, because the mixture's species count may have changed.
  //   Ours are not new[] memory -- they alias rows of the DualViews -- so
  //   that loop must not run on them.  FixEmitFaceFile has no
  //   realloc_nspecies() hook the way FixEmitFace does, so do the equivalent
  //   here instead:
  //     - pick up the new species count (init() sets the same value again)
  //     - resize the DualViews and re-aim every Task pointer at their rows
  //     - zero ntask, which makes those two base loops no-ops.  Nothing else
  //       in init() reads ntask, and FixEmit::create_tasks() -- called at the
  //       end of init() -- zeroes it itself and rebuilds the whole list, so
  //       this changes nothing for the host.

  nspecies = particle->mixture[imix]->nspecies;
  realloc_species_views();
  ntask = 0;

  FixEmitFaceFile::init();

  // create_tasks() ran inside init() and wrote the task values, plus the
  //   per-task fraction/cummulative/vscale/ntargetsp rows, on the host

  k_tasks.modify_host();
  if (perspecies) k_ntargetsp.modify_host();
  k_vscale.modify_host();
  k_cummulative.modify_host();
  k_fraction.modify_host();

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif

  // domain->axisymmetric is not reachable from a device kernel

  axisymmetric = domain->axisymmetric;

  // mixture species indices, the only mixture-wide array the kernels need.
  // vscale/fraction/cummulative are per-task for this style, since they are
  //   interpolated from the file per face, so there is no mixture-wide copy

  k_mspecies = DAT::tdual_int_1d("emit/face/file:mspecies",nspecies);
  d_mspecies = k_mspecies.view_device();

  auto h_mspecies = k_mspecies.view_host();
  for (int isp = 0; isp < nspecies; isp++)
    h_mspecies(isp) = particle->mixture[imix]->species[isp];

  k_mspecies.modify_host();
  // the region is fixed for the run; flatten it once here rather than
  //   rebuilding and re-uploading the token stream in perform_task()

  flatten_region();

}

/* ----------------------------------------------------------------------
   create tasks for all grid cells
   the interpolation of the file mesh onto each cell face is host-only setup:
     it happens here, once per grid change, never per step.  the task values
     it writes land directly in the host half of the DualViews (see
     realloc_species_views()), so all this override has to do is mark them
     dirty
------------------------------------------------------------------------- */

void FixEmitFaceFileKokkos::create_tasks()
{
  k_tasks.sync_host();
  if (perspecies) k_ntargetsp.sync_host();
  k_vscale.sync_host();
  k_cummulative.sync_host();
  k_fraction.sync_host();

  FixEmitFaceFile::create_tasks();

  k_tasks.modify_host();
  if (perspecies) k_ntargetsp.modify_host();
  k_vscale.modify_host();
  k_cummulative.modify_host();
  k_fraction.modify_host();
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

void FixEmitFaceFileKokkos::flatten_region()
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

/* ----------------------------------------------------------------------
   insert particles in grid cells with faces touching the inflow boundary
   always two-pass: count kernel -> offset scan -> generate kernel ->
     compaction
------------------------------------------------------------------------- */

void FixEmitFaceFileKokkos::perform_task()
{
  // the non-Kokkos perform_task_*() read update->dt into a local, leaving the
  //   class member dt (frozen at init) for subsonic_inflow().  keep that
  //   split so a run with fix dt/reset behaves identically

  dt_step = update->dt;

  // face geometry is fix-wide here, not per task as in fix emit/face, so
  //   hoist it into locals the compaction lambda can capture by value
  //   (a device lambda may not capture this)

  auto l_dimension = this->dimension;
  auto l_ndim = this->ndim;
  auto l_pdim = this->pdim;
  auto l_qdim = this->qdim;
  const double l_normal_ndim = normal[ndim];

  // if subsonic, re-compute particle inflow counts for each task
  // also computes current per-task temp_thermal and vstream

  if (subsonic) subsonic_inflow();

  // insert particles for each task = cell/face pair
  // ntarget/ninsert is either perspecies or for all species

  // copy needed task data to device

  if (perspecies) k_ntargetsp.sync_device();
  else k_tasks.sync_device();

  auto ninsert_dim1 = perspecies ? nspecies : 1;
  if (d_ninsert.extent(0) < ntask * ninsert_dim1)
    d_ninsert = DAT::t_int_1d("emit/face/file:ninsert", ntask * ninsert_dim1);

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitFaceFile_ninsert>(0,ntask),*this);
  copymode = 0;

  int ncands;
  d_task2cand = offset_scan(d_ninsert, ncands);

  if (ncands == 0) return;

  // for one particle:
  //   x = random position on subset of face that overlaps with file grid
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
    d_x = DAT::t_float_2d("emit/face/file:x", ncands, l_dimension);

  if (d_task.extent(0) < ncands) {
    d_beta_un  = DAT::t_float_1d("emit/face/file:beta_un", ncands);
    d_theta    = DAT::t_float_1d("emit/face/file:theta", ncands);
    d_vr       = DAT::t_float_1d("emit/face/file:vr", ncands);
    d_erot     = DAT::t_float_1d("emit/face/file:erot", ncands);
    d_evib     = DAT::t_float_1d("emit/face/file:evib", ncands);
    d_dtremain = DAT::t_float_1d("emit/face/file:dtremain", ncands);
    d_id       = DAT::t_int_1d("emit/face/file:id", ncands);
    d_isp      = DAT::t_int_1d("emit/face/file:isp", ncands);
    d_task     = DAT::t_int_1d("emit/face/file:task", ncands);
    d_keep     = DAT::t_int_1d("emit/face/file:keep", ncands);
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
  // fraction/cummulative/vscale are per-task for this style, so all of them
  //   ride along with the tasks rather than coming from the mixture

  k_tasks.sync_device();
  if (perspecies) k_ntargetsp.sync_device();
  k_vscale.sync_device();
  k_cummulative.sync_device();

  auto ld_tasks = d_tasks;
  auto ld_vscale = d_vscale;

  k_mspecies.sync_device();
  auto ld_mspecies = d_mspecies;

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
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixEmitFaceFile_perform_task>(0,ntask),*this,nsingle_reduce);
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

    auto i = ld_task(cand);
    Task task_i = ld_tasks(i);

    const int pcell = task_i.pcell;
    double *vstream = task_i.vstream;

    auto isp = ld_isp(cand);
    auto vscale_val = ld_vscale(i, isp);
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
    v[l_ndim] = beta_un*vscale_val*l_normal_ndim + vstream[l_ndim];
    v[l_pdim] = vr * sin(theta) + vstream[l_pdim];
    v[l_qdim] = vr * cos(theta) + vstream[l_qdim];

    auto inew = ld_cands2new(cand);
    auto ilocal = nlocal_before + inew;

    ParticleKokkos::add_particle_kokkos(ld_particles,ilocal,
        id,ispecies,pcell,x,v,erot,evib);

    ld_particles(ilocal).flag = PINSERT;
    ld_particles(ilocal).dtremain = dtremain;
  });

  particleKK->nlocal = nlocal_before + nnew;
  particleKK->modify(SPARTA_NS::Device, PARTICLE_MASK);
  particleKK->zero_custom_kokkos(nlocal_before,particleKK->nlocal);

  // custom per-particle attributes are still a host-side callback

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
      auto temp_elec = tasks[task].temp_elec;
      auto vstream = tasks[task].vstream;

      auto inew = h_cands2new(cand);
      auto ilocal = nlocal_before + inew;

      modify->update_custom(ilocal,temp_thermal,temp_rot,temp_vib,
                            temp_elec,vstream);
    }
  }
}

/* ----------------------------------------------------------------------
   # of particles to insert for each task
   fix emit/face/file has no "n" and no "modulate" option, so unlike
     fix emit/face this is just ntarget + a uniform deviate
   this is the first of the two passes: every count is drawn here, before
     any particle is generated, which is what perform_task_twopass() on the
     host mirrors
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixEmitFaceFileKokkos::operator()(TagFixEmitFaceFile_ninsert, const int &i) const
{
  rand_type rand_gen = rand_pool.get_state();

  if (perspecies) {
    for (int isp = 0; isp < nspecies; isp++) {
      const double ntarget = d_ntargetsp(i,isp) + rand_gen.drand();
      d_ninsert(i * nspecies + isp) = static_cast<int> (ntarget);
    }
  } else {
    const double ntarget = d_tasks(i).ntarget + rand_gen.drand();
    d_ninsert(i) = static_cast<int> (ntarget);
  }

  rand_pool.free_state(rand_gen);
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixEmitFaceFileKokkos::operator()(TagFixEmitFaceFile_perform_task,
                                       const int &i, int &nsingle_reduce) const
{
  rand_type rand_gen = rand_pool.get_state();

  Task task_i = d_tasks(i);

  double *lo = task_i.lo;
  double *hi = task_i.hi;
  double *vstream = task_i.vstream;

  const double temp_rot = task_i.temp_rot;
  const double temp_vib = task_i.temp_vib;

  // normal is fix-wide for this style, not per task

  const double indot = vstream[0]*normal[0] + vstream[1]*normal[1] +
    vstream[2]*normal[2];

  if (perspecies) {
    for (int isp = 0; isp < nspecies; isp++) {
      const int ispecies = d_mspecies[isp];

      // per-task vscale, not the mixture's: the file can set a per-face
      //   temperature, which interpolate() folded into the task

      const double vscale_val = d_vscale(i,isp);
      const double scosine = indot / vscale_val;

      const int ninsert = d_ninsert(i * nspecies + isp);
      const int start = d_task2cand(i * nspecies + isp);

      int nactual = 0;
      for (int m = 0; m < ninsert; m++) {
        const int cand = start + m;

        double x[3];
        x[0] = lo[0] + rand_gen.drand() * (hi[0]-lo[0]);
        if (axisymmetric)
          x[1] = sqrt(lo[1]*lo[1] +
                      rand_gen.drand() * (hi[1]*hi[1]-lo[1]*lo[1]));
        else x[1] = lo[1] + rand_gen.drand() * (hi[1]-lo[1]);
        if (dimension == 3) x[2] = lo[2] + rand_gen.drand() * (hi[2]-lo[2]);
        else x[2] = 0.0;

        // region_flag must be tested first: with no region there is no token
        //   stream, and region_match_kk() of an empty stream rejects
        //   everything

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

        d_beta_un(cand) = beta_un;

        d_theta(cand) = MY_2PI * rand_gen.drand();
        d_vr(cand) = vscale_val * sqrt(-log(rand_gen.drand()));
        d_erot(cand) = particle_kk_copy.obj.erot(ispecies,temp_rot,rand_gen);
        d_evib(cand) = particle_kk_copy.obj.evib(ispecies,temp_vib,rand_gen);
        d_id(cand) = MAXSMALLINT*rand_gen.drand();
        d_dtremain(cand) = dt_step * rand_gen.drand();
      }

      nsingle_reduce += nactual;
    }

  } else {
    const int ninsert = d_ninsert(i);
    const int start = d_task2cand(i);

    int nactual = 0;
    for (int m = 0; m < ninsert; m++) {
      const int cand = start + m;

      // per-task cummulative, not the mixture's: the file can set per-face
      //   species fractions, which interpolate() folded into the task

      const double rn = rand_gen.drand();
      int isp = 0;
      while (d_cummulative(i,isp) < rn) isp++;

      const int ispecies = d_mspecies[isp];
      const double vscale_val = d_vscale(i,isp);
      const double scosine = indot / vscale_val;

      double x[3];
      x[0] = lo[0] + rand_gen.drand() * (hi[0]-lo[0]);
      if (axisymmetric)
        x[1] = sqrt(lo[1]*lo[1] +
                    rand_gen.drand() * (hi[1]*hi[1]-lo[1]*lo[1]));
      else x[1] = lo[1] + rand_gen.drand() * (hi[1]-lo[1]);
      if (dimension == 3) x[2] = lo[2] + rand_gen.drand() * (hi[2]-lo[2]);
      else x[2] = 0.0;

      // region_flag must be tested first: with no region there is no token
      //   stream, and region_match_kk() of an empty stream rejects everything

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

      d_beta_un(cand) = beta_un;

      d_theta(cand) = MY_2PI * rand_gen.drand();
      d_vr(cand) = vscale_val * sqrt(-log(rand_gen.drand()));
      d_erot(cand) = particle_kk_copy.obj.erot(ispecies,temp_rot,rand_gen);
      d_evib(cand) = particle_kk_copy.obj.evib(ispecies,temp_vib,rand_gen);
      d_id(cand) = MAXSMALLINT*rand_gen.drand();
      d_dtremain(cand) = dt_step * rand_gen.drand();
    }

    nsingle_reduce += nactual;
  }

  rand_pool.free_state(rand_gen);
}

/* ----------------------------------------------------------------------
   recalculate task properties based on subsonic BC
   subsonic is enabled by a "press" column in the input file, so it is a
     property of the file, not of a keyword
------------------------------------------------------------------------- */

void FixEmitFaceFileKokkos::subsonic_inflow()
{
  // for grid cells that are part of tasks:
  // calculate local nrho, vstream, and thermal temperature
  // if needed sort particles for grid cells with tasks

  subsonic_sort();
  subsonic_grid();

  // recalculate particle insertion counts for each task
  // vscale here is recomputed from the per-task temp_thermal, exactly as the
  //   non-Kokkos subsonic_inflow() does -- it does NOT read tasks[i].vscale

  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,SPECIES_MASK);
  d_species_all = particle_kk->k_species.view_device();

  GridKokkos* grid_kk = (GridKokkos*) grid;
  grid_kk->sync(Device,CINFO_MASK);
  d_cinfo = grid_kk->k_cinfo.view_device();

  k_tasks.sync_device();
  if (perspecies) k_ntargetsp.sync_device();
  k_mspecies.sync_device();
  k_fraction.sync_device();

  boltz = update->boltz;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitFaceFile_subsonic_inflow>(0,ntask),*this);
  copymode = 0;

  k_tasks.modify_device();
  if (perspecies) k_ntargetsp.modify_device();

  // release references to reduce memory use

  d_species_all = t_species_1d();
  d_cinfo = {};
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixEmitFaceFileKokkos::operator()(TagFixEmitFaceFile_subsonic_inflow,
                                       const int &i) const
{
  double *vstream = d_tasks(i).vstream;
  const double indot = vstream[0]*normal[0] + vstream[1]*normal[1] +
    vstream[2]*normal[2];

  const double area = d_tasks(i).area;
  const double nrho = d_tasks(i).nrho;
  const double temp_thermal = d_tasks(i).temp_thermal;
  const int icell = d_tasks(i).icell;

  // fraction is per task here, since the file can set per-face species
  //   fractions -- fix emit/face uses the mixture-wide vector instead

  double ntarget = 0.0;
  for (int isp = 0; isp < nspecies; isp++) {
    const double mass = d_species_all[d_mspecies[isp]].mass;
    const double vscale = sqrt(2.0 * boltz * temp_thermal / mass);
    double ntargetsp = mol_inflow_kokkos(indot,vscale,d_fraction(i,isp));
    ntargetsp *= nrho*area*dt / fnum;
    ntargetsp /= d_cinfo[icell].weight;
    ntarget += ntargetsp;
    if (perspecies) d_ntargetsp(i,isp) = ntargetsp;
  }

  d_tasks(i).ntarget = ntarget;
  if (ntarget >= MAXSMALLINT)
    Kokkos::abort("Fix emit/face/file subsonic insertion count "
                  "exceeds 32-bit int");
}

/* ----------------------------------------------------------------------
   sort particles into grid cells on device
   the non-Kokkos FixEmitFaceFile::subsonic_sort() builds its own per-cell
     linked list (Grid::ChildInfo first/count + Particle::next) for the
     "active" cells only.  on device we instead reuse the compressed per-cell
     particle lists ParticleKokkos::sort_kokkos() builds for collisions
     (GridKokkos d_plist/d_cellcount), which cover every cell.  the moment
     sums below are the only consumer, and they are per-cell, so covering
     extra cells costs nothing but the sort itself.
   this is also why the host-side activecell[] bookkeeping and the
     active_current flag are simply not used by the Kokkos path
------------------------------------------------------------------------- */

void FixEmitFaceFileKokkos::subsonic_sort()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;

  // sorted_kk mirrors the host Particle::sorted flag the non-Kokkos path
  //   tests before calling subsonic_sort().  Record it BEFORE sorting: when
  //   the non-Kokkos code has to build the list itself it pushes each
  //   particle on the head, so the list comes out in DEcreasing particle
  //   index; an already-sorted list (built by Particle::sort()'s reverse
  //   loop) comes out in INcreasing index.  d_plist is always increasing,
  //   so remember which order to walk it in.

  plist_descending = !particle_kk->sorted_kk;
  if (!particle_kk->sorted_kk) particle_kk->sort_kokkos();
}

/* ----------------------------------------------------------------------
   compute number density, thermal temperature, stream velocity
   only for grid cells associated with a task
   first compute for grid cells, then adjust due to boundary conditions
------------------------------------------------------------------------- */

void FixEmitFaceFileKokkos::subsonic_grid()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species_all = particle_kk->k_species.view_device();

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

  // only track max thermal temp until the one-time warning has fired
  // avoids a per-step device->host fence once subsonic_warning is set

  if (!subsonic_warning) {
    if (d_tempmax.data() == nullptr)
      d_tempmax = DAT::t_float_scalar("emit/face/file:tempmax");
    Kokkos::deep_copy(d_tempmax,0.0);
  }

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixEmitFaceFile_subsonic_grid>(0,ntask),*this);
  copymode = 0;

  k_tasks.modify_device();
  if (subsonic_style == PONLY) k_vscale.modify_device();

  // release references to reduce memory use

  d_particles = t_particle_1d();
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

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void FixEmitFaceFileKokkos::operator()(TagFixEmitFaceFile_subsonic_grid,
                                       const int &i) const
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
    const int ispecies = d_particles[ip].ispecies;
    const double mass = d_species_all[ispecies].mass;
    const double *v = d_particles[ip].v;
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

  double *vstream = d_tasks(i).vstream;
  if (np) {
    vstream[0] = mv[0] / masstot;
    vstream[1] = mv[1] / masstot;
    vstream[2] = mv[2] / masstot;
  } else vstream[0] = vstream[1] = vstream[2] = 0.0;

  // press is the file-interpolated pressure for THIS task, which is what
  //   distinguishes this style from fix emit/face's global psubsonic

  const double press = d_tasks(i).press;

  double temp_thermal_cell;

  if (subsonic_style == PTBOTH) {
    d_tasks(i).nrho = press / (boltz * d_tasks(i).temp_thermal);
    temp_thermal_cell = d_tasks(i).temp_thermal;

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
      (press - press_cell) / (soundspeed_cell*soundspeed_cell);
    temp_thermal_cell = press / (boltz * d_tasks(i).nrho);
    if (!subsonic_warning && temp_thermal_cell > TEMPLIMIT)
      Kokkos::atomic_max(&d_tempmax(),temp_thermal_cell);

    // the non-Kokkos code guards this update with massrho_cell*soundspeed
    //   > 0.0 as well as np, unlike fix emit/face.  keep the extra guard

    if (np && massrho_cell*soundspeed_cell > 0.0) {
      const double sign = normal[ndim];
      vstream[ndim] += sign *
        (press - press_cell) / (massrho_cell*soundspeed_cell);
    }

    for (int m = 0; m < nspecies; m++) {
      const int ispecies = d_mspecies[m];
      d_vscale(i,m) = sqrt(2.0 * boltz * temp_thermal_cell /
                           d_species_all[ispecies].mass);
    }
  }

  d_tasks(i).temp_thermal = temp_thermal_cell;
  d_tasks(i).temp_rot = d_tasks(i).temp_vib = d_tasks(i).temp_elec =
    temp_thermal_cell;
}

/* ----------------------------------------------------------------------
   grow task list
------------------------------------------------------------------------- */

void FixEmitFaceFileKokkos::grow_task()
{
  ntaskmax += DELTATASK;

  k_tasks.sync_host();
  k_tasks.modify_host(); // force resize on host
  memoryKK->grow_kokkos(k_tasks,tasks,ntaskmax,"emit/face/file:tasks");
  d_tasks = k_tasks.view_device();

  // allocate vectors in each new task or set to NULL
  // ntargetsp is only used for perspecies, exactly as in the non-Kokkos code

  if (perspecies) {
    k_ntargetsp.sync_host();
    k_ntargetsp.modify_host(); // force resize on host
    k_ntargetsp.resize(ntaskmax,nspecies);
    d_ntargetsp = k_ntargetsp.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].ntargetsp = &k_ntargetsp.view_host()(i,0);
  } else {
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].ntargetsp = NULL;
  }

  // fraction/cummulative/vscale are per-task for every run of this style,
  //   subsonic or not, because the file can vary them face by face.  so
  //   unlike fix emit/face they are always allocated, and the Task pointers
  //   the host interpolate() writes through are aimed at their host rows

  k_vscale.sync_host();
  k_vscale.modify_host();
  k_vscale.resize(ntaskmax,nspecies);
  d_vscale = k_vscale.view_device();

  k_cummulative.sync_host();
  k_cummulative.modify_host();
  k_cummulative.resize(ntaskmax,nspecies);
  d_cummulative = k_cummulative.view_device();

  k_fraction.sync_host();
  k_fraction.modify_host();
  k_fraction.resize(ntaskmax,nspecies);
  d_fraction = k_fraction.view_device();

  for (int i = 0; i < ntaskmax; i++) {
    tasks[i].vscale = &k_vscale.view_host()(i,0);
    tasks[i].cummulative = &k_cummulative.view_host()(i,0);
    tasks[i].fraction = &k_fraction.view_host()(i,0);
  }
}

/* ----------------------------------------------------------------------
   (re)allocate the per-task species arrays and re-aim the Task pointers
   called from init(), in place of the realloc_nspecies() hook that
     FixEmitFace has and FixEmitFaceFile does not: the mixture's species
     count may have changed since the last run, and the rows have to be the
     right width before init() -> create_tasks() -> interpolate() starts
     writing through the Task pointers
   the old contents are dropped, exactly as the base class's delete[]/new[]
     pair drops them; create_tasks() rewrites every task from the file mesh
------------------------------------------------------------------------- */

void FixEmitFaceFileKokkos::realloc_species_views()
{
  if (perspecies) {
    k_ntargetsp = DAT::tdual_float_2d_lr("emit/face/file:ntargetsp",
                                         ntaskmax,nspecies);
    d_ntargetsp = k_ntargetsp.view_device();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].ntargetsp = &k_ntargetsp.view_host()(i,0);
  } else {
    k_ntargetsp = DAT::tdual_float_2d_lr();
    d_ntargetsp = DAT::t_float_2d_lr();
    for (int i = 0; i < ntaskmax; i++)
      tasks[i].ntargetsp = NULL;
  }

  k_vscale = DAT::tdual_float_2d_lr("emit/face/file:vscale",
                                    ntaskmax,nspecies);
  k_cummulative = DAT::tdual_float_2d_lr("emit/face/file:cummulative",
                                         ntaskmax,nspecies);
  k_fraction = DAT::tdual_float_2d_lr("emit/face/file:fraction",
                                      ntaskmax,nspecies);

  d_vscale = k_vscale.view_device();
  d_cummulative = k_cummulative.view_device();
  d_fraction = k_fraction.view_device();

  for (int i = 0; i < ntaskmax; i++) {
    tasks[i].vscale = &k_vscale.view_host()(i,0);
    tasks[i].cummulative = &k_cummulative.view_host()(i,0);
    tasks[i].fraction = &k_fraction.view_host()(i,0);
  }
}
