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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "create_particles_kokkos.h"
#include "update.h"
#include "particle_kokkos.h"
#include "mixture.h"
#include "grid_kokkos.h"
#include "modify.h"
#include "comm.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos_type.h"
#include "sparta_masks.h"
#include "kokkos.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // same as Grid
enum{INT,DOUBLE};                       // several files

#define MAXATTEMPT 1024      // max attempts to insert a particle into cut/split cell
#define EPSZERO 1.0e-14

CreateParticlesKokkos::CreateParticlesKokkos(SPARTA* spa):
  CreateParticles(spa)
{
}

void CreateParticlesKokkos::create_local(bigint np)
{
  int dimension = domain->dimension;

  int me = comm->me;
  RanKnuth *random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;
  GridKokkos* grid_kk = ((GridKokkos*)grid);
  grid_kk->sync(Host,CINFO_MASK|CELL_MASK|SINFO_MASK);
  int nglocal = grid->nlocal;

  // flowvol = total weighted flow volume of all cells
  //   skip cells inside surfs and split cells
  //   skip cells outside defined region
  // insertvol = subset of flowvol for cells eligible for insertion
  //   insertvol = flowvol if cutflag = 1
  //   insertvol < flowvol possible if cutflag = 0 (no cut cells)

  double flowvolme = 0.0;
  double insertvolme = 0.0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].type == INSIDE) continue;
    if (cells[icell].nsplit > 1) continue;
    if (region && region->bboxflag &&
        outside_region(dimension,cells[icell].lo,cells[icell].hi))
      continue;

    flowvolme += cinfo[icell].volume / cinfo[icell].weight;
    if (!cutflag && cells[icell].nsurf) continue;
    insertvolme += cinfo[icell].volume / cinfo[icell].weight;
  }

  // calculate total Np if not set explicitly
  // based on total flowvol and mixture density

  if (np == 0) {
    double flowvol;
    MPI_Allreduce(&flowvolme,&flowvol,1,MPI_DOUBLE,MPI_SUM,world);
    np = particle->mixture[imix]->nrho * flowvol / update->fnum;
  }

  // gather cummulative insertion volumes across all procs

  double volupto;
  MPI_Scan(&insertvolme,&volupto,1,MPI_DOUBLE,MPI_SUM,world);

  double *vols;
  int nprocs = comm->nprocs;
  memory->create(vols,nprocs,"create_particles:vols");
  MPI_Allgather(&volupto,1,MPI_DOUBLE,vols,1,MPI_DOUBLE,world);

  // gathered Scan results not guaranteed to be monotonically increasing
  //   can cause epsilon mis-counts for huge particle counts
  //   so enforce monotonic increase by brute force

  for (int i = 1; i < nprocs; i++)
    if (vols[i] != vols[i-1] &&
        fabs(vols[i]-vols[i-1])/vols[nprocs-1] < EPSZERO)
      vols[i] = vols[i-1];

  // nme = # of particles for me to create
  // based on fraction of insertvol I own
  // loop over procs insures sum of nme = Np

  bigint nstart,nstop;
  if (me > 0) nstart = static_cast<bigint> (np * (vols[me-1]/vols[nprocs-1]));
  else nstart = 0;
  nstop = static_cast<bigint> (np * (vols[me]/vols[nprocs-1]));
  bigint nme = nstop-nstart;

  memory->destroy(vols);

  // nfix_update_custom = # of fixes with update_custom() method

  modify->list_init_fixes();
  int nfix_update_custom = modify->n_update_custom;

  // loop over cells I own
  // only add particles to cells eligible for insertion
  // ntarget = floating point # of particles to create in one cell
  // npercell = integer # of particles to create in one cell
  // basing ntarget on accumulated volume and nprev insures Nme total creations
  // particle species = random value based on mixture fractions
  // particle velocity = stream velocity + thermal velocity

  double nrho = particle->mixture[imix]->nrho;
  int *species = particle->mixture[imix]->species;
  double *cummulative = particle->mixture[imix]->cummulative;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;
  int nspecies = particle->mixture[imix]->nspecies;
  double temp_thermal = particle->mixture[imix]->temp_thermal;
  double temp_rot = particle->mixture[imix]->temp_rot;
  double temp_vib = particle->mixture[imix]->temp_vib;

  int npercell,ncreate,isp,ispecies,id,pflag,subcell;
  double x[3],v[3],xcell[3],vstream_var[3];
  double ntarget,scale,rn,vn,vr,theta1,theta2,erot,evib;
  double *lo,*hi;

  double *cummulative_custom = new double[nspecies];

  double tempscale = 1.0;
  double sqrttempscale = 1.0;

  double volsum = 0.0;
  bigint nprev = 0;

  // first pass, just calculate # of particles to create
  // ncreate_values[icell] = # of particles to create in ICELL

  Kokkos::View<int*, DeviceType> d_ncreate_values("create_particles:ncreate", nglocal);
  auto h_ncreate_values = Kokkos::create_mirror_view(d_ncreate_values);

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].type == INSIDE) continue;
    if (cells[icell].nsplit > 1) continue;
    if (cinfo[icell].volume == 0.0) continue;
    if (region && region->bboxflag &&
        outside_region(dimension,cells[icell].lo,cells[icell].hi))
      continue;
    if (!cutflag && cells[icell].nsurf) continue;

    volsum += cinfo[icell].volume / cinfo[icell].weight;

    double ntarget = nme * volsum/insertvolme - nprev;
    auto npercell = static_cast<int> (ntarget);

    if (random->uniform() < ntarget-npercell) npercell++;
    auto ncreate = npercell;

    lo = cells[icell].lo;
    hi = cells[icell].hi;

    if (nrho_flag) {
      if (nrho_var_flag) scale = nrho_variable(lo,hi);
      else scale = nrho_custom[icell] / nrho;
      if (scale < 0.0)
        error->one(FLERR,"Variable/custom density scale factor < 0.0");
      ntarget *= scale;
      ncreate = static_cast<int> (ntarget);
      if (random->uniform() < ntarget-ncreate) ncreate++;
    }

    h_ncreate_values[icell] = ncreate;

    // increment count without effect of density variation
    // so that target insertion count is undisturbed

    nprev += npercell;
  }
  Kokkos::deep_copy(d_ncreate_values, h_ncreate_values);

  // second pass, create particles using ncreate_values

  double *vstream_update_custom = vstream;

  int ncands;
  auto d_cells2cands = offset_scan(d_ncreate_values, ncands);
  auto h_cells2cands = Kokkos::create_mirror_view(d_cells2cands);
  Kokkos::deep_copy(h_cells2cands, d_cells2cands);

  Kokkos::View<int*, DeviceType> d_keep("cand_keep", ncands);
  Kokkos::View<int*, DeviceType> d_isp("cand_isp", ncands);
  Kokkos::View<double*[3], DeviceType> d_x("cand_x", ncands);
  Kokkos::View<double*, DeviceType> d_erot("cand_erot", ncands);
  Kokkos::View<double*, DeviceType> d_evib("cand_evib", ncands);
  Kokkos::View<int*, DeviceType> d_id("cand_id", ncands);
  Kokkos::View<double*[3], DeviceType> d_v("cand_v", ncands);
  auto h_keep = Kokkos::create_mirror_view(d_keep);
  auto h_isp = Kokkos::create_mirror_view(d_isp);
  auto h_x = Kokkos::create_mirror_view(d_x);
  auto h_erot = Kokkos::create_mirror_view(d_erot);
  auto h_evib = Kokkos::create_mirror_view(d_evib);
  auto h_id = Kokkos::create_mirror_view(d_id);
  auto h_v = Kokkos::create_mirror_view(d_v);

  for (int icell = 0; icell < nglocal; icell++) {
    if (cinfo[icell].type == INSIDE) continue;
    if (cells[icell].nsplit > 1) continue;
    if (cinfo[icell].volume == 0.0) continue;
    if (region && region->bboxflag &&
        outside_region(dimension,cells[icell].lo,cells[icell].hi))
      continue;
    if (!cutflag && cells[icell].nsurf) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;

    auto ncreate = h_ncreate_values(icell);

    if (fractions_custom_flag)
      fractions_to_cummulative(nspecies,fractions_custom[icell],cummulative_custom);

    // if surfs in cell, use xcell for all created particle attempts

    if (cells[icell].nsurf)
      pflag = grid->point_outside_surfs(icell,xcell);

    for (int m = 0; m < ncreate; m++) {

      // generate random position X for new particle

      double x[3];
      x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
      x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
      x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
      if (dimension == 2) x[2] = 0.0;

      // if surfs, check if random position is in flow region
      // if subcell, also check if in correct subcell
      // if not, attempt new insertion up to MAXATTEMPT times

      if (cells[icell].nsurf && pflag) {
        int nattempt = 0;
        while (nattempt < MAXATTEMPT) {
          if (grid->outside_surfs(icell,x,xcell)) {
            if (cells[icell].nsplit == 1) break;
            int splitcell = sinfo[cells[icell].isplit].icell;
            if (dimension == 2) subcell = update->split2d(splitcell,x);
            else subcell = update->split3d(splitcell,x);
            if (subcell == icell) break;
          }

          nattempt++;

          x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
          x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
          x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
          if (dimension == 2) x[2] = 0.0;
        }

        // particle insertion was unsuccessful

        if (nattempt >= MAXATTEMPT) continue;
      }

      // if region defined, skip if particle outside region

      if (region && !region->match(x)) continue;

      // insertion of particle at position X is accepted
      // calculate all other particle properties

      rn = random->uniform();

      if (species_flag) {
        if (species_var_flag) {
          isp = species_variable(x) - 1;
          if (isp < 0 || isp >= nspecies) continue;
          ispecies = species[isp];
        } else {
          isp = 0;
          while (cummulative_custom[isp] < rn) isp++;
          ispecies = species[isp];
        }
      } else {
        isp = 0;
        while (cummulative[isp] < rn) isp++;
        ispecies = species[isp];
      }

      h_keep(cand) = 1;
      h_isp(cand) = isp;
      for (int d = 0; d < 3; ++d) h_x(cand, d) = x[d];

      if (temp_flag) {
        if (temp_var_flag) tempscale = temperature_variable(x);
        else tempscale = temp_custom[icell] / temp_thermal;
        if (tempscale < 0.0)
          error->one(FLERR,"Variable/custom temperature scale factor < 0.0");
        sqrttempscale = sqrt(tempscale);
      }

      auto vn = vscale[isp] * sqrttempscale * sqrt(-log(random->uniform()));
      auto vr = vscale[isp] * sqrttempscale * sqrt(-log(random->uniform()));
      auto theta1 = MY_2PI * random->uniform();
      auto theta2 = MY_2PI * random->uniform();

      //these two functions also use random variables...
      h_erot(cand) = particle->erot(ispecies,temp_rot*tempscale,random);
      h_evib(cand) = particle->evib(ispecies,temp_vib*tempscale,random);

      h_id(cand) = MAXSMALLINT*random->uniform();

      double v[3];
      if (vstream_flag) {
        if (vstream_var_flag) {
          vstream_variable(x,vstream,vstream_var);
          v[0] = vstream_var[0] + vn*cos(theta1);
          v[1] = vstream_var[1] + vr*cos(theta2);
          v[2] = vstream_var[2] + vr*sin(theta2);
          vstream_update_custom = vstream_var;
        } else {
          v[0] = vstream_custom[icell][0] + vn*cos(theta1);
          v[1] = vstream_custom[icell][1] + vr*cos(theta2);
          v[2] = vstream_custom[icell][2] + vr*sin(theta2);
          vstream_update_custom = vstream_custom[icell];
        }
      } else {
        v[0] = vstream[0] + vn*cos(theta1);
        v[1] = vstream[1] + vr*cos(theta2);
        v[2] = vstream[2] + vr*sin(theta2);
      }
      for (int d = 0; d < 3; ++d) h_v(cand, d) = v[d];
    }
  }

  Kokkos::deep_copy(d_keep, h_keep);
  Kokkos::deep_copy(d_isp, h_isp);
  Kokkos::deep_copy(d_x, h_x);
  Kokkos::deep_copy(d_erot, h_erot);
  Kokkos::deep_copy(d_evib, h_evib);
  Kokkos::deep_copy(d_id, h_id);
  Kokkos::deep_copy(d_v, h_v);

  int nnew;
  auto d_cands2new = offset_scan(d_keep, nnew);

  auto particleKK = dynamic_cast<ParticleKokkos*>(particle);
  particleKK->grow(nnew);
  particleKK->sync(Device, PARTICLE_MASK | SPECIES_MASK);
  auto d_particles = particleKK->k_particles.d_view;
  Kokkos::View<int*, DeviceType> d_species("species", nspecies);
  auto h_species = Kokkos::create_mirror_view(d_species);
  for (int i = 0; i < nspecies; ++i) h_species(i) = species[i];
  Kokkos::deep_copy(d_species, h_species);
  auto nlocal_before = particleKK->nlocal;

  Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int icell) {
    auto ncreate = d_ncreate_values(icell);
    for (int m = 0; m < ncreate; m++) {
      auto cand = d_cells2cands(icell) + m;
      if (!d_keep(cand)) continue;
      auto inew = d_cands2new(cand) + nlocal_before;
      auto id = d_id(cand);
      auto ispecies = d_species(d_isp(cand));
      double x[3],v[3];
      for (int d = 0; d < 3; ++d) x[d] = d_x(cand, d);
      for (int d = 0; d < 3; ++d) v[d] = d_v(cand, d);
      auto erot = d_erot(cand);
      auto evib = d_evib(cand);
      ParticleKokkos::add_particle_kokkos(d_particles,inew,id,ispecies,icell,x,v,erot,evib);
    }
  });
  particleKK->modify(Device,PARTICLE_MASK);
  particleKK->sync(Host,PARTICLE_MASK);
  particleKK->nlocal += nnew;

  auto h_cands2new = Kokkos::create_mirror_view(d_cands2new);
  Kokkos::deep_copy(h_cands2new, d_cands2new);

  for (int icell = 0; icell < nglocal; icell++) {
    auto ncreate = h_ncreate_values(icell);
    for (int m = 0; m < ncreate; m++) {
      auto cand = h_cells2cands(icell) + m;
      if (!h_keep(cand)) continue;
      auto inew = h_cands2new(cand) + nlocal_before;

      // tempscale and vstream_update_custom are set appropriately
      // if using per-grid variables or per-grid custom attributes

      if (nfix_update_custom)
        modify->update_custom(particle->nlocal-1,tempscale*temp_thermal,
                              tempscale*temp_rot,tempscale*temp_vib,
                              vstream_update_custom);
    }
  }

  // clean up

  delete [] cummulative_custom;
  delete random;
}

