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
  GridKokkos* grid_kk = ((GridKokkos*)grid);
  grid_kk->sync(Host,CINFO_MASK|CELL_MASK);
  int nglocal = grid->nlocal;

  // volme = volume of grid cells I own that are OUTSIDE surfs
  // skip cells entirely outside region
  // Nme = # of particles I will create
  // MPI_Scan() logic insures sum of nme = Np

  double *lo,*hi;
  double volone;

  double volme = 0.0;
  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].type != OUTSIDE) continue;
    lo = cells[i].lo;
    hi = cells[i].hi;
    if (region && region->bboxflag && outside_region(dimension,lo,hi))
      continue;

    if (dimension == 3) volone = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else if (domain->axisymmetric)
      volone = (hi[0]-lo[0]) * (hi[1]*hi[1]-lo[1]*lo[1])*MY_PI;
    else volone = (hi[0]-lo[0]) * (hi[1]-lo[1]);
    volme += volone / cinfo[i].weight;
  }

  double volupto;
  MPI_Scan(&volme,&volupto,1,MPI_DOUBLE,MPI_SUM,world);

  double *vols;
  int nprocs = comm->nprocs;
  memory->create(vols,nprocs,"create_particles:vols");
  MPI_Allgather(&volupto,1,MPI_DOUBLE,vols,1,MPI_DOUBLE,world);

  // nme = # of particles for me to create
  // gathered Scan results not guaranteed to be monotonically increasing
  // can cause epsilon mis-counts for huge particle counts
  // enforce that by brute force

  for (int i = 1; i < nprocs; i++)
    if (vols[i] != vols[i-1] &&
        fabs(vols[i]-vols[i-1])/vols[nprocs-1] < EPSZERO)
      vols[i] = vols[i-1];

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
  // only add particles to OUTSIDE cells
  // skip cells entirely outside region
  // ntarget = floating point # of particles to create in one cell
  // npercell = integer # of particles to create in one cell
  // basing ntarget on accumulated volume and nprev insures Nme total creations
  // particle species = random value based on mixture fractions
  // particle velocity = stream velocity + thermal velocity

  int *species = particle->mixture[imix]->species;
  double *cummulative = particle->mixture[imix]->cummulative;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;
  int nspecies = particle->mixture[imix]->nspecies;
  double temp_thermal = particle->mixture[imix]->temp_thermal;
  double temp_rot = particle->mixture[imix]->temp_rot;
  double temp_vib = particle->mixture[imix]->temp_vib;

  double tempscale = 1.0;
  double sqrttempscale = 1.0;

  double volsum = 0.0;
  bigint nprev = 0;

  double vstream_variable[3];

  Kokkos::View<int*, DeviceType> d_npercell("npercell", nglocal);
  auto h_npercell = Kokkos::create_mirror_view(d_npercell);

  for (int i = 0; i < nglocal; i++) {
    if (cinfo[i].type != OUTSIDE) continue;
    lo = cells[i].lo;
    hi = cells[i].hi;
    if (region && region->bboxflag && outside_region(dimension,lo,hi))
      continue;

    if (dimension == 3) volone = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else if (domain->axisymmetric)
      volone = (hi[0]-lo[0]) * (hi[1]*hi[1]-lo[1]*lo[1])*MY_PI;
    else volone = (hi[0]-lo[0]) * (hi[1]-lo[1]);
    volsum += volone / cinfo[i].weight;

    double ntarget = nme * volsum/volme - nprev;
    auto npercell = static_cast<int> (ntarget);
    if (random->uniform() < ntarget-npercell) npercell++;
    auto ncreate = npercell;

    if (densflag) {
      auto scale = density_variable(lo,hi);
      ntarget *= scale;
      ncreate = static_cast<int> (ntarget);
      if (random->uniform() < ntarget-ncreate) ncreate++;
    }

    if (ncreate < 0) ncreate = 0;
    h_npercell(i) = ncreate;

    // increment count without effect of density variation
    // so that target insertion count is undisturbed

    nprev += npercell;
  }
  Kokkos::deep_copy(d_npercell, h_npercell);

  int ncands;
  auto d_cells2cands = offset_scan(d_npercell, ncands);
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

  for (int i = 0; i < nglocal; i++) {
    auto ncreate = h_npercell(i);
    lo = cells[i].lo;
    hi = cells[i].hi;

    for (int m = 0; m < ncreate; m++) {
      auto cand = h_cells2cands(i) + m;
      auto rn = random->uniform();
      int isp = 0;
      while (cummulative[isp] < rn) isp++;
      auto ispecies = species[isp];

      double x[3];
      x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
      x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
      x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
      if (dimension == 2) x[2] = 0.0;

      if (region && !region->match(x)) continue;
      if (speciesflag) {
        isp = species_variable(x) - 1;
        if (isp < 0 || isp >= nspecies) continue;
        ispecies = species[isp];
      }

      h_keep(cand) = 1;
      h_isp(cand) = isp;
      for (int d = 0; d < 3; ++d) h_x(cand, d) = x[d];

      if (tempflag) {
        tempscale = temperature_variable(x);
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
      if (velflag) {
        velocity_variable(x,vstream,vstream_variable);
        v[0] = vstream_variable[0] + vn*cos(theta1);
        v[1] = vstream_variable[1] + vr*cos(theta2);
        v[2] = vstream_variable[2] + vr*sin(theta2);
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

  Kokkos::parallel_for(nglocal, SPARTA_LAMBDA(int i) {
    auto ncreate = d_npercell(i);
    for (int m = 0; m < ncreate; m++) {
      auto cand = d_cells2cands(i) + m;
      if (!d_keep(cand)) continue;
      auto inew = d_cands2new(cand) + nlocal_before;
      auto id = d_id(cand);
      auto ispecies = d_species(d_isp(cand));
      double x[3],v[3];
      for (int d = 0; d < 3; ++d) x[d] = d_x(cand, d);
      for (int d = 0; d < 3; ++d) v[d] = d_v(cand, d);
      auto erot = d_erot(cand);
      auto evib = d_evib(cand);
      ParticleKokkos::add_particle_kokkos(d_particles,inew,id,ispecies,i,x,v,erot,evib);
    }
  });
  particleKK->modify(Device,PARTICLE_MASK);
  particleKK->sync(Host,PARTICLE_MASK);
  particleKK->nlocal += nnew;

  auto h_cands2new = Kokkos::create_mirror_view(d_cands2new);
  Kokkos::deep_copy(h_cands2new, d_cands2new);

  for (int i = 0; i < nglocal; i++) {
    auto ncreate = h_npercell(i);
    for (int m = 0; m < ncreate; m++) {
      auto cand = h_cells2cands(i) + m;
      if (!h_keep(cand)) continue;
      auto inew = h_cands2new(cand) + nlocal_before;
      if (nfix_update_custom)
        modify->update_custom(inew,temp_thermal,temp_rot,temp_vib,vstream);
    }
  }

  delete random;
}

