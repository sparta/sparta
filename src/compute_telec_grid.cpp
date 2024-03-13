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

#include "compute_telec_grid.h"
#include "error.h"
#include "grid.h"
#include "memory.h"
#include "mixture.h"
#include "particle.h"
#include "update.h"

using namespace SPARTA_NS;

ComputeTelecGrid::ComputeTelecGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute telec/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute telec/grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute telec/grid mixture ID does not exist");

  int iarg = 4;
  while (iarg < narg) {
    error->all(FLERR,"Illegal compute telec/grid command");
  }

  ngroup = particle->mixture[imix]->ngroup;
  mixspecies = particle->mixture[imix]->nspecies;
  nspecies = particle->nspecies;

  ntally = 2*mixspecies;

  per_grid_flag = 1;
  size_per_grid_cols = 0;
  post_process_grid_flag = 1;

  particle->mixture[imix]->init();
  int *groupsize = particle->mixture[imix]->groupsize;
  int **groupspecies = particle->mixture[imix]->groupspecies;

  int nmax = 0;
  for (int i = 0; i < ngroup; i++)
    nmax = MAX(nmax,groupsize[i]);

  nmap = new int[ngroup];
  memory->create(map,ngroup,2*nmax,"telec/grid:map");
  tspecies = new double[nspecies];
  s2t = new int[nspecies];
  t2s = new int[ntally];

  for (int isp = 0; isp < nspecies; isp++) s2t[isp] = -1;

  int itally = 0;
  for (int igroup = 0; igroup < ngroup; igroup++) {
    nmap[igroup] = 2*groupsize[igroup];
    for (int n = 0; n < groupsize[igroup]; n++) {
      map[igroup][2*n] = itally;
      map[igroup][2*n+1] = itally+1;
      s2t[groupspecies[igroup][n]] = itally;
      t2s[itally] = groupspecies[igroup][n];
      t2s[itally+1] = groupspecies[igroup][n];
      itally += 2;
    }
  }

  nglocal = 0;
  vector_grid = NULL;
  tally = NULL;
}

ComputeTelecGrid::~ComputeTelecGrid()
{
  if (copymode) return;

  delete [] nmap;
  memory->destroy(map);

  delete [] tspecies;
  delete [] s2t;
  delete [] t2s;

  memory->destroy(vector_grid);
  memory->destroy(tally);
}

void ComputeTelecGrid::init() {
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute telec/grid "
               "mixture has changed");
  if (mixspecies != particle->mixture[imix]->nspecies)
    error->all(FLERR,"Number of species in compute telec/grid "
               "mixture has changed");
  if (nspecies != particle->nspecies)
    error->all(FLERR,"Number of total species in compute telec/grid "
               "has changed");

  reallocate();
}

void ComputeTelecGrid::compute_per_grid() {
  invoked_per_grid = update->ntimestep;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;

  int *s2g = particle->mixture[imix]->species2group;
  double boltz = update->boltz;
  int nlocal = particle->nlocal;

  int i,j,ispecies,igroup,icell,imode,nmode;

  // zero all accumulators - could do this with memset()

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntally; j++)
      tally[i][j] = 0.0;

  // loop over all particles, skip species not in mixture group
  // tally vib eng and count for species
  for (i = 0; i < nlocal; i++) {
      ispecies = particles[i].ispecies;
      igroup = s2g[ispecies];
      if (igroup < 0) continue;
      icell = particles[i].icell;
      if (!(cinfo[icell].mask & groupbit)) continue;

      j = s2t[ispecies];
      char data_name[] = "eelec";
      int eelec_index = particle->find_custom(data_name);
      double *eelecs = NULL;
      if (eelec_index >= 0) {
        eelecs = particle->edvec[particle->ewhich[eelec_index]];
        tally[icell][j] += eelecs[i];
        tally[icell][j+1] += 1.0;
      }
    }
}

double ComputeTelecGrid::elec_energy(int isp, double temp_elec) {
    Particle::Species species = particle->species[isp];

    double* state_probabilities = particle->electronic_distribution_func(isp, temp_elec);

    double total_energy = 0.0;
    for (int i = 0; i < species.elecdat->nelecstate; ++i) {
      total_energy += state_probabilities[i]*species.elecdat->states[i].temp*update->boltz;
    }

    return total_energy;
}

void ComputeTelecGrid::post_process_grid(int index, int nsample,
  double** etally, int* emap,
  double* vec, int nstride)
{
  int lo = 0;
  int hi = nglocal;
  int k = 0;

  if (!etally) {
    nsample = 1;
    etally = tally;
    emap = map[index];
    vec = vector_grid;
    nstride = 1;
  }

  Particle::Species *species = particle->species;

  int eelec,count,ispecies,isp;
  double first_elec_eng, t_elec, degen0, degen1, numer, denom;
  double boltz = update->boltz;

  int nsp = nmap[index] / 2;

  for (int icell = lo; icell < hi; icell++) {
    eelec = emap[0];
    count = eelec+1;
    for (isp = 0; isp < nsp; isp++) {
      ispecies = t2s[eelec-emap[0]];
      if (species[ispecies].elecdat == NULL ||
          etally[icell][eelec] == 0.0) {
        tspecies[isp] = 0.0;
        eelec += 2;
        count = eelec+1;
        continue;
      }
      // We calculate a first guess at the temp assuming
      // all the electronic energy is stored in the first
      // excited state
      first_elec_eng = species[ispecies].elecdat->states[1].temp*boltz;
      degen0 = species[ispecies].elecdat->states[0].degen;
      degen1 = species[ispecies].elecdat->states[1].degen;
      t_elec = first_elec_eng / (boltz*(
          - log( etally[icell][eelec]*degen0 /
              (etally[icell][count]*first_elec_eng*degen1)
             ))
        );

      // Bisection method to find T accurate to 1%
      double target_energy_per_part = etally[icell][eelec]/etally[icell][count];

      // Find initial bounds based on our first guess
      double T_low = 0.9*t_elec;
      while (elec_energy(isp, T_low) > target_energy_per_part) {
        T_low /= 2.0;
      }

      double T_high = t_elec*1.1;
      while (elec_energy(isp, T_high) < target_energy_per_part) {
        T_high *= 2.0;
      }

      // Bisect
      double T_mid = t_elec;
      double e_mid = elec_energy( isp, T_mid );
      while ((T_high - T_low) > 0.01) {
        if (e_mid > target_energy_per_part) {
          T_high = T_mid;
        } else {
          T_low = T_mid;
        }
        T_mid = (T_high - T_low)/2.0 + T_low;
        e_mid = elec_energy(isp, T_mid);
      }

      tspecies[isp] = T_mid;
      eelec += 2;
      count = eelec+1;
    }

    numer = denom = 0.0;
    count = emap[1];
    for (isp = 0; isp < nsp; isp++) {
      numer += tspecies[isp]*etally[icell][count];
      denom += etally[icell][count];
      count += 2;
    }
    vec[k] = numer/denom;
    k += nstride;
  }
}

int ComputeTelecGrid::query_tally_grid(int index, double**& array, int*& cols) {
  index--;
  array = tally;
  cols = map[index];
  return nmap[index];
}

void ComputeTelecGrid::reallocate() {
  if (grid->nlocal == nglocal) return;

  memory->destroy(vector_grid);
  memory->destroy(tally);
  nglocal = grid->nlocal;
  memory->create(vector_grid,nglocal,"telec/grid:vector_grid");
  memory->create(tally,nglocal,ntally,"telec/grid:tally");

}

bigint ComputeTelecGrid::memory_usage() {
  bigint bytes;
  bytes = nglocal * sizeof(double);
  bytes = ntally*nglocal * sizeof(double);
  return bytes;
}
