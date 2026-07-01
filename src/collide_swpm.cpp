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
#include "math_extra.h"
#include "math_eigen.h"
#include "math_eigen_impl.h"
#include "string.h"
#include "collide.h"
#include "particle.h"
#include "mixture.h"
#include "update.h"
#include "grid.h"
#include "comm.h"
#include "react.h"
#include "modify.h"
#include "fix.h"
#include "fix_ambipolar.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "memory.h"
#include "error.h"
#include <algorithm>

using namespace SPARTA_NS;

enum{ENERGY,HEAT,STRESS};   // particle reduction choices
enum{BINARY,WEIGHT,OCTREE}; // grouping choices

#define DELTADELETE 1024
#define BIG 1.0e20
#define SMALL 1.0e-16

/* ----------------------------------------------------------------------
   Stochastic weighted particle method algorithm
------------------------------------------------------------------------- */

template < int NEARCP, int GASTALLY >
void Collide::collisions_one_stochastic_weighting()
{
  int i,j,n,ip,np;
  int nattempt;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart,*lpart,*mpart;

  // loop over cells I own

  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::ChildCell *cells = grid->cells;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  double isw;
  double *stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 1) continue;

    volume = cinfo[icell].volume / cinfo[icell].weight;
    if (volume == 0.0) error->one(FLERR,"Collision cell volume is zero");

    // setup particle list for this cell

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(plist);
      memory->create(plist,npmax,"collide:plist");
    }

    // build particle list and find maximum particle weight
    // particle weights are relative to update->fnum

    ip = cinfo[icell].first;
    n = 0;
    max_stochastic_weight = 0.0;
    while (ip >= 0) {
      plist[n++] = ip;

      ipart = &particles[ip];
      isw = stochastic_weights[ip];
      max_stochastic_weight = MAX(max_stochastic_weight,isw);

      // error->one, not error->all: this code executes independently on each
      // rank, so a collective error here would hang the other MPI ranks

      if (isw != isw) error->one(FLERR,"Particle has NaN weight");
      if (isw <= 0.0) error->one(FLERR,"Particle has negative or zero weight");
      ip = next[ip];
    }

    // synchronize np with actual plist size
    np = n;

    // attempt = exact collision attempt count for all particles in cell
    // nattempt = rounded attempt with RN
    // if no attempts, continue to next grid cell

    if (np >= Ncmin && Ncmin > 0.0) pre_wtf = 0.0;
    else pre_wtf = 1.0;

    attempt = attempt_collision(icell,np,volume);
    nattempt = static_cast<int> (attempt);

    if (!nattempt) continue;
    nattempt_one += nattempt;

    for (int iattempt = 0; iattempt < nattempt; iattempt++) {

      i = np * random->uniform();
      j = np * random->uniform();
      while (i == j) j = np * random->uniform();
      ipart = &particles[plist[i]];
      jpart = &particles[plist[j]];

      if (!test_collision(icell,0,0,ipart,jpart)) continue;

      // split particles

      split(ipart,jpart,kpart,lpart);

      // add new particles to particle list

      if (kpart != NULL && lpart != NULL) {
        if (np+2 >= npmax) {
          npmax += DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
        }
        plist[np++] = particle->nlocal-2;
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
        stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];
      } else if (kpart != NULL || lpart != NULL) {
        if (np+1 >= npmax) {
          npmax += DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
        }
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
        stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];
      }

      // since ipart and jpart have same weight, do not need
      // ... to account for weight during collision itself
      // also the splits are all handled beforehand

      // check that particles have the same stochastic weight
      int i_index = ipart - particle->particles;
      int j_index = jpart - particle->particles;
      double isw = stochastic_weights[i_index];
      double jsw = stochastic_weights[j_index];
      if (isw != jsw)
        error->one(FLERR,"Particles must have same stochastic weight before collision in stochastic weighting");

      mpart = NULL; // dummy particle
      setup_collision(ipart,jpart);
      perform_collision(ipart,jpart,mpart);
      ncollide_one++;

    } // end attempt loop
  } // loop for cells

  // perform particle reduction if enabled
  if (reduceflag) group();

  return;
}

/* ----------------------------------------------------------------------
   Splits particles and generates two new particles (for SWPM)
------------------------------------------------------------------------- */

void Collide::split(Particle::OnePart *&ip, Particle::OnePart *&jp,
                    Particle::OnePart *&kp, Particle::OnePart *&lp)
{
  double xk[3],vk[3];
  double xl[3],vl[3];
  double erotk, erotl;
  int ks, ls;
  int kcell, lcell;

  // checks if particles properly deleted

  kp = NULL;
  lp = NULL;

  // weight transfer function is assumed to be
  // ... MIN(ip->stochastic_weights,jp->stochastic_weights)/(1 + pre_wtf * wtf)

  double *stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];
  int i_index = ip - particle->particles;
  int j_index = jp - particle->particles;
  double isw = stochastic_weights[i_index];
  double jsw = stochastic_weights[j_index];
  double Gwtf, ksw, lsw;

  if (isw <= 0.0 || jsw <= 0.0)
    error->one(FLERR,"Zero or negative weight before split");

  // particle ip has larger weight

  if(isw >= jsw) {
    Gwtf = jsw/(1.0+pre_wtf*wtf);
    ksw  = isw-Gwtf;
    lsw  = jsw-Gwtf;

    ks = ip->ispecies;
    ls = jp->ispecies;

    kcell = ip->icell;
    lcell = jp->icell;

    memcpy(xk,ip->x,3*sizeof(double));
    memcpy(vk,ip->v,3*sizeof(double));
    memcpy(xl,jp->x,3*sizeof(double));
    memcpy(vl,jp->v,3*sizeof(double));

    erotk = ip->erot;
    erotl = jp->erot;

  // particle jp has larger weight

  } else {
    Gwtf = isw/(1.0+pre_wtf*wtf);
    ksw  = jsw-Gwtf;
    lsw  = isw-Gwtf;

    ks = jp->ispecies;
    ls = ip->ispecies;

    kcell = jp->icell;
    lcell = ip->icell;

    memcpy(xk,jp->x,3*sizeof(double));
    memcpy(vk,jp->v,3*sizeof(double));
    memcpy(xl,ip->x,3*sizeof(double));
    memcpy(vl,ip->v,3*sizeof(double));

    erotk = jp->erot;
    erotl = ip->erot;
  }

  // update stochastic weights in custom array (relative to fnum)
  // ip and jp both keep the transferred weight Gwtf so the equal-weight
  // collision performed on them afterward is valid
  // recalculate indices after potential reallocation in add_particle

  i_index = ip - particle->particles;
  j_index = jp - particle->particles;
  stochastic_weights[i_index] = Gwtf;
  stochastic_weights[j_index] = Gwtf;

  // Gwtf should never be negative or zero

  if (Gwtf <= 0.0)
    error->one(FLERR,"Negative weight assigned after split");

  if (Gwtf > 0.0 && pre_wtf > 0.0)
    if (ksw <= 0.0 || lsw <= 0.0)
      error->one(FLERR,"Zero or negative weight after split");

  // gk is always the bigger of the two

  if(ksw > 0) {
    int id = MAXSMALLINT*random->uniform();
    Particle::OnePart *particles = particle->particles;
    int reallocflag = particle->add_particle(id,ks,kcell,xk,vk,erotk,0.0);
    if (reallocflag) {
      ip = particle->particles + (ip - particles);
      jp = particle->particles + (jp - particles);
    }
    kp = &particle->particles[particle->nlocal-1];
    // add_particle may have grown (reallocated) the custom arrays, so
    // re-acquire the stochastic weight array before writing into it
    stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];
    stochastic_weights[particle->nlocal-1] = ksw;
  }

  if (kp && ksw <= 0.0)
    error->one(FLERR,"New particle [k] has bad weight");

  // there should never be case where you add particle "l" if
  // ... you did not add particle "k"

  if(lsw > 0) {
    if(ksw <= 0) error->one(FLERR,"Bad addition to particle list");
    int id = MAXSMALLINT*random->uniform();
    Particle::OnePart *particles = particle->particles;
    int reallocflag = particle->add_particle(id,ls,lcell,xl,vl,erotl,0.0);
    if (reallocflag) {
      ip = particle->particles + (ip - particles);
      jp = particle->particles + (jp - particles);
      kp = particle->particles + (kp - particles);
    }
    lp = &particle->particles[particle->nlocal-1];
    // re-acquire after possible reallocation in add_particle (see above)
    stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];
    stochastic_weights[particle->nlocal-1] = lsw;
  }

  if (lp && lsw <= 0.0)
    error->one(FLERR,"New particle [l] has bad weight");
}

/* ----------------------------------------------------------------------
   Reorder plist depending on grouping strategy used and prepare for
   grouping and particle reduction
------------------------------------------------------------------------- */

void Collide::group()
{
  int n,nold,np,ip;
  double isw;

  Particle::OnePart *ipart;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::ChildCell *cells = grid->cells;
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  double *stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];

  double swmean, swvar, swstd;
  double d1, d2;
  double lLim, uLim;
  int np_med,np_small;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;

    if (np <= Ncmax) continue;

    // create particle list
    // negative weights ones are the tiny ones remove before

    ip = cinfo[icell].first;
    n = 0;
    while (ip >= 0) {
      ipart = &particles[ip];
      isw = stochastic_weights[ip];
      if(isw > 0) plist[n++] = ip;
      ip = next[ip];
    }

    int group_size_buffer = 0;

    while (n > Ncmax) {
      nold = n;

      // route to appropriate grouping function based on group_type

      if (group_type == BINARY) {
        group_bt(0,n,group_size_buffer);
      } else if (group_type == OCTREE) {
        group_octree(0,n,group_size_buffer);
      } else if (group_type == WEIGHT) {

        // find mean / standard deviation of weight
        // plist[0..n-1] holds only positive-weight particles (built above)

        swmean = swvar = 0.0;
        for (int i = 0; i < n; i++) {
          isw = stochastic_weights[plist[i]];

          // Incremental variance

          d1 = isw - swmean;
          swmean += d1/(i+1.0);
          swvar += i/(i+1.0)*d1*d1;
        }
        swstd = sqrt(swvar/n);

        // weight limits to separate particles
        // uLim must be inclusive below so the degenerate case of all-equal
        // weights (swstd = 0, lLim = uLim = swmean) still classifies every
        // particle as medium; otherwise no particle is ever reduced

        lLim = MAX(swmean-1.25*swstd,0);
        uLim = swmean+2.0*swstd;

        // reorder plist: small weighted particles first

        np_small = 0;
        for (int i = 0; i < n; i++) {
          isw = stochastic_weights[plist[i]];
          if (isw < lLim) {
            std::swap(plist[i],plist[np_small]);
            np_small++;
          }
        }

        // then medium weighted particles, starting at index np_small
        // particles above uLim stay at the end and are not reduced

        np_med = np_small;
        for (int i = np_small; i < n; i++) {
          isw = stochastic_weights[plist[i]];
          if (isw >= lLim && isw <= uLim) {
            std::swap(plist[np_med],plist[i]);
            np_med++;
          }
        }

        // ignore the very large weighted particles
        // group and reduce using binary tree strategy
        // (WEIGHT preprocessing followed by BINARY grouping)

        group_bt(0, np_small, group_size_buffer);
        group_bt(np_small, np_med, group_size_buffer);

      }

      // recreate particle list after reduction

      ip = cinfo[icell].first;
      n = 0;
      while (ip >= 0) {
        ipart = &particles[ip];
        isw = stochastic_weights[ip];
        if(isw > 0) plist[n++] = ip;
        ip = next[ip];
      }

      // if no particles reduced, increase group size
      
      if (n == nold) group_size_buffer += 2;
      if (group_size_buffer > n) break;

    } // while loop for n > ncmax
  }// loop for cells

  return;
}

/* ----------------------------------------------------------------------
   Octree grouping strategy (placeholder for future implementation)
------------------------------------------------------------------------- */

void Collide::group_octree(int istart, int iend, int group_size_buffer)
{
  // TODO: implement octree grouping strategy
  error->all(FLERR,"Octree grouping strategy not yet implemented");
}

/* ----------------------------------------------------------------------
   Recursivley divides particles using the binary tree strategy
------------------------------------------------------------------------- */
void Collide::group_bt(int istart, int iend, int group_size_buffer)
{
  Particle::OnePart *ipart;
  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  double *stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];

  // ignore groups which have too few particles

  // loops don't include iend
  int np = iend-istart;
  if (np <= Ngmin) return;

  // accumulate group totals: weight, mass, momentum, rotational energy
  // first pass forms the mean velocity used to build central moments below

  double gsum, msum, Erot;
  double mV[3];
  gsum = msum = Erot = 0.0;
  for (int i = 0; i < 3; i++) mV[i] = 0.0;

  int ispecies;
  double mass, psw, pmsw;
  for (int p = istart; p < iend; p++) {
    ipart = &particles[plist[p]];
    ispecies = ipart->ispecies;
    mass = species[ispecies].mass;

    psw = stochastic_weights[plist[p]];
    pmsw = psw * mass;
    gsum += psw;
    msum += pmsw;
    Erot += psw*ipart->erot;
    for (int i = 0; i < 3; i++) mV[i] += pmsw*ipart->v[i];
  }

  // mean (mass-weighted) velocity

  double V[3];
  for (int i = 0; i < 3; i++) V[i] = mV[i]/msum;

  // second pass: accumulate the stress tensor and heat flux directly from
  // velocity deviations (v - V).  computing these central moments with
  // nested loops avoids the catastrophic cancellation of the raw-moment
  // form (e.g. mVV - mV*mV/msum), which could yield negative temperatures.

  double pij[3][3], q[3];
  for (int i = 0; i < 3; i++) {
    q[i] = 0.0;
    for (int j = 0; j < 3; j++) pij[i][j] = 0.0;
  }

  for (int p = istart; p < iend; p++) {
    ipart = &particles[plist[p]];
    ispecies = ipart->ispecies;
    mass = species[ispecies].mass;

    psw = stochastic_weights[plist[p]];
    pmsw = psw * mass;

    double dv[3];
    for (int i = 0; i < 3; i++) dv[i] = ipart->v[i] - V[i];
    double dvsq = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];

    for (int i = 0; i < 3; i++) {
      // heat flux: q_i = 1/2 sum w*m*(v_i - V_i)*|v - V|^2
      q[i] += 0.5*pmsw*dv[i]*dvsq;
      // stress tensor: P_ij = sum w*m*(v_i - V_i)*(v_j - V_j)
      for (int j = 0; j < 3; j++) pij[i][j] += pmsw*dv[i]*dv[j];
    }
  }

  // if group is small enough, merge the particles

  if (np <= Ngmax+group_size_buffer) {

    // temperature (becomes the velocity variance after the scaling below)
    double T = (pij[0][0] + pij[1][1] + pij[2][2])/
      (3.0 * gsum * update->boltz);

    // q (heat flux) and pij (stress tensor) were accumulated as central
    // moments in the second pass above

    // scale values to be consistent with definitions in
    // .. stochastic numerics book

    T *= update->boltz/mass;
    for (int i = 0; i < 3; i++) {
      q[i] /= mass;
      for(int j = 0; j < 3; j++) pij[i][j] /= mass;
    }

    // reduce based on type
    if (reduction_type == ENERGY) {
      reduce_energy(istart, iend, gsum, V, T, Erot);
    } else if (reduction_type == HEAT) {
      reduce_heat(istart, iend, gsum, V, T, Erot, q);
    } else if (reduction_type == STRESS) {
      reduce_stress(istart, iend, gsum, V, T, Erot, q, pij);
    }

  // group still too large so divide further

  } else {

    // Compute covariance matrix

    double Rij[3][3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        Rij[i][j] = pij[i][j]/gsum;

    // Find eigenpairs

    double eval[3], evec[3][3];
    int ierror = MathEigen::jacobi3(Rij,eval,evec);

    // Find largest eigenpair
    // maxeval initialized to -1 so maxevec is always set, even in the
    // degenerate case of all-zero eigenvalues (all velocities identical)

    double maxeval;
    double maxevec[3]; // normal of splitting plane

    maxeval = -1.0;
    for (int i = 0; i < 3; i++) {
      if (std::abs(eval[i]) > maxeval) {
        maxeval = std::abs(eval[i]);
        for (int j = 0; j < 3; j++)
          maxevec[j] = evec[j][i];
      }
    }

    // Sort plist based on signed distance from center
    // signed distance = dot(particle_velocity, maxevec) - center

    double center = V[0]*maxevec[0] + V[1]*maxevec[1] + V[2]*maxevec[2];
    
    // Sort plist[istart..iend-1] based on signed distance
    std::sort(plist + istart, plist + iend, 
      [&](int i, int j) {
        Particle::OnePart *pi = &particles[i];
        Particle::OnePart *pj = &particles[j];
        double dist_i = MathExtra::dot3(pi->v, maxevec) - center;
        double dist_j = MathExtra::dot3(pj->v, maxevec) - center;
        return dist_i < dist_j;
      });

    // Find split point: particles with negative signed distance
    int cp = 0; // count of particles with negative signed distance
    for (int i = istart; i < iend; i++) {
      ipart = &particles[plist[i]];
      double dist = MathExtra::dot3(ipart->v, maxevec) - center;
      if (dist < 0.0) cp++;
      else break; // since sorted, can break early
    }

    // if every particle landed on one side of the plane (e.g. all velocities
    // identical so all distances are zero), split the group in half instead;
    // recursing on the unchanged range would never terminate

    if (cp == 0 || cp == np) cp = np/2;

    if(cp > Ngmin) group_bt(istart,istart+cp,group_size_buffer);
    if(np-cp > Ngmin) group_bt(istart+cp,iend,group_size_buffer);
  }

  return;
}

/* ----------------------------------------------------------------------
   explicit instantiations of the SWPM collision template for the
   NEARCP/GASTALLY combinations dispatched from collide.cpp
   (the template is defined here, in a separate translation unit)
------------------------------------------------------------------------- */

namespace SPARTA_NS {
template void Collide::collisions_one_stochastic_weighting<0,0>();
template void Collide::collisions_one_stochastic_weighting<0,1>();
template void Collide::collisions_one_stochastic_weighting<1,0>();
template void Collide::collisions_one_stochastic_weighting<1,1>();
}
