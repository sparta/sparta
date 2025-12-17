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
      isw = stochastic_weights[ip] * ipart->weight;
      max_stochastic_weight = MAX(max_stochastic_weight,isw);

      if (isw != isw) error->all(FLERR,"Particle has NaN weight");
      if (isw <= 0.0) error->all(FLERR,"Particle has negative or zero weight");
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
      } else if (kpart != NULL || lpart != NULL) {
        if (np+1 >= npmax) {
          npmax += DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
        }
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
      }

      // since ipart and jpart have same weight, do not need
      // ... to account for weight during collision itself
      // also the splits are all handled beforehand

      // check that particles have the same weight (weight * stochastic_weights)
      int i_index = ipart - particle->particles;
      int j_index = jpart - particle->particles;
      double isw = stochastic_weights[i_index] * ipart->weight;
      double jsw = stochastic_weights[j_index] * jpart->weight;
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
  double isw = stochastic_weights[i_index] * ip->weight;
  double jsw = stochastic_weights[j_index] * jp->weight;
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

  // update weights in custom array only
  // recalculate indices after potential reallocation in add_particle

  i_index = ip - particle->particles;
  j_index = jp - particle->particles;
  stochastic_weights[i_index] = Gwtf / ip->weight;
  stochastic_weights[j_index] = Gwtf / jp->weight;

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
    stochastic_weights[particle->nlocal-1] = ksw / update->fnum;
  }

  if (kp) {
    if (kp->weight <= 0.0) {
      printf("ksw: %2.3e; kp->weight: %2.3e\n", ksw,kp->weight);
      error->one(FLERR,"New particle [k] has bad weight");
    }
  }

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
    stochastic_weights[particle->nlocal-1] = lsw / update->fnum;
  }

  if (lp) {
    if (lp->weight <= 0.0) {
      printf("lsw: %2.3e; lp->weight: %2.3e\n", lsw,lp->weight);
      error->one(FLERR,"New particle [l] has bad weight");
    }
  }
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
      isw = stochastic_weights[ip] * ipart->weight;
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

        n = 0.0;
        swmean = swvar = 0.0;
        for (int i = 0; i < np; i++) {
          ipart = &particles[plist[i]];
          isw = stochastic_weights[plist[i]] * ipart->weight;

          // Incremental variance

          if(isw > 0) {
            n++;
            d1 = isw - swmean;
            swmean += (d1/n);
            swvar += (n-1.0)/n*d1*d1;
          }
          ip = next[ip];
        }
        swstd = sqrt(swvar/n);

        // weight limits to separate particles

        lLim = MAX(swmean-1.25*swstd,0);
        uLim = swmean+2.0*swstd;

        // recreate particle list and omit large weighted particles

        // first place all small weighted particles to front
        int np_small = 0; // index of center particle
        ip = cinfo[icell].first;
        for (int i = 0; i < np; i++) {
          ipart = &particles[plist[i]];
          isw = stochastic_weights[plist[i]] * ipart->weight;
          if(isw > 0 && isw < lLim) {
            std::swap(plist[ip],plist[np_small]);
            np_small++;
          }
          ip = next[ip];
        }

        // then place all medium weighted particles starting from
        // .., last index (np_small)

        np_med = np_small;
        for (int i = np_small; i < np; i++) {
          ipart = &particles[plist[i]];
          isw = stochastic_weights[plist[i]] * particles[plist[i]].weight;
          if (isw >= lLim && isw < uLim) {
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
        isw = stochastic_weights[ip] * particles[ip].weight;
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

  // ignore groups which have too few particles

  // loops don't include iend
  int np = iend-istart;
  if (np <= Ngmin) return;

  // zero out all values
  double gsum, msum;
  double mV[3];
  double pij[3][3];
  gsum = msum = 0.0;
  for (int i = 0; i < 3; i++) {
    mV[i] = 0.0;
    for (int j = 0; j < 3; j++) pij[i][j] = 0.0;
  }

  // find total weight, mass, momentum, and rotational energy

  int ispecies;
	double mass, psw, pmsw, vp[3];
  double Erot;
  for (int p = istart; p < iend; p++) {
    ipart = &particles[plist[p]];
    ispecies = ipart->ispecies;
    mass = species[ispecies].mass;

    psw = ipart->weight;
    pmsw = psw * mass;
    memcpy(vp, ipart->v, 3*sizeof(double));
   	gsum += psw;
    msum += pmsw;
    Erot += psw*ipart->erot;
    for (int i = 0; i < 3; i++) mV[i] += (pmsw*vp[i]);
  }

  // mean velocity
	double V[3];
  for (int i = 0; i < 3; i++) V[i] = mV[i]/msum;

  // stress tensor

  for (int p = istart; p < iend; p++) {
    ipart = &particles[plist[p]];
    ispecies = ipart->ispecies;
    mass = species[ispecies].mass;

    psw = ipart->weight;
    pmsw = psw * mass;
    memcpy(vp, ipart->v, 3*sizeof(double));
   	gsum += psw;
    msum += pmsw;
    Erot += psw*ipart->erot;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        pij[i][j] += (pmsw*(vp[i]-V[i])*(vp[j]-V[j]));
  }

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      pij[i][j] /= msum;

  // if group is small enough, merge the particles

  if (np <= Ngmax+group_size_buffer) {

    // temperature
    double T = (pij[0][0] + pij[1][1] + pij[2][2])/
      (3.0 * gsum * update->boltz);

    // heat flux
    double Vsq = V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
    double h,h1,h2,q[3];
    int i1,i2;
    for (int i = 0; i < 3; i++) {
      if (i == 0) {
        i1 = 1;
        i2 = 2;
      } else if (i == 1) {
        i1 = 2;
        i2 = 0;
      } else {
        i1 = 0;
        i2 = 1;
      }

      h  = mVVV[i][i] - 3.0*mV[i]*mVV[i][i]/msum +
           2.0*mV[i]*mV[i]*mV[i]/msum/msum;
      h1 = mVVV[i][i1] - 2.0*mVV[i][i1]*mV[i1]/msum -
           mV[i]*mVV[i1][i1]/msum + 2.0*mV[i]*mV[i1]*mV[i1]/msum/msum;
      h2 = mVVV[i][i2] - 2.0*mVV[i][i2]*mV[i2]/msum -
           mV[i]*mVV[i2][i2]/msum + 2.0*mV[i]*mV[i2]*mV[i2]/msum/msum;
      q[i] = (h + h1 + h2) * 0.5;
    }

    // scale values to be consistent with definitions in
    // .. stochastic numerics book

    T *= update->boltz/mass;
    for (int i = 0; i < 3; i++) {
      q[i] /= mass;
      for(int j = 0; j < 3; j++) pij[i][j] /= mass;
    }

    // reduce based on type
    if (reduction_type == ENERGY) {
      reduce(istart, iend, gsum, V, T, Erot);
    } else if (reduction_type == HEAT) {
      reduce(istart, iend, gsum, V, T, Erot, q);
    } else if (reduction_type == STRESS) {
      reduce(istart, iend, gsum, V, T, Erot, q, pij);
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

    double maxeval;
    double maxevec[3]; // normal of splitting plane

    maxeval = 0;
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

    if(cp > Ngmin) group_bt(istart,istart+cp,group_size_buffer);
    if(np-cp > Ngmin) group_bt(istart+cp,iend,group_size_buffer);
  }

  return;
}
