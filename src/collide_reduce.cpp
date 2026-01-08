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

using namespace SPARTA_NS;

#define DELTADELETE 1024
#define SMALL 1.0e-16

/* ----------------------------------------------------------------------
   Merge particles using energy scheme
------------------------------------------------------------------------- */
void Collide::reduce(int istart, int iend, 
                     double rho, double *V, double T, double Erot)
{

  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  int np = iend-istart;
  int ip, jp;
  ip = np * random->uniform() + istart;
  jp = np * random->uniform() + istart;
  while (ip == jp) jp = np * random->uniform() + istart;

  ipart = &particles[plist[ip]];
  jpart = &particles[plist[jp]];

  // find direction of velocity wrt CoM frame

  double theta = 2.0 * 3.14159 * random->uniform();
  double phi = acos(1.0 - 2.0 * random->uniform());
  double uvec[3];
  uvec[0] = sin(phi) * cos(theta);
  uvec[1] = sin(phi) * sin(theta);
  uvec[2] = cos(phi);

  // set reduced particle velocities

  double sqT = sqrt(3.0*T);
  for (int d = 0; d < 3; d++) {
    ipart->v[d] = V[d] + sqT*uvec[d];
    jpart->v[d] = V[d] - sqT*uvec[d];
  }

  // set reduced particle rotational energies

  ipart->erot = Erot/(rho*0.5)*0.5;
  jpart->erot = Erot/(rho*0.5)*0.5;

  // set reduced particle weights

  ipart->weight = rho*0.5;
  jpart->weight = rho*0.5;

  // delete other particles

  for (int i = istart; i < iend; i++) {
    if (i == ip || i == jp) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[plist[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = plist[i];
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using heat flux scheme
------------------------------------------------------------------------- */
void Collide::reduce(int istart, int iend,
                     double rho, double *V, double T, double Erot, double *q)
{

  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  int np = iend-istart;
  int ip, jp;
  ip = np * random->uniform() + istart;
  jp = np * random->uniform() + istart;
  while (ip == jp) jp = np * random->uniform() + istart;

  ipart = &particles[plist[ip]];
  jpart = &particles[plist[jp]];

  // precompute

  double sqT = sqrt(3.0*T);
  double qmag = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
  double qge = qmag / (rho * pow(sqT,3.0));
  double itheta = qge + sqrt(1.0 + qge*qge);
  double alpha = sqT*itheta;
  double beta = sqT/itheta;

  // find direction of velocity wrt CoM frame

  double uvec[3];
  if (qmag < SMALL) {
    for (int d = 0; d < 3; d++) {
      double A = sqrt(-log(random->uniform()));
      double phi = 6.283185308 * random->uniform();
      if (random->uniform() < 0.5) uvec[d] = A * cos(phi);
      else uvec[d] = A * sin(phi);
    }
  } else 
    for (int d = 0; d < 3; d++) uvec[d] = q[d]/qmag;

  // set reduced particle velocities

  for (int d = 0; d < 3; d++) {
    ipart->v[d] = V[d] + alpha*uvec[d];
    jpart->v[d] = V[d] - beta*uvec[d];
  }

  // set reduced particle weights

  double isw = rho/(1.0+itheta*itheta);
  double jsw = rho - isw;

  // set reduced particle rotational energies

  ipart->erot = Erot/isw*0.5;
  jpart->erot = Erot/jsw*0.5;

  ipart->weight = isw;
  jpart->weight = jsw;

  // delete other particles
  for (int i = istart; i < iend; i++) {
    if (i == ip || i == jp) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[plist[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = plist[i];
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using stress scheme
------------------------------------------------------------------------- */
void Collide::reduce(int istart, int iend,
                     double rho, double *V, double T, double Erot,
                     double *q, double pij[3][3])
{

  // find eigenpairs of stress tensor

  double eval[3], evec[3][3];
  int ierror = MathEigen::jacobi3(pij,eval,evec);

  // find number of non-zero eigenvalues

  int nK = 0;
  for (int i = 0; i < 3; i++) {
    if (fabs(eval[i]) >= SMALL && eval[i] > 0) {
      eval[nK] = eval[i];
      for (int d = 0; d < 3; d++) evec[nK][d] = evec[i][d];
      nK++;
    }
  }

  // iterate through nK pairs

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  double qli, itheta;
  double isw, jsw;
  double uvec[3];

  for (int iK = 0; iK < nK; iK++) {

    // reduced particles chosen as first two

    ipart = &particles[plist[2*iK+istart]];
    jpart = &particles[plist[2*iK+1+istart]];

    qli = evec[0][iK]*q[0] + evec[1][iK]*q[1] + evec[2][iK]*q[2];
    if (qli < 0)
      for (int d = 0; d < 3; d++) evec[d][iK] *= -1.0;
    qli = fabs(qli);

    itheta = sqrt(rho) * qli / (sqrt(nK) * pow(eval[iK],1.5))
      + sqrt(1.0 + (rho*qli*qli)/(nK*pow(eval[iK],3.0)));

    // set reduced particle velocities

    for (int d = 0; d < 3; d++) {
      ipart->v[d] = V[d] + itheta*sqrt(nK*eval[iK]/rho)*evec[d][iK];
      jpart->v[d] = V[d] - 1.0/itheta*sqrt(nK*eval[iK]/rho)*evec[d][iK];
    }

    // set reduced particle weights

    isw = rho/(nK*(1.0+itheta*itheta));
    jsw = rho/nK - isw;

    // set reduced particle rotational energies

    ipart->erot = Erot/isw*0.5/nK;
    jpart->erot = Erot/jsw*0.5/nK;

    ipart->weight = isw;
    jpart->weight = jsw;

  } // end nK
  
  // delete other particles
  for (int i = istart; i < iend; i++) {
    if (i < 2*nK) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[plist[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = plist[i];
  }

  return;
}

