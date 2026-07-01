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
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

#define DELTADELETE 1024
#define SMALL 1.0e-16

/* ----------------------------------------------------------------------
   flag plist entry idx for deletion: mark its stochastic weight and queue it
------------------------------------------------------------------------- */
void Collide::reduce_delete(int idx, double *stochastic_weights)
{
  if (ndelete == maxdelete) {
    maxdelete += DELTADELETE;
    memory->grow(dellist,maxdelete,"collide:dellist");
  }
  stochastic_weights[idx] = -1.0;
  dellist[ndelete++] = idx;
}

/* ----------------------------------------------------------------------
   Merge particles using energy scheme
------------------------------------------------------------------------- */
void Collide::reduce_energy(int istart, int iend,
                            double rho, double *V, double T, double Erot)
{

  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;
  double *stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];

  int np = iend-istart;
  int ip, jp;
  ip = np * random->uniform() + istart;
  jp = np * random->uniform() + istart;
  while (ip == jp) jp = np * random->uniform() + istart;

  ipart = &particles[plist[ip]];
  jpart = &particles[plist[jp]];

  // find direction of velocity wrt CoM frame

  double theta = 2.0 * MY_PI * random->uniform();
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

  // set reduced particle stochastic weights (relative to fnum)

  stochastic_weights[plist[ip]] = rho*0.5;
  stochastic_weights[plist[jp]] = rho*0.5;

  // delete other particles

  for (int i = istart; i < iend; i++) {
    if (i == ip || i == jp) continue;
    reduce_delete(plist[i],stochastic_weights);
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using heat flux scheme
------------------------------------------------------------------------- */
void Collide::reduce_heat(int istart, int iend,
                          double rho, double *V, double T, double Erot,
                          double *q)
{

  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;
  double *stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];

  int np = iend-istart;
  int ip, jp;
  ip = np * random->uniform() + istart;
  jp = np * random->uniform() + istart;
  while (ip == jp) jp = np * random->uniform() + istart;

  ipart = &particles[plist[ip]];
  jpart = &particles[plist[jp]];

  // precompute
  // if the group has (near) zero thermal energy, qge below would be 0/0;
  // fall back to the energy scheme (theta = 1, equal weights), which
  // handles T = 0 correctly

  double sqT = sqrt(3.0*T);
  double qmag = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
  double qge;
  if (sqT < SMALL) qge = 0.0;
  else qge = qmag / (rho * pow(sqT,3.0));
  double itheta = qge + sqrt(1.0 + qge*qge);
  double alpha = sqT*itheta;
  double beta = sqT/itheta;

  // find direction of velocity wrt CoM frame
  // for zero heat flux any direction conserves the tallied moments, but the
  // vector must have unit length or the thermal energy sqT placed along it
  // is not conserved; sample uniformly on the unit sphere

  double uvec[3];
  if (qmag < SMALL) {
    double theta = 2.0 * MY_PI * random->uniform();
    double phi = acos(1.0 - 2.0 * random->uniform());
    uvec[0] = sin(phi) * cos(theta);
    uvec[1] = sin(phi) * sin(theta);
    uvec[2] = cos(phi);
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

  stochastic_weights[plist[ip]] = isw;
  stochastic_weights[plist[jp]] = jsw;

  // delete other particles
  for (int i = istart; i < iend; i++) {
    if (i == ip || i == jp) continue;
    reduce_delete(plist[i],stochastic_weights);
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using stress scheme
------------------------------------------------------------------------- */
void Collide::reduce_stress(int istart, int iend,
                            double rho, double *V, double T, double Erot,
                            double *q, double pij[3][3])
{

  // find eigenpairs of stress tensor

  double eval[3], evec[3][3];
  int ierror = MathEigen::jacobi3(pij,eval,evec);

  // find number of non-zero eigenvalues
  // jacobi3 returns eigenvectors as columns (evec[d][i] = component d of
  // eigenvector i), so compact columns, not rows

  int nK = 0;
  for (int i = 0; i < 3; i++) {
    if (fabs(eval[i]) >= SMALL && eval[i] > 0) {
      eval[nK] = eval[i];
      for (int d = 0; d < 3; d++) evec[d][nK] = evec[d][i];
      nK++;
    }
  }

  // degenerate group: all eigenvalues (near) zero, i.e. no thermal spread.
  // fall back to a 2-particle energy-style merge; deleting every particle in
  // the group (empty loop below) would destroy its mass, momentum and energy

  if (nK == 0) {
    reduce_energy(istart,iend,rho,V,T,Erot);
    return;
  }

  // stress reduction can emit up to 2*nK = 6 survivor particles.  groups are
  // guaranteed to have more than 6 particles by the collide_modify reduce
  // command (which requires Ngmin >= 6 for stress), so 2*nK never exceeds the
  // group size below.

  // iterate through nK pairs

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;
  double *stochastic_weights = particle->edvec[particle->ewhich[index_stochastic_weight]];

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

    stochastic_weights[plist[2*iK+istart]] = isw;
    stochastic_weights[plist[2*iK+1+istart]] = jsw;

  } // end nK

  // delete other particles
  // the 2*nK reduced survivors occupy plist[istart .. istart+2*nK-1]
  for (int i = istart + 2*nK; i < iend; i++)
    reduce_delete(plist[i],stochastic_weights);

  return;
}

