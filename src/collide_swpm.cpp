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

enum{ENERGY,HEAT,STRESS};   // particle reduction choices
enum{BINARY,WEIGHT}; // grouping choices

#define DELTADELETE 1024

#define BIG 1.0e20
#define SMALL 1.0e-16

/* ----------------------------------------------------------------------
   Stochastic weighted algorithm
------------------------------------------------------------------------- */

void Collide::collisions_one_sw()
{
  int i,j,n,ip,np,newp;
  int nattempt;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart,*lpart,*mpart;

  // loop over cells I own

  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::ChildCell *cells = grid->cells;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  double isw;

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
      memory->create(pL,npmax,"collide:pL");
      memory->create(pLU,npmax,"collide:pLU");
    }

    // build particle list and find maximum particle weight
    // particle weights are relative to update->fnum

    ip = cinfo[icell].first;
    n = 0;
    sweight_max = 0.0;
    while (ip >= 0) {
      plist[n++] = ip;

      ipart = &particles[ip];
      isw = ipart->weight;

      sweight_max = MAX(sweight_max,isw);

      if (isw != isw) error->all(FLERR,"Particle has NaN weight");
      if (isw <= 0.0) error->all(FLERR,"Particle has negative or zero weight");
      ip = next[ip];
    }
    sweight_max *= update->fnum;

    // attempt = exact collision attempt count for all particles in cell
    // nattempt = rounded attempt with RN
    // if no attempts, continue to next grid cell

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

      if (np >= Ncmin && Ncmin > 0.0) pre_wtf = 0.0;
      else pre_wtf = 1.0;

      newp = split(plist[i],plist[j],ipart,jpart,kpart,lpart);

      // add new particles to particle list

      if (newp > 1) {
        if (np+2 >= npmax) {
          npmax += DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
          memory->grow(pL,npmax,"collide:pL");
          memory->grow(pLU,npmax,"collide:pLU");
        }
        plist[np++] = particle->nlocal-2;
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
      } else if (newp > 0) {
        if (np+1 >= npmax) {
          npmax += DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
          memory->grow(pL,npmax,"collide:pL");
          memory->grow(pLU,npmax,"collide:pLU");
        }
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
      }

      // since ipart and jpart have same weight, do not need
      // ... to account for weight during collision itself
      // also the splits are all handled beforehand

      mpart = NULL; // dummy particle
      setup_collision(ipart,jpart);
      perform_collision(ipart,jpart,mpart);
      ncollide_one++;

    } // end attempt loop
  } // loop for cells

  // remove tiny weighted particles

  remove_tiny();

  return;
}

/* ----------------------------------------------------------------------
   Splits particles and generates two new particles (for SWPM)
------------------------------------------------------------------------- */

int Collide::split(int i, int j,
                   Particle::OnePart *&ip, Particle::OnePart *&jp,
                   Particle::OnePart *&kp, Particle::OnePart *&lp)
{
  double xk[3],vk[3];
  double xl[3],vl[3];
  double erotk, erotl;
  int ks, ls;
  int kcell, lcell;

  // checks if particles properly deleted

  int id;
  int reallocflag;
  Particle::OnePart *particles = particle->particles;

  kp = NULL;
  lp = NULL;

  // weight transfer function is assumed to be
  // ... MIN(ip->sweight,jp->sweight)/(1 + pre_wtf * wtf)

  double isw = ip->weight;
  double jsw = jp->weight;
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

  // update weights

  ip->weight = Gwtf;
  jp->weight = Gwtf;

  // Gwtf should never be negative or zero

  if (Gwtf <= 0.0)
    error->one(FLERR,"Negative weight assigned after split");

  if (Gwtf > 0.0 && pre_wtf > 0.0)
    if (ksw <= 0.0 || lsw <= 0.0)
      error->one(FLERR,"Zero or negative weight after split");

  // number of new particles

  int newp = 0;

  // gk is always the bigger of the two

  if(ksw > 0) {
    id = MAXSMALLINT*random->uniform();
    reallocflag = particle->add_particle(id,ks,kcell,xk,vk,erotk,0.0);
    if (reallocflag) {
      particles = particle->particles;
      ip = particle->particles + (ip - particles);
      jp = particle->particles + (jp - particles);
    }
    kp = &particle->particles[particle->nlocal-1];
    kp->weight = ksw;
    newp++;
  }

  // there should never be case where you add particle "l" if
  // ... you did not add particle "k"

  if(lsw > 0) {
    if(ksw <= 0) error->one(FLERR,"Bad addition to particle list");
    id = MAXSMALLINT*random->uniform();
    reallocflag = particle->add_particle(id,ls,lcell,xl,vl,erotl,0.0);
    if (reallocflag) {
      particles = particle->particles;
      ip = particle->particles + (ip - particles);
      jp = particle->particles + (jp - particles);
      kp = particle->particles + (kp - particles);
    }
    lp = &particle->particles[particle->nlocal-1];
    lp->weight = lsw;
    newp++;
  }

  return newp;
}

/* ----------------------------------------------------------------------
   Reorder plist depending on grouping strategy used and prepare for
   grouping and particle reduction
------------------------------------------------------------------------- */

void Collide::group_reduce()
{
  int n,nold,np,ip;
  double isw;

  Particle::OnePart *ipart;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::ChildCell *cells = grid->cells;
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  double swmean, swvar, swstd;
  double d1, d2;
  double lLim, uLim;
  int npL, npLU;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;

    if (np <= Ncmax) continue;

    // create particle list

    ip = cinfo[icell].first;
    n = 0;
    while (ip >= 0) {
      ipart = &particles[ip];
      isw = ipart->weight;
      if(isw > 0) plist[n++] = ip;
      ip = next[ip];
    }

    gbuf = 0;

    while (n > Ncmax) {
      nold = n;

      // seems to be more stable than weighted

      if (group_type == BINARY) {
        group_bt(plist,n);
      } else if (group_type == WEIGHT) {

        // find mean / standard deviation of weight

        ip = cinfo[icell].first;
        n = 0;
        swmean = swvar = 0.0;
        while (ip >= 0) {
          ipart = &particles[ip];
          isw = ipart->weight;

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

        ip = cinfo[icell].first;
        npL = npLU = 0;
        while (ip >= 0) {
          ipart = &particles[ip];
          isw = ipart->weight;
          if(isw > 0 && isw < lLim) pL[npL++] = ip;
          else if(isw >= lLim && isw < uLim) pLU[npLU++] = ip;
          ip = next[ip];
        }

        // can reuse binary tree division here

        group_bt(pL,  npL);
        group_bt(pLU, npLU);

      }

      // recreate particle list after reduction

      ip = cinfo[icell].first;
      n = 0;
      while (ip >= 0) {
        ipart = &particles[ip];
        isw = ipart->weight;
        if(isw > 0) plist[n++] = ip;
        ip = next[ip];
      }

      // if no particles reduced, increase group size
      
      if (n == nold) gbuf += 2;
      if (gbuf > n) break;

    } // while loop for n > ncmax
  }// loop for cells

  return;
}

/* ----------------------------------------------------------------------
   Recursivley divides particles using the binary tree strategy
------------------------------------------------------------------------- */
void Collide::group_bt(int *plist_leaf, int np)
{
  Particle::OnePart *ipart;
  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;

  // ignore groups which have too few particles

  if (np <= Ngmin) return;

  // compute stress tensor since it's needed for
  // .. further dividing and reduction

  double gsum, msum, mV[3], mVV[3][3], mVVV[3][3];
  gsum = msum = 0.0;
  for (int i = 0; i < 3; i++) {
    mV[i] = 0.0;
    for (int j = 0; j < 3; j++) {
      mVV[i][j] = 0.0;
      mVVV[i][j] = 0.0;
    }
  }

  // find maximum particle weight

  int ispecies;
	double mass, psw, pmsw, vp[3];
  double Erot;
  for (int p = 0; p < np; p++) {
    ipart = &particles[plist_leaf[p]];
    ispecies = ipart->ispecies;
    mass = species[ispecies].mass;

    psw = ipart->weight;
    pmsw = psw * mass;
    memcpy(vp, ipart->v, 3*sizeof(double));
   	gsum += psw;
    msum += pmsw;
    Erot += psw*ipart->erot;
    for (int i = 0; i < 3; i++) {
      mV[i] += (pmsw*vp[i]);
      for (int j = 0; j < 3; j++) {
        mVV[i][j] += (pmsw*vp[i]*vp[j]);
        mVVV[i][j] += (pmsw*vp[i]*vp[j]*vp[j]);
      }
    }
  }

  // mean velocity

	double V[3];
  for (int i = 0; i < 3; i++) V[i] = mV[i]/msum;

  // stress tensor

  double pij[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      pij[i][j] = mVV[i][j] - mV[i]*mV[j]/msum;

  // if group is small enough, merge the particles

  if (np <= Ngmax+gbuf) {

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
      reduce(plist_leaf, np, gsum, V, T, Erot);
    } else if (reduction_type == HEAT) {
      reduce(plist_leaf, np, gsum, V, T, Erot, q);
    } else if (reduction_type == STRESS) {
      reduce(plist_leaf, np, gsum, V, T, Erot, q, pij);
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
        for (int j = 0; j < 3; j++) {
          maxevec[j] = evec[j][i];  
        }
      }
    }

    // Separate based on particle velocity

    double center = V[0]*maxevec[0] + V[1]*maxevec[1] + V[2]*maxevec[2];
    int pid, pidL[np], pidR[np];
    int npL, npR;
    npL = npR = 0;
    for (int i = 0; i < np; i++) {
      pid = plist_leaf[i];
      ipart = &particles[pid];
      if (MathExtra::dot3(ipart->v,maxevec) < center)
        pidL[npL++] = pid;
      else
        pidR[npR++] = pid;
    }

    if(npL > Ngmin) group_bt(pidL,npL);
    if(npR > Ngmin) group_bt(pidR,npR);
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using energy scheme
------------------------------------------------------------------------- */
void Collide::reduce(int *pleaf, int np,
                     double rho, double *V, double T, double Erot)
{

  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  int ip, jp;
  ip = np * random->uniform();
  jp = np * random->uniform();
  while (ip == jp) jp = np * random->uniform();

  ipart = &particles[pleaf[ip]];
  jpart = &particles[pleaf[jp]];

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

  for (int i = 0; i < np; i++) {
    if (i == ip || i == jp) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[pleaf[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = pleaf[i];
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using heat flux scheme
------------------------------------------------------------------------- */
void Collide::reduce(int *pleaf, int np,
                     double rho, double *V, double T, double Erot, double *q)
{

  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  int ip, jp;
  ip = np * random->uniform();
  jp = np * random->uniform();
  while (ip == jp) jp = np * random->uniform();

  ipart = &particles[pleaf[ip]];
  jpart = &particles[pleaf[jp]];

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
  for (int i = 0; i < np; i++) {
    if (i == ip || i == jp) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[pleaf[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = pleaf[i];
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using stress scheme
------------------------------------------------------------------------- */
void Collide::reduce(int *pleaf, int np,
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

    ipart = &particles[pleaf[2*iK]];
    jpart = &particles[pleaf[2*iK+1]];

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
  for (int i = 0; i < np; i++) {
    if (i < 2*nK) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[pleaf[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = pleaf[i];
  }

  return;
}

/* ----------------------------------------------------------------------
   Delete tiny weighted particles
------------------------------------------------------------------------- */
void Collide::remove_tiny()
{
  int np, ip, n;
  double isw, sw_mean;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int *next = particle->next;

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;

    ip = cinfo[icell].first;
    n = 0;
    sw_mean = 0.0;
    while (ip >= 0) {
      ipart = &particles[ip];
      isw = ipart->weight;
      if (isw > 0) sw_mean += isw;
      ip = next[ip];
      n++;
    }

    sw_mean /= n;

    // delete tiny weights

    ip = cinfo[icell].first;
    while (ip >= 0) {
      if (isw < sw_mean*1e-5) {
        if (ndelete == maxdelete) {
          maxdelete += DELTADELETE;
          memory->grow(dellist,maxdelete,"collide:dellist");
        }
        ipart = &particles[ip];
        ipart->weight = -1.0;
        dellist[ndelete++] = ip;
      }
      ip = next[ip];
    }
  } // loop for cells

  return;
}

/* ----------------------------------------------------------------------
   DEBUG : Checks particles
------------------------------------------------------------------------- */

/*template <int ORD> void Collide::check_particles()
{
  //if(ORD==0) printf("before collide\n");
  //else if(ORD==1) printf("after collide\n");
  //else if(ORD==2) printf("after sort\n");
  //else if(ORD==3) printf("after reduce\n");
  //else if(ORD==4) printf("after compress\n");

  int i,j,m,n,ip,np;
  int nattempt,reactflag;
  double attempt,volume;
  Particle::OnePart *ipart;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  double rho;
  double E;
  double V[3], vsq;
  E = rho = 0.0;
  V[0] = V[1] = V[2] = 0.0;
  for (int icell = 0; icell < nglocal; icell++) {
    // create particle list
    ip = cinfo[icell].first;
    while (ip >= 0) {
      ipart = &particles[ip];
      if(ipart->weight > 0) {
        rho += ipart->weight;
        vsq = ipart->v[0]*ipart->v[0]+
              ipart->v[1]*ipart->v[1]+
              ipart->v[2]*ipart->v[2];
        for (int d = 0; d < 3; d++) V[d] += ipart->weight*ipart->v[d];
        E += ipart->weight*vsq;
      }
      ip = next[ip];
    }
  }

  if(update->ntimestep == 1) {
    rho0 = rho;
    for(int d = 0; d < 3; d++) V0[d] = V[d];
    E0 = E;
  }

  if(update->ntimestep % 1000 == 0)
    printf("drho: %2.8e; dE: %2.8e; dV: %2.8e,%2.8e,%2.8e\n",
          1-rho/rho0,
          1-E/E0,
          1-V[0]/V0[0],
          1-V[1]/V0[1],
          1-V[2]/V0[2]);
  return;
}*/
