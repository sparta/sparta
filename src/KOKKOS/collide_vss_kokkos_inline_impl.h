/* ----------------------------------------------------------------------
   Inline device function bodies shared between collide_vss_kokkos.cpp
   and collide_swpm_kokkos.cpp.  Include after collide_vss_kokkos.h.
------------------------------------------------------------------------- */

#ifndef COLLIDE_VSS_KOKKOS_INLINE_IMPL_H
#define COLLIDE_VSS_KOKKOS_INLINE_IMPL_H

#include "math_const.h"

#ifndef EPSZERO
#define EPSZERO 1.0e-14
#endif
#ifndef BIG
#define BIG 1.0e20
#endif

// energy-mode / relaxation-mode enums, mirrored from
// collide_vss_kokkos.cpp so the shared inline bodies below compile in
// any translation unit that includes this header
namespace SPARTA_NS {
enum{NONE,DISCRETE,SMOOTH};
enum{CONSTANT,VARIABLE};
}

namespace SPARTA_NS {

// bodies below use bare MathConst symbols (MY_PI, ...), matching the
// resolution in collide_vss_kokkos.cpp (which has `using namespace MathConst`)
using namespace MathConst;

/* ----------------------------------------------------------------------
   determine if collision actually occurs
   1 = yes, 0 = no
   update vremax either way
   ijsw = SWPM stochastic-weight rejection factor MAX(isw,jsw)/max_sw;
          1.0 for non-SWPM collisions
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::test_collision_kokkos(int icell, int igroup, int jgroup,
                                     Particle::OnePart *ip, Particle::OnePart *jp,
                                     struct State &precoln, rand_type &rand_gen,
                                     double ijsw) const
{
  double *vi = ip->v;
  double *vj = jp->v;
  int ispecies = ip->ispecies;
  int jspecies = jp->ispecies;
  double du  = vi[0] - vj[0];
  double dv  = vi[1] - vj[1];
  double dw  = vi[2] - vj[2];
  double vr2 = du*du + dv*dv + dw*dw;

  // prevent division by zero

  if (vr2 < EPSZERO && d_params(ispecies,jspecies).omega >= 1.0)
    return 0;

  double vro  = pow(vr2,1.0-d_params(ispecies,jspecies).omega);

  // although the vremax is calculated for the group,
  // the individual collisions calculated species dependent vre

  double vre = vro*d_prefactor(ispecies,jspecies);
  d_vremax(icell,igroup,jgroup) = MAX(vre,d_vremax(icell,igroup,jgroup));
  // ijsw accounts for the SWPM stochastic-weight rejection factor
  // (MAX(isw,jsw)/max_stochastic_weight); it is 1.0 for non-SWPM collisions
  if (vre/d_vremax(icell,igroup,jgroup)*ijsw < rand_gen.drand()) return 0;
  precoln.vr2 = vr2;
  return 1;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::setup_collision_kokkos(Particle::OnePart *ip, Particle::OnePart *jp,
                                       struct State &precoln, struct State &postcoln) const
{
  int isp = ip->ispecies;
  int jsp = jp->ispecies;

  precoln.vr = sqrt(precoln.vr2);

  precoln.ave_rotdof = 0.5 * (d_species[isp].rotdof + d_species[jsp].rotdof);
  precoln.ave_vibdof = 0.5 * (d_species[isp].vibdof + d_species[jsp].vibdof);
  precoln.ave_dof = (precoln.ave_rotdof  + precoln.ave_vibdof)/2.;

  precoln.imass = d_species[isp].mass;
  precoln.jmass = d_species[jsp].mass;

  precoln.etrans = 0.5 * d_params(isp,jsp).mr * precoln.vr2;
  precoln.erot = ip->erot + jp->erot;
  precoln.evib = ip->evib + jp->evib;

  precoln.eint   = precoln.erot + precoln.evib;
  precoln.etotal = precoln.etrans + precoln.eint;

  // COM velocity calculated using reactant masses

  double divisor = 1.0 / (d_species[isp].mass + d_species[jsp].mass);
  double *vi = ip->v;
  double *vj = jp->v;
  precoln.ucmf = ((d_species[isp].mass*vi[0])+(d_species[jsp].mass*vj[0]))*divisor;
  precoln.vcmf = ((d_species[isp].mass*vi[1])+(d_species[jsp].mass*vj[1]))*divisor;
  precoln.wcmf = ((d_species[isp].mass*vi[2])+(d_species[jsp].mass*vj[2]))*divisor;

  postcoln.etrans = precoln.etrans;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  postcoln.eint = 0.0;
  postcoln.etotal = precoln.etotal;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::find_nn(rand_type &rand_gen, int i, int np, int icell) const
{
  int jneigh;
  double dx,dy,dz,rsq;
  double *xj;

  // if np = 2, just return J = non-I particle
  // np is never < 2

  if (np == 2) return (i+1) % 2;

  Particle::OnePart *ipart,*jpart;

  // thresh = distance particle I moves in this timestep

  ipart = &d_particles[d_plist(icell,i)];
  double *vi = ipart->v;
  double *xi = ipart->x;
  double threshsq =  dt*dt * (vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2]);
  double minrsq = BIG;

  // nlimit = max # of J candidates to consider

  int nlimit = MIN(nearlimit,np-1);
  int count = 0;

  // pick a random starting J
  // jneigh = collision partner when exit loop
  //   set to initial J as default in case no Nlimit J meets criteria

  int j = np * rand_gen.drand();
  while (i == j) j = np * rand_gen.drand();
  jneigh = j;

  while (count < nlimit) {
    count++;

    // skip this J if I,J last collided with each other

    if (d_nn_last_partner(icell,i) == j+1 && d_nn_last_partner(icell,j) == i+1) {
      j++;
      if (j == np) j = 0;
      continue;
    }

    // rsq = squared distance between particles I and J

    jpart = &d_particles[d_plist(icell,j)];
    xj = jpart->x;
    dx = xi[0] - xj[0];
    dy = xi[1] - xj[1];
    dz = xi[2] - xj[2];
    rsq = dx*dx + dy*dy + dz*dz;

    if (rsq > 0.0) {
      if (rsq <= threshsq) {
        jneigh = j;
        break;
      }
      if (rsq < minrsq) {
        minrsq = rsq;
        jneigh = j;
      }
    }
    j++;
    if (j == np) j = 0;
  }

  return jneigh;
}


/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
int CollideVSSKokkos::perform_collision_kokkos(Particle::OnePart *&ip,
                                  Particle::OnePart *&jp,
                                  Particle::OnePart *&kp,
                                  struct State &precoln, struct State &postcoln, rand_type &rand_gen,
                                  Particle::OnePart *&p3, int &recomb_species, double &recomb_density,
                                  int &index_kpart) const
{
  int reaction,kspecies;
  double x[3],v[3];

  // if gas-phase chemistry defined, attempt and perform reaction
  // if a 3rd particle is created, its kspecies >= 0 is returned
  // if 2nd particle is removed, its jspecies is set to -1
  // reaction = 0 if no reaction occurs
  // reaction = 1 to N for which reaction occurs
  // reaction is returned to caller

  if (react_defined)
    reaction = react_kk_copy.obj.attempt_kk(ip,jp,
                                             precoln.etrans,precoln.erot,
                                             precoln.evib,postcoln.etotal,kspecies,
                                             recomb_species,recomb_density,d_species);
  else reaction = 0;

  // just collision, no reaction

  if (!reaction) {
    if (precoln.ave_dof > 0.0) EEXCHANGE_NonReactingEDisposal(ip,jp,precoln,postcoln,rand_gen);
    SCATTER_TwoBodyScattering(ip,jp,precoln,postcoln,rand_gen);
    return reaction;
  }

  // reaction took place
  // repartition energy and perform velocity scattering for I,J,K particles
  // reaction may have changed species of I,J particles
  // J,K particles may have been removed or created by reaction

  kp = NULL;

  // add 3rd K particle if reaction created it
  // index of new K particle = nlocal-1
  // if add_particle() performs a realloc:
  //   make copy of x,v

  if (kspecies >= 0) {
    int id = MAXSMALLINT*rand_gen.drand();

    memcpy(x,ip->x,3*sizeof(double));
    memcpy(v,ip->v,3*sizeof(double));
    index_kpart = Kokkos::atomic_fetch_add(&d_nlocal(),1);
    int reallocflag =
      ParticleKokkos::add_particle_kokkos(d_particles,index_kpart,id,kspecies,ip->icell,x,v,0.0,0.0);
    if (reallocflag) {
      d_retry() = 1;
      d_part_grow() = 1;
      return 0;
    }

    kp = &d_particles[index_kpart];
    EEXCHANGE_ReactingEDisposal(ip,jp,kp,precoln,postcoln,rand_gen);
    SCATTER_ThreeBodyScattering(ip,jp,kp,precoln,postcoln,rand_gen);

  // remove 2nd J particle if recombination reaction removed it
  // p3 is 3rd particle participating in energy exchange

  } else if (jp->ispecies < 0) {
    double *vi = ip->v;
    double *vj = jp->v;

    const double divisor = 1.0 / (precoln.imass + precoln.jmass);
    const double ucmf = ((precoln.imass*vi[0]) + (precoln.jmass*vj[0])) * divisor;
    const double vcmf = ((precoln.imass*vi[1]) + (precoln.jmass*vj[1])) * divisor;
    const double wcmf = ((precoln.imass*vi[2]) + (precoln.jmass*vj[2])) * divisor;

    vi[0] = ucmf;
    vi[1] = vcmf;
    vi[2] = wcmf;

    jp = NULL;

    // account for 3rd body energy via another call to setup_collision()
    // set precoln.vr2 = relative velocity between ip and 3rd body p3

    const double *vp3 = p3->v;
    const double du  = vi[0] - vp3[0];
    const double dv  = vi[1] - vp3[1];
    const double dw  = vi[2] - vp3[2];
    const double vr2 = du*du + dv*dv + dw*dw;
    precoln.vr2 = vr2;

    // save postcoln.etotal from previous setup_collision()
    // add 3rd body internal energy to it
    // ip internal energy is already included in postcoln.etotal

    double partial_energy =  postcoln.etotal + p3->erot + p3->evib;
    ip->erot = 0.0;
    ip->evib = 0.0;
    p3->erot = 0.0;
    p3->evib = 0.0;

    // 2nd call to setup_collision() sets new postcoln.etotal
    // then add saved partial_energy to it

    setup_collision_kokkos(ip,p3,precoln,postcoln);
    postcoln.etotal += partial_energy;

    if (precoln.ave_dof > 0.0) EEXCHANGE_ReactingEDisposal(ip,p3,jp,precoln,postcoln,rand_gen);
    SCATTER_TwoBodyScattering(ip,p3,precoln,postcoln,rand_gen);

  } else {
    EEXCHANGE_ReactingEDisposal(ip,jp,kp,precoln,postcoln,rand_gen);
    SCATTER_TwoBodyScattering(ip,jp,precoln,postcoln,rand_gen);
  }

  return reaction;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::SCATTER_TwoBodyScattering(Particle::OnePart *ip,
                                                 Particle::OnePart *jp,
                                                 struct State &precoln, struct State &postcoln,
                                                 rand_type &rand_gen) const
{
  double ua,vb,wc;
  double vrc[3];

  double *vi = ip->v;
  double *vj = jp->v;
  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  double mass_i = d_species[isp].mass;
  double mass_j = d_species[jsp].mass;

  double alpha_r = 1.0 / d_params(isp,jsp).alpha;

  double eps = rand_gen.drand() * 2*MY_PI;
  if (fabs(alpha_r - 1.0) < 0.001) {
    double vr = sqrt(2.0 * postcoln.etrans / d_params(isp,jsp).mr);
    double cosX = 2.0*rand_gen.drand() - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0 * postcoln.etrans) / (d_params(isp,jsp).mr * precoln.vr2));
    double cosX = 2.0*pow(rand_gen.drand(),alpha_r) - 1.0;
    double sinX = sqrt(1.0 - cosX*cosX);
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.0e-6) {
      ua = scale * ( cosX*vrc[0] + sinX*d*sin(eps) );
      vb = scale * ( cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) -
                                         vrc[0]*vrc[1]*sin(eps))/d );
      wc = scale * ( cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) +
                                         vrc[0]*vrc[2]*sin(eps))/d );
    } else {
      ua = scale * ( cosX*vrc[0] );
      vb = scale * ( sinX*vrc[0]*cos(eps) );
      wc = scale * ( sinX*vrc[0]*sin(eps) );
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_i + mass_j);
  vi[0] = precoln.ucmf + (mass_j*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_j*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_j*divisor)*wc;
  vj[0] = precoln.ucmf - (mass_i*divisor)*ua;
  vj[1] = precoln.vcmf - (mass_i*divisor)*vb;
  vj[2] = precoln.wcmf - (mass_i*divisor)*wc;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::EEXCHANGE_NonReactingEDisposal(Particle::OnePart *ip,
                                                      Particle::OnePart *jp,
                                                      struct State &precoln, struct State &postcoln,
                                                      rand_type &rand_gen) const
{
  double State_prob,Fraction_Rot,Fraction_Vib,E_Dispose;
  int i,rotdof,vibdof,max_level,ivib;

  Particle::OnePart *p;

  double AdjustFactor = 0.99999999;
  postcoln.erot = 0.0;
  postcoln.evib = 0.0;
  double pevib = 0.0;

  // handle each kind of energy disposal for non-reacting reactants

  if (precoln.ave_dof == 0) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;

  } else {
    E_Dispose = precoln.etrans;

    for (i = 0; i < 2; i++) {
      if (i == 0) p = ip;
      else p = jp;

      int sp = p->ispecies;
      rotdof = d_species[sp].rotdof;
      double rotn_phi = d_species[sp].rotrel;

      if (rotdof) {
        if (relaxflag == VARIABLE) rotn_phi = rotrel(sp,E_Dispose+p->erot);
        if (rotn_phi >= rand_gen.drand()) {
          if (rotstyle == NONE) {
            p->erot = 0.0 ;

          } else if (rotstyle != NONE && rotdof == 2) {
            E_Dispose += p->erot;
            Fraction_Rot =
              1- pow(rand_gen.drand(),
                     (1/(2.5-d_params(ip->ispecies,jp->ispecies).omega)));
            p->erot = Fraction_Rot * E_Dispose;
            E_Dispose -= p->erot;
          } else {
            E_Dispose += p->erot;
            p->erot = E_Dispose *
              sample_bl(rand_gen,0.5*d_species[sp].rotdof-1.0,
                        1.5-d_params(ip->ispecies,jp->ispecies).omega);
            E_Dispose -= p->erot;
          }
        }
      }
      postcoln.erot += p->erot;

      vibdof = d_species[sp].vibdof;
      double vibn_phi = d_species[sp].vibrel[0];

      if (vibdof) {
        if (relaxflag == VARIABLE) vibn_phi = vibrel(sp,E_Dispose+p->evib);
        if (vibn_phi >= rand_gen.drand()) {
          if (vibstyle == NONE) {
            p->evib = 0.0;

          } else if (vibdof == 2) {
            if (vibstyle == SMOOTH) {
              E_Dispose += p->evib;
              Fraction_Vib =
                1.0 - pow(rand_gen.drand(),(1.0/(2.5-d_params(ip->ispecies,jp->ispecies).omega)));
              p->evib= Fraction_Vib * E_Dispose;
              E_Dispose -= p->evib;

            } else if (vibstyle == DISCRETE) {
              E_Dispose += p->evib;
              max_level = static_cast<int>
                (E_Dispose / (boltz * d_species[sp].vibtemp[0]));
              do {
                ivib = static_cast<int>
                  (rand_gen.drand()*(max_level+AdjustFactor));
                p->evib = ivib * boltz * d_species[sp].vibtemp[0];
                State_prob = pow((1.0 - p->evib / E_Dispose),
                                 (1.5 - d_params(ip->ispecies,jp->ispecies).omega));
              } while (State_prob < rand_gen.drand());
              E_Dispose -= p->evib;
            }
          } else if (vibdof > 2) {
            if (vibstyle == SMOOTH) {
              E_Dispose += p->evib;
              p->evib = E_Dispose *
                sample_bl(rand_gen,0.5*d_species[sp].vibdof-1.0,
                          1.5-d_params(ip->ispecies,jp->ispecies).omega);
              E_Dispose -= p->evib;

            } else if (vibstyle == DISCRETE) {
              p->evib = 0.0;

              int nmode = d_species[sp].nvibmode;
              const auto &d_vibmode = k_eiarray.view_device()[d_ewhich[index_vibmode]].k_view.view_device();
              int pindex = p - d_particles.data();

              for (int imode = 0; imode < nmode; imode++) {
                ivib = d_vibmode(pindex,imode);
                E_Dispose += ivib * boltz *
                  d_species[sp].vibtemp[imode];
                max_level = static_cast<int>
                  (E_Dispose / (boltz * d_species[sp].vibtemp[imode]));

                do {
                  ivib = static_cast<int>
                    (rand_gen.drand()*(max_level+AdjustFactor));
                  pevib = ivib * boltz * d_species[sp].vibtemp[imode];
                  State_prob = pow((1.0 - pevib / E_Dispose),
                                   (1.5 - d_params(ip->ispecies,jp->ispecies).omega));
                } while (State_prob < rand_gen.drand());

                d_vibmode(pindex,imode) = ivib;
                p->evib += pevib;
                E_Dispose -= pevib;
              }
            }
          } // end of vibstyle/vibdof if
        }
        postcoln.evib += p->evib;
      } // end of vibdof if
    }
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::SCATTER_ThreeBodyScattering(Particle::OnePart *ip,
                                                   Particle::OnePart *jp,
                                                   Particle::OnePart *kp,
                                                   struct State &precoln, struct State &postcoln,
                                                   rand_type &rand_gen) const
{
  double vrc[3],ua,vb,wc;

  int isp = ip->ispecies;
  int jsp = jp->ispecies;
  int ksp = kp->ispecies;
  double mass_i = d_species[isp].mass;
  double mass_j = d_species[jsp].mass;
  double mass_k = d_species[ksp].mass;
  double mass_ij = mass_i + mass_j;
  double *vi = ip->v;
  double *vj = jp->v;
  double *vk = kp->v;

  double alpha_r = 1.0 / d_params(isp,jsp).alpha;
  double mr = mass_ij * mass_k / (mass_ij + mass_k);
  postcoln.eint = ip->erot + jp->erot + ip->evib + jp->evib
                + kp->erot + kp->evib;

  double cosX = 2.0*pow(rand_gen.drand(), alpha_r) - 1.0;
  double sinX = sqrt(1.0 - cosX*cosX);
  double eps = rand_gen.drand() * 2*MY_PI;

  if (fabs(alpha_r - 1.0) < 0.001) {
    double vr = sqrt(2*postcoln.etrans/mr);
    ua = vr*cosX;
    vb = vr*sinX*cos(eps);
    wc = vr*sinX*sin(eps);
  } else {
    double scale = sqrt((2.0*postcoln.etrans) / (mr*precoln.vr2));
    vrc[0] = vi[0]-vj[0];
    vrc[1] = vi[1]-vj[1];
    vrc[2] = vi[2]-vj[2];
    double d = sqrt(vrc[1]*vrc[1]+vrc[2]*vrc[2]);
    if (d > 1.E-6 ) {
      ua = scale * (cosX*vrc[0] + sinX*d*sin(eps));
      vb = scale * (cosX*vrc[1] + sinX*(precoln.vr*vrc[2]*cos(eps) -
                                        vrc[0]*vrc[1]*sin(eps))/d);
      wc = scale * (cosX*vrc[2] - sinX*(precoln.vr*vrc[1]*cos(eps) +
                                        vrc[0]*vrc[2]*sin(eps))/d);
    } else {
      ua = scale * cosX*vrc[0];
      vb = scale * sinX*vrc[0]*cos(eps);
      wc = scale * sinX*vrc[0]*sin(eps);
    }
  }

  // new velocities for the products

  double divisor = 1.0 / (mass_ij + mass_k);
  vi[0] = precoln.ucmf + (mass_k*divisor)*ua;
  vi[1] = precoln.vcmf + (mass_k*divisor)*vb;
  vi[2] = precoln.wcmf + (mass_k*divisor)*wc;
  vk[0] = precoln.ucmf - (mass_ij*divisor)*ua;
  vk[1] = precoln.vcmf - (mass_ij*divisor)*vb;
  vk[2] = precoln.wcmf - (mass_ij*divisor)*wc;
  vj[0] = vi[0];
  vj[1] = vi[1];
  vj[2] = vi[2];
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void CollideVSSKokkos::EEXCHANGE_ReactingEDisposal(Particle::OnePart *ip,
                                                   Particle::OnePart *jp,
                                                   Particle::OnePart *kp,
                                                   struct State &precoln, struct State &postcoln,
                                                   rand_type &rand_gen) const
{
  double State_prob,Fraction_Rot,Fraction_Vib;
  int i,numspecies,rotdof,vibdof,max_level,ivib;
  double aveomega,pevib;

  Particle::OnePart *p;
  double AdjustFactor = 0.99999999;

  if (!kp) {
    ip->erot = 0.0;
    jp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    numspecies = 2;
    aveomega = d_params(ip->ispecies,jp->ispecies).omega;
  } else {
    ip->erot = 0.0;
    jp->erot = 0.0;
    kp->erot = 0.0;
    ip->evib = 0.0;
    jp->evib = 0.0;
    kp->evib = 0.0;
    numspecies = 3;
    aveomega = (d_params(ip->ispecies,ip->ispecies).omega + d_params(jp->ispecies,jp->ispecies).omega +
                d_params(kp->ispecies,kp->ispecies).omega)/3.0;
  }

  // handle each kind of energy disposal for non-reacting reactants
  // clean up memory for the products

  double E_Dispose = postcoln.etotal;

  for (i = 0; i < numspecies; i++) {
    if (i == 0) p = ip;
    else if (i == 1) p = jp;
    else p = kp;

    int sp = p->ispecies;
    rotdof = d_species[sp].rotdof;

    if (rotdof) {
      if (rotstyle == NONE) {
        p->erot = 0.0 ;
      } else if (rotdof == 2) {
        Fraction_Rot =
          1- pow(rand_gen.drand(),(1/(2.5-aveomega)));
        p->erot = Fraction_Rot * E_Dispose;
        E_Dispose -= p->erot;

      } else if (rotdof > 2) {
        p->erot = E_Dispose *
          sample_bl(rand_gen,0.5*d_species[sp].rotdof-1.0,
                    1.5-aveomega);
        E_Dispose -= p->erot;
      }
    }

    vibdof = d_species[sp].vibdof;

    if (vibdof) {
      if (vibstyle == NONE) {
        p->evib = 0.0;
      } else if (vibdof == 2 && vibstyle == DISCRETE) {
        max_level = static_cast<int>
          (E_Dispose / (boltz * d_species[sp].vibtemp[0]));
        do {
          ivib = static_cast<int>
            (rand_gen.drand()*(max_level+AdjustFactor));
          p->evib = (double)
            (ivib * boltz * d_species[sp].vibtemp[0]);
          State_prob = pow((1.0 - p->evib / E_Dispose),
                           (1.5 - aveomega));
        } while (State_prob < rand_gen.drand());
        E_Dispose -= p->evib;

      } else if (vibdof == 2 && vibstyle == SMOOTH) {
        Fraction_Vib =
          1.0 - pow(rand_gen.drand(),(1.0 / (2.5-aveomega)));
        p->evib = Fraction_Vib * E_Dispose;
        E_Dispose -= p->evib;

      } else if (vibdof > 2 && vibstyle == SMOOTH) {
        p->evib = E_Dispose *
          sample_bl(rand_gen,0.5*d_species[sp].vibdof-1.0,
                    1.5-aveomega);
        E_Dispose -= p->evib;
      } else if (vibdof > 2 && vibstyle == DISCRETE) {
        p->evib = 0.0;

        int nmode = d_species[sp].nvibmode;
        const auto &d_vibmode = k_eiarray.view_device()[d_ewhich[index_vibmode]].k_view.view_device();
        int pindex = p - d_particles.data();

        for (int imode = 0; imode < nmode; imode++) {
          ivib = d_vibmode(pindex,imode);
          E_Dispose += ivib * boltz *
          d_species[sp].vibtemp[imode];
          max_level = static_cast<int>
          (E_Dispose / (boltz * d_species[sp].vibtemp[imode]));
          do {
            ivib = static_cast<int>
            (rand_gen.drand()*(max_level+AdjustFactor));
            pevib = ivib * boltz * d_species[sp].vibtemp[imode];
            State_prob = pow((1.0 - pevib / E_Dispose),
                             (1.5 - aveomega));
          } while (State_prob < rand_gen.drand());

          d_vibmode(pindex,imode) = ivib;
          p->evib += pevib;
          E_Dispose -= pevib;
        }
      }
    }
  }

  // compute post-collision internal energies

  postcoln.erot = ip->erot + jp->erot;
  postcoln.evib = ip->evib + jp->evib;

  if (kp) {
    postcoln.erot += kp->erot;
    postcoln.evib += kp->evib;
  }

  // compute portion of energy left over for scattering

  postcoln.eint = postcoln.erot + postcoln.evib;
  postcoln.etrans = E_Dispose;
}

/* ---------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::sample_bl(rand_type &rand_gen, double Exp_1, double Exp_2) const
{
  double Exp_s = Exp_1 + Exp_2;
  double x,y;
  do {
    x = rand_gen.drand();
    y = pow(x*Exp_s/Exp_1, Exp_1)*pow((1.0-x)*Exp_s/Exp_2, Exp_2);
  } while (y < rand_gen.drand());
  return x;
}

/* ----------------------------------------------------------------------
   compute a variable rotational relaxation parameter
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::rotrel(int isp, double Ec) const
{
  // Because we are only relaxing one of the particles in each call, we only
  //  include its DoF, consistent with Bird 2013 (3.32)

  double Tr = Ec /(boltz * (2.5-d_params(isp,isp).omega + d_species[isp].rotdof/2.0));
  double rotphi = (1.0+d_params(isp,isp).rotc2/sqrt(Tr) + d_params(isp,isp).rotc3/Tr)
                / d_params(isp,isp).rotc1;
  return rotphi;
}

/* ----------------------------------------------------------------------
   compute a variable vibrational relaxation parameter
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double CollideVSSKokkos::vibrel(int isp, double Ec) const
{
  double Tr = Ec /(boltz * (3.5-d_params(isp,isp).omega));
  double omega = d_params(isp,isp).omega;
  double vibphi = 1.0 / (d_params(isp,isp).vibc1/pow(Tr,omega) *
                         exp(d_params(isp,isp).vibc2/pow(Tr,1.0/3.0)));
  return vibphi;
}


} // namespace SPARTA_NS

#endif // COLLIDE_VSS_KOKKOS_INLINE_IMPL_H
