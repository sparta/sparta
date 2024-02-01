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

#include "mpi.h"
#include "math.h"
#include "ctype.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "surf_react_adsorb.h"
#include "input.h"
#include "update.h"
#include "comm.h"
#include "surf.h"
#include "grid.h"
#include "domain.h"
#include "particle.h"
#include "collide.h"
#include "collide_vss.h"
#include "surf_collide.h"
#include "style_surf_collide.h"
#include "compute_react_surf.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{INT,DOUBLE};                      // several files
enum{FACE,SURF};
enum{PSFACE,PSLINE,PSTRI};

#define DELTA_TALLY 1024
#define DELTA_PART 8                   // make this bigger once debugged

// GS chemistry

enum{DISSOCIATION,EXCHANGE,RECOMBINATION,AA,DA,LH1,LH3,CD,ER,CI};
enum{NOMODEL,SPECULAR,DIFFUSE,ADIABATIC,CLL,TD,IMPULSIVE,MAXMODELS};

#define MAXREACTANT_GS 5
#define MAXPRODUCT_GS 5
#define MAXCOEFF_GS 4

// PS chemistry

enum{DS,LH2,LH4,SB};
//enum{GRID,LINE,TRI};
enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain

#define MAXREACTANT_PS 4
#define MAXPRODUCT_PS 4
#define MAXCOEFF_PS 3

// both kinds of chemistry

enum{GS,PS,GSPS};
enum{SIMPLE,ARRHENIUS};                       // type of reaction
enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};           // several files
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files

#define MAXLINE 1024
#define DELTALIST 10

// Syntax: surf_react ID adsorb (gs or ps or gs/ps) filename1 {filename2}
//                    Nsync (face or surf) Tsurf max_cover O CO CO_a CO_b ...

/* ---------------------------------------------------------------------- */

SurfReactAdsorb::SurfReactAdsorb(SPARTA *sparta, int narg, char **arg) :
  SurfReact(sparta, narg, arg)
{
  if (surf->distributed || surf->implicit)
    error->all(FLERR,
               "Cannot yet use surf_react adsorb with distributed or "
               "implicit surf elements");

  me = comm->me;
  nprocs = comm->nprocs;

  // 1st arg: gas chemistry or surf chemistry or both

  if (narg < 3) error->all(FLERR,"Illegal surf_react adsorb command");

  gsflag = psflag = 0;
  if (strcmp(arg[2],"gs") == 0) gsflag = 1;
  else if (strcmp(arg[2],"ps") == 0) psflag = 1;
  else if (strcmp(arg[2],"gs/ps") == 0) gsflag = psflag = 1;
  else error->all(FLERR,"Illegal surf_react adsorb command");

  int iarg,gs_filearg,ps_filearg;

  if (gsflag && psflag) {
    gs_filearg = 3;
    ps_filearg = 4;
    iarg = 5;
    if (narg < 10) error->all(FLERR,"Illegal surf_react adsorb command");
  } else {
    if (gsflag) gs_filearg = 3;
    else if (psflag) ps_filearg = 3;
    iarg = 4;
    if (narg < 9) error->all(FLERR,"Illegal surf_react adsorb command");
  }

  if (strcmp(arg[iarg],"nsync") != 0)
    error->all(FLERR,"Illegal surf_react adsorb command");
  nsync = input->numeric(FLERR,arg[iarg+1]);
  if (nsync < 1) error->all(FLERR,"Illegal surf_react adsorb command");

  if (strcmp(arg[iarg+2],"face") == 0) mode = FACE;
  else if (strcmp(arg[iarg+2],"surf") == 0) mode = SURF;
  else error->all(FLERR,"Illegal surf_react adsorb command");

  if (mode == SURF && surf->nsurf == 0)
    error->all(FLERR,"Cannot use surf_react adsorb when no surfs exist");

  twall = input->numeric(FLERR,arg[iarg+3]);
  max_cover = input->numeric(FLERR,arg[iarg+4]);

  // species_surf = list of surface species IDs

  iarg += 5;
  species_surf = new char*[narg-iarg];
  nspecies_surf = 0;

  while (iarg < narg) {
    int isp = particle->find_species(arg[iarg]);
    if (isp < 0) error->all(FLERR,"Surf_react adsorb species is not defined");
    int n = strlen(arg[iarg]) + 1;
    species_surf[nspecies_surf] = new char[n];
    strcpy(species_surf[nspecies_surf],arg[iarg]);
    nspecies_surf++;
    iarg++;
  }

  // initialize reaction data structs

  nlist_gs = maxlist_gs = 0;
  rlist_gs = NULL;
  reactions_gs = NULL;
  indices_gs = NULL;

  nlist_ps = maxlist_ps = 0;
  rlist_ps = NULL;
  reactions_ps_list = NULL;
  nactive_ps = 0;
  n_PS_react = 0;

  // initialize PS added particle data structs

  mypart = NULL;
  allpart = NULL;
  maxmypart = maxallpart = 0;

  memory->create(recvcounts,comm->nprocs,"sr_adsorb:recvcounts");
  memory->create(displs,comm->nprocs,"sr_adsorb:displs");

  // list of surface collision models

  cmodels = new SurfCollide*[MAXMODELS];
  for (int i = 0; i < MAXMODELS; i++) cmodels[i] = NULL;

  // read the file defining GS and/or PS reactions

  if (gsflag) readfile_gs(arg[gs_filearg]);
  if (psflag) readfile_ps(arg[ps_filearg]);

  // reaction tally setup

  nlist = 0;
  if (gsflag) nlist += nlist_gs;
  if (psflag) nlist += nlist_ps;
  tally_single = new int[nlist];
  tally_total = new int[nlist];
  tally_single_all = new int[nlist];
  tally_total_all = new int[nlist];

  size_vector = 2 + 2*nlist;

  nsingle = ntotal = 0;

  // initialize RN generator

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  // create and initialize per-face or custom per-surf attributes

  if (mode == FACE) create_per_face_state();
  if (mode == SURF) create_per_surf_state();
}

/* ---------------------------------------------------------------------- */

SurfReactAdsorb::~SurfReactAdsorb()
{
  delete random;

  // surface species

  for (int i = 0; i < nspecies_surf; i++) delete [] species_surf[i];
  delete [] species_surf;

  // GS chemistry

  if (gsflag) {
    for (int i = 0; i < maxlist_gs; i++) {
      for (int j = 0; j < rlist_gs[i].nreactant; j++) {
        delete [] rlist_gs[i].id_reactants[j];
        delete [] rlist_gs[i].state_reactants[j];
      }
      for (int j = 0; j < rlist_gs[i].nproduct; j++) {
        delete [] rlist_gs[i].id_products[j];
        delete [] rlist_gs[i].state_products[j];
      }
      delete [] rlist_gs[i].id;
      delete [] rlist_gs[i].id_reactants;
      delete [] rlist_gs[i].id_products;
      delete [] rlist_gs[i].state_reactants;
      delete [] rlist_gs[i].state_products;
      delete [] rlist_gs[i].part_reactants;
      delete [] rlist_gs[i].part_products;
      delete [] rlist_gs[i].stoich_reactants;
      delete [] rlist_gs[i].stoich_products;
      delete [] rlist_gs[i].reactants;
      delete [] rlist_gs[i].products;
      delete [] rlist_gs[i].reactants_ad_index;
      delete [] rlist_gs[i].products_ad_index;
      delete [] rlist_gs[i].coeff;
      delete [] rlist_gs[i].cmodel_ip_flags;
      delete [] rlist_gs[i].cmodel_ip_coeffs;
      delete [] rlist_gs[i].cmodel_jp_flags;
      delete [] rlist_gs[i].cmodel_jp_coeffs;
    }

    memory->destroy(rlist_gs);
    memory->destroy(reactions_gs);
    memory->destroy(indices_gs);
  }

  // PS chemistry

  if (psflag) {
    for (int i = 0; i < maxlist_ps; i++) {
      for (int j = 0; j < rlist_ps[i].nreactant; j++) {
        delete [] rlist_ps[i].id_reactants[j];
        delete [] rlist_ps[i].state_reactants[j];
      }
      for (int j = 0; j < rlist_ps[i].nproduct; j++) {
        delete [] rlist_ps[i].id_products[j];
        delete [] rlist_ps[i].state_products[j];
      }
      delete [] rlist_ps[i].id_reactants;
      delete [] rlist_ps[i].id_products;
      delete [] rlist_ps[i].state_reactants;
      delete [] rlist_ps[i].state_products;
      delete [] rlist_ps[i].part_reactants;
      delete [] rlist_ps[i].part_products;
      delete [] rlist_ps[i].stoich_reactants;
      delete [] rlist_ps[i].stoich_products;
      delete [] rlist_ps[i].reactants;
      delete [] rlist_ps[i].products;
      delete [] rlist_ps[i].reactants_ad_index;
      delete [] rlist_ps[i].products_ad_index;
      delete [] rlist_ps[i].coeff;
      delete [] rlist_ps[i].id;
      delete [] rlist_ps[i].cmodel_ip_flags;
      delete [] rlist_ps[i].cmodel_ip_coeffs;
      delete [] rlist_ps[i].cmodel_jp_flags;
      delete [] rlist_ps[i].cmodel_jp_coeffs;
    }
    memory->destroy(rlist_ps);
    memory->destroy(reactions_ps_list);

    // added PS particles

    memory->sfree(mypart);
    memory->sfree(allpart);
  }

  // parallel comm

  memory->destroy(recvcounts);
  memory->destroy(displs);

  // surface collision models

  for (int i = 0; i < MAXMODELS; i++)  delete cmodels[i];
  delete [] cmodels;

  // deallocate per-face state data

  if (mode == FACE) {
    memory->destroy(face_species_state);
    memory->destroy(face_total_state);
    memory->destroy(face_area);
    memory->destroy(face_weight);
    if (psflag) memory->destroy(face_tau);
  }

  // delete custom per-surf state data owned by Surf class
  // this must be done exactly once by first_owner
  // in case multiple surf react/adsorb instances are used

  if (mode == SURF && first_owner) {
    surf->remove_custom(nstick_species_custom);
    surf->remove_custom(nstick_total_custom);
    surf->remove_custom(area_custom);
    surf->remove_custom(weight_custom);
    if (psflag) surf->remove_custom(tau_custom);
  }

  // delete local face and per-surf data

  if (mode == FACE) {
    memory->destroy(face_species_delta);
    memory->destroy(face_sum_delta);
    memory->destroy(face_norm);
  }

  if (mode == SURF) {
    memory->destroy(surf_species_delta);
    memory->destroy(mark);
    memory->destroy(tally2surf);
    memory->destroy(intally);
    memory->destroy(outtally);
    memory->destroy(incollate);
    memory->destroy(outcollate);
  }
}

/* ----------------------------------------------------------------------
   create and initialize per-face data structs
------------------------------------------------------------------------- */

void SurfReactAdsorb::create_per_face_state()
{
  // every proc allocates state

  nface = 2 * domain->dimension;
  memory->create(face_species_state,nface,nspecies_surf,"face_species_state");
  memory->create(face_total_state,nface,"face_total_state");
  memory->create(face_area,nface,"face_area");
  memory->create(face_weight,nface,"face_weight");

  // local delta and norm storage

  memory->create(face_species_delta,nface,nspecies_surf,"face_species_delta");
  memory->create(face_sum_delta,nface,nspecies_surf,"face_sum_delta");
  memory->create(face_norm,nface,3,"face_norm");

  // initialize 4 state quantities and face normals
  // initialize face_species_delta for first time

  for (int iface = 0; iface < nface; iface++) {
    face_total_state[iface] = 0;
    for (int isp = 0; isp < nspecies_surf; isp++) {
      face_species_state[iface][isp] = 0;
      face_species_delta[iface][isp] = 0;
    }

    face_norm[iface][0] = face_norm[iface][1] = face_norm[iface][2] = 0.0;
    if (domain->dimension == 2) {
      if (iface < 2) face_area[iface] = domain->prd[1];
      else if (iface < 4) face_area[iface] = domain->prd[0];
      if (iface == 0) face_norm[iface][0] = 1.0;
      else if (iface == 1) face_norm[iface][0] = -1.0;
      else if (iface == 2) face_norm[iface][1] = 1.0;
      else if (iface == 3) face_norm[iface][1] = -1.0;
    } else if (domain->dimension == 3) {
      if (iface < 2) face_area[iface] = domain->prd[1]*domain->prd[2];
      else if (iface < 4) face_area[iface] = domain->prd[0]*domain->prd[2];
      else if (iface < 6) face_area[iface] = domain->prd[0]*domain->prd[1];
      if (iface == 0) face_norm[iface][0] = 1.0;
      else if (iface == 1) face_norm[iface][0] = -1.0;
      else if (iface == 2) face_norm[iface][1] = 1.0;
      else if (iface == 3) face_norm[iface][1] = -1.0;
      else if (iface == 4) face_norm[iface][2] = 1.0;
      else if (iface == 5) face_norm[iface][2] = -1.0;
    }
    face_weight[iface] = 1.0;
  }

  // set ptrs used by react() and PS_react() to per-face data structs

  species_state = face_species_state;
  total_state = face_total_state;
  area = face_area;
  weight = face_weight;
  species_delta = face_species_delta;
}

/* ----------------------------------------------------------------------
   create and initialize custom per-surf-element data structs
   this must be done exactly once even if multiple SRA instances are used
------------------------------------------------------------------------- */

void SurfReactAdsorb::create_per_surf_state()
{
  // SRA instance with first_owner = 1 allocates the custom per-surf data
  // add_custom() intializes all state data to zero
  // custom indices are for any type of persurf vec or array
  // direct indices are for specific types of vecs or arrays
  // cannot perform add_custom() for tau now, b/c nactive_ps is not yet known
  //   do this later in init

  if (surf->find_custom((char *) "nstick") < 0) {
    first_owner = 1;
    nstick_species_custom = surf->add_custom((char *) "nstick_species",
                                             INT,nspecies_surf);
    nstick_total_custom = surf->add_custom((char *) "nstick_total",INT,0);
    area_custom = surf->add_custom((char *) "area",DOUBLE,0);
    weight_custom = surf->add_custom((char *) "weight",DOUBLE,0);
    if (psflag) tau_custom = -1;

  } else {
    first_owner = 0;
    nstick_species_custom = surf->find_custom((char *) "nstick_species");
    nstick_total_custom = surf->find_custom((char *) "nstick_total");
    area_custom = surf->find_custom((char *) "area");
    weight_custom = surf->find_custom((char *) "weight");
    if (psflag) tau_custom = surf->find_custom((char *) "tau");
  }

  int nstick_species_direct = surf->ewhich[nstick_species_custom];
  int nstick_total_direct = surf->ewhich[nstick_total_custom];
  int area_direct = surf->ewhich[area_custom];
  int weight_direct = surf->ewhich[weight_custom];

  surf_species_state = surf->eiarray[nstick_species_direct];
  surf_total_state = surf->eivec[nstick_total_direct];
  surf_area = surf->edvec[area_direct];
  surf_weight = surf->edvec[weight_direct];

  // local delta storage
  // initialize surf_species_delta for first time

  int nlocal = surf->nlocal;
  memory->create(surf_species_delta,nlocal,nspecies_surf,
                 "react/adsorb:surf_species_delta");

  for (int isurf = 0; isurf < nlocal; isurf++)
    for (int isp = 0; isp < nspecies_surf; isp++)
      surf_species_delta[isurf][isp] = 0;

  // set ptrs used by react() and PS_react() to per-surf data structs

  species_state = surf_species_state;
  total_state = surf_total_state;
  area = surf_area;
  weight = surf_weight;
  species_delta = surf_species_delta;

  // allocate data structs for periodic sync

  int nown = surf->nown;

  memory->create(mark,nlocal,"react/adsorb:mark");
  memory->create(tally2surf,nlocal,"react/adsorb:tally2surf");
  memory->create(intally,nlocal,nspecies_surf,"react/adsorb:intally");
  memory->create(outtally,nlocal,nspecies_surf,"react/adsorb:outtally");
  memory->create(outcollate,nown,nspecies_surf,"react/adsorb:outcollate");

  maxtally = 0;
  incollate = NULL;

  // clear mark vector

  for (int isurf = 0; isurf < nlocal; isurf++) mark[isurf] = 0;
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorb::init()
{
  SurfReact::init();

  // this_index = index of this instance of surf react/adsorb
  //              in Surf list of reaction models

  this_index = surf->find_react(id);

  // initialize GS and PS models
  // for PS, this sets nactive_ps

  if (gsflag) init_reactions_gs();
  if (psflag) init_reactions_ps();

  // initialze tau only for PS models

  if (psflag) {
    if (mode == FACE) {
      nface = domain->dimension*2;
      memory->create(face_tau,nface,nactive_ps,"face_tau");
      for (int iface = 0; iface < nface; iface++) {
        for (int isp = 0; isp < nactive_ps; isp++)
          face_tau[iface][isp] = 0;
      }
      tau = face_tau;

    } else if (mode == SURF) {
      if (tau_custom == -1)
        tau_custom = surf->add_custom((char *) "tau",DOUBLE,nactive_ps);
      int itau = surf->ewhich[tau_custom];
      tau = surf->edarray[itau];
    }
  }

  // NOTE: should check that surf count has not changed since constructor
  //       b/c have lots of internal surf arrays
  //       else wait to allocate them until 1st init, then check
  //       count has not changed on subsequent init()

  // one-time initialize of custom per-surf attributes for area and weight
  // counts were initialized to zero by Surf class when custom vecs/arrays added
  // only set for surf elements assigned to this command ID
  // b/c there may be multiple instances of this command, each for different surfs
  // cannot do until now b/c surf->sr[isr] not set until run is performed

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nlocal = surf->nlocal;

  int isr;

  if (domain->dimension == 2) {
    for (int isurf = 0; isurf < nlocal; isurf++) {
      isr = lines[isurf].isr;
      if (surf->sr[isr] != this) return;
      area[isurf] = surf->line_size(&lines[isurf]);
      weight[isurf] = 1.0;
      if (psflag)
        for (int isp = 0; isp < nactive_ps; isp++) tau[isurf][isp] = 0.0;
    }
  } else {
    double tmp;
    for (int isurf = 0; isurf < nlocal; isurf++) {
      isr = tris[isurf].isr;
      if (surf->sr[isr] != this) return;
      area[isurf] = surf->tri_size(&tris[isurf],tmp);
      weight[isurf] = 1.0;
      if (psflag)
        for (int isp = 0; isp < nactive_ps; isp++) tau[isurf][isp] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   GS = gas/surface chemistry reactions
   select surface reaction to perform for particle with ptr IP on surface
   isurf < 0 for faces indexed from -1 to -6
   isurf >= 0 for surf element indexed from 0 to Nsurf-1
   norm = unit normal vector to face or line or tri
   ip = incident particle
   return 0 = no reaction
   return 1 = destroy reaction
   return 2 = create reaction
   if create, add particle and return ptr JP
------------------------------------------------------------------------- */

int SurfReactAdsorb::react(Particle::OnePart *&ip, int isurf, double *norm,
                           Particle::OnePart *&jp, int &velreset)
{
  // just return if GS model not defined

  if (!gsflag) return 0;

  // error checks

  if (isurf < 0 && mode == SURF)
    error->one(FLERR,"Surf_react adsorb surf used with box faces");
  if (isurf >= 0 && mode == FACE)
    error->one(FLERR,"Surf_react adsorb face used with surface elements");

  // convert face index from negative value to 0 to 5 inclusive

  if (mode == FACE) isurf = -(isurf+1);

  // n = # of possible reactions for particle IP

  int *list = reactions_gs[ip->ispecies].list;
  int n = reactions_gs[ip->ispecies].n;
  if (n == 0) return 0;

  double fnum = update->fnum;
  long int maxstick = ceil(max_cover*area[isurf] / (fnum*weight[isurf]));
  double factor = fnum * weight[isurf] / area[isurf];
  double ms_inv = factor / max_cover;

  // loop over possible reactions for this species

  Particle::Species *species = particle->species;

  OneReaction_GS *r;
  double prob_value[n], sum_prob = 0.0;
  double scatter_prob = 0.0, correction = 1.0;
  //int check_ads = 0, ads_index = -1;

  int coeff_val = 1;

  for (int i = 0; i < n; i++) {
    r = &rlist_gs[list[i]];

    if (r->style == ARRHENIUS) coeff_val = 3;

    switch (r->type) {
    case DISSOCIATION:
      {
        //prob_value[i] = r->coeff[0];
        prob_value[i] = r->k_react;
        break;
      }

    case EXCHANGE:
      {
        //prob_value[i] = r->coeff[0];
        prob_value[i] = r->k_react;
        break;
      }

    case RECOMBINATION:
      {
        //prob_value[i] = r->coeff[0];
        prob_value[i] = r->k_react;
        break;
      }

    case AA:
      {
        //check_ads = 1;
        //ads_index = i;
        double surf_cover = total_state[isurf] * ms_inv;
        double S_theta = 0.0;

        if (r->kisliuk_flag)
        {
          double K_ads = r->kisliuk_coeff[0] * pow(twall,r->kisliuk_coeff[1]) *
          exp(-r->kisliuk_coeff[2]/twall);
          if (surf_cover < 1)
            S_theta = pow((1 - surf_cover) /
                          (1 - surf_cover +
                           K_ads*surf_cover),r->coeff[coeff_val]);
        }
        else
        {
          S_theta = pow((1-surf_cover),r->coeff[coeff_val]);
        }
        //scatter_prob = 1 - r->k_react*S_theta;
        //prob_value[i] = 1.0;
        prob_value[i] = r->k_react*S_theta;
        break;
      }

    case DA:
      {
        double surf_cover = total_state[isurf] * ms_inv;
        double S_theta = 0.0;

        if (r->kisliuk_flag) {
          double K_ads = r->kisliuk_coeff[0] * pow(twall,r->kisliuk_coeff[1]) *
            exp(-r->kisliuk_coeff[2]/twall);
          if (surf_cover < 1)
            S_theta = pow((1 - surf_cover) /
                          (1 - surf_cover +
                           K_ads*surf_cover),r->coeff[coeff_val]);
        } else {
          S_theta = pow((1-surf_cover),r->coeff[coeff_val]);
        }

        prob_value[i] = r->k_react*S_theta;

        /*
        if (r->state_products[1][0] == 's') {
          double K_ads2 = r->coeff[6] * pow(twall,r->coeff[7]) *
            exp(-r->coeff[8]/twall);
          double S_ratio2 = (1 - surf_cover) /
               (1 - surf_cover + K_ads*surf_cover);
          prob_value[i] *= (S_ratio2);
        }
        */

        break;
      }

    case LH1:
      {
        double surf_cover = total_state[isurf] * ms_inv;
        double S_theta = 0.0;

        if (r->kisliuk_flag) {
          double K_ads = r->kisliuk_coeff[0] * pow(twall,r->kisliuk_coeff[1]) *
            exp(-r->kisliuk_coeff[2]/twall);
          if (surf_cover < 1)
            S_theta = pow((1 - surf_cover) /
                          (1 - surf_cover +
                           K_ads*surf_cover),r->coeff[coeff_val]);
        } else {
          S_theta = pow((1-surf_cover),r->coeff[coeff_val]);
        }

        prob_value[i] = r->k_react*S_theta;
        break;
      }

    case LH3:
      {
        double surf_cover = total_state[isurf] * ms_inv;
        double S_theta = 0.0;

        if (r->kisliuk_flag) {
          double K_ads = r->kisliuk_coeff[0] * pow(twall,r->kisliuk_coeff[1]) *
            exp(-r->kisliuk_coeff[2]/twall);
          if (surf_cover < 1)
            S_theta = pow((1 - surf_cover) /
                          (1 - surf_cover +
                           K_ads*surf_cover),r->coeff[coeff_val]);
        } else {
          S_theta = pow((1-surf_cover),r->coeff[coeff_val]);
        }

        prob_value[i] = r->k_react*S_theta;
        break;
      }

    case CD:
      {
        double surf_cover = total_state[isurf] * ms_inv;
        double S_theta = 0.0;

        if (r->kisliuk_flag) {
          double K_ads = r->kisliuk_coeff[0] * pow(twall,r->kisliuk_coeff[1]) *
            exp(-r->kisliuk_coeff[2]/twall);
          if (surf_cover < 1)
            S_theta = pow((1 - surf_cover) /
                          (1 - surf_cover +
                           K_ads*surf_cover),r->coeff[coeff_val]);
        } else {
          S_theta = pow((1-surf_cover),r->coeff[coeff_val]);
        }

        prob_value[i] = r->k_react*S_theta;
        break;
      }

    case ER:
      {
        double dot = MathExtra::dot3(ip->v,norm);
        dot = 2.0;

        if (r->nreactant == 1) {
          prob_value[i] = 2.0 * r->k_react *
            (maxstick - total_state[isurf]) * ms_inv / fabs(dot);
        } else {
          prob_value[i] = 2.0 * r->k_react / fabs(dot);
        }
        break;

      }

    case CI:
      {
        prob_value[i] = r->k_react;
        if (r->energy_flag) {
          double *v = ip->v;
          double dot = MathExtra::dot3(v,norm);
          double vmag_sq = MathExtra::lensq3(v);
          double E_i = 0.5 * species[ip->ispecies].mass * vmag_sq;
          double cos_theta = abs(dot) / sqrt(vmag_sq);
          prob_value[i] *= pow(E_i,r->energy_coeff[0]) *
          pow(cos_theta,r->energy_coeff[1]);
        }
        break;
      }

    /*
    case CI2:
      {
        double *v = ip->v;
        double dot = MathExtra::dot3(v,norm);
        double vmag_sq = MathExtra::lensq3(v);
        double E_i = 0.5 * species[ip->ispecies].mass * vmag_sq;
        double cos_theta = dot/sqrt(vmag_sq);
        prob_value[i] = r->k_react * pow(E_i,r->coeff[3]) *
          pow(cos_theta,r->coeff[4]);
        break;
      }
      */
    }

    for (int j = 1; j < r->nreactant; j++) {
      if (r->state_reactants[j][0] == 's') {
        if (r->part_reactants[j] == 0) {
          prob_value[i] *=
            stoich_pow(total_state[isurf],
                       r->stoich_reactants[j]) *
            pow(ms_inv,r->stoich_reactants[j]);
        } else {
          prob_value[i] *=
            stoich_pow(species_state[isurf][r->reactants_ad_index[j]],
                       r->stoich_reactants[j]) *
            pow(ms_inv,r->stoich_reactants[j]);
        }
      }
    }

    sum_prob += prob_value[i];
  }

  // NOTE: scatter_prob is always zero?
  /*
  if (check_ads) {
    if (sum_prob > (1-scatter_prob))
      correction = (1-scatter_prob) / sum_prob;
  } else if (sum_prob > 1.0) correction = 1.0/sum_prob;
  */

  if (sum_prob > 1.0) correction = 1.0/sum_prob;
  else scatter_prob = 1.0 - sum_prob;

  //if (sum_prob > (1.0-scatter_prob))
  //  correction = (1.0-scatter_prob) / sum_prob;

  // probablity to compare to reaction probability

  double react_prob = scatter_prob;
  double random_prob = random->uniform();

  if (react_prob > random_prob) return 0;
  else {
    // NOTE: at this point is it guaranteed a reaction will take place?

    if (mode == SURF) mark[isurf] = 1;

    // NOTE: why summing to previous react_prob?

    for (int i = 0; i < n; i++) {
      r = &rlist_gs[list[i]];
      react_prob += prob_value[i] * correction;
      if (react_prob <= random_prob) continue;

      // perform the reaction and return
      // if dissociation or CI2 performs a realloc:
      //   make copy of x,v, then repoint ip to new particles data struct

      nsingle++;
      tally_single[list[i]]++;

      for (int j = 0; j < r->nreactant; j++) {
        if (r->part_reactants[j] == 1) {
          switch(r->state_reactants[j][0]) {
          case 's':
            {
              species_delta[isurf][r->reactants_ad_index[j]] -=
                r->stoich_reactants[j];
              break;
            }

          case 'g': {}
          case 'b': {}
          }
        }
      }

      for (int j = 0; j < r->nproduct; j++) {
        if (r->part_products[j] == 1) {
          switch(r->state_products[j][0]) {
          case 's':
            {
              species_delta[isurf][r->products_ad_index[j]] +=
                r->stoich_products[j];
              break;
            }

          case 'g': {}
          case 'b': {}
          }
        }
      }

      // for each reaction, post-reaction velocities must be set
      // if NOMODEL:
      //   leave velreset = 0 (value passed by caller)
      //   SC instance associated with the surf/face sets vels after return
      // else:
      //   set velreset = 1
      //   SC style created when GS file was read does velocity reset

      switch (r->type) {

      case DISSOCIATION:
        {
          double x[3],v[3];
          ip->ispecies = r->products[0];
          int id = MAXSMALLINT*random->uniform();
          memcpy(x,ip->x,3*sizeof(double));
          memcpy(v,ip->v,3*sizeof(double));
          Particle::OnePart *particles = particle->particles;

          int jp_species;
          if (r->stoich_products[0] == 2) jp_species = r->products[0];
          else jp_species = r->products[1];

          int reallocflag =
            particle->add_particle(id,jp_species,ip->icell,x,v,0.0,0.0);
          if (reallocflag) ip = particle->particles + (ip - particles);
          jp = &particle->particles[particle->nlocal-1];
          return (list[i] + 1);
         break;
        }

      case EXCHANGE:
        {
          ip->ispecies = r->products[0];
          return (list[i] + 1);
          break;
        }

      case RECOMBINATION:
        {
          ip = NULL;
          return (list[i] + 1);
          break;
        }

      case AA:
        {
          ip = NULL;
          return (list[i] + 1);
          break;
        }

      case DA:
        {
          if (r->nprod_g == 0) ip = NULL;
          else {
            int n = 1;
            for (int j = 1; j < r->nproduct; j++) {
              if (r->state_products[j][0] == 'g') {
                if (n == 1) {
                  n++;
                  ip->ispecies = r->products[j];
                  if (r->cmodel_ip != NOMODEL)
                    cmodels[r->cmodel_ip]->
                      wrapper(ip,norm,r->cmodel_ip_flags,r->cmodel_ip_coeffs);

                  if (r->stoich_products[j] == 2) {
                    double x[3],v[3];
                    memcpy(x,ip->x,3*sizeof(double));
                    memcpy(v,ip->v,3*sizeof(double));

                    int id = MAXSMALLINT*random->uniform();
                    Particle::OnePart *particles = particle->particles;

                    int reallocflag =
                      particle->add_particle(id,r->products[j],ip->icell,
                                             x,v,0.0,0.0);
                    if (reallocflag) ip = particle->particles + (ip - particles);
                    jp = &particle->particles[particle->nlocal-1];

                    if (r->cmodel_ip != NOMODEL)
                      cmodels[r->cmodel_ip]->
                        wrapper(jp,norm,r->cmodel_ip_flags,r->cmodel_ip_coeffs);
                  }

                } else {
                    double x[3],v[3];
                    memcpy(x,ip->x,3*sizeof(double));
                    memcpy(v,ip->v,3*sizeof(double));

                    int id = MAXSMALLINT*random->uniform();
                    Particle::OnePart *particles = particle->particles;

                    int reallocflag =
                      particle->add_particle(id,r->products[j],ip->icell,
                                             x,v,0.0,0.0);
                    if (reallocflag) ip = particle->particles + (ip - particles);
                    jp = &particle->particles[particle->nlocal-1];

                    if (r->cmodel_jp != NOMODEL)
                      cmodels[r->cmodel_jp]->
                        wrapper(jp,norm,r->cmodel_jp_flags,r->cmodel_jp_coeffs);
                }
              }
            }
          }

          if (r->cmodel_ip != NOMODEL) velreset = 1;
          return (list[i] + 1);
          break;
        }

      case LH1:
        {
          ip->ispecies = r->products[0];
          if (r->cmodel_ip != NOMODEL)
            cmodels[r->cmodel_ip]->
              wrapper(ip,norm,r->cmodel_ip_flags,r->cmodel_ip_coeffs);

          if (r->cmodel_ip != NOMODEL) velreset = 1;
          return (list[i] + 1);
          break;
        }

      case LH3:
        {
          ip = NULL;
          return (list[i] + 1);
          break;
        }

      case CD:
        {
          ip = NULL;
          return (list[i] + 1);
          break;
        }

      case ER:
        {
          ip->ispecies = r->products[0];
          if (r->cmodel_ip != NOMODEL)
            cmodels[r->cmodel_ip]->
              wrapper(ip,norm,r->cmodel_ip_flags,r->cmodel_ip_coeffs);

          if (r->cmodel_ip != NOMODEL) velreset = 1;
          return (list[i] + 1);
          break;
        }

      case CI:
        {
          ip->ispecies = r->products[0];

          if (r->cmodel_ip != NOMODEL)
            cmodels[r->cmodel_ip]->
              wrapper(ip,norm,r->cmodel_ip_flags,r->cmodel_ip_coeffs);

          if (r->nprod_g_tot == 2) {
            double x[3],v[3];
            memcpy(x,ip->x,3*sizeof(double));
            memcpy(v,ip->v,3*sizeof(double));

            int id = MAXSMALLINT*random->uniform();
            Particle::OnePart *particles = particle->particles;

            if (r->stoich_products[0] == 2) {
              int reallocflag =
                particle->add_particle(id,r->products[0],ip->icell,x,v,0.0,0.0);
              if (reallocflag) ip = particle->particles + (ip - particles);
              jp = &particle->particles[particle->nlocal-1];

              if (r->cmodel_ip != NOMODEL)
                cmodels[r->cmodel_ip]->wrapper(jp,norm,r->cmodel_ip_flags,
                                               r->cmodel_ip_coeffs);
            } else {
              int reallocflag =
                particle->add_particle(id,r->products[1],ip->icell,x,v,0.0,0.0);
              if (reallocflag) ip = particle->particles + (ip - particles);
              jp = &particle->particles[particle->nlocal-1];

              if (r->cmodel_jp != NOMODEL)
                cmodels[r->cmodel_jp]->wrapper(jp,norm,r->cmodel_jp_flags,
                                               r->cmodel_jp_coeffs);
            }
          }

          if (r->cmodel_ip != NOMODEL) velreset = 1;
          return (list[i] + 1);
          break;
        }
      }
    }
  }

  // no reaction

  return 0;
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorb::tally_update()
{
  // tally only the gas phase reactions

  int nsingle_gs = nsingle;
  ntotal += nsingle;
  for (int i = 0; i < nlist_gs; i++) tally_total[i] += tally_single[i];

  // sync surface state across all procs
  // only done once every Nsync steps

  if (update->ntimestep % nsync) return;

  // perform on-surface chemistry for PS
  // will insert new particles desorbing from faces/surfs as needed
  // first sync gas/surf chem changes to surf states since last sync

  if (psflag) {
    if (mode == FACE) update_state_face();
    else if (mode == SURF) update_state_surf();

    PS_chemistry();
  }

  // update the state of all faces or surf elements
  // if no PS chemistry, syncs gas/surf chem changes to surf states
  //   since last sync
  // if yes PS chemistry, syncs PS chem changes to surf states

  if (mode == FACE) update_state_face();
  else if (mode == SURF) update_state_surf();

  // tally only the surf phase reactions

  ntotal += nsingle - nsingle_gs;
  for (int i = nlist_gs; i < nlist; i++) tally_total[i] += tally_single[i];
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorb::PS_chemistry()
{
  // zero the mypart vector of particles this proc is adding

  npart = 0;

  // extract list of Computes which tally reactions from Update
  // used by PS_react() to tally on-surface reactions

  if (mode == FACE) {
    ncompute_tally = update->nboundary_tally;
    clist_active = update->blist_active;
  } else if (mode == SURF) {
    ncompute_tally = update->nsurf_tally;
    clist_active = update->slist_active;
  }

  // for box faces: a single proc updates all faces
  // for surf elements: each proc updates every Pth surf it owns
  // only call PS_react() if this instance of surf react/adsorb matches
  //   the reaction model assigned to surf or face

  if (mode == FACE) {
    if (me == 0) {
      for (int iface = 0; iface < nface; iface++)
        if (domain->surf_react[iface] == this_index)
          PS_react(PSFACE,iface,face_norm[iface]);
    }

  } else if (mode == SURF) {
    int nsurf = surf->nsurf;
    if (domain->dimension == 2) {
      for (int m = me; m < nsurf; m += nprocs)
        if (surf->lines[m].isr == this_index)
          PS_react(PSLINE,m,surf->lines[m].norm);
    } else {
      for (int m = me; m < nsurf; m += nprocs)
        if (surf->tris[m].isr == this_index)
          PS_react(PSTRI,m,surf->tris[m].norm);
    }
  }

  // allpart = nall-length vector of particles all procs are adding
  // accumulate via Allgatherv

  int nall;
  MPI_Allreduce(&npart,&nall,1,MPI_INT,MPI_SUM,world);

  if (nall > maxallpart) {
    while (maxallpart < nall) maxallpart += DELTA_PART;
    memory->sfree(allpart);
    allpart = (AddParticle *)
      memory->smalloc(maxallpart*sizeof(AddParticle),"sr_adsorb:allpart");
  }

  int nsend = npart*sizeof(AddParticle);
  MPI_Allgather(&nsend,1,MPI_INT,recvcounts,1,MPI_INT,world);
  displs[0] = 0;
  for (int i = 1; i < nprocs; i++) displs[i] = displs[i-1] + recvcounts[i-1];

  MPI_Allgatherv(mypart,nsend,MPI_CHAR,allpart,recvcounts,displs,MPI_CHAR,world);

  // loop over all particles
  // check if inside a child cell I own via id_find_child()
  //   if child cells is a split cell, find subcell via update->split()
  // if not, skip the particle, another proc will add it
  // if yes, add it to particle list using values in allpart
  // dtremain must be set separately
  // grid->hash must be filled to use grid->id_find_child()
  // NOTE: does this need logic for handling split cells ?

  Grid::ChildCell *cells = grid->cells;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  int dimension = domain->dimension;

  int icell;
  double *x;
  AddParticle *p;

  for (int i = 0; i < nall; i++) {
    p = &allpart[i];
    x = p->x;

    icell = grid->id_find_child(0,0,boxlo,boxhi,x);
    if (icell < 0) continue;
    if (icell >= grid->nlocal) continue;
    if (cells[icell].nsplit > 1) {
      if (dimension == 3) icell = update->split3d(icell,x);
      else icell = update->split2d(icell,x);
    }

    particle->add_particle(p->id,p->ispecies,icell,p->x,p->v,p->erot,p->evib);
    particle->particles[particle->nlocal-1].dtremain = p->dtremain;
  }
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorb::update_state_face()
{
  int i,j;

  // sum perspecies deltas across all procs

  MPI_Allreduce(&species_delta[0][0],&face_sum_delta[0][0],
                nface*nspecies_surf,MPI_INT,MPI_SUM,world);

  // new perspecies state = old perspecies state + summed delta
  // insure no counts < 0
  // new total state = sum of new perspecies state over species
  // re-initialize species_delta, i.e. face_species_delta

  for (i = 0; i < nface; i++) {
    total_state[i] = 0;
    for (j = 0; j < nspecies_surf; j++) {
      species_state[i][j] += face_sum_delta[i][j];
      species_state[i][j] = MAX(0,species_state[i][j]);
      total_state[i] += species_state[i][j];
      species_delta[i][j] = 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorb::update_state_surf()
{
  int i,j,m,isr;

  // incollate = array of deltas for surfs I marked
  // ntally = # of surfs I marked = # of rows in incollate
  // tally2surf = global surf index (1 to Nsurf) for each row of array
  // re-initialize non-zero species_delta values, i.e. surf_species_delta

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nlocal = surf->nlocal;

  int ntally = 0;

  if (domain->dimension == 2) {
    for (int i = 0; i < nlocal; i++) {
      isr = lines[i].isr;
      if (surf->sr[isr] != this) continue;
      if (!mark[i]) continue;

      if (ntally == maxtally) {
        maxtally += DELTA_TALLY;
        memory->grow(tally2surf,maxtally,"react/adsorb:tally2surf");
        memory->grow(incollate,maxtally,nspecies_surf,
                     "react/adsorb:incollate");
      }

      tally2surf[ntally] = i+1;
      for (j = 0; j < nspecies_surf; j++) {
        incollate[ntally][j] = species_delta[i][j];
        species_delta[i][j] = 0;
      }
      ntally++;
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      isr = tris[i].isr;
      if (surf->sr[isr] != this) continue;
      if (!mark[i]) continue;

      if (ntally == maxtally) {
        maxtally += DELTA_TALLY;
        memory->grow(tally2surf,maxtally,"react/adsorb:tally2surf");
        memory->grow(incollate,maxtally,nspecies_surf,
                     "react/adsorb:incollate");
      }

      tally2surf[ntally] = i+1;
      for (j = 0; j < nspecies_surf; j++) {
        incollate[ntally][j] = species_delta[i][j];
        species_delta[i][j] = 0;
      }
      ntally++;
    }
  }

  // perform the collate
  // outcollate = values only for owned surfs

  surf->collate_array(ntally,nspecies_surf,tally2surf,incollate,outcollate);

  // MPI_Allreduce of outcollate so all procs know new state of all surfs
  // NOTE: inefficient, but necessary for non-distributed surfs?
  //       since each proc may perform a collision with any surf

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < nspecies_surf; j++)
      intally[i][j] = 0;

  m = 0;
  for (i = me; i < nlocal; i += nprocs) {
    for (j = 0; j < nspecies_surf; j++)
      intally[i][j] = static_cast<int> (outcollate[m][j]);
    m++;
  }

  MPI_Allreduce(&intally[0][0],&outtally[0][0],nlocal*nspecies_surf,
                MPI_INT,MPI_SUM,world);

  // new perspecies state = old perspecies state + summed delta (outtally)
  // insure no counts < 0
  // new total state = sum of new perspecies state over species

  for (i = 0; i < nlocal; i++) {
    total_state[i] = 0;
    for (j = 0; j < nspecies_surf; j++) {
      species_state[i][j] += outtally[i][j];
      species_state[i][j] = MAX(0,species_state[i][j]);
      total_state[i] += species_state[i][j];
    }
  }

  // clear mark vector

  for (int isurf = 0; isurf < nlocal; isurf++) mark[isurf] = 0;
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorb::init_reactions_gs()
{
  // convert species IDs to species indices
  // flag reactions as active/inactive depending on whether all species exist

  for (int m = 0; m < nlist_gs; m++) {
    OneReaction_GS *r = &rlist_gs[m];
    r->active = 1;
    for (int i = 0; i < r->nreactant; i++) {
      r->reactants[i] = particle->find_species(r->id_reactants[i]);
      if (r->reactants[i] < 0) {
        r->active = 0;
        break;
      }
      if (r->state_reactants[i][0] == 's') {
        r->reactants_ad_index[i] = find_surf_species(r->id_reactants[i]);
        if (r->reactants_ad_index[i] < 0) {
          r->active = 0;
          break;
        }
      }
      else r->reactants_ad_index[i] = -1;
    }

    for (int i = 0; i < r->nproduct; i++) {
      r->products[i] = particle->find_species(r->id_products[i]);
      if (r->products[i] < 0) {
        r->active = 0;
        break;
      }
      if (r->state_products[i][0] == 's') {
        r->products_ad_index[i] = find_surf_species(r->id_products[i]);
        if (r->products_ad_index[i] < 0) {
          r->active = 0;
          break;}
      }
      else r->products_ad_index[i] = -1;
    }
  }

  // count possible reactions for each species

  memory->destroy(reactions_gs);
  int nspecies = particle->nspecies;
  reactions_gs = memory->create(reactions_gs,nspecies,
                             "surf_adsorb:reactions");

  for (int i = 0; i < nspecies; i++) reactions_gs[i].n = 0;

  int n = 0;
  for (int m = 0; m < nlist_gs; m++) {
    OneReaction_GS *r = &rlist_gs[m];
    if (!r->active) continue;
    int i = r->reactants[0];
    reactions_gs[i].n++;
    n++;
  }

  // allocate indices_gs = entire list of reactions for all I species

  memory->destroy(indices_gs);
  memory->create(indices_gs,n,"surf_adsorb:indices_gs");

  // reactions_gs[i].list = offset into full indices_gs vector

  int offset = 0;
  for (int i = 0; i < nspecies; i++) {
    reactions_gs[i].list = &indices_gs[offset];
    offset += reactions_gs[i].n;
  }

  // reactions_gs[i].list = indices_gs of possible reactions for each species

  for (int i = 0; i < nspecies; i++) reactions_gs[i].n = 0;

  for (int m = 0; m < nlist_gs; m++) {
    OneReaction_GS *r = &rlist_gs[m];
    if (!r->active) continue;
    int i = r->reactants[0];
    reactions_gs[i].list[reactions_gs[i].n++] = m;
  }

  // check that summed reaction probabilities for each species <= 1.0

//  double sum;
//  for (int i = 0; i < nspecies; i++) {
//    sum = 0.0;
//    for (int j = 0; j < reactions_gs[i].n; j++)
//      sum += rlist_gs[reactions_gs[i].list[j]].coeff[0];
//    if (sum > 1.0)
//      error->all(FLERR,"Surface reaction probability for a species > 1.0");
//  }
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorb::readfile_gs(char *fname)
{
  int n,n1,n2,eof;
  char line1[MAXLINE],line2[MAXLINE];
  char copy1[MAXLINE],copy2[MAXLINE],copy3[MAXLINE],copy4[MAXLINE];
  char *word;
  OneReaction_GS *r;

  // proc 0 opens file

  if (me == 0) {
    fp = fopen(fname,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open reaction file %s",fname);
      error->one(FLERR,str);
    }
  }

  // read reactions one at a time and store their info in rlist_gs

  while (1) {
    if (me == 0) eof = readone(line1,line2,n1,n2);
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;

    MPI_Bcast(&n1,1,MPI_INT,0,world);
    MPI_Bcast(&n2,1,MPI_INT,0,world);
    MPI_Bcast(line1,n1,MPI_CHAR,0,world);
    MPI_Bcast(line2,n2,MPI_CHAR,0,world);

    if (nlist_gs == maxlist_gs) {
      maxlist_gs += DELTALIST;
      rlist_gs = (OneReaction_GS *)
        memory->srealloc(rlist_gs,maxlist_gs*sizeof(OneReaction_GS),
                         "surf_react/GS:rlist_gs");

      for (int i = nlist_gs; i < maxlist_gs; i++) {
        r = &rlist_gs[i];
        r->nreactant = r->nproduct = 0;
        r->nprod_g = r->nprod_g_tot = 0;
        r->id = NULL;
        r->id_reactants = new char*[MAXREACTANT_GS];
        r->id_products = new char*[MAXPRODUCT_GS];
        r->state_reactants = new char*[MAXREACTANT_GS];
        r->state_products = new char*[MAXPRODUCT_GS];
        r->part_reactants = new int[MAXREACTANT_GS];
        r->part_products = new int[MAXPRODUCT_GS];
        r->stoich_reactants = new int[MAXREACTANT_GS];
        r->stoich_products = new int[MAXPRODUCT_GS];
        r->reactants = new int[MAXREACTANT_GS];
        r->products = new int[MAXPRODUCT_GS];
        r->reactants_ad_index = new int[MAXREACTANT_GS];
        r->products_ad_index = new int[MAXPRODUCT_GS];
        r->coeff = new double[MAXCOEFF_GS];
        r->cmodel_ip = NOMODEL;
        r->cmodel_ip_flags = NULL;
        r->cmodel_ip_coeffs = NULL;
        r->cmodel_jp = NOMODEL;
        r->cmodel_jp_flags = NULL;
        r->cmodel_jp_coeffs = NULL;
      }
    }

    strcpy(copy1,line1);
    strcpy(copy2,line2);

    r = &rlist_gs[nlist_gs];

    int side = 0;
    int species = 1;
    int start = 0;

    // process 1st line of reaction

    n = strlen(line1) - 1;
    r->id = new char[n+1];
    strncpy(r->id,line1,n);
    r->id[n] = '\0';

    word = strtok(line1," \t\n");

    while (1) {
      if (!word) {
        if (side == 0) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        if (species) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        break;
      }

      if (species) {
        species = 0;
        if (side == 0) {
          if (r->nreactant == MAXREACTANT_GS) {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Too many reactants in a reaction formula");
            }
          n = strlen(word) + 1;
          start = 0;
          r->part_reactants[r->nreactant] = 1;
          r->stoich_reactants[r->nreactant] = 1;
          if (word[n-2] != ')') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Specify the state of the reactants");
          }
          if (word[n-3] == 'c') {
            r->part_reactants[r->nreactant] = 0;
            n--;
          }
          /*
          if (r->nreactant == 0 && word[n-3] != 'g') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"The first reactant must be gas phase");
          }
          */
          if (r->nreactant != 0 && word[n-3] == 'g') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Only one gas phase reactant can be present");
          }
          if (r->part_reactants[r->nreactant] == 0 && word[n-3] == 'g') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Gas phase reactants cannot be catalytic");
          }
          if (r->part_reactants[r->nreactant] == 0 && word[n-3] == 'b') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Bulk phase reactants cannot be catalytic");
          }
          if (isdigit(word[0])) {
            r->stoich_reactants[r->nreactant] = atoi(word);
            start=1;
          }
          r->id_reactants[r->nreactant] = new char[n-start-3];
          strncpy(r->id_reactants[r->nreactant],&(word[start]),n-start-4);
          r->id_reactants[r->nreactant][n-start-4] = '\0';
          r->state_reactants[r->nreactant] = new char[1]();
          strncpy(r->state_reactants[r->nreactant],1+strstr(word,"("),1);
          r->nreactant++;

        } else {
          if (r->nproduct == MAXPRODUCT_GS) {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Too many products in a reaction formula");
          }

          n = strlen(word) + 1;
          start = 0;
          r->part_products[r->nproduct] = 1;
          r->stoich_products[r->nproduct] = 1;
          if (word[n-2] != ')') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Specify the state of the products");
          }
          if (word[n-3] == 'c') {
            r->part_products[r->nproduct] = 0;
            n--;
          }
          if (r->part_products[r->nproduct] == 0 && word[n-3] == 'g') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Gas phase products cannot be catalytic");
          }
          if (r->part_products[r->nproduct] == 0 && word[n-3] == 'b') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Bulk phase products cannot be catalytic");
          }
          if (isdigit(word[0])) {
            r->stoich_products[r->nproduct] = atoi(word);
            start=1;
          }
          r->id_products[r->nproduct] = new char[n-start-3]();
          strncpy(r->id_products[r->nproduct],&(word[start]),n-start-4);
          r->id_products[r->nproduct][n-start-4] = '\0';
          r->state_products[r->nproduct] = new char[1]();
          strncpy(r->state_products[r->nproduct],1+strstr(word,"("),1);

          if(r->state_products[r->nproduct][0] == 'g')
            {
              r->nprod_g++;
              r->nprod_g_tot += r->stoich_products[r->nproduct];
            }

          r->nproduct++;
        }

      } else {
        species = 1;
        if (strcmp(word,"+") == 0) {
          word = strtok(NULL," \t\n");
          continue;
        }
        if (strcmp(word,"-->") != 0) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        side = 1;
      }

      word = strtok(NULL," \t\n");
    }

    // replace single NULL product with no products

    if (r->nproduct == 1 && strcmp(r->id_products[0],"NULL") == 0) {
      delete [] r->id_products[0];
      r->id_products[0] = NULL;
      r->nproduct = 0;
    }

    // process 2nd line of reaction

    word = strtok(line2," \t\n");
    if (!word) {
      print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction type in file");
    }

    r->ncoeff = 0;
    if (strcmp(word,"D") == 0 || strcmp(word,"d") == 0) {
      r->type = DISSOCIATION;
    } else if (strcmp(word,"E") == 0 || strcmp(word,"e") == 0) {
      r->type = EXCHANGE;
    } else if (strcmp(word,"R") == 0 || strcmp(word,"r") == 0) {
      r->type = RECOMBINATION;
    } else if (strcmp(word,"AA") == 0 || strcmp(word,"aa") == 0) {
      r->type = AA;
      r->ncoeff = 1;
    } else if (strcmp(word,"DA")==0 || strcmp(word,"da")==0) {
      r->type = DA;
      r->ncoeff = 1;
    } else if (strcmp(word,"LH1") == 0 || strcmp(word,"lh1") == 0) {
      r->type = LH1;
      r->ncoeff = 1;
    } else if (strcmp(word,"LH3") == 0 || strcmp(word,"lh3") == 0) {
      r->type = LH3;
      r->ncoeff = 1;
    } else if (strcmp(word,"CD") == 0 || strcmp(word,"cd") == 0) {
      r->type = CD;
      r->ncoeff = 1;
    } else if (strcmp(word,"ER") == 0 || strcmp(word,"er") == 0) {
      r->type = ER;
    } else if (strcmp(word,"CI") == 0 || strcmp(word,"ci") == 0) {
      r->type = CI;
    } else {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction type in file");
    }

    word = strtok(NULL," \t\n");
    if (!word) {
      print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction type in file");
    }

    if (word[0] == 'S' || word[0] == 's') {
      r->style = SIMPLE; r->ncoeff += 1;
    } else if (word[0] == 'A' || word[0] == 'a') {
      r->style = ARRHENIUS; r->ncoeff += 3;
    } else {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction type in file");
    }

    for (int i = 0; i < r->ncoeff; i++) {
      word = strtok(NULL," \t\n");
      if (!word) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction coefficients in file");
      }
      r->coeff[i] = input->numeric(FLERR,word);
    }

    r->kisliuk_flag = 0;
    r->energy_flag = 0;

    word = strtok(NULL," \t\n");
    while (word != NULL) {
      if (strcmp(word,"kisliuk") == 0) {
        r->kisliuk_flag = 1;
        for (int i = 0; i < 3; i++) {
          word = strtok(NULL," \t\n");
          if (!word) {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Invalid reaction coefficients in file");
          }
          r->kisliuk_coeff[i] = input->numeric(FLERR,word);
        }
      } else if (strcmp(word,"energy") == 0)
      {
        if (r->type != CI)
        {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Energy option can only be used to define "
                     "the reaction rate constant in CI reaction");
        }

        r->energy_flag = 1;
        for (int i = 0; i < 2; i++) {
          word = strtok(NULL," \t\n");
          if (!word) {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Invalid reaction coefficients in file");
          }
          r->energy_coeff[i] = input->numeric(FLERR,word);
        }
      } else {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction type in file");
      }
      word = strtok(NULL," \t\n");
    }

    /*
    word = strtok(NULL," \t\n");
    if (word) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Too many coefficients in a reaction formula");
    }
    */

    // ERROR CHECKS
    // check that reactant/product counts are consistent with type
    /*
    int n_reactant_exp, n_product_exp;

    if (r->type == DISSOCIATION) {
      n_reactant_exp = 1;
      n_product_exp = 2;
    }
    */

    if (r->state_reactants[0][0] != 'g') {
      print_reaction(copy1,copy2);
      error->all(FLERR,"The first reactant must be gas phase");
    }

    if (r->nprod_g_tot > 2) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Number of gas phase products cannot be greater than 2");
    }

    /*
    if (r->type == DISSOCIATION) {
      if (r->nreactant != 1 || r->nproduct != 2) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction type in file");
      }
    } else if (r->type == EXCHANGE) {
      if (r->nreactant != 1 || r->nproduct != 1) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction type in file");
      }
    } else if (r->type == RECOMBINATION) {
      if (r->nreactant != 1 || r->nproduct != 0) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction type in file");
      }
      */
    //} else

    switch (r->type) {

    case DISSOCIATION:
      {
        /*
          if (r->nreactant != 1 || r->nproduct != 2)
          {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction type in file");
          }
        */
        if (r->kisliuk_flag)
        {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Kisliuk option can only be used to define the reaction rate constant in adsorption mediated reaction such as AA, DA, LH1, LH3, and CD");
        }
        break;
      }

    case EXCHANGE:
      {
        /*
          if (r->nreactant != 1 || r->nproduct != 1)
          {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction type in file");
          }
        */
        if (r->kisliuk_flag)
        {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Kisliuk option can only be used to define the reaction rate constant in adsorption mediated reaction such as AA, DA, LH1, LH3, and CD");
        }
        break;
      }

    case RECOMBINATION:
      {
        /*
          if (r->nreactant != 1 || r->nproduct != 0)
          {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction type in file");
          }
        */
        if (r->kisliuk_flag)
        {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Kisliuk option can only be used to define the reaction rate constant in adsorption mediated reaction such as AA, DA, LH1, LH3, and CD");
        }
        break;
      }

    case AA:
      {
        if (r->state_products[0][0] != 's') {
          print_reaction(copy1,copy2);
          error->all(FLERR,
                     "First product must be surface phase in AA reaction");
        }
        break;
      }

    case DA:
      {
        if (r->state_products[0][0] == 'g') {
          print_reaction(copy1,copy2);
          error->all(FLERR,
                     "First product must be surface or bulk phase in DA reaction");
        }
        /*
        for (int i=2; i < r->nproduct; i++) {
          if (r->state_products[i][0] == 'g') {
            print_reaction(copy1,copy2);
            error->all(FLERR,
                       "Gas phase species must be second product "
                       "in DA reaction");
          }
        }
        */
        break;
      }

    case LH1:
      {
        if (r->state_products[0][0] != 'g') {
          print_reaction(copy1,copy2);
          error->all(FLERR,
                     "First product must be gas phase in LH1 reaction");
        }
        break;
      }

    case LH3:
      {
        if (r->state_products[0][0] != 's') {
          print_reaction(copy1,copy2);
          error->all(FLERR,
                     "First product must be surface phase in LH3 reaction");
        }
        break;
      }

    case CD:
      {
        if (r->state_products[0][0] != 'b') {
          print_reaction(copy1,copy2);
          error->all(FLERR,"First product must be bulk phase in CD reaction");
        }
        break;
      }

    case ER:
      {
        if (r->state_products[0][0] != 'g') {
          print_reaction(copy1,copy2);
          error->all(FLERR,"First product must be gas phase in ER reaction");
        }
        if (r->kisliuk_flag)
        {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Kisliuk option can only be used to define the reaction rate constant in adsorption mediated reaction such as AA, DA, LH1, LH3, and CD");
        }
        break;
      }

    case CI:
      {
        if (r->state_products[0][0] != 'g') {
          print_reaction(copy1,copy2);
          error->all(FLERR,"First product must be gas phase in CI reaction");
        }

        for (int i=2; i < r->nproduct; i++) {
          if (r->state_products[i][0] == 'g'){
            print_reaction(copy1,copy2);
            error->all(FLERR,"Gas phase species must be first or "
                       "second product in CI reaction");
          }
        }

        if (r->kisliuk_flag)
        {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Kisliuk option can only be used to define the reaction rate constant in adsorption mediated reaction such as AA, DA, LH1, LH3, and CD");
        }
        break;
      }
    }

    r->k_react = r->coeff[0];
    if (r->style == ARRHENIUS)
      r->k_react = r->k_react * pow(twall,r->coeff[1]) *
        exp(-r->coeff[2]/(twall));

    // process 3rd line of reaction
    // NOTE: RIGHT HERE

    // nextra = # of extra lines to read in this reaction: 0,1,2
    // please add some code that computes nextra based on gas species count

    int nextra = r->nprod_g;

    // NOTE: END of ADDED CODE

    if (nextra == 0) {
      nlist_gs++;
      continue;
    }

    if (me == 0) eof = readextra(nextra,line1,line2,n1,n2);
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) error->all(FLERR,"Missing line(s) for collision model to use");

    MPI_Bcast(&n1,1,MPI_INT,0,world);
    MPI_Bcast(line1,n1,MPI_CHAR,0,world);
    if (nextra == 2) {
      MPI_Bcast(&n2,1,MPI_INT,0,world);
      MPI_Bcast(line2,n2,MPI_CHAR,0,world);
    }

    // process reactions for particles IP and JP, if they exist

    for (int m = 1; m <= nextra; m++) {
      int nwords;
      if (m == 1) nwords = input->count_words(line1);
      if (m == 2) nwords = input->count_words(line2);

      if (nwords == 0) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction collide style in file");
      }

      // convert line into char **words
      // leave space for leading ID arg = "react_adsorb"
      // doesn't matter that all saved models will have same ID

      nwords++;
      char **words = new char*[nwords];
      words[0] = (char *) "react_adsorb";

      int narg = 1;
      if (m == 1) words[narg] = strtok(line1," \t\n");
      if (m == 2) words[narg] = strtok(line2," \t\n");
      narg++;

      while (narg < nwords) {
        words[narg] = strtok(NULL," \t\n");
        narg++;
      }

      // create an instance of a SurfCollide model
      // supported models are in enum above and in this if/then/else block

      int model,nflags,ncoeffs;
      int *flags = NULL;
      double *coeffs = NULL;
      SurfCollide *sc;

      if (strcmp(words[1],"none") == 0) {
        model = NOMODEL;
        nflags = ncoeffs = 0;
        sc = NULL;
      } else if (strcmp(words[1],"specular") == 0) {
        model = SPECULAR;
        nflags = 1;
        ncoeffs = 0;
        sc = new SurfCollideSpecular(sparta,nwords,words);
      } else if (strcmp(words[1],"diffuse") == 0) {
        model = DIFFUSE;
        nflags = 0;
        ncoeffs = 2;
        sc = new SurfCollideDiffuse(sparta,nwords,words);
      } else if (strcmp(words[1],"adiabatic") == 0) {
        model = ADIABATIC;
        nflags = ncoeffs = 0;
        sc = new SurfCollideAdiabatic(sparta,nwords,words);
      } else if (strcmp(words[1],"cll") == 0) {
        model = CLL;
        nflags = 1;
        ncoeffs = 5;
        sc = new SurfCollideCLL(sparta,nwords,words);
      } else if (strcmp(words[1],"td") == 0) {
        model = TD;
        nflags = 3;
        ncoeffs = 8;
        sc = new SurfCollideTD(sparta,nwords,words);
      } else if (strcmp(words[1],"impulsive") == 0) {
        model = IMPULSIVE;
        nflags = 4;
        ncoeffs = 11;
        sc = new SurfCollideImpulsive(sparta,nwords,words);
      } else {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction collide style in file");
      }

      delete [] words;

      // use parsing from SurfCollide instance to fill flags, coeffs

      if (nflags) flags = new int[nflags];
      if (ncoeffs) coeffs = new double[ncoeffs];
      if (model != NOMODEL) sc->flags_and_coeffs(flags,coeffs);

      // if first SurfCollide instance of a particular style, save in cmodels
      // else delete it

      if (cmodels[model] == NULL) cmodels[model] = sc;
      else delete sc;

      // populate reaction-specific collision model settings

      if (m == 1) {
        r->cmodel_ip = model;
        r->cmodel_ip_flags = flags;
        r->cmodel_ip_coeffs = coeffs;
      } else if (m == 2) {
        r->cmodel_jp = model;
        r->cmodel_jp_flags = flags;
        r->cmodel_jp_coeffs = coeffs;
      }
    }

    // if nextra == 2, check that IP and JP models are both defined or neither

    if (nextra == 2)
      if ((r->cmodel_ip == NOMODEL && r->cmodel_jp != NOMODEL) ||
          (r->cmodel_ip != NOMODEL && r->cmodel_jp == NOMODEL)) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Both or neither collide styles must be defined");
      }

    // increment reaction count

    nlist_gs++;
  }

  // close reaction file

  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

char *SurfReactAdsorb::reactionID(int m)
{
  if (m < nlist_gs) return rlist_gs[m].id;
  return rlist_ps[m-nlist_gs].id;
}

/* ---------------------------------------------------------------------- */

int SurfReactAdsorb::match_reactant(char *species, int m)
{
  for (int i = 0; i < rlist_gs[m].nreactant; i++)
    if (strcmp(species,rlist_gs[m].id_reactants[i]) == 0) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

int SurfReactAdsorb::match_product(char *species, int m)
{
  for (int i = 0; i < rlist_gs[m].nproduct; i++)
    if (strcmp(species,rlist_gs[m].id_products[i]) == 0) return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorb::init_reactions_ps()
{
  // convert species IDs to species indices_ps
  // flag reactions as active/inactive_ps depending on whether all species exist

  for (int m = 0; m < nlist_ps; m++) {
    OneReaction_PS *r = &rlist_ps[m];
    r->active = 1;
    for (int i = 0; i < r->nreactant; i++) {
      r->reactants[i] = particle->find_species(r->id_reactants[i]);
      if (r->reactants[i] < 0) {
        r->active = 0;
        break;
      }
      if (r->state_reactants[i][0] == 's') {
        r->reactants_ad_index[i] = find_surf_species(r->id_reactants[i]);
        if ( r->reactants_ad_index[i] < 0) {
        r->active = 0;
        break;
        }
    }
    else r->reactants_ad_index[i] = -1; // SGK_change
    }
    for (int i = 0; i < r->nproduct; i++) {
      r->products[i] = particle->find_species(r->id_products[i]);
      if (r->products[i] < 0) {
        r->active = 0;
        break;
      }
      if (r->state_products[i][0] == 's') {
        r->products_ad_index[i] = find_surf_species(r->id_products[i]);
        if ( r->products_ad_index[i] < 0) {
        r->active = 0;
        break;
        }
      }
      else r->products_ad_index[i] = -1;
    }
  }

  // count possible reactions for each species

  nactive_ps = 0;
  for (int m = 0; m < nlist_ps; m++) {
    OneReaction_PS *r = &rlist_ps[m];
    if (!r->active) continue;
    nactive_ps++;
  }

  memory->destroy(reactions_ps_list);
  memory->create(reactions_ps_list,nactive_ps,"surf_adsorb:reactions_ps_list");
  int n = 0;

  for (int m = 0; m < nlist_ps; m++) {
    OneReaction_PS *r = &rlist_ps[m];
    if (!r->active) continue;
    reactions_ps_list[n] = m;
    n++;
  }
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorb::readfile_ps(char *fname)
{
  int n,n1,n2,eof;
  char line1[MAXLINE],line2[MAXLINE];
  char copy1[MAXLINE],copy2[MAXLINE];
  char *word;
  OneReaction_PS *r;

  // proc 0 opens file

  if (me == 0) {
    fp = fopen(fname,"r");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open reaction file %s",fname);
      error->one(FLERR,str);
    }
  }

  // read reactions one at a time and store their info in rlist_ps

  while (1) {
    if (me == 0) eof = readone(line1,line2,n1,n2);
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;

    MPI_Bcast(&n1,1,MPI_INT,0,world);
    MPI_Bcast(&n2,1,MPI_INT,0,world);
    MPI_Bcast(line1,n1,MPI_CHAR,0,world);
    MPI_Bcast(line2,n2,MPI_CHAR,0,world);

    if (nlist_ps == maxlist_ps) {
      maxlist_ps += DELTALIST;
      rlist_ps = (OneReaction_PS *)
        memory->srealloc(rlist_ps,maxlist_ps*sizeof(OneReaction_PS),
                         "react/tce:rlist_ps");
      for (int i = nlist_ps; i < maxlist_ps; i++) {
        r = &rlist_ps[i];
        r->nreactant = r->nproduct = 0;
        r->nprod_g = r->nprod_g_tot = 0;
        r->id = NULL;
        r->id_reactants = new char*[MAXREACTANT_PS];
        r->id_products = new char*[MAXPRODUCT_PS];
        r->state_reactants = new char*[MAXREACTANT_PS];
        r->state_products = new char*[MAXPRODUCT_PS];
        r->part_reactants = new int[MAXREACTANT_PS];
        r->part_products = new int[MAXPRODUCT_PS];
        r->stoich_reactants = new int[MAXREACTANT_PS];
        r->stoich_products = new int[MAXPRODUCT_PS];
        r->reactants = new int[MAXREACTANT_PS];
        r->products = new int[MAXPRODUCT_PS];
        r->reactants_ad_index = new int[MAXREACTANT_PS];
        r->products_ad_index = new int[MAXPRODUCT_PS];
        r->coeff = new double[MAXCOEFF_PS];
        r->cmodel_ip = NOMODEL;
        r->cmodel_ip_flags = NULL;
        r->cmodel_ip_coeffs = NULL;
        r->cmodel_jp = NOMODEL;
        r->cmodel_jp_flags = NULL;
        r->cmodel_jp_coeffs = NULL;
      }
    }

    strcpy(copy1,line1);
    strcpy(copy2,line2);

    r = &rlist_ps[nlist_ps];
    r->index = n_PS_react;
    n_PS_react++;

    int side = 0;
    int species = 1;
    int start = 0;

    n = strlen(line1) - 1;
    r->id = new char[n+1];
    strncpy(r->id,line1,n);
    r->id[n] = '\0';

    word = strtok(line1," \t\n");

    while (1) {
      if (!word) {
        if (side == 0) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        if (species) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        break;
      }
      if (species) {
        species = 0;
        if (side == 0) {
          if (r->nreactant == MAXREACTANT_PS) {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Too many reactants in a reaction formula");
          }
          n = strlen(word) + 1;
          start = 0;
          r->part_reactants[r->nreactant] = 1;
          r->stoich_reactants[r->nreactant] = 1;
          if (word[n-2] != ')') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Specify the state of the reactants");
          }
          if (word[n-3] == 'c') {
            r->part_reactants[r->nreactant] = 0;
            n--;
          }
          if (word[n-3] != 's' && word[n-3] != 'b') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Only adsorbed and bulk species can be "
                       "reactants of an PS reaction");
          }
          if (isdigit(word[0])) {
            r->stoich_reactants[r->nreactant] = atoi(word);  //word[0] - '0';
            start=1;
          }
          r->id_reactants[r->nreactant] = new char[n-start-3]();
          strncpy(r->id_reactants[r->nreactant],&(word[start]),n-start-4);
          r->state_reactants[r->nreactant] = new char[1]();
          strncpy(r->state_reactants[r->nreactant],1+strstr(word,"("),1);
          r->nreactant++;
        } else {
          if (r->nproduct == MAXPRODUCT_PS) {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Too many products in a reaction formula");
          }

          n = strlen(word) + 1;
          start = 0;
          r->part_products[r->nproduct] = 1;
          r->stoich_products[r->nproduct] = 1;
          if (word[n-2] != ')') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Specify the state of the products");
          }
          if (word[n-3] == 'c') {
            r->part_products[r->nproduct] = 0;
            n--;
          }
          if (isdigit(word[0])) {
            r->stoich_products[r->nproduct] = atoi(word);  //word[0] - '0';
            start=1;
          }
          r->id_products[r->nproduct] = new char[n-start-3]();
          strncpy(r->id_products[r->nproduct],&(word[start]),n-start-4);
          r->state_products[r->nproduct] = new char[1]();
          strncpy(r->state_products[r->nproduct],1+strstr(word,"("),1);

          if(r->state_products[r->nproduct][0] == 'g')
            {
              r->nprod_g++;
              r->nprod_g_tot += r->stoich_products[r->nproduct];
            }

          r->nproduct++;
        }
      } else {
        species = 1;
        if (strcmp(word,"+") == 0) {
          word = strtok(NULL," \t\n");
          continue;
        }
        if (strcmp(word,"-->") != 0) {
          print_reaction(copy1,copy2);
          error->all(FLERR,"Invalid reaction formula in file");
        }
        side = 1;
      }
      word = strtok(NULL," \t\n");
    }

    // process 2nd line of reaction

    word = strtok(line2," \t\n");
    if (!word) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction type in file");
    }

    r->ncoeff = 0;

    if (strcmp(word,"DS")==0 || strcmp(word,"Ds")==0 ||
        strcmp(word,"dS")==0|| strcmp(word,"ds")==0) {
      r->type = DS;
    } else if (strcmp(word,"LH2") == 0 || strcmp(word,"Lh2") == 0 ||
               strcmp(word,"lH2") == 0 || strcmp(word,"lh2") == 0) {
      r->type = LH2;
    } else if (strcmp(word,"LH4") == 0 || strcmp(word,"Lh4") == 0 ||
               strcmp(word,"lH4") == 0 || strcmp(word,"lh4") == 0) {
      r->type = LH4;
    } else if (strcmp(word,"SB")==0 || strcmp(word,"Sb")==0 ||
               strcmp(word,"sB")==0|| strcmp(word,"sb")==0) {
      r->type = SB;
    } else {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction type in file");
    }

    word = strtok(NULL," \t\n");
    if (!word) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction style in file");
    }
    if (word[0] == 'A' || word[0] == 'a') {r->style = ARRHENIUS;r->ncoeff += 3;}
    else if (word[0] == 'S' || word[0] == 's') {r->style = SIMPLE;r->ncoeff += 1;}
    else {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Invalid reaction style in file");
    }

    for (int i = 0; i < r->ncoeff; i++) {
      word = strtok(NULL," \t\n");
      if (!word) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction coefficients in file");
      }
      r->coeff[i] = input->numeric(FLERR,word);
    }

    word = strtok(NULL," \t\n");
    if (word) {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Too many coefficients in a reaction formula");
    }

    // ERROR CHECKS
    // check that reactant/product counts are consistent with type

    if (r->nprod_g_tot > 2)
    {
      print_reaction(copy1,copy2);
      error->all(FLERR,"Number of gas phase products cannot be greater than 2");
    }

    switch (r->type) {

      case DS:
        {
          if (r->state_reactants[0][0] != 's') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"First reactant must be surface phase "
                       "in DS reaction");
            }

          if (r->state_products[0][0] != 'g') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"First product must be gas phase in DS reaction");
            }

          for (int i=1; i < r->nproduct; i++) {
            if (r->state_products[i][0] == 'g')
              {
                print_reaction(copy1,copy2);
                error->all(FLERR,"Fas phase species must be "
                           "first product in DS reaction");
              }
          }
          break;
        }

    case LH2:
      {
        if (r->state_reactants[0][0] != 's') {
          print_reaction(copy1,copy2);
          error->all(FLERR,"First reactant must be surface phase in "
                     "LH2 reaction");
        }

        if (r->state_products[0][0] != 'g') {
          print_reaction(copy1,copy2);
          error->all(FLERR,"First product must be gas phase in "
                     "LH2 reaction");
        }

        for (int i=1; i < r->nproduct; i++) {
          if (r->state_products[i][0] == 'g')
            {
              print_reaction(copy1,copy2);
              error->all(FLERR,"Gas phase species must be "
                         "first product in LH2 reaction");
            }
        }
        break;
      }

    case LH4:
      {
        if (r->state_reactants[0][0] != 's') {
          print_reaction(copy1,copy2);
          error->all(FLERR,"First reactant must be surface phase in "
                     "LH4 reaction");
        }

        if (r->state_products[0][0] != 's') {
          print_reaction(copy1,copy2);
          error->all(FLERR,"First product must be surface phase in "
                     "LH4 reaction");
        }
        break;
      }

    case SB:
      {
        if (r->state_reactants[0][0] != 'b') {
          print_reaction(copy1,copy2);
          error->all(FLERR,"First reactant must be bulk phase in SB reaction");
        }

        if (r->state_products[0][0] != 'g') {
          print_reaction(copy1,copy2);
          error->all(FLERR,"First product must be gas phase in SB reaction");
        }

        for (int i=1; i < r->nproduct; i++) {
          if (r->state_products[i][0] == 'g') {
            print_reaction(copy1,copy2);
            error->all(FLERR,"Gas phase species must be "
                       "first product in SB reaction");
          }
        }

        break;
      }
    }

    r->k_react = r->coeff[0];
    if (r->style == ARRHENIUS) r->k_react = r->k_react * pow(twall,r->coeff[1]) *
                                 exp(-r->coeff[2]/(twall));
    //nlist_ps++;

    // process 3rd line of reaction
    // NOTE: RIGHT HERE

    // nextra = # of extra lines to read in this reaction: 0,1,2
    // please add some code that computes nextra based on gas species count

    int nextra = r->nprod_g;

    // NOTE: END of ADDED CODE

    if (nextra == 0) {
      nlist_ps++;
      continue;
    }

    if (me == 0) eof = readextra(nextra,line1,line2,n1,n2);
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) error->all(FLERR,"Missing line(s) for collision model to use");

    MPI_Bcast(&n1,1,MPI_INT,0,world);
    MPI_Bcast(line1,n1,MPI_CHAR,0,world);
    if (nextra == 2) {
      MPI_Bcast(&n2,1,MPI_INT,0,world);
      MPI_Bcast(line2,n2,MPI_CHAR,0,world);
    }

    // process reactions for particles IP and JP, if they exist

    for (int m = 1; m <= nextra; m++) {
      int nwords;
      if (m == 1) nwords = input->count_words(line1);
      if (m == 2) nwords = input->count_words(line2);

      if (nwords == 0) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction collide style in file");
      }

      // convert line into char **words
      // leave space for leading ID arg = "react_adsorb"
      // doesn't matter that all saved models will have same ID

      nwords++;
      char **words = new char*[nwords];
      words[0] = (char *) "react_adsorb";

      int narg = 1;
      if (m == 1) words[narg] = strtok(line1," \t\n");
      if (m == 2) words[narg] = strtok(line2," \t\n");
      narg++;

      while (narg < nwords) {
        words[narg] = strtok(NULL," \t\n");
        narg++;
      }

      // create an instance of a SurfCollide model
      // supported models are in enum above and in this if/then/else block

      int model,nflags,ncoeffs;
      int *flags = NULL;
      double *coeffs = NULL;
      SurfCollide *sc;

      if (strcmp(words[1],"none") == 0) {
        model = NOMODEL;
        nflags = ncoeffs = 0;
        sc = NULL;
      } else if (strcmp(words[1],"specular") == 0) {
        model = SPECULAR;
        nflags = 1;
        ncoeffs = 0;
        sc = new SurfCollideSpecular(sparta,nwords,words);
      } else if (strcmp(words[1],"diffuse") == 0) {
        model = DIFFUSE;
        nflags = 0;
        ncoeffs = 2;
      } else if (strcmp(words[1],"adiabatic") == 0) {
        model = ADIABATIC;
        nflags = ncoeffs = 0;
        sc = new SurfCollideAdiabatic(sparta,nwords,words);
      } else if (strcmp(words[1],"cll") == 0) {
        model = CLL;
        nflags = 1;
        ncoeffs = 5;
        sc = new SurfCollideCLL(sparta,nwords,words);
      } else if (strcmp(words[1],"td") == 0) {
        model = TD;
        nflags = 3;
        ncoeffs = 8;
        sc = new SurfCollideTD(sparta,nwords,words);
      } else if (strcmp(words[1],"impulsive") == 0) {
        model = IMPULSIVE;
        nflags = 4;
        ncoeffs = 11;
        sc = new SurfCollideImpulsive(sparta,nwords,words);
      } else {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Invalid reaction collide style in file");
      }

      delete [] words;

      // use parsing from SurfCollide instance to fill flags, coeffs

      if (nflags) flags = new int[nflags];
      if (ncoeffs) coeffs = new double[ncoeffs];
      if (model != NOMODEL) sc->flags_and_coeffs(flags,coeffs);

      // if first SurfCollide instance of a particular style, save in cmodels
      // else delete it

      if (cmodels[model] == NULL) cmodels[model] = sc;
      else delete sc;

      // populate reaction-specific collision model settings

      if (m == 1) {
        r->cmodel_ip = model;
        r->cmodel_ip_flags = flags;
        r->cmodel_ip_coeffs = coeffs;
      } else if (m == 2) {
        r->cmodel_jp = model;
        r->cmodel_jp_flags = flags;
        r->cmodel_jp_coeffs = coeffs;
      }
    }

    // if nextra == 2, check that IP and JP models are both defined or neither

    if (nextra == 2)
      if ((r->cmodel_ip == NOMODEL && r->cmodel_jp != NOMODEL) ||
          (r->cmodel_ip != NOMODEL && r->cmodel_jp == NOMODEL)) {
        print_reaction(copy1,copy2);
        error->all(FLERR,"Both or neither collide styles must be defined");
      }

    // increment reaction count

    nlist_ps++;
  }

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   perform on-surf reactions for one face or one surf element
   modePS = PSFACE, PSLINE, PSTRI
   isurf = 0 to 5 for box faces
   isurf >= 0 for line or tri indexed from 0 to Nsurf-1
   invoked once per Nsync steps
------------------------------------------------------------------------- */

void SurfReactAdsorb::PS_react(int modePS, int isurf, double *norm)
{
  // mark this surface element since performing on-surf chemistry

  if (mode == SURF) mark[isurf] = 1;

  // use these 5 data structs for either faces or surface elements
  // in either case, can be indexed by isurf
  // NOTE: not seeing where species_state or total_state is used in code below?

  // int **species_delta;       // change in perspecies count since last sync
  // int **species_state;       // perspecies count at last sync
  // int *total_state;          // total count at last sync
  // double *area;              // area of surf
  // double *weight;            // weight of surf

  if (nactive_ps == 0) return;

  //1 Particle::Species *species = particle->species;
  //1 Particle::OnePart *particles;
  //1 particles = particle->particles;
  // int nlocal = particle->nlocal;

  // line or tri data

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  double fnum = update->fnum;
  double factor = fnum * weight[isurf] / area[isurf];
  double ms_inv = factor/max_cover;

  Particle::OnePart *p;
  int id,isc;

  double nu_react[nactive_ps];
  OneReaction_PS *r;
  int rxn_occur[nactive_ps];

  for (int i = 0; i < nactive_ps; i++) {
    r = &rlist_ps[reactions_ps_list[i]];
    //int react_num = r->index;
    rxn_occur[i] = 1;

    for (int j=0; j<r->nreactant; j++) {
      if (r->state_reactants[j][0] == 's') {
        if (species_state[isurf][r->reactants_ad_index[j]] < r->stoich_reactants[j]) rxn_occur[i] = 0;
      }
    }
    //if (rxn_occur[i]) tau[isurf][react_num] += update->dt;
    if (rxn_occur[i]) tau[isurf][i] += update->dt*nsync;
  }

  while (1) {
    long int sum_nu_tau = 0;
    long int nu_tau[nactive_ps];

    for (int i = 0; i < nactive_ps; i++) {
      nu_react[i] = 0.0;
      nu_tau[i] = 0;
      if (rxn_occur[i]) {
        r = &rlist_ps[reactions_ps_list[i]];
        //int react_num = r->index;


        nu_react[i] = r->k_react;

        if (r->type == SB) {
          double surf_cover = total_state[isurf] * ms_inv;
          nu_react[i] *= pow((1-surf_cover),r->stoich_reactants[0]);
        } else {
          int factor_pow = -1;
          for (int j=0; j<r->nreactant; j++) {
          if (r->state_reactants[j][0] == 's') {
            factor_pow += r->stoich_reactants[j];
            if (r->part_reactants[j] == 0) {
              nu_react[i] *= stoich_pow(total_state[isurf],r->stoich_reactants[j]);
            } else {
              nu_react[i] *= stoich_pow(species_state[isurf][r->reactants_ad_index[j]],r->stoich_reactants[j]);
            }
          }
        }


        nu_react[i] *= pow(ms_inv,factor_pow);
        }

        /*
        for (int j=0; j<r->nreactant; j++) {
          nu_react[i] *=
            stoich_pow(species_state[isurf][r->reactants_ad_index[j]],
                       r->stoich_reactants[j]);
          factor_pow += r->stoich_reactants[j];
        }
        */

        //nu_tau[i] = MAX(floor(nu_react[i] * tau[isurf][react_num]),0);
        nu_tau[i] = MAX(floor(nu_react[i] * tau[isurf][i]),0);
        sum_nu_tau += nu_tau[i];
      }
    }


    if (sum_nu_tau == 0) break;

    double sum_inv = 1.0/sum_nu_tau;
    double random_prob = random->uniform();
    double react_prob = 0.0;
    int check_break = 0;

    int m,ireaction;

    for (int i = 0; i < nactive_ps; i++) {
      react_prob += nu_tau[i]*sum_inv;
      if (react_prob > random_prob) {
        check_break++;

        r = &rlist_ps[reactions_ps_list[i]];
        //int react_num = r->index;

        // if computes which tally on-surface reactions exist:
        //    invoke them here so can be tallied on a per-surf basis
        //    Update::run() does same thing for gas/surf reactions

        nsingle++;
        ireaction = nlist_gs + reactions_ps_list[i];
        tally_single[ireaction]++;
        if (ncompute_tally)
          for (m = 0; m < ncompute_tally; m++)
            clist_active[m]->surf_tally(isurf,-1,ireaction,NULL,NULL,NULL);

        // update tau

        double t = -log(random->uniform())/nu_react[i];
        //tau[isurf][react_num] -= t;
        tau[isurf][i] -= t;

        for (int j=0;j<r->nreactant;j++) {
          if (r->part_reactants[j] == 1) {
            switch(r->state_reactants[j][0]) {
            case 's':
              {
                species_delta[isurf][r->reactants_ad_index[j]] -=
                  r->stoich_reactants[j];
              }
            case 'g': {}
            case 'b': {}
            }
          }
        }

        for (int j=0;j<r->nproduct;j++) {
          if (r->part_products[j] == 1) {
            switch(r->state_products[j][0]) {
            case 's':
              {
                species_delta[isurf][r->products_ad_index[j]] +=
                  r->stoich_products[j];
              }
            case 'g': {}
            case 'b': {}
            }
          }
        }

        // for each reaction, post-reaction velocities must be set
        // if NOMODEL then SC instance associated with surf/face sets vels
        // else SC style created when PS file was read sets velocities
        // calls to add_particle(), followed by add_particle_mine()
        //   for add_particle() use dummy icell = 0, will be reset later
        //   add_particle_mine() copies new particle to mypart list,
        //     then removes it from Particle class
        //   concatenated mypart list is processed in PS_chemistry()
        //   added particle's grid cells are identified by owning procs

        switch (r->type) {

        case DS:
          {
            double x[3],v[3];

            int id = MAXSMALLINT*random->uniform();
            random_point(isurf,x);
            v[0] = v[1] = v[2] = 0.0;

            particle->add_particle(id,r->products[0],0,x,v,0.0,0.0);
            p = &particle->particles[particle->nlocal-1];
            p->dtremain = update->dt*random->uniform();

            if (r->cmodel_ip != NOMODEL)
              cmodels[r->cmodel_ip]->wrapper(p,norm,r->cmodel_ip_flags,
                                             r->cmodel_ip_coeffs);
            else {
              if (modePS == PSFACE) isc = domain->surf_collide[isurf];
              else if (modePS == PSLINE) isc = lines[isurf].isc;
              else if (modePS == PSTRI) isc = tris[isurf].isc;
              surf->sc[isc]->wrapper(p,norm,NULL,NULL);
            }

            add_particle_mine(p);
            particle->nlocal--;

            break;
          }

        case LH2:
          {
            double x[3],v[3];

            int id = MAXSMALLINT*random->uniform();
            random_point(isurf,x);
            v[0] = v[1] = v[2] = 0.0;

            particle->add_particle(id,r->products[0],0,x,v,0.0,0.0);
            p = &particle->particles[particle->nlocal-1];
            p->dtremain = update->dt*random->uniform();

            if (r->cmodel_ip != NOMODEL)
              cmodels[r->cmodel_ip]->wrapper(p,norm,r->cmodel_ip_flags,
                                             r->cmodel_ip_coeffs);
            else {
              if (modePS == PSFACE) isc = domain->surf_collide[isurf];
              else if (modePS == PSLINE) isc = lines[isurf].isc;
              else if (modePS == PSTRI) isc = tris[isurf].isc;
              surf->sc[isc]->wrapper(p,norm,NULL,NULL);
            }

            add_particle_mine(p);
            particle->nlocal--;

            break;
          }

        case LH4:
          {
            break;
          }

        case SB:
          {
            double x[3],v[3];

            int id = MAXSMALLINT*random->uniform();
            random_point(isurf,x);
            v[0] = v[1] = v[2] = 0.0;

            particle->add_particle(id,r->products[0],0,x,v,0.0,0.0);
            p = &particle->particles[particle->nlocal-1];
            p->dtremain = update->dt*random->uniform();

            if (r->cmodel_ip != NOMODEL)
              cmodels[r->cmodel_ip]->wrapper(p,norm,r->cmodel_ip_flags,
                                             r->cmodel_ip_coeffs);
            else {
              if (modePS == PSFACE) isc = domain->surf_collide[isurf];
              else if (modePS == PSLINE) isc = lines[isurf].isc;
              else if (modePS == PSTRI) isc = tris[isurf].isc;
              surf->sc[isc]->wrapper(p,norm,NULL,NULL);
            }

            add_particle_mine(p);
            particle->nlocal--;

            break;
          }
        }

        break;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   add new particle P to mypart list
------------------------------------------------------------------------- */

void SurfReactAdsorb::add_particle_mine(Particle::OnePart *p)
{
  if (npart == maxmypart) {
    maxmypart += DELTA_PART;
    mypart = (AddParticle *)
      memory->srealloc(mypart,maxmypart*sizeof(AddParticle),"sr_adsorb:mypart");
  }

  mypart[npart].id = p->id;
  mypart[npart].ispecies = p->ispecies;
  memcpy(mypart[npart].x,p->x,3*sizeof(double));
  memcpy(mypart[npart].v,p->v,3*sizeof(double));
  mypart[npart].erot = p->erot;
  mypart[npart].evib = p->evib;
  mypart[npart].dtremain = p->dtremain;
  npart++;
}

/* ---------------------------------------------------------------------- */

/*
void SurfReactAdsorb::energy_barrier_scatter(Particle::OnePart *p, double *norm,
                                             double barrier_cos_pow,
                                             double sigma1, double sigma2)
{
  Particle::Species *species = particle->species;
  double tangent1[3],tangent2[3];
  int ispecies = p->ispecies;

  double *v = p->v;
  double mass = species[ispecies].mass;
  double E_i = 0.5 * mass * MathExtra::lensq3(v);

  double E_t = update->boltz*(twall+sigma2) + sigma1*E_i;
  double E_n = E_t + update->boltz*twall*0.5*(barrier_cos_pow-1);
  double vrm_n = sqrt(2.0*E_n / mass);
  double vrm_t = sqrt(2.0*E_t / mass);
  double vperp = vrm_n * sqrt(-log(random->uniform()));

  double theta = MY_2PI * random->uniform();
  double vtangent = vrm_t * sqrt(-log(random->uniform()));
  double vtan1 = vtangent * sin(theta);
  double vtan2 = vtangent * cos(theta);

  double dot = MathExtra::dot3(v,norm);

  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  if (MathExtra::lensq3(tangent1) == 0.0) {
    tangent2[0] = random->uniform();
    tangent2[1] = random->uniform();
    tangent2[2] = random->uniform();
    MathExtra::cross3(norm,tangent2,tangent1);
  }

  MathExtra::norm3(tangent1);
  MathExtra::cross3(norm,tangent1,tangent2);

  v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
  v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
  v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

  p->erot = particle->erot(ispecies,twall,random);
  p->evib = particle->evib(ispecies,twall,random);
}
*/

/* ---------------------------------------------------------------------- */

/*
void SurfReactAdsorb::non_thermal_scatter(Particle::OnePart *p, double *norm,
                                          double NT_alpha, double NT_u0_a,
                                          double NT_u0_b, double NT_barrier)
{
  Particle::Species *species = particle->species;
  double tangent1[3],tangent2[3];
  int ispecies = p->ispecies;

  double *v = p->v;
  double mass = species[ispecies].mass;

  double dot = MathExtra::dot3(v,norm);

  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  if (MathExtra::lensq3(tangent1) == 0.0) {
    tangent2[0] = random->uniform();
    tangent2[1] = random->uniform();
    tangent2[2] = random->uniform();
    MathExtra::cross3(norm,tangent2,tangent1);
  }

  MathExtra::norm3(tangent1);
  MathExtra::cross3(norm,tangent1,tangent2);

  double NT_u0 = NT_u0_a*twall + NT_u0_b;
  double NT_alpha_sq = NT_alpha * NT_alpha;

  double vrm_n = sqrt(2.0*update->boltz * (twall + NT_barrier) / mass);
  double vrm_t = sqrt(2.0*update->boltz * twall / mass);

  double NT_vf_max = 0.5 * (NT_u0 + sqrt(NT_u0*NT_u0 + 6*NT_alpha_sq));
  double NT_f_max = NT_vf_max*NT_vf_max*NT_vf_max *
    exp(-(NT_vf_max - NT_u0)*(NT_vf_max - NT_u0)/(NT_alpha_sq));

  double P = 0, NT_vf_mag;
  while (random->uniform() > P) {
    NT_vf_mag = NT_vf_max + 3 * NT_alpha * ( 2 * random->uniform() - 1 );
    P = NT_vf_mag*NT_vf_mag*NT_vf_mag/(NT_f_max) *
      exp(-(NT_vf_mag - NT_u0)*(NT_vf_mag - NT_u0)/(NT_alpha_sq));
  }

  double NT_phi = MY_2PI * random->uniform();
  double NT_theta = atan2(vrm_t * sqrt(-log(random->uniform())),vrm_n *
                          sqrt(-log(random->uniform())));

  double vperp = NT_vf_mag * cos(NT_theta);
  double vtan1 = NT_vf_mag * sin(NT_theta) * cos(NT_phi);
  double vtan2 = NT_vf_mag * sin(NT_theta) * sin(NT_phi);

  v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
  v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
  v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

  p->erot = particle->erot(ispecies,twall,random);
  p->evib = particle->evib(ispecies,twall,random);
}
*/

/* ---------------------------------------------------------------------- */

/*
void SurfReactAdsorb::cll(Particle::OnePart *p, double *norm, double acc_n,
                          double acc_t, double eccen)
{
  // cll reflection
  // vrm = most probable speed of species, eqns (4.1) and (4.7)
  // vperp = velocity component perpendicular to surface along norm, eqn (12.3)
  // vtan12 = 2 velocity components tangential to surface
  // tangent1 = component of particle v tangential to surface,
  //   check if tangent1 = 0 (normal collision), set randomly
  // tangent2 = norm x tangent1 = orthogonal tangential direction
  // tangent12 are both unit vectors

  Particle::Species *species = particle->species;
  double tangent1[3],tangent2[3];
  int ispecies = p->ispecies;

  double *v = p->v;
  double dot = MathExtra::dot3(v,norm);
  double tan = sqrt(MathExtra::lensq3(v) - dot*dot);

  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  if (MathExtra::lensq3(tangent1) == 0.0) {
    tangent2[0] = random->uniform();
    tangent2[1] = random->uniform();
    tangent2[2] = random->uniform();
    MathExtra::cross3(norm,tangent2,tangent1);
  }

  MathExtra::norm3(tangent1);
  MathExtra::cross3(norm,tangent1,tangent2);

  double tan1 = MathExtra::dot3(v,tangent1);
  double vrm = sqrt(2.0*update->boltz * twall / species[ispecies].mass);

  // CLL model normal velocity

  double r_1 = sqrt(-acc_n*log(random->uniform()));
  double theta_1 = MY_2PI * random->uniform();
  double dot_norm = fabs(dot/vrm) * sqrt(1-acc_n);
  double vperp = vrm * sqrt( r_1*r_1 + dot_norm*dot_norm +
                             2*r_1*dot_norm*cos(theta_1) );

  // CLL model tangential velocities

  double r_2 = sqrt(-acc_t*log(random->uniform()));
  double theta_2 = MY_2PI * random->uniform();
  double vtangent = fabs(tan/vrm) * sqrt(1-acc_t);
  double vtan1 = vrm * (vtangent + r_2*cos(theta_2));
  double vtan2 = vrm * r_2 * sin(theta_2);

  int pflag = 0;
  if (eccen >= 0 && eccen < 1) pflag = 1;

  if (pflag) {
    double tan2 = MathExtra::dot3(v,tangent2);
    double theta_i, phi_i, psi_i, theta_f, phi_f, psi_f, cos_beta;

    theta_i = acos(dot/sqrt(MathExtra::lensq3(v)));
    psi_i = acos(dot*dot/MathExtra::lensq3(v));
    phi_i = atan2(tan2,tan1);

    double v_mag = sqrt(vperp*vperp + vtan1*vtan1 + vtan2*vtan2);

    double P = 0;
    while (random->uniform() > P) {
      phi_f = MY_2PI*random->uniform();
      psi_f = acos(1-random->uniform());
      cos_beta =  cos(psi_i)*cos(psi_f) + sin(psi_i)*sin(psi_f) *
        cos(phi_i - phi_f);
      P = (1-eccen)/(1-eccen*cos_beta);
    }

    theta_f = acos(sqrt(cos(psi_f)));

    vperp = v_mag * cos(theta_f);
    vtan1 = v_mag * sin(theta_f) * cos(phi_f);
    vtan2 = v_mag * sin(theta_f) * sin(phi_f);
  }

  v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
  v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
  v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

  p->erot = particle->erot(ispecies,twall,random);
  p->evib = particle->evib(ispecies,twall,random);
}
*/

/* ----------------------------------------------------------------------
   pick a random point X on a face or surf
   for face: isurf = 0 to 5 inclusive for which face
   for surf: isurf = index to line or triangle
   return X
------------------------------------------------------------------------- */

void SurfReactAdsorb::random_point(int isurf, double *x)
{
  if (mode == FACE) {
    double *lo = domain->boxlo;
    double *hi = domain->boxhi;
    double rand1 = random->uniform();
    double rand2 = random->uniform();

    switch (isurf) {
    case XLO:
      {
        x[0] = lo[0];
        x[1] = lo[1] + rand1 * (hi[1] - lo[1]);
        x[2] = lo[2] + rand2 * (hi[2] - lo[2]);
        break;
      }

    case XHI:
      {
        x[0] = hi[0];
        x[1] = lo[1] + rand1 * (hi[1] - lo[1]);
        x[2] = lo[2] + rand2 * (hi[2] - lo[2]);
        break;
      }

    case YLO:
      {
        x[0] = lo[0] + rand2 * (hi[0] - lo[0]);
        x[1] = lo[1];
        x[2] = lo[2] + rand1 * (hi[2] - lo[2]);
        break;
      }

    case YHI:
      {
        x[0] = lo[0] + rand2 * (hi[0] - lo[0]);
        x[1] = hi[1];
        x[2] = lo[2] + rand1 * (hi[2] - lo[2]);
        break;
      }

    case ZLO:
      {
        x[0] = lo[0] + rand1 * (hi[0] - lo[0]);
        x[1] = lo[1] + rand2 * (hi[1] - lo[1]);
        x[2] = lo[2];
        break;
      }

    case ZHI:
      {
        x[0] = lo[0] + rand1 * (hi[0] - lo[0]);
        x[1] = lo[1] + rand2 * (hi[1] - lo[1]);
        x[2] = hi[2];
        break;
      }
    }

  } else if (mode == SURF) {
    if (domain->dimension == 2) {
      Surf::Line *lines = surf->lines;
      double *p1,*p2;
      double rand = random->uniform();

      p1 = lines[isurf].p1;
      p2 = lines[isurf].p2;

      x[0] = p1[0] + rand * (p2[0] - p1[0]);
      x[1] = p1[1] + rand * (p2[1] - p1[1]);
      x[2] = 0.0;

    } else if (domain->dimension == 3) {
      // NOTE: to avoid sqrt() could use 2 uniform RNs: r1,r2
      // if r1+r2 > 1 then r1 = 1-r1, r2 = 1-r2
      // x[i] = p1[i] + r1*(p2[i]-p1[i]) + r2*(p3[i]-p1[i])

      Surf::Tri *tris = surf->tris;
      double *p1,*p2,*p3;
      double rand1 = sqrt(random->uniform());
      double rand2 = random->uniform();
      double factor1 = 1-rand1;
      double factor2 = rand1*(1-rand2);
      double factor3 = rand1*rand2;

      p1 = tris[isurf].p1;
      p2 = tris[isurf].p2;
      p3 = tris[isurf].p3;

      x[0] = factor1 * p1[0] + factor2 * p2[0] + factor3 * p3[0];
      x[1] = factor1 * p1[1] + factor2 * p2[1] + factor3 * p3[1];
      x[2] = factor1 * p1[2] + factor2 * p2[2] + factor3 * p3[2];
    }
  }

  /*
  switch(element) {

  case GRID: {
    Grid::ChildCell *cells = grid->cells;
    double rand1 = random->uniform();
    double rand2 = random->uniform();
    double *lo = cells[ielem].lo;
    double *hi = cells[ielem].hi;

    double d_beam = 1.5e-3;
    double theta_beam = 45 * MY_PI /180;

    double rand_r = sqrt(random->uniform());
    double rand_angle = MY_2PI * random->uniform();

    double x_strike = 0.0;
    double y_strike = 0.0;
    double z_strike = 0.0;

    x[0] = x_strike ;
    x[1] = y_strike + 0.5 * d_beam / cos(theta_beam) * rand_r * cos(rand_angle);
    x[2] = z_strike + 0.5 * d_beam * rand_r * sin(rand_angle);
    break;
  }

  case LINE: {
    Surf::Line *lines = surf->lines;
    double *p1,*p2;
    double rand = random->uniform();

    p1 = lines[ielem].p1;
    p2 = lines[ielem].p2;

    x[0] = p1[0] + rand * (p2[0] - p1[0]);
    x[1] = p1[1] + rand * (p2[1] - p1[1]);
    x[2] = 0.0;
    break;
  }

  case TRI: {
    Surf::Tri *tris = surf->tris;
    double *p1,*p2,*p3;
    double rand1 = sqrt(random->uniform());
    double rand2 = random->uniform();
    double factor1 = 1-rand1;
    double factor2 = rand1*(1-rand2);
    double factor3 = rand1*rand2;

    p1 = tris[ielem].p1;
    p2 = tris[ielem].p2;
    p3 = tris[ielem].p3;

    x[0] = factor1 * p1[0] + factor2 * p2[0] + factor3 * p3[0];
    x[1] = factor1 * p1[1] + factor2 * p2[1] + factor3 * p3[1];
    x[2] = factor1 * p1[2] + factor2 * p2[2] + factor3 * p3[2];

    break;
  }
  }
  */
}

/* ---------------------------------------------------------------------- */

/*
int SurfReactAdsorb::find_cell(int isurf, double *x)
{
  int value = -1;
  switch(element) {

  case GRID: {
    value = ielem;
    break;
  }

  case LINE: {
    Surf::Line *lines = surf->lines;
    Grid::ChildCell *cells = grid->cells;
    for (int icell=0; icell<lines[ielem].ncell; icell++) {
      Grid::ChildCell *cell = &cells[lines[ielem].cell_list[icell]];
      if (x[0] <= cell->hi[0] && x[0] >= cell->lo[0] &&
          x[1] <= cell->hi[1] && x[1] >= cell->lo[1]) {
        value = icell;
        break;
      }
    }
    if (value == -1)
      error->all(FLERR,"Cell corresponding to the surface element was not found");
    break;
  }

  case TRI: {
    Surf::Tri *tris = surf->tris;
    Grid::ChildCell *cells = grid->cells;
    for (int icell=0; icell<tris[ielem].ncell; icell++) {
      Grid::ChildCell *cell = &cells[tris[ielem].cell_list[icell]];
      if (x[0] <= cell->hi[0] && x[0] >= cell->lo[0] &&
          x[1] <= cell->hi[1] && x[1] >= cell->lo[1] &&
          x[2] <= cell->hi[2] && x[2] >= cell->lo[2])  {
        value = icell;
        break;
      }
    }
    if (value == -1)
      error->all(FLERR,"Cell corresponding to the surface element is not found");
    break;
  }
  }

  return value;
}
*/

/* ---------------------------------------------------------------------- */

double SurfReactAdsorb::stoich_pow(int base, int pow)
{
  double value = 0.0;
  switch (pow) {

  case 0: {
    value = 1.0;
    break;
  }

  case 1: {
    if (base >= pow) value = double(base);
    break;
  }

  case 2: {
    if (base >= pow) value = 0.5*base*(base-1);
    break;
  }

  case 3: {
    if (base >= pow) value = 0.5*THIRD*base*(base-1)*(base-2);
    break;
  }

  case 4: {
    if (base >= pow) value = 0.125*THIRD*base*(base-1)*(base-2)*(base-3);
    break;
  }

  case 5: {
    if (base >= pow) value = 0.025*THIRD*base*(base-1)*(base-2)*(base-3)*(base-4);
    break;
  }

  case 6: {
    if (base >= pow)
      value = 0.0125*THIRD*THIRD*base*(base-1)*(base-2)*(base-3)*
        (base-4)*(base-5);
    break;
  }
  }

  return value;
}

/* ----------------------------------------------------------------------
   return index of ID in list of species IDs
   return -1 if not found
------------------------------------------------------------------------- */

int SurfReactAdsorb::find_surf_species(char *id)
{
  for (int i = 0; i < nspecies_surf; i++)
    if (strcmp(id,species_surf[i]) == 0) return i;
  return -1;
}

/* ----------------------------------------------------------------------
   read one reaction from file
   just first 2 lines, any extra are read by readextra()
   reaction = 2 lines of length n1 and n2
   return 1 if end-of-file, else return 0
------------------------------------------------------------------------- */

int SurfReactAdsorb::readone(char *line1, char *line2, int &n1, int &n2)
{
  char *eof;
  while ((eof = fgets(line1,MAXLINE,fp))) {
    int pre = strspn(line1," \t\n");
    if (pre == strlen(line1) || line1[pre] == '#') continue;
    eof = fgets(line2,MAXLINE,fp);
    if (!eof) break;
    n1 = strlen(line1) + 1;
    n2 = strlen(line2) + 1;
    return 0;
  }

  return 1;
}

/* ----------------------------------------------------------------------
   read nextra lines of a reaction, for surf collide models
   nextra = # of extra lines = 1 or 2
   line1 and line2 are extra lines of length n1 and n2
   return 1 if end-of-file, else return 0
------------------------------------------------------------------------- */

int SurfReactAdsorb::readextra(int nextra, char *line1, char *line2,
                               int &n1, int &n2)
{
  char *eof = fgets(line1,MAXLINE,fp);
  if (!eof) return 1;
  n1 = strlen(line1) + 1;

  if (nextra == 2) {
    eof = fgets(line2,MAXLINE,fp);
    if (!eof) return 1;
    n2 = strlen(line2) + 1;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorb::print_reaction(char *line1, char *line2)
{
  if (me) return;
  printf("Bad reaction format:\n");
  printf("%s%s",line1,line2);
};
