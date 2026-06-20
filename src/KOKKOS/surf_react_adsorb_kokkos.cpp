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
#include "string.h"
#include "surf_react_adsorb_kokkos.h"
#include "input.h"
#include "update.h"
#include "comm.h"
#include "domain.h"
#include "particle.h"
#include "error.h"
#include "particle_kokkos.h"
#include "sparta_masks.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

SurfReactAdsorbKokkos::SurfReactAdsorbKokkos(SPARTA *sparta, int narg, char **arg) :
  SurfReactAdsorb(sparta,narg,arg),
  rand_pool(12345 + comm->me
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  kokkosable = 1;

  d_scalars = DAT::t_int_1d("surf_react_adsorb:scalars",nlist_gs+1);
  d_nsingle = Kokkos::subview(d_scalars,0);
  d_tally_single = Kokkos::subview(d_scalars,std::make_pair(1,nlist_gs+1));

  h_scalars = HAT::t_int_1d("surf_react_adsorb:scalars_mirror",nlist_gs+1);
  h_nsingle = Kokkos::subview(h_scalars,0);
  h_tally_single = Kokkos::subview(h_scalars,std::make_pair(1,nlist_gs+1));

  random_backup = NULL;
}

SurfReactAdsorbKokkos::SurfReactAdsorbKokkos(SPARTA *sparta) :
  SurfReactAdsorb(sparta),
  rand_pool(12345
#ifdef SPARTA_KOKKOS_EXACT
            , sparta
#endif
            )
{
  copy = 1;
}

/* ---------------------------------------------------------------------- */

SurfReactAdsorbKokkos::~SurfReactAdsorbKokkos()
{
  if (copy) return;

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.destroy();
  if (random_backup) delete random_backup;
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorbKokkos::init()
{
  SurfReactAdsorb::init();

  // Kokkos GS adsorb currently supports a restricted feature set;
  //   error clearly at init rather than silently producing wrong results

  if (!gsflag)
    error->all(FLERR,"Kokkos surf_react adsorb requires gas-surface (gs) chemistry");
  if (psflag)
    error->all(FLERR,"Kokkos surf_react adsorb does not yet support on-surface (ps) chemistry");
  if (mode != SRA_KK::FACE)
    error->all(FLERR,"Kokkos surf_react adsorb only supports the box-face (face) option");

  for (int i = 0; i < nlist_gs; i++) {
    OneReaction_GS *r = &rlist_gs[i];
    if (!r->active) continue;
    // post-reaction collision model (cmodel) scatter on device currently
    //   supports NOMODEL and SPECULAR (no RNG); RNG-based cmodels deferred

    if ((r->cmodel_ip != SRA_KK::NOMODEL && r->cmodel_ip != SRA_KK::SPECULAR) ||
        (r->cmodel_jp != SRA_KK::NOMODEL && r->cmodel_jp != SRA_KK::SPECULAR))
      error->all(FLERR,"Kokkos surf_react adsorb does not yet support reactions with "
                 "a diffuse/cll/td/adiabatic/impulsive post-reaction collision model");
  }

  Kokkos::deep_copy(d_scalars,0);

  init_reactions_gs_kokkos();

#ifdef SPARTA_KOKKOS_EXACT
  rand_pool.init(random);
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorbKokkos::init_reactions_gs_kokkos()
{
  int nspecies = particle->nspecies;

  // per-species reaction lists

  int nmax = 0;
  d_reactions_n = DAT::t_int_1d("surf_react_adsorb:reactions_n",nspecies);
  auto h_reactions_n = Kokkos::create_mirror_view(d_reactions_n);
  for (int i = 0; i < nspecies; i++) {
    int n = reactions_gs[i].n;
    h_reactions_n(i) = n;
    nmax = MAX(nmax,n);
  }
  if (nmax > SRA_KK_MAXPERSPECIES)
    error->all(FLERR,"Too many Kokkos surf_react adsorb reactions per species");

  d_list = DAT::t_int_2d("surf_react_adsorb:list",nspecies,MAX(nmax,1));
  auto h_list = Kokkos::create_mirror_view(d_list);
  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < reactions_gs[i].n; j++)
      h_list(i,j) = reactions_gs[i].list[j];

  // flattened per-reaction tables

  int nr = MAX(nlist_gs,1);
  d_type = DAT::t_int_1d("sra:type",nr);
  d_style = DAT::t_int_1d("sra:style",nr);
  d_kreact = DAT::t_float_1d("sra:kreact",nr);
  d_kisliuk_flag = DAT::t_int_1d("sra:kflag",nr);
  d_kisliuk = DAT::t_float_2d("sra:kisliuk",nr,3);
  d_energy_flag = DAT::t_int_1d("sra:eflag",nr);
  d_energy = DAT::t_float_2d("sra:energy",nr,2);
  d_coeff = DAT::t_float_2d("sra:coeff",nr,SRA_KK_MAXCOEFF);
  d_nreactant = DAT::t_int_1d("sra:nreactant",nr);
  d_nproduct = DAT::t_int_1d("sra:nproduct",nr);
  d_nprod_g = DAT::t_int_1d("sra:nprod_g",nr);
  d_nprod_g_tot = DAT::t_int_1d("sra:nprod_g_tot",nr);
  d_cmodel_ip = DAT::t_int_1d("sra:cmodel_ip",nr);
  d_cmodel_jp = DAT::t_int_1d("sra:cmodel_jp",nr);
  d_rstate = DAT::t_int_2d("sra:rstate",nr,SRA_KK_MAXREACTANT);
  d_rpart = DAT::t_int_2d("sra:rpart",nr,SRA_KK_MAXREACTANT);
  d_rstoich = DAT::t_int_2d("sra:rstoich",nr,SRA_KK_MAXREACTANT);
  d_rad = DAT::t_int_2d("sra:rad",nr,SRA_KK_MAXREACTANT);
  d_pstate = DAT::t_int_2d("sra:pstate",nr,SRA_KK_MAXPRODUCT);
  d_ppart = DAT::t_int_2d("sra:ppart",nr,SRA_KK_MAXPRODUCT);
  d_pstoich = DAT::t_int_2d("sra:pstoich",nr,SRA_KK_MAXPRODUCT);
  d_pad = DAT::t_int_2d("sra:pad",nr,SRA_KK_MAXPRODUCT);
  d_products = DAT::t_int_2d("sra:products",nr,SRA_KK_MAXPRODUCT);

  auto h_type = Kokkos::create_mirror_view(d_type);
  auto h_style = Kokkos::create_mirror_view(d_style);
  auto h_kreact = Kokkos::create_mirror_view(d_kreact);
  auto h_kflag = Kokkos::create_mirror_view(d_kisliuk_flag);
  auto h_kisliuk = Kokkos::create_mirror_view(d_kisliuk);
  auto h_eflag = Kokkos::create_mirror_view(d_energy_flag);
  auto h_energy = Kokkos::create_mirror_view(d_energy);
  auto h_coeff = Kokkos::create_mirror_view(d_coeff);
  auto h_nreactant = Kokkos::create_mirror_view(d_nreactant);
  auto h_nproduct = Kokkos::create_mirror_view(d_nproduct);
  auto h_nprod_g = Kokkos::create_mirror_view(d_nprod_g);
  auto h_nprod_g_tot = Kokkos::create_mirror_view(d_nprod_g_tot);
  auto h_cmodel_ip = Kokkos::create_mirror_view(d_cmodel_ip);
  auto h_cmodel_jp = Kokkos::create_mirror_view(d_cmodel_jp);
  auto h_rstate = Kokkos::create_mirror_view(d_rstate);
  auto h_rpart = Kokkos::create_mirror_view(d_rpart);
  auto h_rstoich = Kokkos::create_mirror_view(d_rstoich);
  auto h_rad = Kokkos::create_mirror_view(d_rad);
  auto h_pstate = Kokkos::create_mirror_view(d_pstate);
  auto h_ppart = Kokkos::create_mirror_view(d_ppart);
  auto h_pstoich = Kokkos::create_mirror_view(d_pstoich);
  auto h_pad = Kokkos::create_mirror_view(d_pad);
  auto h_products = Kokkos::create_mirror_view(d_products);

  for (int i = 0; i < nlist_gs; i++) {
    OneReaction_GS *r = &rlist_gs[i];
    h_type(i) = r->type;
    h_style(i) = r->style;
    h_kreact(i) = r->k_react;
    h_kflag(i) = r->kisliuk_flag;
    for (int k = 0; k < 3; k++) h_kisliuk(i,k) = r->kisliuk_coeff[k];
    h_eflag(i) = r->energy_flag;
    for (int k = 0; k < 2; k++) h_energy(i,k) = r->energy_coeff[k];
    for (int k = 0; k < SRA_KK_MAXCOEFF; k++)
      h_coeff(i,k) = (k < r->ncoeff) ? r->coeff[k] : 0.0;
    h_nreactant(i) = r->nreactant;
    h_nproduct(i) = r->nproduct;
    h_nprod_g(i) = r->nprod_g;
    h_nprod_g_tot(i) = r->nprod_g_tot;
    h_cmodel_ip(i) = r->cmodel_ip;
    h_cmodel_jp(i) = r->cmodel_jp;
    for (int k = 0; k < r->nreactant && k < SRA_KK_MAXREACTANT; k++) {
      h_rstate(i,k) = r->state_reactants[k][0];
      h_rpart(i,k) = r->part_reactants[k];
      h_rstoich(i,k) = r->stoich_reactants[k];
      h_rad(i,k) = r->reactants_ad_index[k];
    }
    for (int k = 0; k < r->nproduct && k < SRA_KK_MAXPRODUCT; k++) {
      h_pstate(i,k) = r->state_products[k][0];
      h_ppart(i,k) = r->part_products[k];
      h_pstoich(i,k) = r->stoich_products[k];
      h_pad(i,k) = r->products_ad_index[k];
      h_products(i,k) = r->products[k];
    }
  }

  Kokkos::deep_copy(d_reactions_n,h_reactions_n);
  Kokkos::deep_copy(d_list,h_list);
  Kokkos::deep_copy(d_type,h_type);
  Kokkos::deep_copy(d_style,h_style);
  Kokkos::deep_copy(d_kreact,h_kreact);
  Kokkos::deep_copy(d_kisliuk_flag,h_kflag);
  Kokkos::deep_copy(d_kisliuk,h_kisliuk);
  Kokkos::deep_copy(d_energy_flag,h_eflag);
  Kokkos::deep_copy(d_energy,h_energy);
  Kokkos::deep_copy(d_coeff,h_coeff);
  Kokkos::deep_copy(d_nreactant,h_nreactant);
  Kokkos::deep_copy(d_nproduct,h_nproduct);
  Kokkos::deep_copy(d_nprod_g,h_nprod_g);
  Kokkos::deep_copy(d_nprod_g_tot,h_nprod_g_tot);
  Kokkos::deep_copy(d_cmodel_ip,h_cmodel_ip);
  Kokkos::deep_copy(d_cmodel_jp,h_cmodel_jp);
  Kokkos::deep_copy(d_rstate,h_rstate);
  Kokkos::deep_copy(d_rpart,h_rpart);
  Kokkos::deep_copy(d_rstoich,h_rstoich);
  Kokkos::deep_copy(d_rad,h_rad);
  Kokkos::deep_copy(d_pstate,h_pstate);
  Kokkos::deep_copy(d_ppart,h_ppart);
  Kokkos::deep_copy(d_pstoich,h_pstoich);
  Kokkos::deep_copy(d_pad,h_pad);
  Kokkos::deep_copy(d_products,h_products);

  // per-face state device storage (SRA_KK::FACE mode)

  d_total_state = DAT::t_int_1d("sra:total_state",nface);
  d_area = DAT::t_float_1d("sra:area",nface);
  d_weight = DAT::t_float_1d("sra:weight",nface);
  d_species_state = DAT::t_int_2d("sra:species_state",nface,nspecies_surf);

  k_species_delta = DAT::tdual_int_2d("sra:species_delta",nface,nspecies_surf);
  d_species_delta = k_species_delta.view_device();
  Kokkos::deep_copy(d_species_delta,0);
}

/* ----------------------------------------------------------------------
   sync per-face state host->device and refresh particle/scalar views
   called each step from the surf collide pre_collide
------------------------------------------------------------------------- */

void SurfReactAdsorbKokkos::pre_react()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  particle_kk->sync(Device,PARTICLE_MASK|SPECIES_MASK);
  d_particles = particle_kk->k_particles.view_device();
  d_species = particle_kk->k_species.view_device();

  fnum_ = update->fnum;

  // copy current per-face state (changes only at sync) host->device

  auto h_total = Kokkos::create_mirror_view(d_total_state);
  auto h_area = Kokkos::create_mirror_view(d_area);
  auto h_weight = Kokkos::create_mirror_view(d_weight);
  auto h_sstate = Kokkos::create_mirror_view(d_species_state);
  for (int i = 0; i < nface; i++) {
    h_total(i) = total_state[i];
    h_area(i) = area[i];
    h_weight(i) = weight[i];
    for (int j = 0; j < nspecies_surf; j++)
      h_sstate(i,j) = species_state[i][j];
  }
  Kokkos::deep_copy(d_total_state,h_total);
  Kokkos::deep_copy(d_area,h_area);
  Kokkos::deep_copy(d_weight,h_weight);
  Kokkos::deep_copy(d_species_state,h_sstate);
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorbKokkos::tally_reset()
{
  SurfReact::tally_reset();
  Kokkos::deep_copy(d_scalars,0);
}

/* ----------------------------------------------------------------------
   bring device tallies + per-face deltas to host, then run the host
   state-sync logic (MPI reduce + per-face state update), then re-zero
------------------------------------------------------------------------- */

void SurfReactAdsorbKokkos::tally_update()
{
  // device -> host: reaction counts

  Kokkos::deep_copy(h_scalars,d_scalars);
  nsingle = h_nsingle();
  for (int i = 0; i < nlist_gs; i++) tally_single[i] = h_tally_single[i];

  // device -> host: per-face perspecies deltas accumulated since last sync

  k_species_delta.modify_device();
  k_species_delta.sync_host();
  auto h_delta = k_species_delta.view_host();
  for (int i = 0; i < nface; i++)
    for (int j = 0; j < nspecies_surf; j++)
      species_delta[i][j] = h_delta(i,j);

  // host logic: accumulate tallies and (every nsync) MPI-sync per-face state;
  //   update_state_face() also re-zeros host species_delta

  SurfReactAdsorb::tally_update();

  // mirror re-zeroed host deltas back to device (only changed on a sync step)

  if (update->ntimestep % nsync == 0) {
    for (int i = 0; i < nface; i++)
      for (int j = 0; j < nspecies_surf; j++)
        h_delta(i,j) = species_delta[i][j];
    k_species_delta.modify_host();
    k_species_delta.sync_device();
    Kokkos::deep_copy(d_scalars,0);
  }
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorbKokkos::backup()
{
  ParticleKokkos* particle_kk = (ParticleKokkos*) particle;
  d_particles = particle_kk->k_particles.view_device();

#ifdef SPARTA_KOKKOS_EXACT
  if (!random_backup)
    random_backup = new RanKnuth(12345 + comm->me);
  memcpy(random_backup,random,sizeof(RanKnuth));
#endif
}

/* ---------------------------------------------------------------------- */

void SurfReactAdsorbKokkos::restore()
{
#ifdef SPARTA_KOKKOS_EXACT
  memcpy(random,random_backup,sizeof(RanKnuth));
#endif
}
