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

#include "string.h"
#include "react_tce_qk_kokkos.h"
#include "particle.h"
#include "collide.h"
#include "update.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ReactTCEQKKokkos::ReactTCEQKKokkos(SPARTA *sparta, int narg, char **arg) :
  ReactBirdKokkos(sparta, narg, arg) {}

/* ---------------------------------------------------------------------- */

void ReactTCEQKKokkos::init()
{
  if (!collide || (strcmp(collide->style,"vss") != 0 &&
                   strcmp(collide->style,"vss/kk") != 0))
    error->all(FLERR,"React tce/qk can only be used with collide vss");

  ReactBirdKokkos::init();

  // these two restrictions mirror the host (react_tce_qk.cpp:44-52) exactly --
  //   same condition, same message.  They are NOT Kokkos limitations:
  //   react tce/qk does not support recombination or compute_chem_rates on
  //   the CPU either, so lifting them here would diverge from the host.
  //   Only react tce supports both, and react_tce_kokkos.h ports that.

  for (int i = 0; i < nlist; i++)
    if (rlist[i].active && rlist[i].type == RECOMBINATION)
      error->all(FLERR,"React tce/qk does not currently support recombination reactions");
  if (computeChemRates)
    error->all(FLERR,"React tce/qk does not currently support the "
               "'react_modify compute_chem_rates' option");

  boltz = update->boltz;

  // flatten VSS omega for all species pairs to device

  int nspecies = particle->nspecies;
  d_omega = DAT::t_float_2d("react/tce/qk:omega",nspecies,nspecies);
  auto h_omega = Kokkos::create_mirror_view(d_omega);
  for (int i = 0; i < nspecies; i++)
    for (int j = 0; j < nspecies; j++)
      h_omega(i,j) = collide->extract(i,j,"omega");
  Kokkos::deep_copy(d_omega,h_omega);
}
