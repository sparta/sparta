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

#include "stdlib.h"
#include "string.h"
#include "scale_particles.h"
#include "update.h"
#include "particle.h"
#include "grid.h"
#include "mixture.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

ScaleParticles::ScaleParticles(SPARTA *sparta) : Pointers(sparta) {}

/* ---------------------------------------------------------------------- */

void ScaleParticles::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot scale particles before grid is defined");

  if (narg != 2) error->all(FLERR,"Illegal scale_particles command");

  int imix = particle->find_mixture(arg[0]);
  if (imix < 0) error->all(FLERR,"Scale_particles mixture ID does not exist");
  particle->mixture[imix]->init();

  double factor = atof(arg[1]);
  if (factor < 0.0) error->all(FLERR,"Illegal scale_particles command");

  // RNG for cloning/deletion

  RanKnuth *random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);

  // nbefore = current total # of particles

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  bigint nbefore;
  bigint nme = particle->nlocal;
  MPI_Allreduce(&nme,&nbefore,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // scale particles in mixture
  // clone or delete randomly based on factor

  Particle::OnePart *particles = particle->particles;
  int nbytes = sizeof(Particle::OnePart);
  int *s2g = particle->mixture[imix]->species2group;
  int ncustom = particle->ncustom;

  int m,ispecies,igroup,nclone,flag;
  double fraction;

  int nlocal = particle->nlocal;
  int nlocal_original = nlocal;
  int i = 0;

  while (i < nlocal_original) {

    // skip if particle species not in mixture

    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;

    // factor < 1.0 is candidate for deletion
    // if deleted and particle that takes its place is cloned (Nloc > Norig)
    //   then skip it via i++, else will examine it on next iteration

    if (factor < 1.0) {
      if (random->uniform() > factor) {
        memcpy(&particles[i],&particles[nlocal-1],nbytes);
        if (ncustom) particle->copy_custom(i,nlocal-1);
        if (nlocal > nlocal_original) i++;
        else nlocal_original--;
        particle->nlocal--;
        nlocal--;
      } else i++;
      continue;
    }

    // factor > 1.0 is candidate for cloning
    // create Nclone new particles each with unique ID

    nclone = static_cast<int> (factor);
    fraction = factor - nclone;
    nclone--;
    if (random->uniform() < fraction) nclone++;

    for (m = 0; m < nclone; m++) {
      flag = particle->clone_particle(i);
      if (flag) particles = particle->particles;
      nlocal++;
      particles[nlocal-1].id = MAXSMALLINT*random->uniform();
    }
    i++;
  }

  // nafter = new total # of particles

  bigint nafter;
  nme = particle->nlocal;
  MPI_Allreduce(&nme,&nafter,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // clean up

  delete random;

  // print stats

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Scaled particles\n");
      fprintf(screen,"  before = " BIGINT_FORMAT
              ", after = " BIGINT_FORMAT "\n",
              nbefore,nafter);
      fprintf(screen,"  CPU time = %g secs\n",time2-time1);
    }
    if (logfile) {
      fprintf(logfile,"Scaled particles\n");
      fprintf(logfile,"  before = " BIGINT_FORMAT
              ", after = " BIGINT_FORMAT "\n",
              nbefore,nafter);
      fprintf(logfile,"  CPU time = %g secs\n",time2-time1);
    }
  }
}
