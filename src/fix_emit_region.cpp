/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "fix_emit_region.h"
#include "update.h"
#include "domain.h"
#include "grid.h"
#include "particle.h"
#include "mixture.h"
#include "comm.h"
#include "modify.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // same as Grid
enum{MAXWELL,CHAPMAN};

/* ---------------------------------------------------------------------- */

FixEmitRegion::FixEmitRegion(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;

  if (narg < 5) error->all(FLERR,"Illegal fix emit/region command");
                  
  imix = particle->find_mixture(arg[2]);
  if (imix < 0) error->all(FLERR,"Fix emit/region mixture ID does not exist");

  int igroup = grid->find_group(arg[3]);
  if (igroup < 0) error->all(FLERR,"Fix emit/region group ID does not exist");
  groupbit = grid->bitmask[igroup];

  if (strcmp(arg[4],"mb") == 0) mode = MAXWELL;
  else if (strcmp(arg[4],"ce") == 0) mode = CHAPMAN;
  else error->all(FLERR,"Illegal fix emit/region command");

  if (mode == CHAPMAN)
    error->all(FLERR,"Fix emit/region does not yet support Chapman-Enskog distribution");

  // optional args

  np = 0;
  nevery = 1;
  
  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix emit/region command");
      np = ATOBIGINT(arg[iarg+1]);
      if (np < 0) error->all(FLERR,"Illegal emit/region command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nevery") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix emit/region command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix emit/region command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix emit/region command");
  }

  // RNG

  int me = comm->me;
  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,me,100);

  // counters common to all emit styles for output from fix

  nsingle = ntotal = 0;
}

/* ---------------------------------------------------------------------- */

FixEmitRegion::~FixEmitRegion()
{
  if (copymode) return;

  delete random;
}

/* ---------------------------------------------------------------------- */

int FixEmitRegion::setmask()
{
  int mask = 0;
  mask |= START_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEmitRegion::init()
{
  particle->exist = 1;
  ntotal = 0;

  // voltotal = weighted volume of insertion region

  int dimension = domain->dimension;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  double *lo,*hi;
  double volone;

  double volme = 0.0;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsurf) continue;
    if (cinfo[icell].type == INSIDE) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;

    if (dimension == 3) volone = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else if (domain->axisymmetric)
      volone = (hi[0]-lo[0]) * (hi[1]*hi[1]-lo[1]*lo[1])*MY_PI;
    else volone = (hi[0]-lo[0]) * (hi[1]-lo[1]);
    volme += volone / cinfo[icell].weight;
  }

  MPI_Allreduce(&volme,&voltotal,1,MPI_DOUBLE,MPI_SUM,world);

  // calculate Np if not set explicity

  if (np == 0) np = particle->mixture[imix]->nrho * voltotal / update->fnum;
}

/* ----------------------------------------------------------------------
   perform particle deletion/creation for cells I own
   participating cells are those in group and flow volume, with no surfs
------------------------------------------------------------------------- */

void FixEmitRegion::start_of_step()
{
  if (update->ntimestep % nevery) return;

  // delete all current particles from cell group
  // flag their cell as -1, then invoke Particle::compress_rebalance()
  
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  int nplocal = particle->nlocal;

  for (int ip = 0; ip < nplocal; ip++) {
    int icell = particles[ip].icell;
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsurf) continue;
    if (cinfo[icell].type == INSIDE) continue;
    
    particles[ip].icell = -1;
  }

  particle->compress_rebalance();

  // nfix_update_custom = # of fixes with update_custom() method

  particle->error_custom();
  modify->list_init_fixes();
  int nfix_update_custom = modify->n_update_custom;

  // add new particles to cell group
  // particle properties based on mode = MAXWELL or CHAPMAN
  // ntarget = floating point # of particles to create in one cell
  // npercell = integer # of particles to create in one cell

  int dimension = domain->dimension;
  int *species = particle->mixture[imix]->species;
  double *cummulative = particle->mixture[imix]->cummulative;
  double *vstream = particle->mixture[imix]->vstream;
  double *vscale = particle->mixture[imix]->vscale;
  double temp_thermal = particle->mixture[imix]->temp_thermal;
  double temp_rot = particle->mixture[imix]->temp_rot;
  double temp_vib = particle->mixture[imix]->temp_vib;
  int nglocal = grid->nlocal;

  double ntarget,rn,theta1,theta2,volone,vn,vr,erot,evib;
  int npercell,isp,ispecies,id;
  double x[3],v[3];
  double *lo,*hi;
  
  nsingle = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsurf) continue;
    if (cinfo[icell].type == INSIDE) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;

    if (dimension == 3) volone = (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]);
    else if (domain->axisymmetric)
      volone = (hi[0]-lo[0]) * (hi[1]*hi[1]-lo[1]*lo[1])*MY_PI;
    else volone = (hi[0]-lo[0]) * (hi[1]-lo[1]);

    ntarget = np * volone / cinfo[icell].weight / voltotal;
    npercell = static_cast<int> (ntarget);
    if (random->uniform() < ntarget-npercell) npercell++;
    nsingle += npercell;

    // generate particles from Maxwell-Boltzmann distribution
    
    if (mode == MAXWELL) {

      for (int m = 0; m < npercell; m++) {
        rn = random->uniform();

        isp = 0;
        while (cummulative[isp] < rn) isp++;
        ispecies = species[isp];

        x[0] = lo[0] + random->uniform() * (hi[0]-lo[0]);
        x[1] = lo[1] + random->uniform() * (hi[1]-lo[1]);
        x[2] = lo[2] + random->uniform() * (hi[2]-lo[2]);
        if (dimension == 2) x[2] = 0.0;

        vn = vscale[isp] * sqrt(-log(random->uniform()));
        vr = vscale[isp] * sqrt(-log(random->uniform()));
        theta1 = MY_2PI * random->uniform();
        theta2 = MY_2PI * random->uniform();

        v[0] = vstream[0] + vn*cos(theta1);
        v[1] = vstream[1] + vr*cos(theta2);
        v[2] = vstream[2] + vr*sin(theta2);


        erot = particle->erot(ispecies,temp_rot,random);
        evib = particle->evib(ispecies,temp_vib,random);
      
        id = MAXSMALLINT*random->uniform();

        particle->add_particle(id,ispecies,icell,x,v,erot,evib);
        
        if (nfix_update_custom)
          modify->update_custom(particle->nlocal-1,temp_thermal,
                                temp_rot,temp_vib,vstream);
      }
    }
        
    // generate particles from Chapman-Enskog distribution
    // NOTE: add C-E code here
      
    else if (mode == CHAPMAN) {

    }
  }
  
  ntotal += nsingle;
}

/* ----------------------------------------------------------------------
   return one-step or total count of particle insertions
------------------------------------------------------------------------- */

double FixEmitRegion::compute_vector(int i)
{
  double one,all;

  if (i == 0) one = nsingle;
  else one = ntotal;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}
