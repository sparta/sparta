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

#include "string.h"
#include "compute_tvib_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{NONE,COUNT,MASSWT,DOF};

/* ---------------------------------------------------------------------- */

ComputeTvibGrid::ComputeTvibGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute tvib/grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute tvib/grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute tvib/grid mixture ID does not exist");

  // optional args

  modeflag = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mode") == 0) {
      modeflag = 2;
      iarg++;
    } else error->all(FLERR,"Illegal compute tvib/grid command");
  }

  // convert modeflag = 0,1,2 and error check
  // possible vibstyle settings = NONE,DISCRETE,SMOOTH
  // per-particle fix vibmode vectors may exist or not
  //   can only exist if DISCRETE, but may not exist if DISCRETE
  // user mode setting can be on or off
  //   if off and no vib vectors exist -> modeflag = 0
  //   if off and vib vectors exist -> modeflag = 1
  //   if on, vib vectors must exist -> modeflag = 2

  if (modeflag == 0 && particle->find_custom((char *) "vibmode") >= 0)
    modeflag = 1;

  if (modeflag == 2 && particle->find_custom((char *) "vibmode") < 0)
    error->all(FLERR,"Cannot use compute tvib/grid mode without "
               "fix vibmode defined");

  // maxmode = max # of vib modes for any DISCRETE species
  // not used for modeflag = 0

  if (modeflag != 0) {
    maxmode = particle->maxvibmode;
    if (maxmode == 0)
      error->all(FLERR,"No species in compute tvib/grid has vibrational modes");
  }

  // ngroup = # of groups in mixture
  // ntally = # of tally values per grid cell in tally array
  // size_per_grid_cols = # of output columns from this compute

  ngroup = particle->mixture[imix]->ngroup;
  mixspecies = particle->mixture[imix]->nspecies;
  nspecies = particle->nspecies;

  if (modeflag == 0) ntally = 2*mixspecies;
  else ntally = 2*mixspecies*maxmode;

  per_grid_flag = 1;
  if (modeflag != 2) size_per_grid_cols = ngroup;
  else size_per_grid_cols = ngroup * maxmode;
  post_process_grid_flag = 1;

  // allocate and initialize nmap,map
  //   nmax = max # of species in any group
  // for modeflag = 0: also tspecies,s2t,t2s
  // for modeflag = 1 and 2: also tspecies_mode,s2t_mode,t2s_mode
  // must first set groupsize, groupspecies by mixture->init()
  // 2 tally quantities per species and per mode for every group

  particle->mixture[imix]->init();
  int *groupsize = particle->mixture[imix]->groupsize;
  int **groupspecies = particle->mixture[imix]->groupspecies;

  int nmax = 0;
  for (int i = 0; i < ngroup; i++)
    nmax = MAX(nmax,groupsize[i]);

  // modeflag = 0
  // Ngroup outputs, 2*Nspecies tallies per output

  if (modeflag == 0) {
    nmap = new int[ngroup];
    memory->create(map,ngroup,2*nmax,"tvib/grid:map");

    tspecies = new double[nspecies];
    s2t = new int[nspecies];
    t2s = new int[ntally];

    for (int isp = 0; isp < nspecies; isp++) s2t[isp] = -1;

    int itally = 0;
    for (int igroup = 0; igroup < ngroup; igroup++) {
      nmap[igroup] = 2*groupsize[igroup];
      for (int n = 0; n < groupsize[igroup]; n++) {
        map[igroup][2*n] = itally;
        map[igroup][2*n+1] = itally+1;
        s2t[groupspecies[igroup][n]] = itally;
        t2s[itally] = groupspecies[igroup][n];
        t2s[itally+1] = groupspecies[igroup][n];
        itally += 2;
      }
    }

  // modeflag = 1
  // Ngroup outputs, 2*Nspecies*Nmode tallies per output

  } else if (modeflag == 1) {
    nmap = new int[ngroup];
    memory->create(map,ngroup,2*nmax*maxmode,"tvib/grid:map");

    memory->create(tspecies_mode,nspecies,maxmode,"tvib/grid:tspecies_mode");
    memory->create(s2t_mode,nspecies,maxmode,"tvib/grid:s2t_mode");
    t2s_mode = new int[ntally];

    for (int isp = 0; isp < nspecies; isp++)
      for (int imode = 0; imode < maxmode; imode++)
        s2t_mode[isp][imode] = -1;

    int itally = 0;
    for (int igroup = 0; igroup < ngroup; igroup++) {
      nmap[igroup] = 2*groupsize[igroup]*maxmode;
      for (int n = 0; n < groupsize[igroup]; n++) {
        for (int imode = 0; imode < maxmode; imode++) {
          map[igroup][2*n*maxmode + 2*imode] = itally;
          map[igroup][2*n*maxmode + 2*imode + 1] = itally+1;
          s2t_mode[groupspecies[igroup][n]][imode] = itally;
          t2s_mode[itally] = groupspecies[igroup][n];
          t2s_mode[itally+1] = groupspecies[igroup][n];
          itally += 2;
        }
      }
    }

  // modeflag = 2
  // Ngroup*Nmode outputs, 2*Nspecies*Nmode tallies per output
  // tally count and list of tally columns are duplicated
  //   Nmode times for each group

  } else if (modeflag == 2) {
    nmap = new int[ngroup*maxmode];
    memory->create(map,ngroup*maxmode,2*nmax*maxmode,"tvib/grid:map");

    memory->create(tspecies_mode,nspecies,maxmode,"tvib/grid:tspecies_mode");
    memory->create(s2t_mode,nspecies,maxmode,"tvib/grid:s2t_mode");
    t2s_mode = new int[ntally];

    for (int isp = 0; isp < nspecies; isp++)
      for (int imode = 0; imode < maxmode; imode++)
        s2t_mode[isp][imode] = -1;

    int itally = 0;
    for (int igroup = 0; igroup < ngroup; igroup++) {
      nmap[igroup*maxmode] = 2*groupsize[igroup]*maxmode;
      for (int n = 0; n < groupsize[igroup]; n++) {
        for (int imode = 0; imode < maxmode; imode++) {
          map[igroup*maxmode][2*n*maxmode + 2*imode] = itally;
          map[igroup*maxmode][2*n*maxmode + 2*imode + 1] = itally+1;
          s2t_mode[groupspecies[igroup][n]][imode] = itally;
          t2s_mode[itally] = groupspecies[igroup][n];
          t2s_mode[itally+1] = groupspecies[igroup][n];
          itally += 2;
        }
      }
    }

    // copy each nmap[] value and map[][*] row into all other modes

    for (int igroup = 0; igroup < ngroup; igroup++) {
      for (int imode = 1; imode < maxmode; imode++)
        nmap[igroup*maxmode + imode] = nmap[igroup*maxmode];

      int ncol = 2*groupsize[igroup]*maxmode;
      for (int imode = 1; imode < maxmode; imode++)
        for (int icol = 0; icol < ncol; icol++)
          map[igroup*maxmode + imode][icol] = map[igroup*maxmode][icol];
    }
  }

  nglocal = 0;
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeTvibGrid::~ComputeTvibGrid()
{
  if (copymode) return;

  delete [] nmap;
  memory->destroy(map);

  if (modeflag == 0) {
    delete [] tspecies;
    delete [] s2t;
    delete [] t2s;
  } else {
    memory->destroy(tspecies_mode);
    memory->destroy(s2t_mode);
    delete [] t2s_mode;
  }

  memory->destroy(vector_grid);
  memory->destroy(tally);
}

/* ---------------------------------------------------------------------- */

void ComputeTvibGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute tvib/grid "
               "mixture has changed");
  if (mixspecies != particle->mixture[imix]->nspecies)
    error->all(FLERR,"Number of species in compute tvib/grid "
               "mixture has changed");
  if (nspecies != particle->nspecies)
    error->all(FLERR,"Number of total species in compute tvib/grid "
               "has changed");

  if (modeflag > 0 && particle->find_custom((char *) "vibmode") < 0)
    error->all(FLERR,"Cannot use compute tvib/grid mode without "
               "fix vibmode defined");

  int *groupsize = particle->mixture[imix]->groupsize;
  for (int i = 0; i < ngroup; i++)
    if ((modeflag == 0 && 2*groupsize[i] != nmap[i]) ||
        (modeflag > 0  && 2*groupsize[i]*maxmode != nmap[i]))
      error->all(FLERR,"Number of species in compute tvib/grid "
                 "group has changed");

  if (modeflag > 0) index_vibmode = particle->find_custom((char *) "vibmode");

  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeTvibGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  int *s2g = particle->mixture[imix]->species2group;
  double boltz = update->boltz;
  int nlocal = particle->nlocal;

  int i,j,ispecies,igroup,icell,imode,nmode;

  // zero all accumulators - could do this with memset()

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntally; j++)
      tally[i][j] = 0.0;

  // loop over all particles, skip species not in mixture group
  // mode = 0: tally vib eng and count for each species
  // mode >= 1: tally vib level and count for each species and each vib mode

  if (modeflag == 0) {
    for (i = 0; i < nlocal; i++) {
      ispecies = particles[i].ispecies;
      if (!species[ispecies].vibdof) continue;
      igroup = s2g[ispecies];
      if (igroup < 0) continue;
      icell = particles[i].icell;
      if (!(cinfo[icell].mask & groupbit)) continue;

      j = s2t[ispecies];
      tally[icell][j] += particles[i].evib;
      tally[icell][j+1] += 1.0;
    }

  } else if (modeflag >= 1) {
    int **vibmode =
      particle->eiarray[particle->ewhich[index_vibmode]];

    for (i = 0; i < nlocal; i++) {
      ispecies = particles[i].ispecies;
      if (!species[ispecies].vibdof) continue;
      igroup = s2g[ispecies];
      if (igroup < 0) continue;
      icell = particles[i].icell;
      if (!(cinfo[icell].mask & groupbit)) continue;

      // tally only the modes this species has

      nmode = particle->species[ispecies].nvibmode;
      for (imode = 0; imode < nmode; imode++) {
        j = s2t_mode[ispecies][imode];
        if (nmode > 1) tally[icell][j] += vibmode[i][imode];
        else tally[icell][j] +=
               particles[i].evib / (boltz*species[ispecies].vibtemp[0]);
        tally[icell][j+1] += 1.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   query info about internal tally array for this compute
   index = which column of output (0 for vec, 1 to N for array)
   return # of tally quantities for this index
   also return array = ptr to tally array
   also return cols = ptr to list of columns in tally for this index
------------------------------------------------------------------------- */

int ComputeTvibGrid::query_tally_grid(int index, double **&array, int *&cols)
{
  index--;
  array = tally;
  cols = map[index];
  return nmap[index];
}

/* ----------------------------------------------------------------------
   tally accumulated info to compute final normalized values
   index = which column of output (0 for vec, 1 to N for array)
   for etally = NULL:
     use internal tallied info for single timestep, set nsample = 1
     compute values for all grid cells
       store results in vector_grid with nstride = 1 (single col of array_grid)
   for etally = ptr to caller array:
     use external tallied info for many timesteps
     nsample = additional normalization factor used by some values
     emap = list of etally columns to use, # of columns determined by index
     store results in caller's vec, spaced by nstride
   if norm = 0.0, set result to 0.0 directly so do not divide by 0.0
------------------------------------------------------------------------- */

void ComputeTvibGrid::post_process_grid(int index, int nsample,
                                        double **etally, int *emap,
                                        double *vec, int nstride)
{
  index--;

  int lo = 0;
  int hi = nglocal;
  int k = 0;

  if (!etally) {
    nsample = 1;
    etally = tally;
    emap = map[index];
    vec = vector_grid;
    nstride = 1;
  }

  // compute normalized single Tgroup value for each grid cell
  // compute Ibar and Tspecies for each species in the requested group
  // ditto for individual vibrational modes if modeflag = 1 or 2
  // loop over species/modes in group to compute normalized Tgroup

  Particle::Species *species = particle->species;

  int isp,evib,count,ispecies,imode;
  double theta,ibar,numer,denom;
  double boltz = update->boltz;

  // modeflag = 0, no vib modes exist
  // nsp = # of species in the group
  // inputs: 2*nsp tallies
  // output: Tgroup = weighted sum over all Tsp for species in group

  if (modeflag == 0) {
    int nsp = nmap[index] / 2;

    for (int icell = lo; icell < hi; icell++) {
      evib = emap[0];
      count = evib+1;
      for (isp = 0; isp < nsp; isp++) {
        ispecies = t2s[evib-emap[0]];
        theta = species[ispecies].vibtemp[0];
        if (theta == 0.0 || etally[icell][count] == 0.0) {
          tspecies[isp] = 0.0;
          evib += 2;
          count = evib+1;
          continue;
        }
        ibar = etally[icell][evib] / (etally[icell][count] * boltz * theta);
        if (ibar == 0.0) {
          tspecies[isp] = 0.0;
          evib += 2;
          count = evib+1;
          continue;
        }
        tspecies[isp] = theta / (log(1.0 + 1.0/ibar));
        //denom = boltz * etally[icell][count] * ibar * log(1.0 + 1.0/ibar);
        //tspecies[isp] = etally[icell][evib] / denom;
        evib += 2;
        count = evib+1;
      }

      // loop over species in group to accumulate numerator & denominator

      numer = denom = 0.0;
      count = emap[1];
      for (isp = 0; isp < nsp; isp++) {
        numer += tspecies[isp]*etally[icell][count];
        denom += etally[icell][count];
        count += 2;
      }

      if (denom == 0.0) vec[k] = 0.0;
      else vec[k] = numer/denom;
      k += nstride;
    }

  // modeflag = 1, vib modes exist
  // Tgroup = weighted sum over all Tsp and modes for species in group
  // nsp = # of species in the group
  // maxmode = max # of modes for any species (unused values are zero)
  // inputs: 2*nsp*maxmode tallies

  } else if (modeflag == 1) {
    int nsp = nmap[index] / maxmode / 2;
    int **vibmode =
      particle->eiarray[particle->ewhich[index_vibmode]];

    for (int icell = lo; icell < hi; icell++) {
      evib = emap[0];
      count = evib+1;
      for (isp = 0; isp < nsp; isp++) {
        ispecies = t2s_mode[evib-emap[0]];
        for (imode = 0; imode < maxmode; imode++) {
          theta = species[ispecies].vibtemp[imode];
          if (theta == 0.0 || etally[icell][count] == 0.0) {
            tspecies_mode[isp][imode] = 0.0;
            evib += 2;
            count = evib+1;
            continue;
          }
          ibar = etally[icell][evib] / etally[icell][count];
          if (ibar == 0.0) {
            tspecies_mode[isp][imode] = 0.0;
            evib += 2;
            count = evib+1;
            continue;
          }
          tspecies_mode[isp][imode] = theta / (log(1.0 + 1.0/ibar));
          //denom = boltz * etally[icell][count] * ibar * log(1.0 + 1.0/ibar);
          //tspecies_mode[isp][imode] = etally[icell][evib] / denom;
          evib += 2;
          count = evib+1;
        }
      }

      // loop over species in group and all their modes
      // to accumulate numerator & denominator

      numer = denom = 0.0;
      count = emap[1];
      for (isp = 0; isp < nsp; isp++) {
        for (imode = 0; imode < maxmode; imode++) {
          numer += tspecies_mode[isp][imode]*etally[icell][count];
          denom += etally[icell][count];
          count += 2;
        }
      }

      if (denom == 0.0) vec[k] = 0.0;
      else vec[k] = numer/denom;
      k += nstride;
    }

  // modeflag = 2, vib modes exist
  // Tgroup = weighted sum over all Tsp and single mode for species in group
  // nsp = # of species in the group
  // imode = single mode correpsonding to caller index
  // inputs: 2*nsp tallies strided by maxmode

  } else if (modeflag == 2) {
    int nsp = nmap[index] / maxmode / 2;
    imode = index % maxmode;
    int **vibmode =
      particle->eiarray[particle->ewhich[index_vibmode]];

    for (int icell = lo; icell < hi; icell++) {
      evib = emap[2*imode];
      count = evib+1;
      for (isp = 0; isp < nsp; isp++) {
        ispecies = t2s_mode[evib-emap[0]];
        theta = species[ispecies].vibtemp[imode];
        if (theta == 0.0 || etally[icell][count] == 0.0) {
          tspecies_mode[isp][imode] = 0.0;
          evib += 2*maxmode;
          count = evib+1;
          continue;
        }
        ibar = etally[icell][evib] / etally[icell][count];
        if (ibar == 0.0) {
          tspecies_mode[isp][imode] = 0.0;
          evib += 2*maxmode;
          count = evib+1;
          continue;
        }
        tspecies_mode[isp][imode] = theta / (log(1.0 + 1.0/ibar));
        //denom = boltz * etally[icell][count] * ibar * log(1.0 + 1.0/ibar);
        //tspecies_mode[isp][imode] = etally[icell][evib] / denom;
        evib += 2*maxmode;
        count = evib+1;
      }

      // loop over species in group and single mode for each species
      // to accumulate numerator & denominator

      numer = denom = 0.0;
      count = emap[2*imode+1];
      for (isp = 0; isp < nsp; isp++) {
        numer += tspecies_mode[isp][imode]*etally[icell][count];
        denom += etally[icell][count];
        count += 2*maxmode;
      }

      if (denom == 0.0) vec[k] = 0.0;
      else vec[k] = numer/denom;
      k += nstride;
    }
  }
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeTvibGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memory->destroy(vector_grid);
  memory->destroy(tally);
  nglocal = grid->nlocal;
  memory->create(vector_grid,nglocal,"tvib/grid:vector_grid");
  memory->create(tally,nglocal,ntally,"tvib/grid:tally");
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based data
------------------------------------------------------------------------- */

bigint ComputeTvibGrid::memory_usage()
{
  bigint bytes;
  bytes = nglocal * sizeof(double);
  bytes = ntally*nglocal * sizeof(double);
  return bytes;
}
