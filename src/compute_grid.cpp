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
#include "compute_grid.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace SPARTA_NS;

// user keywords

// SWS - adding the new output variable: NUMWI
enum{NUM,NUMWI,NRHO,NFRAC,MASS,MASSRHO,MASSFRAC,
     U,V,W,USQ,VSQ,WSQ,KE,TEMPERATURE,EROT,TROT,EVIB,TVIB,
     PXRHO,PYRHO,PZRHO,KERHO};

// internal accumulators
// SWS - add the keyword for the new accumulators: COUNT_WI and CELLCOUNTWI.
enum{COUNT,COUNT_WI,MASSSUM,MVX,MVY,MVZ,MVXSQ,MVYSQ,MVZSQ,MVSQ,
     ENGROT,ENGVIB,DOFROT,DOFVIB,CELLCOUNT,CELLMASS,CELLCOUNTWI,LASTSIZE};

// max # of quantities to accumulate for any user value

#define MAXACCUMULATE 2

/* ---------------------------------------------------------------------- */

ComputeGrid::ComputeGrid(SPARTA *sparta, int narg, char **arg) :
  Compute(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal compute grid command");

  int igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Compute grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  imix = particle->find_mixture(arg[3]);
  if (imix < 0) error->all(FLERR,"Compute grid mixture ID does not exist");
  ngroup = particle->mixture[imix]->ngroup;

  nvalue = narg - 4;
  value = new int[nvalue];
  cellcountwi = 0;  // SWS
  tvib_flag = 0;

  npergroup = cellmass = cellcount = 0;
  unique = new int[LASTSIZE];
  nmap = new int[nvalue];
  memory->create(map,ngroup*nvalue,MAXACCUMULATE,"grid:map");
  for (int i = 0; i < nvalue; i++) nmap[i] = 0;

  int ivalue = 0;
  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"n") == 0) {
      value[ivalue] = NUM;
      set_map(ivalue,COUNT);
    } else if (strcmp(arg[iarg],"sumwi") == 0) {  // SWS - wi summation output
      value[ivalue] = NUMWI;
      set_map(ivalue,COUNT_WI);
    } else if (strcmp(arg[iarg],"nrho") == 0) {
      value[ivalue] = NRHO;
      set_map(ivalue,COUNT_WI);  // SWS
      cellcountwi = 1;
    } else if (strcmp(arg[iarg],"nfrac") == 0) {
      value[ivalue] = NFRAC;
      set_map(ivalue,COUNT_WI);  // SWS
      set_map(ivalue,CELLCOUNTWI);  // SWS
      cellcountwi = 1;  // SWS
      cellcount = 1;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      value[ivalue] = MASS;
      set_map(ivalue,MASSSUM);
      set_map(ivalue,COUNT_WI); // SWS
      cellcountwi = 1;  // SWS
    } else if (strcmp(arg[iarg],"massrho") == 0) {
      value[ivalue] = MASSRHO;
      set_map(ivalue,MASSSUM);
    } else if (strcmp(arg[iarg],"massfrac") == 0) {
      value[ivalue] = MASSFRAC;
      set_map(ivalue,MASSSUM);
      set_map(ivalue,CELLMASS);
      cellmass = 1;
    } else if (strcmp(arg[iarg],"u") == 0) {
      value[ivalue] = U;
      set_map(ivalue,MVX);
      set_map(ivalue,MASSSUM);
    } else if (strcmp(arg[iarg],"v") == 0) {
      value[ivalue] = V;
      set_map(ivalue,MVY);
      set_map(ivalue,MASSSUM);
    } else if (strcmp(arg[iarg],"w") == 0) {
      value[ivalue] = W;
      set_map(ivalue,MVZ);
      set_map(ivalue,MASSSUM);
    } else if (strcmp(arg[iarg],"usq") == 0) {
      value[ivalue] = USQ;
      set_map(ivalue,MVXSQ);
      set_map(ivalue,MASSSUM);
    } else if (strcmp(arg[iarg],"vsq") == 0) {
      value[ivalue] = VSQ;
      set_map(ivalue,MVYSQ);
      set_map(ivalue,MASSSUM);
    } else if (strcmp(arg[iarg],"wsq") == 0) {
      value[ivalue] = WSQ;
      set_map(ivalue,MVZSQ);
      set_map(ivalue,MASSSUM);
    } else if (strcmp(arg[iarg],"ke") == 0) {
      value[ivalue] = KE;
      set_map(ivalue,MVSQ);
      set_map(ivalue,COUNT_WI);  // SWS
      cellcountwi = 1;  // SWS
    } else if (strcmp(arg[iarg],"temp") == 0) {
      value[ivalue] = TEMPERATURE;
      set_map(ivalue,MVSQ);
      set_map(ivalue,COUNT_WI);  // SWS
      cellcountwi = 1;  // SWS
    } else if (strcmp(arg[iarg],"erot") == 0) {
      value[ivalue] = EROT;
      set_map(ivalue,ENGROT);
      set_map(ivalue,COUNT_WI);  // SWS
      cellcountwi = 1;  // SWS
    } else if (strcmp(arg[iarg],"trot") == 0) {
      value[ivalue] = TROT;
      set_map(ivalue,ENGROT);
      set_map(ivalue,DOFROT);
    } else if (strcmp(arg[iarg],"evib") == 0) {
      value[ivalue] = EVIB;
      set_map(ivalue,ENGVIB);
      set_map(ivalue,COUNT_WI);  // SWS
      cellcountwi = 1;  // SWS
    } else if (strcmp(arg[iarg],"tvib") == 0) {
      value[ivalue] = TVIB;
      set_map(ivalue,ENGVIB);
      set_map(ivalue,DOFVIB);
      tvib_flag = 1;
    } else if (strcmp(arg[iarg],"pxrho") == 0) {
      value[ivalue] = PXRHO;
      set_map(ivalue,MVX);
    } else if (strcmp(arg[iarg],"pyrho") == 0) {
      value[ivalue] = PYRHO;
      set_map(ivalue,MVY);
    } else if (strcmp(arg[iarg],"pzrho") == 0) {
      value[ivalue] = PZRHO;
      set_map(ivalue,MVZ);
    } else if (strcmp(arg[iarg],"kerho") == 0) {
      value[ivalue] = KERHO;
      set_map(ivalue,MVSQ);
    } else error->all(FLERR,"Illegal compute grid command");

    ivalue++;
    iarg++;
  }

  // ntotal = total # of columns in tally array
  // reset_map() adjusts indices in initial map() using final npergroup
  // also adds columns to tally array for CELLCOUNT/CELLMASS

  ntotal = ngroup*npergroup;
  reset_map();

  per_grid_flag = 1;
  size_per_grid_cols = ngroup*nvalue;
  post_process_grid_flag = 1;

  nglocal = 0;
  vector_grid = NULL;
  tally = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeGrid::~ComputeGrid()
{
  if (copymode) return;

  delete [] value;
  delete [] unique;

  delete [] nmap;
  memory->destroy(map);

  memory->destroy(vector_grid);
  memory->destroy(tally);
}

/* ---------------------------------------------------------------------- */

void ComputeGrid::init()
{
  if (ngroup != particle->mixture[imix]->ngroup)
    error->all(FLERR,"Number of groups in compute grid mixture has changed");

  if (tvib_flag && particle->find_custom((char *) "vibmode") >= 0)
    if (comm->me == 0)
      error->warning(FLERR,"Using compute grid tvib with fix vibmode may give "
                     "incorrect temperature, use compute tvib/grid instead");

  eprefactor = 0.5*update->mvv2e;
  tprefactor = update->mvv2e / (3.0*update->boltz);
  rvprefactor = 2.0*update->mvv2e / update->boltz;

  reallocate();
}

/* ---------------------------------------------------------------------- */

void ComputeGrid::compute_per_grid()
{
  invoked_per_grid = update->ntimestep;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::Species *species = particle->species;
  Particle::OnePart *particles = particle->particles;
  int *s2g = particle->mixture[imix]->species2group;
  int nlocal = particle->nlocal;

  int i,j,k,m,ispecies,igroup,icell;
  double mass;
  double *v,*vec;
  double specwt;  // SWS

  // zero all accumulators - could do this with memset()

  for (i = 0; i < nglocal; i++)
    for (j = 0; j < ntotal; j++)
      tally[i][j] = 0.0;

  // loop over all particles, skip species not in mixture group
  // skip cells not in grid group
  // perform all tallies needed for each particle
  // depends on its species group and the user-requested values

  for (i = 0; i < nlocal; i++) {
    ispecies = particles[i].ispecies;
    igroup = s2g[ispecies];
    if (igroup < 0) continue;
    icell = particles[i].icell;
    if (!(cinfo[icell].mask & groupbit)) continue;

    mass = species[ispecies].mass;
    v = particles[i].v;

    vec = tally[icell];
    // ========================================================================
    // SWS - Accumulate the number of physical particles: sum of the weight, and 
    // the corresponding mass: sum of the product mass * weight.
    // ========================================================================
    // Baseline code:
    //~ if (cellmass) vec[cellmass] += mass;
    //~ if (cellcount) vec[cellcount] += 1.0;
    // Modified code:
    specwt = species[ispecies].specwt;
    if (cellcountwi) vec[cellcountwi] += specwt;
    if (cellmass) vec[cellmass] += mass*specwt;
    if (cellcount) vec[cellcount] += 1.0;

    // loop has all possible values particle needs to accumulate
    // subset defined by user values are indexed by accumulate vector

    k = igroup*npergroup;

    for (m = 0; m < npergroup; m++) {
      switch (unique[m]) {
      case COUNT:
        vec[k++] += 1.0;
        break;
      case COUNT_WI:         // SWS
        vec[k++] += specwt;  // SWS
        break;               // SWS
      case MASSSUM:
        vec[k++] += mass*specwt;  // SWS
        break;
      case MVX:
        vec[k++] += mass*specwt*v[0];  // SWS
        break;
      case MVY:
        vec[k++] += mass*specwt*v[1];  // SWS
        break;
      case MVZ:
        vec[k++] += mass*specwt*v[2];  // SWS
        break;
      case MVXSQ:
        vec[k++] += mass*specwt*v[0]*v[0];  // SWS
        break;
      case MVYSQ:
        vec[k++] += mass*specwt*v[1]*v[1];  // SWS
        break;
      case MVZSQ:
        vec[k++] += mass*specwt*v[2]*v[2];  // SWS
        break;
      case MVSQ:
        vec[k++] += mass*specwt * (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);  // SWS
        break;
      case ENGROT:
        vec[k++] += specwt * particles[i].erot;  // SWS
        break;
      case ENGVIB:
        vec[k++] += specwt * particles[i].evib;  // SWS
        break;
      case DOFROT:
        vec[k++] += specwt * species[ispecies].rotdof;  // SWS
        break;
      case DOFVIB:
        vec[k++] += specwt * species[ispecies].vibdof;  // SWS
        break;
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

int ComputeGrid::query_tally_grid(int index, double **&array, int *&cols)
{
  index--;
  int ivalue = index % nvalue;
  array = tally;
  cols = map[index];
  return nmap[ivalue];
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

void ComputeGrid::post_process_grid(int index, int nsample,
                                    double **etally, int *emap,
                                    double *vec, int nstride)
{
  index--;
  int ivalue = index % nvalue;

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

  // compute normalized final value for each grid cell

  switch (value[ivalue]) {

  case NUM:
    {
      int count = emap[0];
      for (int icell = lo; icell < hi; icell++) {
        vec[k] = etally[icell][count] / nsample;
        k += nstride;
      }
      break;
    }
  case NUMWI:  // SWS
    {
      int count_wi = emap[0];
      for (int icell = lo; icell < hi; icell++) {
        vec[k] = etally[icell][count_wi] / nsample;
        k += nstride;
      }
      break;
    }
  case MASS:
    {
      double norm;
      int mass = emap[0];
      int count_wi = emap[1];  // SWS
      for (int icell = lo; icell < hi; icell++) {
        // Baseline code:
        //~ norm = etally[icell][count];
        // SWS - modifified code:
        norm = etally[icell][count_wi];
        if (norm == 0.0) vec[k] = 0.0;
        else vec[k] = etally[icell][mass] / norm;
        k += nstride;
      }
      break;
    }

  case NRHO:
    {
      double wt;
      double fnum = update->fnum;
      Grid::ChildInfo *cinfo = grid->cinfo;

      double norm;
      int count_wi = emap[0]; // SWS
      for (int icell = lo; icell < hi; icell++) {
        norm = cinfo[icell].volume;
        if (norm == 0.0) vec[k] = 0.0;
        else {
          wt = fnum * cinfo[icell].weight / norm;     
          vec[k] = wt * etally[icell][count_wi] / nsample;  // SWS
        }
        k += nstride;
      }
      break;
    }

  case MASSRHO:
    {
      double wt;
      double fnum = update->fnum;
      Grid::ChildInfo *cinfo = grid->cinfo;

      double norm;
      int mass = emap[0];
      for (int icell = lo; icell < hi; icell++) {
        norm = cinfo[icell].volume;
        if (norm == 0.0) vec[k] = 0.0;
        else {
          wt = fnum * cinfo[icell].weight / norm;
          vec[k] = wt * etally[icell][mass] / nsample;
        }
        k += nstride;
      }
      break;
    }

  case NFRAC:
  case MASSFRAC:
    {
      double norm;
      int count_or_mass = emap[0];
      int cell_count_or_mass = emap[1];
      for (int icell = lo; icell < hi; icell++) {
        norm = etally[icell][cell_count_or_mass];
        if (norm == 0.0) vec[k] = 0.0;
        else vec[k] = etally[icell][count_or_mass] / norm;
        k += nstride;
      }
      break;
    }

  case U:
  case V:
  case W:
  case USQ:
  case VSQ:
  case WSQ:
    {
      double norm;
      int velocity = emap[0];
      int mass = emap[1];
      for (int icell = lo; icell < hi; icell++) {
        norm = etally[icell][mass];
        if (norm == 0.0) vec[k] = 0.0;
        else vec[k] = etally[icell][velocity] / norm;
        k += nstride;
      }
      break;
    }

  case KE:
    {
      double norm;
      int mvsq = emap[0];
      int count_wi = emap[1];  // SWS
      for (int icell = lo; icell < hi; icell++) {
        norm = etally[icell][count_wi];  // SWS
        if (norm == 0.0) vec[k] = 0.0;
        else vec[k] = eprefactor * etally[icell][mvsq] / norm;
        k += nstride;
      }
      break;
    }

  case TEMPERATURE:
    {
      double norm;
      int mvsq = emap[0];
      int count_wi = emap[1];  // SWS
      for (int icell = lo; icell < hi; icell++) {
        norm = etally[icell][count_wi];  // SWS
        if (norm == 0.0) vec[k] = 0.0;
        else vec[k] = tprefactor * etally[icell][mvsq] / norm;
        k += nstride;
      }
      break;
    }

  case EROT:
  case EVIB:
    {
      double norm;
      int eng = emap[0];
      int count_wi = emap[1];  // SWS
      for (int icell = lo; icell < hi; icell++) {
        norm = etally[icell][count_wi];  // SWS
        if (norm == 0.0) vec[k] = 0.0;
        else vec[k] = etally[icell][eng] / norm;
        k += nstride;
      }
      break;
    }

  case TROT:
  case TVIB:
    {
      double norm;
      int eng = emap[0];
      int dof = emap[1];
      for (int icell = lo; icell < hi; icell++) {
        norm = etally[icell][dof];
        if (norm == 0.0) vec[k] = 0.0;
        else vec[k] = rvprefactor * etally[icell][eng] / norm;
        k += nstride;
      }
      break;
    }

  case PXRHO:
  case PYRHO:
  case PZRHO: 
    {
      double wt;
      double fnum = update->fnum;
      Grid::ChildInfo *cinfo = grid->cinfo;

      double norm;
      int mom = emap[0];
      for (int icell = lo; icell < hi; icell++) {
        norm = cinfo[icell].volume;
        if (norm == 0.0) vec[k] = 0.0;
        else {
          wt = fnum * cinfo[icell].weight / norm;
          vec[k] = wt * etally[icell][mom] / nsample;
        }
        k += nstride;
      }
      break;
    }

  case KERHO: 
    {
      double wt;
      double fnum = update->fnum;
      Grid::ChildInfo *cinfo = grid->cinfo;

      double norm;
      int ke = emap[0];
      for (int icell = lo; icell < hi; icell++) {
        norm = cinfo[icell].volume;
        if (norm == 0.0) vec[k] = 0.0;
        else {
          wt = fnum * cinfo[icell].weight / norm;
          vec[k] = eprefactor * wt * etally[icell][ke] / nsample;
        }
        k += nstride;
      }
      break;
    }
  }
}

/* ----------------------------------------------------------------------
   add a tally quantity to all groups for ivalue
   also add it to unique list if first time this name is used
   name = name of tally quantity from enum{} at top of file
   nmap[i] = # of tally quantities for user value I
   map[i][k] = index of Kth tally quantity for output value I
   npergroup = length of unique list
   do not add CELLCOUNT/CELLMASS to unique, since not a pergroup tally
------------------------------------------------------------------------- */

void ComputeGrid::set_map(int ivalue, int name)
{
  // index = loc of name in current unique list if there, else npergroup

  int index = 0;
  for (index = 0; index < npergroup; index++)
    if (unique[index] == name) break;

  // if name = CELLCOUNT/CELLMASS, just set index to name for now
  // if name is not already in unique, add it and increment npergroup

  if (name == CELLCOUNT || name == CELLMASS || name == CELLCOUNTWI) index = name; // SWS

  else if (index == npergroup) {
    index = npergroup;
    unique[npergroup++] = name;
  }

  // add index to map and nmap for all groups
  // will add group offset in reset_map()

  for (int igroup = 0; igroup < ngroup; igroup++)
    map[igroup*nvalue+ivalue][nmap[ivalue]] = index;
  nmap[ivalue]++;
}

/* ----------------------------------------------------------------------
   increment ntotal = # of tally quantities for CELLCOUNT/CELLMASS
   reset map indices to reflect final npergroup = unique quantities/group
------------------------------------------------------------------------- */

void ComputeGrid::reset_map()
{
  if (cellcount) cellcount = ntotal++;
  if (cellmass) cellmass = ntotal++;
  if (cellcountwi) cellcountwi = ntotal++;  // SWS

  for (int i = 0; i < ngroup*nvalue; i++) {
    int igroup = i / nvalue;
    int ivalue = i % nvalue;
    for (int k = 0; k < nmap[ivalue]; k++) {
      if (map[i][k] == CELLCOUNT) map[i][k] = cellcount;
      else if (map[i][k] == CELLMASS) map[i][k] = cellmass;
      else if (map[i][k] == CELLCOUNTWI) map[i][k] = cellcountwi;  // SWS
      else map[i][k] += igroup*npergroup;
    }
  }
}

/* ----------------------------------------------------------------------
   reallocate data storage if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void ComputeGrid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memory->destroy(vector_grid);
  memory->destroy(tally);
  nglocal = grid->nlocal;
  memory->create(vector_grid,nglocal,"grid:vector_grid");
  memory->create(tally,nglocal,ntotal,"grid:tally");
}

/* ----------------------------------------------------------------------
   memory usage of local grid-based data
------------------------------------------------------------------------- */

bigint ComputeGrid::memory_usage()
{
  bigint bytes;
  bytes = nglocal * sizeof(double);
  bytes = ntotal*nglocal * sizeof(double);
  return bytes;
}
