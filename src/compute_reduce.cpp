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
#include "stdlib.h"
#include "compute_reduce.h"
#include "update.h"
#include "domain.h"
#include "particle.h"
#include "mixture.h"
#include "grid.h"
#include "surf.h"
#include "modify.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{SUM,SUMSQ,MINN,MAXX,AVE,AVESQ,SUMAREA,AVEAREA};
enum{X,V,KE,EROT,EVIB,COMPUTE,FIX,VARIABLE};
enum{PARTICLE,GRID,SURF};

#define INVOKED_PER_PARTICLE 8
#define INVOKED_PER_GRID 16
#define INVOKED_PER_SURF 32

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

ComputeReduce::ComputeReduce(SPARTA *spa, int narg, char **arg) :
  Compute(spa, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute reduce command");

  if (strcmp(arg[2],"sum") == 0) mode = SUM;
  else if (strcmp(arg[2],"sumsq") == 0) mode = SUMSQ;
  else if (strcmp(arg[2],"min") == 0) mode = MINN;
  else if (strcmp(arg[2],"max") == 0) mode = MAXX;
  else if (strcmp(arg[2],"ave") == 0) mode = AVE;
  else if (strcmp(arg[2],"avesq") == 0) mode = AVESQ;
  else if (strcmp(arg[2],"sum-area") == 0) mode = SUMAREA;
  else if (strcmp(arg[2],"ave-area") == 0) mode = AVEAREA;
  else error->all(FLERR,"Illegal compute reduce command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // expand args if any have wildcard character "*"

  int expand = 0;
  char **earg;
  int nargnew = input->expand_args(narg-3,&arg[3],1,earg);

  if (earg != &arg[3]) expand = 1;
  arg = earg;

  // parse remaining values until one isn't recognized

  which = new int[nargnew];
  argindex = new int[nargnew];
  flavor = new int[nargnew];
  ids = new char*[nargnew];
  value2index = new int[nargnew];
  nvalues = 0;

  int iarg = 0;
  while (iarg < nargnew) {
    ids[nvalues] = NULL;

    if (strcmp(arg[iarg],"x") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"y") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"z") == 0) {
      which[nvalues] = X;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"vx") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"ke") == 0) {
      which[nvalues] = KE;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"erot") == 0) {
      which[nvalues] = EROT;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"evib") == 0) {
      which[nvalues] = EVIB;
      argindex[nvalues++] = 0;

    } else if (strncmp(arg[iarg],"c_",2) == 0 ||
               strncmp(arg[iarg],"f_",2) == 0 ||
               strncmp(arg[iarg],"v_",2) == 0) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal compute reduce command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else break;

    iarg++;
  }

  // optional args

  replace = new int[nvalues];
  for (int i = 0; i < nvalues; i++) replace[i] = -1;
  subsetID = NULL;

  while (iarg < nargnew) {
    if (strcmp(arg[iarg],"replace") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal compute reduce command");
      if (mode != MINN && mode != MAXX)
        error->all(FLERR,"Compute reduce replace requires min or max mode");
      int col1 = atoi(arg[iarg+1]) - 1;
      int col2 = atoi(arg[iarg+2]) - 1;
      if (col1 < 0 || col1 >= nvalues || col2 < 0 || col2 >= nvalues)
        error->all(FLERR,"Illegal compute reduce command");
      if (col1 == col2) error->all(FLERR,"Illegal compute reduce command");
      if (replace[col1] >= 0 || replace[col2] >= 0)
        error->all(FLERR,"Invalid replace values in compute reduce");
      replace[col1] = col2;
      iarg += 3;
    } else if (strcmp(arg[iarg],"subset") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute reduce command");
      int n = strlen(arg[iarg+1]) + 1;
      subsetID = new char[n];
      strcpy(subsetID,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute reduce command");
  }

  // delete replace if not set

  int flag = 0;
  for (int i = 0; i < nvalues; i++)
    if (replace[i] >= 0) flag = 1;
  if (!flag) {
    delete [] replace;
    replace = NULL;
  }

  // if wildcard expansion occurred, free earg memory from expand_args()

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete [] earg[i];
    memory->sfree(earg);
  }

  // setup and error check
  // cannot use per-surf compute since data not yet summed across surfs

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == X || which[i] == V) flavor[i] = PARTICLE;

    else if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute reduce does not exist");
      if (modify->compute[icompute]->per_particle_flag) {
        flavor[i] = PARTICLE;
        if (argindex[i] == 0 &&
            modify->compute[icompute]->size_per_particle_cols != 0)
          error->all(FLERR,"Compute reduce compute does not "
                     "calculate a per-particle vector");
        if (argindex[i] &&
            modify->compute[icompute]->size_per_particle_cols == 0)
          error->all(FLERR,"Compute reduce compute does not "
                     "calculate a per-particle array");
        if (argindex[i] &&
            argindex[i] > modify->compute[icompute]->size_per_particle_cols)
          error->all(FLERR,
                     "Compute reduce compute array is accessed out-of-range");
      } else if (modify->compute[icompute]->per_grid_flag) {
        flavor[i] = GRID;
        if (argindex[i] == 0 &&
            modify->compute[icompute]->size_per_grid_cols != 0)
          error->all(FLERR,"Compute reduce compute does not "
                     "calculate a per-grid vector");
        if (argindex[i] && modify->compute[icompute]->size_per_grid_cols == 0)
          error->all(FLERR,"Compute reduce compute does not "
                     "calculate a per-grid array");
        if (argindex[i] &&
            argindex[i] > modify->compute[icompute]->size_per_grid_cols)
          error->all(FLERR,
                     "Compute reduce compute array is accessed out-of-range");
      } else error->all(FLERR,"Compute reduce compute calculates "
                        "global or surf values");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute reduce does not exist");
      if (modify->fix[ifix]->per_particle_flag) {
        flavor[i] = PARTICLE;
        if (argindex[i] == 0 &&
            modify->fix[ifix]->size_per_particle_cols != 0)
          error->all(FLERR,"Compute reduce fix does not "
                     "calculate a per-particle vector");
        if (argindex[i] && modify->fix[ifix]->size_per_particle_cols == 0)
          error->all(FLERR,"Compute reduce fix does not "
                     "calculate a per-particle array");
        if (argindex[i] &&
            argindex[i] > modify->fix[ifix]->size_per_particle_cols)
          error->all(FLERR,"Compute reduce fix array is accessed out-of-range");
      } else if (modify->fix[ifix]->per_grid_flag) {
        flavor[i] = GRID;
        if (argindex[i] == 0 &&
            modify->fix[ifix]->size_per_grid_cols != 0)
          error->all(FLERR,"Compute reduce fix does not "
                     "calculate a per-grid vector");
        if (argindex[i] && modify->fix[ifix]->size_per_grid_cols == 0)
          error->all(FLERR,"Compute reduce fix does not "
                     "calculate a per-grid array");
        if (argindex[i] &&
            argindex[i] > modify->fix[ifix]->size_per_grid_cols)
          error->all(FLERR,"Compute reduce fix array is accessed out-of-range");
      } else if (modify->fix[ifix]->per_surf_flag) {
        flavor[i] = SURF;
        if (argindex[i] == 0 &&
            modify->fix[ifix]->size_per_surf_cols != 0)
          error->all(FLERR,"Compute reduce fix does not "
                     "calculate a per-surf vector");
        if (argindex[i] && modify->fix[ifix]->size_per_surf_cols == 0)
          error->all(FLERR,"Compute reduce fix does not "
                     "calculate a per-surf array");
        if (argindex[i] &&
            argindex[i] > modify->fix[ifix]->size_per_surf_cols)
          error->all(FLERR,"Compute reduce fix array is accessed out-of-range");
      } else error->all(FLERR,"Compute reduce fix calculates global values");

    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for compute reduce does not exist");
      if (input->variable->particle_style(ivariable)) flavor[i] = PARTICLE;
      else if (input->variable->grid_style(ivariable)) flavor[i] = GRID;
      else
        error->all(FLERR,"Compute reduce variable is not "
                   "particle-style or grid-style variable");
    }

    // require all values have same flavor

    if (i && flavor[i] != flavor[i-1])
      error->all(FLERR,"Compute reduce inputs must be all "
                 "particle, grid, or surf values");
  }

  // require per-surf values for SUMAREA or AVEAREA

  if ((mode == SUMAREA || mode == AVEAREA) && (flavor[0] != SURF))
    error->all(FLERR,"Compure reduce sum-area/ave-area require surf values");

  // this compute produces either a scalar or vector

  if (nvalues == 1) {
    scalar_flag = 1;
    vector = onevec = NULL;
    indices = owner = NULL;
  } else {
    vector_flag = 1;
    size_vector = nvalues;
    vector = new double[size_vector];
    onevec = new double[size_vector];
    indices = new int[size_vector];
    owner = new int[size_vector];
  }

  maxparticle = maxgrid = 0;
  varparticle = vargrid = NULL;
  areasurf = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeReduce::~ComputeReduce()
{
  delete [] which;
  delete [] argindex;
  delete [] flavor;
  for (int m = 0; m < nvalues; m++) delete [] ids[m];
  delete [] ids;
  delete [] value2index;
  delete [] replace;

  delete [] vector;
  delete [] onevec;
  delete [] indices;
  delete [] owner;

  if (subsetID) delete [] subsetID;

  memory->destroy(varparticle);
  memory->destroy(vargrid);
  memory->destroy(areasurf);
}

/* ---------------------------------------------------------------------- */

void ComputeReduce::init()
{
  // if subsetID is defined, identify mixture or grid/surf group

  if (subsetID) {
    if (flavor[0] == PARTICLE) {
      int imix = particle->find_mixture(subsetID);
      if (imix < 0)
        error->all(FLERR,"Compute reduce particle mixture ID does not exist");
      s2g = particle->mixture[imix]->species2group;
    } else if (flavor[0] == GRID) {
      int igroup = grid->find_group(subsetID);
      if (igroup < 0)
        error->all(FLERR,"Compute reduce grid group ID does not exist");
      gridgroupbit = grid->bitmask[igroup];
    } else if (flavor[0] == SURF) {
      int igroup = surf->find_group(subsetID);
      if (igroup < 0)
        error->all(FLERR,"Compute reduce surf group ID does not exist");
      surfgroupbit = surf->bitmask[igroup];
    }
  }

  // flavor = SURF, pre-compute area of all surfs
  // if mode != SUMAREA or AVEAREA, area not computed, just set to 1.0

  if (flavor[0] == SURF) area_total = area_per_surf();

  // set indices of all computes,fixes,variables

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for compute reduce does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for compute reduce does not exist");
      value2index[m] = ifix;

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for compute reduce does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }
}

/* ---------------------------------------------------------------------- */

double ComputeReduce::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double one = compute_one(0,-1);

  if (mode == SUM || mode == SUMSQ || mode == SUMAREA) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  } else if (mode == MINN) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_MIN,world);
  } else if (mode == MAXX) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_MAX,world);
  } else if (mode == AVE || mode == AVESQ) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
    bigint n = count_included();
    if (n) scalar /= n;
  } else if (mode == AVEAREA) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
    if (area_total > 0.0) scalar /= area_total;
  }

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeReduce::compute_vector()
{
  invoked_vector = update->ntimestep;

  for (int m = 0; m < nvalues; m++)
    if (!replace || replace[m] < 0) {
      onevec[m] = compute_one(m,-1);
      indices[m] = index;
    }

  if (mode == SUM || mode == SUMSQ || mode == SUMAREA) {
    for (int m = 0; m < nvalues; m++)
      MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);

  } else if (mode == MINN) {
    if (!replace) {
      for (int m = 0; m < nvalues; m++)
        MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_MIN,world);

    } else {
      for (int m = 0; m < nvalues; m++)
        if (replace[m] < 0) {
          pairme.value = onevec[m];
          pairme.proc = me;
          MPI_Allreduce(&pairme,&pairall,1,MPI_DOUBLE_INT,MPI_MINLOC,world);
          vector[m] = pairall.value;
          owner[m] = pairall.proc;
        }
      for (int m = 0; m < nvalues; m++)
        if (replace[m] >= 0) {
          if (me == owner[replace[m]])
            vector[m] = compute_one(m,indices[replace[m]]);
          MPI_Bcast(&vector[m],1,MPI_DOUBLE,owner[replace[m]],world);
        }
    }

  } else if (mode == MAXX) {
    if (!replace) {
      for (int m = 0; m < nvalues; m++)
        MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_MAX,world);

    } else {
      for (int m = 0; m < nvalues; m++)
        if (replace[m] < 0) {
          pairme.value = onevec[m];
          pairme.proc = me;
          MPI_Allreduce(&pairme,&pairall,1,MPI_DOUBLE_INT,MPI_MAXLOC,world);
          vector[m] = pairall.value;
          owner[m] = pairall.proc;
        }
      for (int m = 0; m < nvalues; m++)
        if (replace[m] >= 0) {
          if (me == owner[replace[m]])
            vector[m] = compute_one(m,indices[replace[m]]);
          MPI_Bcast(&vector[m],1,MPI_DOUBLE,owner[replace[m]],world);
        }
    }

  } else if (mode == AVE || mode == AVESQ) {
    bigint n = count_included();
    for (int m = 0; m < nvalues; m++) {
      MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);
      if (n) vector[m] /= n;
    }

  } else if (mode == AVEAREA) {
    for (int m = 0; m < nvalues; m++) {
      MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);
      if (area_total > 0.0) vector[m] /= area_total;
    }
  }
}

/* ----------------------------------------------------------------------
   calculate reduced value for one input M and return it
   if flag = -1:
     sum/min/max/ave all values in vector
     if mode = MIN or MAX, also set index to which vector value wins
   if flag >= 0: simply return vector[flag]
------------------------------------------------------------------------- */

double ComputeReduce::compute_one(int m, int flag)
{
  int i;

  // invoke the appropriate attribute,compute,fix,variable
  // for flag = -1, compute scalar quantity by scanning over atom properties
  // only include atoms in group for atom properties and per-atom quantities

  index = -1;
  int vidx = value2index[m];
  int aidx = argindex[m];

  double one = 0.0;
  if (mode == MINN) one = BIG;
  else if (mode == MAXX) one = -BIG;

  if (which[m] == X) {
    Particle::OnePart *particles = particle->particles;
    int n = particle->nlocal;
    if (flag < 0) {
      for (i = 0; i < n; i++) {
        if (subsetID && s2g[particles[i].ispecies] < 0) continue;
        combine(one,particles[i].x[aidx],i);
      }
    } else one = particles[flag].x[aidx];

  } else if (which[m] == V) {
    Particle::OnePart *particles = particle->particles;
    int n = particle->nlocal;
    if (flag < 0) {
      for (i = 0; i < n; i++) {
        if (subsetID && s2g[particles[i].ispecies] < 0) continue;
        combine(one,particles[i].v[aidx],i);
      }
    } else one = particles[flag].v[aidx];

  } else if (which[m] == KE) {
    Particle::OnePart *particles = particle->particles;
    Particle::Species *species = particle->species;
    double mvv2e = update->mvv2e;
    int n = particle->nlocal;
    Particle::OnePart *p;
    double *v;
    double ke;

    if (flag < 0) {
      for (i = 0; i < n; i++) {
        if (subsetID && s2g[particles[i].ispecies] < 0) continue;
        p = &particles[i];
        v = p->v;
        ke = mvv2e * 0.5 * species[p->ispecies].mass *
          (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        combine(one,ke,i);
      }
    } else {
      p = &particles[flag];
      v = p->v;
      one = mvv2e * 0.5 * species[p->ispecies].mass *
        (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }

  } else if (which[m] == EROT) {
    Particle::OnePart *particles = particle->particles;
    int n = particle->nlocal;
    if (flag < 0) {
      for (i = 0; i < n; i++) {
        if (subsetID && s2g[particles[i].ispecies] < 0) continue;
        combine(one,particles[i].erot,i);
      }
    } else one = particles[flag].erot;

  } else if (which[m] == EVIB) {
    Particle::OnePart *particles = particle->particles;
    int n = particle->nlocal;
    if (flag < 0) {
      for (i = 0; i < n; i++) {
        if (subsetID && s2g[particles[i].ispecies] < 0) continue;
        combine(one,particles[i].evib,i);
      }
    } else one = particles[flag].evib;

  // invoke compute if not previously invoked
  // for per-grid compute, invoke post_process() method if necessary

  } else if (which[m] == COMPUTE) {
    Compute *c = modify->compute[vidx];

    if (flavor[m] == PARTICLE) {
      if (!(c->invoked_flag & INVOKED_PER_PARTICLE)) {
        c->compute_per_particle();
        c->invoked_flag |= INVOKED_PER_PARTICLE;
      }

      if (aidx == 0) {
        double *cvec = c->vector_particle;
        Particle::OnePart *particles = particle->particles;
        int n = particle->nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++) {
            if (subsetID && s2g[particles[i].ispecies] < 0) continue;
            combine(one,cvec[i],i);
          }
        } else one = cvec[flag];
      } else {
        double **carray = c->array_particle;
        Particle::OnePart *particles = particle->particles;
        int n = particle->nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++) {
            if (subsetID && s2g[particles[i].ispecies] < 0) continue;
            combine(one,carray[i][aidxm1],i);
          }
        } else one = carray[flag][aidxm1];
      }

    } else if (flavor[m] == GRID) {
      if (!(c->invoked_flag & INVOKED_PER_GRID)) {
        c->compute_per_grid();
        c->invoked_flag |= INVOKED_PER_GRID;
      }

      if (c->post_process_grid_flag)
        c->post_process_grid(aidx,1,NULL,NULL,NULL,1);
      else if (c->post_process_isurf_grid_flag)
        c->post_process_isurf_grid();

      if (aidx == 0 || c->post_process_grid_flag) {
        double *cvec = c->vector_grid;
        Grid::ChildInfo *cinfo = grid->cinfo;
        int n = grid->nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++) {
            if (subsetID && !(cinfo[i].mask & gridgroupbit)) continue;
            combine(one,cvec[i],i);
          }
        } else one = cvec[flag];
      } else {
        double **carray = c->array_grid;
        Grid::ChildInfo *cinfo = grid->cinfo;
        int n = grid->nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++) {
            if (subsetID && !(cinfo[i].mask & gridgroupbit)) continue;
            combine(one,carray[i][aidxm1],i);
          }
        } else one = carray[flag][aidxm1];
      }
    }

  // access fix fields, check if fix frequency is a match

  } else if (which[m] == FIX) {
    Fix *fix = modify->fix[vidx];

    if (flavor[m] == PARTICLE) {
      if (update->ntimestep % modify->fix[vidx]->per_particle_freq)
        error->all(FLERR,"Fix used in compute reduce not "
                   "computed at compatible time");
      if (aidx == 0) {
        double *fvec = fix->vector_particle;
        Particle::OnePart *particles = particle->particles;
        int n = particle->nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++) {
            if (subsetID && s2g[particles[i].ispecies] < 0) continue;
            combine(one,fvec[i],i);
          }
        } else one = fvec[flag];
      } else {
        double **farray = fix->array_particle;
        Particle::OnePart *particles = particle->particles;
        int n = particle->nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++) {
            if (subsetID && s2g[particles[i].ispecies] < 0) continue;
            combine(one,farray[i][aidxm1],i);
          }
        } else one = farray[flag][aidxm1];
      }

    } else if (flavor[m] == GRID) {
      if (update->ntimestep % modify->fix[vidx]->per_grid_freq)
        error->all(FLERR,"Fix used in compute reduce not "
                   "computed at compatible time");
      if (aidx == 0) {
        double *fvec = fix->vector_grid;
        Grid::ChildInfo *cinfo = grid->cinfo;
        int n = grid->nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++) {
            if (subsetID && !(cinfo[i].mask & gridgroupbit)) continue;
            combine(one,fvec[i],i);
          }
        } else one = fvec[flag];
      } else {
        double **farray = fix->array_grid;
        Grid::ChildInfo *cinfo = grid->cinfo;
        int n = grid->nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++) {
            if (subsetID && !(cinfo[i].mask & gridgroupbit)) continue;
            combine(one,farray[i][aidxm1],i);
          }
        } else one = farray[flag][aidxm1];
      }

    } else if (flavor[m] == SURF) {
      if (update->ntimestep % modify->fix[vidx]->per_surf_freq)
        error->all(FLERR,"Fix used in compute reduce not "
                   "computed at compatible time");
      if (aidx == 0) {
        double *fvec = fix->vector_surf;

        int dimension = domain->dimension;
        Surf::Line *lines = surf->lines;
        Surf::Line *mylines = surf->mylines;
        Surf::Tri *tris = surf->tris;
        Surf::Tri *mytris = surf->mytris;
        int distributed = surf->distributed;
        int n = surf->nown;

        if (flag < 0) {

          // mask/group logic for selecting surfs matches DumpSurf

          for (i = 0; i < n; i++) {
            if (subsetID) {
              if (dimension == 2) {
                if (!distributed) {
                  if (!(lines[me+i*nprocs].mask & surfgroupbit)) continue;
                } else {
                  if (!(mylines[i].mask & surfgroupbit)) continue;
                }
              } else {
                if (!distributed) {
                  if (!(tris[me+i*nprocs].mask & surfgroupbit)) continue;
                } else {
                  if (!(mytris[i].mask & surfgroupbit)) continue;
                }
              }
            }

            combine(one,areasurf[i]*fvec[i],i);
          }
        } else one = fvec[flag];

      } else {
        double **farray = fix->array_surf;

        int dimension = domain->dimension;
        Surf::Line *lines = surf->lines;
        Surf::Line *mylines = surf->mylines;
        Surf::Tri *tris = surf->tris;
        Surf::Tri *mytris = surf->mytris;
        int distributed = surf->distributed;
        int n = surf->nown;

        int aidxm1 = aidx - 1;
        if (flag < 0) {

          // mask/group logic for selecting surfs matches DumpSurf

          for (i = 0; i < n; i++) {
            if (subsetID) {
              if (dimension == 2) {
                if (!distributed) {
                  if (!(lines[me+i*nprocs].mask & surfgroupbit)) continue;
                } else {
                  if (!(mylines[i].mask & surfgroupbit)) continue;
                }
              } else {
                if (!distributed) {
                  if (!(tris[me+i*nprocs].mask & surfgroupbit)) continue;
                } else {
                  if (!(mytris[i].mask & surfgroupbit)) continue;
                }
              }
            }

            combine(one,areasurf[i]*farray[i][aidxm1],i);
          }
        } else one = farray[flag][aidxm1];
      }
    }

  // evaluate particle-style or grid-style variable

  } else if (which[m] == VARIABLE) {
    if (flavor[m] == PARTICLE) {
      Particle::OnePart *particles = particle->particles;
      int n = particle->nlocal;
      if (n > maxparticle) {
        maxparticle = particle->maxlocal;
        memory->destroy(varparticle);
        memory->create(varparticle,maxparticle,"reduce:varparticle");
      }

      input->variable->compute_particle(vidx,varparticle,1,0);
      if (flag < 0) {
        for (i = 0; i < n; i++) {
          if (subsetID && s2g[particles[i].ispecies] < 0) continue;
          combine(one,varparticle[i],i);
        }
      } else one = varparticle[flag];

    } else if (flavor[m] == GRID) {
      Grid::ChildInfo *cinfo = grid->cinfo;
      int n = grid->nlocal;
      if (n > maxgrid) {
        maxgrid = grid->maxlocal;
        memory->destroy(vargrid);
        memory->create(vargrid,maxgrid,"reduce:vargrid");
      }

      input->variable->compute_grid(vidx,vargrid,1,0);
      if (flag < 0) {
        for (i = 0; i < n; i++) {
          if (subsetID && !(cinfo[i].mask & gridgroupbit)) continue;
          combine(one,vargrid[i],i);
        }
      } else one = vargrid[flag];
    }
  }

  return one;
}

/* ---------------------------------------------------------------------- */

bigint ComputeReduce::count_included()
{
  bigint ncount,ncountall;

  if (flavor[0] == PARTICLE) {
    if (!subsetID) ncount = particle->nlocal;
    else {
      Particle::OnePart *particles = particle->particles;
      int n = particle->nlocal;
      ncount = 0;
      for (int i = 0; i < n; i++)
        if (s2g[particles[i].ispecies] >= 0) ncount++;
    }

  } else if (flavor[0] == GRID) {
    if (!subsetID) ncount = grid->nlocal;
    else {
      Grid::ChildInfo *cinfo = grid->cinfo;
      int n = grid->nlocal;
      ncount = 0;
      for (int i = 0; i < n; i++)
        if (cinfo[i].mask & gridgroupbit) ncount++;
    }

  } else if (flavor[0] == SURF) {
    if (!subsetID) ncount = surf->nown;
    else {
      int dimension = domain->dimension;
      Surf::Line *lines = surf->lines;
      Surf::Line *mylines = surf->mylines;
      Surf::Tri *tris = surf->tris;
      Surf::Tri *mytris = surf->mytris;
      int distributed = surf->distributed;
      int n = surf->nown;

      ncount = 0;
      for (int i = 0; i < n; i++)
        if (dimension == 2) {
          if (!distributed) {
            if (lines[me+i*nprocs].mask & surfgroupbit) ncount++;
          } else {
            if (mylines[i].mask & surfgroupbit) ncount++;
          }
        } else {
          if (!distributed) {
            if (tris[me+i*nprocs].mask & surfgroupbit) ncount++;
          } else {
            if (mytris[i].mask & surfgroupbit) ncount++;
          }
        }
    }
  }

  MPI_Allreduce(&ncount,&ncountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  return ncountall;
}

/* ---------------------------------------------------------------------- */

double ComputeReduce::area_per_surf()
{
  double area_mine,area_all,tmp;

  int dimension = domain->dimension;
  Surf::Line *lines = surf->lines;
  Surf::Line *mylines = surf->mylines;
  Surf::Tri *tris = surf->tris;
  Surf::Tri *mytris = surf->mytris;
  int distributed = surf->distributed;
  int n = surf->nown;

  memory->destroy(areasurf);
  memory->create(areasurf,n,"reduce:areasurf");

  area_mine = 0.0;
  for (int i = 0; i < n; i++) {
    if (dimension == 2) {
      if (!distributed) {
        if (subsetID && !(lines[me+i*nprocs].mask & surfgroupbit)) continue;
        if (mode == SUMAREA || mode == AVEAREA)
          areasurf[i] = surf->line_size(&lines[me+i*nprocs]);
        else areasurf[i] = 1.0;
      } else {
        if (subsetID && !(mylines[i].mask & surfgroupbit)) continue;
        if (mode == SUMAREA || mode == AVEAREA)
          areasurf[i] = surf->line_size(&mylines[i]);
        else areasurf[i] = 1.0;
      }
    } else {
      if (!distributed) {
        if (subsetID && !(tris[me+i*nprocs].mask & surfgroupbit)) continue;
        if (mode == SUMAREA || mode == AVEAREA)
          areasurf[i] = surf->tri_size(&tris[me+i*nprocs],tmp);
        else areasurf[i] = 1.0;
      } else {
        if (subsetID && !(mytris[i].mask & surfgroupbit)) continue;
        if (mode == SUMAREA || mode == AVEAREA)
          areasurf[i] = surf->tri_size(&mytris[i],tmp);
        else areasurf[i] = 1.0;
      }
    }


    area_mine += areasurf[i];
  }

  MPI_Allreduce(&area_mine,&area_all,1,MPI_DOUBLE,MPI_SUM,world);
  return area_all;
}

/* ----------------------------------------------------------------------
   combine two values according to reduction mode
   for MIN/MAX, also update index with winner
------------------------------------------------------------------------- */

void ComputeReduce::combine(double &one, double two, int i)
{
  if (mode == SUM || mode == AVE) one += two;
  else if (mode == SUMSQ || mode == AVESQ) one += two*two;
  else if (mode == MINN) {
    if (two < one) {
      one = two;
      index = i;
    }
  } else if (mode == MAXX) {
    if (two > one) {
      one = two;
      index = i;
    }
  } else if (mode == SUMAREA || mode == AVEAREA) one += two;
}

/* ----------------------------------------------------------------------
   memory usage of varatom
------------------------------------------------------------------------- */

bigint ComputeReduce::memory_usage()
{
  bigint bytes = maxparticle * sizeof(double);
  bytes += maxgrid * sizeof(double);
  return bytes;
}
