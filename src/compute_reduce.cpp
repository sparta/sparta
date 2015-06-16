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

#include "string.h"
#include "stdlib.h"
#include "compute_reduce.h"
#include "update.h"
#include "particle.h"
#include "grid.h"
#include "surf.h"
#include "modify.h"
#include "fix.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{SUM,MINN,MAXX,AVE};
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
  else if (strcmp(arg[2],"min") == 0) mode = MINN;
  else if (strcmp(arg[2],"max") == 0) mode = MAXX;
  else if (strcmp(arg[2],"ave") == 0) mode = AVE;
  else error->all(FLERR,"Illegal compute reduce command");

  MPI_Comm_rank(world,&me);

  // parse remaining values until one isn't recognized

  which = new int[narg-3];
  argindex = new int[narg-3];
  flavor = new int[narg-3];
  ids = new char*[narg-3];
  value2index = new int[narg-3];
  nvalues = 0;

  int iarg = 3;
  while (iarg < narg) {
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

  while (iarg < narg) {
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
      if (input->variable->particle_style(ivariable) == 0)
        error->all(FLERR,
                   "Compute reduce variable is not particle-style variable");
      flavor[i] = PARTICLE;
    }
  }

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

  maxparticle = 0;
  varparticle = NULL;
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

  memory->destroy(varparticle);
}

/* ---------------------------------------------------------------------- */

void ComputeReduce::init()
{
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

  if (mode == SUM) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  } else if (mode == MINN) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_MIN,world);
  } else if (mode == MAXX) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_MAX,world);
  } else if (mode == AVE) {
    MPI_Allreduce(&one,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
    bigint n = count(0);
    if (n) scalar /= n;
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

  if (mode == SUM) {
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

  } else if (mode == AVE) {
    for (int m = 0; m < nvalues; m++) {
      MPI_Allreduce(&onevec[m],&vector[m],1,MPI_DOUBLE,MPI_SUM,world);
      bigint n = count(m);
      if (n) vector[m] /= n;
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

  double one;
  if (mode == SUM) one = 0.0;
  else if (mode == MINN) one = BIG;
  else if (mode == MAXX) one = -BIG;
  else if (mode == AVE) one = 0.0;

  if (which[m] == X) {
    Particle::OnePart *particles = particle->particles;
    int n = particle->nlocal;
    if (flag < 0) {
      for (i = 0; i < n; i++)
        combine(one,particles[i].x[aidx],i);
    } else one = particles[flag].x[aidx];
  } else if (which[m] == V) {
    Particle::OnePart *particles = particle->particles;
    int n = particle->nlocal;
    if (flag < 0) {
      for (i = 0; i < n; i++)
        combine(one,particles[i].v[aidx],i);
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
      for (i = 0; i < n; i++)
        combine(one,particles[i].erot,i);
    } else one = particles[flag].erot;
  } else if (which[m] == EVIB) {
    Particle::OnePart *particles = particle->particles;
    int n = particle->nlocal;
    if (flag < 0) {
      for (i = 0; i < n; i++)
        combine(one,particles[i].evib,i);
    } else one = particles[flag].evib;

  // invoke compute if not previously invoked

  } else if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[vidx];

    if (flavor[m] == PARTICLE) {
      if (!(compute->invoked_flag & INVOKED_PER_PARTICLE)) {
        compute->compute_per_particle();
        compute->invoked_flag |= INVOKED_PER_PARTICLE;
      }

      if (aidx == 0) {
        double *cvec = compute->vector_particle;
        int n = particle->nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            combine(one,cvec[i],i);
        } else one = cvec[flag];
      } else {
        double **carray = compute->array_particle;
        int n = particle->nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            combine(one,carray[i][aidxm1],i);
        } else one = carray[flag][aidxm1];
      }

    } else if (flavor[m] == GRID) {
      if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
        compute->compute_per_grid();
        compute->invoked_flag |= INVOKED_PER_GRID;
      }

      if (aidx == 0) {
        double *cvec = compute->vector_grid;
        int n = grid->nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            combine(one,cvec[i],i);
        } else one = cvec[flag];
      } else {
        double **carray = compute->array_grid;
        int n = grid->nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            combine(one,carray[i][aidxm1],i);
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
        int n = particle->nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            combine(one,fvec[i],i);
        } else one = fvec[flag];
      } else {
        double **farray = fix->array_particle;
        int n = particle->nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            combine(one,farray[i][aidxm1],i);
        } else one = farray[flag][aidxm1];
      }

    } else if (flavor[m] == GRID) {
      if (update->ntimestep % modify->fix[vidx]->per_grid_freq)
        error->all(FLERR,"Fix used in compute reduce not "
                   "computed at compatible time");
      if (aidx == 0) {
        double *fvec = fix->vector_grid;
        int n = grid->nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            combine(one,fvec[i],i);
        } else one = fvec[flag];
      } else {
        double **farray = fix->array_grid;
        int n = grid->nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            combine(one,farray[i][aidxm1],i);
        } else one = farray[flag][aidxm1];
      }

    } else if (flavor[m] == SURF) {
      if (update->ntimestep % modify->fix[vidx]->per_surf_freq)
        error->all(FLERR,"Fix used in compute reduce not "
                   "computed at compatible time");
      if (aidx == 0) {
        double *fvec = fix->vector_surf;
        int n = surf->nlocal;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            combine(one,fvec[i],i);
        } else one = fvec[flag];
      } else {
        double **farray = fix->array_surf;
        int n = surf->nlocal;
        int aidxm1 = aidx - 1;
        if (flag < 0) {
          for (i = 0; i < n; i++)
            combine(one,farray[i][aidxm1],i);
        } else one = farray[flag][aidxm1];
      }
    }

  // evaluate particle-style variable

  } else if (which[m] == VARIABLE) {
    int n = particle->nlocal;
    if (n > maxparticle) {
      maxparticle = particle->maxlocal;
      memory->destroy(varparticle);
      memory->create(varparticle,maxparticle,"reduce:varparticle");
    }

    input->variable->compute_particle(vidx,varparticle,1,0);
    if (flag < 0) {
      for (i = 0; i < n; i++)
        combine(one,varparticle[i],i);
    } else one = varparticle[flag];
  }

  return one;
}

/* ---------------------------------------------------------------------- */

bigint ComputeReduce::count(int m)
{
  bigint ncount,ncountall;

  int vidx = value2index[m];

  if (which[m] == X || which[m] == V) {
    ncount = particle->nlocal;
    MPI_Allreduce(&ncount,&ncountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  } else if (which[m] == KE || which[m] == EROT || which[m] == EVIB) {
    ncount = particle->nlocal;
    MPI_Allreduce(&ncount,&ncountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  } else if (which[m] == COMPUTE) {
    Compute *compute = modify->compute[vidx];
    if (flavor[m] == PARTICLE) {
      ncount = particle->nlocal;
      MPI_Allreduce(&ncount,&ncountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    } else if (flavor[m] == GRID) {
      ncount = grid->nlocal;
      MPI_Allreduce(&ncount,&ncountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    }
  } else if (which[m] == FIX) {
    Fix *fix = modify->fix[vidx];
    if (flavor[m] == PARTICLE) {
      ncount = particle->nlocal;
      MPI_Allreduce(&ncount,&ncountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    } else if (flavor[m] == GRID) {
      ncount = grid->nlocal;
      MPI_Allreduce(&ncount,&ncountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    } else if (flavor[m] == SURF) {
      ncount = surf->nlocal;
      MPI_Allreduce(&ncount,&ncountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    }
  } else if (which[m] == VARIABLE) {
    ncount = particle->nlocal;
    MPI_Allreduce(&ncount,&ncountall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  }

  return ncountall;
}

/* ----------------------------------------------------------------------
   combine two values according to reduction mode
   for MIN/MAX, also update index with winner
------------------------------------------------------------------------- */

void ComputeReduce::combine(double &one, double two, int i)
{
  if (mode == SUM || mode == AVE) one += two;
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
  }
}

/* ----------------------------------------------------------------------
   memory usage of varatom
------------------------------------------------------------------------- */

bigint ComputeReduce::memory_usage()
{
  bigint bytes = maxparticle * sizeof(double);
  return bytes;
}
