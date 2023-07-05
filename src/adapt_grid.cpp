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
#include "adapt_grid.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "comm.h"
#include "grid.h"
#include "surf.h"
#include "collide.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "fix_ave_grid.h"
#include "cut3d.h"
#include "cut2d.h"
#include "output.h"
#include "dump.h"
#include "write_grid.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "hashlittle.h"
#include "my_page.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathExtra;
using namespace MathConst;

enum{NONE,REFINE,COARSEN};              // also in FixAdapt
enum{PARTICLE,SURF,VALUE,COMPUTE,FIX,RANDOM};
enum{REGION_ALL,REGION_ONE,REGION_CENTER};
enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};   // several files
enum{LESS,MORE};
enum{SUM,MINIMUM,MAXIMUM};
enum{ONE,RUNNING};                      // also in FixAveGrid

#define INVOKED_PER_GRID 16
#define DELTA_LIST 1024
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

AdaptGrid::AdaptGrid(SPARTA *sparta) : Pointers(sparta)
{
  me = comm->me;
  nprocs = comm->nprocs;

  valueID = NULL;
  file = NULL;
}

/* ---------------------------------------------------------------------- */

AdaptGrid::~AdaptGrid()
{
  delete [] valueID;
  delete [] file;
}

/* ---------------------------------------------------------------------- */

void AdaptGrid::command(int narg, char **arg)
{
  if (!grid->exist)
    error->all(FLERR,"Cannot adapt grid when grid is not defined");

  if (narg < 1) error->all(FLERR,"Illegal adapt_grid command");

  // process command-line args

  mode = 0;
  process_args(narg,arg);
  check_args(-1);

  // perform adaptation

  if (me == 0) {
    if (screen) fprintf(screen,"Adapting grid ...\n");
    if (logfile) fprintf(logfile,"Adapting grid ...\n");
  }

  MPI_Barrier(world);
  double time1 = MPI_Wtime();

  // invoke init()
  // so all grid cell info including collide & fixes is ready to migrate

  sparta->init();
  grid->remove_ghosts();

  // iterate over refinement and coarsening actions

  bigint nrefine_total = 0;
  bigint ncoarsen_total = 0;
  bigint nrefine = 0;
  bigint ncoarsen = 0;

  for (int iter = 0; iter < niterate; iter++) {
    setup(iter);

    if (action1 == REFINE || action2 == REFINE) grid->maxlevel++;

    if (action1 == REFINE) nrefine = refine();
    else if (action1 == COARSEN) ncoarsen = coarsen();

    if (action2 == REFINE) nrefine = refine();
    else if (action2 == COARSEN) ncoarsen = coarsen();

    grid->set_maxlevel();
    grid->rehash();

    nrefine_total += nrefine;
    ncoarsen_total += ncoarsen;

    cleanup();

    if (nrefine == 0 && ncoarsen == 0) break;
  }

  // if no refine or coarsen, just reghost/reneighbor and return

  if (nrefine_total == 0 && ncoarsen_total == 0) {
    grid->acquire_ghosts();
    grid->find_neighbors();

    if (me == 0) {
      if (screen) fprintf(screen,"  no grid adaptation performed\n");
      if (logfile) fprintf(logfile,"  no grid adaptation performed\n");
    }
    return;
  }

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // reset all attributes of adapted grid
  // same steps as in create_grid

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->find_neighbors();
  grid->check_uniform();
  comm->reset_neighbors();

  if (surf->exist) {
    grid->set_inout();
    grid->type_check();
  }

  // if not before first run,
  // notify all classes that store per-grid data that grid may have changed

  if (update->first_update) grid->notify_changed();

  // if explicit distributed surfs
  // set redistribute timestep and clear custom status flags

  if (surf->distributed && !surf->implicit) {
    surf->localghost_changed_step = update->ntimestep;
    for (int i = 0; i < surf->ncustom; i++)
      surf->estatus[i] = 0;
  }

  // write out new grid file

  if (file) write_file();

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // stats

  double time_total = time3-time1;

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  " BIGINT_FORMAT " cells refined, " BIGINT_FORMAT
              " cells coarsened\n",nrefine_total,ncoarsen_total);
      fprintf(screen,"  adapted to " BIGINT_FORMAT " grid cells\n",grid->ncell);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  adapt/redo percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
    }

    if (logfile) {
      fprintf(logfile,"  " BIGINT_FORMAT " cells refined, " BIGINT_FORMAT
              " cells coarsened\n",nrefine_total,ncoarsen_total);
      fprintf(logfile,"  adapted to " BIGINT_FORMAT " grid cells\n",grid->ncell);
      fprintf(logfile,"  CPU time = %g secs\n",time_total);
      fprintf(logfile,"  adapt/redo percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
    }
  }
}

/* ----------------------------------------------------------------------
   process command args for both adapt_grid and fix adapt
------------------------------------------------------------------------- */

void AdaptGrid::process_args(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal adapt command");

  int igroup = grid->find_group(arg[0]);
  if (igroup < 0) error->all(FLERR,"Adapt_grid group ID does not exist");
  groupbit = grid->bitmask[igroup];

  // define action(s)

  if (strcmp(arg[1],"refine") == 0) action1 = REFINE;
  else if (strcmp(arg[1],"coarsen") == 0) action1 = COARSEN;
  else error->all(FLERR,"Illegal adapt command");
  int iarg = 2;

  if (strcmp(arg[2],"refine") == 0) {
    action2 = REFINE;
    iarg = 3;
  } else if (strcmp(arg[2],"coarsen") == 0) {
    action2 = COARSEN;
    iarg = 3;
  } else action2 = NONE;

  if (action1 == action2) error->all(FLERR,"Illegal adapt command");

  // define style

  if (strcmp(arg[iarg],"particle") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal adapt command");
      style = PARTICLE;
      rcount = input->numeric(FLERR,arg[iarg+1]);
      ccount = input->numeric(FLERR,arg[iarg+2]);
      iarg += 3;

  } else if (strcmp(arg[iarg],"surf") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal adapt command");
      style = SURF;
      int igroup = surf->find_group(arg[iarg+1]);
      if (igroup < 0)
        error->all(FLERR,"Adapt command surface group does not exist");
      sgroupbit = surf->bitmask[igroup];
      surfsize = input->numeric(FLERR,arg[iarg+2]);
      iarg += 3;

  } else if (strcmp(arg[iarg],"value") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal adapt command");
      style = VALUE;
      if (strncmp(arg[iarg+1],"c_",2) == 0) valuewhich = COMPUTE;
      else if (strncmp(arg[iarg+1],"f_",2) == 0) valuewhich = FIX;
      else error->all(FLERR,"Illegal adapt command");

      int n = strlen(arg[iarg+1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg+1][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal adapt command");
        valindex = atoi(ptr+1);
        *ptr = '\0';
      } else valindex = 0;
      n = strlen(suffix) + 1;
      valueID = new char[n];
      strcpy(valueID,suffix);
      delete [] suffix;

      rvalue = input->numeric(FLERR,arg[iarg+2]);
      cvalue = input->numeric(FLERR,arg[iarg+3]);
      iarg += 4;

  } else if (strcmp(arg[iarg],"random") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal adapt command");
      style = RANDOM;
      rfrac = input->numeric(FLERR,arg[iarg+1]);
      cfrac = input->numeric(FLERR,arg[iarg+2]);
      if (rfrac < 0.0 || rfrac > 1.0 || cfrac < 0.0 || cfrac > 1.0)
        error->all(FLERR,"Illegal adapt command");
      iarg += 3;

  } else error->all(FLERR,"Illegal adapt command");

  // optional args

  niterate = 1;
  minlevel = 1;
  maxlevel = 0;
  rdecide = MORE;
  cdecide = LESS;
  combine = SUM;
  nx = ny = nz = 2;
  if (domain->dimension == 2) nz = 1;
  region = NULL;
  sdir[0] = sdir[1] = sdir[2] = 0.0;
  file = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"iterate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal adapt command");
      niterate = input->inumeric(FLERR,arg[iarg+1]);
      if (mode) error->all(FLERR,"Illegal adapt command");
      if (niterate <= 0) error->all(FLERR,"Illegal adapt command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"maxlevel") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal adapt command");
      maxlevel = input->inumeric(FLERR,arg[iarg+1]);
      if (maxlevel < 0) error->all(FLERR,"Illegal adapt command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"minlevel") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal adapt command");
      minlevel = input->inumeric(FLERR,arg[iarg+1]);
      if (minlevel < 1) error->all(FLERR,"Illegal adapt command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"thresh") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal adapt command");
      if (strcmp(arg[iarg+1],"less") == 0) rdecide = LESS;
      else if (strcmp(arg[iarg+1],"more") == 0) rdecide = MORE;
      else error->all(FLERR,"Illegal adapt command");
      if (strcmp(arg[iarg+2],"less") == 0) cdecide = LESS;
      else if (strcmp(arg[iarg+2],"more") == 0) cdecide = MORE;
      else error->all(FLERR,"Illegal adapt command");
      iarg += 3;
    } else if (strcmp(arg[iarg],"combine") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal adapt command");
      if (strcmp(arg[iarg+1],"sum") == 0) combine = SUM;
      else if (strcmp(arg[iarg+1],"min") == 0) combine = MINIMUM;
      else if (strcmp(arg[iarg+1],"max") == 0) combine = MAXIMUM;
      else error->all(FLERR,"Illegal adapt command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"cells") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal adapt command");
      nx = input->inumeric(FLERR,arg[iarg+1]);
      ny = input->inumeric(FLERR,arg[iarg+2]);
      nz = input->inumeric(FLERR,arg[iarg+3]);
      if (nx < 1 || ny < 1 || nz < 1)
        error->all(FLERR,"Illegal adapt command");
      if (domain->dimension == 2 && nz != 1)
        error->all(FLERR,"Adapt cells nz must be 1 for a 2d simulation");
      if (nx < 1 || ny < 1 || nz < 1)
        error->all(FLERR,"Adapt cells nx,ny,nz cannot be < 1");
      if (nx == 1 && ny == 1 && nz == 1)
        error->all(FLERR,"Adapt cells nx,ny,nz cannot all be one");
      iarg += 4;

    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal adapt command");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Adapt region ID does not exist");
      region = domain->regions[iregion];
      if (strcmp(arg[iarg+2],"all") == 0) regstyle = REGION_ALL;
      else if (strcmp(arg[iarg+2],"one") == 0) regstyle = REGION_ONE;
      else if (strcmp(arg[iarg+2],"center") == 0) regstyle = REGION_CENTER;
      else error->all(FLERR,"Illegal group command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"dir") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal adapt command");
      sdir[0] = input->numeric(FLERR,arg[iarg+1]);
      sdir[1] = input->numeric(FLERR,arg[iarg+2]);
      sdir[2] = input->numeric(FLERR,arg[iarg+3]);
      if (domain->dimension == 2 && sdir[2] != 0.0)
        error->all(FLERR,"Illegal adapt command");
      iarg += 4;

    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal adapt command");
      int n = strlen(arg[iarg+1]) + 1;
      file = new char[n];
      strcpy(file,arg[iarg+1]);
      iarg += 2;

    } else error->all(FLERR,"Illegal adapt command");
  }

  // use maxlevel and nx,ny,nz to set params for Grid::plevels
  // if necessary, limit maxlevel to what is allowed by cellint size

  Grid::ParentLevel *plevels = grid->plevels;

  int maxlevel_request = maxlevel;
  if (maxlevel == 0) maxlevel = grid->plevel_limit;
  int level = MIN(maxlevel,grid->maxlevel);
  int newbits = grid->id_bits(nx,ny,nz);

  while (level < maxlevel) {
    if (plevels[level-1].nbits + plevels[level-1].newbits +
        newbits > 8*sizeof(cellint)) {
      maxlevel = level;
      break;
    }

    plevels[level].nbits = plevels[level-1].nbits + plevels[level-1].newbits;
    plevels[level].newbits = newbits;
    plevels[level].nx = nx;
    plevels[level].ny = ny;
    plevels[level].nz = nz;
    plevels[level].nxyz = (bigint) nx * ny * nz;
    level++;
  }

  if (maxlevel_request && maxlevel < maxlevel_request && me == 0) {
    char str[128];
    sprintf(str,"Reduced maxlevel because it induces "
            "cell IDs that exceed %d bits",(int) sizeof(cellint)*8);
    error->warning(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   error check on value compute/fix for both adapt_grid and fix adapt
   nevery is passed only from fix adapt (mode = 1)
------------------------------------------------------------------------- */

void AdaptGrid::check_args(int nevery)
{
  // if fix ave/grid used require that:
  //   (1) fix adapt Nevery is multiple of fix ave Nfreq
  //   (2) fix ave/grid is defined before fix adapt (checked in fix adapt)
  // if any fix ave/grid defined require that:
  //   (3) fix ave/grid is not a running ave
  // this insures fix ave/grid values will be up-to-date before
  //   adaptation changes the grid
  // NOTE: at some point could make fix ave/grid do full interpolation
  //       for new grid cells some of these restrictions are not needed

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"ave/grid") == 0) {
      if (((FixAveGrid *) modify->fix[i])->ave == RUNNING)
        error->all(FLERR,"Adapt command does not yet allow use of "
                   "fix ave/grid ave running");
    }
  }

  if (style != VALUE) return;

  if (valuewhich == COMPUTE) {
    icompute = modify->find_compute(valueID);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for adapt does not exist");
    if (modify->compute[icompute]->per_grid_flag == 0)
      error->all(FLERR,
                 "Adapt compute does not calculate per-grid values");
    if (valindex == 0 && modify->compute[icompute]->size_per_grid_cols != 0)
      error->all(FLERR,"Adapt compute does not calculate a per-grid vector");
    if (valindex && modify->compute[icompute]->size_per_grid_cols == 0)
      error->all(FLERR,"Adapt compute does not calculate a per-grid array");
    if (valindex && valindex > modify->compute[icompute]->size_per_grid_cols)
      error->all(FLERR,"Adapt compute array is accessed out-of-range");

  } else if (valuewhich == FIX) {
    ifix = modify->find_fix(valueID);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for adapt does not exist");
    if (modify->fix[ifix]->per_grid_flag == 0)
      error->all(FLERR,"Adapt fix does not calculate per-grid values");
    if (valindex == 0 && modify->fix[ifix]->size_per_grid_cols != 0)
      error->all(FLERR,"Adapt fix does not calculate a per-grid vector");
    if (valindex && modify->fix[ifix]->size_per_grid_cols == 0)
      error->all(FLERR,"Adapt fix does not calculate a per-grid array");
    if (valindex && valindex > modify->fix[ifix]->size_per_grid_cols)
      error->all(FLERR,"Adapt fix array is accessed out-of-range");

    // do not allow adapt_grid to use niterate > 1 with fix ave/grid
    // this is b/c the fix's values will not be valid on 2nd iteration
    // OK on 1st iteration, b/c altered cells are not eligble for action2
    // NOTE: at some point could make fix ave/grid do full interpolation
    //       for new grid cells so this restriction is not needed

    if (mode == 0) {
      if (update->ntimestep % modify->fix[ifix]->per_grid_freq)
        error->all(FLERR,"Fix for adapt not computed at compatible time");
      if (niterate > 1)
        error->all(FLERR,"Adapt_grid can not perform multiple iterations "
                   "if using fix as a value");
    } else {
      if (nevery % modify->fix[ifix]->per_grid_freq)
        error->all(FLERR,"Fix for adapt not computed at compatible time");
    }
  }
}

/* ----------------------------------------------------------------------
   setup for both adapt_grid and fix adapt
------------------------------------------------------------------------- */

void AdaptGrid::setup(int iter)
{
  // create RNG for style = RANDOM

  if (style == RANDOM) {
    random = new RanKnuth(update->ranmaster->uniform());
    double seed = update->ranmaster->uniform();
    random->reset(seed,comm->me,100);
  } else random = NULL;

  // list of new cell indices for one refined cell
  // refinement may occur between minlevel and maxlevel-1 inclusive

  Grid::ParentLevel *plevels = grid->plevels;

  int nmax = 0;
  for (int i = minlevel; i < maxlevel; i++)
    nmax = MAX(nmax,plevels[i].nx * plevels[i].ny * plevels[i].nz);
  childlist = new int[nmax];

  // rlist and clist for refine/coarsen

  rlist = NULL;
  clist = NULL;
  alist = NULL;
  cnummax = anummax = 0;

  // for cut/split of new cells by surfaces

  cut3d = NULL;
  cut2d = NULL;
  if (surf->exist) {
    if (domain->dimension == 3) cut3d = new Cut3d(sparta);
    else cut2d = new Cut2d(sparta,domain->axisymmetric);
  }

  // allocate rhash for new parent cell IDs created by refining
  // allocate chash for new child cell IDs created by coarsening

  rhash = new MyHash();
  chash = new MyHash();
}

/* ----------------------------------------------------------------------
   perform refinement
------------------------------------------------------------------------- */

bigint AdaptGrid::refine()
{
  grid->rehash();
  particle->sort();

  candidates_refine();

  if (style == PARTICLE) refine_particle();
  else if (style == SURF) refine_surf();
  else if (style == VALUE) refine_value();
  else if (style == RANDOM) refine_random();

  bigint nme = perform_refine();

  // nrefine = # of refined cells across all processors

  bigint nrefine;
  MPI_Allreduce(&nme,&nrefine,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // if any refinement
  // compress() removes child cells that became parents

  if (nrefine) {
    if (collide) collide->adapt_grid();
    grid->compress();

    // reallocate per grid cell arrays in per grid computes
    // also unset their invoked_flag so that if needed on this timestep
    //   by any other caller, it will be invoked again with changed grid

    Compute **compute = modify->compute;
    for (int i = 0; i < modify->ncompute; i++)
      if (compute[i]->per_grid_flag) {
        compute[i]->reallocate();
        compute[i]->invoked_flag = 0;
      }
  }

  return nrefine;
}

/* ----------------------------------------------------------------------
   perform coarsening
------------------------------------------------------------------------- */

bigint AdaptGrid::coarsen()
{
  grid->rehash();
  particle->sort();

  candidates_coarsen();

  if (style == PARTICLE) coarsen_particle();
  else if (style == SURF) coarsen_surf();
  else if (style == VALUE) coarsen_value();
  else if (style == RANDOM) coarsen_random();

  particle_surf_comm();
  bigint nme = perform_coarsen();

  // ncoarsen = # of coarsened cells across all processors

  bigint ncoarsen;
  MPI_Allreduce(&nme,&ncoarsen,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // if any coarsening
  // compress() removes child cells that vanished
  // surf->compress() removes surfs no longer referenced by owned cells
  //   can be needed when owned cells are removed by coarsening
  //   no need to call in refine() since new child cells own same surfs

  if (ncoarsen) {
    if (collide) collide->adapt_grid();
    grid->compress();
    if (particle->exist) particle->compress_rebalance();
    if (surf->distributed) surf->compress_explicit();

    // reallocate per grid cell arrays in per grid computes
    // also unset their invoked_flag so that if needed on this timestep
    //   by any other caller, it will be invoked again with changed grid

    Compute **compute = modify->compute;
    for (int i = 0; i < modify->ncompute; i++)
      if (compute[i]->per_grid_flag) {
        compute[i]->reallocate();
        compute[i]->invoked_flag = 0;
      }
  }

  return ncoarsen;
}

/* ----------------------------------------------------------------------
   create list of child cells I will possibly refine
   check grid group, split status, maxlevel, INSIDE status, region
   check with chash skips a child cell created by coarsening
------------------------------------------------------------------------- */

void AdaptGrid::candidates_refine()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  memory->create(rlist,nglocal,"adapt_grid:rlist");
  rnum = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;
    if (cells[icell].level >= maxlevel) continue;
    if (cinfo[icell].type == INSIDE) continue;
    if (chash->find(cells[icell].id) != chash->end()) continue;
    if (region && !region_check(cells[icell].lo,cells[icell].hi)) continue;
    rlist[rnum++] = icell;
  }
}

/* ----------------------------------------------------------------------
   refine based on particle count
------------------------------------------------------------------------- */

void AdaptGrid::refine_particle()
{
  int icell,np,nsplit,jcell;
  int *csubs;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;

  int n = 0;
  for (int i = 0; i < rnum; i++) {
    icell = rlist[i];

    if (cells[icell].nsplit == 1) np = cinfo[icell].count;
    else {
      np = 0;
      nsplit = cells[icell].nsplit;
      csubs = sinfo[cells[icell].isplit].csubs;
      for (int j = 0; j < nsplit; j++) {
        jcell = csubs[j];
        np += cinfo[jcell].count;
      }
    }

    if (np > rcount) rlist[n++] = icell;
  }

  rnum = n;
}

/* ----------------------------------------------------------------------
   refine if cell contains eligible surfs
   eligible = part of sgroup and norm opposed to sdir (or sdir = 0,0,0)
   do not create cells smaller in any dimension than surfsize
------------------------------------------------------------------------- */

void AdaptGrid::refine_surf()
{
  int j,m,icell,flag,nsurf,plevel;
  surfint *csurfs;
  double *norm,*lo,*hi;

  int dim = domain->dimension;
  Grid::ChildCell *cells = grid->cells;
  Grid::ParentLevel *plevels = grid->plevels;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  int n = 0;
  for (int i = 0; i < rnum; i++) {
    icell = rlist[i];
    if (!cells[icell].nsurf) continue;
    nsurf = cells[icell].nsurf;
    csurfs = cells[icell].csurfs;
    for (j = 0; j < nsurf; j++) {
      m = csurfs[j];
      if (dim == 2) {
        if (!(lines[m].mask & sgroupbit)) continue;
      } else {
        if (!(tris[m].mask & sgroupbit)) continue;
      }
      if (dim == 2) norm = lines[m].norm;
      else norm = tris[m].norm;
      if (MathExtra::dot3(norm,sdir) <= 0.0) break;
    }
    if (j == nsurf) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    plevel = cells[icell].level;
    flag = 1;
    if (fabs(hi[0]-lo[0])/plevels[plevel].nx < surfsize) flag = 0;
    if (fabs(hi[1]-lo[1])/plevels[plevel].ny < surfsize) flag = 0;
    if (dim == 3 && fabs(hi[2]-lo[2])/plevels[plevel].nz < surfsize) flag = 0;
    if (flag) rlist[n++] = icell;
  }
  rnum = n;
}

/* ----------------------------------------------------------------------
   refine based on compute or fix value
------------------------------------------------------------------------- */

void AdaptGrid::refine_value()
{
  int icell,nsplit,jcell;
  double value;
  int *csubs;

  // invoke compute each time refinement is done
  // grid could have changed from previous refinement or coarsening

  if (valuewhich == COMPUTE) {
    compute = modify->compute[icompute];
    compute->compute_per_grid();
    if (compute->post_process_grid_flag)
      compute->post_process_grid(valindex,1,NULL,NULL,NULL,1);
  } else if (valuewhich == FIX) fix = modify->fix[ifix];

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  int n = 0;
  for (int i = 0; i < rnum; i++) {
    icell = rlist[i];

    if (cells[icell].nsplit == 1) {
      if (valuewhich == COMPUTE) value = value_compute(icell);
      else if (valuewhich == FIX) value = value_fix(icell);
    } else {
      nsplit = cells[icell].nsplit;
      csubs = sinfo[cells[icell].isplit].csubs;
      if (combine == SUM) value = 0.0;
      else if (combine == MINIMUM) value = BIG;
      else if (combine == MAXIMUM) value = -BIG;
      for (int j = 0; j < nsplit; j++) {
        jcell = csubs[j];
        if (combine == SUM) {
          if (valuewhich == COMPUTE) value += value_compute(jcell);
          else if (valuewhich == FIX) value += value_fix(jcell);
        } else if (combine == MINIMUM) {
          if (valuewhich == COMPUTE)
            value = MIN(value,value_compute(jcell));
          else if (valuewhich == FIX)
            value = MIN(value,value_fix(jcell));
        } else if (combine == MAXIMUM) {
          if (valuewhich == COMPUTE)
            value = MAX(value,value_compute(jcell));
          else if (valuewhich == FIX)
            value = MAX(value,value_fix(jcell));
        }
      }
    }

    if (rdecide == LESS) {
      if (value < rvalue) rlist[n++] = icell;
    } else {
      if (value > rvalue) rlist[n++] = icell;
    }
  }

  rnum = n;
}

/* ----------------------------------------------------------------------
   refine based on random #
   useful for debugging
------------------------------------------------------------------------- */

void AdaptGrid::refine_random()
{
  int icell;

  int n = 0;
  for (int i = 0; i < rnum; i++) {
    icell = rlist[i];
    if (random->uniform() >= rfrac) continue;
    rlist[n++] = icell;
  }

  rnum = n;
}

/* ----------------------------------------------------------------------
   perform refinement of child cells in rlist
------------------------------------------------------------------------- */

int AdaptGrid::perform_refine()
{
  int icell,jcell;
  int nglocal,nglocalprev;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int ilist = 0; ilist < rnum; ilist++) {
    icell = rlist[ilist];

    // add cell ID being refined to rhash

    (*rhash)[cells[icell].id] = 0;

    // flag icell for deletion
    // don't need to flag sub cells, since grid compress will do that

    cells[icell].proc = -1;

    // add Nx by Ny by Nz child cells to replace icell
    // refine_cell() adds all new unsplit/split cells to grid hash

    nglocalprev = grid->nlocal;
    grid->refine_cell(icell,childlist,cut2d,cut3d);
    cells = grid->cells;
    cinfo = grid->cinfo;
    nglocal = grid->nlocal;

    // set each new child cell group mask to same value as parent
    // for unsplit and split and sub cells

    for (jcell = nglocalprev; jcell < nglocal; jcell++)
      cinfo[jcell].mask = cinfo[icell].mask;
  }

  return rnum;
}

/* ----------------------------------------------------------------------
   build clist = list of candidate parent cells this proc may coarsen
   consists of 3 kinds of parent cells
   (a) ones where this proc owns all its child cells and they exist
   (b) ones where multiple procs own all its child cells and they exist
       the parent is assigned to this proc via random hash()
   (c) ones where all its child cells do not exist
       parent is assigned same as in (b)
   for a,b: clist->flag = 1 (possible coarsening)
   for c: clist->flag = 0 (no coarsening will be done)
   final decision on which to coarsen depends on child cell attributes
------------------------------------------------------------------------- */

void AdaptGrid::candidates_coarsen()
{
  int m,n,proc,level,nxyz,nchild;
  cellint parentID;
  double lo[3],hi[3];

  // for style = VALUE, invoke compute each time coarsening is done
  // grid could have changed from previous refinement or coarsening

  if (style == VALUE) {
    if (valuewhich == COMPUTE) {
      compute = modify->compute[icompute];
      compute->compute_per_grid();
      if (compute->post_process_grid_flag)
        compute->post_process_grid(valindex,1,NULL,NULL,NULL,1);
    } else if (valuewhich == FIX) fix = modify->fix[ifix];
  }

  // scan my child cells to identify possible parents cells to coarsen
  // exclude sub-cells and parents which do not meet level/mask/region criteria
  // phash = unique parent cells of my children
  //         key = parentID, value = # of child cells I own
  // active = 1 for participating child cells

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  int *active;
  memory->create(active,nglocal,"adapt_grid:active");
  MyHash *phash = new MyHash();

  for (int icell = 0; icell < nglocal; icell++) {
    active[icell] = 0;

    // exclude child cells from consideration
    // check with rhash skips a child cell created by refining

    if (cells[icell].nsplit < 1) continue;
    if (cells[icell].level <= minlevel) continue;
    if (!(cinfo[icell].mask & groupbit)) continue;
    parentID = grid->id_coarsen(cells[icell].id,cells[icell].level);
    if (rhash->find(parentID) != rhash->end()) continue;
    if (region) {
      grid->id_lohi(parentID,cells[icell].level-1,boxlo,boxhi,lo,hi);
      if (!region_check(lo,hi)) continue;
    }

    active[icell] = 1;
    if (phash->find(parentID) != phash->end()) (*phash)[parentID]++;
    else (*phash)[parentID] = 1;
  }

  // create 1st portion of clist = data struct for parentIDs assigned to me
  //   1st portion = I own all children of the parent, so I own the parent
  // parent cells where children are owned by multiple procs
  //   these are assigned randomly to rendezvous procs via hashlittle()
  //   the ones assined to this proc will become 2nd portion of clist
  // inbuf = datums to send to owners of 2nd portion parents, one datum per child cell
  // clhash = list of parentIDs assigned to me, only for 1st portion at this point
  //          key = parentID, value = index into clist

  Grid::ParentLevel *plevels = grid->plevels;
  MyHash *clhash = new MyHash();

  int *proclist;
  memory->create(proclist,nglocal,"adapt_grid:proclist");
  Rvous1 *inbuf = (Rvous1 *) memory->smalloc((bigint) nglocal*sizeof(Rvous1),
                                             "adapt_grid:inbuf");
  cnum = 0;
  int nsend = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!active[icell]) continue;
    level = cells[icell].level;
    parentID = grid->id_coarsen(cells[icell].id,level);
    nxyz = plevels[level-1].nxyz;
    nchild = (*phash)[parentID];

    // this proc owns all children of parentID
    // create a clist entry, and add parentID to clhash

    if (nchild == nxyz) {
      if (clhash->find(parentID) == clhash->end()) {
        if (cnum == cnummax) {
          cnummax += DELTA_LIST;
          clist = (CList *) memory->srealloc(clist,cnummax*sizeof(CList),
                                             "adapt_grid:clist");
        }
        (*clhash)[parentID] = cnum;
        m = cnum++;
        clist[m].parentID = parentID;
        clist[m].plevel = level-1;
        clist[m].nchild = nxyz;
        clist[m].nexist = 0;
        clist[m].proc = new int[nxyz];
        clist[m].index = new int[nxyz];
        clist[m].value = new double[nxyz];
      } else m = (*clhash)[parentID];

      n = clist[m].nexist;
      clist[m].proc[n] = me;
      clist[m].index[n] = icell;
      if (style == PARTICLE) clist[m].value[n] = coarsen_particle_cell(icell);
      else if (style == SURF) clist[m].value[n] = coarsen_surf_cell(icell);
      else if (style == VALUE) clist[m].value[n] = coarsen_value_cell(icell);
      else if (style == RANDOM) clist[m].value[n] = 0.0;
      clist[m].nexist++;

    // this proc does not own all children of parentID
    // add cell info to inbuf to perform rendezvous comm with

    } else {
      proclist[nsend] = hashlittle(&parentID,sizeof(cellint),0) % nprocs;
      inbuf[nsend].parentID = parentID;
      inbuf[nsend].plevel = level-1;
      inbuf[nsend].proc = me;
      inbuf[nsend].icell = icell;
      inbuf[nsend].ichild = grid->id_ichild(cells[icell].id,parentID,level-1) - 1;
      if (style == PARTICLE) inbuf[nsend].value = coarsen_particle_cell(icell);
      else if (style == SURF) inbuf[nsend].value = coarsen_surf_cell(icell);
      else if (style == VALUE) inbuf[nsend].value = coarsen_value_cell(icell);
      else if (style == RANDOM) inbuf[nsend].value = 0.0;
      nsend++;
    }
  }

  memory->destroy(active);

  // perform rendezvous communication to acquire inbuf data
  // each Rvous proc receives child cell data from SPARTA decomp for
  //   all children of parent cells assigned to it
  // callback() method is NULL, just receive data via MPI_All2allv()

  char *buf;
  int nreturn = comm->rendezvous(1,nsend,(char *) inbuf,sizeof(Rvous1),
                                 0,proclist,NULL,0,buf,0,this,0);

  Rvous1 *outbuf = (Rvous1 *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // scan received outbuf to add child cell data to 2nd portion of clist
  // if new parentIDs to clhash

  for (int i = 0; i < nreturn; i++) {
    parentID = outbuf[i].parentID;
    if (clhash->find(parentID) == clhash->end()) {
      (*clhash)[parentID] = cnum;
      if (cnum == cnummax) {
        cnummax += DELTA_LIST;
        clist = (CList *) memory->srealloc(clist,cnummax*sizeof(CList),
                                           "adapt_grid:clist");
      }
      clist[cnum].parentID = parentID;
      clist[cnum].plevel = outbuf[i].plevel;
      nchild = plevels[outbuf[i].plevel].nxyz;
      clist[cnum].nchild = nchild;
      clist[cnum].nexist = 0;
      clist[cnum].proc = new int[nchild];
      clist[cnum].index = new int[nchild];
      clist[cnum].value = new double[nchild];
      m = cnum++;
    } else m = (*clhash)[parentID];

    n = outbuf[i].ichild;
    clist[m].proc[n] = outbuf[i].proc;
    clist[m].index[n] = outbuf[i].icell;
    clist[m].value[n] = outbuf[i].value;
    clist[m].nexist++;
  }

  // can now determine whether all the child cells of 2nd portion parent cells exist
  // flag = 1 if all children of a parentID exist
  // flag = 0 if they do not

  for (int i = 0; i < cnum; i++) {
    if (clist[i].nexist == clist[i].nchild) clist[i].flag = 1;
    else clist[i].flag = 0;
  }

  // clean up

  memory->sfree(outbuf);
  delete phash;
  delete clhash;
}

/* ----------------------------------------------------------------------
   return particle count for one child cell
------------------------------------------------------------------------- */

double AdaptGrid::coarsen_particle_cell(int icell)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;

  int np = 0;

  if (cells[icell].nsplit == 1) np = cinfo[icell].count;
  else {
    int jcell;
    int nsplit = cells[icell].nsplit;
    int *csubs = sinfo[cells[icell].isplit].csubs;
    for (int j = 0; j < nsplit; j++) {
      jcell = csubs[j];
      np += cinfo[jcell].count;
    }
  }

  return (double) np;
}

/* ----------------------------------------------------------------------
   return 1 if any surf in child cell meets coarsening criterion
   else return 0
------------------------------------------------------------------------- */

double AdaptGrid::coarsen_surf_cell(int icell)
{
  double *norm;

  int dim = domain->dimension;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  Grid::ChildCell *cells = grid->cells;
  int nsurf = cells[icell].nsurf;
  surfint *csurfs = cells[icell].csurfs;

  int anysurf = 0;

  for (int i = 0; i < nsurf; i++) {
    if (dim == 2) {
      if (!(lines[i].mask & sgroupbit)) continue;
    } else {
      if (!(tris[i].mask & sgroupbit)) continue;
    }
    if (dim == 2) norm = lines[csurfs[i]].norm;
    else norm = tris[csurfs[i]].norm;
    if (MathExtra::dot3(norm,sdir) < 0.0) {
      anysurf = 1;
      break;
    }
  }

  return (double) anysurf;
}

/* ----------------------------------------------------------------------
   return compute or fix value for child cell
   if child cell is split cell, accumulate value over sub cells
------------------------------------------------------------------------- */

double AdaptGrid::coarsen_value_cell(int icell)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  if (cells[icell].nsplit == 1) {
    if (valuewhich == COMPUTE) return value_compute(icell);
    else if (valuewhich == FIX) return value_fix(icell);
  }

  int *csubs = sinfo[cells[icell].isplit].csubs;
  int nsplit = cells[icell].nsplit;

  int jcell;
  double value = 0.0;

  for (int i = 0; i < nsplit; i++) {
    jcell = csubs[i];
    if (combine == SUM) {
      if (valuewhich == COMPUTE) value += value_compute(jcell);
      else if (valuewhich == FIX) value += value_fix(jcell);
    } else if (combine == MINIMUM) {
      if (valuewhich == COMPUTE)
        value = MIN(value,value_compute(jcell));
      else if (valuewhich == FIX)
        value = MIN(value,value_fix(jcell));
    } else if (combine == MAXIMUM) {
      if (valuewhich == COMPUTE)
        value = MAX(value,value_compute(jcell));
      else if (valuewhich == FIX)
        value = MAX(value,value_fix(jcell));
    }
  }

  return value;
}

/* ----------------------------------------------------------------------
   coarsen clist cells based on particle count
------------------------------------------------------------------------- */

void AdaptGrid::coarsen_particle()
{
  int m,nchild,np;
  double *values;

  for (int i = 0; i < cnum; i++) {
    if (clist[i].flag == 0) continue;

    values = clist[i].value;
    nchild = clist[i].nchild;

    np = 0;
    for (m = 0; m < nchild; m++) np += static_cast<int> (values[m]);
    if (np < ccount) clist[i].flag = 1;
    else clist[i].flag = 0;
  }
}

/* ----------------------------------------------------------------------
   coarsen based on surfs in cell
   if any child cell smaller in any dimension than surfsize, coarsen to parent
   only use surfs with normals opposing vstream
------------------------------------------------------------------------- */

void AdaptGrid::coarsen_surf()
{
  int m,nchild,anysurf;
  double *values;

  for (int i = 0; i < cnum; i++) {
    if (clist[i].flag == 0) continue;

    values = clist[i].value;
    nchild = clist[i].nchild;

    anysurf = 0;
    for (m = 0; m < nchild; m++) anysurf += static_cast<int> (values[m]);
    if (anysurf) clist[i].flag = 1;
    else clist[i].flag = 0;
  }
}

/* ----------------------------------------------------------------------
   coarsen based on compute or fix value
------------------------------------------------------------------------- */

void AdaptGrid::coarsen_value()
{
  int m,nchild;
  double onevalue,allvalues;
  double *values;

  for (int i = 0; i < cnum; i++) {
    if (clist[i].flag == 0) continue;

    values = clist[i].value;
    nchild = clist[i].nchild;

    if (combine == SUM) allvalues = 0.0;
    else if (combine == MINIMUM) allvalues = BIG;
    else if (combine == MAXIMUM) allvalues = -BIG;

    for (m = 0; m < nchild; m++) {
      onevalue = values[m];
      if (combine == SUM) allvalues += onevalue;
      else if (combine == MINIMUM)
        allvalues = MIN(allvalues,onevalue);
      else if (combine == MAXIMUM)
        allvalues = MAX(allvalues,onevalue);
    }

    if (cdecide == LESS) {
      if (allvalues < cvalue) clist[i].flag = 1;
      else clist[i].flag = 0;
    } else {
      if (allvalues > cvalue) clist[i].flag = 1;
      else clist[i].flag = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   coarsen based on random #
   useful for debugging
------------------------------------------------------------------------- */

void AdaptGrid::coarsen_random()
{
  // loop over my coarsen list
  // skip unflagged tasks being considered by another proc
  // decide whether to coarsen or not, flag accordingly

  for (int i = 0; i < cnum; i++) {
    if (clist[i].flag == 0) continue;
    clist[i].flag = 0;
    if (random->uniform() >= cfrac) continue;
    clist[i].flag = 1;
  }
}

/* ----------------------------------------------------------------------
   now have clist with flag=1 for each parent cell that will be coarsened
   (1) 1st portion of clist can be coarsened by this proc (owns all child cells)
   (2) assign 2nd portion of clist to procs that own center child of each parent
   (3) rendezvous comm that owning proc to owners of all child cells in clist
   (4) child cells can then communicate their particles and surfs to new parent owner
       this step is performed by comm->send_cells_adapt()
   (5) create alist with info for only parent cell that will be coarsened
       using data received from step (4)
------------------------------------------------------------------------- */

void AdaptGrid::particle_surf_comm()
{
  int j,m,plevel,ihalf,jhalf,khalf,ichild,nchild,owner;
  int icell,jcell,np,nsplit;
  cellint parentID;
  int *csubs;

  Grid::ParentLevel *plevels = grid->plevels;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;

  // count child cells for parents flagged for coarsening
  // children will be notified via rendezvous comm, owned by this proc or others

  int nsend = 0;
  for (int i = 0; i < cnum; i++) {
    if (!clist[i].flag) continue;
    nsend += clist[i].nchild;
  }

  // assign clist parent cells to an owning proc
  // owning proc = one that owns center cell of child sub-grid
  // fill inbuf with info to send to each child cell of coarsened parents

  int *proclist;
  memory->create(proclist,nsend,"adapt_grid:proclist");
  Rvous2 *inbuf = (Rvous2 *) memory->smalloc((bigint) nsend*sizeof(Rvous2),
                                             "adapt_grid:inbuf");
  nsend = 0;

  for (int i = 0; i < cnum; i++) {
    if (!clist[i].flag) continue;
    plevel = clist[i].plevel;
    ihalf = plevels[plevel].nx / 2;
    jhalf = plevels[plevel].ny / 2;
    khalf = plevels[plevel].nz / 2;
    ichild = khalf*plevels[plevel].ny*plevels[plevel].nx +
      jhalf*plevels[plevel].nx + ihalf;
    owner = clist[i].proc[ichild];
    nchild = clist[i].nchild;

    for (m = 0; m < nchild; m++) {
      proclist[nsend] = clist[i].proc[m];
      inbuf[nsend].parentID = clist[i].parentID;
      inbuf[nsend].owner = owner;
      inbuf[nsend].icell = clist[i].index[m];
      inbuf[nsend].ichild = m;
      nsend++;
    }
  }

  // perform rendezvous communication to acquire inbuf data
  // each SPARTA proc receives child cell data from Rvous decomp for
  //   all children it owns in parent cells that are being coarsened
  // callback() method is NULL, just receive data via MPI_All2allv()

  char *buf;
  int nreturn = comm->rendezvous(1,nsend,(char *) inbuf,sizeof(Rvous2),
                                 0,proclist,NULL,0,buf,0,this,0);

  Rvous2 *outbuf = (Rvous2 *) buf;

  memory->destroy(proclist);
  memory->sfree(inbuf);

  // use outbuf to fill SendAdapt data struct for my requested child cells
  // flag my child cells for deletion

  nsend = nreturn;
  memory->create(proclist,nsend,"adapt_grid:proclist");
  SendAdapt *sadapt = (SendAdapt *) memory->smalloc(nsend*sizeof(SendAdapt),
                                                    "adapt_grid:sadapt");
  for (int i = 0; i < nreturn; i++) {
    icell = outbuf[i].icell;
    cells[icell].proc = -1;

    proclist[i] = outbuf[i].owner;
    sadapt[i].parentID = outbuf[i].parentID;
    sadapt[i].plevel = cells[icell].level - 1;
    sadapt[i].owner = outbuf[i].owner;
    sadapt[i].proc = me;
    sadapt[i].icell = icell;
    sadapt[i].type = cinfo[icell].type;
    sadapt[i].ichild = outbuf[i].ichild;
    sadapt[i].nsurf = cells[icell].nsurf;

    nsplit = cells[icell].nsplit;
    if (nsplit == 1) sadapt[i].np = cinfo[icell].count;
    else {
      csubs = sinfo[cells[icell].isplit].csubs;
      np = 0;
      for (int j = 0; j < nsplit; j++) {
        jcell = csubs[j];
        np += cinfo[jcell].count;
      }
      sadapt[i].np = np;
    }
  }

  memory->sfree(outbuf);

  // use Comm::send_cells_adapt()
  //   sends SendAdapt + surfs/particles for each cell to parent cell owner
  // it invokes Grid::pack_one_adapt()
  //   if this proc is parent owner, no surfs/particles are sent to self

  int nrecv = comm->send_cells_adapt(nsend,proclist,(char *) sadapt,&spbuf);

  memory->destroy(proclist);
  memory->sfree(sadapt);

  // create alist = list of parent cells I will coarsen in perform_coarsen()

  MyHash *alhash = new MyHash();
  SendAdapt *s;

  int nbytes_total = sizeof(Particle::OnePart) + particle->sizeof_custom();

  int dim = domain->dimension;
  int distributed = surf->distributed;

  alist = NULL;
  anum = anummax = 0;

  char *ptr = spbuf;

  for (int i = 0; i < nrecv; i++) {
    s = (SendAdapt *) ptr;
    parentID = s->parentID;

    if (alhash->find(parentID) == alhash->end()) {
      if (anum == anummax) {
        anummax += DELTA_LIST;
        alist = (ActionList *) memory->srealloc(alist,anummax*sizeof(ActionList),
                                                "adapt_grid:alist");
      }
      (*alhash)[parentID] = anum;
      m = anum++;
      alist[m].parentID = s->parentID;
      alist[m].plevel = s->plevel;
      alist[m].anyinside = 0;
      nchild = plevels[s->plevel].nxyz;
      alist[m].nchild = nchild;
      alist[m].index = new int[nchild];
      alist[m].nsurf = new int[nchild];
      alist[m].np = new int[nchild];
      alist[m].surfs = new void*[nchild];
      alist[m].particles = new char*[nchild];
    } else m = (*alhash)[parentID];

    if (s->type == INSIDE) alist[m].anyinside = 1;

    ichild = s->ichild;
    if (s->proc == me) alist[m].index[ichild] = s->icell;
    else alist[m].index[ichild] = -1;

    alist[m].np[ichild] = s->nsurf;
    alist[m].nsurf[ichild] = s->np;

    ptr += sizeof(SendAdapt);
    ptr = ROUNDUP(ptr);

    if (s->proc == me) continue;

    // surfs are packed differently for distributed vs non-distributed
    // non = indices, dist = line/tri data

    alist[m].nsurf[ichild] = s->nsurf;
    alist[m].surfs[ichild] = (void *) ptr;
    if (!distributed) ptr += s->nsurf * sizeof(surfint);
    else if (dim == 2) ptr += s->nsurf * sizeof(Surf::Line);
    else ptr += s->nsurf * sizeof(Surf::Tri);
    ptr = ROUNDUP(ptr);

    alist[m].np[ichild] = s->np;
    alist[m].particles[ichild] = (char *) ptr;
    ptr += s->np * nbytes_total;
    ptr = ROUNDUP(ptr);
  }

  // clean up
  // spbuf is onwed by Comm, so no need to delete in AdaptGrid

  delete alhash;
}

/* ----------------------------------------------------------------------
   perform coarsening of flagged parent cells in clist
------------------------------------------------------------------------- */

int AdaptGrid::perform_coarsen()
{
  int i,m,icell,nchild,newcell,mask;
  int plevel,nsplit,jcell,ip;
  cellint parentID;
  double plo[3],phi[3];
  int *csubs;

  Grid::ChildCell *cells;
  Grid::ChildInfo *cinfo;;
  Grid::SplitInfo *sinfo;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  // coarsening with distributed surfs requires use of hash
  // used by grid->coarsen_cell() method to add new unique surfs

  if (surf->distributed) surf->rehash();

  // loop over coarsening action list

  for (i = 0; i < anum; i++) {
    parentID = alist[i].parentID;
    plevel = alist[i].plevel;
    grid->id_lohi(parentID,plevel,boxlo,boxhi,plo,phi);
    nchild = alist[i].nchild;

    // coarsen parentID to become a new child cell

    grid->coarsen_cell(parentID,plevel,plo,phi,nchild,
                       alist[i].index,alist[i].nsurf,alist[i].np,
                       alist[i].surfs,alist[i].particles,cut2d,cut3d);

    cells = grid->cells;
    cinfo = grid->cinfo;
    sinfo = grid->sinfo;
    newcell = grid->nlocal - 1;

    // if new child has no surfs and any of its children was INSIDE
    // then type of new child cell = INSIDE

    if (cells[newcell].nsurf == 0 && alist[i].anyinside)
      cinfo[newcell].type = INSIDE;

    // set group mask of new child and its sub-cells to groupbit + all
    // NOTE: could set mask to intersection of all old children
    //       to attempt to preserve other group settings

    mask = groupbit | 1;
    cinfo[newcell].mask = mask;

    if (cells[newcell].nsplit > 1) {
      sinfo = grid->sinfo;
      nsplit = cells[newcell].nsplit;
      csubs = sinfo[cells[newcell].isplit].csubs;
      for (int j = 0; j < nsplit; j++) {
        jcell = csubs[j];
        cinfo[jcell].mask = mask;
      }
    }

    // add ID of new coarsened cell to chash

    (*chash)[parentID] = 0;
  }

  // done with surf hash

  if (surf->distributed) {
    surf->hash->clear();
    surf->hashfilled = 0;
  }

  // return # of cells this proc coarsened

  return anum;
}

/* ----------------------------------------------------------------------
   extract a value for icell,valindex from a compute
------------------------------------------------------------------------- */

double AdaptGrid::value_compute(int icell)
{
  double value;

  if (valindex == 0 || compute->post_process_grid_flag)
    value = compute->vector_grid[icell];
  else value = compute->array_grid[icell][valindex-1];

  return value;
}

/* ----------------------------------------------------------------------
   extract a value for icell,valindex from a fix
------------------------------------------------------------------------- */

double AdaptGrid::value_fix(int icell)
{
  double value;

  if (valindex == 0) value = fix->vector_grid[icell];
  else value = fix->array_grid[icell][valindex-1];

  return value;
}

/* ----------------------------------------------------------------------
   repoint particles in all refined grid cells to new icell index
   particles are still sorted
   so count/first values in cinfo are still valid
------------------------------------------------------------------------- */

int AdaptGrid::region_check(double *lo, double *hi)
{
  double pt[3];

  if (domain->dimension == 2) {
    if (regstyle == REGION_ALL) {
      pt[0] = lo[0]; pt[1] = lo[1]; pt[2] = 0.0;
      if (!region->match(pt)) return 0;
      pt[0] = hi[0]; pt[1] = lo[1]; pt[2] = 0.0;
      if (!region->match(pt)) return 0;
      pt[0] = lo[0]; pt[1] = hi[1]; pt[2] = 0.0;
      if (!region->match(pt)) return 0;
      pt[0] = hi[0]; pt[1] = hi[1]; pt[2] = 0.0;
      if (!region->match(pt)) return 0;
      return 1;
    } else if (regstyle == REGION_ONE) {
      pt[0] = lo[0]; pt[1] = lo[1]; pt[2] = 0.0;
      if (region->match(pt)) return 1;
      pt[0] = hi[0]; pt[1] = lo[1]; pt[2] = 0.0;
      if (region->match(pt)) return 1;
      pt[0] = lo[0]; pt[1] = hi[1]; pt[2] = 0.0;
      if (region->match(pt)) return 1;
      pt[0] = hi[0]; pt[1] = hi[1]; pt[2] = 0.0;
      if (region->match(pt)) return 1;
      return 0;
    } else if (regstyle == REGION_CENTER) {
      pt[0] = 0.5 * (lo[0] + hi[0]);
      pt[1] = 0.5 * (lo[1] + hi[1]);
      pt[2] = 0.0;
      if (region->match(pt)) return 1;
      return 0;
    }
  }

  if (domain->dimension == 3) {
    if (regstyle == REGION_ALL) {
      pt[0] = lo[0]; pt[1] = lo[1]; pt[2] = lo[2];
      if (!region->match(pt)) return 0;
      pt[0] = hi[0]; pt[1] = lo[1]; pt[2] = lo[2];
      if (!region->match(pt)) return 0;
      pt[0] = lo[0]; pt[1] = hi[1]; pt[2] = lo[2];
      if (!region->match(pt)) return 0;
      pt[0] = hi[0]; pt[1] = hi[1]; pt[2] = lo[2];
      if (!region->match(pt)) return 0;
      pt[0] = lo[0]; pt[1] = lo[1]; pt[2] = hi[2];
      if (!region->match(pt)) return 0;
      pt[0] = hi[0]; pt[1] = lo[1]; pt[2] = hi[2];
      if (!region->match(pt)) return 0;
      pt[0] = lo[0]; pt[1] = hi[1]; pt[2] = hi[2];
      if (!region->match(pt)) return 0;
      pt[0] = hi[0]; pt[1] = hi[1]; pt[2] = hi[2];
      if (!region->match(pt)) return 0;
      return 1;
    } else if (regstyle == REGION_ONE) {
      pt[0] = lo[0]; pt[1] = lo[1]; pt[2] = lo[2];
      if (region->match(pt)) return 1;
      pt[0] = hi[0]; pt[1] = lo[1]; pt[2] = lo[2];
      if (region->match(pt)) return 1;
      pt[0] = lo[0]; pt[1] = hi[1]; pt[2] = lo[2];
      if (region->match(pt)) return 1;
      pt[0] = hi[0]; pt[1] = hi[1]; pt[2] = lo[2];
      if (region->match(pt)) return 1;
      pt[0] = lo[0]; pt[1] = lo[1]; pt[2] = hi[2];
      if (region->match(pt)) return 1;
      pt[0] = hi[0]; pt[1] = lo[1]; pt[2] = hi[2];
      if (region->match(pt)) return 1;
      pt[0] = lo[0]; pt[1] = hi[1]; pt[2] = hi[2];
      if (region->match(pt)) return 1;
      pt[0] = hi[0]; pt[1] = hi[1]; pt[2] = hi[2];
      if (region->match(pt)) return 1;
      return 0;
    } else if (regstyle == REGION_CENTER) {
      pt[0] = 0.5 * (lo[0] + hi[0]);
      pt[1] = 0.5 * (lo[1] + hi[1]);
      pt[2] = 0.5 * (lo[2] + hi[2]);
      if (region->match(pt)) return 1;
      return 0;
    }
  }

  return 0;
}

/* ----------------------------------------------------------------------
   cleanup memory use
------------------------------------------------------------------------- */

void AdaptGrid::cleanup()
{
  delete random;
  delete [] childlist;

  memory->destroy(rlist);
  if (clist) {
    for (int i = 0; i < cnum; i++) {
      delete [] clist[i].proc;
      delete [] clist[i].index;
      delete [] clist[i].value;
    }
    memory->sfree(clist);
  }
  if (alist) {
    for (int i = 0; i < anum; i++) {
      delete [] alist[i].index;
      delete [] alist[i].nsurf;
      delete [] alist[i].np;
      delete [] alist[i].surfs;
      delete [] alist[i].particles;
    }
    memory->sfree(alist);
  }

  if (domain->dimension == 3) delete cut3d;
  else delete cut2d;

  delete rhash;
  delete chash;
}

/* ----------------------------------------------------------------------
   write out new grid file via WriteGrid
------------------------------------------------------------------------- */

void AdaptGrid::write_file()
{
  WriteGrid *wg = new WriteGrid(sparta);
  wg->silent = 1;

  int narg = 1;
  char **args = new char*[narg];

  char *expandfile = NULL;
  if (strchr(file,'*')) {
    expandfile = new char[strlen(file) + 16];
    char *ptr = strchr(file,'*');
    *ptr = '\0';
    sprintf(expandfile,"%s" BIGINT_FORMAT "%s",file,update->ntimestep,ptr+1);
    *ptr = '*';
    args[0] = expandfile;
  } else args[0] = file;

  wg->command(narg,args);

  // NOTE: could persist WriteGrid instance for fix adapt

  delete [] expandfile;
  delete [] args;
  delete wg;
}
