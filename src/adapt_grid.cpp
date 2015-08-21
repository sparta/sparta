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
#include "random_park.h"
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
#define DELTA_NEW 1024
#define DELTA_SEND 1024
#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

AdaptGrid::AdaptGrid(SPARTA *sparta) : Pointers(sparta)
{
  me = comm->me;
  nprocs = comm->nprocs;

  valueID = NULL;
  maxnew = 0;
  newcells = NULL;
  file = NULL;
}

/* ---------------------------------------------------------------------- */

AdaptGrid::~AdaptGrid()
{
  delete [] valueID;
  memory->destroy(newcells);
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

  /*
  printf("PRE ADAPT\n",grid->nlocal);
  for (int i = 0; i < grid->nlocal; i++) {
    Grid::ChildCell *g = &grid->cells[i];
    if (g->nsplit <= 1)
    printf("ICELL %d id %d iparent %d proc %d nsplit %d isplit %d lo %g %g "
           "hi %g %g inout %d parts %d %d\n",
           i,g->id,g->iparent,g->proc,g->nsplit,g->isplit,
           g->lo[0],
           g->lo[1],
           g->hi[0],
           g->hi[1],
           grid->cinfo[i].type,grid->cinfo[i].first,grid->cinfo[i].count);
    else
    printf("ICELLSPLIT %d id %d iparent %d proc %d nsplit %d isplit %d "
           "spicell %d csubs %d %d lo %g %g "
           "hi %g %g inout %d parts %d %d\n",
           i,g->id,g->iparent,g->proc,g->nsplit,g->isplit,
           grid->sinfo[g->isplit].icell,
           grid->sinfo[g->isplit].csubs[0],grid->sinfo[g->isplit].csubs[1],
           g->lo[0],
           g->lo[1],
           g->hi[0],
           g->hi[1],
           grid->cinfo[i].type,grid->cinfo[i].first,grid->cinfo[i].count);
  }
  */

  // iterate over refinement and coarsening actions
  // use of shrinking pstop prevents any new parent created by refinement
  //   ever being considered for coarsening
  // use of chash prevents any new child created by coarsening
  //   ever being considered for refinement

  bigint nrefine_total = 0;
  bigint ncoarsen_total = 0;
  int nrefine = 0;
  int ncoarsen = 0;

  int pstop = grid->nparent;

  for (int iter = 0; iter < niterate; iter++) {
    setup(iter);

    if (action1 == REFINE) nrefine = refine();
    else if (action1 == COARSEN) ncoarsen = coarsen(pstop);

    if (action2 == REFINE) nrefine = refine();
    else if (action2 == COARSEN) ncoarsen = coarsen(pstop);

    if (nrefine == 0 && ncoarsen == 0) break;
    nrefine_total += nrefine;
    ncoarsen_total += ncoarsen;
    pstop -= ncoarsen;

    // NOTE: if mode = value, need to trigger computes to realloc to new grid
    //       for subsequent iterations - done below at very end

    cleanup();
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

  // DEBUG

  /*
  char str[32];
  printf("POST GATHER %d: %d\n",comm->me,grid->nparent);
  for (int i = 0; i < grid->nparent; i++) {
    Grid::ParentCell *p = &grid->pcells[i];
    grid->id_num2str(p->id,str);
    printf("PCELL %d: %d id %d %s iparent %d level %d "
           "nxyz %d %d %d lo %g %g %g "
           "hi %g %g %g\n",
           comm->me,i,p->id,str,p->iparent,p->level,
           p->nx,p->ny,p->nz,
           p->lo[0],
           p->lo[1],
           p->lo[2],
           p->hi[0],
           p->hi[1],
           p->hi[2]);
  }

  printf("POST ADAPT %d: %d\n",comm->me,grid->nlocal);
  for (int i = 0; i < grid->nlocal; i++) {
    Grid::ChildCell *g = &grid->cells[i];
    if (g->nsplit <= 1)
    printf("ICELL %d: %d id %d iparent %d proc %d nsplit %d isplit %d "
           "lo %g %g %g "
           "hi %g %g %g parts %d %d\n",
           comm->me,i,g->id,g->iparent,g->proc,g->nsplit,g->isplit,
           g->lo[0],
           g->lo[1],
           g->lo[2],
           g->hi[0],
           g->hi[1],
           g->hi[2],grid->cinfo[i].first,grid->cinfo[i].count);
    else
    printf("ICELLSPLIT %d id %d iparent %d proc %d nsplit %d isplit %d "
           "spicell %d csubs %d %d lo %g %g "
           "hi %g %g inout %d\n",
           i,g->id,g->iparent,g->proc,g->nsplit,g->isplit,
           grid->sinfo[g->isplit].icell,
           grid->sinfo[g->isplit].csubs[0],grid->sinfo[g->isplit].csubs[1],
           g->lo[0],
           g->lo[1],
           g->hi[0],
           g->hi[1],
           grid->cinfo[i].type);
  }
  */

  MPI_Barrier(world);
  double time2 = MPI_Wtime();

  // reset all attributes of adapted grid
  // same steps as in create_grid

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->find_neighbors();
  grid->check_uniform();
  comm->reset_neighbors();

  // DEBUG

  /*
  printf("PRE INOUT %d: %d\n",comm->me,grid->nlocal);
  Grid::ParentCell *pcells = grid->pcells;
  for (int i = 0; i < grid->nlocal; i++) {
    Grid::ChildCell *g = &grid->cells[i];
    printf("ICELL %d: %d id %d pid %d lo %g %g "
           "hi %g %g type %d corners %d %d %d %d\n",
           comm->me,i,g->id,pcells[g->iparent].id,
           g->lo[0],
           g->lo[1],
           g->hi[0],
           g->hi[1],
           grid->cinfo[i].type,
           grid->cinfo[i].corner[0],
           grid->cinfo[i].corner[1],
           grid->cinfo[i].corner[2],
           grid->cinfo[i].corner[3]);
  }
  */

  if (surf->exist) {
    grid->set_inout();
    grid->type_check();
  }

  // final update of any per grid fixes for all new child cells
  
  if (modify->n_pergrid) add_grid_fixes();

  // reallocate per grid arrays in per grid dumps

  for (int i = 0; i < output->ndump; i++)
    output->dump[i]->reset_grid();

  // write out new parent grid file

  if (file) write_file();

  /*
  printf("POST INOUT %d: %d\n",comm->me,grid->nlocal);
  Grid::ParentCell *pcells = grid->pcells;
  for (int i = 0; i < grid->nlocal; i++) {
    Grid::ChildCell *g = &grid->cells[i];
    printf("ICELL %d: %d id %d pid %d lo %g %g "
           "hi %g %g type %d corners %d %d %d %d vol %g\n",
           comm->me,i,g->id,pcells[g->iparent].id,
           g->lo[0],
           g->lo[1],
           g->hi[0],
           g->hi[1],
           grid->cinfo[i].type,
           grid->cinfo[i].corner[0],
           grid->cinfo[i].corner[1],
           grid->cinfo[i].corner[2],
           grid->cinfo[i].corner[3],grid->cinfo[i].volume);
  }
  */

  MPI_Barrier(world);
  double time3 = MPI_Wtime();

  // stats

  double time_total = time3-time1;

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  " BIGINT_FORMAT " cells refined, " BIGINT_FORMAT 
              " cells coarsened\n",nrefine_total,ncoarsen_total);
      fprintf(screen,"  adapted to " BIGINT_FORMAT " child grid cells, "
              "%d parent cells\n",grid->ncell,grid->nparent);
      fprintf(screen,"  CPU time = %g secs\n",time_total);
      fprintf(screen,"  adapt/redo percent = %g %g\n",
              100.0*(time2-time1)/time_total,100.0*(time3-time2)/time_total);
    }

    if (logfile) {
      fprintf(logfile,"  " BIGINT_FORMAT " cells refined, " BIGINT_FORMAT 
              " cells coarsened\n",nrefine_total,ncoarsen_total);
      fprintf(logfile,"  adapted to " BIGINT_FORMAT " child grid cells, "
              "%d parent cells\n",grid->ncell,grid->nparent);
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
  if (narg < 2) error->all(FLERR,"Illegal adapt command");

  // define action(s)

  if (strcmp(arg[0],"refine") == 0) action1 = REFINE;
  else if (strcmp(arg[0],"coarsen") == 0) action1 = COARSEN;
  else error->all(FLERR,"Illegal adapt command");
  int iarg = 1;
  
  if (strcmp(arg[1],"refine") == 0) {
    action2 = REFINE;
    iarg = 2;
  } else if (strcmp(arg[1],"coarsen") == 0) {
    action2 = COARSEN;
    iarg = 2;
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
        error->all(FLERR,"Adapt cells nz must be 1 for 2d simulation");
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
}

/* ----------------------------------------------------------------------
   error check on value compute/fix for both adapt_grid and fix adapt
   nevery is passed only from fix adapt (mode = 1)
------------------------------------------------------------------------- */

void AdaptGrid::check_args(int nevery)
{
  // for every fix ave/grid defined, require that:
  // (1) fix adapt Nevery is multiple of fix ave Nfreq
  // (2) fix ave/grid is not a running ave
  // (3) fix ave/grid is defined before this fix (checked in fix adapt)
  // this insures fix ave/grid values will be up-to-date before
  //   this adaptation changes the grid
  // NOTE: at some point could make fix ave/grid do full interpolation
  //       for new grid cells so this restriction is not needed

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
        error->all(FLERR,"Adapt can not perform multiple iterations "
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
  // zero new child cell counter on first iteration
  // used in add_grid_fixes after all new child cells have been added

  if (iter == 0) nnew = 0;

  // create RNG for style = RANDOM

  if (style == RANDOM) {
    random = new RanPark(update->ranmaster->uniform());
    double seed = update->ranmaster->uniform();
    random->reset(seed,comm->me,100);
  } else random = NULL;

  // list of new cell indices for one refined cell

  childlist = new int[nx*ny*nz];

  // for cut/split of new cells by surfaces

  cut3d = NULL;
  cut2d = NULL;
  if (surf->exist) {
    if (domain->dimension == 3) cut3d = new Cut3d(sparta);
    else cut2d = new Cut2d(sparta,domain->axisymmetric);
  }

  // data structs

  rlist = NULL;
  ctask = NULL;
  sadapt = NULL;
  delparent = NULL;

  // allocate chash for child cell IDs created by coarsening a parent

#ifdef SPARTA_MAP
  chash = new std::map<cellint,int>();
#elif defined SPARTA_UNORDERED_MAP
  chash = new std::unordered_map<cellint,int>();
#else
  chash = new std::tr1::unordered_map<cellint,int>();
#endif
}

/* ----------------------------------------------------------------------
   perform refinement
------------------------------------------------------------------------- */

int AdaptGrid::refine()
{
  // NOTE: why is grid not really hashed at this point even though it says so?
  // WHEN done with adapt, do I need to mark it hashed or unhashed?
  // NOTE: does fix balance always do a sort?  yes it does

  //if (!grid->hashfilled) grid->rehash();
  //if (particle->exist && !particle->sorted) particle->sort();      

  grid->rehash();
  particle->sort();      

  candidates_refine();

  if (style == PARTICLE) refine_particle();
  else if (style == SURF) refine_surf();
  else if (style == VALUE) refine_value();
  else if (style == RANDOM) refine_random();

  int delta = perform_refine();

  // nrefine = # of new parents across all processors

  int nrefine;
  MPI_Allreduce(&delta,&nrefine,1,MPI_INT,MPI_SUM,world);

  if (nrefine) {
    gather_parents_refine(delta,nrefine);
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

int AdaptGrid::coarsen(int pstop)
{
  // NOTE: why is grid not really hashed at this point even though it says so?
  // WHEN done with adapt, do I need to mark it hashed or unhashed?
  // NOTE: does fix balance always do a sort?  yes it does

  //if (!grid->hashfilled) grid->rehash();
  //if (particle->exist && !particle->sorted) particle->sort();      

  grid->rehash();
  particle->sort();      

  assign_parents_coarsen(pstop);
  candidates_coarsen(pstop);

  if (style == PARTICLE) coarsen_particle();
  else if (style == SURF) coarsen_surf();
  else if (style == VALUE) coarsen_value();
  else if (style == RANDOM) coarsen_random();

  int delta = perform_coarsen();

  // ncoarsen = # of removed parents across all processors

  int ncoarsen;
  MPI_Allreduce(&delta,&ncoarsen,1,MPI_INT,MPI_SUM,world);

  if (ncoarsen) {
    gather_parents_coarsen(delta,ncoarsen);
    if (collide) collide->adapt_grid();
    grid->compress();
    if (particle->exist && replyany) particle->compress_rebalance();

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
   observe maxlevel
   chash check skips child cells created by coarsening a parent cell
------------------------------------------------------------------------- */

void AdaptGrid::candidates_refine()
{
  Grid::ParentCell *pcells = grid->pcells;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  memory->create(rlist,nglocal,"adapt_grid:rlist");
  rnum = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (maxlevel && pcells[cells[icell].iparent].level >= maxlevel-1) continue;
    if (cinfo[icell].type == INSIDE) continue;
    if (region && !region_check(cells[icell].lo,cells[icell].hi)) continue;
    if (chash->find(cells[icell].id) != chash->end()) continue;
    rlist[rnum++] = icell;
  }

  //printf("RNUM %d\n",rnum);
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
      for (int i = 0; i < nsplit; i++) {
        jcell = csubs[i];
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
  int m,icell,flag,nsurf;
  int *csurfs;
  double *norm,*lo,*hi;

  int dim = domain->dimension;
  Grid::ChildCell *cells = grid->cells;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  int n = 0;
  for (int i = 0; i < rnum; i++) {
    icell = rlist[i];
    if (!cells[icell].nsurf) continue;
    nsurf = cells[icell].nsurf;
    csurfs = cells[icell].csurfs;
    for (m = 0; m < nsurf; m++) {
      if (dim == 2) {
        if (!(lines[m].mask & sgroupbit)) continue;
      } else {
        if (!(tris[m].mask & sgroupbit)) continue;
      }
      if (dim == 2) norm = lines[csurfs[m]].norm;
      else norm = tris[csurfs[m]].norm;
      if (MathExtra::dot3(norm,sdir) < 0.0) break;
    }
    if (m == nsurf) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    flag = 1;
    if (fabs(hi[0]-lo[0])/nx < surfsize) flag = 0;
    if (fabs(hi[1]-lo[1])/ny < surfsize) flag = 0;
    if (dim == 3 && fabs(hi[2]-lo[2])/nz < surfsize) flag = 0;
    if (flag) rlist[n++] = icell;
  }
  rnum = n;

  //printf("RNUM2 %d: %d\n",me,rnum);
}

/* ----------------------------------------------------------------------
   refine based on compute or fix value
   // NOTE: insure fix ave/grid comes before fix adapt at end of step
------------------------------------------------------------------------- */

void AdaptGrid::refine_value()
{
  int icell,nsplit,jcell;
  double value;
  int *csubs;

  // NOTE: always invoke compute - is there someway to check more carefully?
  // needed if iterating, needed if refine follows coarsen
  
  if (valuewhich == COMPUTE) {
    compute = modify->compute[icompute];
    compute->compute_per_grid();
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
    }

    //printf("AAA %d %d: %g %g\n",
    //       icell,cells[icell].id,value,rvalue);

    if (rdecide == LESS) {
      if (value < rvalue) rlist[n++] = icell;
    } else {
      //printf("AAA %d %d: %g %g\n",icell,cells[icell].id,value,rthresh);
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
    //if (grid->cells[icell].id != 6) continue;
    rlist[n++] = icell;
  }
  rnum = n;

  //printf("RNUM2 %d: %d\n",me,rnum);
}

/* ----------------------------------------------------------------------
   final add of new child cells to fixes
   allow for each cell ID to no longer exist or be for a parent cell
     due to subsequent iterations
   NOTE: first stage must happen earlier so fix memory for new cell is valid
------------------------------------------------------------------------- */

void AdaptGrid::add_grid_fixes()
{
#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  grid->rehash();

  int icell;
  for (int i = 0; i < nnew; i++) {
    if (hash->find(newcells[i]) == hash->end()) continue;
    icell = (*hash)[newcells[i]];
    if (icell < 0) continue;
    modify->add_grid_one(icell-1,1);
  }
}

/* ----------------------------------------------------------------------
   perform refinement of child cells in rlist
------------------------------------------------------------------------- */

int AdaptGrid::perform_refine()
{
  int icell,jcell,iparent;
  int nglocal,nglocalprev;
  
  Grid::ParentCell *pcells = grid->pcells;
  Grid::ChildCell *cells = grid->cells;

  int nxyz = nx*ny*nz;
  int newparent = 0;
  
  for (int ilist = 0; ilist < rnum; ilist++) {
    icell = rlist[ilist];

    // flag icell for deletion
    // refine_cell() will remove it from grid hash
    // don't need to flag sub cells, since grid compress will do that

    cells[icell].proc = -1;

    // add new parent cell at end of my pcells

    grid->add_parent_cell(cells[icell].id,cells[icell].iparent,
                          nx,ny,nz,cells[icell].lo,cells[icell].hi);
    newparent++;
    pcells = grid->pcells;
    iparent = grid->nparent - 1;

    // add Nx by Ny by Nz child cells to replace icell
    // refine_cell() adds all new unsplit/split cells to grid hash

    nglocalprev = grid->nlocal;
    grid->refine_cell(icell,iparent,nx,ny,nz,childlist,cut2d,cut3d);
    cells = grid->cells;
    nglocal = grid->nlocal;

    // add each new child cell to newcells list
    // only for unsplit/split cells, not sub cells

    if (nnew + nxyz > maxnew) {
      while (nnew + nxyz > maxnew) maxnew += DELTA_NEW;
      memory->grow(newcells,maxnew,"adapt_grid:newcells");
    }

    for (jcell = nglocalprev; jcell < nglocal; jcell++) {
      if (cells[jcell].nsplit < 1) continue;
      newcells[nnew++] = cells[jcell].id;
    }
  }

  // return # of new parents on this proc

  return newparent;
}

/* ----------------------------------------------------------------------
   gather refined parent cells from all procs, on each proc
   must end up with same list, in same order, on every proc
   must repoint child cells to new parent
   also point child cells that became parent to new parent
     this is so can reassign particles to correct new child cells
------------------------------------------------------------------------- */

void AdaptGrid::gather_parents_refine(int delta, int nrefine)
{
  int i,j,m,icell,nx,ny,nz,nsplit;
  cellint id;
  Grid::ParentCell *p;

  // pcell_mine = copy of parent cells I added via refinement

  Grid::ParentCell *pcells = grid->pcells;
  int psize = sizeof(Grid::ParentCell);
  int nprev = grid->nparent - delta;

  Grid::ParentCell *pcells_mine = (Grid::ParentCell *)
    memory->smalloc(delta*psize,"adapt_grid:pcell_mine");
  memcpy(pcells_mine,&pcells[nprev],delta*psize);

  // perform Allgatherv of added parent cells from all procs
  // append to end of original pcells list

  int *recvcounts,*displs;
  memory->create(recvcounts,nprocs,"adapt_grid:recvcounts");
  memory->create(displs,nprocs,"adapt_grid:displs");

  int nsend = delta * psize;
  MPI_Allgather(&nsend,1,MPI_INT,recvcounts,1,MPI_INT,world);
  displs[0] = 0;
  for (i = 1; i < nprocs; i++) displs[i] = displs[i-1] + recvcounts[i-1];

  grid->grow_pcells(nrefine-delta);
  pcells = grid->pcells;

  MPI_Allgatherv(pcells_mine,nsend,MPI_CHAR,
                 &pcells[nprev],recvcounts,displs,MPI_CHAR,world);
  grid->nparent = nprev + nrefine;
  
  // loop over all added pcells
  // repoint any child cells I own to new parent

  int myfirst = nprev + displs[me]/psize;
  int mylast = myfirst + recvcounts[me]/psize;

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  int nparent = grid->nparent;
  int *csubs;

  //printf("GATHER: %d %d %d\n",nprev,delta,nrefine);
  //printf("GATHER: myfirst mylast %d %d\n",myfirst,mylast);

  for (i = nprev; i < nparent; i++) {
    p = &pcells[i];
    (*hash)[p->id] = -(i+1);
    pcells[p->iparent].grandparent = 1;

    // not a new parent cell I contributed, just continue

    if (i < myfirst || i >= mylast) continue;

    // this is new parent cell I contributed
    // repoint new child cells of this parent to new iparent value
    // do this for sub cells of new child cells as well

    nx = p->nx;
    ny = p->ny;
    nz = p->nz;

    m = 0;
    for (int iz = 0; iz < nz; iz++)
      for (int iy = 0; iy < ny; iy++)
        for (int ix = 0; ix < nx; ix++) {
          m++;
          id = p->id | ((cellint) m << p->nbits);
          icell = (*hash)[id] - 1;
          cells[icell].iparent = i;

          if (cells[icell].nsplit > 1) {
            nsplit = cells[icell].nsplit;
            csubs = sinfo[cells[icell].isplit].csubs;
            for (j = 0; j < nsplit; j++) cells[csubs[j]].iparent = i;
          }
        }
  }

  // clean up

  memory->sfree(pcells_mine);
  memory->destroy(recvcounts);
  memory->destroy(displs);
}

/* ----------------------------------------------------------------------
   assign each parent to an owning proc
   for non-clumped distribution use rendezvous proc = round robin
   for clumped distribution assign to proc owning central child cell
     userendezvous comm to figure that out
   observe minlevel
------------------------------------------------------------------------- */

void AdaptGrid::assign_parents_coarsen(int pstop)
{
  int i,m,n,icell,iparent,nxyz;
  cellint id;

  Grid::ParentCell *pcells = grid->pcells;
  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  // pcount = number of children I own for each parent
  // only tabulate for eligible parents < pstop
  // only tally unsplit/split child cells, not sub cells

  memory->create(pcount,pstop,"adapt_grid:pcount");
  for (i = 0; i < pstop; i++) pcount[i] = 0;

  for (icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit < 1) continue;
    iparent = cells[icell].iparent;
    if (iparent >= pstop) continue;
    pcount[iparent]++;
  }

  memory->create(powner,pstop,"adapt_grid:powner");
  for (i = 0; i < pstop; i++) powner[i] = -1;

  for (i = 0; i < pstop; i++) {
    if (pcells[i].grandparent) continue;
    if (pcells[i].level < minlevel) continue;
    if (region && !region_check(pcells[i].lo,pcells[i].hi)) continue;

    nxyz = pcells[i].nx * pcells[i].ny * pcells[i].nz;
    if (pcount[i] == nxyz) powner[i] = me;
    else powner[i] = i % nprocs;
  }

  cnummax = 0;
  for (i = 0; i < pstop; i++) 
    if (powner[i] == me) cnummax++;

  if (nprocs == 1 || !grid->clumped) return;

  // to insure grid decomposition stays clumped,
  //   need to insure any coarsened parents are assigned to a proc
  //   that owns at least one of its children
  // redezvous procs do not do that
  // so use rendezvous procs to choose an owner and communicate
  //   the owner to all the child cells of the coarsened parent

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  // send 4 values per child cell to rendezvous proc
  // sending proc, local icell, iparent, ichild (0 to Nxyz-1)

  int *sbuf = NULL;
  int *procsend = NULL;
  int sendmax = 0;
  int nsend = 0;

  n = 0;
  for (i = 0; i < pstop; i++) {
    if (pcount[i] == 0) continue;
    if (powner[i] < 0 || powner[i] == me) continue;
    nxyz = pcells[i].nx * pcells[i].ny * pcells[i].nz;

    for (m = 0; m < nxyz; m++) {
      id = pcells[i].id | ((m+1) << pcells[i].nbits);
      if (hash->find(id) == hash->end()) continue;
      icell = (*hash)[id] - 1;
      
      if (nsend == sendmax) {
        sendmax += DELTA_SEND;
        memory->grow(sbuf,4*sendmax,"adapt_grid:sbuf");
        memory->grow(procsend,sendmax,"adapt_grid:procsend");
      }
      
      procsend[nsend] = powner[i];
      sbuf[n] = me;
      sbuf[n+1] = icell;
      sbuf[n+2] = i;
      sbuf[n+3] = m;
      nsend++;
      n += 4;
    }
  }

  // perform irregular comm

  int nsize = 4*sizeof(int);
  int *rbuf;
  int nrecv = comm->irregular_uniform(nsend,procsend,(char *) sbuf,
                                      nsize,(char **) &rbuf);

  // scan received child cells,
  // looking for one which determines owner of parent cell
  // NOTE: make it central child of parent

  m = 0;
  for (i = 0; i < nrecv; i++) {
    if (rbuf[m+3] == 0) powner[rbuf[m+2]] = rbuf[m];
    m += 4;
  }

  // send new owning proc for parent of child cells back to senders

  memory->grow(sbuf,2*nrecv,"adapt_grid:sbuf");
  memory->grow(procsend,nrecv,"adapt_grid:procsend");
  nsend = 0;

  m = n = 0;
  for (i = 0; i < nrecv; i++) {
    procsend[nsend] = rbuf[m];
    sbuf[n] = rbuf[m+1];
    //    printf("AAA %d %d %d: %d %d %d %d\n",
    //       me,i,nrecv,rbuf[m],rbuf[m+1],rbuf[m+2],rbuf[m+3]);
    sbuf[n+1] = powner[rbuf[m+2]];
    nsend++;
    n += 2;
    m += 4;
  }

  // perform irregular comm

  nsize = 2*sizeof(int);
  nrecv = comm->irregular_uniform(nsend,procsend,(char *) sbuf,
                                  nsize,(char **) &rbuf);

  // scan received child cells and set powner of their parent

  m = 0;
  for (i = 0; i < nrecv; i++) {
    iparent = cells[rbuf[m]].iparent;
    powner[iparent] = rbuf[m+1];
    m += 2;
  }

  // udpate cnummax to reflect new parent owners

  cnummax = 0;
  for (i = 0; i < pstop; i++) 
    if (powner[i] == me) cnummax++;

  // clean up

  memory->destroy(sbuf);
  memory->destroy(procsend);
}

/* ----------------------------------------------------------------------
   create list of parent cells I will coarsen
------------------------------------------------------------------------- */

void AdaptGrid::candidates_coarsen(int pstop)
{
  int i,m,icell,jcell,iparent,nxyz,oproc;
  cellint id;
  double value;
  int *proc,*index;

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  // for style = VALUE, setup compute or fix values
  // NOTE: always invoke compute - is there someway to check more carefully?
  // needed if iterating, needed if refine follows coarsen

  if (style == VALUE) {
    if (valuewhich == COMPUTE) {
      compute = modify->compute[icompute];
      compute->compute_per_grid();
    } else if (valuewhich == FIX) fix = modify->fix[ifix];
  }

  Grid::ParentCell *pcells = grid->pcells;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;

  // create ctask entries = list of parents I will possibly coarsen
  // create sadapt = list of child cells I will send to other procs

  ctask = (CTask *) memory->smalloc(cnummax*sizeof(CTask),"adapt_grid:ctask");
  cnum = 0;

  int *parent2task;
  memory->create(parent2task,pstop,"adapt_grid:parent2task");
  for (i = 0; i < pstop; i++) parent2task[i] = -1;

  int *procsend = NULL;
  int sendmax = 0;
  nsend = 0;

  for (i = 0; i < pstop; i++) {

    if (powner[i] < 0) continue;
    nxyz = pcells[i].nx * pcells[i].ny * pcells[i].nz;

    // add to clist, b/c I am owner

    if (powner[i] == me) {
      //printf("AAA %d: %d %d\n",me,i,pcells[i].id);
      ctask[cnum].iparent = i;
      ctask[cnum].id = pcells[i].id;
      ctask[cnum].nchild = nxyz;
      ctask[cnum].ncomplete = pcount[i];
      proc = ctask[cnum].proc = new int[nxyz];
      index = ctask[cnum].index = new int[nxyz];
      ctask[cnum].recv = new int[nxyz];
      parent2task[i] = cnum;
      cnum++;

      if (pcount[i]) {
        for (m = 0; m < nxyz; m++) {
          id = pcells[i].id | ((m+1) << pcells[i].nbits);
          if (hash->find(id) == hash->end()) proc[m] = index[m] = -1;
          else {
            icell = (*hash)[id] - 1;
            proc[m] = me;
            index[m] = icell;
          }
        }
      }

    // another proc is the owner
    // fill sadapt with info on child cells I own to send to that proc

    } else if (pcount[i]) {
      //printf("CCC %d: %d %d\n",me,i,pcells[i].id);
      oproc = powner[i];
      for (m = 0; m < nxyz; m++) {
        id = pcells[i].id | ((m+1) << pcells[i].nbits);
        if (hash->find(id) == hash->end()) continue;
        icell = (*hash)[id] - 1;

        if (nsend == sendmax) {
          sendmax += DELTA_SEND;
          sadapt = (SendAdapt *) 
            memory->srealloc(sadapt,sendmax*sizeof(SendAdapt),
                             "adapt_grid:sadapt");
          memory->grow(procsend,sendmax,"adapt_grid:procsend");
        }

        procsend[nsend] = oproc;
        sadapt[nsend].proc = me;
        sadapt[nsend].icell = icell;
        sadapt[nsend].ichild = m;
        sadapt[nsend].iparent = i;
        sadapt[nsend].type = cinfo[icell].type;
        sadapt[nsend].nsurf = cells[icell].nsurf;
        if (cells[icell].nsplit == 1) {
          sadapt[nsend].np = cinfo[icell].count;
          if (style == VALUE) {
            if (valuewhich == COMPUTE) value = value_compute(icell);
            else if (valuewhich == FIX) value = value_fix(icell);
            sadapt[nsend].value = value;
          } else sadapt[nsend].value = 0.0;

        } else {
          int nsplit = cells[icell].nsplit;
          int isplit = cells[icell].isplit;
          int *csubs = sinfo[isplit].csubs;
          int np = 0;
          if (combine == SUM) value = 0.0;
          else if (combine == MINIMUM) value = BIG;
          else if (combine == MAXIMUM) value = -BIG;
          for (int i = 0; i < nsplit; i++) {
            jcell = csubs[i];
            np += cinfo[jcell].count;
            if (style == VALUE) {
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
          sadapt[nsend].np = np;
          if (style == VALUE) sadapt[nsend].value = value;
          else sadapt[nsend].value = 0.0;
        }

        nsend++;
      }
    }
  }

  /*
  printf("NCOMPLETE %d:",me);
  for (i = 0; i < cnum; i++) printf(" %d",ctask[i].ncomplete);
  printf("\n");
  */

  // clean up

  memory->destroy(pcount);
  memory->destroy(powner);

  // return if no one has any child cell info to send to another proc

  int sendany;
  MPI_Allreduce(&nsend,&sendany,1,MPI_INT,MPI_SUM,world);
  if (!sendany) {
    for (i = 0; i < cnum; i++) {
      if (ctask[i].ncomplete) ctask[i].flag = 1;
      else ctask[i].flag = 0;
    }
    memory->destroy(parent2task);
    memory->destroy(procsend);
    sa_header = NULL;
    sa_csurfs = NULL;
    sa_particles = NULL;
    nrecv = 0;

    //printf("CNUM1a %d: %d %d\n",me,cnum,nrecv);

    return;
  }

  // perform irregular comm to rendezvous procs
  // scan received buf with nrecv grid cells to fill in ctask fields

  char *buf;
  nrecv = comm->send_cells_adapt(nsend,procsend,(char *) sadapt,&buf);

  //printf("NSEND me %d: nsend %d nrecv %d \n",me,nsend,nrecv);

  sa_header = (SendAdapt **) 
    memory->smalloc(nrecv*sizeof(SendAdapt *),"adapt_grid::sa_header");
  sa_csurfs = (int **) 
    memory->smalloc(nrecv*sizeof(int *),"adapt_grid::sa_csurfs");
  sa_particles = (char **) 
    memory->smalloc(nrecv*sizeof(char *),"adapt_grid::sa_particles");
  int nbytes_total = sizeof(Particle::OnePart) + particle->sizeof_custom();

  int inum;

  char *ptr = buf;
  for (i = 0; i < nrecv; i++) {
    sa_header[i] = (SendAdapt *) ptr;
    ptr += sizeof(SendAdapt);
    ptr = ROUNDUP(ptr);

    sa_csurfs[i] = (int *) ptr;
    ptr += sa_header[i]->nsurf * sizeof(int);
    ptr = ROUNDUP(ptr);

    sa_particles[i] = (char *) ptr;
    ptr += sa_header[i]->np * nbytes_total;

    iparent = sa_header[i]->iparent;
    inum = parent2task[iparent];
    if (inum < 0) error->one(FLERR,"Parent cell has no coarsen task");

    m = sa_header[i]->ichild;
    ctask[inum].proc[m] = sa_header[i]->proc;
    ctask[inum].index[m] = sa_header[i]->icell;
    ctask[inum].recv[m] = i;
    ctask[inum].ncomplete++;

    /*
    printf("BUF %d %d inum %d: proc %d icell %d ichild %d iparentID %d\n",
           me,i,inum,
           sa_header[i]->proc,
           sa_header[i]->icell,
           sa_header[i]->ichild,
           pcells[sa_header[i]->iparent].id);
    */
  }

  // flag = 1 if my tasks now have all necessary child cell info
  // flag = 0 if they do not, b/c some other proc owns all child cells
  // error if ncompete is non-zero but not equal to nxyz

  // DEBUG

  /*
  for (i = 0; i < cnum; i++) {
    //if (!ctask[i].flag) continue;
    if (ctask[i].ncomplete == 4)
      printf("CAND %d: %d: parent %d: %d %d %d %d: %d %d %d %d\n",
             me,i,pcells[ctask[i].iparent].id,
             ctask[i].proc[0],
             ctask[i].proc[1],
             ctask[i].proc[2],
             ctask[i].proc[3],
             ctask[i].index[0],
             ctask[i].index[1],
             ctask[i].index[2],
             ctask[i].index[3]);
    else
      printf("CAND %d: %d: parent %d: ncomplete %d\n",
             me,i,pcells[ctask[i].iparent].id,ctask[i].ncomplete);
  }
  */

  for (i = 0; i < cnum; i++) {
    if (ctask[i].ncomplete == ctask[i].nchild) ctask[i].flag = 1;
    else if (ctask[i].ncomplete == 0) ctask[i].flag = 0;
    else error->one(FLERR,"Parent cell to coarsen has incomplete info");
  }

  // clean up

  //printf("DESTROY AAA %d\n",me);
  memory->destroy(parent2task);
  memory->destroy(procsend);

  //printf("CNUM1b %d: %d %d\n",me,cnum,nrecv);
}

/* ----------------------------------------------------------------------
   coarsen based on particle count
------------------------------------------------------------------------- */

void AdaptGrid::coarsen_particle()
{
  int m,iparent,nchild,np,icell,nsplit,jcell;
  int *proc,*index,*recv,*csubs;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;

  for (int i = 0; i < cnum; i++) {
    if (ctask[i].flag == 0) continue;
    ctask[i].flag = 0;

    iparent = ctask[i].iparent;
    nchild = ctask[i].nchild;
    proc = ctask[i].proc;
    index = ctask[i].index;
    recv = ctask[i].recv;

    np = 0;

    for (m = 0; m < nchild; m++) {
      if (proc[m] == me) {
        icell = index[m];
        if (cells[icell].nsplit == 1) np += cinfo[icell].count;
        else {
          nsplit = cells[icell].nsplit;
          csubs = sinfo[cells[icell].isplit].csubs;
          for (int i = 0; i < nsplit; i++) {
            jcell = csubs[i];
            np += cinfo[jcell].count;
          }
        }
      } else np += sa_header[recv[m]]->np;
    }

    if (np < ccount) ctask[i].flag = 1;
  }
}

/* ----------------------------------------------------------------------
   coarsen based on surfs in cell
   if any child cell smaller in any dimension than surfsize, coarsen to parent
   only use surfs with normals opposing vstream
------------------------------------------------------------------------- */

void AdaptGrid::coarsen_surf()
{
  int i,j,m,iparent,nchild,icell,nsurf,flag;
  int *proc,*index,*recv,*csurfs;
  double *lo,*hi,*norm;

  int dim = domain->dimension;
  Grid::ParentCell *pcells = grid->pcells;
  Grid::ChildCell *cells = grid->cells;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  for (i = 0; i < cnum; i++) {
    if (ctask[i].flag == 0) continue;
    ctask[i].flag = 0;

    iparent = ctask[i].iparent;
    nchild = ctask[i].nchild;
    proc = ctask[i].proc;
    index = ctask[i].index;
    recv = ctask[i].recv;

    lo = pcells[iparent].lo;
    hi = pcells[iparent].hi;
    flag = 0;
    if (fabs(hi[0]-lo[0])/pcells[iparent].nx < surfsize) flag = 1;
    if (fabs(hi[1]-lo[1])/pcells[iparent].ny < surfsize) flag = 1;
    if (dim == 3 && fabs(hi[2]-lo[2])/pcells[iparent].nz < surfsize) flag = 1;
    if (!flag) continue;

    flag = 0;
    for (m = 0; m < nchild; m++) {
      if (proc[m] == me) {
        icell = index[m];
        if (!cells[icell].nsurf) continue;
        nsurf = cells[icell].nsurf;
        csurfs = cells[icell].csurfs;
        for (j = 0; j < nsurf; j++) {
          if (dim == 2) {
            if (!(lines[j].mask & sgroupbit)) continue;
          } else {
            if (!(tris[j].mask & sgroupbit)) continue;
          }
          if (dim == 2) norm = lines[csurfs[j]].norm;
          else norm = tris[csurfs[j]].norm;
          if (MathExtra::dot3(norm,sdir) < 0.0) break;
        }
        if (j < nsurf) {
          flag = 1;
          break;
        }
      } else {
	nsurf = sa_header[recv[m]]->nsurf;
        if (!nsurf) continue;
	csurfs = sa_csurfs[recv[m]];
        for (j = 0; j < nsurf; j++) {
          if (dim == 2) {
            if (!(lines[j].mask & sgroupbit)) continue;
          } else {
            if (!(tris[j].mask & sgroupbit)) continue;
          }
          if (dim == 2) norm = lines[csurfs[j]].norm;
          else norm = tris[csurfs[j]].norm;
          if (MathExtra::dot3(norm,sdir) < 0.0) break;
        }
        if (j < nsurf) {
          flag = 1;
          break;
        }
      }
    }

    if (flag) ctask[i].flag = 1;
  }
}

/* ----------------------------------------------------------------------
   coarsen based on compute or fix value
   NOTE: insure fix ave/grid comes before fix adapt at end of step
------------------------------------------------------------------------- */

void AdaptGrid::coarsen_value()
{
  int m,iparent,nchild,icell,jcell,nsplit;
  int *proc,*index,*recv,*csubs;
  double value;

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  for (int i = 0; i < cnum; i++) {
    if (ctask[i].flag == 0) continue;
    ctask[i].flag = 0;

    iparent = ctask[i].iparent;
    nchild = ctask[i].nchild;
    proc = ctask[i].proc;
    index = ctask[i].index;
    recv = ctask[i].recv;

    if (combine == SUM) value = 0.0;
    else if (combine == MINIMUM) value = BIG;
    else if (combine == MAXIMUM) value = -BIG;

    for (m = 0; m < nchild; m++) {
      if (proc[m] == me) {
        icell = index[m];
        if (cells[icell].nsplit == 1) {
          if (combine == SUM) {
            if (valuewhich == COMPUTE) value += value_compute(icell);
            else if (valuewhich == FIX) value += value_fix(icell);
          } else if (combine == MINIMUM) {
            if (valuewhich == COMPUTE) 
              value = MIN(value,value_compute(icell));
            else if (valuewhich == FIX) 
              value = MIN(value,value_fix(icell));
          } else if (combine == MAXIMUM) {
            if (valuewhich == COMPUTE) 
              value = MAX(value,value_compute(icell));
            else if (valuewhich == FIX) 
              value = MAX(value,value_fix(icell));
          }
        } else {

          nsplit = cells[icell].nsplit;
          csubs = sinfo[cells[icell].isplit].csubs;
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
        }
      } else {
        if (combine == SUM) value += sa_header[recv[m]]->value;
        else if (combine == MINIMUM) 
          value = MIN(value,sa_header[recv[m]]->value);
        else if (combine == MAXIMUM) 
          value = MAX(value,sa_header[recv[m]]->value);
      }
    }
    
    //printf("BBB %d %d: %g %g\n",
    //      i,grid->pcells[iparent].id,value,cvalue);

    if (cdecide == LESS) {
      if (value < cvalue) ctask[i].flag = 1;
    } else {
      if (value > cvalue) ctask[i].flag = 1;
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
    if (ctask[i].flag == 0) continue;
    ctask[i].flag = 0;
    if (random->uniform() >= cfrac) continue;
    ctask[i].flag = 1;
  }
}

/* ----------------------------------------------------------------------
   perform coarsening of parent cells in ctask
------------------------------------------------------------------------- */

int AdaptGrid::perform_coarsen()
{
  int i,m,icell,iparent,nchild;
  int nsplit,jcell,ip;
  int *proc,*index,*recv,*csubs;

  Grid::ChildCell *cells = grid->cells;

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  // loop over coarsening tasks

  int *procreply,*cellreply;
  memory->create(procreply,nrecv,"adapt_grid:procreply");
  memory->create(cellreply,nrecv,"adapt_grid:cellreply");
  int nreply = 0;

  int removeparent = 0;
  memory->create(delparent,cnum,"adapt_grid:delparent");

  // DEBUG

  /*
  for (i = 0; i < cnum; i++) {
    if (!ctask[i].flag) continue;
    printf("CT %d: %d: parent %d\n",me,i,pcells[ctask[i].iparent].id);
  }
  */

  for (i = 0; i < cnum; i++) {

    // skip if not selected
    // no need to notify owners of children of this parent

    if (!ctask[i].flag) continue;

    iparent = ctask[i].iparent;
    nchild = ctask[i].nchild;
    proc = ctask[i].proc;
    index = ctask[i].index;
    recv = ctask[i].recv;
    delparent[removeparent++] = iparent;

    // flag children for deletion
    // for mine, set proc = -1 and remove from grid hash
    // for others, add to reply list

    for (m = 0; m < nchild; m++) {
      //if (me == 0 && pcells[iparent].id == 3)
      //  printf("DELETING CHILD m %d index %d proc %d cellID %d\n",
      //         m,index[m],proc[m],cells[index[m]].id);
      if (proc[m] == me) {
	icell = index[m];
	cells[icell].proc = -1;
	hash->erase(cells[icell].id);
      } else {
        procreply[nreply] = proc[m];
        cellreply[nreply] = index[m];
        nreply++;
      }
    }

    // coarsen iparent to a new child cell
    // coarsen_cell() will add new child to grid hash
    // also add it to chash to prevent it from being refined

    grid->coarsen_cell(iparent,nchild,proc,index,recv,this,cut2d,cut3d);
    cells = grid->cells;
    (*chash)[ctask[i].id] = 0;

    // add new child cell to newcells list

    if (nnew == maxnew) {
      maxnew += DELTA_NEW;
      memory->grow(newcells,maxnew,"adapt_grid:newcells");
    }
    newcells[nnew++] = ctask[i].id;
  }

  // return if no one has any child cell info to reply to another proc

  MPI_Allreduce(&nreply,&replyany,1,MPI_INT,MPI_SUM,world);
  if (!replyany) {
    memory->sfree(sa_header);
    memory->sfree(sa_csurfs);
    memory->sfree(sa_particles);
    memory->destroy(procreply);
    memory->destroy(cellreply);
    nrecv = 0;

    //printf("CNUM3 %d: %d %d\n",me,removeparent,nrecv);

    return removeparent;
  }

  // perform irregular comm from rendezvous procs
  // recv list of my cells coarsened into new child cell on another proc
  // scan received buf to mark my grid cells for deletion 
  //   and remove from grid hash
  // don't need to flag sub cells, since grid compress will do that
  // also mark particles of cell for deletion
  // if cell is split, mark particles in all its sub cells for deletion

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;

  char *buf;
  nrecv = comm->irregular_uniform(nreply,procreply,
                                  (char *) cellreply,sizeof(int),&buf);
  int *ibuf = (int *) buf;

  for (i = 0; i < nrecv; i++) {
    icell = ibuf[i];
    cells[icell].proc = -1;
    hash->erase(cells[icell].id);

    if (cells[icell].nsplit == 1) {
      ip = cinfo[icell].first;
      while (ip >= 0) {
        particles[ip].icell = -1;
        ip = next[ip];
      }
    } else {
      nsplit = cells[icell].nsplit;
      csubs = sinfo[cells[icell].isplit].csubs;
      for (m = 0; m < nsplit; m++) {
        jcell = csubs[m];
        ip = cinfo[jcell].first;
        while (ip >= 0) {
          particles[ip].icell = -1;
          ip = next[ip];
        }
      }
    }
  }

  // clean up

  memory->sfree(sa_header);
  memory->sfree(sa_csurfs);
  memory->sfree(sa_particles);
  memory->destroy(procreply);
  memory->destroy(cellreply);

  // return # of parents removed by this proc

  //printf("CNUM3 %d: %d %d\n",me,removeparent,nrecv);

  return removeparent;
}

/* ----------------------------------------------------------------------
   gather deleted parent cells from all procs, on each proc
   must end up with same list, in same order, on every proc
   must repoint child cells to new parent
   also point child cells that became parent to new parent
     this is so can reassign particles to correct new child cells
------------------------------------------------------------------------- */

void AdaptGrid::gather_parents_coarsen(int delta, int nrefine)
{
  int i,m,iparent,nsplit;
  cellint tmp;
  int *csubs;

#ifdef SPARTA_MAP
  std::map<cellint,int> *hash = grid->hash;
#elif defined SPARTA_UNORDERED_MAP
  std::unordered_map<cellint,int> *hash = grid->hash;
#else
  std::tr1::unordered_map<cellint,int> *hash = grid->hash;
#endif

  // perform Allgatherv of removed parent cells from all procs

  int *recvcounts,*displs;
  memory->create(recvcounts,nprocs,"adapt_grid:recvcounts");
  memory->create(displs,nprocs,"adapt_grid:displs");

  MPI_Allgather(&delta,1,MPI_INT,recvcounts,1,MPI_INT,world);
  displs[0] = 0;
  for (i = 1; i < nprocs; i++) displs[i] = displs[i-1] + recvcounts[i-1];

  int deltaall;
  MPI_Allreduce(&delta,&deltaall,1,MPI_INT,MPI_SUM,world);
  int *delparentall;
  memory->create(delparentall,deltaall,"adapt_grid:delparentall");

  MPI_Allgatherv(delparent,delta,MPI_INT,
                 delparentall,recvcounts,displs,MPI_INT,world);

  // flag all parent cells to remove with id = 0
  // remove parent cell from hash

  Grid::ParentCell *pcells = grid->pcells;

  for (i = 0; i < deltaall; i++) {
    iparent = delparentall[i];
    hash->erase(pcells[iparent].id);
    pcells[iparent].id = 0;
  }

  // remove all flagged parent cells
  // shift cells that are kept downward to preserve original ordering
  // reset hash to new location in pcells

  int psize = sizeof(Grid::ParentCell);

  int ncurrent = grid->nparent;
  int nparent = 1;
  for (i = 1; i < ncurrent; i++) {
    if (pcells[i].id) {
      if (i > nparent) memcpy(&pcells[nparent],&pcells[i],psize);
      (*hash)[pcells[i].id] = -(nparent+1);
      nparent++;
    }
  }

  grid->nparent = nparent;

  // loop over all parent cells
  // reset their iparent field and grandparent flag

  for (i = 1; i < nparent; i++) {
    iparent = grid->id_find_parent(pcells[i].id,tmp);
    pcells[i].iparent = iparent;
    pcells[i].grandparent = 0;
    pcells[pcells[i].iparent].grandparent = 1;
  }

  // DEBUG

  /*
  char str[32];
  printf("POST COARSEN PARENT %d: %d\n",comm->me,grid->nparent);
  for (int i = 0; i < grid->nparent; i++) {
    Grid::ParentCell *p = &grid->pcells[i];
    grid->id_num2str(p->id,str);
    printf("PCELL %d: %d id %d %s iparent %d %d level %d "
           "nxyz %d %d %d lo %g %g %g "
           "hi %g %g %g\n",
           comm->me,i,p->id,str,p->iparent,p->grandparent,p->level,
           p->nx,p->ny,p->nz,
           p->lo[0],
           p->lo[1],
           p->lo[2],
           p->hi[0],
           p->hi[1],
           p->hi[2]);
  }


  printf("POST COARSEN CHILD %d: %d\n",comm->me,grid->nlocal);
  for (int i = 0; i < grid->nlocal; i++) {
    Grid::ChildCell *g = &grid->cells[i];
    if (g->nsplit <= 1) {
      grid->id_num2str(g->id,str);
      printf("ICELL %d: %d id %d iparent %d proc %d\n",
             comm->me,i,g->id,g->iparent,g->proc);
    } else
      printf("ICELLSPLIT %d id %d iparent %d proc %d nsplit %d isplit %d "
             "spicell %d csubs %d %d lo %g %g "
             "hi %g %g inout %d\n",
             i,g->id,g->iparent,g->proc,g->nsplit,g->isplit,
             grid->sinfo[g->isplit].icell,
             grid->sinfo[g->isplit].csubs[0],grid->sinfo[g->isplit].csubs[1],
             g->lo[0],
             g->lo[1],
             g->hi[0],
             g->hi[1],
             grid->cinfo[i].type);
  }
  */

  // loop over all child cells and reset their iparent field
  // skip cells with proc < 0 that are marked for deletion
  // skip sub cells since parent no longer exists for ones about to be deleted
  // set non-deleted sub cells via parent split cell

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  int nglocal = grid->nlocal;

  for (i = 0; i < nglocal; i++) {
    if (cells[i].proc < 0) continue;
    if (cells[i].nsplit <= 0) continue;

    //printf("FIND PARENT %d %d: %d iparent %d\n",me,i,cells[i].id,
    //       cells[i].iparent);
    iparent = grid->id_find_parent(cells[i].id,tmp);
    //printf("  found PARENT %d %d: %d %d\n",me,i,cells[i].id,iparent);
    cells[i].iparent = iparent;
    if (cells[i].nsplit > 1) {
      nsplit = cells[i].nsplit;
      csubs = sinfo[cells[i].isplit].csubs;
      for (m = 0; m < nsplit; m++) cells[csubs[m]].iparent = iparent;
    }
  }

  // clean up

  memory->destroy(recvcounts);
  memory->destroy(displs);
  memory->sfree(delparentall);
}


/* ----------------------------------------------------------------------
   extract a value for icell,valindex from a compute
------------------------------------------------------------------------- */

double AdaptGrid::value_compute(int icell)
{
  double value;

  if (valindex == 0) {
    if (compute->post_process_grid_flag) 
      compute->post_process_grid(NULL,NULL,icell,0,&value,1);
    else value = compute->vector_grid[icell];
    
  } else {
    if (compute->post_process_grid_flag)
      compute->post_process_grid(NULL,NULL,icell,valindex,&value,1);
    else value = compute->array_grid[icell][valindex-1];
  }

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

  if (domain->dimension == 3) delete cut3d;
  else delete cut2d;

  memory->destroy(rlist);
  if (ctask) {
    for (int i = 0; i < cnum; i++) {
      delete [] ctask[i].proc;
      delete [] ctask[i].index;
      delete [] ctask[i].recv;
    }
    memory->sfree(ctask);
  }
  memory->sfree(sadapt);
  memory->sfree(delparent);

  delete chash;
}

/* ----------------------------------------------------------------------
   write out new parent grid file via WriteGrid
------------------------------------------------------------------------- */

void AdaptGrid::write_file()
{
  WriteGrid *wg = new WriteGrid(sparta);
  wg->silent = 1;

  int narg = 2;
  char **args = new char*[narg];
  args[0] = (char *) "parent";

  char *expandfile = NULL;
  if (strchr(file,'*')) {
    expandfile = new char[strlen(file) + 16];
    char *ptr = strchr(file,'*');
    *ptr = '\0';
    sprintf(expandfile,"%s" BIGINT_FORMAT "%s",file,update->ntimestep,ptr+1);
    *ptr = '*';
    args[1] = expandfile;
  } else args[1] = file;

  wg->command(narg,args);

  // NOTE: could persist wg for fix adapt

  delete [] expandfile;
  delete [] args;
  delete wg;
}
