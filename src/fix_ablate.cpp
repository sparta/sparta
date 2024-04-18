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

#include "spatype.h"
#include "stdlib.h"
#include "string.h"
#include "fix_ablate.h"
#include "update.h"
#include "grid.h"
#include "domain.h"
#include "decrement_lookup_table.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "output.h"
#include "input.h"
#include "variable.h"
#include "dump.h"
#include "marching_squares.h"
#include "marching_cubes.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{COMPUTE,FIX,VARIABLE,RANDOM};
enum{CVALUE,CDELTA};

#define INVOKED_PER_GRID 16
#define DELTAGRID 1024            // must be bigger than split cells per cell
#define DELTASEND 1024
#define EPSILON 1.0e-4            // this is on a scale of 0 to 255

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Update

// remove if fix particles-inside-surfs issue
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};  // several files

// NOTES
// should I store one value or 8 per cell
// how to preserve svalues and set type of new surfs,
//   maybe need local per-cell type
// after create new surfs, need to assign group, sc, sr
// how to impose sgroup like ReadIsurf does after surfs created
// how to prevent adaptation of any cells in the implicit grid group(s)?
//   just testing for surfs is not enough, since it may have no surfs
// need to update neigh corner points even if neigh cell is not in group?
// do a run-time test for gridcut = 0.0 ?
// worry about array_grid having updated values for sub cells?  line in store()
// can I output both per-grid and global values?

/* ---------------------------------------------------------------------- */

FixAblate::FixAblate(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  MPI_Comm_rank(world,&me);

  if (narg < 6) error->all(FLERR,"Illegal fix ablate command");

  igroup = grid->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Could not find fix ablate group ID");
  groupbit = grid->bitmask[igroup];

  nevery = atoi(arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix ablate command");

  idsource = NULL;

  scale = atof(arg[4]);
  if (scale < 0.0) error->all(FLERR,"Illegal fix ablate command");

  int iarg = 6;
  if ((strncmp(arg[5],"c_",2) == 0) || (strncmp(arg[5],"f_",2) == 0)) {
    if (arg[5][0] == 'c') which = COMPUTE;
    else if (arg[5][0] == 'f') which = FIX;

    int n = strlen(arg[5]);
    char *suffix = new char[n];
    strcpy(suffix,&arg[5][2]);

    char *ptr = strchr(suffix,'[');
    if (ptr) {
      if (suffix[strlen(suffix)-1] != ']')
        error->all(FLERR,"Illegal fix ablate command");
      argindex = atoi(ptr+1);
      *ptr = '\0';
    } else argindex = 0;

    n = strlen(suffix) + 1;
    idsource = new char[n];
    strcpy(idsource,suffix);
    delete [] suffix;

  } else if (strncmp(arg[5],"v_",2) == 0) {
    which = VARIABLE;

    int n = strlen(arg[5]);
    char *idsource = new char[n];
    strcpy(idsource,&arg[5][2]);

  } else if (strcmp(arg[5],"random") == 0) {
    iarg++; // one additional input
    if (narg < 7) error->all(FLERR,"Illegal fix ablate command");
    which = RANDOM;
    maxrandom = atoi(arg[6]);

  } else error->all(FLERR,"Illegal fix ablate command");

  process_args(narg-iarg,&arg[iarg]);

  // error check

  if (which == COMPUTE) {
    icompute = modify->find_compute(idsource);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for fix ablate does not exist");
    if (modify->compute[icompute]->per_grid_flag == 0)
      error->all(FLERR,
                 "Fix ablate compute does not calculate per-grid values");
    if (modify->compute[icompute]->post_process_isurf_grid_flag == 0)
      error->all(FLERR,
                 "Fix ablate compute does not calculate isurf per-grid values");
    if (argindex == 0 &&
        modify->compute[icompute]->size_per_grid_cols != 0)
      error->all(FLERR,"Fix ablate compute does not "
                 "calculate per-grid vector");
    if (argindex && modify->compute[icompute]->size_per_grid_cols == 0)
      error->all(FLERR,"Fix ablate compute does not "
                 "calculate per-grid array");
    if (argindex && argindex > modify->compute[icompute]->size_per_grid_cols)
      error->all(FLERR,"Fix ablate compute array is accessed out-of-range");

  } else if (which == FIX) {
    ifix = modify->find_fix(idsource);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for fix ablate does not exist");
    if (modify->fix[ifix]->per_grid_flag == 0)
      error->all(FLERR,"Fix ablate fix does not calculate per-grid values");
    if (argindex == 0 && modify->fix[ifix]->size_per_grid_cols != 0)
      error->all(FLERR,
                 "Fix ablate fix does not calculate per-grid vector");
    if (argindex && modify->fix[ifix]->size_per_grid_cols == 0)
      error->all(FLERR,
                 "Fix ablate fix does not calculate per-grid array");
    if (argindex && argindex > modify->fix[ifix]->size_per_grid_cols)
      error->all(FLERR,"Fix ablate fix array is accessed out-of-range");
    if (nevery % modify->fix[ifix]->per_grid_freq)
      error->all(FLERR,
                 "Fix for fix ablate not computed at compatible time");

  } else if (which == VARIABLE) {
    ivariable = input->variable->find(idsource);
    if (ivariable < 0)
      error->all(FLERR,"Could not find fix ablate variable name");
    if (input->variable->grid_style(ivariable) == 0)
      error->all(FLERR,"Fix ablate variable is not grid-style variable");
  }

  // this fix produces a per-grid array and a scalar

  dim = domain->dimension;

  per_grid_flag = 1;
  if (dim == 2) size_per_grid_cols = 4;
  else size_per_grid_cols = 8;
  per_grid_freq = 1;
  gridmigrate = 1;

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  sum_delta = 0.0;
  ndelete = 0;

  storeflag = adjacentflag = 0;
  array_grid = cvalues = NULL;
  cavalues = NULL;
  tvalues = NULL;
  ncorner = size_per_grid_cols;

  if(dim == 2) nadj = 4;
  else nadj = 6;

  // local storage

  ixyz = NULL;
  mcflags = NULL;
  celldelta = NULL;
  cdelta = NULL;
  cdelta_ghost = NULL;
  cadelta = NULL;
  cadelta_ghost = NULL;
  numsend = NULL;
  maxgrid = maxghost = 0;

  proclist = NULL;
  locallist = NULL;
  maxsend = 0;

  sbuf = NULL;
  maxbuf = 0;

  vbuf = NULL;
  maxvar = 0;

  ms = NULL;
  mc = NULL;

  // RNG for random decrements
  // for now, use same RNG on every proc
  // uncomment two lines if want to change that
  // b/c set_delta_random() is decrementing the same no matter who owns a cell

  random = NULL;
  if (which == RANDOM) {
    random = new RanKnuth(update->ranmaster->uniform());
    //double seed = update->ranmaster->uniform();
    //random->reset(seed,comm->me,100);
  }

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  if (nevery) {
    bigint nvalid = (update->ntimestep/nevery)*nevery + nevery;
    modify->addstep_compute_all(nvalid);
  }
}



/* ---------------------------------------------------------------------- */

FixAblate::~FixAblate()
{
  delete [] idsource;
  memory->destroy(cvalues);
  memory->destroy(cavalues);
  memory->destroy(tvalues);

  memory->destroy(ixyz);
  memory->destroy(mcflags);
  memory->destroy(celldelta);
  memory->destroy(cdelta);
  memory->destroy(cdelta_ghost);
  memory->destroy(cadelta);
  memory->destroy(cadelta_ghost);
  memory->destroy(numsend);

  memory->destroy(proclist);
  memory->destroy(locallist);

  memory->destroy(sbuf);
  memory->destroy(vbuf);

  delete ms;
  delete mc;

  delete random;
}

/* ---------------------------------------------------------------------- */

int FixAblate::setmask()
{
  int mask = 0;
  if (nevery) mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   store grid corner point and type values in cvalues and tvalues
   then create implicit surfaces
   called by ReadIsurf when corner point grid is read in
------------------------------------------------------------------------- */
void FixAblate::store_corners(int nx_caller, int ny_caller, int nz_caller,
                              double *cornerlo_caller, double *xyzsize_caller,
                              double **cvalues_caller, int *tvalues_caller,
                              double thresh_caller, char *sgroupID, int pushflag)
{
  storeflag = 1;

  nx = nx_caller;
  ny = ny_caller;
  nz = nz_caller;
  cornerlo[0] = cornerlo_caller[0];
  cornerlo[1] = cornerlo_caller[1];
  cornerlo[2] = cornerlo_caller[2];
  xyzsize[0] = xyzsize_caller[0];
  xyzsize[1] = xyzsize_caller[1];
  xyzsize[2] = xyzsize_caller[2];
  thresh = thresh_caller;

  tvalues_flag = 0;
  if (tvalues_caller) tvalues_flag = 1;

  if (sgroupID) {
    int sgroup = surf->find_group(sgroupID);
    if (sgroup < 0) sgroup = surf->add_group(sgroupID);
    sgroupbit = surf->bitmask[sgroup];
  } else sgroupbit = 0;

  // allocate per-grid cell data storage

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  nglocal = grid->nlocal;

  grow_percell(0);

  // copy caller values into local values of FixAblate

  for (int icell = 0; icell < nglocal; icell++) {
    for (int m = 0; m < ncorner; m++)
      cvalues[icell][m] = cvalues_caller[icell][m];
    if (tvalues_flag) tvalues[icell] = tvalues_caller[icell];
  }

  // set ix,iy,iz indices from 1 to Nxyz for each of my owned grid cells
  // same logic as ReadIsurf::create_hash()

  for (int i = 0; i < nglocal; i++)
    ixyz[i][0] = ixyz[i][1] = ixyz[i][2] = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ixyz[icell][0] =
      static_cast<int> ((cells[icell].lo[0]-cornerlo[0]) / xyzsize[0] + 0.5) + 1;
    ixyz[icell][1] =
      static_cast<int> ((cells[icell].lo[1]-cornerlo[1]) / xyzsize[1] + 0.5) + 1;
    ixyz[icell][2] =
      static_cast<int> ((cells[icell].lo[2]-cornerlo[2]) / xyzsize[2] + 0.5) + 1;
  }

  // push corner pt values that are fully external/internal to 0 or 255

  if (pushflag) push_lohi();

  if (adjacentflag) epsilon_adjust_adjacent();
  else epsilon_adjust();

  // create marching squares/cubes classes, now that have group & threshold

  if (dim == 2) ms = new MarchingSquares(sparta,igroup,thresh);
  else mc = new MarchingCubes(sparta,igroup,thresh);

  // create implicit surfaces

  create_surfs(1);
}

/* ----------------------------------------------------------------------
   store grid corner point and type values in cvalues and tvalues
   then create implicit surfaces
   called by ReadIsurf when corner point grid is read in
------------------------------------------------------------------------- */

void FixAblate::store_corners(int nx_caller, int ny_caller, int nz_caller,
                              double *cornerlo_caller, double *xyzsize_caller,
                              double ***cavalues_caller, int *tvalues_caller,
                              double thresh_caller, char *sgroupID, int pushflag)
{
  storeflag = 1;
  adjacentflag = 1;

  nx = nx_caller;
  ny = ny_caller;
  nz = nz_caller;
  cornerlo[0] = cornerlo_caller[0];
  cornerlo[1] = cornerlo_caller[1];
  cornerlo[2] = cornerlo_caller[2];
  xyzsize[0] = xyzsize_caller[0];
  xyzsize[1] = xyzsize_caller[1];
  xyzsize[2] = xyzsize_caller[2];
  thresh = thresh_caller;

  tvalues_flag = 0;
  if (tvalues_caller) tvalues_flag = 1;

  if (sgroupID) {
    int sgroup = surf->find_group(sgroupID);
    if (sgroup < 0) sgroup = surf->add_group(sgroupID);
    sgroupbit = surf->bitmask[sgroup];
  } else sgroupbit = 0;

  // allocate per-grid cell data storage

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  nglocal = grid->nlocal;

  grow_percell(0);

  // copy caller values into local values of FixAblate

  for (int icell = 0; icell < nglocal; icell++) {
    for (int m = 0; m < ncorner; m++)
      for(int n = 0; n < nadj; n++)
        cavalues[icell][m][n] = cavalues_caller[icell][m][n];
    if (tvalues_flag) tvalues[icell] = tvalues_caller[icell];
  }

  // set ix,iy,iz indices from 1 to Nxyz for each of my owned grid cells
  // same logic as ReadIsurf::create_hash()

  for (int i = 0; i < nglocal; i++)
    ixyz[i][0] = ixyz[i][1] = ixyz[i][2] = 0;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ixyz[icell][0] =
      static_cast<int> ((cells[icell].lo[0]-cornerlo[0]) / xyzsize[0] + 0.5) + 1;
    ixyz[icell][1] =
      static_cast<int> ((cells[icell].lo[1]-cornerlo[1]) / xyzsize[1] + 0.5) + 1;
    ixyz[icell][2] =
      static_cast<int> ((cells[icell].lo[2]-cornerlo[2]) / xyzsize[2] + 0.5) + 1;
  }

  if (adjacentflag) epsilon_adjust_adjacent();
  else epsilon_adjust();

  // create marching squares/cubes classes, now that have group & threshold

  if (dim == 2) ms = new MarchingSquares(sparta,igroup,thresh);
  else mc = new MarchingCubes(sparta,igroup,thresh);

  // create implicit surfaces

  create_surfs(1);
}

/* ---------------------------------------------------------------------- */

void FixAblate::init()
{
  if (!storeflag)
    error->all(FLERR,"Fix ablate corner point values not stored");

  if (which == COMPUTE) {
    icompute = modify->find_compute(idsource);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for fix ablate does not exist");
  } else if (which == FIX) {
    ifix = modify->find_fix(idsource);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for fix ablate does not exist");
  } else if (which == VARIABLE) {
    ivariable = input->variable->find(idsource);
    if (ivariable < 0)
      error->all(FLERR,"Variable ID for fix ablate does not exist");
  }

  // reallocate per-grid data if necessary

  nglocal = grid->nlocal;
  grow_percell(0);
}

/* ---------------------------------------------------------------------- */

void FixAblate::end_of_step()
{
  // set per-cell delta vector randomly or from compute/fix source

  //if (which == RANDOM) set_delta_random();
  if (which == RANDOM) set_delta_uniform();
  else set_delta();

  // find decrement at each corner -> sync decrements -> update corners

  if (distributedflag) {
    decrement_distributed_outside();
    sync_distributed_outside();
    if (adjacentflag) {
      if(dim == 2) decrement_adjacent_distributed_inside2d();
      else decrement_adjacent_distributed_inside3d();
      sync_adjacent_distributed_inside();
    } else {
      decrement_distributed_inside();
      sync_distributed_inside();
    }
  } else {
    if (adjacentflag) {
      decrement_adjacent();
      sync_adjacent();
    } else {
      decrement();
      sync();
    }
  }

  // handle small corner values

  if (adjacentflag) epsilon_adjust_adjacent();
  else epsilon_adjust();

  // re-create implicit surfs

  create_surfs(0);
}

/* ---------------------------------------------------------------------- */

void FixAblate::create_surfs(int outflag)
{
  // DEBUG
  // store copy of last ablation's per-cell MC flags before a new ablation

  int **mcflags_old = mcflags;
  memory->create(mcflags,maxgrid,4,"ablate:mcflags");
  for (int i = 0; i < maxgrid; i++)
    mcflags[i][0] = mcflags[i][1] = mcflags[i][2] = mcflags[i][3] = -1;

  // sort existing particles since may be clearing split cells

  if (!particle->sorted) particle->sort();

  // reassign particles in sub cells to all be in parent split cell

  if (grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->combine_split_cell_particles(icell,1);
  }

  // call clear_surf before create new surfs, so cell/corner flags are all set

  grid->unset_neighbors();
  grid->remove_ghosts();
  grid->clear_surf();
  surf->clear_implicit();

  // perform Marching Squares/Cubes to create new implicit surfs
  // cvalues = corner point values
  // tvalues = surf type for surfs in each grid cell

  if (adjacentflag) {
    if (dim == 2) ms->invoke(cavalues,tvalues);
    else mc->invoke(cavalues,tvalues,mcflags);
  } else {
    if (dim == 2) ms->invoke(cvalues,tvalues);
    else mc->invoke(cvalues,tvalues,mcflags);
  }

  // set surf->nsurf and surf->nown

  surf->nown = surf->nlocal;
  bigint nlocal = surf->nlocal;
  MPI_Allreduce(&nlocal,&surf->nsurf,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  // output extent of implicit surfs, some may be tiny

  if (outflag) {
    if (dim == 2) surf->output_extent(0);
    else surf->output_extent(0);
  }

  // compute normals of new surfs

  if (dim == 2) surf->compute_line_normal(0);
  else surf->compute_tri_normal(0);

  // MC->cleanup() checks for consistent triangles on grid cell faces
  // needs to come after normals are computed
  // it requires neighbor indices and ghost cell info
  // so first acquire ghosts (which will also grab surfs),
  //   then remove ghost surfs and ghost grid cells again

  if (dim == 3) {
    grid->acquire_ghosts(0);
    grid->reset_neighbors();
    mc->cleanup();
    surf->remove_ghosts();
    grid->unset_neighbors();
    grid->remove_ghosts();
  }

  // assign optional surf group to masks of new surfs

  if (sgroupbit) {
    int nsurf = surf->nlocal;
    if (dim == 3) {
      Surf::Tri *tris = surf->tris;
      for (int i = 0; i < nsurf; i++) tris[i].mask |= sgroupbit;
    } else {
      Surf::Line *lines = surf->lines;
      for (int i = 0; i < nsurf; i++) lines[i].mask |= sgroupbit;
    }
  }

  // assign surf collision/reaction models to newly created surfs
  // this assignment can be made in input script via surf_modify
  //   after implicit surfs are created
  // for active ablation, must be re-assigned at every ablation atep
  // for now just assume all surfs are assigned to first collide/react model
  // NOTE: need a more flexible way to do this

  int nslocal = surf->nlocal;

  if (dim == 2) {
    Surf::Line *lines = surf->lines;
    if (surf->nsc)
      for (int i = 0; i < nslocal; i++)
        lines[i].isc = 0;
    if (surf->nsr)
      for (int i = 0; i < nslocal; i++)
        lines[i].isr = 0;
  } else {
    Surf::Tri *tris = surf->tris;
    if (surf->nsc)
      for (int i = 0; i < nslocal; i++)
        tris[i].isc = 0;
    if (surf->nsr)
      for (int i = 0; i < nslocal; i++)
        tris[i].isr = 0;
  }

  // watertight check can be done before surfs are mapped to grid cells

  if (dim == 2) surf->check_watertight_2d();
  else surf->check_watertight_3d();

  // if no surfs created, use clear_surf to set all celltypes = OUTSIDE

  if (surf->nsurf == 0) {
    surf->exist = 0;
    grid->clear_surf();
  }

  // -----------------------
  // map surfs to grid cells
  // -----------------------

  // surfs are already assigned to grid cells
  // create split cells due to new surfs

  grid->surf2grid_implicit(1,outflag);

  // re-setup grid ghosts and neighbors

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  // flag cells and corners as OUTSIDE or INSIDE

  grid->set_inout();
  grid->type_check(outflag);

  // reassign particles in a split cell to sub cell owner
  // particles are unsorted afterwards, within new sub cells

  if (grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->assign_split_cell_particles(icell);
    particle->sorted = 0;
  }

  // notify all classes that store per-grid data that grid has changed

  grid->notify_changed();

  // ------------------------------------------------------------------------
  // DEBUG - should not have to do any of this once marching cubes is perfect
  // only necessary for 3d

  if (dim == 2) {
    memory->destroy(mcflags_old);
    return;
  }

  // DEBUG - if this line is uncommented, code will do delete no particles
  //         eventually this should work

  // if (dim == 3) {
  //   memory->destroy(mcflags_old);
  //   return;
  // }

  // DEBUG - remove all particles
  // if these lines are uncommented, all particles are wiped out

  // particle->nlocal = 0;
  // memory->destroy(mcflags_old);
  // return;

  // DEBUG - remove only the particles that are inside the surfs
  //         after ablation
  // similar code as in fix grid/check

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  Particle::OnePart *particles = particle->particles;
  int pnlocal = particle->nlocal;

  int ncount;
  int icell,splitcell,subcell,pflag;
  double *x;
  double xcell[3];

  ncount = 0;
  for (int i = 0; i < pnlocal; i++) {
    particles[i].flag = PKEEP;
    icell = particles[i].icell;
    if (cells[icell].nsurf == 0) continue;

    x = particles[i].x;

    // check that particle is outside surfs
    // if no xcell found, cannot check

    pflag = grid->point_outside_surfs(icell,xcell);
    if (!pflag) continue;
    pflag = grid->outside_surfs(icell,x,xcell);

    // check that particle is in correct split subcell

    if (pflag && cells[icell].nsplit <= 0) {
      splitcell = sinfo[cells[icell].isplit].icell;
      if (dim == 2) subcell = update->split2d(splitcell,x);
      else subcell = update->split3d(splitcell,x);
      if (subcell != icell) pflag = 0;
    }

    // discard the particle if either test failed

    if (!pflag) {
      particles[i].flag = PDISCARD;
      // DEBUG - print message about MC flags for cell of deleted particle
      /*
      printf("INSIDE PART: me %d id %d coords %g %g %g "
             "cellID %d celltype %d nsplit %d MCflags old %d %d %d %d "
             "MCflags now %d %d %d %d corners %g %g %g %g %g %g %g %g\n",
             comm->me,particles[i].id,x[0],x[1],x[2],
             cells[mcell].id,cinfo[mcell].type,cells[mcell].nsplit,
             mcflags_old[mcell][0],
             mcflags_old[mcell][1],
             mcflags_old[mcell][2],
             mcflags_old[mcell][3],
             mcflags[mcell][0],
             mcflags[mcell][1],
             mcflags[mcell][2],
             mcflags[mcell][3],
             cvalues[mcell][0],
             cvalues[mcell][1],
             cvalues[mcell][2],
             cvalues[mcell][3],
             cvalues[mcell][4],
             cvalues[mcell][5],
             cvalues[mcell][6],
             cvalues[mcell][7]);
      */
      ncount++;
    }
  }

  memory->destroy(mcflags_old);

  // compress out the deleted particles
  // NOTE: if end up keeping this section, need logic for custom particle vectors
  //       see Particle::compress_rebalance()

  int nbytes = sizeof(Particle::OnePart);

  int i = 0;
  while (i < pnlocal) {
    if (particles[i].flag == PDISCARD) {
      memcpy(&particles[i],&particles[pnlocal-1],nbytes);
      pnlocal--;
    } else i++;
  }

  MPI_Allreduce(&ncount,&ndelete,1,MPI_INT,MPI_SUM,world);

  particle->nlocal = pnlocal;
  particle->sorted = 0;
}

/* ----------------------------------------------------------------------
   set per-cell delta vector randomly
   celldelta = random integer between 0 and maxrandom
   scale = fraction of cells that are decremented
------------------------------------------------------------------------- */

void FixAblate::set_delta_random()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // enforce same decrement no matter who owns which cells
  // NOTE: could change this at some point, use differnet RNG for each proc

  if (!grid->hashfilled) grid->rehash();
  Grid::MyHash *hash = grid->hash;
  cellint cellID;
  int rn2,icell;
  double rn1;
  for (int i = 0; i < grid->ncell; i++) {
    rn1 = random->uniform();
    rn2 = static_cast<int> (random->uniform()*maxrandom) + 1.0;
    cellID = i+1;
    if (hash->find(cellID) == hash->end()) continue;
    icell = (*hash)[cellID];
    if (icell >= nglocal) continue;     // ghost cell

    if (rn1 > scale) celldelta[icell] = 0.0;
    else celldelta[icell] = rn2;
  }

  // total decrement for output

  double sum = 0.0;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;
    sum += celldelta[icell];
  }

  MPI_Allreduce(&sum,&sum_delta,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   set per-cell delta vector randomly
   celldelta = random integer between 0 and maxrandom
   scale = fraction of cells that are decremented
------------------------------------------------------------------------- */

void FixAblate::set_delta_uniform()
{
  int nin;
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // enforce same decrement no matter who owns which cells
  // NOTE: could change this at some point, use differnet RNG for each proc

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    nin = 0;
    for (int i = 0; i < ncorner; i++) {
      if (adjacentflag) {
        if (cavalues[icell][i][0] > thresh) nin++;
      } else {
        if (cvalues[icell][i] > thresh) nin++;
      }
    }
    if (nin == 0 || nin == ncorner) celldelta[icell] = 0.0;
    else celldelta[icell] = maxrandom*scale;
  }

  // total decrement for output

  double sum = 0.0;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;
    sum += celldelta[icell];
  }

  MPI_Allreduce(&sum,&sum_delta,1,MPI_DOUBLE,MPI_SUM,world);
}


/* ----------------------------------------------------------------------
   set per-cell delta vector from compute/fix/variable source
   celldelta = nevery * scale * source-value
   NOTE: how does this work for split cells? should only do parent split?
------------------------------------------------------------------------- */

void FixAblate::set_delta()
{
  int i;

  double prefactor = nevery*scale;

  // compute/fix may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  if (which == COMPUTE) {
    Compute *c = modify->compute[icompute];

    if (!(c->invoked_flag & INVOKED_PER_GRID)) {
      c->compute_per_grid();
      c->invoked_flag |= INVOKED_PER_GRID;
    }
    c->post_process_isurf_grid();

    if (argindex == 0) {
      double *cvec = c->vector_grid;
      for (i = 0; i < nglocal; i++)
        celldelta[i] = prefactor * cvec[i];
    } else {
      double **carray = c->array_grid;
      int im1 = argindex - 1;
      for (i = 0; i < nglocal; i++)
        celldelta[i] = prefactor * carray[i][im1];
    }

  } else if (which == FIX) {
    Fix *f = modify->fix[ifix];

    if (argindex == 0) {
      double *fvec = f->vector_grid;
      for (i = 0; i < nglocal; i++) {
        celldelta[i] = prefactor * fvec[i];
      }
    } else {
      double **farray = f->array_grid;
      int im1 = argindex - 1;
      for (i = 0; i < nglocal; i++)
        celldelta[i] = prefactor * farray[i][im1];
    }

  } else if (which == VARIABLE) {
    if (nglocal > maxvar) {
      maxvar = grid->maxlocal;
      memory->destroy(vbuf);
      memory->create(vbuf,maxvar,"ablate:vbuf");
    }

    input->variable->compute_grid(ivariable,vbuf,1,0);
    for (i = 0; i < nglocal; i++)
      celldelta[i] = prefactor * vbuf[i];
  }

  // NOTE: this does not get invoked on step 100,
  //   b/c needs to also be done in constructor
  //   ditto for fix adapt?
  //   they need nextvalid() methods like fix_ave_time
  //   or do it how output calcs next_stats for next thermo step

  modify->addstep_compute(update->ntimestep + nevery);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  double sum = 0.0;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;
    sum += celldelta[icell];
  }

  MPI_Allreduce(&sum,&sum_delta,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   decrement corner points of each owned grid cell
   skip cells not in group, with no surfs, and sub-cells
   algorithm:
     no corner pt value can be < 0.0
     decrement smallest corner pt by full delta
     if cannot, decrement to 0.0, decrement next smallest by remainder, etc
------------------------------------------------------------------------- */

void FixAblate::decrement()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,imin;
  double minvalue,total;
  double *corners;

  // total = full amount to decrement from cell
  // cdelta[icell] = amount to decrement from each corner point of icell

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++) cdelta[icell][i] = 0.0;

    total = celldelta[icell];
    corners = cvalues[icell];
    while (total > 0.0) {
      imin = -1;
      minvalue = 256.0;
      for (i = 0; i < ncorner; i++) {
        if (corners[i] > 0.0 && corners[i] < minvalue &&
            cdelta[icell][i] == 0.0) {
          imin = i;
          minvalue = corners[i];
        }
      }
      if (imin == -1) break;
      if (total < corners[imin]) {
        cdelta[icell][imin] += total;
        total = 0.0;
      } else {
        cdelta[icell][imin] = corners[imin];
        total -= corners[imin];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   find how much to decrement outside corner points
------------------------------------------------------------------------- */

void FixAblate::decrement_distributed_inside()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,j;
  int Nneigh, *neighbors, neigh1, neigh2;  
  double total_remain;
  double *corners;

  // to avoid couting same decrement multiple times in different cells

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++) cdelta[icell][i] = 0.0;

    // reference should be same as decrement_distributed_outside since
    // only outside points were updated

    if (dim == 2) setup_distributed2d(icell);
    else setup_distributed3d(icell);

    for (i = 0; i < ncorner; i++) {
      if (cvalues[icell][i] < 0) {

        // count which neighbors are inside and connected
        // record how much corner value to pass on if decrement was too large

        Nneigh = 0;
        if (dim == 2) {
          if (i == 0 || i == 3) {
            neigh1 = 1;
            neigh2 = 2;
          } else {
            neigh1 = 0;
            neigh2 = 3;
          }

          if (refcorners[neigh1] == 1) Nneigh++;
          if (refcorners[neigh2] == 1) Nneigh++;

          total_remain = fabs(cvalues[icell][i]) / Nneigh;
          if (refcorners[neigh1] == 1)
            cdelta[icell][neigh1] += total_remain;
          if (refcorners[neigh2] == 1)
            cdelta[icell][neigh2] += total_remain;

        } else {
          neighbors = corner_neighbor[i];
          for (j = 0; j < 3; j++)
            if (refcorners[neighbors[j]] == 1) Nneigh++;

          total_remain = fabs(cvalues[icell][i]) / Nneigh;
          for (j = 0; j < 3; j++)
            if(refcorners[neighbors[j]] == 1)
              cdelta[icell][neighbors[j]] += total_remain;
        } 

        // set negative corner value to zero
        cvalues[icell][i] = 0.0;

      } // end if for negative cvalues
    } // end corner
  } // end cells
}

/* ----------------------------------------------------------------------
   determine decrement for each corner point of each owned grid cell
   skip cells not in group, with no surfs, and sub-cells
   algorithm:
     no corner pt value can be < 0.0
     decrement smallest corner pt by full delta
     if cannot, decrement to 0.0, decrement next smallest by remainder, etc
------------------------------------------------------------------------- */

void FixAblate::decrement_adjacent()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,j,imin,jmin;
  double minvalue,total;
  double cmax[ncorner];

  // total = full amount to decrement from cell
  // cdelta[icell] = amount to decrement from each corner point of icell

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++) 
      for (j = 0; j < nadj; j++)
        cadelta[icell][i][j] = 0.0;

    total = celldelta[icell];
    while (total > 0.0) {

      // find max inner value
      for (i = 0; i < ncorner; i++) {
        cmax[i] = 0.0;
        for (j = 0; j < nadj; j++)
          cmax[i] = MAX(cmax[i],cavalues[icell][i][j]);
      }

      imin = -1;
      minvalue = 256.0;
      for (i = 0; i < ncorner; i++) {
        if (cmax[i] > 0.0 &&
            cadelta[icell][i][0] == 0.0 &&
            cmax[i] < minvalue) {
          imin = i;
          minvalue = cmax[i];
        }
      }

      if (imin == -1) break;
      if (total < cmax[imin]) {
        for (j = 0; j < nadj; j++)
          cadelta[icell][imin][j] += total;
        total = 0.0;
      } else {
        for (j = 0; j < nadj; j++)
          cadelta[icell][imin][j] = cmax[imin];
        total -= cmax[imin];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   find how much to decrement outside corner points
------------------------------------------------------------------------- */

void FixAblate::decrement_adjacent_distributed_inside2d()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,j,k,ind;
  int Nneigh, neigh1, neigh2;  
  double total_remain;

  // to avoid couting same decrement multiple times in different cells

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++)
        for (j = 0; j < nadj; j++) cadelta[icell][i][j] = 0.0;

    // reference should be same as decrement_distributed_outside since
    // only outside points were updated

    setup_distributed2d(icell);

    for (i = 0; i < ncorner; i++) {
      for (j = 0; j < nadj; j++) {
        if (cavalues[icell][i][j] < 0) {

          // count which neighbors are inside and connected
          // record how much corner value to pass on if decrement was too large

          Nneigh = 0;
          if (i == 0 || i == 3) {
            neigh1 = 1;
            neigh2 = 2;
          } else {
            neigh1 = 0;
            neigh2 = 3;
          }

          if (refcorners[neigh1] == 1) Nneigh++;
          if (refcorners[neigh2] == 1) Nneigh++;

          total_remain = fabs(cavalues[icell][i][j]) / Nneigh;

          if (refcorners[neigh1] == 1) {
            ind = find_inner(j,neigh1);
            cadelta[icell][neigh1][ind] += total_remain;
          }

          if (refcorners[neigh2] == 1) {
            ind = find_inner(j,neigh2);
            cadelta[icell][neigh2][ind] += total_remain;
          }

          // set negative inner value to zero
          cavalues[icell][i][j] = 0.0;

        } // end if for negative cvalues
      } // end inner
    } // end corner
  } // end cells
}

/* ----------------------------------------------------------------------
   find how much to decrement outside corner points
------------------------------------------------------------------------- */

void FixAblate::decrement_distributed_outside()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,j,Nout;
  double total,perout;

  // total = full amount to decrement from cell
  // cdelta[icell] = amount to decrement from each corner point of icell

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++) {
      if (adjacentflag)
        for (j = 0; j < nadj; j++) cadelta[icell][i][j] = 0.0;
      else
        cdelta[icell][i] = 0.0;
    }

    if (dim == 2) Nout = setup_distributed2d(icell);
    else Nout = setup_distributed3d(icell);
    /*printf("icell: %i; Nout: %i; ref: %i,%i,%i,%i,%i,%i,%i,%i\n",
      icell, Nout,
      refcorners[0],refcorners[1],refcorners[2],refcorners[3],
      refcorners[4],refcorners[5],refcorners[6],refcorners[7]);*/

    if (Nout == 0) continue; // all zero

    total = celldelta[icell];
    perout = total/Nout;
    //printf("celldelta: %4.2e; perout: %4.3e; Nout: %i\n", 
    //  total, perout, Nout);

    for (i = 0; i < ncorner; i++) {
      if (refcorners[i] != -1) continue;
      if (adjacentflag)
        for (j = 0; j < nadj; j++) cadelta[icell][i][j] += perout;
      else
        cdelta[icell][i] += perout;
    }

    /*for (i = 0; i < ncorner; i++) {
        for (j = 0; j < nadj; j++) 
          if(cadelta[icell][i][j] > 0.0)
            printf("cd[%i][%i][%i]: %4.2e\n", icell, i, j, cadelta[icell][i][j]);
    }*/

  } // end cells

  //error->one(FLERR,"ck dec out");
}

/* ----------------------------------------------------------------------
   sync and update corner values. some values may be negative. this is
   ok and will be handled in proceeding routines
------------------------------------------------------------------------- */

void FixAblate::sync_distributed_outside()
{
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,icorner,jcorner;
  int icell,jcell;
  double total[6];

  comm_neigh_corners(CDELTA);

  // perform update of corner pts for all my owned grid cells
  //   using contributions from all cells that share the corner point
  // insure order of numeric operations will give exact same answer
  //   for all Ncorner duplicates of a corner point (stored by other cells)

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    // loop over corner points

    for (i = 0; i < ncorner; i++) {

      // ixyz first = offset from icell of lower left cell of 2x2x2 stencil
      //              that shares the Ith corner point

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      // loop over 2x2x2 stencil of cells that share the corner point
      // also works for 2d, since izfirst = 0

      for (j = 0; j < 6; j++) total[j] = 0.0;
      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            // check if neighbor cell is within bounds of ablate grid

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            // jcell = local index of (jx,jy,jz) neighbor cell of icell

            jcell = walk_to_neigh(icell,jx,jy,jz);

            // update total with one corner point of jcell
            // jcorner descends from ncorner

            if (adjacentflag) {
              for (j = 0; j < nadj; j++) {
                if (jcell < nglocal) total[j] += cadelta[jcell][jcorner][j];
                else total[j] += cadelta_ghost[jcell-nglocal][jcorner][j];
              }
            } else {
              if (jcell < nglocal) total[0] += cdelta[jcell][jcorner];
              else total[0] += cdelta_ghost[jcell-nglocal][jcorner];
            }

          } // end jx
        } // end jy
      } // end jz

      /*for (j = 0; j < nadj; j++) {
        if (total[j] > 0)
          printf("total[%i][%i][%i]: %4.3e\n",
            icell, i, j, total[j]);
      }*/

      //for (j = 0; j < nadj; j++)
      //    printf("before - cval[%i][%i][%i]: %4.3e\n",
      //      icell, i, j, cavalues[icell][i][j]);

      if (adjacentflag)
        for (j = 0; j < nadj; j++) cavalues[icell][i][j] -= total[j];
      else
        cvalues[icell][i] -= total[0];

      //for (j = 0; j < nadj; j++)
      //    printf("after - cval[%i][%i][%i]: %4.3e\n",
      //      icell, i, j, cavalues[icell][i][j]);

    } // end corners
  } // end cells
  //error->one(FLERR,"ck sync out");
}


/* ----------------------------------------------------------------------
   find how much to decrement outside corner points
------------------------------------------------------------------------- */

void FixAblate::decrement_adjacent_distributed_inside3d()
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int i,j,k;
  int i_inner, i_neighbor_corner;
  int Nneigh, *ineighbors, *neighbors, oinner, kcorner;
  double total_remain;

  // to avoid couting same decrement multiple times in different cells

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (i = 0; i < ncorner; i++) {
      for (j = 0; j < nadj; j++) {
        cadelta[icell][i][j] = 0.0;
      }
    }

    // reference should be same as decrement_distributed_outside since
    // only outside points were updated

    setup_distributed3d(icell);

    for (i = 0; i < ncorner; i++) {

      // find neighboring corners and which are inside
      neighbors = corner_neighbor[i];
      Nneigh = 0;
      for (k = 0; k < 3; k++)
        if (refcorners[neighbors[k]] == 1) Nneigh++;

      //if (Nneigh > 0)
      //  printf("[%i][%i]: %i\n", icell, i, Nneigh);

      // three inner neighbors to test
      ineighbors = inner_neighbor[i];

      for (k = 0; k < 3; k++) {
        i_inner = ineighbors[k];
        i_neighbor_corner = neighbors[k];

        /*if (i_inner == 0) oinner = 1;
        else if (i_inner == 1) oinner = 0;
        else if (i_inner == 2) oinner = 3;
        else if (i_inner == 3) oinner = 2;
        else if (i_inner == 4) oinner = 5;
        else if (i_inner == 5) oinner = 4;*/
        //printf("%i: [%i] -> [%i] : inner: %i\n",
        //  icell, i, i_neighbor_corner, i_inner);

        if (cavalues[icell][i][i_inner] < 0) {

          // is the neighbor corner inside?
          if (refcorners[i_neighbor_corner] == 1) {

            // which inner index for connected corner point
            if (i_inner == 0) oinner = 1;
            else if (i_inner == 1) oinner = 0;
            else if (i_inner == 2) oinner = 3;
            else if (i_inner == 3) oinner = 2;
            else if (i_inner == 4) oinner = 5;
            else if (i_inner == 5) oinner = 4;

            //printf("%i: [%i][%i] -> [%i][%i]\n",
            //  icell, i, i_inner, i_neighbor_corner, oinner);

            total_remain = fabs(cavalues[icell][i][i_inner]) / Nneigh;
            cadelta[icell][i_neighbor_corner][oinner] += total_remain;
            //printf("cd: %4.3e\n", cadelta[icell][i_neighbor_corner][oinner]);
          }
        } // end if for negative cvalues
      } // end  k - neighbors

      // zero out negative values
      for (j = 0; j < nadj; j++)
        if (cavalues[icell][i][j] < 0) cavalues[icell][i][j] = 0.0;

    } // end corner

    /*for (i = 0; i < ncorner; i++) {
    for (j = 0; j < nadj; j++) {
      if(cadelta[icell][i][j] > 0.0)
        printf("cell: %i -- [%i][%i]: %4.3e\n",
          icell, i, j, cadelta[icell][i][j]);
    }
    }*/

  } // end cells
  //error->one(FLERR,"ck sync out");
}

/* ----------------------------------------------------------------------
   sync and update corner values. some values may be negative. this is
   ok and will be handled in proceeding routines
------------------------------------------------------------------------- */

void FixAblate::sync_adjacent_distributed_inside()
{
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,icorner,jcorner;
  int icell,jcell;
  int ndelta[nadj];
  double total[nadj];

  comm_neigh_corners(CDELTA);

  // perform update of corner pts for all my owned grid cells
  //   using contributions from all cells that share the corner point
  // insure order of numeric operations will give exact same answer
  //   for all Ncorner duplicates of a corner point (stored by other cells)

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    // loop over corner points

    for (i = 0; i < ncorner; i++) {

      // ixyz first = offset from icell of lower left cell of 2x2x2 stencil
      //              that shares the Ith corner point

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      // loop over 2x2x2 stencil of cells that share the corner point
      // also works for 2d, since izfirst = 0

      for (j = 0; j < nadj; j++) total[j] = 0.0;
      for (j = 0; j < nadj; j++) ndelta[j] = 0;
      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            // check if neighbor cell is within bounds of ablate grid

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            // jcell = local index of (jx,jy,jz) neighbor cell of icell

            jcell = walk_to_neigh(icell,jx,jy,jz);

            // update total with one corner point of jcell
            // jcorner descends from ncorner

            if (jcell < nglocal) {
              for (j = 0; j < nadj; j++) {
                if (cadelta[jcell][jcorner][j] > 0) {
                  total[j] += cadelta[jcell][jcorner][j];
                  ndelta[j] = ndelta[j] + 1;
                }
              }
            } else {
              for (j = 0; j < nadj; j++) {
                if (cadelta_ghost[jcell-nglocal][jcorner][j] > 0) {
                  total[j] += cadelta_ghost[jcell-nglocal][jcorner][j];
                  ndelta[j] = ndelta[j] + 1;
                }
              }
            }

          } // jz
        } // jy
      } // jx

      /*for (j = 0; j < nadj; j++)
        if (total[j] > 0)
          printf("final total[%i][%i][%i]: %4.3e -- %i\n",
            icell, i, j, total[j], ndelta[j]);*/

      for (j = 0; j < nadj; j++)
        if(ndelta[j] > 0) total[j] /= ndelta[j];

      /*if (icell <= 7) {
        for (j = 0; j < nadj; j++)
          printf("init cval[%i][%i][%i]: %4.3e -- %4.3e\n",
            icell, i, j, cavalues[icell][i][j], total[j]);
      }*/

      for (j = 0; j < nadj; j++) {
        cavalues[icell][i][j] -= total[j];
        if (cavalues[icell][i][j] < 0) cavalues[icell][i][j] = 0.0;
      }

      /*if (icell <= 7) {
        for (j = 0; j < nadj; j++)
          printf("final cval[%i][%i][%i]: %4.3e\n",
            icell, i, j, cavalues[icell][i][j]);
      }*/

    } // end corners
  } // end cells


  // check consistent
  /*int allin;
  int bad;
  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    // loop over corner points
 
    for (i = 0; i < ncorner; i++) {
      if(cavalues[icell][i][0] < thresh) allin = 0;
      else allin = 1;
      bad = 0;
      for (j = 0; j < nadj; j++) {
        if(cavalues[icell][i][j] > thresh && !allin) bad = 1;
        if(cavalues[icell][i][j] < thresh && allin) bad = 1;
      }
      if(bad) {
        for (j = 0; j < nadj; j++) {
          printf("cval: %4.3e\n", cavalues[icell][i][j]);
        }
        error->one(FLERR,"bad");
      }
    }
  }*/

  //error->one(FLERR,"ck sync adj in");
}


/* ----------------------------------------------------------------------
   lookup table for finding corresponding inner index for "end point"
   of a cell edge
------------------------------------------------------------------------- */

int FixAblate::find_inner(int istart, int iend)
{
  int ind = -1;

  if (dim == 2) {
    if (istart == 0) {
      if (iend == 1) ind = 0;
      else if (iend == 2) ind = 2;

    } else if (istart == 1) {
      if (iend == 0) ind = 1;
      else if (iend == 3) ind = 2;

    } else if (istart == 2) {
      if (iend == 0) ind = 3;
      else if (iend == 3) ind = 0;

    } else if (istart == 3) {
      if (iend == 1) ind = 3;
      else if (iend == 2) ind = 1;

    }
  } else {
    if (istart == 0) {
      if (iend == 1) ind = 0;
      else if (iend == 2) ind = 2;
      else if (iend == 4) ind = 4;

    } else if (istart == 1) {
      if (iend == 0) ind = 1;
      else if (iend == 3) ind = 2;
      else if (iend == 5) ind = 4;

    } else if (istart == 2) {
      if (iend == 0) ind = 3;
      else if (iend == 3) ind = 0;
      else if (iend == 6) ind = 4;

    } else if (istart == 3) {
      if (iend == 1) ind = 3;
      else if (iend == 2) ind = 1;
      else if (iend == 7) ind = 4;

    } else if (istart == 4) {
      if (iend == 5) ind = 0;
      else if (iend == 6) ind = 2;
      else if (iend == 0) ind = 5;

    } else if (istart == 5) {
      if (iend == 4) ind = 1;
      else if (iend == 7) ind = 2;
      else if (iend == 1) ind = 5;

    } else if (istart == 6) {
      if (iend == 4) ind = 3;
      else if (iend == 7) ind = 0;
      else if (iend == 2) ind = 5;

    } else if (istart == 7) {
      if (iend == 5) ind = 3;
      else if (iend == 6) ind = 1;
      else if (iend == 3) ind = 5;

    }
  }

  return ind;
}

/* ----------------------------------------------------------------------
   sync all copies of corner points values for all owned grid cells
   algorithm:
     comm my cdelta values that are shared by neighbor
     each corner point is shared by N cells, less on borders
     dsum = sum of decrements to that point by all N cells
     newvalue = MAX(oldvalue-dsum,0)
   all N copies of corner pt are set to newvalue
     in numerically consistent manner (same order of operations)
------------------------------------------------------------------------- */

void FixAblate::sync()
{
  int i,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner;
  int icell,jcell;
  double total;

  comm_neigh_corners(CDELTA);

  // perform update of corner pts for all my owned grid cells
  //   using contributions from all cells that share the corner point
  // insure order of numeric operations will give exact same answer
  //   for all Ncorner duplicates of a corner point (stored by other cells)

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    // loop over corner points

    for (i = 0; i < ncorner; i++) {

      // ixyz first = offset from icell of lower left cell of 2x2x2 stencil
      //              that shares the Ith corner point

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      // loop over 2x2x2 stencil of cells that share the corner point
      // also works for 2d, since izfirst = 0

      total = 0.0;
      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            // check if neighbor cell is within bounds of ablate grid

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            // jcell = local index of (jx,jy,jz) neighbor cell of icell

            jcell = walk_to_neigh(icell,jx,jy,jz);

            // update total with one corner point of jcell
            // jcorner descends from ncorner

            if (jcell < nglocal) total += cdelta[jcell][jcorner];
            else total += cdelta_ghost[jcell-nglocal][jcorner];
          }
        }
      }

      //if (total > 0)
      //  printf("final total[%i][%i]: %4.3e\n",
      //    icell, i, total);

      if (total > cvalues[icell][i]) cvalues[icell][i] = 0.0;
      else cvalues[icell][i] -= total;
    }
  }
}

/* ----------------------------------------------------------------------
   sync and update corner values. some values may be negative. this is
   ok and will be handled in proceeding routines
------------------------------------------------------------------------- */

void FixAblate::sync_distributed_inside()
{
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,icorner,jcorner;
  int icell,jcell;
  int ndelta;
  double total;

  comm_neigh_corners(CDELTA);

  // perform update of corner pts for all my owned grid cells
  //   using contributions from all cells that share the corner point
  // insure order of numeric operations will give exact same answer
  //   for all Ncorner duplicates of a corner point (stored by other cells)

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    // loop over corner points

    for (i = 0; i < ncorner; i++) {

      // ixyz first = offset from icell of lower left cell of 2x2x2 stencil
      //              that shares the Ith corner point

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      // loop over 2x2x2 stencil of cells that share the corner point
      // also works for 2d, since izfirst = 0

      total = 0.0;
      ndelta = 0;
      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            // check if neighbor cell is within bounds of ablate grid

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            // jcell = local index of (jx,jy,jz) neighbor cell of icell

            jcell = walk_to_neigh(icell,jx,jy,jz);

            // update total with one corner point of jcell
            // jcorner descends from ncorner

            if (jcell < nglocal) {
              if (cdelta[jcell][jcorner] > 0) {
                total += cdelta[jcell][jcorner];
                ndelta++;
              }
            } else {
              if (cdelta_ghost[jcell-nglocal][jcorner] > 0) {
                total += cdelta_ghost[jcell-nglocal][jcorner];
                ndelta++;
              }
            }

          }
        }
      }

      // accounts for double counting in sync_outside
      if(ndelta > 0) total /= (ndelta*0.5);

      cvalues[icell][i] -= total;
      if (cvalues[icell][i] < 0) error->all(FLERR,"negative corner value");

    } // end corners
  } // end cells
}

/* ----------------------------------------------------------------------
   sync all copies of corner points values for all owned grid cells
   algorithm for 2D
------------------------------------------------------------------------- */

int FixAblate::setup_distributed2d(int icell)
{
  for (int i = 0; i < ncorner; i++)
    refcorners[i] = 0;

  // mark inside corners first
  int Nin = 0;

  if(adjacentflag) {
    if (cavalues[icell][0][0] > thresh) {
      refcorners[0] = 1;
      Nin++;
    }

    if (cavalues[icell][1][0] > thresh) {
      refcorners[1] = 1;
      Nin++;
    }

    if (cavalues[icell][2][0] > thresh) {
      refcorners[3] = 1;
      Nin++;
    }

    if (cavalues[icell][3][0] > thresh) {
      refcorners[2] = 1;
      Nin++;
    }

  } else {
    if (cvalues[icell][0] > thresh) {
      refcorners[0] = 1;
      Nin++;
    }

    if (cvalues[icell][1] > thresh) {
      refcorners[1] = 1;
      Nin++;
    }

    if (cvalues[icell][2] > thresh) {
      refcorners[3] = 1;
      Nin++;
    }

    if (cvalues[icell][3] > thresh) {
      refcorners[2] = 1;
      Nin++;
    }
  }

  // built-in lookup table

  // 0 inside points -> 1 case
  int Nout;
  if (Nin == 0) {
    Nout = 4;
    refcorners[0] = refcorners[1] = refcorners[2] = refcorners[3] = -1.0;
  // 1 inside points -> 4 cases
  } else if (Nin == 1) {
    Nout = 2;
    if (refcorners[0] == 1) refcorners[2] = -1;
    else if (refcorners[1] == 1) refcorners[3] = -1;
    else if (refcorners[2] == 1) refcorners[0] = -1;
    else refcorners[1] = -1;
  }
  else if (Nin == 2) Nout = 2;
  else if (Nin == 3) Nout = 1;
  else Nout = 0;

  return Nout;
}


/* ----------------------------------------------------------------------
   sync all copies of corner points values for all owned grid cells
   algorithm 
------------------------------------------------------------------------- */

int FixAblate::setup_distributed3d(int icell)
{
  int i, bit[8], which;
  if (adjacentflag) {
    for (i = 0; i < ncorner; i++)
      bit[i] = cavalues[icell][i][0] <= thresh ? 0 : 1;
  } else {
    for (i = 0; i < ncorner; i++)
      bit[i] = cvalues[icell][i] <= thresh ? 0 : 1;
  }

  which = (bit[6] << 7) + (bit[7] << 6) + (bit[5] << 5) + (bit[4] << 4) +
    (bit[2] << 3) + (bit[3] << 2) + (bit[1] << 1) + bit[0];

  // find number of nonzero corner points and index of minimum corner
  // note: Nin + Nout does not always equal ncorner
  int Nin = Ninside[which];
  int Nout = Noutside[which];

  int *incorners, *outcorners;
  incorners = inside[which];
  outcorners = outside[which];

  // record which are inside and outside

  int icorner;
  for (i = 0; i < ncorner; i++)
    refcorners[i] = 0;

  for (i = 0; i < Nin; i++) {
    if (incorners[i] == 2) refcorners[3] = 1;
    else if (incorners[i] == 3) refcorners[2] = 1;
    else if (incorners[i] == 6) refcorners[7] = 1;
    else if (incorners[i] == 7) refcorners[6] = 1;
    else refcorners[incorners[i]] = 1;
  }

  for (i = 0; i < Nout; i++) {
    if (outcorners[i] == 2) refcorners[3] = -1;
    else if (outcorners[i] == 3) refcorners[2] = -1;
    else if (outcorners[i] == 6) refcorners[7] = -1;
    else if (outcorners[i] == 7) refcorners[6] = -1;
    else refcorners[outcorners[i]] = -1;
  }

  return Nout;
}

/* ----------------------------------------------------------------------
   sync all copies of adjacent values at each corner point
   for all owned grid cells ( very similar to sync() )
------------------------------------------------------------------------- */
void FixAblate::sync_adjacent()
{
  int i,j,ix,iy,iz,jx,jy,jz,ixfirst,iyfirst,izfirst,jcorner;
  int icell,jcell;
  int numin;
  double total[nadj];

  comm_neigh_corners(CDELTA);

  // perform update of corner pts for all my owned grid cells
  //   using contributions from all cells that share the corner point
  // insure order of numeric operations will give exact same answer
  //   for all Ncorner duplicates of a corner point (stored by other cells)

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    // loop over corner points

    for (i = 0; i < ncorner; i++) {

      /*-----------------------------------------------------------*/

      // ixyz first = offset from icell of lower left cell of 2x2x2 stencil
      //              that shares the Ith corner point

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      // loop over 2x2x2 stencil of cells that share the corner point
      // also works for 2d, since izfirst = 0

      for (j = 0; j < nadj; j++) total[j] = 0.0;
      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            // check if neighbor cell is within bounds of ablate grid

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            // jcell = local index of (jx,jy,jz) neighbor cell of icell

            jcell = walk_to_neigh(icell,jx,jy,jz);

            // update total with one corner point of jcell
            // jcorner descends from ncorner

            for (j = 0; j < nadj; j++) {
              if (jcell < nglocal) total[j] += cadelta[jcell][jcorner][j];
              else total[j] += cadelta_ghost[jcell-nglocal][jcorner][j];
            }
          }
        }
      }

      /*-----------------------------------------------------------*/

      // now decrement corners

      for(j = 0; j < nadj; j++) {
        cavalues[icell][i][j] -= total[j];
        if (cavalues[icell][i][j] < 0) cavalues[icell][i][j] = 0.0;
      }

      // check for consistency

      //numin = 0;
      //for(j = 0; j < nadj; j++)
      //  if (cavalues[icell][i][j] < thresh) numin++;
      //if (numin != 0 && numin != nadj) {
      //  for(j = 0; j < nadj; j++)
      //    printf("cval: %4.3e\n", cavalues[icell][i][j]);
      //  error->one(FLERR,"bad inner values - sync");
      //}

    }
  }
}

/* ----------------------------------------------------------------------
   adjust corner point values by epsilon of too close to threshold
   to avoid creating tiny or zero-size surface elements
------------------------------------------------------------------------- */

void FixAblate::epsilon_adjust()
{
  // insure no corner point is within EPSILON of threshold
  // if so, set it to threshold - EPSILON

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (int i = 0; i < ncorner; i++) {
      if (cvalues[icell][i] < cbufmin && cvalues[icell][i] >= thresh)
        cvalues[icell][i] = cbufmax;
      else if (cvalues[icell][i] > cbufmax && cvalues[icell][i] < thresh)
        cvalues[icell][i] = cbufmax;
    } // end corner
  } // end cells
}

/* ----------------------------------------------------------------------
   adjust corner point values by epsilon of too close to threshold
   to avoid creating tiny or zero-size surface elements
------------------------------------------------------------------------- */

void FixAblate::epsilon_adjust_adjacent()
{
  int allin, mixflag;
  double maxval;

  // insure no corner point is within EPSILON of threshold
  // if so, set it to threshold - EPSILON

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    for (int i = 0; i < ncorner; i++) {

      // check all in or out
      if (cavalues[icell][i][0] > thresh) allin = 1;
      else allin = 0;

      mixflag = 0;
      for (int j = 0; j < nadj; j++) {
        if (cavalues[icell][i][j] < cbufmin && allin) mixflag = 1;
        if (cavalues[icell][i][j] > thresh && !allin) mixflag = 1;
      }

      // if mixflag = 1, inner indices in disagreement in terms of side
      if (mixflag) {
        // find largest value from cbufmax
        maxval = -1.0;
        for (int j = 0; j < nadj; j++)
          maxval = MAX(maxval,cavalues[icell][i][j]-cbufmax);

        if (maxval < 0) error->one(FLERR,"negative maxval");
        for (int j = 0; j < nadj; j++) {
          cavalues[icell][i][j] -= maxval;
          cavalues[icell][i][j] = MAX(cavalues[icell][i][j],0.0);
        }               

      } else {
        for (int j = 0; j < nadj; j++) {
          if (cavalues[icell][i][j] < cbufmin &&
              cavalues[icell][i][j] >= thresh)
            cavalues[icell][i][j] = cbufmax;
          else if (cavalues[icell][i][j] > cbufmax &&
                   cavalues[icell][i][j] < thresh)
            cavalues[icell][i][j] = cbufmax;
        }
      }

    } // end corner
  } // end cells

  /*int bad;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    // loop over corner points
 
    for (int i = 0; i < ncorner; i++) {
      if(cavalues[icell][i][0] < thresh) allin = 0;
      else allin = 1;
      bad = 0;
      for (int j = 0; j < nadj; j++) {
        if(cavalues[icell][i][j] > thresh && !allin) bad = 1;
        if(cavalues[icell][i][j] < thresh && allin) bad = 1;
      }
      if(bad) {
        for (int j = 0; j < nadj; j++) {
          printf("cval: %4.3e\n", cavalues[icell][i][j]);
        }
        error->one(FLERR,"bad");
      }
    }
  }*/
}

/* ----------------------------------------------------------------------
   push corner points value to 0 or 255
     if all surrounding neighs are below or above threshold
     do this for all N copies of an affected corner point
   algorithm:
     comm my cdelta values that are shared by neighbor
     each corner point is shared by N cells, less on borders
     dsum = sum of decrements to that point by all N cells
     newvalue = MAX(oldvalue-dsum,0)
   all N copies of corner pt are set to newvalue
     in numerically consistent manner (same order of operations)
------------------------------------------------------------------------- */

void FixAblate::push_lohi()
{
  int i,ix,iy,iz,ixfirst,iyfirst,izfirst,jx,jy,jz;
  int icell,jcell,jcorner,pushflag;

  comm_neigh_corners(CVALUE);

  // perform push of corner pt values for all my owned grid cells
  //   by checking corner pt values of all cells that share same corner pt
  // if all surrounding corner pts are > threshold, push corner pt -> 255
  // if all surrounding corner pts are < threshold, push corner pt -> 0

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int plo = 0;
  int phi = 0;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    // loop over corner points

    for (i = 0; i < ncorner; i++) {

      // flag = -1 if corner pt value < threshold, +1 if > threshold

      if (cvalues[icell][i] < thresh) pushflag = -1;
      else if (cvalues[icell][i] > thresh) pushflag = 1;
      else continue;

      // ixyz first = offset from icell of lower left cell of 2x2x2 stencil
      //              that shares the Ith corner point

      ixfirst = (i % 2) - 1;
      iyfirst = (i/2 % 2) - 1;
      if (dim == 2) izfirst = 0;
      else izfirst = (i / 4) - 1;

      // loop over 2x2x2 stencil of cells that share the corner point
      // also works for 2d, since izfirst = 0

      jcorner = ncorner;

      for (jz = izfirst; jz <= izfirst+1; jz++) {
        for (jy = iyfirst; jy <= iyfirst+1; jy++) {
          for (jx = ixfirst; jx <= ixfirst+1; jx++) {
            jcorner--;

            // check if neighbor cell is within bounds of ablate grid

            if (ix+jx < 1 || ix+jx > nx) continue;
            if (iy+jy < 1 || iy+jy > ny) continue;
            if (iz+jz < 1 || iz+jz > nz) continue;

            // jcell = local index of (jx,jy,jz) neighbor cell of icell

            jcell = walk_to_neigh(icell,jx,jy,jz);

            // set pushflag to 0 if jcorner pt of jcell is not
            //   on same side of threshold as icorner or icell

            if (jcell < nglocal) {
              if ((pushflag == -1 && cvalues[jcell][jcorner] > thresh) ||
                  (pushflag == 1 && cvalues[jcell][jcorner] < thresh))
                pushflag = 0;
            } else {
              if ((pushflag == -1 && cdelta_ghost[jcell-nglocal][jcorner] >
                   thresh) ||
                  (pushflag == 1 && cdelta_ghost[jcell-nglocal][jcorner] <
                   thresh))
                pushflag = 0;
            }
          }
        }
      }

      // DEBUG OFF
      if (pushflag == -1) cvalues[icell][i] = 0;
      else if (pushflag == 1) cvalues[icell][i] = 255;

      if (pushflag == -1) plo++;
      else if (pushflag == 1) phi++;
    }
  }

  bigint bplo = plo;
  bigint bphi = phi;
  bigint ploall,phiall;
  MPI_Allreduce(&bplo,&ploall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&bphi,&phiall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);

  if (me == 0) {
    if (screen)
      fprintf(screen,"  " BIGINT_FORMAT " " BIGINT_FORMAT
              " pushed corner pt values\n",ploall,phiall);
    if (logfile)
      fprintf(logfile,"  " BIGINT_FORMAT " " BIGINT_FORMAT
              " pushed corner pt values\n",ploall,phiall);
  }
}

/* ----------------------------------------------------------------------
   comm my cdelta values that are shared by neighbor cells
   each corner point is shared by N cells, less on borders
   done via irregular comm
------------------------------------------------------------------------- */

void FixAblate::comm_neigh_corners(int which)
{
  int i,j,k,m,n,ix,iy,iz,jx,jy,jz;
  int icell,ifirst,jcell,proc,ilocal;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // make list of datums to send to neighbor procs
  // 8 or 26 cells surrounding icell need icell's cdelta info
  // but only if they are owned by a neighbor proc
  // insure icell is only sent once to same neighbor proc
  // also set proclist and locallist for each sent datum

  int nsend = 0;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];
    ifirst = nsend;

    // loop over 3x3x3 stencil of neighbor cells centered on icell

    for (jz = -1; jz <= 1; jz++) {
      for (jy = -1; jy <= 1; jy++) {
        for (jx = -1; jx <= 1; jx++) {

          // skip neigh = self

          if (jx == 0 && jy == 0 && jz == 0) continue;

          // check if neighbor cell is within bounds of ablate grid

          if (ix+jx < 1 || ix+jx > nx) continue;
          if (iy+jy < 1 || iy+jy > ny) continue;
          if (iz+jz < 1 || iz+jz > nz) continue;

          // jcell = local index of (jx,jy,jz) neighbor cell of icell

          jcell = walk_to_neigh(icell,jx,jy,jz);

          // add a send list entry of icell to proc != me if haven't already

          proc = cells[jcell].proc;
          if (proc != me) {
            for (j = ifirst; j < nsend; j++)
              if (proc == proclist[j]) break;
            if (j == nsend) {
              if (nsend == maxsend) grow_send();
              proclist[nsend] = proc;
              locallist[nsend++] = cells[icell].id;
            }
          }
        }
      }
    }

    // # of neighbor procs to send icell to

    numsend[icell] = nsend - ifirst;
  }

  // realloc sbuf if necessary
  // ncomm = ilocal + Ncorner values

  int ncomm;
  if(adjacentflag) ncomm = 1 + ncorner*nadj;
  else ncomm = 1 + ncorner;

  if (nsend*ncomm > maxbuf) {
    memory->destroy(sbuf);
    maxbuf = nsend*ncomm;
    memory->create(sbuf,maxbuf,"ablate:sbuf");
  }

  // pack datums to send
  // datum = ilocal of neigh cell on other proc + Ncorner values

  nsend = 0;
  m = 0;

  for (icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    n = numsend[icell];
    for (i = 0; i < n; i++) {
      sbuf[m++] = ubuf(locallist[nsend]).d;

      if (adjacentflag) {
        if (which == CDELTA) {
          for (j = 0; j < ncorner; j++)
            for (k = 0; k < nadj; k++)
              sbuf[m++] = cadelta[icell][j][k];
        } else if (which == CVALUE) {
          for (j = 0; j < ncorner; j++)
            for (k = 0; k < nadj; k++)
             sbuf[m++] = cavalues[icell][j][k];
        }
      } else {
        if (which == CDELTA) {
          for (j = 0; j < ncorner; j++)
            sbuf[m++] = cdelta[icell][j];
        } else if (which == CVALUE) {
          for (j = 0; j < ncorner; j++)
            sbuf[m++] = cvalues[icell][j];
        }
      }

      nsend++;
    }
  }

  // perform irregular neighbor comm
  // Comm class manages rbuf memory

  double *rbuf;
  int nrecv = comm->irregular_uniform_neighs(nsend,proclist,(char *) sbuf,
                                             ncomm*sizeof(double),
                                             (char **) &rbuf);

  // realloc cdelta_ghost if necessary

  if (grid->nghost > maxghost) {
    if (adjacentflag) {
      memory->destroy(cadelta_ghost);
      maxghost = grid->nghost;
      memory->create(cadelta_ghost,maxghost,ncorner,nadj,"ablate:cadelta_ghost");
    } else {
      memory->destroy(cdelta_ghost);
      maxghost = grid->nghost;
      memory->create(cdelta_ghost,maxghost,ncorner,"ablate:cdelta_ghost");
    }
  }

  // unpack received data into cdelta_ghost = ghost cell corner points

  // NOTE: need to check if hashfilled
  cellint cellID;
  Grid::MyHash *hash = grid->hash;

  m = 0;
  for (i = 0; i < nrecv; i++) {
    cellID = (cellint) ubuf(rbuf[m++]).u;
    ilocal = (*hash)[cellID];
    icell = ilocal - nglocal;
    if (adjacentflag)
      for (j = 0; j < ncorner; j++)
        for (k = 0; k < nadj; k++)
          cadelta_ghost[icell][j][k] = rbuf[m++];
    else
      for (j = 0; j < ncorner; j++)
        cdelta_ghost[icell][j] = rbuf[m++];
  }
}

/* ----------------------------------------------------------------------
   walk to neighbor of icell, offset by (jx,jy,jz)
   walk first by x, then by y, last by z
   return jcell = local index of neighbor cell
------------------------------------------------------------------------- */

int FixAblate::walk_to_neigh(int icell, int jx, int jy, int jz)
{
  Grid::ChildCell *cells = grid->cells;

  int jcell = icell;

  if (jx < 0) {
    if (grid->neigh_decode(cells[jcell].nmask,XLO) != NCHILD)
      error->one(FLERR,"Fix ablate walk to neighbor cell failed");
    jcell = cells[jcell].neigh[0];
  } else if (jx > 0) {
    if (grid->neigh_decode(cells[jcell].nmask,XHI) != NCHILD)
      error->one(FLERR,"Fix ablate walk to neighbor cell failed");
    jcell = cells[jcell].neigh[1];
  }

  if (jy < 0) {
    if (grid->neigh_decode(cells[jcell].nmask,YLO) != NCHILD)
      error->one(FLERR,"Fix ablate walk to neighbor cell failed");
    jcell = cells[jcell].neigh[2];
  } else if (jy > 0) {
    if (grid->neigh_decode(cells[jcell].nmask,YHI) != NCHILD)
      error->one(FLERR,"Fix ablate walk to neighbor cell failed");
    jcell = cells[jcell].neigh[3];
  }

  if (jz < 0) {
    if (grid->neigh_decode(cells[jcell].nmask,ZLO) != NCHILD)
      error->one(FLERR,"Fix ablate walk to neighbor cell failed");
    jcell = cells[jcell].neigh[4];
  } else if (jz > 0) {
    if (grid->neigh_decode(cells[jcell].nmask,ZHI) != NCHILD)
      error->one(FLERR,"Fix ablate walk to neighbor cell failed");
    jcell = cells[jcell].neigh[5];
  }

  return jcell;
}

/* ----------------------------------------------------------------------
   pack icell values for per-cell arrays into buf
   if icell is a split cell, also pack all sub cell values
   return byte count of amount packed
   if memflag, only return count, do not fill buf
------------------------------------------------------------------------- */

int FixAblate::pack_grid_one(int icell, char *buf, int memflag)
{
  char *ptr = buf;
  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  if (!adjacentflag) {
    if (memflag) memcpy(ptr,cvalues[icell],ncorner*sizeof(double));
    ptr += ncorner*sizeof(double);
  } else {
    for(int j = 0; j < ncorner; j++) {
      if (memflag) memcpy(ptr,cavalues[icell][j],nadj*sizeof(double));
      ptr += nadj*sizeof(double);
    }
  }

  if (tvalues_flag) {
    if (memflag) {
      double *dbuf = (double *) ptr;
      dbuf[0] = tvalues[icell];
    }
    ptr += sizeof(double);
  }

  if (memflag) {
    double *dbuf = (double *) ptr;
    dbuf[0] = ixyz[icell][0];
    dbuf[1] = ixyz[icell][1];
    dbuf[2] = ixyz[icell][2];
  }
  ptr += 3*sizeof(double);

  // DEBUG

  if (memflag) {
    double *dbuf = (double *) ptr;
    dbuf[0] = mcflags[icell][0];
    dbuf[1] = mcflags[icell][1];
    dbuf[2] = mcflags[icell][2];
    dbuf[3] = mcflags[icell][3];
  }
  ptr += 4*sizeof(double);

  if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    for (int i = 0; i < nsplit; i++) {
      int jcell = sinfo[isplit].csubs[i];
      
      if (!adjacentflag) {
        if (memflag) memcpy(ptr,cvalues[jcell],ncorner*sizeof(double));
        ptr += ncorner*sizeof(double);
      } else {
        for(int j = 0; j < ncorner; j++) {
          if (memflag) memcpy(ptr,cavalues[jcell][j],nadj*sizeof(double));
          ptr += nadj*sizeof(double);
        }
      }

    }
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   unpack icell values for per-cell array from buf
   return byte count of amount unpacked
------------------------------------------------------------------------- */

int FixAblate::unpack_grid_one(int icell, char *buf)
{
  char *ptr = buf;
  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;

  grow_percell(1);

  if (!adjacentflag) {
    memcpy(cvalues[icell],ptr,ncorner*sizeof(double));
    ptr += ncorner*sizeof(double);
  } else {
    for(int j = 0; j < ncorner; j++) {
      memcpy(cavalues[icell][j],ptr,nadj*sizeof(double));
      ptr += nadj*sizeof(double);
    }
  }

  if (tvalues_flag) {
    double *dbuf = (double *) ptr;
    tvalues[icell] = static_cast<int> (dbuf[0]);
    ptr += sizeof(double);
  }

  double *dbuf = (double *) ptr;
  ixyz[icell][0] = static_cast<int> (dbuf[0]);
  ixyz[icell][1] = static_cast<int> (dbuf[1]);
  ixyz[icell][2] = static_cast<int> (dbuf[2]);
  ptr += 3*sizeof(double);

  dbuf = (double *) ptr;
  mcflags[icell][0] = static_cast<int> (dbuf[0]);
  mcflags[icell][1] = static_cast<int> (dbuf[1]);
  mcflags[icell][2] = static_cast<int> (dbuf[2]);
  mcflags[icell][3] = static_cast<int> (dbuf[3]);
  ptr += 4*sizeof(double);

  nglocal++;

 if (cells[icell].nsplit > 1) {
    int isplit = cells[icell].isplit;
    int nsplit = cells[icell].nsplit;
    grow_percell(nsplit);
    for (int i = 0; i < nsplit; i++) {
      int jcell = sinfo[isplit].csubs[i];

      if (!adjacentflag) {
        memcpy(cvalues[jcell],ptr,ncorner*sizeof(double));
        ptr += ncorner*sizeof(double);
      } else {
        for(int j = 0; j < ncorner; j++) {
          memcpy(cavalues[jcell][j],ptr,nadj*sizeof(double));
          ptr += nadj*sizeof(double);
        }
      }

    }
    nglocal += nsplit;
  }

  return ptr-buf;
}

/* ----------------------------------------------------------------------
   copy per-cell info from Icell to Jcell
   called whenever a grid cell is removed from this processor's list
   caller checks that Icell != Jcell
------------------------------------------------------------------------- */

void FixAblate::copy_grid_one(int icell, int jcell)
{
  if (!adjacentflag) {
    memcpy(cvalues[jcell],cvalues[icell],ncorner*sizeof(double));
  } else {
    for (int j = 0; j < ncorner; j++)
      memcpy(cavalues[jcell][j],cavalues[icell][j],nadj*sizeof(double));
  }

  if (tvalues_flag) tvalues[jcell] = tvalues[icell];

  ixyz[jcell][0] = ixyz[icell][0];
  ixyz[jcell][1] = ixyz[icell][1];
  ixyz[jcell][2] = ixyz[icell][2];

  mcflags[jcell][0] = mcflags[icell][0];
  mcflags[jcell][1] = mcflags[icell][1];
  mcflags[jcell][2] = mcflags[icell][2];
  mcflags[jcell][3] = mcflags[icell][3];
}

/* ----------------------------------------------------------------------
   add a grid cell
   called when a grid cell is added to this processor's list
   initialize values to 0.0
------------------------------------------------------------------------- */

void FixAblate::add_grid_one()
{
  grow_percell(1);

  if (!adjacentflag) {
    for (int i = 0; i < ncorner; i++) cvalues[nglocal][i] = 0.0;
  } else {
    for (int i = 0; i < ncorner; i++)
      for (int j = 0; j < nadj; j++)
        cavalues[nglocal][i][j] = 0.0;
  }

  if (tvalues_flag) tvalues[nglocal] = 0;
  ixyz[nglocal][0] = 0;
  ixyz[nglocal][1] = 0;
  ixyz[nglocal][2] = 0;

  mcflags[nglocal][0] = -1;
  mcflags[nglocal][1] = -1;
  mcflags[nglocal][2] = -1;
  mcflags[nglocal][3] = -1;

  nglocal++;
}

/* ----------------------------------------------------------------------
   reset final grid cell count after grid cell removals
------------------------------------------------------------------------- */

void FixAblate::reset_grid_count(int nlocal)
{
  nglocal = nlocal;
}

/* ----------------------------------------------------------------------
   insure per-cell arrays are allocated long enough for Nnew cells
------------------------------------------------------------------------- */

void FixAblate::grow_percell(int nnew)
{
  if (nglocal+nnew < maxgrid) return;
  if (nnew == 0) maxgrid = nglocal;
  else maxgrid += DELTAGRID;
  if (adjacentflag) memory->grow(cavalues,maxgrid,ncorner,nadj,"ablate:cavalues");
  else memory->grow(cvalues,maxgrid,ncorner,"ablate:cvalues");
  if (tvalues_flag) memory->grow(tvalues,maxgrid,"ablate:tvalues");
  memory->grow(ixyz,maxgrid,3,"ablate:ixyz");
  memory->grow(mcflags,maxgrid,4,"ablate:mcflags");
  memory->grow(celldelta,maxgrid,"ablate:celldelta");
  if (adjacentflag) memory->grow(cadelta,maxgrid,ncorner,nadj,"ablate:cadelta");
  else memory->grow(cdelta,maxgrid,ncorner,"ablate:cdelta");
  memory->grow(numsend,maxgrid,"ablate:numsend");

  array_grid = cvalues;
}

/* ----------------------------------------------------------------------
   reallocate send vectors
------------------------------------------------------------------------- */

void FixAblate::grow_send()
{
  maxsend += DELTASEND;
  memory->grow(proclist,maxsend,"ablate:proclist");
  memory->grow(locallist,maxsend,"ablate:locallist");
}

/* ----------------------------------------------------------------------
   output sum of grid cell corner point values
   assume boundary corner points have value = 0.0
   NOTE: else would have to apply duplication weights to each of 4/8 corner pts
------------------------------------------------------------------------- */

double FixAblate::compute_scalar()
{
  int ix,iy,iz;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;

  double sum = 0.0;
  double cavg;
  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;

    ix = ixyz[icell][0];
    iy = ixyz[icell][1];
    iz = ixyz[icell][2];

    if (dim == 2 && (ix == 0 || iy == 0)) continue;
    if (dim == 3 && (ix == 0 || iy == 0 || iz == 0)) continue;

    if (!adjacentflag)
      sum += cvalues[icell][0];
    else {
      cavg = 0.0;
      for (int j = 0; j < nadj; j++) cavg += cavalues[icell][0][j];
      sum += cavg/nadj;
    }

  }

  double sumall;
  MPI_Allreduce(&sum,&sumall,1,MPI_DOUBLE,MPI_SUM,world);
  return sumall;
}

/* ----------------------------------------------------------------------
   vector outputs
   1 = last ablation decrement
   2 = # of deleted inside particles at last ablation
------------------------------------------------------------------------- */

double FixAblate::compute_vector(int i)
{
  if (i == 0) return sum_delta;
  if (i == 1) return 1.0*ndelete;
  return 0.0;
}

/* ----------------------------------------------------------------------
   process command line args
------------------------------------------------------------------------- */

void FixAblate::process_args(int narg, char **arg)
{
  distributedflag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if(strcmp(arg[iarg],"distributed") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid read_isurf command");
      if (strcmp(arg[iarg+1],"no") == 0) distributedflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) distributedflag = 1;
      else error->all(FLERR,"Illegal fix_ablate command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix_ablate command");
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double FixAblate::memory_usage()
{
  double bytes = 0.0;
  bytes += maxgrid*ncorner * sizeof(double);   // cvalues
  if (tvalues_flag) bytes += maxgrid * sizeof(int);   // tvalues
  bytes += maxgrid*3 * sizeof(int);            // ixyz
  // NOTE: add for mcflags if keep
  bytes += maxgrid * sizeof(double);           // celldelta
  bytes += maxgrid*ncorner * sizeof(double);   // cdelta
  bytes += maxghost*ncorner * sizeof(double);  // cdelta_ghost
  bytes += 3*maxsend * sizeof(int);            // proclist,locallist,numsend
  bytes += maxbuf * sizeof(double);            // sbuf
  return bytes;
}
