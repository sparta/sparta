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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include <array>
#include <map>
#include <algorithm>
#include "fix_rigid.h"
#include "update.h"
#include "domain.h"
#include "surf.h"
#include "grid.h"
#include "particle.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "compute_surf.h"
#include "input.h"
#include "geometry.h"
#include "cut2d.h"
#include "cut3d.h"
#include "math_extra.h"
#include "math_eigen.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

static constexpr double EPSILON = 1.0e-7;

#define INVOKED_PER_SURF 32
#define MAXLINE 1024
#define EPSSURF 1.0e-4          // same as Grid
#define EPSENCLOSED 1.0e-8      // min enclosed area/volume, relative
#define BIG 1.0e20
#define DELTA_MODIFY 1024

enum{INT,DOUBLE};                      // several files

enum{OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN};    // same as Update

// cell types, same as Grid
// renamed to avoid clash with surf collision side enum above

enum{CELLUNKNOWN,CELLOUTSIDE,CELLINSIDE,CELLOVERLAP};

enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain

enum{CUTCELL,INCREMENTAL};          // remap modes

// reasons incremental_recut() requests a full grid re-map

enum{FALLBACK_NONE,FALLBACK_NOPREV,FALLBACK_SPLIT,FALLBACK_SURFMAX};

enum{LINEAR,HERTZ};             // push-off force laws
enum{EULER,RICHARDSON};         // quaternion rotation update schemes

// local box/box overlap test, touching counts as overlap

static inline int box_overlap(double *alo, double *ahi,
                              double *blo, double *bhi)
{
  if (ahi[0] < blo[0] || alo[0] > bhi[0]) return 0;
  if (ahi[1] < blo[1] || alo[1] > bhi[1]) return 0;
  if (ahi[2] < blo[2] || alo[2] > bhi[2]) return 0;
  return 1;
}

/* ---------------------------------------------------------------------- */

FixRigid::FixRigid(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix rigid command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 22;
  global_freq = 1;
  nevery = 1;

  // gridmigrate insures grid_changed() is invoked when grid cells
  // are rebuilt or migrated, so incremental re-cut data can be reset

  gridmigrate = 1;

  if (!surf->exist) error->all(FLERR,"Fix rigid requires surf elements exist");
  kokkosable = 0;
  if (domain->axisymmetric)
    error->all(FLERR,"Fix rigid cannot be used with axisymmetric domains");
  if (surf->implicit)
    error->all(FLERR,"Fix rigid cannot be used with implicit surfs");

  // for distributed surfs, each proc owns a subset of the surfs;
  // the fix gathers a replicated copy of its (compact) body elements
  //   in gather_body() and maintains local Surf copies of them


  igroup = surf->find_group(arg[2]);
  if (igroup < 0) error->all(FLERR,"Fix rigid surf group ID does not exist");
  groupbit = surf->bitmask[igroup];
  
  int n = strlen(arg[3]) + 1;
  csurfID = new char[n];
  strcpy(csurfID,arg[3]);

  n = modify->find_compute(csurfID);
  if (n < 0) error->all(FLERR,"Fix rigid compute ID does not exist");

  // parse body params

  dim = domain->dimension;
  infile = NULL;
  slist = NULL;
  displace = NULL;

  forceinfile = 0;
  fcm_infile[0] = fcm_infile[1] = fcm_infile[2] = 0.0;
  torque_infile[0] = torque_infile[1] = torque_infile[2] = 0.0;

  int iarg = 4;
  if (strcmp(arg[iarg],"body") == 0) {
    if (iarg+22 > narg) error->all(FLERR,"Fix rigid body args not valid");
    massflag = comflag = vcomflag = moiflag = angmomflag = 0;
    int jarg = iarg+1;

    // NVALUE = # of args each keyword consumes, including the keyword
    // a keyword must not read past the 22 args this style is defined to
    //   take, which the check above showed are present

#define BODY_ARGS(nvalue)                                               \
    if (jarg+(nvalue) > iarg+22)                                        \
      error->all(FLERR,"Fix rigid body args not valid");

    while (jarg < iarg+22) {
      if (strcmp(arg[jarg],"mass") == 0) {
        BODY_ARGS(2);
	massflag = 1;
	massbody = input->numeric(FLERR,arg[jarg+1]);
	jarg += 2;
      } else if (strcmp(arg[jarg],"com") == 0) {
        BODY_ARGS(4);
	comflag = 1;
	xcm[0] = input->numeric(FLERR,arg[jarg+1]);
	xcm[1] = input->numeric(FLERR,arg[jarg+2]);
	xcm[2] = input->numeric(FLERR,arg[jarg+3]);
	jarg += 4;
      } else if (strcmp(arg[jarg],"moi") == 0) {
        BODY_ARGS(7);
	moiflag = 1;
	moi[0] = input->numeric(FLERR,arg[jarg+1]);
	moi[1] = input->numeric(FLERR,arg[jarg+2]);
	moi[2] = input->numeric(FLERR,arg[jarg+3]);
	moi[3] = input->numeric(FLERR,arg[jarg+4]);
	moi[4] = input->numeric(FLERR,arg[jarg+5]);
	moi[5] = input->numeric(FLERR,arg[jarg+6]);
	jarg += 7;
      } else if (strcmp(arg[jarg],"vcom") == 0) {
        BODY_ARGS(4);
	vcomflag = 1;
	vcm[0] = input->numeric(FLERR,arg[jarg+1]);
	vcm[1] = input->numeric(FLERR,arg[jarg+2]);
	vcm[2] = input->numeric(FLERR,arg[jarg+3]);
	jarg += 4;
      } else if (strcmp(arg[jarg],"angmom") == 0) {
        BODY_ARGS(4);
	angmomflag = 1;
	angmom[0] = input->numeric(FLERR,arg[jarg+1]);
	angmom[1] = input->numeric(FLERR,arg[jarg+2]);
	angmom[2] = input->numeric(FLERR,arg[jarg+3]);
	jarg += 4;
      } else
	error->all(FLERR,"Fix rigid body keyword not recognized");

    }

#undef BODY_ARGS
    if (!massflag || !comflag || !moiflag || !vcomflag || !angmomflag)
      error->all(FLERR,"Fix rigid body args not valid");
    iarg += 22;
    
  } else if (strcmp(arg[iarg],"infile") == 0) {
    if (iarg+2 > narg) error->all(FLERR,"Fix rigid infile args not valid");
    int n = strlen(arg[iarg+1]) + 1;
    infile = new char[n];
    strcpy(infile,arg[iarg+1]);
    read_infile(infile);
    iarg += 2;
    
  } else error->all(FLERR,"Fix rigid define style not recognized");

  // optional args

  outfile = NULL;
  outevery = 0;
  remapmode = INCREMENTAL;
  rotstyle = EULER;
  pushflag = 0;
  pushboundflag = 0;
  pushstyle = LINEAR;
  gammapush = 0.0;
  fext[0] = fext[1] = fext[2] = 0.0;
  int pushstyleflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"push") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Fix rigid body args not valid");
      pushflag = 1;
      kpush = input->numeric(FLERR,arg[iarg+1]);
      pushcutoff = input->numeric(FLERR,arg[iarg+2]);
      if (kpush < 0.0 || pushcutoff <= 0.0)
        error->all(FLERR,"Fix rigid body args not valid");
      iarg += 3;
    } else if (strcmp(arg[iarg],"pushbound") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Fix rigid body args not valid");
      if (strcmp(arg[iarg+1],"yes") == 0) pushboundflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pushboundflag = 0;
      else error->all(FLERR,"Fix rigid body args not valid");
      iarg += 2;
    } else if (strcmp(arg[iarg],"pushstyle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Fix rigid body args not valid");
      pushstyleflag = 1;
      if (strcmp(arg[iarg+1],"linear") == 0) pushstyle = LINEAR;
      else if (strcmp(arg[iarg+1],"hertz") == 0) pushstyle = HERTZ;
      else error->all(FLERR,"Fix rigid body args not valid");
      iarg += 2;
    } else if (strcmp(arg[iarg],"pushdamp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Fix rigid body args not valid");
      pushstyleflag = 1;
      gammapush = input->numeric(FLERR,arg[iarg+1]);
      if (gammapush < 0.0)
        error->all(FLERR,"Fix rigid body args not valid");
      iarg += 2;
    } else if (strcmp(arg[iarg],"remap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Fix rigid body args not valid");
      if (strcmp(arg[iarg+1],"cutcell") == 0) remapmode = CUTCELL;
      else if (strcmp(arg[iarg+1],"incremental") == 0)
        remapmode = INCREMENTAL;
      else error->all(FLERR,"Fix rigid body args not valid");
      iarg += 2;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Fix rigid body args not valid");
      if (strcmp(arg[iarg+1],"euler") == 0) rotstyle = EULER;
      else if (strcmp(arg[iarg+1],"richardson") == 0) rotstyle = RICHARDSON;
      else error->all(FLERR,"Fix rigid body args not valid");
      iarg += 2;
    } else if (strcmp(arg[iarg],"force") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Fix rigid body args not valid");
      fext[0] = input->numeric(FLERR,arg[iarg+1]);
      fext[1] = input->numeric(FLERR,arg[iarg+2]);
      fext[2] = input->numeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"outfile") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Fix rigid body args not valid");
      int n = strlen(arg[iarg+1]) + 1;
      outfile = new char[n];
      strcpy(outfile,arg[iarg+1]);
      outevery = input->inumeric(FLERR,arg[iarg+2]);
      if (outevery <= 0) error->all(FLERR,"Fix rigid body args not valid");
      iarg += 3;
    } else error->all(FLERR,"Fix rigid body args not valid");
  }

  if ((pushboundflag || pushstyleflag) && !pushflag)
    error->all(FLERR,"Fix rigid pushbound, pushstyle, and pushdamp "
               "require push keyword");

  if (massbody <= 0.0)
    error->all(FLERR,"Fix rigid body mass must be positive");

  // for 2d, insure all body params are consistent with in-plane motion

  if (dim == 2) {
    if (xcm[2] != 0.0 || vcm[2] != 0.0)
      error->all(FLERR,"Fix rigid z components of com and vcom "
                 "must be zero for 2d");
    if (angmom[0] != 0.0 || angmom[1] != 0.0)
      error->all(FLERR,"Fix rigid x,y components of angmom "
                 "must be zero for 2d");
    if (moi[4] != 0.0 || moi[5] != 0.0)
      error->all(FLERR,"Fix rigid ixz,iyz components of moi "
                 "must be zero for 2d");
    if (fext[2] != 0.0)
      error->all(FLERR,"Fix rigid z component of force must be zero for 2d");
  }

  // setup the rigid body

  setup_body();

  // restore the force/torque of the step before a continuation, which
  //   setup_body() zeroed; the body is moved by them on the first step

  if (forceinfile) {
    for (int j = 0; j < 3; j++) {
      fcm[j] = fcm_infile[j];
      torque[j] = torque_infile[j];
    }
  }

  // irigid = per-surf flags, indexed by local surf index
  // -1 for static surfs, else index into slist of body surfs
  // used by Update::build_rigidmap() to detect moving surfs
  // non-distributed only: every proc stores all surfs, length = nlocal
  // for distributed surfs the local surf list changes as the body
  //   moves, so build_rigidmap() maps by global surf ID instead

  irigid = NULL;
  nsurfall = surf->nlocal;

  if (!surf->distributed) {
    int nslocal = surf->nlocal;
    memory->create(irigid,nslocal,"fix_rigid:irigid");
    for (int i = 0; i < nslocal; i++) irigid[i] = -1;
    for (int i = 0; i < nsurf; i++) irigid[slist[i]] = i;
  }

  // remap data structs
  // body surfs are cut/split into grid cells by the normal surf
  //   pipeline (done at read_surf time), so no special setup is needed
  //   here beyond the swept-assignment and incremental work buffers

  ndeleted = 0;
  ndeleted_all = 0;
  ndelvalid = -1;

  // elemlo/elemhi are allocated by setup_body() above

  nmodified = maxmodified = 0;
  modified = NULL;
  nsurf_saved = NULL;
  csurfs_saved = NULL;
  cpage = NULL;             // allocated in setup(), needs all bodies' sizes

  pbodyflag = 0;
  noldinside = maxoldinside = 0;
  oldinside = NULL;
  nrcand = maxrcand = 0;
  rcand = NULL;
  newlist = NULL;
  newmap = NULL;
  reclist = NULL;
  maxreclist = 0;
  maxnewlist = 0;
  cut2d = NULL;
  cut3d = NULL;

  pushbinstart = NULL;
  pushbinlist = NULL;
  pushstamp = NULL;
  pushstampcur = 0;
  ftbuf_mine = ftbuf_all = NULL;
  tqpush[0] = tqpush[1] = tqpush[2] = 0.0;
  warnfallback = 0;
  warndelete = 0;
  ndelrun = 0;

  swstamp = NULL;
  swhead = NULL;
  maxswcell = 0;
  swcur = 0;
  swcells = NULL;
  nswcell = maxswcells = 0;
  entnext = NULL;
  entelem = NULL;
  nent = maxent = 0;

  // for incremental mode: cutters and work bufs for re-cutting cells

  // work bufs are sized per run in setup(), since global surfmax can be
  //   changed between runs, after this fix is defined

  if (remapmode == INCREMENTAL) {
    if (dim == 2) cut2d = new Cut2d(sparta,0);
    else cut3d = new Cut3d(sparta);
  }
}

/* ---------------------------------------------------------------------- */
 
FixRigid::~FixRigid()
{
  delete [] csurfID;
  delete [] infile;
  delete [] outfile;
  memory->destroy(slist);
  memory->destroy(displace);
  memory->destroy(irigid);
  memory->destroy(elemlo);
  memory->destroy(elemhi);
  memory->destroy(bodypt);
  memory->destroy(bodynorm);
  memory->destroy(sids);
  memory->destroy(bodymask);
  memory->destroy(bodytype);
  memory->destroy(bodytrans);
  memory->destroy(bodyisc);
  memory->destroy(bodyisr);
  memory->destroy(lblist);
  memory->destroy(copy_index);
  memory->destroy(copy_elem);
  memory->destroy(olist_own);
  memory->destroy(olist_elem);
  memory->destroy(modified);
  memory->destroy(nsurf_saved);
  memory->sfree(csurfs_saved);
  delete cpage;

  // registry csurfs lists may still be installed in live grid cells,
  //   e.g. this fix is unfixed between runs after incremental re-cuts
  // copy each such list into grid-owned page storage before freeing it
  // Modify is destroyed before Grid at program teardown, so grid is valid

  copy_registry_to_grid();
  free_registry();
  memory->destroy(oldinside);
  memory->destroy(rcand);
  memory->destroy(newlist);
  memory->destroy(newmap);
  memory->destroy(reclist);
  delete cut2d;
  delete cut3d;

  memory->destroy(pushbinstart);
  memory->destroy(pushbinlist);
  memory->destroy(pushstamp);
  memory->destroy(ftbuf_mine);
  memory->destroy(ftbuf_all);
  memory->destroy(swstamp);
  memory->destroy(swhead);
  memory->destroy(swcells);
  memory->destroy(entnext);
  memory->destroy(entelem);
}

/* ---------------------------------------------------------------------- */

int FixRigid::setmask()
{
  int mask = 0;
  mask |= START_OF_STEP;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRigid::init()
{
  // check that global rigid flag is set

  if (update->rigidflag == 0)
    error->all(FLERR,"Cannot use fix rigid unless global rigid is set");

  // with the KOKKOS package active, only the rigid/kk variant brackets
  //   its host-side work with the required device transfers

  if (sparta->kokkos && !kokkosable)
    error->all(FLERR,"Must use fix rigid/kk with KOKKOS");

  // the recoil correction of collisions is exact down to a body as
  //   light as one simulation particle; below that the exact result
  //   would carry the particle past the wall within the step, which the
  //   one-step coupling cannot represent, so the correction is inactive

  double mmax = 0.0;
  for (int i = 0; i < particle->nspecies; i++)
    mmax = MAX(mmax,particle->species[i].mass);
  if (massbody < update->fnum*mmax && comm->me == 0)
    error->warning(FLERR,"Fix rigid body mass is less than the mass of a "
                   "simulation particle, collisions are not corrected "
                   "for body recoil");

  // check that specified compute is valid for use with fix rigid
  // NOTE: check that it operates on same surf group ?

  int n = modify->find_compute(csurfID);
  if (n < 0) error->all(FLERR,"Could not find fix rigid compute ID");
  if (strcmp(modify->compute[n]->style,"surf") != 0 &&
      strcmp(modify->compute[n]->style,"surf/kk") != 0)
    error->all(FLERR,"Fix rigid compute is not style surf");
  csurf = (ComputeSurf *) modify->compute[n];
  if (csurf->per_surf_flag == 0)
    error->all(FLERR,"Fix rigid compute does not compute per-surf info");
  if (csurf->size_per_surf_cols != 6 || !csurf->force_torque_colcheck())
    error->all(FLERR,"Fix rigid compute must tally exactly "
               "fx fy fz tx ty tz for a single group");

  // insure the compute tallies on the first step of the next run
  // end_of_step() extends this to every step of the run

  csurf->addstep(update->ntimestep+1);

  // body surfs cannot be transparent or have surface reactions assigned
  // all body surfs must be in the surf group tallied by the compute
  // attributes come from the replicated body table, valid for both
  //   non-distributed and distributed surfs

  int cbit = csurf->surf_groupbit();

  for (int i = 0; i < nsurf; i++) {
    if (bodytrans[i])
      error->all(FLERR,"Fix rigid body surfs cannot be transparent");
    if (bodyisr[i] >= 0)
      error->all(FLERR,"Fix rigid body surfs cannot have surface reactions");
    if (!(bodymask[i] & cbit))
      error->all(FLERR,"Fix rigid compute surf group does not include "
                 "all body surfs");
  }

  // smallest grid cell edge length, for motion-rate warnings

  Grid::ChildCell *cells = grid->cells;
  int nglocal = grid->nlocal;

  double mine = BIG;
  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    mine = MIN(mine,cells[icell].hi[0]-cells[icell].lo[0]);
    mine = MIN(mine,cells[icell].hi[1]-cells[icell].lo[1]);
    if (dim == 3) mine = MIN(mine,cells[icell].hi[2]-cells[icell].lo[2]);
  }
  MPI_Allreduce(&mine,&mincellsize,1,MPI_DOUBLE,MPI_MIN,world);

  // re-enable single-shot warnings for this run

  warnrotate = warntranslate = warnexit = warnfallback = 0;
  warndelete = 0;
  ndelrun = 0;

  // surfs cannot change once a fix rigid is defined:
  //   removal invalidates the body element table; a change to the
  //   body group invalidates the body definition
  // surfs appended after the fix was defined are allowed: grow irigid
  //   and flag them static
  // Update::init() clamps its rigidmap scan to nsurfall, so it is
  //   correct even though it runs before this method

  if (surf->count_group(igroup) != nsurf)
    error->all(FLERR,"Fix rigid body surf group was changed "
               "after fix rigid was defined");

  if (!surf->distributed) {
    if (surf->nlocal < nsurfall)
      error->all(FLERR,"Surfs were removed after fix rigid was defined");
    if (surf->nlocal > nsurfall) {
      memory->grow(irigid,surf->nlocal,"fix_rigid:irigid");
      for (int i = nsurfall; i < surf->nlocal; i++) irigid[i] = -1;
      nsurfall = surf->nlocal;
    }
  }

  // each fix rigid defines its own body: no surf can be in two bodies
  // each fix rigid must have its own compute: the fix resets the
  //   compute's COM to its own body COM every step, so a shared compute
  //   would tally torques about the wrong body's COM

  for (int ifix = 0; ifix < modify->nfix; ifix++) {
    if (modify->fix[ifix] == this) continue;
    if (!fix_rigid_style(modify->fix[ifix]->style)) continue;
    FixRigid *other = (FixRigid *) modify->fix[ifix];
    for (int i = 0; i < nsurf; i++)
      if (other->body_elem(sids[i]) >= 0)
        error->all(FLERR,"Surf element is in more than one fix rigid body");
    if (strcmp(other->csurfID,csurfID) == 0)
      error->all(FLERR,"Two fix rigid commands cannot use the same compute");
  }

  // fix rigid must be defined before fixes which change the grid,
  // so its end_of_step() restores overlaid grid cells before they run

  int myindex = modify->find_fix(id);
  for (int ifix = 0; ifix < myindex; ifix++)
    if (strncmp(modify->fix[ifix]->style,"balance",7) == 0 ||
        strncmp(modify->fix[ifix]->style,"adapt",5) == 0)
      error->all(FLERR,
                 "Fix rigid must be defined before fix balance or fix adapt");
}

/* ----------------------------------------------------------------------
   called at start of each run, after grid and particles are setup
------------------------------------------------------------------------- */

void FixRigid::setup()
{
  FixRigid **flist = update->fixrigidlist;
  int nb = update->nfixrigid;

  // work buffer for the fused force/torque Allreduce over all bodies
  // the # of rigid fixes can change between runs, so (re)allocate

  memory->destroy(ftbuf_mine);
  memory->destroy(ftbuf_all);
  memory->create(ftbuf_mine,6*nb,"fix_rigid:ftbuf_mine");
  memory->create(ftbuf_all,6*nb,"fix_rigid:ftbuf_all");

  // page for merged csurfs lists built by swept_assign_all each step
  // owned by the last-defined rigid fix, which performs the single
  //   cell pass for all bodies; a merged list can hold one cell's
  //   current surfs plus the swept surfs of every body

  delete cpage;
  cpage = NULL;
  if (flist[nb-1] == this) {
    int nsurftotal = 0;
    for (int m = 0; m < nb; m++) nsurftotal += flist[m]->nsurf;
    int maxchunk = grid->maxsurfpercell + nsurftotal;
    cpage = new MyPage<surfint>(maxchunk,MAX(65536,4*maxchunk));
    if (cpage->errorflag)
      error->all(FLERR,"Fix rigid could not allocate collision-list page");
  }

  // distributed surfs: insure local copies of body surfs exist and
  //   rigidmap covers them, before the first step's swept assignment
  // rigidmap also spans the ghost surfs acquired for this run, and
  //   per-surf computes must size their arrays for any appended copies

  if (surf->distributed) {
    ensure_local_copies();
    update->build_rigidmap();
    surfs_changed();
  }

  // work bufs for incremental re-cutting of one cell, sized for the
  //   current global surfmax, which can change between runs:
  // newlist/newmap = the cell's new surf list and its split map, both
  //   capped at maxsurfpercell by the cut routines
  // reclist = candidate surfs, the cell's current static surfs plus
  //   every element of every body

  if (remapmode == INCREMENTAL) {
    if (grid->maxsurfpercell > maxnewlist) {
      maxnewlist = grid->maxsurfpercell;
      memory->destroy(newlist);
      memory->destroy(newmap);
      memory->create(newlist,maxnewlist,"fix_rigid:newlist");
      memory->create(newmap,maxnewlist,"fix_rigid:newmap");
    }

    int nsurftotal = 0;
    for (int m = 0; m < nb; m++) nsurftotal += flist[m]->nsurf;
    int n = grid->maxsurfpercell + nsurftotal;
    if (n > maxreclist) {
      maxreclist = n;
      memory->destroy(reclist);
      memory->create(reclist,maxreclist,"fix_rigid:reclist");
    }
  }

  // bin static surfs for push-off candidate pruning

  if (pushflag) push_bins();

  // delete any particles inside the body
  // create_particles marks the body's cells INSIDE via the surf pipeline
  //   and normally avoids them, but this is a safety net for any that
  //   end up inside, e.g. via an emit region overlapping the body

  if (particle->exist) ndeleted += remove_inside_particles(0);

  // for incremental remap: grid state is now consistent with the
  //   body at its current position

  if (remapmode == INCREMENTAL) {
    body_bbox(0);
    for (int j = 0; j < 3; j++) {
      pbodylo[j] = bbodylo[j];
      pbodyhi[j] = bbodyhi[j];
    }
    pbodyflag = 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixRigid::start_of_step()
{
  // reset COM used by compute surf for torque tallies to current COM
  // torque is thus about the start-of-step COM,
  //   consistent with the accuracy of the time integration below

  csurf->set_com(xcm);

  // body inverse mass and inertia for collision recoil this step,
  //   from the start-of-step axes before they are advanced below

  set_recoil();

  // time integrate from current position to end-of-step position
  // velocity Verlet: this is the first half kick and the drift,
  //   the second half kick is applied in end_of_step() once the force
  //   and torque of this step are known
  // fcm,torque = from particle collisions and push-off contacts during
  //   the previous step, i.e. the force at the start of this step,
  //   plus the constant external force
  // vcm/angmom/omega are thus half-step values during the step: the
  //   body moves, and particles collide with it, at these velocities,
  //   which is second-order accurate and exact for a constant force
  // xcmnew/quatnew/exyz_space = end-of-step values

  double dt = update->dt;
  double dtfhalf = 0.5 * dt / massbody;
  double dthalf = 0.5 * dt;

  vcm[0] += dtfhalf * (fcm[0] + fext[0]);
  vcm[1] += dtfhalf * (fcm[1] + fext[1]);
  vcm[2] += dtfhalf * (fcm[2] + fext[2]);

  // drift xcm by full step with the half-step velocity
  // store as xcmnew so have start/stop position for this timestep

  xcmnew[0] = xcm[0] + dt * vcm[0];
  xcmnew[1] = xcm[1] + dt * vcm[1];
  xcmnew[2] = xcm[2] + dt * vcm[2];

  // half kick of angular momentum in spatial frame

  angmom[0] += dthalf * torque[0];
  angmom[1] += dthalf * torque[1];
  angmom[2] += dthalf * torque[2];

  // compute new omega from new angmom, both in spatial frame

  MathExtra::angmom_to_omega(angmom,ex_space,ey_space,ez_space,inertia,omega);

  // for 2d, insure COM stays in plane and rotation is about z axis
  // guards against small numeric drift in principal axes

  if (dim == 2) {
    xcmnew[2] = 0.0;
    omega[0] = 0.0;
    omega[1] = 0.0;
  }

  // update quaternion by full step using new omega in spatial frame
  // store as quatnew so have start/stop orientation for this timestep
  // rotate euler (default): omega is held constant over the step, so
  //   dq/dt = 1/2 omega q integrates exactly to the rotation by
  //   angle |omega|*dt about omega; the moving-surf collision tests
  //   assume this same rotation, so the end-of-step geometry the
  //   particles were reflected from is exactly the one installed
  // rotate richardson: LAMMPS-style Richardson iteration which
  //   re-evaluates omega at the half step from the (constant over the
  //   step) angular momentum; useful for rotation-dominated bodies

  if (rotstyle == RICHARDSON) {
    quatnew[0] = quat[0];
    quatnew[1] = quat[1];
    quatnew[2] = quat[2];
    quatnew[3] = quat[3];
    MathExtra::richardson(quatnew,angmom,omega,inertia,dthalf);
  } else {
    double wmag = MathExtra::len3(omega);
    if (wmag > 0.0) {
      double axis[3],qrot[4];
      axis[0] = omega[0]/wmag;
      axis[1] = omega[1]/wmag;
      axis[2] = omega[2]/wmag;
      MathExtra::axisangle_to_quat(axis,wmag*dt,qrot);
      MathExtra::quatquat(qrot,quat,quatnew);
      MathExtra::qnormalize(quatnew);
    } else {
      quatnew[0] = quat[0];
      quatnew[1] = quat[1];
      quatnew[2] = quat[2];
      quatnew[3] = quat[3];
    }
  }
  MathExtra::q_to_exyz(quatnew,ex_space,ey_space,ez_space);

  // warn once per run if body motion in a single step is too large
  // rotation > 0.1 radian degrades the chord approximation used for
  //   collisions of particles with rotating surfs
  // max surf pt displacement > smallest grid cell degrades the
  //   accuracy of surf assignment to grid cells for cutcell remapping

  if (!warnrotate && MathExtra::len3(omega)*dt > 0.1) {
    warnrotate = 1;
    if (comm->me == 0)
      error->warning(FLERR,"Fix rigid body rotation per timestep exceeds "
                     "0.1 radian, collision accuracy degrades");
  }

  if (!warntranslate) {
    double dispmax =
      (MathExtra::len3(vcm) + MathExtra::len3(omega)*rmaxbody) * dt;
    if (dispmax > mincellsize) {
      warntranslate = 1;
      if (comm->me == 0)
        error->warning(FLERR,"Fix rigid body moves more than a grid cell "
                       "per timestep, cell assignment accuracy degrades");
    }
  }

  // augment collision lists of all cells any body sweeps through during
  //   the step, so particles in the swept paths are tested against the
  //   moving surfs and reflected rather than overtaken and later deleted
  // one pass over grid cells for all bodies, by the last-defined fix:
  //   start_of_step runs fixes in definition order, so when the last
  //   fix runs every body's end-of-step pose is known

  if (update->fixrigidlist[update->nfixrigid-1] == this) swept_assign_all();
}

/* ---------------------------------------------------------------------- */

void FixRigid::end_of_step()
{
  FixRigid **flist = update->fixrigidlist;
  int nb = update->nfixrigid;

  // the first-defined rigid fix coordinates two all-body operations,
  //   before any fix reads per-cell surf lists or per-surf tallies:
  // (1) undo the swept collision-list augmentation from start_of_step;
  //     the last-defined fix installed it and owns the bookkeeping
  // (2) sum per-surf force/torque tallies to fcm/torque of every body,
  //     fused into a single Allreduce of 6 values per body
  // end_of_step runs fixes in definition order, so the first fix runs
  //   before any other fix's end_of_step touches the grid

  if (flist[0] == this) {

    flist[nb-1]->swept_restore();

    // sum per-surf force/torque to each body's fcm/torque
    // read the compute's RAW local tally rows: values are fully
    //   normalized at tally time, and each row's surf ID maps to a
    //   body element via the body's ID table, so a local sum plus the
    //   single fused Allreduce below is exactly the collated result
    // this avoids Surf::collate_array entirely, whose reduce path is
    //   an Allreduce over ALL global surfs per compute per step, and
    //   avoids any scan over the surf list: cost is O(local tallies)
    // identical for non-distributed and distributed surfs

    for (int i = 0; i < 6*nb; i++) ftbuf_mine[i] = 0.0;

    surfint *t2s;
    for (int m = 0; m < nb; m++) {
      FixRigid *f = flist[m];

      ComputeSurf *cs = f->csurf;
      if (!(cs->invoked_flag & INVOKED_PER_SURF)) {
        cs->compute_per_surf();
        cs->invoked_flag |= INVOKED_PER_SURF;
      }

      int ntally = cs->tallyinfo(t2s);
      double **tally = cs->tally_array();

      double *ft = &ftbuf_mine[6*m];
      for (int i = 0; i < ntally; i++) {
        if (f->body_elem(t2s[i]) < 0) continue;
        for (int j = 0; j < 6; j++) ft[j] += tally[i][j];
      }

      // insure the compute tallies on the next step

      cs->addstep(update->ntimestep+1);
    }

    MPI_Allreduce(ftbuf_mine,ftbuf_all,6*nb,MPI_DOUBLE,MPI_SUM,world);

    for (int m = 0; m < nb; m++) {
      FixRigid *f = flist[m];
      f->fcm[0] = ftbuf_all[6*m];
      f->fcm[1] = ftbuf_all[6*m+1];
      f->fcm[2] = ftbuf_all[6*m+2];
      f->torque[0] = ftbuf_all[6*m+3];
      f->torque[1] = ftbuf_all[6*m+4];
      f->torque[2] = ftbuf_all[6*m+5];
    }
  }

  // for incremental remap: record cells interior to the body
  //   before its surfs move to their end-of-step positions

  if (remapmode == INCREMENTAL) record_oldinside();

  // reset xcm/quat to new xcm/quat calculated in start_of_step()

  xcm[0] = xcmnew[0];
  xcm[1] = xcmnew[1];
  xcm[2] = xcmnew[2];

  quat[0] = quatnew[0];
  quat[1] = quatnew[1];
  quat[2] = quatnew[2];
  quat[3] = quatnew[3];

  // enforce2d on all body properties
  // NOTE: should we also enforce this in start_of_step() for xcmnew,quatnew,omega ?
  
  if (dim == 2) {
    xcm[2] = 0.0;
    vcm[2] = 0.0;
    fcm[2] = 0.0;
    torque[0] = 0.0;
    torque[1] = 0.0;
    angmom[0] = 0.0;
    angmom[1] = 0.0;
    omega[0] = 0.0;
    omega[1] = 0.0;
    // what about quat for 2d rotations ?
  }

  // regenerate the replicated body geometry from the new pose:
  //   corner pts from displace rotated to the space frame + new COM,
  //   normals recomputed from the corner pts
  // then write it into the Surf copies the mover and cut pipeline read
  // matvec() converts displace vector from body frame to space frame

  {
    double z[3],delta[3],delta12[3],delta13[3];
    z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;

    for (int i = 0; i < nsurf; i++) {
      for (int j = 0; j < dim; j++) {
        MathExtra::matvec(ex_space,ey_space,ez_space,displace[i][j],delta);
        if (dim == 2) delta[2] = 0.0;
        MathExtra::add3(xcm,delta,bodypt[i][j]);
      }

      if (dim == 2) {
        MathExtra::sub3(bodypt[i][1],bodypt[i][0],delta);
        MathExtra::cross3(z,delta,bodynorm[i]);
        MathExtra::norm3(bodynorm[i]);
        bodynorm[i][2] = 0.0;
      } else {
        MathExtra::sub3(bodypt[i][1],bodypt[i][0],delta12);
        MathExtra::sub3(bodypt[i][2],bodypt[i][0],delta13);
        MathExtra::cross3(delta12,delta13,bodynorm[i]);
        MathExtra::norm3(bodynorm[i]);
      }
    }
  }

  update_surf_copies();

  // bbox around body elements at their new positions

  body_bbox(0);

  // push-off forces are computed for all bodies at once by the
  //   last-defined fix in its end_of_step below, after every body has
  //   committed its end-of-step geometry, so that body-body contact
  //   forces can be applied equal-and-opposite to both bodies

  // error if body now extends beyond a periodic boundary,
  //   b/c body coords are not wrapped across periodic boundaries
  // body is allowed to exit thru non-periodic boundaries
  // test the true body extent, not the eps-inflated bbox

  int outflag = 0;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  int *bflag = domain->bflag;
  double eps = bboxeps;

  // warn once per run if body is entirely outside the simulation box,
  //   b/c it no longer interacts with any particles

  if (!warnexit) {
    if (bbodyhi[0] < boxlo[0] || bbodylo[0] > boxhi[0] ||
        bbodyhi[1] < boxlo[1] || bbodylo[1] > boxhi[1] ||
        (dim == 3 &&
         (bbodyhi[2] < boxlo[2] || bbodylo[2] > boxhi[2]))) {
      warnexit = 1;
      if (comm->me == 0)
        error->warning(FLERR,"Fix rigid body has exited the simulation box "
                       "and no longer interacts with particles");
    }
  }

  if (bflag[0] == PERIODIC && bbodylo[0]+eps < boxlo[0]) outflag = 1;
  if (bflag[1] == PERIODIC && bbodyhi[0]-eps > boxhi[0]) outflag = 1;
  if (bflag[2] == PERIODIC && bbodylo[1]+eps < boxlo[1]) outflag = 1;
  if (bflag[3] == PERIODIC && bbodyhi[1]-eps > boxhi[1]) outflag = 1;
  if (dim == 3) {
    if (bflag[4] == PERIODIC && bbodylo[2]+eps < boxlo[2]) outflag = 1;
    if (bflag[5] == PERIODIC && bbodyhi[2]-eps > boxhi[2]) outflag = 1;
  }

  if (outflag)
    error->all(FLERR,"Fix rigid body moved beyond a periodic boundary");

  // the last-defined rigid fix coordinates the all-body end-of-step
  //   work, after every body has moved to its end-of-step position;
  //   end_of_step runs fixes in definition order, so when the last fix
  //   runs, all earlier bodies are already moved:
  // (1) push-off forces for all bodies; body-body contact forces are
  //     applied equal-and-opposite to both bodies of a contact, so
  //     body-body interactions conserve momentum; reactions are applied
  //     even to bodies without the push keyword
  // (2) re-map body surfs to grid cells: cut/split cells and
  //     INSIDE/OUTSIDE typing from the new body positions; if every
  //     body is incremental, attempt the cheap incremental re-cut of
  //     only the affected cells, else do an exact full grid re-map;
  //     the fallback decision is per-proc but a full re-map is
  //     collective, so all procs must agree via Allreduce
  // (3) remove particles inside any body in one fused pass over
  //     particles, with split-cell reassignment only after a full
  //     re-map; no reduction here, deletion counts stay per-proc and
  //     are reduced lazily by compute_scalar()

  if (flist[nb-1] == this) {

    int anypush = 0;
    for (int m = 0; m < nb; m++) {
      FixRigid *f = flist[m];
      f->fpush[0] = f->fpush[1] = f->fpush[2] = 0.0;
      f->tqpush[0] = f->tqpush[1] = f->tqpush[2] = 0.0;
      if (f->pushflag) anypush = 1;
    }
    if (anypush) {
      for (int m = 0; m < nb; m++)
        if (flist[m]->pushflag) flist[m]->push_off();

      // for distributed surfs the static-contact contributions are
      //   disjoint per-proc partial sums (each proc handles the static
      //   surfs it owns): merge with one Allreduce for all bodies
      // for non-distributed surfs every proc computed identical totals

      if (surf->distributed) {
        for (int m = 0; m < nb; m++) {
          FixRigid *f = flist[m];
          ftbuf_mine[6*m]   = f->fpush[0];
          ftbuf_mine[6*m+1] = f->fpush[1];
          ftbuf_mine[6*m+2] = f->fpush[2];
          ftbuf_mine[6*m+3] = f->tqpush[0];
          ftbuf_mine[6*m+4] = f->tqpush[1];
          ftbuf_mine[6*m+5] = f->tqpush[2];
        }
        MPI_Allreduce(ftbuf_mine,ftbuf_all,6*nb,MPI_DOUBLE,MPI_SUM,world);
        for (int m = 0; m < nb; m++) {
          FixRigid *f = flist[m];
          f->fpush[0] = ftbuf_all[6*m];
          f->fpush[1] = ftbuf_all[6*m+1];
          f->fpush[2] = ftbuf_all[6*m+2];
          f->tqpush[0] = ftbuf_all[6*m+3];
          f->tqpush[1] = ftbuf_all[6*m+4];
          f->tqpush[2] = ftbuf_all[6*m+5];
        }
      }

      for (int m = 0; m < nb; m++) {
        FixRigid *f = flist[m];
        f->fcm[0] += f->fpush[0];
        f->fcm[1] += f->fpush[1];
        f->fcm[2] += f->fpush[2];
        f->torque[0] += f->tqpush[0];
        f->torque[1] += f->tqpush[1];
        f->torque[2] += f->tqpush[2];
      }
    }

    // second half kick of velocity Verlet for every body, now that
    //   its end-of-step force and torque are complete: vcm/angmom/omega
    //   become the velocities at the end of the step, synchronized
    //   with xcm/quat, as reported by the fix and written to outfile

    for (int m = 0; m < nb; m++) flist[m]->final_kick();

    int all_incremental = 1;
    for (int m = 0; m < nb; m++)
      if (flist[m]->remapmode != INCREMENTAL) all_incremental = 0;

    // incremental_recut() returns a reason code > 0 if a full re-map is
    //   required; all procs must agree, so reduce the max
    // warn once per run when the fallback occurs, since a fallback on
    //   every step silently costs as much as remap cutcell

    int fallback = 1;
    if (all_incremental) {
      int fallmine = incremental_recut();
      MPI_Allreduce(&fallmine,&fallback,1,MPI_INT,MPI_MAX,world);
      if (fallback && !warnfallback) {
        warnfallback = 1;
        if (comm->me == 0) {
          const char *why;
          if (fallback == FALLBACK_SPLIT)
            why = "a split cell is in the re-cut region";
          else if (fallback == FALLBACK_SURFMAX)
            why = "a cell would exceed global surfmax";
          else why = "no previous body position is known";
          char str[256];
          snprintf(str,sizeof(str),"Fix rigid incremental remap fell back "
                   "to a full grid re-map because %s",why);
          error->warning(FLERR,str);
        }
      }
    }
    if (fallback) grid_rebuild();

    if (particle->exist) remove_inside_all(fallback);

    // advance each incremental body's previous-region bookkeeping

    for (int m = 0; m < nb; m++) {
      FixRigid *f = flist[m];
      if (f->remapmode != INCREMENTAL) continue;
      for (int j = 0; j < 3; j++) {
        f->pbodylo[j] = f->bbodylo[j];
        f->pbodyhi[j] = f->bbodyhi[j];
      }
      f->pbodyflag = 1;
    }
  }

  // write body state to output file every outevery steps
  // file is compatible with the infile option for run continuation

  if (outfile && update->ntimestep % outevery == 0) write_outfile();
}

/* ----------------------------------------------------------------------
   write current rigid body attributes to output file
   format matches what the infile option reads
   moi is written in the space frame for the current body orientation
------------------------------------------------------------------------- */

void FixRigid::write_outfile()
{
  if (comm->me) return;

  // reconstruct space-frame moi from principal moments and current axes
  // I_space = sum over K of inertia[K] e_K outer-product e_K

  double ispace[6];
  ispace[0] = inertia[0]*ex_space[0]*ex_space[0] +
    inertia[1]*ey_space[0]*ey_space[0] + inertia[2]*ez_space[0]*ez_space[0];
  ispace[1] = inertia[0]*ex_space[1]*ex_space[1] +
    inertia[1]*ey_space[1]*ey_space[1] + inertia[2]*ez_space[1]*ez_space[1];
  ispace[2] = inertia[0]*ex_space[2]*ex_space[2] +
    inertia[1]*ey_space[2]*ey_space[2] + inertia[2]*ez_space[2]*ez_space[2];
  ispace[3] = inertia[0]*ex_space[0]*ex_space[1] +
    inertia[1]*ey_space[0]*ey_space[1] + inertia[2]*ez_space[0]*ez_space[1];
  ispace[4] = inertia[0]*ex_space[0]*ex_space[2] +
    inertia[1]*ey_space[0]*ey_space[2] + inertia[2]*ez_space[0]*ez_space[2];
  ispace[5] = inertia[0]*ex_space[1]*ex_space[2] +
    inertia[1]*ey_space[1]*ey_space[2] + inertia[2]*ez_space[1]*ez_space[2];

  FILE *fp = fopen(outfile,"w");
  if (fp == nullptr) error->one(FLERR,"Cannot open fix rigid outfile");

  fprintf(fp,"# rigid body state from fix %s rigid at timestep " BIGINT_FORMAT
          "\n",id,update->ntimestep);
  // fcm/torque are written after the 16 body params, so that a
  //   continuation run resumes with the force and torque which would
  //   have moved the body on the next step

  fprintf(fp,"# mtotal xcm ycm zcm ixx iyy izz ixy ixz iyz "
          "vxcm vycm vzcm lx ly lz fx fy fz tx ty tz\n");
  fprintf(fp,"%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g "
          "%.15g %.15g %.15g %.15g %.15g %.15g "
          "%.15g %.15g %.15g %.15g %.15g %.15g\n",
          massbody,xcm[0],xcm[1],xcm[2],
          ispace[0],ispace[1],ispace[2],ispace[3],ispace[4],ispace[5],
          vcm[0],vcm[1],vcm[2],angmom[0],angmom[1],angmom[2],
          fcm[0],fcm[1],fcm[2],torque[0],torque[1],torque[2]);

  fclose(fp);
}

/* ----------------------------------------------------------------------
   one-time initialization of rigid body attributes from file
------------------------------------------------------------------------- */

void FixRigid::read_infile(char *filename)
{
  // open file and read first non-empty, non-comment line
  // only done by proc 0
  
  if (comm->me == 0) {
    char *start;
    char line[MAXLINE];
    FILE *fp = fopen(filename,"r");
    if (fp == nullptr)
      error->one(FLERR,"Cannot open fix rigid infile");
    while (true) {
      char *eof = fgets(line,MAXLINE,fp);
      if (eof == nullptr) error->one(FLERR,"Unexpected end of fix rigid infile");
      start = &line[strspn(line," \t\n\v\f\r")];
      if (*start != '\0' && *start != '#') break;
    }

    // check that line has correct number of words
    
    // 16 params, optionally followed by the force and torque which
    //   act on the body during the first step of a continuation run

    int nwords = input->count_words(line);
    if (nwords != 16 && nwords != 22)
      error->one(FLERR,"Incorrect rigid body format in fix rigid infile");
    if (nwords == 22) forceinfile = 1;

    // convert each word to a rigid body param
    // totalmass, xcm, moi, vcm, angmom

    massbody = atof(strtok(line," \t\n\r\f"));
    xcm[0] = atof(strtok(NULL," \t\n\r\f"));
    xcm[1] = atof(strtok(NULL," \t\n\r\f"));
    xcm[2] = atof(strtok(NULL," \t\n\r\f"));
    moi[0] = atof(strtok(NULL," \t\n\r\f"));
    moi[1] = atof(strtok(NULL," \t\n\r\f"));
    moi[2] = atof(strtok(NULL," \t\n\r\f"));
    moi[3] = atof(strtok(NULL," \t\n\r\f"));
    moi[4] = atof(strtok(NULL," \t\n\r\f"));
    moi[5] = atof(strtok(NULL," \t\n\r\f"));
    vcm[0] = atof(strtok(NULL," \t\n\r\f"));
    vcm[1] = atof(strtok(NULL," \t\n\r\f"));
    vcm[2] = atof(strtok(NULL," \t\n\r\f"));
    angmom[0] = atof(strtok(NULL," \t\n\r\f"));
    angmom[1] = atof(strtok(NULL," \t\n\r\f"));
    angmom[2] = atof(strtok(NULL," \t\n\r\f"));

    if (forceinfile) {
      for (int j = 0; j < 3; j++)
        fcm_infile[j] = atof(strtok(NULL," \t\n\r\f"));
      for (int j = 0; j < 3; j++)
        torque_infile[j] = atof(strtok(NULL," \t\n\r\f"));
    }

    fclose(fp);
  }

  // broadcast result of file read to all procs
    
  MPI_Bcast(&massbody,1,MPI_DOUBLE,0,world);
  MPI_Bcast(xcm,3,MPI_DOUBLE,0,world);
  MPI_Bcast(moi,6,MPI_DOUBLE,0,world);
  MPI_Bcast(vcm,3,MPI_DOUBLE,0,world);
  MPI_Bcast(angmom,3,MPI_DOUBLE,0,world);
  MPI_Bcast(&forceinfile,1,MPI_INT,0,world);
  MPI_Bcast(fcm_infile,3,MPI_DOUBLE,0,world);
  MPI_Bcast(torque_infile,3,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   build the replicated table of body elements: geometry + attributes
   every proc stores every body element, the authoritative source for
     all body computations; bodies are compact so this is small even
     when the full surf collection is distributed
   non-distributed surfs: filled directly from the local surf arrays,
     which hold all surfs on every proc
   distributed surfs: each proc contributes the body elements it owns,
     gathered to all procs and sorted by surf ID so every proc builds
     the identical table
------------------------------------------------------------------------- */

// per-element record exchanged between procs for distributed surfs

struct BodyElemRecord {
  double pts[9];
  surfint id;
  int type,mask,trans,isc,isr;
};

static bool body_record_cmp(const BodyElemRecord &a, const BodyElemRecord &b)
{
  return a.id < b.id;
}

void FixRigid::gather_body()
{
  int i,j,k;

  // nsurf = # of lines/tris in rigid body

  bigint bnsurf = surf->count_group(igroup);
  if (bnsurf > MAXSMALLINT) error->all(FLERR,"Too many surfs in rigid body");
  nsurf = bnsurf;
  if (nsurf == 0) error->all(FLERR,"Fix rigid body has no surface elements");

  memory->create(bodypt,nsurf,dim,3,"fix_rigid:bodypt");
  memory->create(bodynorm,nsurf,3,"fix_rigid:bodynorm");
  memory->create(sids,nsurf,"fix_rigid:sids");
  memory->create(bodymask,nsurf,"fix_rigid:bodymask");
  memory->create(bodytype,nsurf,"fix_rigid:bodytype");
  memory->create(bodytrans,nsurf,"fix_rigid:bodytrans");
  memory->create(bodyisc,nsurf,"fix_rigid:bodyisc");
  memory->create(bodyisr,nsurf,"fix_rigid:bodyisr");
  memory->create(lblist,nsurf,"fix_rigid:lblist");

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  if (!surf->distributed) {

    // slist = list of local surf indices in the body group

    memory->create(slist,nsurf,"fix_rigid:slist");

    int nlocal = surf->nlocal;
    int n = 0;

    for (i = 0; i < nlocal; i++) {
      int mask = (dim == 2) ? lines[i].mask : tris[i].mask;
      if (!(mask & groupbit)) continue;
      slist[n] = i;
      if (dim == 2) {
        sids[n] = lines[i].id;
        bodymask[n] = lines[i].mask;
        bodytype[n] = lines[i].type;
        bodytrans[n] = lines[i].transparent;
        bodyisc[n] = lines[i].isc;
        bodyisr[n] = lines[i].isr;
        memcpy(bodypt[n][0],lines[i].p1,3*sizeof(double));
        memcpy(bodypt[n][1],lines[i].p2,3*sizeof(double));
        memcpy(bodynorm[n],lines[i].norm,3*sizeof(double));
      } else {
        sids[n] = tris[i].id;
        bodymask[n] = tris[i].mask;
        bodytype[n] = tris[i].type;
        bodytrans[n] = tris[i].transparent;
        bodyisc[n] = tris[i].isc;
        bodyisr[n] = tris[i].isr;
        memcpy(bodypt[n][0],tris[i].p1,3*sizeof(double));
        memcpy(bodypt[n][1],tris[i].p2,3*sizeof(double));
        memcpy(bodypt[n][2],tris[i].p3,3*sizeof(double));
        memcpy(bodynorm[n],tris[i].norm,3*sizeof(double));
      }
      n++;
    }

  } else {

    // pack the body elements this proc owns, gather to all procs,
    //   sort by surf ID so the table is identical everywhere

    slist = NULL;

    Surf::Line *mylines = surf->mylines;
    Surf::Tri *mytris = surf->mytris;
    int nown = surf->nown;

    int nmine = 0;
    for (i = 0; i < nown; i++) {
      int mask = (dim == 2) ? mylines[i].mask : mytris[i].mask;
      if (mask & groupbit) nmine++;
    }

    BodyElemRecord *mine = new BodyElemRecord[MAX(nmine,1)];

    int n = 0;
    for (i = 0; i < nown; i++) {
      int mask = (dim == 2) ? mylines[i].mask : mytris[i].mask;
      if (!(mask & groupbit)) continue;
      BodyElemRecord &r = mine[n++];
      if (dim == 2) {
        r.id = mylines[i].id;
        r.type = mylines[i].type;
        r.mask = mylines[i].mask;
        r.trans = mylines[i].transparent;
        r.isc = mylines[i].isc;
        r.isr = mylines[i].isr;
        memcpy(&r.pts[0],mylines[i].p1,3*sizeof(double));
        memcpy(&r.pts[3],mylines[i].p2,3*sizeof(double));
      } else {
        r.id = mytris[i].id;
        r.type = mytris[i].type;
        r.mask = mytris[i].mask;
        r.trans = mytris[i].transparent;
        r.isc = mytris[i].isc;
        r.isr = mytris[i].isr;
        memcpy(&r.pts[0],mytris[i].p1,3*sizeof(double));
        memcpy(&r.pts[3],mytris[i].p2,3*sizeof(double));
        memcpy(&r.pts[6],mytris[i].p3,3*sizeof(double));
      }
    }

    int nprocs = comm->nprocs;
    int *counts = new int[nprocs];
    int *displs = new int[nprocs];
    int nbytes = nmine * (int) sizeof(BodyElemRecord);
    MPI_Allgather(&nbytes,1,MPI_INT,counts,1,MPI_INT,world);
    displs[0] = 0;
    for (i = 1; i < nprocs; i++) displs[i] = displs[i-1] + counts[i-1];
    bigint btotal = (bigint) displs[nprocs-1] + counts[nprocs-1];
    if (btotal != (bigint) nsurf * sizeof(BodyElemRecord))
      error->all(FLERR,"Fix rigid body element gather is inconsistent");

    BodyElemRecord *all = new BodyElemRecord[nsurf];
    MPI_Allgatherv(mine,nbytes,MPI_BYTE,all,counts,displs,MPI_BYTE,world);
    std::sort(all,all+nsurf,body_record_cmp);

    double d12[3],d13[3];
    double z[3] = {0.0,0.0,1.0};

    for (i = 0; i < nsurf; i++) {
      BodyElemRecord &r = all[i];
      sids[i] = r.id;
      bodymask[i] = r.mask;
      bodytype[i] = r.type;
      bodytrans[i] = r.trans;
      bodyisc[i] = r.isc;
      bodyisr[i] = r.isr;
      memcpy(bodypt[i][0],&r.pts[0],3*sizeof(double));
      memcpy(bodypt[i][1],&r.pts[3],3*sizeof(double));
      if (dim == 3) memcpy(bodypt[i][2],&r.pts[6],3*sizeof(double));

      // recompute the outward normal the same way Surf does

      if (dim == 2) {
        MathExtra::sub3(bodypt[i][1],bodypt[i][0],d12);
        MathExtra::cross3(z,d12,bodynorm[i]);
        MathExtra::norm3(bodynorm[i]);
        bodynorm[i][2] = 0.0;
      } else {
        MathExtra::sub3(bodypt[i][1],bodypt[i][0],d12);
        MathExtra::sub3(bodypt[i][2],bodypt[i][0],d13);
        MathExtra::cross3(d12,d13,bodynorm[i]);
        MathExtra::norm3(bodynorm[i]);
      }
    }

    delete [] mine;
    delete [] all;
    delete [] counts;
    delete [] displs;
  }

  // idmap = global surf ID -> body element index

  idmap.clear();
  for (i = 0; i < nsurf; i++) idmap[sids[i]] = i;

  // lblist = local surf index of each body element on this proc
  // for distributed surfs, ensure_local_copies() fills lblist and the
  //   list of all local copies at setup and after every surf change
  // olist = owned-array index of the body elements this proc owns

  ncopy = maxcopy = 0;
  copy_index = copy_elem = NULL;

  for (i = 0; i < nsurf; i++) lblist[i] = -1;
  int nslocal = surf->nlocal;
  for (i = 0; i < nslocal; i++) {
    surfint id = (dim == 2) ? lines[i].id : tris[i].id;
    k = body_elem(id);
    if (k >= 0) lblist[k] = i;
  }

  nolist = 0;
  olist_own = olist_elem = NULL;
  if (surf->distributed) {
    Surf::Line *mylines = surf->mylines;
    Surf::Tri *mytris = surf->mytris;
    int nown = surf->nown;
    memory->create(olist_own,nsurf,"fix_rigid:olist_own");
    memory->create(olist_elem,nsurf,"fix_rigid:olist_elem");
    for (i = 0; i < nown; i++) {
      surfint id = (dim == 2) ? mylines[i].id : mytris[i].id;
      k = body_elem(id);
      if (k < 0) continue;
      olist_own[nolist] = i;
      olist_elem[nolist] = k;
      nolist++;
    }
  }
}

/* ----------------------------------------------------------------------
   distributed surfs: insure this proc's local (non-ghost) surf arrays
     contain a copy of every body element, at its current position
   the swept collision lists and the particle mover reference body
     surfs by local index, and a fast body can sweep into cells on a
     proc whose local arrays do not yet hold its surfs
   copies must be in the local range: owned cells may only reference
     local surfs (Surf::compress_explicit relies on it)
   ghost surfs follow the local range in the same array, so appending a
     local copy while ghosts exist requires re-packing the ghosts and
     re-indexing the csurfs lists of ghost cells; a ghost copy of a
     promoted element is dropped, its references map to the local copy
   the surf hash is empty outside of Grid::acquire_ghosts(), so it
     needs no maintenance here
   called at setup and from grid_changed() after any grid/surf change;
     refreshes lblist and the list of all local copies
   caller is responsible for update->build_rigidmap() and, if surfs
     were appended, surfs_changed()
------------------------------------------------------------------------- */

void FixRigid::ensure_local_copies()
{
  int i,j,k,m;

  if (!surf->distributed) return;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nslocal = surf->nlocal;
  int nsghost = surf->nghost;

  // every copy of a body element in the local range

  for (k = 0; k < nsurf; k++) lblist[k] = -1;
  ncopy = 0;

  for (i = 0; i < nslocal; i++) {
    surfint id = (dim == 2) ? lines[i].id : tris[i].id;
    k = body_elem(id);
    if (k < 0) continue;
    if (lblist[k] < 0) lblist[k] = i;
    if (ncopy == maxcopy) {
      maxcopy += DELTA_MODIFY;
      memory->grow(copy_index,maxcopy,"fix_rigid:copy_index");
      memory->grow(copy_elem,maxcopy,"fix_rigid:copy_elem");
    }
    copy_index[ncopy] = i;
    copy_elem[ncopy] = k;
    ncopy++;
  }

  int nmissing = 0;
  for (k = 0; k < nsurf; k++)
    if (lblist[k] < 0) nmissing++;
  if (!nmissing) return;

  // save the ghost entries, then truncate the ghost range

  Surf::Line *glines = NULL;
  Surf::Tri *gtris = NULL;
  int *gmap = new int[MAX(nsghost,1)];

  if (nsghost) {
    if (dim == 2) {
      glines = new Surf::Line[nsghost];
      memcpy(glines,&lines[nslocal],nsghost*sizeof(Surf::Line));
    } else {
      gtris = new Surf::Tri[nsghost];
      memcpy(gtris,&tris[nslocal],nsghost*sizeof(Surf::Tri));
    }
  }
  surf->remove_ghosts();

  // append a local copy of each missing element from the body table

  for (k = 0; k < nsurf; k++) {
    if (lblist[k] >= 0) continue;
    if (dim == 2) {
      Surf::Line line;
      memset(&line,0,sizeof(Surf::Line));
      line.id = sids[k];
      line.type = bodytype[k];
      line.mask = bodymask[k];
      line.transparent = bodytrans[k];
      line.isc = bodyisc[k];
      line.isr = bodyisr[k];
      memcpy(line.p1,bodypt[k][0],3*sizeof(double));
      memcpy(line.p2,bodypt[k][1],3*sizeof(double));
      memcpy(line.norm,bodynorm[k],3*sizeof(double));
      surf->add_line_copy(1,&line);
    } else {
      Surf::Tri tri;
      memset(&tri,0,sizeof(Surf::Tri));
      tri.id = sids[k];
      tri.type = bodytype[k];
      tri.mask = bodymask[k];
      tri.transparent = bodytrans[k];
      tri.isc = bodyisc[k];
      tri.isr = bodyisr[k];
      memcpy(tri.p1,bodypt[k][0],3*sizeof(double));
      memcpy(tri.p2,bodypt[k][1],3*sizeof(double));
      memcpy(tri.p3,bodypt[k][2],3*sizeof(double));
      memcpy(tri.norm,bodynorm[k],3*sizeof(double));
      surf->add_tri_copy(1,&tri);
    }
    lblist[k] = surf->nlocal - 1;
    if (ncopy == maxcopy) {
      maxcopy += DELTA_MODIFY;
      memory->grow(copy_index,maxcopy,"fix_rigid:copy_index");
      memory->grow(copy_elem,maxcopy,"fix_rigid:copy_elem");
    }
    copy_index[ncopy] = lblist[k];
    copy_elem[ncopy] = k;
    ncopy++;
  }

  // re-append the saved ghosts after the enlarged local range
  // gmap = new index of each old ghost, a promoted element's ghost copy
  //   maps to its new local copy

  for (m = 0; m < nsghost; m++) {
    surfint id = (dim == 2) ? glines[m].id : gtris[m].id;
    k = body_elem(id);
    if (k >= 0) {
      gmap[m] = lblist[k];
      continue;
    }
    if (dim == 2) surf->add_line_copy(0,&glines[m]);
    else surf->add_tri_copy(0,&gtris[m]);
    gmap[m] = surf->nlocal + surf->nghost - 1;
  }

  // re-index ghost-range entries in the csurfs lists of ghost cells
  // sub cells share the list of their split cell, so visit each once

  if (nsghost) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    int ngtotal = grid->nlocal + grid->nghost;

    for (int icell = nglocal; icell < ngtotal; icell++) {
      if (cells[icell].nsplit <= 0) continue;
      if (cells[icell].nsurf <= 0) continue;
      surfint *csurfs = cells[icell].csurfs;
      int n = cells[icell].nsurf;
      for (j = 0; j < n; j++)
        if (csurfs[j] >= nslocal) csurfs[j] = gmap[csurfs[j]-nslocal];
    }
  }

  delete [] glines;
  delete [] gtris;
  delete [] gmap;
}

/* ----------------------------------------------------------------------
   notify per-surf computes that the local surf arrays changed, so they
     re-size per-surf storage (e.g. ComputeSurf normflux) and refresh
     cached surf pointers
   same action Grid::notify_changed() takes for computes
------------------------------------------------------------------------- */

void FixRigid::surfs_changed()
{
  Compute **compute = modify->compute;
  for (int i = 0; i < modify->ncompute; i++)
    if (compute[i]->per_surf_flag) compute[i]->reallocate();
}

/* ----------------------------------------------------------------------
   write the current replicated body geometry (bodypt/bodynorm) into
     the Surf storage the particle mover and cut pipeline read:
   non-distributed: the local copies every proc stores (via slist)
   distributed: every local copy on this proc (via the copy list) and
     the owned copies (via olist), so a later re-map redistributes
     current coords
------------------------------------------------------------------------- */

void FixRigid::update_surf_copies()
{
  int i,index;

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  if (!surf->distributed) {
    for (i = 0; i < nsurf; i++) {
      index = slist[i];
      if (dim == 2) {
        memcpy(lines[index].p1,bodypt[i][0],3*sizeof(double));
        memcpy(lines[index].p2,bodypt[i][1],3*sizeof(double));
        memcpy(lines[index].norm,bodynorm[i],3*sizeof(double));
      } else {
        memcpy(tris[index].p1,bodypt[i][0],3*sizeof(double));
        memcpy(tris[index].p2,bodypt[i][1],3*sizeof(double));
        memcpy(tris[index].p3,bodypt[i][2],3*sizeof(double));
        memcpy(tris[index].norm,bodynorm[i],3*sizeof(double));
      }
    }
    return;
  }

  for (int m = 0; m < ncopy; m++) {
    index = copy_index[m];
    i = copy_elem[m];
    if (dim == 2) {
      memcpy(lines[index].p1,bodypt[i][0],3*sizeof(double));
      memcpy(lines[index].p2,bodypt[i][1],3*sizeof(double));
      memcpy(lines[index].norm,bodynorm[i],3*sizeof(double));
    } else {
      memcpy(tris[index].p1,bodypt[i][0],3*sizeof(double));
      memcpy(tris[index].p2,bodypt[i][1],3*sizeof(double));
      memcpy(tris[index].p3,bodypt[i][2],3*sizeof(double));
      memcpy(tris[index].norm,bodynorm[i],3*sizeof(double));
    }
  }

  Surf::Line *mylines = surf->mylines;
  Surf::Tri *mytris = surf->mytris;

  for (int m = 0; m < nolist; m++) {
    int iown = olist_own[m];
    i = olist_elem[m];
    if (dim == 2) {
      memcpy(mylines[iown].p1,bodypt[i][0],3*sizeof(double));
      memcpy(mylines[iown].p2,bodypt[i][1],3*sizeof(double));
      memcpy(mylines[iown].norm,bodynorm[i],3*sizeof(double));
    } else {
      memcpy(mytris[iown].p1,bodypt[i][0],3*sizeof(double));
      memcpy(mytris[iown].p2,bodypt[i][1],3*sizeof(double));
      memcpy(mytris[iown].p3,bodypt[i][2],3*sizeof(double));
      memcpy(mytris[iown].norm,bodynorm[i],3*sizeof(double));
    }
  }
}

/* ----------------------------------------------------------------------
   one-time initialization of rigid body attributes
------------------------------------------------------------------------- */

void FixRigid::setup_body()
{
  // build the replicated body element table: geometry + attributes

  gather_body();

  // insure body surfs form a closed (watertight) object
  //   which encloses a non-zero area or volume

  check_watertight();
  check_enclosed();

  // tensor = inertia tensor in space frame
		    
  double tensor[3][3],evectors[3][3];

  tensor[0][0] = moi[0];
  tensor[1][1] = moi[1];
  tensor[2][2] = moi[2];
  tensor[1][2] = tensor[2][1] = moi[5];
  tensor[0][2] = tensor[2][0] = moi[4];
  tensor[0][1] = tensor[1][0] = moi[3];

  // diagonalize the inertia tensor to create body frame
  
  int ierror = MathEigen::jacobi3(tensor,inertia,evectors,1);
  if (ierror) error->all(FLERR,"Insufficient Jacobi rotations for rigid body");

  ex_space[0] = evectors[0][0];
  ex_space[1] = evectors[1][0];
  ex_space[2] = evectors[2][0];
  ey_space[0] = evectors[0][1];
  ey_space[1] = evectors[1][1];
  ey_space[2] = evectors[2][1];
  ez_space[0] = evectors[0][2];
  ez_space[1] = evectors[1][2];
  ez_space[2] = evectors[2][2];

  // for 2d, insure the principal axis aligned with z is in the 3rd slot
  // the z axis is a principal axis b/c ixz = iyz = 0 is enforced for 2d

  if (dim == 2) {
    if (fabs(ez_space[2]) < 1.0-EPSILON) {
      if (fabs(ey_space[2]) > 1.0-EPSILON) {
        std::swap(inertia[1],inertia[2]);
        std::swap(ey_space[0],ez_space[0]);
        std::swap(ey_space[1],ez_space[1]);
        std::swap(ey_space[2],ez_space[2]);
      } else if (fabs(ex_space[2]) > 1.0-EPSILON) {
        std::swap(inertia[0],inertia[2]);
        std::swap(ex_space[0],ez_space[0]);
        std::swap(ex_space[1],ez_space[1]);
        std::swap(ex_space[2],ez_space[2]);
      } else
        error->all(FLERR,"Fix rigid 2d body inertia tensor has "
                   "no principal axis along z");
    }
  }

  // if any principal moment < scaled EPSILON, set to 0.0

  double max;
  max = MAX(inertia[0],inertia[1]);
  max = MAX(max,inertia[2]);

  if (inertia[0] < EPSILON*max) inertia[0] = 0.0;
  if (inertia[1] < EPSILON*max) inertia[1] = 0.0;
  if (inertia[2] < EPSILON*max) inertia[2] = 0.0;

  // validity checks on principal moments of inertia
  // for 2d only the moment about the z axis matters
  // for 3d all must be positive and satisfy the triangle inequality,
  //   else the moi settings are not those of a physical rigid body
  // jacobi3() sorted the moments in increasing order

  if (dim == 2) {
    if (inertia[2] <= 0.0)
      error->all(FLERR,
                 "Fix rigid moment of inertia about z axis must be positive");
  } else {
    if (inertia[0] <= 0.0 || inertia[1] <= 0.0 || inertia[2] <= 0.0)
      error->all(FLERR,
                 "Fix rigid principal moments of inertia must be positive");
    if (inertia[0] + inertia[1] < (1.0-EPSILON)*inertia[2])
      error->all(FLERR,"Fix rigid moments of inertia do not satisfy "
                 "the triangle inequality");
  }

  // enforce 3 evectors as a right-handed coordinate system
  // flip 3rd vector if needed

  double cross[3];
  MathExtra::cross3(ex_space,ey_space,cross);
  if (MathExtra::dot3(cross,ez_space) < 0.0) MathExtra::negate3(ez_space);

  // create initial quaternion
  
  MathExtra::exyz_to_q(ex_space,ey_space,ez_space,quat);

  set_recoil();

  // set displacement for each end/corner point in each line/tri
  // delta = vector from COM to end/corner point in space frame
  // displace = delta rotated to be in basis of principal axes, i.e. in body frame
  // corner pts come from the replicated body geometry built by gather_body()

  double delta[3];

  memory->create(displace,nsurf,dim,3,"fix_rigid:displace");

  for (int i = 0; i < nsurf; i++)
    for (int j = 0; j < dim; j++) {
      delta[0] = bodypt[i][j][0] - xcm[0];
      delta[1] = bodypt[i][j][1] - xcm[1];
      if (dim == 3) delta[2] = bodypt[i][j][2] - xcm[2];
      else delta[2] = 0.0;
      MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
                                  delta,&displace[i][j][0]);
    }

  // initial omega, consistent with initial angmom

  MathExtra::angmom_to_omega(angmom,ex_space,ey_space,ez_space,inertia,omega);

  // rmaxbody = max distance of any body corner pt from the COM

  rmaxbody = 0.0;
  for (int i = 0; i < nsurf; i++)
    for (int j = 0; j < dim; j++)
      rmaxbody = MAX(rmaxbody,MathExtra::len3(&displace[i][j][0]));

  // per-element bbox arrays for swept collision assignment
  //   and push-off candidate pruning

  memory->create(elemlo,nsurf,3,"fix_rigid:elemlo");
  memory->create(elemhi,nsurf,3,"fix_rigid:elemhi");

  // zero body force/torque in case accessed via compute_vector() on step 0

  fcm[0] = fcm[1] = fcm[2] = 0.0;
  torque[0] = torque[1] = torque[2] = 0.0;
  fpush[0] = fpush[1] = fpush[2] = 0.0;
}

/* ----------------------------------------------------------------------
   bin static surfs for push-off candidate pruning
   built once per run in setup(); static surfs never move during a run
   each static surf is added to every bin its bbox overlaps (CSR layout);
     a query gathers the bins overlapping the pushcutoff-inflated body
     bbox and dedups multi-bin surfs with a visit stamp
   bin edge lengths are at least pushcutoff, at most 64 bins per dim
------------------------------------------------------------------------- */

void FixRigid::push_bins()
{
  int i,k,m,ibx,iby,ibz;
  int blo[3],bhi[3];

  // non-distributed: bin the static surfs of the local arrays, which
  //   hold all surfs on every proc (identical bins everywhere)
  // distributed: bin the static surfs this proc OWNS; each surf is
  //   binned on exactly one proc, so per-proc push contributions are
  //   disjoint partial sums
  // static = surf not in any rigid body

  int distributed = surf->distributed;

  Surf::Line *lines;
  Surf::Tri *tris;
  int nslocal;
  if (!distributed) {
    lines = surf->lines;
    tris = surf->tris;
    nslocal = surf->nlocal;
  } else {
    lines = surf->mylines;
    tris = surf->mytris;
    nslocal = surf->nown;
  }

  FixRigid **flist = update->fixrigidlist;
  int nb = update->nfixrigid;
  int *rigidmap = update->rigidmap;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;

  for (k = 0; k < 3; k++) {
    pushbinlo[k] = boxlo[k];
    double len = boxhi[k] - boxlo[k];
    int n = (int) (len/pushcutoff);
    n = MAX(n,1);
    n = MIN(n,64);
    if (dim == 2 && k == 2) n = 1;
    pushnbin[k] = n;
    pushbininv[k] = n/len;
  }
  int nbins = pushnbin[0]*pushnbin[1]*pushnbin[2];

  memory->destroy(pushbinstart);
  memory->destroy(pushbinlist);
  memory->destroy(pushstamp);
  memory->create(pushbinstart,nbins+1,"fix_rigid:pushbinstart");
  memory->create(pushstamp,nslocal,"fix_rigid:pushstamp");
  for (i = 0; i < nslocal; i++) pushstamp[i] = 0;
  pushstampcur = 0;

  // two passes: count entries per bin, then fill

  double slo[3],shi[3];

  for (int pass = 0; pass < 2; pass++) {
    if (pass == 0)
      for (i = 0; i <= nbins; i++) pushbinstart[i] = 0;

    for (m = 0; m < nslocal; m++) {

      // skip surfs belonging to any rigid body

      if (!distributed) {
        if (rigidmap[m] >= 0) continue;
      } else {
        surfint id = (dim == 2) ? lines[m].id : tris[m].id;
        int inbody = 0;
        for (k = 0; k < nb; k++)
          if (flist[k]->body_elem(id) >= 0) { inbody = 1; break; }
        if (inbody) continue;
      }

      if (dim == 2) {
        for (k = 0; k < 2; k++) {
          slo[k] = MIN(lines[m].p1[k],lines[m].p2[k]);
          shi[k] = MAX(lines[m].p1[k],lines[m].p2[k]);
        }
        slo[2] = shi[2] = 0.0;
      } else {
        for (k = 0; k < 3; k++) {
          slo[k] = MIN(tris[m].p1[k],MIN(tris[m].p2[k],tris[m].p3[k]));
          shi[k] = MAX(tris[m].p1[k],MAX(tris[m].p2[k],tris[m].p3[k]));
        }
      }

      for (k = 0; k < 3; k++) {
        blo[k] = (int) ((slo[k]-pushbinlo[k]) * pushbininv[k]);
        bhi[k] = (int) ((shi[k]-pushbinlo[k]) * pushbininv[k]);
        blo[k] = MAX(0,MIN(blo[k],pushnbin[k]-1));
        bhi[k] = MAX(0,MIN(bhi[k],pushnbin[k]-1));
      }

      for (ibz = blo[2]; ibz <= bhi[2]; ibz++)
        for (iby = blo[1]; iby <= bhi[1]; iby++)
          for (ibx = blo[0]; ibx <= bhi[0]; ibx++) {
            int ibin = (ibz*pushnbin[1] + iby)*pushnbin[0] + ibx;
            if (pass == 0) pushbinstart[ibin+1]++;
            else pushbinlist[pushbinstart[ibin]++] = m;
          }
    }

    if (pass == 0) {
      for (i = 0; i < nbins; i++) pushbinstart[i+1] += pushbinstart[i];
      memory->create(pushbinlist,pushbinstart[nbins],
                     "fix_rigid:pushbinlist");
    } else {
      // filling advanced the starts by one bin: shift them back
      for (i = nbins; i > 0; i--) pushbinstart[i] = pushbinstart[i-1];
      pushbinstart[0] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   contact forces between all corner pts of this body and one source
     element with corner pts p1,p2 (p3 for 3d) and outward normal norm
   for each body corner pt within pushcutoff of the element, apply a
     repulsive force along the element outward normal, with overlap
     delta = pushcutoff - dist:
     linear spring F = kpush * delta, or
     Hertzian contact F = kpush * delta^3/2 (smooth onset, standard
     model for elastic contact of spherical particulates)
   if gammapush > 0, a dashpot term F -= gammapush * d(delta)/dt is
     added (the DEM spring-dashpot pair), computed from the normal
     approach rate of the corner pt relative to the source surface;
     the total contact force is clamped at zero, so the dashpot never
     produces adhesion as a contact ends
   src = the rigid body the element belongs to, or NULL if static
   if src is set, the reaction force -F is applied to src at the same
     contact point, so body-body contacts conserve momentum exactly
------------------------------------------------------------------------- */

void FixRigid::push_contact(double *p1, double *p2, double *p3,
                            double *norm, FixRigid *src)
{
  int i,j;
  double dsq,d,scale;
  double **pts;
  double fone[3],rdelta[3],tq[3];

  int npoint = dim;     // 2 corner pts per line, 3 per tri
  double cutsq = pushcutoff*pushcutoff;

  for (i = 0; i < nsurf; i++) {
    pts = bodypt[i];

    for (j = 0; j < npoint; j++) {
      if (dim == 2)
        dsq = Geometry::distsq_point_line(pts[j],p1,p2);
      else
        dsq = Geometry::distsq_point_tri(pts[j],p1,p2,p3,norm);
      if (dsq >= cutsq) continue;

      d = sqrt(dsq);
      if (pushstyle == LINEAR) scale = kpush * (pushcutoff-d);
      else scale = kpush * (pushcutoff-d) * sqrt(pushcutoff-d);

      // dashpot: damp by the normal approach rate of the corner pt
      //   relative to the source surface,
      //   which moves if it belongs to another rigid body

      if (gammapush > 0.0) {
        double vpt[3],vsrc[3],rd[3];
        MathExtra::sub3(pts[j],xcm,rd);
        MathExtra::cross3(omega,rd,vpt);
        MathExtra::add3(vcm,vpt,vpt);
        if (src) {
          MathExtra::sub3(pts[j],src->xcm,rd);
          MathExtra::cross3(src->omega,rd,vsrc);
          MathExtra::add3(src->vcm,vsrc,vsrc);
          MathExtra::sub3(vpt,vsrc,vpt);
        }
        scale -= gammapush * MathExtra::dot3(vpt,norm);
        if (scale < 0.0) scale = 0.0;
      }

      fone[0] = scale*norm[0];
      fone[1] = scale*norm[1];
      fone[2] = scale*norm[2];

      fpush[0] += fone[0];
      fpush[1] += fone[1];
      fpush[2] += fone[2];
      MathExtra::sub3(pts[j],xcm,rdelta);
      MathExtra::cross3(rdelta,fone,tq);
      tqpush[0] += tq[0];
      tqpush[1] += tq[1];
      tqpush[2] += tq[2];

      // equal-and-opposite reaction on the source body,
      //   applied at the same contact point

      if (src) {
        src->fpush[0] -= fone[0];
        src->fpush[1] -= fone[1];
        src->fpush[2] -= fone[2];
        MathExtra::sub3(pts[j],src->xcm,rdelta);
        MathExtra::cross3(rdelta,fone,tq);
        src->tqpush[0] -= tq[0];
        src->tqpush[1] -= tq[1];
        src->tqpush[2] -= tq[2];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   push-off forces on the body from too-close static surfs, other
     rigid bodies, and (if pushboundflag) non-periodic box boundaries
   called by the last-defined rigid fix for each body with the push
     keyword, after all bodies have committed end-of-step geometry
   static surf candidates come from the bins built by push_bins();
     other bodies are pruned by a body-body bbox test, then per element
   forces accumulate in fpush/torque; the caller adds fpush into fcm
     for the next step's time integration
   non-distributed surfs: computed identically on every proc, so no
     communication is needed; distributed surfs: per-proc partial sums
     which the caller merges with one Allreduce
   NOTE: a corner pt shared by adjacent body elements contributes once
     per element, and a corner close to several source elements
     interacts with each of them, so kpush is a per-contact stiffness
------------------------------------------------------------------------- */

void FixRigid::push_off()
{
  int i,j,m,e;
  double d,scale;
  double **pts;
  double fone[3],rdelta[3],tq[3];
  int blo[3],bhi[3];

  // static surf sources:
  //   non-distributed: local surf arrays hold all surfs on every proc,
  //     every proc computes the identical full contribution
  //   distributed: each proc's bins hold only the static surfs it OWNS,
  //     so contributions are disjoint partial sums, merged by the
  //     coordinating fix with an Allreduce
  // pair and boundary contributions are identical on every proc, so
  //   for distributed surfs only proc 0 computes them before the merge

  int distributed = surf->distributed;
  Surf::Line *lines;
  Surf::Tri *tris;
  if (!distributed) {
    lines = surf->lines;
    tris = surf->tris;
  } else {
    lines = surf->mylines;
    tris = surf->mytris;
  }

  int npoint = dim;     // 2 corner pts per line, 3 per tri

  // cutlo/cuthi = bbox around body inflated by pushcutoff
  // requires body_bbox() was called for current body position

  double cutlo[3],cuthi[3];
  for (j = 0; j < 3; j++) {
    cutlo[j] = bbodylo[j] - pushcutoff;
    cuthi[j] = bbodyhi[j] + pushcutoff;
  }

  // static surf candidates: bins overlapping the inflated body bbox
  // stamp dedups surfs binned into more than one of the bins

  for (j = 0; j < 3; j++) {
    blo[j] = (int) ((cutlo[j]-pushbinlo[j]) * pushbininv[j]);
    bhi[j] = (int) ((cuthi[j]-pushbinlo[j]) * pushbininv[j]);
    blo[j] = MAX(0,MIN(blo[j],pushnbin[j]-1));
    bhi[j] = MAX(0,MIN(bhi[j],pushnbin[j]-1));
  }

  pushstampcur++;

  for (int ibz = blo[2]; ibz <= bhi[2]; ibz++)
    for (int iby = blo[1]; iby <= bhi[1]; iby++)
      for (int ibx = blo[0]; ibx <= bhi[0]; ibx++) {
        int ibin = (ibz*pushnbin[1] + iby)*pushnbin[0] + ibx;
        for (i = pushbinstart[ibin]; i < pushbinstart[ibin+1]; i++) {
          m = pushbinlist[i];
          if (pushstamp[m] == pushstampcur) continue;
          pushstamp[m] = pushstampcur;

          if (dim == 2) {
            if (MAX(lines[m].p1[0],lines[m].p2[0]) < cutlo[0]) continue;
            if (MIN(lines[m].p1[0],lines[m].p2[0]) > cuthi[0]) continue;
            if (MAX(lines[m].p1[1],lines[m].p2[1]) < cutlo[1]) continue;
            if (MIN(lines[m].p1[1],lines[m].p2[1]) > cuthi[1]) continue;
            push_contact(lines[m].p1,lines[m].p2,NULL,lines[m].norm,NULL);
          } else {
            if (MAX(tris[m].p1[0],MAX(tris[m].p2[0],tris[m].p3[0])) <
                cutlo[0]) continue;
            if (MIN(tris[m].p1[0],MIN(tris[m].p2[0],tris[m].p3[0])) >
                cuthi[0]) continue;
            if (MAX(tris[m].p1[1],MAX(tris[m].p2[1],tris[m].p3[1])) <
                cutlo[1]) continue;
            if (MIN(tris[m].p1[1],MIN(tris[m].p2[1],tris[m].p3[1])) >
                cuthi[1]) continue;
            if (MAX(tris[m].p1[2],MAX(tris[m].p2[2],tris[m].p3[2])) <
                cutlo[2]) continue;
            if (MIN(tris[m].p1[2],MIN(tris[m].p2[2],tris[m].p3[2])) >
                cuthi[2]) continue;
            push_contact(tris[m].p1,tris[m].p2,tris[m].p3,
                         tris[m].norm,NULL);
          }
        }
      }

  // other rigid bodies: body-body bbox prefilter, then per-element
  //   bbox tests using the current-position element boxes set by
  //   body_bbox(0) when each body committed its end-of-step geometry
  // element geometry from the source body's replicated bodypt/bodynorm
  // each contact applies equal-and-opposite forces to both bodies

  if (!distributed || comm->me == 0) {
    FixRigid **flist = update->fixrigidlist;
    int nb = update->nfixrigid;

    for (int mb = 0; mb < nb; mb++) {
      FixRigid *g = flist[mb];
      if (g == this) continue;
      if (!box_overlap(cutlo,cuthi,g->bbodylo,g->bbodyhi)) continue;

      for (e = 0; e < g->nsurf; e++) {
        if (!box_overlap(cutlo,cuthi,g->elemlo[e],g->elemhi[e])) continue;
        if (dim == 2)
          push_contact(g->bodypt[e][0],g->bodypt[e][1],NULL,
                       g->bodynorm[e],g);
        else
          push_contact(g->bodypt[e][0],g->bodypt[e][1],g->bodypt[e][2],
                       g->bodynorm[e],g);
      }
    }
  }

  // spring force from non-periodic simulation box boundaries
  // corner pts from the replicated body geometry

  if (pushboundflag && (!distributed || comm->me == 0)) {
    double *boxlo = domain->boxlo;
    double *boxhi = domain->boxhi;
    int *bflag = domain->bflag;

    int nface = 2*dim;
    double fsign[6] = {1.0,-1.0,1.0,-1.0,1.0,-1.0};

    for (i = 0; i < nsurf; i++) {
      pts = bodypt[i];

      for (j = 0; j < npoint; j++) {
        for (int iface = 0; iface < nface; iface++) {
          if (bflag[iface] == PERIODIC) continue;
          int idim = iface/2;
          if (iface % 2 == 0) d = pts[j][idim] - boxlo[idim];
          else d = boxhi[idim] - pts[j][idim];
          if (d >= pushcutoff) continue;

          if (pushstyle == LINEAR) scale = kpush * (pushcutoff-d);
          else scale = kpush * (pushcutoff-d) * sqrt(pushcutoff-d);

          // dashpot vs the static boundary, face normal = fsign*e_idim

          if (gammapush > 0.0) {
            double vpt[3],rd[3];
            MathExtra::sub3(pts[j],xcm,rd);
            MathExtra::cross3(omega,rd,vpt);
            MathExtra::add3(vcm,vpt,vpt);
            scale -= gammapush * fsign[iface] * vpt[idim];
            if (scale < 0.0) scale = 0.0;
          }

          scale *= fsign[iface];
          fone[0] = fone[1] = fone[2] = 0.0;
          fone[idim] = scale;

          fpush[0] += fone[0];
          fpush[1] += fone[1];
          fpush[2] += fone[2];
          MathExtra::sub3(pts[j],xcm,rdelta);
          MathExtra::cross3(rdelta,fone,tq);
          tqpush[0] += tq[0];
          tqpush[1] += tq[1];
          tqpush[2] += tq[2];
        }
      }
    }
  }

  // fpush/tqpush are merged into fcm/torque by the coordinating fix
  //   after all bodies' push-off forces, including reactions from
  //   other bodies' contacts, have been accumulated
}

/* ----------------------------------------------------------------------
   full re-map of surfs to grid cells
   same sequence of operations as in FixMoveSurf::end_of_step()
   called every step (cutcell mode) or on incremental-mode fallback,
     after body surfs move, to cut/split cells and set INSIDE/OUTSIDE
     cell typing from the new body positions
------------------------------------------------------------------------- */

void FixRigid::grid_rebuild()
{
  // sort particles, grid rebuild requires it

  if (particle->exist) particle->sort();

  // assign split cell particles to parent split cell

  grid->unset_neighbors();
  grid->remove_ghosts();

  if (grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->combine_split_cell_particles(icell,1);
  }

  // assign surfs to grid cells

  grid->clear_surf();
  grid->surf2grid(1,0);

  // re-setup owned and ghost cell info

  grid->setup_owned();
  grid->acquire_ghosts();
  grid->reset_neighbors();
  comm->reset_neighbors();

  // flag cells and corners as OUTSIDE or INSIDE

  grid->set_inout();
  grid->type_check(0);

  // notify all classes that store per-grid data that grid may have changed
  // invokes grid_changed() of every rigid fix, which re-establishes the
  //   local body-surf copies and the per-surf rigidmap for the rebuilt
  //   local/ghost surf arrays, before per-surf computes re-size

  grid->notify_changed();
}

/* ----------------------------------------------------------------------
   add every body's surfs to the collision lists (csurfs) of all grid
     cells they sweep through during this step, so particles anywhere in
     a body's swept path are tested against the moving surfs and
     reflected (rather than overtaken by a fast body and deleted)
   called each start_of_step by the last-defined rigid fix, after the
     end-of-step pose of every body is known; a single pass over the
     grid cells handles all bodies
   for each overlapped cell a merged list = its current csurfs plus all
     swept surfs not already present is installed; the original is
     saved for swept_restore(); split cells share the merged list with
     their sub cells, where particles reside
   cut-cell volumes are not changed: this augments collision lists only
------------------------------------------------------------------------- */

void FixRigid::swept_assign_all()
{
  int i,j,m,icell,isub,ncur,isplit,dup,nmerged;
  surfint *merged,*cur;
  FixRigid *f;

  FixRigid **flist = update->fixrigidlist;
  int nb = update->nfixrigid;

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  int ntotal = grid->nlocal + grid->nghost;

  // per-element swept bounding boxes of every body for this step

  for (m = 0; m < nb; m++) flist[m]->body_bbox(1);

  cpage->reset();
  nmodified = 0;

  // phase 1: gather (cell, swept element) entries per body, visiting
  //   only the candidate cells near each body from the box->cell
  //   index, so cost scales with the bodies' swept regions and not
  //   with the number of cells this proc owns
  // entries for one cell are chained; a per-cell stamp detects the
  //   first touch of a cell this step

  if (ntotal > maxswcell) {
    int oldmax = maxswcell;
    maxswcell = ntotal;
    memory->grow(swstamp,maxswcell,"fix_rigid:swstamp");
    memory->grow(swhead,maxswcell,"fix_rigid:swhead");
    for (i = oldmax; i < maxswcell; i++) swstamp[i] = 0;
  }
  swcur++;
  nswcell = 0;
  nent = 0;

  int ncand,icand;
  int *cand;

  for (m = 0; m < nb; m++) {
    f = flist[m];
    ncand = update->rigid_cell_box(f->bbodylo,f->bbodyhi,&cand);

    for (icand = 0; icand < ncand; icand++) {
      icell = cand[icand];
      if (cells[icell].nsplit <= 0) continue;
      if (cells[icell].nsurf < 0) continue;
      if (!box_overlap(cells[icell].lo,cells[icell].hi,
                       f->bbodylo,f->bbodyhi)) continue;

      for (i = 0; i < f->nsurf; i++) {
        if (!box_overlap(cells[icell].lo,cells[icell].hi,
                         f->elemlo[i],f->elemhi[i])) continue;

        if (swstamp[icell] != swcur) {
          swstamp[icell] = swcur;
          swhead[icell] = -1;
          if (nswcell == maxswcells) {
            maxswcells += DELTA_MODIFY;
            memory->grow(swcells,maxswcells,"fix_rigid:swcells");
          }
          swcells[nswcell++] = icell;
        }

        // lblist = local surf index of the element on this proc;
        // for distributed surfs grid_changed() keeps it current

        if (nent == maxent) {
          maxent += DELTA_MODIFY;
          memory->grow(entnext,maxent,"fix_rigid:entnext");
          memory->grow(entelem,maxent,"fix_rigid:entelem");
        }
        entelem[nent] = (surfint) f->lblist[i];
        entnext[nent] = swhead[icell];
        swhead[icell] = nent++;
      }
    }
  }

  // phase 2: for each touched cell, install a merged list =
  //   current csurfs + chained swept elements not already present
  // dedup vs current csurfs skips body surfs the cut pipeline placed
  //   at a body's start-of-step position; chained entries are unique
  //   among themselves (bodies are disjoint, one entry per element)

  for (int ic = 0; ic < nswcell; ic++) {
    icell = swcells[ic];

    ncur = cells[icell].nsurf;
    cur = cells[icell].csurfs;
    merged = cpage->vget();
    for (j = 0; j < ncur; j++) merged[j] = cur[j];
    nmerged = ncur;

    for (int e = swhead[icell]; e >= 0; e = entnext[e]) {
      surfint selem = entelem[e];
      dup = 0;
      for (j = 0; j < ncur; j++)
        if (cur[j] == selem) { dup = 1; break; }
      if (!dup) merged[nmerged++] = selem;
    }

    if (nmerged == ncur) continue;   // all swept surfs already present
    cpage->vgot(nmerged);

    // save cell settings so they can be restored, then override them
    // install the merged list where particles reside: the cell itself
    //   if unsplit, else only its sub cells; the split cell's own list
    //   must keep its original length, since Update::split2d/3d()
    //   index the sinfo csplits array in lockstep with it

    if (nmodified+MAX(cells[icell].nsplit,1) > maxmodified) {
      maxmodified += DELTA_MODIFY;
      memory->grow(modified,maxmodified,"fix_rigid:modified");
      memory->grow(nsurf_saved,maxmodified,"fix_rigid:nsurf_saved");
      csurfs_saved = (surfint **)
        memory->srealloc(csurfs_saved,maxmodified*sizeof(surfint *),
                         "fix_rigid:csurfs_saved");
    }

    if (cells[icell].nsplit == 1) {
      modified[nmodified] = icell;
      nsurf_saved[nmodified] = cells[icell].nsurf;
      csurfs_saved[nmodified] = cells[icell].csurfs;
      nmodified++;
      cells[icell].nsurf = nmerged;
      cells[icell].csurfs = merged;
    } else {
      isplit = cells[icell].isplit;
      for (j = 0; j < cells[icell].nsplit; j++) {
        isub = sinfo[isplit].csubs[j];
        modified[nmodified] = isub;
        nsurf_saved[nmodified] = cells[isub].nsurf;
        csurfs_saved[nmodified] = cells[isub].csurfs;
        nmodified++;
        cells[isub].nsurf = nmerged;
        cells[isub].csurfs = merged;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   restore the csurfs collision lists augmented by swept_assign()
------------------------------------------------------------------------- */

void FixRigid::swept_restore()
{
  Grid::ChildCell *cells = grid->cells;

  for (int m = 0; m < nmodified; m++) {
    int icell = modified[m];
    cells[icell].nsurf = nsurf_saved[m];
    cells[icell].csurfs = csurfs_saved[m];
  }
  nmodified = 0;
}

/* ----------------------------------------------------------------------
   grid cells were rebuilt, adapted, or migrated to other procs
   called via Grid::notify_changed(), after the new owned cells and
     ghost cells (and for distributed surfs, the local/ghost surf
     arrays) are in place, and before per-surf computes re-size
   any merged csurfs lists were discarded by the grid rebuild, and any
     csurfs lists installed by incremental re-cutting were copied into
     grid storage by Grid::compress() or discarded by Grid::clear_surf()
   next re-map re-cuts body surfs into the new grid cells
------------------------------------------------------------------------- */

void FixRigid::grid_changed()
{
  nmodified = 0;
  if (cpage) cpage->reset();
  free_registry();
  update->rigid_bins_clear();

  // distributed surfs: the local surf arrays were rebuilt, so
  //   re-establish this fix's local body-surf copies and the per-surf
  //   rigidmap, which must also span the newly acquired ghost surfs
  // per-surf computes re-size after all fixes are notified

  if (surf->distributed) {
    ensure_local_copies();
    update->build_rigidmap();
  }
}

/* ----------------------------------------------------------------------
   for incremental remap: record cells interior to the body,
     i.e. INSIDE cells with no surfs whose center is within the body
   called before the body surfs move to their end-of-step positions
------------------------------------------------------------------------- */

void FixRigid::record_oldinside()
{
  double ctr[3];
  int icell;

  // bbox around body at its current (pre-move) position

  body_bbox(0);

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  noldinside = 0;

  // candidate cells near the body from the box->cell index,
  //   restricted to owned cells

  int *cand;
  int ncand = update->rigid_cell_box(bbodylo,bbodyhi,&cand);

  for (int ic = 0; ic < ncand; ic++) {
    icell = cand[ic];
    if (icell >= nglocal) continue;
    if (cells[icell].nsplit != 1) continue;
    if (cells[icell].nsurf) continue;
    if (cinfo[icell].type != CELLINSIDE) continue;
    if (!box_overlap(cells[icell].lo,cells[icell].hi,bbodylo,bbodyhi))
      continue;

    ctr[0] = 0.5 * (cells[icell].lo[0] + cells[icell].hi[0]);
    ctr[1] = 0.5 * (cells[icell].lo[1] + cells[icell].hi[1]);
    if (dim == 3) ctr[2] = 0.5 * (cells[icell].lo[2] + cells[icell].hi[2]);
    else ctr[2] = 0.0;
    if (!inside_body(ctr)) continue;

    if (noldinside == maxoldinside) {
      maxoldinside += DELTA_MODIFY;
      memory->grow(oldinside,maxoldinside,"fix_rigid:oldinside");
    }
    oldinside[noldinside++] = icell;
  }
}

/* ----------------------------------------------------------------------
   for incremental remap: re-cut only grid cells near the body
   a cell is re-cut if the set of surfs overlapping it changed,
     or if it is overlapped by a body surf (whose geometry moved)
   candidate surfs for a cell = the static surfs already in its list
     (static surfs never move, so the set overlapping a cell is fixed)
     plus every element of every body, so the cost per cell is
     O(surfs in cell + body surfs) and independent of the total surf
     count; only local surf indices are ever referenced, as required
     for distributed surfs
   cells interior to the body at its old or new position are re-typed
     as INSIDE/OUTSIDE via parity tests, all other cells are untouched
   ghost cell copies of re-cut cells become stale, which is acceptable:
     the ghost cell surf lists the mover consults are re-covered by the
     swept assignment every step, and cell volumes/types of ghost cells
     are not used
   return 0 if done, else a FALLBACK reason code requesting a full grid
     re-map for a structural change:
     a cell would become or stop being a split cell, a cell's surf
     count exceeds maxsurfpercell, or no previous body position is set
------------------------------------------------------------------------- */

int FixRigid::incremental_recut()
{
  int i,n,ncand,icell,nsplitone,xsub,moving,nontrans;
  double vol;
  double xsplit[3],ctr[3],rlo[3],rhi[3];
  double *vols;
  double *clo,*chi;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nglocal = grid->nlocal;
  int maxsurfpercell = grid->maxsurfpercell;
  int *rigidmap = update->rigidmap;

  int ncorner = 4;
  if (dim == 3) ncorner = 8;

  // gather all incremental bodies; every one must have a previous region
  // R = rlo/rhi = union over all incremental bodies of the region each
  //   occupied before and after its move this step

  rlo[0] = rlo[1] = rlo[2] = BIG;
  rhi[0] = rhi[1] = rhi[2] = -BIG;
  int nincr = 0;

  FixRigid **flist = update->fixrigidlist;
  int nb = update->nfixrigid;

  for (int m = 0; m < nb; m++) {
    FixRigid *f = flist[m];
    if (f->remapmode != INCREMENTAL) continue;
    if (!f->pbodyflag) return FALLBACK_NOPREV;
    for (i = 0; i < 3; i++) {
      rlo[i] = MIN(rlo[i],MIN(f->pbodylo[i],f->bbodylo[i]));
      rhi[i] = MAX(rhi[i],MAX(f->pbodyhi[i],f->bbodyhi[i]));
    }
    nincr++;
  }
  if (!nincr) return FALLBACK_NOPREV;

  // collect the owned cells overlapping R from the box->cell index;
  //   the re-cut and re-type passes below iterate only this list

  int *cand;
  int ncells = update->rigid_cell_box(rlo,rhi,&cand);

  nrcand = 0;
  for (int ic = 0; ic < ncells; ic++) {
    icell = cand[ic];
    if (icell >= nglocal) continue;
    if (cells[icell].nsplit <= 0) continue;
    if (!box_overlap(cells[icell].lo,cells[icell].hi,rlo,rhi)) continue;
    if (nrcand == maxrcand) {
      maxrcand += DELTA_MODIFY;
      memory->grow(rcand,maxrcand,"fix_rigid:rcand");
    }
    rcand[nrcand++] = icell;
  }

  // pass 1: re-cut cells in R whose surf overlap changed
  //   or which are overlapped by a moved body surf (from any body)
  // candidate list keeps the cell's static surfs in their current
  //   order, followed by the body elements, so an unchanged cell
  //   yields an identical list and is skipped

  for (int ic = 0; ic < nrcand; ic++) {
    icell = rcand[ic];

    // structural change unsupported: split cells trigger a full re-map

    if (cells[icell].nsplit > 1) return FALLBACK_SPLIT;

    ncand = 0;
    surfint *cur = cells[icell].csurfs;
    for (i = 0; i < cells[icell].nsurf; i++)
      if (rigidmap[cur[i]] < 0) reclist[ncand++] = cur[i];
    for (int m = 0; m < nb; m++) {
      FixRigid *f = flist[m];
      for (i = 0; i < f->nsurf; i++) reclist[ncand++] = f->lblist[i];
    }

    // new list of surfs overlapping this cell

    if (dim == 2)
      n = cut2d->surf2grid_list(cells[icell].id,
                                cells[icell].lo,cells[icell].hi,
                                ncand,reclist,newlist,maxsurfpercell);
    else
      n = cut3d->surf2grid_list(cells[icell].id,
                                cells[icell].lo,cells[icell].hi,
                                ncand,reclist,newlist,maxsurfpercell);
    if (n > maxsurfpercell) return FALLBACK_SURFMAX;

    // order the list by local surf index, as Grid::surf2grid() does
    //   before cutting, so the cut sees surfs in the same order in both
    //   remap modes and lists of unchanged cells compare equal

    std::sort(newlist,newlist+n);

    // skip cell if surf list is unchanged and contains no moving surf
    // a moving surf belongs to any rigid body (via rigidmap)

    moving = 0;
    for (i = 0; i < n; i++)
      if (rigidmap[newlist[i]] >= 0) {
        moving = 1;
        break;
      }

    if (!moving && n == cells[icell].nsurf) {
      if (n == 0) continue;
      if (memcmp(newlist,cells[icell].csurfs,n*sizeof(surfint)) == 0)
        continue;
    }

    clo = cells[icell].lo;
    chi = cells[icell].hi;
    ctr[0] = 0.5 * (clo[0] + chi[0]);
    ctr[1] = 0.5 * (clo[1] + chi[1]);
    if (dim == 3) ctr[2] = 0.5 * (clo[2] + chi[2]);
    else ctr[2] = 0.0;

    if (n == 0) {

      // cell no longer overlaps any surf
      // full flow volume, interior/exterior typing via parity test

      registry_remove(icell);
      cells[icell].nsurf = 0;
      cells[icell].csurfs = NULL;

      if (dim == 3)
        vol = (chi[0]-clo[0]) * (chi[1]-clo[1]) * (chi[2]-clo[2]);
      else vol = (chi[0]-clo[0]) * (chi[1]-clo[1]);
      cinfo[icell].volume = vol;

      if (inside_any_body(ctr)) cinfo[icell].type = CELLINSIDE;
      else cinfo[icell].type = CELLOUTSIDE;
      for (i = 0; i < ncorner; i++)
        cinfo[icell].corner[i] = cinfo[icell].type;

    } else {

      // install new surf list and re-cut the cell

      surfint *list =
        (surfint *) memory->smalloc(n*sizeof(surfint),"fix_rigid:recut");
      memcpy(list,newlist,n*sizeof(surfint));
      registry_replace(icell,list);
      cells[icell].nsurf = n;
      cells[icell].csurfs = list;

      if (dim == 2)
        nsplitone = cut2d->split(cells[icell].id,
                                 cells[icell].lo,cells[icell].hi,
                                 n,list,vols,newmap,
                                 cinfo[icell].corner,xsub,xsplit);
      else
        nsplitone = cut3d->split(cells[icell].id,
                                 cells[icell].lo,cells[icell].hi,
                                 n,list,vols,newmap,
                                 cinfo[icell].corner,xsub,xsplit);

      // cell would become a split cell: fall back to full re-map

      if (nsplitone > 1) return FALLBACK_SPLIT;

      // same as Grid::surf2grid_one(): volume only if corners are known
      // cell is OVERLAP only if it has a non-transparent surf,
      //   else typed via parity test like Grid::set_inout()

      if (cinfo[icell].corner[0] != CELLUNKNOWN)
        cinfo[icell].volume = vols[0];

      nontrans = 0;
      for (i = 0; i < n; i++) {
        int trans = (dim == 2) ? lines[list[i]].transparent :
          tris[list[i]].transparent;
        if (!trans) {
          nontrans = 1;
          break;
        }
      }

      if (nontrans) cinfo[icell].type = CELLOVERLAP;
      else if (inside_any_body(ctr)) cinfo[icell].type = CELLINSIDE;
      else cinfo[icell].type = CELLOUTSIDE;
    }
  }

  // pass 2: cells a body interior moved away from become OUTSIDE
  // process every incremental body's recorded interior cells
  // only cells which are now surf-free and inside no body,
  //   which leaves any static (non-body) INSIDE cells untouched

  for (int mb = 0; mb < nb; mb++) {
    FixRigid *f = flist[mb];
    if (f->remapmode != INCREMENTAL) continue;

    for (int m = 0; m < f->noldinside; m++) {
      icell = f->oldinside[m];
      if (cells[icell].nsurf) continue;

      clo = cells[icell].lo;
      chi = cells[icell].hi;
      ctr[0] = 0.5 * (clo[0] + chi[0]);
      ctr[1] = 0.5 * (clo[1] + chi[1]);
      if (dim == 3) ctr[2] = 0.5 * (clo[2] + chi[2]);
      else ctr[2] = 0.0;
      if (inside_any_body(ctr)) continue;

      cinfo[icell].type = CELLOUTSIDE;
      if (dim == 3)
        cinfo[icell].volume = (chi[0]-clo[0]) * (chi[1]-clo[1]) *
          (chi[2]-clo[2]);
      else cinfo[icell].volume = (chi[0]-clo[0]) * (chi[1]-clo[1]);
      for (i = 0; i < ncorner; i++)
        cinfo[icell].corner[i] = CELLOUTSIDE;
    }
  }

  // pass 3: surf-free cells a body interior moved over become INSIDE
  // catches cells swept over entirely within one step, which never
  //   overlap a body surf at start- or end-of-step positions
  // R covers the swept corridor since it unions old and new positions
  // guard on type != INSIDE leaves static INSIDE cells untouched

  for (int ic = 0; ic < nrcand; ic++) {
    icell = rcand[ic];
    if (cells[icell].nsplit != 1) continue;
    if (cells[icell].nsurf) continue;
    if (cinfo[icell].type == CELLINSIDE) continue;

    clo = cells[icell].lo;
    chi = cells[icell].hi;
    ctr[0] = 0.5 * (clo[0] + chi[0]);
    ctr[1] = 0.5 * (clo[1] + chi[1]);
    if (dim == 3) ctr[2] = 0.5 * (clo[2] + chi[2]);
    else ctr[2] = 0.0;
    if (!inside_any_body(ctr)) continue;

    cinfo[icell].type = CELLINSIDE;
    for (i = 0; i < ncorner; i++)
      cinfo[icell].corner[i] = CELLINSIDE;
  }

  return FALLBACK_NONE;
}

/* ----------------------------------------------------------------------
   return 1 if point x is inside any rigid body, else 0
------------------------------------------------------------------------- */

int FixRigid::inside_any_body(double *x)
{
  for (int m = 0; m < update->nfixrigid; m++)
    if (update->fixrigidlist[m]->inside_body(x)) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   registry of cells whose csurfs lists are allocated by this fix
   grid-owned csurfs lists live in page memory and are never freed
     individually, so lists installed by incremental re-cutting are
     tracked here and freed when replaced or when the grid changes
------------------------------------------------------------------------- */

void FixRigid::registry_replace(int icell, surfint *list)
{
  std::map<int,surfint *>::iterator it = registry.find(icell);
  if (it != registry.end()) {
    memory->sfree(it->second);
    it->second = list;
  } else registry[icell] = list;
}

void FixRigid::registry_remove(int icell)
{
  std::map<int,surfint *>::iterator it = registry.find(icell);
  if (it == registry.end()) return;
  memory->sfree(it->second);
  registry.erase(it);
}

void FixRigid::free_registry()
{
  for (std::map<int,surfint *>::iterator it = registry.begin();
       it != registry.end(); ++it)
    memory->sfree(it->second);
  registry.clear();
}

/* ----------------------------------------------------------------------
   copy every registry list still installed in a live grid cell into
     grid-owned page storage, preserving pointer sharing between a
     split cell and its sub cells, so no cell is left pointing at
     memory this fix is about to free
   used by the destructor: a fix can be unfixed between runs after
     incremental re-cuts installed lists in cells
------------------------------------------------------------------------- */

void FixRigid::copy_registry_to_grid()
{
  if (registry.empty()) return;

  std::map<surfint *,surfint *> replaced;
  Grid::ChildCell *cells = grid->cells;
  int ntotal = grid->nlocal + grid->nghost;

  for (std::map<int,surfint *>::iterator it = registry.begin();
       it != registry.end(); ++it)
    replaced[it->second] = NULL;

  for (int icell = 0; icell < ntotal; icell++) {
    if (cells[icell].nsurf <= 0) continue;
    std::map<surfint *,surfint *>::iterator it =
      replaced.find(cells[icell].csurfs);
    if (it == replaced.end()) continue;
    if (it->second == NULL) {
      surfint *copy = grid->csurfs->get(cells[icell].nsurf);
      memcpy(copy,cells[icell].csurfs,cells[icell].nsurf*sizeof(surfint));
      it->second = copy;
    }
    cells[icell].csurfs = it->second;
  }
}

/* ----------------------------------------------------------------------
   compute per-element bounding boxes (elemlo/elemhi) and the whole-body
     bounding box (bbodylo/bbodyhi)
   sweepflag = 0: boxes bound elements at their current positions
   sweepflag = 1: boxes also bound elements at their end-of-step positions
     (from xcmnew and exyz_space set from quatnew in start_of_step),
     giving the region each element sweeps through during the step
   boxes are inflated by EPSSURF * body extent to avoid round-off misses
------------------------------------------------------------------------- */

void FixRigid::body_bbox(int sweepflag)
{
  int i,j,k;
  double **pts;
  double delta[3],ptnew[3];
  double *lo,*hi;

  int npoint = dim;     // 2 points per line, 3 per tri

  bbodylo[0] = bbodylo[1] = bbodylo[2] = BIG;
  bbodyhi[0] = bbodyhi[1] = bbodyhi[2] = -BIG;

  for (i = 0; i < nsurf; i++) {
    pts = bodypt[i];

    lo = elemlo[i];
    hi = elemhi[i];
    lo[0] = lo[1] = lo[2] = BIG;
    hi[0] = hi[1] = hi[2] = -BIG;

    for (j = 0; j < npoint; j++)
      for (k = 0; k < 3; k++) {
        lo[k] = MIN(lo[k],pts[j][k]);
        hi[k] = MAX(hi[k],pts[j][k]);
      }

    if (sweepflag) {
      for (j = 0; j < npoint; j++) {
        MathExtra::matvec(ex_space,ey_space,ez_space,displace[i][j],delta);
        if (dim == 2) delta[2] = 0.0;
        MathExtra::add3(xcmnew,delta,ptnew);
        for (k = 0; k < 3; k++) {
          lo[k] = MIN(lo[k],ptnew[k]);
          hi[k] = MAX(hi[k],ptnew[k]);
        }
      }
    }

    for (k = 0; k < 3; k++) {
      bbodylo[k] = MIN(bbodylo[k],lo[k]);
      bbodyhi[k] = MAX(bbodyhi[k],hi[k]);
    }
  }

  double eps = EPSSURF * MAX(bbodyhi[0]-bbodylo[0],bbodyhi[1]-bbodylo[1]);
  eps = EPSSURF * MAX(eps/EPSSURF,bbodyhi[2]-bbodylo[2]);
  bboxeps = eps;

  for (i = 0; i < nsurf; i++)
    for (k = 0; k < 3; k++) {
      elemlo[i][k] -= eps;
      elemhi[i][k] += eps;
    }
  for (k = 0; k < 3; k++) {
    bbodylo[k] -= eps;
    bbodyhi[k] += eps;
  }
}

/* ----------------------------------------------------------------------
   determine if point X is inside the closed body via a parity test
   count intersections of segment from X to a point outside the body
     with all body elements: odd = inside, even = outside
   segment direction is oblique to coordinate axes to reduce the chance
     of exactly grazing element edges or vertices
   requires body_bbox() was called to set bbodylo/bbodyhi
------------------------------------------------------------------------- */

int FixRigid::inside_body(double *x)
{
  int hitflag,side;
  double param;
  double xout[3],xc[3];

  double dmax = MAX(bbodyhi[0]-bbodylo[0],bbodyhi[1]-bbodylo[1]);
  dmax = MAX(dmax,bbodyhi[2]-bbodylo[2]);

  xout[0] = bbodyhi[0] + 0.414159*dmax;
  xout[1] = x[1] + 0.271828*dmax;
  if (dim == 3) xout[2] = x[2] + 0.161803*dmax;
  else xout[2] = 0.0;

  int count = 0;
  for (int i = 0; i < nsurf; i++) {
    if (dim == 2)
      hitflag = Geometry::
        line_line_intersect(x,xout,bodypt[i][0],bodypt[i][1],
                            bodynorm[i],xc,param,side);
    else
      hitflag = Geometry::
        line_tri_intersect(x,xout,bodypt[i][0],bodypt[i][1],
                           bodypt[i][2],bodynorm[i],xc,param,side);
    if (hitflag) count++;
  }

  return count % 2;
}

/* ----------------------------------------------------------------------
   remove particles which are inside this body
   also remove all particles in INSIDE cells
   used at setup; each step uses the fused remove_inside_all() instead
   splitflag = 1 if called after a grid rebuild,
     to first reassign particles in split cells to their sub cells
   return # of particles deleted by this proc; counts are summed
     across procs lazily by compute_scalar()
------------------------------------------------------------------------- */

bigint FixRigid::remove_inside_particles(int splitflag)
{
  // reassign particles in split cells to sub cell owner
  // requires sorted particles, done by grid_rebuild()

  if (splitflag && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->assign_split_cell_particles(icell);
  }

  // bbox around body at its current position

  body_bbox(0);

  // flag particles inside the body or in INSIDE cells for deletion

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  int nplocal = particle->nlocal;

  int icell;
  double *x;
  int delflag = 0;

  for (int i = 0; i < nplocal; i++) {
    icell = particles[i].icell;
    if (icell < 0) continue;

    if (cinfo[icell].type == CELLINSIDE) {
      particles[i].icell = -1;
      delflag = 1;
      continue;
    }

    x = particles[i].x;
    if (x[0] < bbodylo[0] || x[0] > bbodyhi[0]) continue;
    if (x[1] < bbodylo[1] || x[1] > bbodyhi[1]) continue;
    if (dim == 3 && (x[2] < bbodylo[2] || x[2] > bbodyhi[2])) continue;

    if (inside_body(x)) {
      particles[i].icell = -1;
      delflag = 1;
    }
  }

  // compress out deleted particles

  int nlocal_old = particle->nlocal;
  if (delflag) particle->compress_rebalance();
  return nlocal_old - particle->nlocal;
}

/* ----------------------------------------------------------------------
   remove particles inside any rigid body, in one pass over particles
   called each step by the last-defined rigid fix, after the grid
     re-map; every body's bbodylo/bbodyhi is current (set when each
     body committed its end-of-step geometry)
   also removes all particles in INSIDE cells
   splitflag = 1 if called after a full grid rebuild,
     to first reassign particles in split cells to their sub cells
   deletions increment the owning body's per-proc ndeleted count with
     no communication; compute_scalar() reduces the counts on demand
------------------------------------------------------------------------- */

void FixRigid::remove_inside_all(int splitflag)
{
  int m;
  double *x;

  FixRigid **flist = update->fixrigidlist;
  int nb = update->nfixrigid;

  // reassign particles in split cells to sub cell owner
  // requires sorted particles, done by grid_rebuild()

  if (splitflag && grid->nsplitlocal) {
    Grid::ChildCell *cells = grid->cells;
    int nglocal = grid->nlocal;
    for (int icell = 0; icell < nglocal; icell++)
      if (cells[icell].nsplit > 1)
        grid->assign_split_cell_particles(icell);
  }

  // flag particles inside any body or in INSIDE cells for deletion
  // attribute each deletion to the body containing the particle;
  //   a particle in an INSIDE cell claimed by no body (e.g. inside
  //   static closed geometry) is attributed to this fix

  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  int nplocal = particle->nlocal;

  int icell;
  int delflag = 0;

  for (int i = 0; i < nplocal; i++) {
    icell = particles[i].icell;
    if (icell < 0) continue;

    x = particles[i].x;
    int inside = (cinfo[icell].type == CELLINSIDE);

    int owner = -1;
    for (m = 0; m < nb; m++) {
      FixRigid *f = flist[m];
      if (x[0] < f->bbodylo[0] || x[0] > f->bbodyhi[0]) continue;
      if (x[1] < f->bbodylo[1] || x[1] > f->bbodyhi[1]) continue;
      if (dim == 3 &&
          (x[2] < f->bbodylo[2] || x[2] > f->bbodyhi[2])) continue;
      if (f->inside_body(x)) {
        owner = m;
        break;
      }
    }

    if (owner < 0 && !inside) continue;

    particles[i].icell = -1;
    delflag = 1;
    if (owner >= 0) {
      flist[owner]->ndeleted++;
      flist[owner]->ndelrun++;
    } else {
      ndeleted++;
      ndelrun++;
    }
  }

  // compress out deleted particles, once for all bodies

  if (delflag) particle->compress_rebalance();

  // warn once per run if any particle was deleted after the setup pass
  // with swept collision coverage a particle in the body's path is
  //   reflected, so this should not happen; if it does, either the body
  //   moves so far in one step that it jumps past a particle, or its
  //   surfs coincide with grid cell boundaries, which makes the cut
  //   cells degenerate until the body moves off the alignment
  // checked once, on the last step of the run, to avoid a collective
  //   on every step

  if (!warndelete && update->ntimestep == update->laststep) {
    bigint mine = 0;
    for (m = 0; m < nb; m++) mine += flist[m]->ndelrun;
    bigint all;
    MPI_Allreduce(&mine,&all,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    for (m = 0; m < nb; m++) flist[m]->warndelete = 1;
    if (all && comm->me == 0) {
      char str[256];
      snprintf(str,sizeof(str),BIGINT_FORMAT " particles were deleted inside "
               "a rigid body during this run.  The body may be moving too "
               "far per timestep, or its surfs may lie exactly on grid cell "
               "boundaries",all);
      error->warning(FLERR,str);
    }
  }
}

/* ----------------------------------------------------------------------
   second half kick of velocity Verlet with the end-of-step force/torque
   omega is recomputed from angmom with the end-of-step axes
------------------------------------------------------------------------- */

void FixRigid::final_kick()
{
  double dt = update->dt;
  double dtfhalf = 0.5 * dt / massbody;
  double dthalf = 0.5 * dt;

  vcm[0] += dtfhalf * (fcm[0] + fext[0]);
  vcm[1] += dtfhalf * (fcm[1] + fext[1]);
  vcm[2] += dtfhalf * (fcm[2] + fext[2]);

  angmom[0] += dthalf * torque[0];
  angmom[1] += dthalf * torque[1];
  angmom[2] += dthalf * torque[2];

  MathExtra::angmom_to_omega(angmom,ex_space,ey_space,ez_space,inertia,omega);

  if (dim == 2) {
    vcm[2] = 0.0;
    angmom[0] = 0.0;
    angmom[1] = 0.0;
    omega[0] = 0.0;
    omega[1] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   set invmass and the space-frame inverse inertia tensor from the
     current principal axes and moments, used by the particle mover to
     correct collisions for the recoil of the finite-mass body
   Iinv = sum over K of (1/inertia[K]) e_K outer-product e_K
   for 2d only rotation about z is possible, so Iinv has only a zz
     component = 1/Izz, which keeps the in-plane response decoupled
------------------------------------------------------------------------- */

void FixRigid::set_recoil()
{
  int i,j,k;
  double *e[3] = {ex_space,ey_space,ez_space};

  invmass = 1.0 / massbody;
  for (k = 0; k < 9; k++) invinertia[k] = 0.0;

  if (dim == 2) {
    double izz = 0.0;
    for (k = 0; k < 3; k++) izz += inertia[k]*e[k][2]*e[k][2];
    invinertia[8] = 1.0 / izz;
  } else {
    for (k = 0; k < 3; k++)
      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
          invinertia[3*i+j] += e[k][i]*e[k][j] / inertia[k];
  }
}

/* ----------------------------------------------------------------------
   check that the body surfs form one or more closed (watertight) objects
   2d: each point must appear exactly as often as the 1st endpoint of a
     line as it does as the 2nd endpoint of a line
   3d: each edge must be traversed the same number of times in each
     direction by the tris that share it
   matching of points is on exact floating point values, the same as
     the watertight checks applied to all surfs by the Surf class
   all procs store all surfs, so the check is identical on every proc
------------------------------------------------------------------------- */

void FixRigid::check_watertight()
{
  int unmatched = 0;

  if (dim == 2) {
    std::map<std::array<double,2>,int> count;
    std::array<double,2> key;

    for (int i = 0; i < nsurf; i++) {
      key[0] = bodypt[i][0][0]; key[1] = bodypt[i][0][1];
      count[key]++;
      key[0] = bodypt[i][1][0]; key[1] = bodypt[i][1][1];
      count[key]--;
    }

    for (std::map<std::array<double,2>,int>::iterator it = count.begin();
         it != count.end(); ++it)
      if (it->second != 0) unmatched++;

  } else {
    std::map<std::array<double,6>,int> count;
    std::array<double,6> key;
    double *pts[4];
    double *a,*b;
    int dir;

    for (int i = 0; i < nsurf; i++) {
      pts[0] = bodypt[i][0]; pts[1] = bodypt[i][1];
      pts[2] = bodypt[i][2]; pts[3] = bodypt[i][0];

      for (int j = 0; j < 3; j++) {
        a = pts[j];
        b = pts[j+1];

        // store edge with endpoints in canonical order
        // count is +1 if traversed in that order, -1 if reversed

        dir = 1;
        if (b[0] < a[0] ||
            (b[0] == a[0] &&
             (b[1] < a[1] || (b[1] == a[1] && b[2] < a[2])))) {
          std::swap(a,b);
          dir = -1;
        }

        key[0] = a[0]; key[1] = a[1]; key[2] = a[2];
        key[3] = b[0]; key[4] = b[1]; key[5] = b[2];
        count[key] += dir;
      }
    }

    for (std::map<std::array<double,6>,int>::iterator it = count.begin();
         it != count.end(); ++it)
      if (it->second != 0) unmatched++;
  }

  if (unmatched) {
    char str[128];
    if (dim == 2)
      sprintf(str,"Fix rigid body is not watertight: "
              "%d unmatched points",unmatched);
    else
      sprintf(str,"Fix rigid body is not watertight: "
              "%d unmatched edges",unmatched);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   check that the body surfs enclose a non-zero area (2d) or volume (3d)
   a zero-thickness body, e.g. a line or tri traversed once in each
     direction, passes the watertight check but has no interior:
     the cut-cell routines cannot mark cells inside/outside it
   measure = signed area via the shoelace sum over the lines, or
     signed volume via the divergence theorem over the tris,
     each computed relative to the centroid of the body points
     so that round-off is set by the body extent, not its position
   only the magnitude is tested, since a watertight body with
     inward normals (a container) is a valid object
   all procs store all surfs, so the check is identical on every proc
------------------------------------------------------------------------- */

void FixRigid::check_enclosed()
{
  int i,j,k;
  double c[3],a[3],b[3],d[3],e[3];

  int npoint = dim;
  double lo[3],hi[3];
  lo[0] = lo[1] = lo[2] = BIG;
  hi[0] = hi[1] = hi[2] = -BIG;
  c[0] = c[1] = c[2] = 0.0;

  for (i = 0; i < nsurf; i++)
    for (j = 0; j < npoint; j++)
      for (k = 0; k < 3; k++) {
        c[k] += bodypt[i][j][k];
        lo[k] = MIN(lo[k],bodypt[i][j][k]);
        hi[k] = MAX(hi[k],bodypt[i][j][k]);
      }
  for (k = 0; k < 3; k++) c[k] /= nsurf*npoint;

  double extent = MAX(hi[0]-lo[0],hi[1]-lo[1]);
  if (dim == 3) extent = MAX(extent,hi[2]-lo[2]);

  double measure = 0.0;

  if (dim == 2) {
    for (i = 0; i < nsurf; i++) {
      MathExtra::sub3(bodypt[i][0],c,a);
      MathExtra::sub3(bodypt[i][1],c,b);
      measure += a[0]*b[1] - a[1]*b[0];
    }
    measure *= 0.5;
  } else {
    for (i = 0; i < nsurf; i++) {
      MathExtra::sub3(bodypt[i][0],c,a);
      MathExtra::sub3(bodypt[i][1],c,b);
      MathExtra::sub3(bodypt[i][2],c,d);
      MathExtra::cross3(b,d,e);
      measure += MathExtra::dot3(a,e);
    }
    measure /= 6.0;
  }

  // scale = extent^dim, a body thinner than EPSENCLOSED of its extent
  //   has no usable interior on any grid that could resolve it

  double scale = extent*extent;
  if (dim == 3) scale *= extent;

  if (fabs(measure) <= EPSENCLOSED*scale) {
    if (dim == 2)
      error->all(FLERR,"Fix rigid body encloses zero area");
    else
      error->all(FLERR,"Fix rigid body encloses zero volume");
  }
}

/* ----------------------------------------------------------------------
   return cummulative count of particles deleted inside the moving body
------------------------------------------------------------------------- */

double FixRigid::compute_scalar()
{
  // ndeleted is a per-proc count, summed across procs on demand
  //   rather than with a collective every step
  // cached per timestep: repeated outputs on one step reduce once
  // like all fix scalar outputs, must be accessed on all procs

  if (ndelvalid != update->ntimestep) {
    MPI_Allreduce(&ndeleted,&ndeleted_all,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
    ndelvalid = update->ntimestep;
  }
  return (double) ndeleted_all;
}

/* ----------------------------------------------------------------------
   return properties of the single rigid body
------------------------------------------------------------------------- */

double FixRigid::compute_vector(int index)
{
  if (index < 3) return xcm[index];
  if (index < 6) return vcm[index-3];
  if (index < 9) return fcm[index-6];
  if (index < 12) return torque[index-9];
  if (index < 15) return omega[index-12];
  if (index < 19) return quat[index-15];
  if (index < 22) return fpush[index-19];

  return 0.0;
}
