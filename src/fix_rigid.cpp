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
#include "random_knuth.h"
#include "random_mars.h"
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
#define BIG 1.0e20
#define DELTA_MODIFY 1024

enum{INT,DOUBLE};                      // several files

enum{OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN};    // same as Update

// cell types, same as Grid
// renamed to avoid clash with surf collision side enum above

enum{CELLUNKNOWN,CELLOUTSIDE,CELLINSIDE,CELLOVERLAP};

enum{PERIODIC,OUTFLOW,REFLECT,SURFACE,AXISYM};  // same as Domain

enum{CUTCELL,INCREMENTAL};          // remap modes

enum{LINEAR,HERTZ};             // push-off force laws

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
  if (domain->axisymmetric)
    error->all(FLERR,"Fix rigid cannot be used with axisymmetric domains");
  if (surf->implicit || surf->distributed)
    error->all(FLERR,"Fix rigid can only be used with explicit non-distributed surf elements");


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

  int iarg = 4;
  if (strcmp(arg[iarg],"body") == 0) {
    if (iarg+22 > narg) error->all(FLERR,"Fix rigid body args not valid");
    massflag = comflag = vcomflag = moiflag = angmomflag = 0;
    int jarg = iarg+1;
    while (jarg < iarg+22) {
      if (strcmp(arg[jarg],"mass") == 0) {
	massflag = 1;
	massbody = input->numeric(FLERR,arg[jarg+1]);
	jarg += 2;
      } else if (strcmp(arg[jarg],"com") == 0) {
	comflag = 1;
	xcm[0] = input->numeric(FLERR,arg[jarg+1]);
	xcm[1] = input->numeric(FLERR,arg[jarg+2]);
	xcm[2] = input->numeric(FLERR,arg[jarg+3]);
	jarg += 4;
      } else if (strcmp(arg[jarg],"moi") == 0) {
	moiflag = 1;
	moi[0] = input->numeric(FLERR,arg[jarg+1]);
	moi[1] = input->numeric(FLERR,arg[jarg+2]);
	moi[2] = input->numeric(FLERR,arg[jarg+3]);
	moi[3] = input->numeric(FLERR,arg[jarg+4]);
	moi[4] = input->numeric(FLERR,arg[jarg+5]);
	moi[5] = input->numeric(FLERR,arg[jarg+6]);
	jarg += 7;
      } else if (strcmp(arg[jarg],"vcom") == 0) {
	vcomflag = 1;
	vcm[0] = input->numeric(FLERR,arg[jarg+1]);
	vcm[1] = input->numeric(FLERR,arg[jarg+2]);
	vcm[2] = input->numeric(FLERR,arg[jarg+3]);
	jarg += 4;
      } else if (strcmp(arg[jarg],"angmom") == 0) {
	angmomflag = 1;
	angmom[0] = input->numeric(FLERR,arg[jarg+1]);
	angmom[1] = input->numeric(FLERR,arg[jarg+2]);
	angmom[2] = input->numeric(FLERR,arg[jarg+3]);
	jarg += 4;
      } else
	error->all(FLERR,"Fix rigid body keyword not recognized");

    }
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

  pseudoflag = 0;
  outfile = NULL;
  outevery = 0;
  remapmode = CUTCELL;
  pushflag = 0;
  pushboundflag = 0;
  pushstyle = LINEAR;
  gammapush = 0.0;
  int pushstyleflag = 0;
  double scale = 1.0;

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
    } else if (strcmp(arg[iarg],"pseudo") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Fix rigid body args not valid");
      pseudoflag = 1;
      nparticle_user = input->inumeric(FLERR,arg[iarg+1]);
      pmass_user = input->numeric(FLERR,arg[iarg+2]);
      frac_user = input->numeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Fix rigid body args not valid");
      scale = input->numeric(FLERR,arg[iarg+1]);
      if (scale <= 0.0)
        error->all(FLERR,"Fix rigid scale factor must be positive");
      iarg += 2;
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
  }

  // RNG for pseudo particles

  random = NULL;
  if (pseudoflag) random = new RanKnuth(update->ranmaster->uniform());

  // apply mass scaling to body params which depend on mass
  // whether defined by fix rigid keywords or read from infile

  if (scale != 1.0) {
    massbody *= scale;
    moi[0] *= scale;
    moi[1] *= scale;
    moi[2] *= scale;
    moi[3] *= scale;
    moi[4] *= scale;
    moi[5] *= scale;
    angmom[0] *= scale;
    angmom[1] *= scale;
    angmom[2] *= scale;
  }
  
  // setup the rigid body

  setup_body();

  // irigid = per-surf flags, indexed by local surf index
  // -1 for static surfs, else index into slist of body surfs
  // accessed by Update::move() to detect moving surfs
  // every proc stores all surfs (non-distributed), so length = nlocal

  int nslocal = surf->nlocal;
  memory->create(irigid,nslocal,"fix_rigid:irigid");

  for (int i = 0; i < nslocal; i++) irigid[i] = -1;
  for (int i = 0; i < nsurf; i++) irigid[slist[i]] = i;

  // remap data structs
  // body surfs are cut/split into grid cells by the normal surf
  //   pipeline (done at read_surf time), so no special setup is needed
  //   here beyond the incremental re-cut work buffers below

  ndeleted = 0;

  pbodyflag = 0;
  noldinside = maxoldinside = 0;
  oldinside = NULL;
  nreg = maxreg = 0;
  regcell = NULL;
  reglist = NULL;
  newlist = NULL;
  newmap = NULL;
  cut2d = NULL;
  cut3d = NULL;

  // for incremental mode: cutters and work bufs for re-cutting cells

  if (remapmode == INCREMENTAL) {
    if (dim == 2) cut2d = new Cut2d(sparta,0);
    else cut3d = new Cut3d(sparta);
    memory->create(newlist,grid->maxsurfpercell,"fix_rigid:newlist");
    memory->create(newmap,grid->maxsurfpercell,"fix_rigid:newmap");
  }
}

/* ---------------------------------------------------------------------- */
 
FixRigid::~FixRigid()
{
  delete [] csurfID;
  delete [] infile;
  delete [] outfile;
  delete random;
  memory->destroy(slist);
  memory->destroy(displace);
  memory->destroy(irigid);

  free_registry();
  memory->destroy(regcell);
  memory->sfree(reglist);
  memory->destroy(oldinside);
  memory->destroy(newlist);
  memory->destroy(newmap);
  delete cut2d;
  delete cut3d;
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

  // check that specified compute is valid for use with fix rigid
  // NOTE: check that it operates on same surf group ?

  int n = modify->find_compute(csurfID);
  if (n < 0) error->all(FLERR,"Could not find fix rigid compute ID");
  if (strcmp(modify->compute[n]->style,"surf") != 0)
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

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int cbit = csurf->surf_groupbit();

  int mask,transparent,isr;
  for (int i = 0; i < nsurf; i++) {
    if (dim == 2) {
      mask = lines[slist[i]].mask;
      transparent = lines[slist[i]].transparent;
      isr = lines[slist[i]].isr;
    } else {
      mask = tris[slist[i]].mask;
      transparent = tris[slist[i]].transparent;
      isr = tris[slist[i]].isr;
    }
    if (transparent)
      error->all(FLERR,"Fix rigid body surfs cannot be transparent");
    if (isr >= 0)
      error->all(FLERR,"Fix rigid body surfs cannot have surface reactions");
    if (!(mask & cbit))
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

  warnrotate = warntranslate = warnexit = 0;

  // each fix rigid defines its own body: no surf can be in two bodies

  for (int ifix = 0; ifix < modify->nfix; ifix++) {
    if (modify->fix[ifix] == this) continue;
    if (strcmp(modify->fix[ifix]->style,"rigid") != 0) continue;
    FixRigid *other = (FixRigid *) modify->fix[ifix];
    for (int i = 0; i < nsurf; i++)
      if (other->irigid[slist[i]] >= 0)
        error->all(FLERR,"Surf element is in more than one fix rigid body");
  }

  // fix rigid must be defined before fixes which change the grid,
  // so its end_of_step() restores overlaid grid cells before they run

  int myindex = modify->find_fix(id);
  for (int ifix = 0; ifix < myindex; ifix++)
    if (strcmp(modify->fix[ifix]->style,"balance") == 0 ||
        strcmp(modify->fix[ifix]->style,"adapt") == 0)
      error->all(FLERR,
                 "Fix rigid must be defined before fix balance or fix adapt");
}

/* ----------------------------------------------------------------------
   called at start of each run, after grid and particles are setup
------------------------------------------------------------------------- */

void FixRigid::setup()
{
  // delete any particles inside the body
  // create_particles marks the body's cells INSIDE via the surf pipeline
  //   and normally avoids them, but this is a safety net for any that
  //   end up inside, e.g. via an emit region overlapping the body

  if (particle->exist) ndeleted += remove_inside_particles(0);

  // for incremental remap: grid state is now consistent with the
  //   body at its current position

  if (remapmode == INCREMENTAL) {
    body_bbox();
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

  // time integrate from current position to end-of-step position
  // use full-step semi-implicit Euler algorithm
  // apply forces and torques accumulated from collisions during last step
  // update vcm/angmom/omega to start-of-step values
  // use them to calculate xcmnew/quatnew and exyz_space for end-of-step values
  
  double dt = update->dt;
  double dtf = dt / massbody;
  double dthalf = 0.5 * dt;
  
  // update vcm by full step

  vcm[0] += dtf * fcm[0];
  vcm[1] += dtf * fcm[1];
  vcm[2] += dtf * fcm[2];

  // update xcm by full step
  // use of new vcm turns Euler into semi-implicit Euler
  // store as xcmnew so have start/stop position for this timestep
  
  xcmnew[0] = xcm[0] + dt * vcm[0];
  xcmnew[1] = xcm[1] + dt * vcm[1];
  xcmnew[2] = xcm[2] + dt * vcm[2];

  // update angular momentum in spatial frame by full step

  angmom[0] += dt * torque[0];
  angmom[1] += dt * torque[1];
  angmom[2] += dt * torque[2];

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
  // use dq/dt = 1/2 omega q
  // store as qusatnew so have start/stop orientation for this timestep
  
  double wq[4];
  MathExtra::vecquat(omega,quat,wq);
  quatnew[0] = quat[0] + dthalf * wq[0];
  quatnew[1] = quat[1] + dthalf * wq[1];
  quatnew[2] = quat[2] + dthalf * wq[2];
  quatnew[3] = quat[3] + dthalf * wq[3];
  MathExtra::qnormalize(quatnew);
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
}

/* ---------------------------------------------------------------------- */

void FixRigid::end_of_step()
{
  // invoke compute surf and extract per-surf force/torque info
  // NOTE: this could access fix ave/surf for f/t from many steps of collisions ?

  if (!pseudoflag) {

    if (!(csurf->invoked_flag & INVOKED_PER_SURF)) {
      csurf->compute_per_surf();
      csurf->invoked_flag |= INVOKED_PER_SURF;
    }

    csurf->post_process_surf();
    double **array = csurf->array_surf;

    // sum per-surf force/torque to body fcm/torque
    // array rows are surfs each proc owns:
    //   for explicit non-distributed surfs, proc owns every Pth surf,
    //   row M of array = surf with local index me + M*nprocs
    // sum rows for owned surfs in body group, then Allreduce for all procs

    Surf::Line *lines = surf->lines;
    Surf::Tri *tris = surf->tris;
    int nslocal = surf->nlocal;
    int me = comm->me;
    int nprocs = comm->nprocs;

    double ft_mine[6],ft_all[6];
    for (int i = 0; i < 6; i++) ft_mine[i] = 0.0;

    int mask;
    int m = 0;
    for (int i = me; i < nslocal; i += nprocs, m++) {
      if (dim == 2) mask = lines[i].mask;
      else mask = tris[i].mask;
      if (!(mask & groupbit)) continue;
      for (int j = 0; j < 6; j++) ft_mine[j] += array[m][j];
    }

    MPI_Allreduce(ft_mine,ft_all,6,MPI_DOUBLE,MPI_SUM,world);

    fcm[0] = ft_all[0];
    fcm[1] = ft_all[1];
    fcm[2] = ft_all[2];
    torque[0] = ft_all[3];
    torque[1] = ft_all[4];
    torque[2] = ft_all[5];

    // insure the compute tallies on the next step

    csurf->addstep(update->ntimestep+1);
  }

  // bounce N pseudo (fictitious) particles off object, see how it moves
  // bbox = epsilon-augmented bbox around current object
  // shoot N particles from random positions from left face to right face of bbox
  // whichever surf it hits first, compute force/torque on object

  if (pseudoflag) {

    fcm[0] = fcm[1] = fcm[2] = 0.0;
    torque[0] = torque[1] = torque[2] = 0.0;

    Surf::Line *lines = surf->lines;
    Surf::Tri *tris = surf->tris;
    Surf::Line *line;
    Surf::Tri *tri;
  
    int index;
    double bboxlo[3],bboxhi[3];

    if (dim == 2) {
      bboxlo[0] = bboxlo[1] = 1.e20;
      bboxhi[0] = bboxhi[1] = -1.e20;

      for (int i = 0; i < nsurf; i++) {
	index = slist[i];
	bboxlo[0] = MIN(bboxlo[0],lines[index].p1[0]);
	bboxlo[0] = MIN(bboxlo[0],lines[index].p2[0]);
	bboxhi[0] = MAX(bboxhi[0],lines[index].p1[0]);
	bboxhi[0] = MAX(bboxhi[0],lines[index].p2[0]);
	bboxlo[1] = MIN(bboxlo[1],lines[index].p1[1]);
	bboxlo[1] = MIN(bboxlo[1],lines[index].p2[1]);
	bboxhi[1] = MAX(bboxhi[1],lines[index].p1[1]);
	bboxhi[1] = MAX(bboxhi[1],lines[index].p2[1]);
	bboxlo[2] = bboxhi[2] = 0.0;
      }
    } else if (dim == 3) {
      bboxlo[0] = bboxlo[1] = bboxlo[2] = 1.e20;
      bboxhi[0] = bboxhi[1] = bboxhi[2] = -1.e20;
      
      for (int i = 0; i < nsurf; i++) {
	index = slist[i];
	bboxlo[0] = MIN(bboxlo[0],tris[index].p1[0]);
	bboxlo[0] = MIN(bboxlo[0],tris[index].p2[0]);
	bboxlo[0] = MIN(bboxlo[0],tris[index].p3[0]);
	bboxhi[0] = MAX(bboxhi[0],tris[index].p1[0]);
	bboxhi[0] = MAX(bboxhi[0],tris[index].p2[0]);
	bboxhi[0] = MAX(bboxhi[0],tris[index].p3[0]);
	bboxlo[1] = MIN(bboxlo[1],tris[index].p1[1]);
	bboxlo[1] = MIN(bboxlo[1],tris[index].p2[1]);
	bboxlo[1] = MIN(bboxlo[1],tris[index].p3[1]);
	bboxhi[1] = MAX(bboxhi[1],tris[index].p1[1]);
	bboxhi[1] = MAX(bboxhi[1],tris[index].p2[1]);
	bboxhi[1] = MAX(bboxhi[1],tris[index].p3[1]);
	bboxlo[2] = MIN(bboxlo[2],tris[index].p1[2]);
	bboxlo[2] = MIN(bboxlo[2],tris[index].p2[2]);
	bboxlo[2] = MIN(bboxlo[2],tris[index].p3[2]);
	bboxhi[2] = MAX(bboxhi[2],tris[index].p1[2]);
	bboxhi[2] = MAX(bboxhi[2],tris[index].p2[2]);
	bboxhi[2] = MAX(bboxhi[2],tris[index].p3[2]);
      }
    }

    // expand bbox by one percent in all dims

    bboxlo[0] -= 0.01 * (bboxhi[0]-bboxlo[0]);
    bboxlo[1] -= 0.01 * (bboxhi[1]-bboxlo[1]);
    bboxlo[2] -= 0.01 * (bboxhi[2]-bboxlo[2]);
    bboxhi[0] += 0.01 * (bboxhi[0]-bboxlo[0]);
    bboxhi[1] += 0.01 * (bboxhi[1]-bboxlo[1]);
    bboxhi[2] += 0.01 * (bboxhi[2]-bboxlo[2]);
    
    int nhits = 0;
    
    if (dim == 2) {
      for (int m = 0; m < nparticle_user; m++) {
	int side,minsurf;
	double param;
	double x[3],xnew[3];
	double xc[3],minxc[3];
      
	// set x,xnew randomly for each of N particles	
	// limit by fraction of bbox face in y

	x[0] = bboxlo[0];
	xnew[0] = bboxhi[0];
	double rn = random->uniform();
	xnew[1] = x[1] = bboxlo[1] + frac_user * rn * (bboxhi[1]-bboxlo[1]);
	xnew[2] = x[2] = 0.0;
      
	// initial vel = +x
	// final vel = reflect off surf

	double vpre[3],vpost[3];
	vpre[0] = 1.0;
	vpre[1] = vpre[2] = 0.0;
      
	int cflag = 0;
	double minparam = 2.0;
      
	for (int i = 0; i < nsurf; i++) {
	  index = slist[i];
	  line = &lines[index];
	  int hitflag = Geometry::
	    line_line_intersect(x,xnew,line->p1,line->p2,
				line->norm,xc,param,side);
	  if (hitflag && param < minparam && side == OUTSIDE) {
	    cflag = 1;
	    minparam = param;
	    minsurf = index;
	    minxc[0] = xc[0];
	    minxc[1] = xc[1];
	    minxc[2] = 0.0;
	  }
	}

	// add force/torque from collision
	
	if (cflag) {
	  nhits++;
	  double pforce[3],rdelta[3],tq[3];
	  
	  vpost[0] = vpre[0]; vpost[1] = vpre[1]; vpost[2] = vpre[2];
	  MathExtra::reflect3(vpost,lines[minsurf].norm);
	  
	  pforce[0] = pforce[1] = pforce[2] = 0.0;
	  MathExtra::axpy3(pmass_user,vpre,pforce);
	  MathExtra::axpy3(-pmass_user,vpost,pforce);
	  fcm[0] += pforce[0];
	  fcm[1] += pforce[1];
	  fcm[2] += pforce[2];
	  
	  MathExtra::sub3(minxc,xcm,rdelta);
	  MathExtra::cross3(rdelta,pforce,tq);
	  torque[0] += tq[0];
	  torque[1] += tq[1];
	  torque[2] += tq[2];
	}
      }
      
    } else if (dim == 3) {
      for (int m = 0; m < nparticle_user; m++) {
	int side,minsurf;
	double param;
	double x[3],xnew[3];
	double xc[3],minxc[3];

	// set x,xnew randomly for each of N particles
	// limit by fraction of bbox face in z
      
	x[0] = bboxlo[0];
	xnew[0] = bboxhi[0];
	double rn = random->uniform();
	xnew[1] = x[1] = bboxlo[1] + rn * (bboxhi[1]-bboxlo[1]);
	rn = random->uniform();
	xnew[2] = x[2] = bboxlo[2] + frac_user * rn * (bboxhi[2]-bboxlo[2]);
	
	// initial vel = +x
	// final vel = reflect off surf
	
	double vpre[3],vpost[3];
	vpre[0] = 1.0;
	vpre[1] = vpre[2] = 0.0;
	
	int cflag = 0;
	double minparam = 2.0;
	
	for (int i = 0; i < nsurf; i++) {
	  index = slist[i];
	  tri = &tris[index];
	  int hitflag = Geometry::
	    line_tri_intersect(x,xnew,tri->p1,tri->p2,tri->p3,
			       tri->norm,xc,param,side);
	  if (hitflag && param < minparam && side == OUTSIDE) {
	    cflag = 1;
	    minparam = param;
	    minsurf = index;
	    minxc[0] = xc[0];
	    minxc[1] = xc[1];
	    minxc[2] = xc[2];
	  }
	}
	
	// add force/torque from collision
	
	if (cflag) {
	  nhits++;
	  double pforce[3],rdelta[3],tq[3];
	  
	  vpost[0] = vpre[0]; vpost[1] = vpre[1]; vpost[2] = vpre[2];
	  MathExtra::reflect3(vpost,tris[minsurf].norm);
	  
	  pforce[0] = pforce[1] = pforce[2] = 0.0;
	  MathExtra::axpy3(pmass_user,vpre,pforce);
	  MathExtra::axpy3(-pmass_user,vpost,pforce);
	  fcm[0] += pforce[0];
	  fcm[1] += pforce[1];
	  fcm[2] += pforce[2];
	  
	  MathExtra::sub3(minxc,xcm,rdelta);
	  MathExtra::cross3(rdelta,pforce,tq);
	  torque[0] += tq[0];
	  torque[1] += tq[1];
	  torque[2] += tq[2];
	}
      }
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

  // update Surf class properties of lines and tris in body
  // line and tri positions and orientations
  // set via Line/Tri end/corner points and norm
  // matvec() converts displace vector from body frame to space frame

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int index;
  
  if (dim == 2) {
    double z[3],delta[3];
    z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;
    
    for (int i = 0; i < nsurf; i++) {
      index = slist[i];

      MathExtra::matvec(ex_space,ey_space,ez_space,displace[i][0],delta);
      delta[2] = 0.0;
      MathExtra::add3(xcm,delta,lines[index].p1);
      MathExtra::matvec(ex_space,ey_space,ez_space,displace[i][1],delta);
      delta[2] = 0.0;
      MathExtra::add3(xcm,delta,lines[index].p2);

      MathExtra::sub3(lines[index].p2,lines[index].p1,delta);
      MathExtra::cross3(z,delta,lines[index].norm);
      MathExtra::norm3(lines[index].norm);
      lines[index].norm[2] = 0.0;
    }

  } else if (dim == 3) {
    double delta[3],delta12[3],delta13[3];

    for (int i = 0; i < nsurf; i++) {
      index = slist[i];

      MathExtra::matvec(ex_space,ey_space,ez_space,displace[i][0],delta);
      MathExtra::add3(xcm,delta,tris[index].p1);
      MathExtra::matvec(ex_space,ey_space,ez_space,displace[i][1],delta);
      MathExtra::add3(xcm,delta,tris[index].p2);
      MathExtra::matvec(ex_space,ey_space,ez_space,displace[i][2],delta);
      MathExtra::add3(xcm,delta,tris[index].p3);

      MathExtra::sub3(tris[index].p2,tris[index].p1,delta12);
      MathExtra::sub3(tris[index].p3,tris[index].p1,delta13);
      MathExtra::cross3(delta12,delta13,tris[index].norm);
      MathExtra::norm3(tris[index].norm);
    }
  }

  // bbox around body elements at their new positions

  body_bbox();

  // push-off forces from static surfs and domain boundaries
  //   which are within pushcutoff of the body
  // computed for the end-of-step geometry, so they are part of the
  //   force/torque applied in the next step's time integration

  fpush[0] = fpush[1] = fpush[2] = 0.0;
  if (pushflag) push_off();

  // error if body now extends beyond a periodic boundary,
  //   b/c body coords are not wrapped across periodic boundaries
  // body is allowed to exit thru non-periodic boundaries

  int outflag = 0;
  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  int *bflag = domain->bflag;

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

  if (bflag[0] == PERIODIC && bbodylo[0] < boxlo[0]) outflag = 1;
  if (bflag[1] == PERIODIC && bbodyhi[0] > boxhi[0]) outflag = 1;
  if (bflag[2] == PERIODIC && bbodylo[1] < boxlo[1]) outflag = 1;
  if (bflag[3] == PERIODIC && bbodyhi[1] > boxhi[1]) outflag = 1;
  if (dim == 3) {
    if (bflag[4] == PERIODIC && bbodylo[2] < boxlo[2]) outflag = 1;
    if (bflag[5] == PERIODIC && bbodyhi[2] > boxhi[2]) outflag = 1;
  }

  if (outflag)
    error->all(FLERR,"Fix rigid body moved beyond a periodic boundary");

  // re-map body surfs to grid cells: cut/split cells and INSIDE/OUTSIDE
  //   cell typing from the new body positions
  // all bodies are re-mapped together by the last-defined rigid fix,
  //   once per step, after every body has moved to its end-of-step
  //   position; end_of_step runs fixes in definition order, so when the
  //   last fix runs, all earlier bodies are already moved
  // if every body is incremental, attempt the cheap incremental re-cut
  //   of only the affected cells; otherwise (any cutcell body) do an
  //   exact full grid re-map, which correctly handles all bodies
  // the incremental fallback decision is per-proc but a full re-map is
  //   collective, so all procs must agree via Allreduce
  // then remove particles inside each body, with split-cell particle
  //   reassignment done once (only needed after a full re-map)

  int ilast = -1;
  int all_incremental = 1;
  for (int ifix = 0; ifix < modify->nfix; ifix++) {
    if (strcmp(modify->fix[ifix]->style,"rigid") != 0) continue;
    ilast = ifix;
    if (((FixRigid *) modify->fix[ifix])->remapmode != INCREMENTAL)
      all_incremental = 0;
  }

  if (modify->fix[ilast] == this) {
    int fallback = 1;
    if (all_incremental) {
      int fallmine = incremental_recut();
      MPI_Allreduce(&fallmine,&fallback,1,MPI_INT,MPI_MAX,world);
    }
    if (fallback) grid_rebuild();

    if (particle->exist) {
      int splitflag = fallback;
      for (int ifix = 0; ifix < modify->nfix; ifix++) {
        if (strcmp(modify->fix[ifix]->style,"rigid") != 0) continue;
        FixRigid *f = (FixRigid *) modify->fix[ifix];
        f->ndeleted += f->remove_inside_particles(splitflag);
        splitflag = 0;
      }
    }

    // advance each incremental body's previous-region bookkeeping

    for (int ifix = 0; ifix < modify->nfix; ifix++) {
      if (strcmp(modify->fix[ifix]->style,"rigid") != 0) continue;
      FixRigid *f = (FixRigid *) modify->fix[ifix];
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
  fprintf(fp,"# mtotal xcm ycm zcm ixx iyy izz ixy ixz iyz "
          "vxcm vycm vzcm lx ly lz\n");
  fprintf(fp,"%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g "
          "%.15g %.15g %.15g %.15g %.15g %.15g\n",
          massbody,xcm[0],xcm[1],xcm[2],
          ispace[0],ispace[1],ispace[2],ispace[3],ispace[4],ispace[5],
          vcm[0],vcm[1],vcm[2],angmom[0],angmom[1],angmom[2]);

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
    
    int ncorrect = 16;
    int nwords = input->count_words(line);
    if (nwords != ncorrect)
      error->all(FLERR,"Incorrect rigid body format in fix rigid infile");

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

    fclose(fp);
  }

  // broadcast result of file read to all procs
    
  MPI_Bcast(&massbody,1,MPI_DOUBLE,0,world);
  MPI_Bcast(xcm,3,MPI_DOUBLE,0,world);
  MPI_Bcast(moi,6,MPI_DOUBLE,0,world);
  MPI_Bcast(vcm,3,MPI_DOUBLE,0,world);
  MPI_Bcast(angmom,3,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   one-time initialization of rigid body attributes
------------------------------------------------------------------------- */

void FixRigid::setup_body()
{
  // nsurf = # of lines/tris in rigid body
  // NOTE: add check that group is a closed object ?

  bigint bnsurf = surf->count_group(igroup);
  if (bnsurf > MAXSMALLINT) error->all(FLERR,"Too many surfs in rigid body");
  nsurf = bnsurf;
  if (nsurf == 0) error->all(FLERR,"Fix rigid body has no surface elements");
  //if (check) error->all(FLERR,"Fix rigid surfs are an invalid object);

  // slist = list of surf indices in surf group

  memory->create(slist,nsurf,"fix_rigid:slist");

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nlocal = surf->nlocal;

  int n = 0;
  if (dim == 2) {
    for (int i = 0; i < nlocal; i++)
      if (lines[i].mask & groupbit)
	slist[n++] = i;
  } else if (dim == 3) {
    for (int i = 0; i < nlocal; i++)
      if (tris[i].mask & groupbit)
	slist[n++] = i;
  }

  // insure body surfs form a closed (watertight) object

  check_watertight();

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

  // set displacement for each end/corner point in each line/tri
  // delta = vector from COM to end/corner point in space frame
  // displace = delta rotated to be in basis of principal axes, i.e. in body frame
  
  double delta[3];
  int index;
  
  if (dim == 2) {
    memory->create(displace,nsurf,2,3,"fix_rigid:displace");
    
    for (int i = 0; i < nsurf; i++) {
      index = slist[i];

      delta[0] = lines[index].p1[0] - xcm[0];
      delta[1] = lines[index].p1[1] - xcm[1];
      delta[2] = 0.0;
      MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
				  delta,&displace[i][0][0]);
      delta[0] = lines[index].p2[0] - xcm[0];
      delta[1] = lines[index].p2[1] - xcm[1];
      delta[2] = 0.0;
      MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
				  delta,&displace[i][1][0]);
    }
    
  } else if (dim == 3) {
    memory->create(displace,nsurf,3,3,"fix_rigid:displace");

    for (int i = 0; i < nsurf; i++) {
      index = slist[i];

      delta[0] = tris[index].p1[0] - xcm[0];
      delta[1] = tris[index].p1[1] - xcm[1];
      delta[2] = tris[index].p1[2] - xcm[2];
      MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
				  delta,&displace[i][0][0]);
      delta[0] = tris[index].p2[0] - xcm[0];
      delta[1] = tris[index].p2[1] - xcm[1];
      delta[2] = tris[index].p2[2] - xcm[2];
      MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
				  delta,&displace[i][1][0]);
      delta[0] = tris[index].p3[0] - xcm[0];
      delta[1] = tris[index].p3[1] - xcm[1];
      delta[2] = tris[index].p3[2] - xcm[2];
      MathExtra::transpose_matvec(ex_space,ey_space,ez_space,
				  delta,&displace[i][2][0]);
    }
  }

  // initial omega, consistent with initial angmom

  MathExtra::angmom_to_omega(angmom,ex_space,ey_space,ez_space,inertia,omega);

  // rmaxbody = max distance of any body corner pt from the COM

  rmaxbody = 0.0;
  for (int i = 0; i < nsurf; i++)
    for (int j = 0; j < dim; j++)
      rmaxbody = MAX(rmaxbody,MathExtra::len3(&displace[i][j][0]));

  // zero body force/torque in case accessed via compute_vector() on step 0

  fcm[0] = fcm[1] = fcm[2] = 0.0;
  torque[0] = torque[1] = torque[2] = 0.0;
  fpush[0] = fpush[1] = fpush[2] = 0.0;
}

/* ----------------------------------------------------------------------
   push-off forces on the body from too-close static surfs
   for each body element corner pt within pushcutoff of a static surf,
     apply a repulsive force in the direction of the static surf
     outward normal, with overlap delta = pushcutoff - dist:
     linear spring F = kpush * delta, or
     Hertzian contact F = kpush * delta^3/2 (smooth onset, standard
     model for elastic contact of spherical particulates)
   if gammapush > 0, a dashpot term F -= gammapush * d(delta)/dt is
     added (the DEM spring-dashpot pair), computed from the normal
     approach rate of the corner pt relative to the source surface,
     including the motion of another rigid body as the source
   the total contact force is clamped at zero, so the dashpot never
     produces adhesion as a contact ends
   if pushboundflag is set, the same spring force is applied by
     non-periodic simulation box boundaries
   forces are added to fcm/torque for the next step's time integration
     and accumulated in fpush for diagnostic output
   computed identically on every proc: all surfs are stored everywhere,
     so no communication is needed and all procs stay in sync
   NOTE: a corner pt shared by adjacent body elements contributes once
     per element, and a corner close to several static elements
     interacts with each of them, so kpush is a per-contact stiffness
------------------------------------------------------------------------- */

void FixRigid::push_off()
{
  int i,j,m,index;
  double dsq,d,scale;
  double *pts[3],*norm;
  double fone[3],rdelta[3],tq[3];

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nslocal = surf->nlocal;

  int npoint = dim;     // 2 corner pts per line, 3 per tri
  double cutsq = pushcutoff*pushcutoff;

  // cutlo/cuthi = bbox around body inflated by pushcutoff
  // requires body_bbox() was called for current body position

  double cutlo[3],cuthi[3];
  for (j = 0; j < 3; j++) {
    cutlo[j] = bbodylo[j] - pushcutoff;
    cuthi[j] = bbodyhi[j] + pushcutoff;
  }

  // loop over static surfs whose bbox overlaps the inflated body bbox

  for (m = 0; m < nslocal; m++) {
    if (irigid[m] >= 0) continue;

    if (dim == 2) {
      if (MAX(lines[m].p1[0],lines[m].p2[0]) < cutlo[0]) continue;
      if (MIN(lines[m].p1[0],lines[m].p2[0]) > cuthi[0]) continue;
      if (MAX(lines[m].p1[1],lines[m].p2[1]) < cutlo[1]) continue;
      if (MIN(lines[m].p1[1],lines[m].p2[1]) > cuthi[1]) continue;
      norm = lines[m].norm;
    } else {
      if (MAX(tris[m].p1[0],MAX(tris[m].p2[0],tris[m].p3[0])) < cutlo[0])
        continue;
      if (MIN(tris[m].p1[0],MIN(tris[m].p2[0],tris[m].p3[0])) > cuthi[0])
        continue;
      if (MAX(tris[m].p1[1],MAX(tris[m].p2[1],tris[m].p3[1])) < cutlo[1])
        continue;
      if (MIN(tris[m].p1[1],MIN(tris[m].p2[1],tris[m].p3[1])) > cuthi[1])
        continue;
      if (MAX(tris[m].p1[2],MAX(tris[m].p2[2],tris[m].p3[2])) < cutlo[2])
        continue;
      if (MIN(tris[m].p1[2],MIN(tris[m].p2[2],tris[m].p3[2])) > cuthi[2])
        continue;
      norm = tris[m].norm;
    }

    // spring force on each body corner pt within pushcutoff

    for (i = 0; i < nsurf; i++) {
      index = slist[i];
      if (dim == 2) {
        pts[0] = lines[index].p1;
        pts[1] = lines[index].p2;
      } else {
        pts[0] = tris[index].p1;
        pts[1] = tris[index].p2;
        pts[2] = tris[index].p3;
      }

      for (j = 0; j < npoint; j++) {
        if (dim == 2)
          dsq = Geometry::distsq_point_line(pts[j],lines[m].p1,lines[m].p2);
        else
          dsq = Geometry::distsq_point_tri(pts[j],tris[m].p1,tris[m].p2,
                                           tris[m].p3,tris[m].norm);
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
          if (update->rigidflag && update->rigidmap &&
              update->rigidmap[m] >= 0) {
            FixRigid *src = update->fixrigidlist[update->rigidmap[m]];
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
        torque[0] += tq[0];
        torque[1] += tq[1];
        torque[2] += tq[2];
      }
    }
  }

  // spring force from non-periodic simulation box boundaries

  if (pushboundflag) {
    double *boxlo = domain->boxlo;
    double *boxhi = domain->boxhi;
    int *bflag = domain->bflag;

    int nface = 2*dim;
    double fsign[6] = {1.0,-1.0,1.0,-1.0,1.0,-1.0};

    for (i = 0; i < nsurf; i++) {
      index = slist[i];
      if (dim == 2) {
        pts[0] = lines[index].p1;
        pts[1] = lines[index].p2;
      } else {
        pts[0] = tris[index].p1;
        pts[1] = tris[index].p2;
        pts[2] = tris[index].p3;
      }

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
          torque[0] += tq[0];
          torque[1] += tq[1];
          torque[2] += tq[2];
        }
      }
    }
  }

  // add push-off force to body force for next step's integration

  fcm[0] += fpush[0];
  fcm[1] += fpush[1];
  fcm[2] += fpush[2];
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

  grid->notify_changed();
}

/* ----------------------------------------------------------------------
   grid cells were rebuilt or migrated to other procs
   next re-map re-cuts body surfs into the new grid cells
------------------------------------------------------------------------- */

void FixRigid::grid_changed()
{
  // csurfs lists installed by incremental re-cutting were discarded
  //   by the grid rebuild, or cell indices changed due to migration

  free_registry();
}

/* ----------------------------------------------------------------------
   for incremental remap: record cells interior to the body,
     i.e. INSIDE cells with no surfs whose center is within the body
   called before the body surfs move to their end-of-step positions
------------------------------------------------------------------------- */

void FixRigid::record_oldinside()
{
  double ctr[3];

  // bbox around body at its current (pre-move) position

  body_bbox();

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int nglocal = grid->nlocal;

  noldinside = 0;

  for (int icell = 0; icell < nglocal; icell++) {
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
   cells interior to the body at its old or new position are re-typed
     as INSIDE/OUTSIDE via parity tests, all other cells are untouched
   ghost cell copies of re-cut cells become stale, which is acceptable:
     the particle mover only consults surf lists of owned cells
   return 1 to request a fallback to a full grid re-map if a
     structural change occurs:
     a cell would become or stop being a split cell, a cell's surf
     count exceeds maxsurfpercell, or no previous body position is set
------------------------------------------------------------------------- */

int FixRigid::incremental_recut()
{
  int i,n,icell,nsplitone,xsub,moving;
  double vol;
  double xsplit[3],ctr[3],rlo[3],rhi[3];
  double *vols;
  double *clo,*chi;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
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

  for (int ifix = 0; ifix < modify->nfix; ifix++) {
    if (strcmp(modify->fix[ifix]->style,"rigid") != 0) continue;
    FixRigid *f = (FixRigid *) modify->fix[ifix];
    if (f->remapmode != INCREMENTAL) continue;
    if (!f->pbodyflag) return 1;
    for (i = 0; i < 3; i++) {
      rlo[i] = MIN(rlo[i],MIN(f->pbodylo[i],f->bbodylo[i]));
      rhi[i] = MAX(rhi[i],MAX(f->pbodyhi[i],f->bbodyhi[i]));
    }
    nincr++;
  }
  if (!nincr) return 1;

  // pass 1: re-cut cells in R whose surf overlap changed
  //   or which are overlapped by a moved body surf (from any body)
  // NOTE: surf lists are compared elementwise, both in cut2d/cut3d
  //   surf index order; lists built by the rendezvous surf2grid
  //   algorithm may be ordered differently, causing a one-time
  //   spurious re-cut of unchanged cells, which is harmless

  for (icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (!box_overlap(cells[icell].lo,cells[icell].hi,rlo,rhi)) continue;

    // structural change unsupported: split cells trigger a full re-map

    if (cells[icell].nsplit > 1) return 1;

    // new list of surfs overlapping this cell

    if (dim == 2)
      n = cut2d->surf2grid(cells[icell].id,cells[icell].lo,cells[icell].hi,
                           newlist,maxsurfpercell);
    else
      n = cut3d->surf2grid(cells[icell].id,cells[icell].lo,cells[icell].hi,
                           newlist,maxsurfpercell);
    if (n > maxsurfpercell) return 1;

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

    if (n == 0) {

      // cell no longer overlaps any surf
      // full flow volume, interior/exterior typing via parity test

      registry_remove(icell);
      cells[icell].nsurf = 0;
      cells[icell].csurfs = NULL;

      clo = cells[icell].lo;
      chi = cells[icell].hi;
      if (dim == 3)
        vol = (chi[0]-clo[0]) * (chi[1]-clo[1]) * (chi[2]-clo[2]);
      else vol = (chi[0]-clo[0]) * (chi[1]-clo[1]);
      cinfo[icell].volume = vol;

      ctr[0] = 0.5 * (clo[0] + chi[0]);
      ctr[1] = 0.5 * (clo[1] + chi[1]);
      if (dim == 3) ctr[2] = 0.5 * (clo[2] + chi[2]);
      else ctr[2] = 0.0;

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

      if (nsplitone > 1) return 1;

      cinfo[icell].volume = vols[0];
      cinfo[icell].type = CELLOVERLAP;
    }
  }

  // pass 2: cells a body interior moved away from become OUTSIDE
  // process every incremental body's recorded interior cells
  // only cells which are now surf-free and inside no body,
  //   which leaves any static (non-body) INSIDE cells untouched

  for (int ifix = 0; ifix < modify->nfix; ifix++) {
    if (strcmp(modify->fix[ifix]->style,"rigid") != 0) continue;
    FixRigid *f = (FixRigid *) modify->fix[ifix];
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

  for (icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit != 1) continue;
    if (cells[icell].nsurf) continue;
    if (cinfo[icell].type == CELLINSIDE) continue;
    if (!box_overlap(cells[icell].lo,cells[icell].hi,rlo,rhi)) continue;

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

  return 0;
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
  for (int i = 0; i < nreg; i++)
    if (regcell[i] == icell) {
      memory->sfree(reglist[i]);
      reglist[i] = list;
      return;
    }

  if (nreg == maxreg) {
    maxreg += DELTA_MODIFY;
    memory->grow(regcell,maxreg,"fix_rigid:regcell");
    reglist = (surfint **)
      memory->srealloc(reglist,maxreg*sizeof(surfint *),
                       "fix_rigid:reglist");
  }

  regcell[nreg] = icell;
  reglist[nreg] = list;
  nreg++;
}

void FixRigid::registry_remove(int icell)
{
  for (int i = 0; i < nreg; i++)
    if (regcell[i] == icell) {
      memory->sfree(reglist[i]);
      regcell[i] = regcell[nreg-1];
      reglist[i] = reglist[nreg-1];
      nreg--;
      return;
    }
}

void FixRigid::free_registry()
{
  for (int i = 0; i < nreg; i++) memory->sfree(reglist[i]);
  nreg = 0;
}

/* ----------------------------------------------------------------------
   compute the bounding box around the whole body at its current position
   box is inflated by EPSSURF * body extent to avoid round-off misses
------------------------------------------------------------------------- */

void FixRigid::body_bbox()
{
  int i,j,k,index;
  double *pts[3];

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  int npoint = dim;     // 2 points per line, 3 per tri

  bbodylo[0] = bbodylo[1] = bbodylo[2] = BIG;
  bbodyhi[0] = bbodyhi[1] = bbodyhi[2] = -BIG;

  for (i = 0; i < nsurf; i++) {
    index = slist[i];
    if (dim == 2) {
      pts[0] = lines[index].p1;
      pts[1] = lines[index].p2;
    } else {
      pts[0] = tris[index].p1;
      pts[1] = tris[index].p2;
      pts[2] = tris[index].p3;
    }

    for (j = 0; j < npoint; j++)
      for (k = 0; k < 3; k++) {
        bbodylo[k] = MIN(bbodylo[k],pts[j][k]);
        bbodyhi[k] = MAX(bbodyhi[k],pts[j][k]);
      }
  }

  double eps = EPSSURF * MAX(bbodyhi[0]-bbodylo[0],bbodyhi[1]-bbodylo[1]);
  eps = EPSSURF * MAX(eps/EPSSURF,bbodyhi[2]-bbodylo[2]);

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
  int index,hitflag,side;
  double param;
  double xout[3],xc[3];

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  double dmax = MAX(bbodyhi[0]-bbodylo[0],bbodyhi[1]-bbodylo[1]);
  dmax = MAX(dmax,bbodyhi[2]-bbodylo[2]);

  xout[0] = bbodyhi[0] + 0.414159*dmax;
  xout[1] = x[1] + 0.271828*dmax;
  if (dim == 3) xout[2] = x[2] + 0.161803*dmax;
  else xout[2] = 0.0;

  int count = 0;
  for (int i = 0; i < nsurf; i++) {
    index = slist[i];
    if (dim == 2)
      hitflag = Geometry::
        line_line_intersect(x,xout,lines[index].p1,lines[index].p2,
                            lines[index].norm,xc,param,side);
    else
      hitflag = Geometry::
        line_tri_intersect(x,xout,tris[index].p1,tris[index].p2,
                           tris[index].p3,tris[index].norm,xc,param,side);
    if (hitflag) count++;
  }

  return count % 2;
}

/* ----------------------------------------------------------------------
   remove particles which are inside the body
   also remove all particles in INSIDE cells (cutcell mode)
   splitflag = 1 if called after a grid rebuild,
     to first reassign particles in split cells to their sub cells
   return # of particles deleted across all procs
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

  body_bbox();

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
  bigint delta = nlocal_old - particle->nlocal;
  bigint nall;
  MPI_Allreduce(&delta,&nall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  return nall;
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

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;

  if (dim == 2) {
    std::map<std::array<double,2>,int> count;
    std::array<double,2> key;

    for (int i = 0; i < nsurf; i++) {
      Surf::Line *l = &lines[slist[i]];
      key[0] = l->p1[0]; key[1] = l->p1[1];
      count[key]++;
      key[0] = l->p2[0]; key[1] = l->p2[1];
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
      Surf::Tri *t = &tris[slist[i]];
      pts[0] = t->p1; pts[1] = t->p2; pts[2] = t->p3; pts[3] = t->p1;

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
   return cummulative count of particles deleted inside the moving body
------------------------------------------------------------------------- */

double FixRigid::compute_scalar()
{
  return (double) ndeleted;
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
