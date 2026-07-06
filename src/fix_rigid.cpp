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

enum{OVERLAY,CUTCELL};          // remap modes

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
  // are rebuilt or migrated, so overlay remap data can be reset

  gridmigrate = 1;

  if (!surf->exist) error->all(FLERR,"Fix rigid requires surf elements exist");
  if (domain->axisymmetric)
    error->all(FLERR,"Fix rigid cannot be used with axisymmetric domains");
  if (surf->implicit || surf->distributed)
    error->all(FLERR,"Fix rigid can only be used with explicit non-distributed surf elements");

  // only a single rigid body is allowed

  for (int ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"rigid") == 0)
      error->all(FLERR,"Only one fix rigid command can be defined");

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
  remapmode = OVERLAY;
  double scale = 1.0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"remap") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Fix rigid body args not valid");
      if (strcmp(arg[iarg+1],"overlay") == 0) remapmode = OVERLAY;
      else if (strcmp(arg[iarg+1],"cutcell") == 0) remapmode = CUTCELL;
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

  // create rigid = custom per-surf vector
  //   unless already exists, due to restart file
  // rigid = index into short list of rigid surfs for mobile surfs
  // rigid = -1 for static surfs

  rigidindex = surf->find_custom((char *) "rigid");
  if (rigidindex < 0) rigidindex = surf->add_custom((char *) "rigid",INT,0);

  int *irigid = surf->eivec[surf->ewhich[rigidindex]];

  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  int nslocal = surf->nlocal;
  
  for (int i = 0; i < nslocal; i++) irigid[i] = -1;
  for (int i = 0; i < nsurf; i++) irigid[slist[i]] = i;

  // remap data structs

  nmodified = maxmodified = 0;
  modified = NULL;
  nsurf_saved = NULL;
  csurfs_saved = NULL;
  cpage = NULL;
  ndeleted = 0;

  // for overlay mode:
  // exclude body surfs from static assignment of surfs to grid cells,
  //   then perform a full re-map of surfs to grid cells
  // cells overlapped by the body then have full flow volume,
  //   are not cut or split by body surfs, and are not marked INSIDE
  // each step, start_of_step() overlays body surfs onto the grid cells
  //   they sweep through, so particles can collide with them

  if (remapmode == OVERLAY) {
    int maxchunk = grid->maxsurfpercell + nsurf;
    cpage = new MyPage<surfint>(maxchunk,MAX(65536,4*maxchunk));
    if (cpage->errorflag)
      error->all(FLERR,"Fix rigid could not allocate overlay page");

    surf->rigidbits |= groupbit;
    grid_rebuild();
  }
}

/* ---------------------------------------------------------------------- */
 
FixRigid::~FixRigid()
{
  // restore any overlaid grid cells and static surf assignment flags
  // NOTE: grid cells are not re-mapped here, so after an unfix the body
  //   surfs do not interact with particles until the next full re-map

  overlay_restore();
  surf->rigidbits &= ~groupbit;

  delete [] csurfID;
  delete [] infile;
  delete [] outfile;
  delete random;
  memory->destroy(slist);
  memory->destroy(displace);
  memory->destroy(elemlo);
  memory->destroy(elemhi);
  memory->destroy(tmplist);
  memory->destroy(modified);
  memory->destroy(nsurf_saved);
  memory->sfree(csurfs_saved);
  delete cpage;
  surf->remove_custom(rigidindex);
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
  // e.g. created by create_particles, which is unaware of the body
  //   when the body's grid cells are not marked INSIDE (overlay mode)

  if (particle->exist) ndeleted += remove_inside_particles(0);
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

  // overlay body surfs onto all grid cells
  //   their swept bbox overlaps at any point during this step

  if (remapmode == OVERLAY) overlay_assign();
}

/* ---------------------------------------------------------------------- */

void FixRigid::end_of_step()
{
  // restore static surf assignment of grid cells overlaid this step
  // must be done before any other operation changes grid or particle data

  if (remapmode == OVERLAY) overlay_restore();

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

  // for cutcell mode: full re-map of surfs to grid cells,
  //   including cut/split cells and INSIDE/OUTSIDE cell typing
  // then remove any particles now inside the moved body

  if (remapmode == CUTCELL) {
    grid_rebuild();
    if (particle->exist) ndeleted += remove_inside_particles(1);
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

  if (dim == 2) {
      if (fabs(ez_space[0]) > EPSILON || fabs(ez_space[1]) > EPSILON) {
        std::swap(inertia[1],inertia[2]);
        std::swap(ey_space[0],ez_space[0]);
        std::swap(ey_space[1],ez_space[1]);
        std::swap(ey_space[2],ez_space[2]);
      }
    }

  // if any principal moment < scaled EPSILON, set to 0.0

  double max;
  max = MAX(inertia[0],inertia[1]);
  max = MAX(max,inertia[2]);

  if (inertia[0] < EPSILON*max) inertia[0] = 0.0;
  if (inertia[1] < EPSILON*max) inertia[1] = 0.0;
  if (inertia[2] < EPSILON*max) inertia[2] = 0.0;

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

  // work arrays for remap of body surfs to grid cells

  memory->create(elemlo,nsurf,3,"fix_rigid:elemlo");
  memory->create(elemhi,nsurf,3,"fix_rigid:elemhi");
  memory->create(tmplist,nsurf,"fix_rigid:tmplist");

  // zero body force/torque in case accessed via compute_vector() on step 0

  fcm[0] = fcm[1] = fcm[2] = 0.0;
  torque[0] = torque[1] = torque[2] = 0.0;
  fpush[0] = fpush[1] = fpush[2] = 0.0;
}

/* ----------------------------------------------------------------------
   full re-map of surfs to grid cells
   same sequence of operations as in FixMoveSurf::end_of_step()
   for overlay mode: called once at fix creation,
     body surfs are excluded via Surf::rigidbits
   for cutcell mode: called every step after body surfs move,
     body surfs cut/split cells and set INSIDE/OUTSIDE cell typing
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
   overlay body surfs onto grid cells for the current timestep
   a body surf is added to the csurfs list of every owned or ghost cell
     which its swept bbox for this step overlaps
   cells settings are saved so overlay_restore() can undo the overlay
------------------------------------------------------------------------- */

void FixRigid::overlay_assign()
{
  int i,j,icell,isub,nstatic,nmov,isplit;
  surfint *merged;

  Grid::ChildCell *cells = grid->cells;
  Grid::SplitInfo *sinfo = grid->sinfo;
  int ntotal = grid->nlocal + grid->nghost;

  // swept bounding boxes of body elements over this timestep

  body_bbox(1);

  // reset page of merged surf lists and overlay restore lists

  cpage->reset();
  nmodified = 0;

  // loop over owned + ghost cells
  // skip sub cells, handled via their split cell
  // skip empty ghost cells, flagged with nsurf = -1

  for (icell = 0; icell < ntotal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    if (cells[icell].nsurf < 0) continue;
    if (!box_overlap(cells[icell].lo,cells[icell].hi,bbodylo,bbodyhi))
      continue;

    // nmov = # of body elements whose swept bbox overlaps this cell

    nmov = 0;
    for (i = 0; i < nsurf; i++)
      if (box_overlap(cells[icell].lo,cells[icell].hi,elemlo[i],elemhi[i]))
        tmplist[nmov++] = (surfint) slist[i];
    if (!nmov) continue;

    // merged surf list = cell's static surfs + overlaid body elements

    nstatic = cells[icell].nsurf;
    merged = cpage->get(nstatic+nmov);
    for (j = 0; j < nstatic; j++) merged[j] = cells[icell].csurfs[j];
    for (j = 0; j < nmov; j++) merged[nstatic+j] = tmplist[j];

    // save cell settings so they can be restored, then override them
    // for a split cell, also override its sub cells,
    //   which share nsurf/csurfs with their split cell

    if (nmodified+cells[icell].nsplit > maxmodified) {
      maxmodified += DELTA_MODIFY;
      memory->grow(modified,maxmodified,"fix_rigid:modified");
      memory->grow(nsurf_saved,maxmodified,"fix_rigid:nsurf_saved");
      csurfs_saved = (surfint **)
        memory->srealloc(csurfs_saved,maxmodified*sizeof(surfint *),
                         "fix_rigid:csurfs_saved");
    }

    modified[nmodified] = icell;
    nsurf_saved[nmodified] = cells[icell].nsurf;
    csurfs_saved[nmodified] = cells[icell].csurfs;
    nmodified++;
    cells[icell].nsurf = nstatic + nmov;
    cells[icell].csurfs = merged;

    if (cells[icell].nsplit > 1) {
      isplit = cells[icell].isplit;
      for (j = 0; j < cells[icell].nsplit; j++) {
        isub = sinfo[isplit].csubs[j];
        modified[nmodified] = isub;
        nsurf_saved[nmodified] = cells[isub].nsurf;
        csurfs_saved[nmodified] = cells[isub].csurfs;
        nmodified++;
        cells[isub].nsurf = nstatic + nmov;
        cells[isub].csurfs = merged;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   restore static surf assignment of all grid cells overlaid this step
------------------------------------------------------------------------- */

void FixRigid::overlay_restore()
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
   grid cells were rebuilt or migrated to other procs
   any overlaid csurfs lists were discarded by the grid rebuild,
     so just invalidate the overlay restore records
   next start_of_step() re-assigns body surfs to grid cells
------------------------------------------------------------------------- */

void FixRigid::grid_changed()
{
  nmodified = 0;
  if (cpage) cpage->reset();
}

/* ----------------------------------------------------------------------
   compute bounding boxes around each body element and around whole body
   sweepflag = 0: boxes bound elements at their current positions
   sweepflag = 1: boxes also bound elements at their end-of-step positions
     computed from xcmnew and current exyz_space (set from quatnew),
     so only valid after time integration in start_of_step()
   boxes are inflated by EPSSURF * body extent to avoid round-off misses
------------------------------------------------------------------------- */

void FixRigid::body_bbox(int sweepflag)
{
  int i,j,k,index;
  double *pts[3];
  double delta[3],ptnew[3];
  double *lo,*hi;

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

    lo = elemlo[i];
    hi = elemhi[i];
    lo[0] = lo[1] = lo[2] = BIG;
    hi[0] = hi[1] = hi[2] = -BIG;

    for (j = 0; j < npoint; j++) {
      for (k = 0; k < 3; k++) {
        lo[k] = MIN(lo[k],pts[j][k]);
        hi[k] = MAX(hi[k],pts[j][k]);
      }
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

  for (i = 0; i < nsurf; i++) {
    for (k = 0; k < 3; k++) {
      elemlo[i][k] -= eps;
      elemhi[i][k] += eps;
    }
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
  bigint delta = nlocal_old - particle->nlocal;
  bigint nall;
  MPI_Allreduce(&delta,&nall,1,MPI_SPARTA_BIGINT,MPI_SUM,world);
  return nall;
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
