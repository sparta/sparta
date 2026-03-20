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
#include "modify.h"
#include "compute.h"
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

// DEBUG
enum{OUTSIDE,INSIDE,ONSURF2OUT,ONSURF2IN};    // same as Update

/* ---------------------------------------------------------------------- */

FixRigid::FixRigid(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix rigid command");

  vector_flag = 1;
  size_vector = 12;
  global_freq = 1;
  nevery = 1;
  
  if (!surf->exist) error->all(FLERR,"Fix rigid requiers surf elements exist");
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
	ixx = input->numeric(FLERR,arg[jarg+1]);
	iyy = input->numeric(FLERR,arg[jarg+2]);
	izz = input->numeric(FLERR,arg[jarg+3]);
	ixy = input->numeric(FLERR,arg[jarg+4]);
	ixz = input->numeric(FLERR,arg[jarg+5]);
	iyz = input->numeric(FLERR,arg[jarg+6]);
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
    
  } else error->all(FLERR,"Fix rigid define not recognized");

  // optional args
  
  nparticleflag = 0;
  pmassflag = 0;
  double scale = 1.0;
  
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nparticle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Fix rigid body args not valid");
      nparticleflag = 1;
      nparticle_user = input->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"pmass") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Fix rigid body args not valid");
      pmassflag = 1;
      pmass_user = input->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Fix rigid body args not valid");
      scale = input->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Fix rigid body args not valid");
  }

  // apply scaling to either body or infile
  // NOTE: how to apply to group of surfs themselves
  //       allow this to be done by read_surf ?  rotate is hard
  
  // setup the rigid body

  setup_body();

  printf("EX %g %g %g\n",ex_space[0],ex_space[1],ex_space[2]);
  printf("EY %g %g %g\n",ey_space[0],ey_space[1],ey_space[2]);
  printf("EZ %g %g %g\n",ez_space[0],ez_space[1],ez_space[2]);
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
  // check that specified compute is valid for use with fix rigid
  // NOTE: how to check it operates on same surf group ?
  
  int n = modify->find_compute(csurfID);
  if (n < 0) error->all(FLERR,"Could not find fix rigid compute ID");
  csurf = modify->compute[n];
  if (strcmp(csurf->style,"surf") != 0)
    error->all(FLERR,"Fix rigid compute is not style surf");
  if (csurf->per_surf_flag == 0)
    error->all(FLERR,"Fix rigid compute does not compute per-surf info");
  if (csurf->size_per_surf_cols != 6)
    error->all(FLERR,"Fix rigid compute must calcualte per-surf array with 6 columns");
}

/* ---------------------------------------------------------------------- */

void FixRigid::start_of_step()
{
  // TODO:
  // update RB props
  // change COM within compute surf - note that COM is not constant thru step
  // remap RB surfs to grid cells they cover during advection

  // time integrate from current position to end-of-step position
  // use full-step semi-implicit Euler algorithm
  // apply forces and torques accumulated from collisions during last step
  // vcm/angmmom/omega are updated now, b/c they are used during
  //   this step to update xnew/quat at end of this step
  // xcmnew/quatnew will become position/orientation at end of this step
  // save old and new xcm/quat to use in particle collision detection
  
  double dt = update->dt;
  double dtf = dt / massbody;
  double dthalf = 0.5 * dt;
  
  // update vcm by full step

  vcm[0] += dtf * fcm[0];
  vcm[1] += dtf * fcm[1];
  vcm[2] += dtf * fcm[2];

  // update xcm by full step
  // using new vcm turns Euler into semi-implicit Euler

  xcmnew[0] = xcm[0] + dt * vcm[0];
  xcmnew[1] = xcm[1] + dt * vcm[1];
  xcmnew[2] = xcm[2] + dt * vcm[2];

  // update angular momentum in spatial frame by full step

  angmom[0] += dt * torque[0];
  angmom[1] += dt * torque[1];
  angmom[2] += dt * torque[2];

  // compute new omega from new angmom, both in spatial frame
  
  MathExtra::angmom_to_omega(angmom,ex_space,ey_space,ez_space,inertia,omega);

  // update quaternion by full step using new omega in spatial frame
  // use dq/dt = 1/2 omega q

  double wq[4];
  MathExtra::vecquat(omega,quat,wq);
  quatnew[0] = quat[0] + dthalf * wq[0];
  quatnew[1] = quat[1] + dthalf * wq[1];
  quatnew[2] = quat[2] + dthalf * wq[2];
  quatnew[3] = quat[3] + dthalf * wq[3];
  MathExtra::qnormalize(quatnew);
  // NOTE: need to also keep old exyz_space ??  for collision detection ??
  MathExtra::q_to_exyz(quatnew,ex_space,ey_space,ez_space);
}

/* ---------------------------------------------------------------------- */

void FixRigid::end_of_step()
{
  // invoke compute surf and extract per-surf force/torque info
  // NOTE: could this access a fix ave/surf for many steps of collisions ?

  /*
  if (!(csurf->invoked_flag & INVOKED_PER_SURF)) {
    csurf->compute_per_surf();
    csurf->invoked_flag |= INVOKED_PER_SURF;
  }

  csurf->post_process_surf();
  double **array = csurf->array_surf;

  // sum per-surf force/torque to body fcm/torque

  fcm[0] = fcm[1] = fcm[2] = 0.0;
  torque[0] = torque[1] = torque[2] = 0.0;

  int index;
  
  for (int i = 0; i < nsurf; i++) {
    index = slist[i];
    fcm[0] += array[index][0];
    fcm[1] += array[index][1];
    fcm[2] += array[index][2];
    torque[0] += array[index][3];
    torque[1] += array[index][4];
    torque[2] += array[index][5];
  }
  */

  // DEBUG: bounce N fictitious particles off object, see how it moves
  // bbox = epsilon-augmented bbox around current object
  // shoot N particles from random positions from left face to right face of bbox
  // whichever surf it hits first, compute force/torque on object

  int nparticle = 10;
  double pmass = 1.0;

  if (nparticleflag) nparticle = nparticle_user;
  if (pmassflag) pmass = pmass_user;

  RanKnuth *random = new RanKnuth(update->ranmaster->uniform());

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

    printf("BBOX x %g %g y %g %g z %g %g\n",
	 bboxlo[0],bboxhi[0],
	 bboxlo[1],bboxhi[1],
	 bboxlo[2],bboxhi[2]);

  bboxlo[0] -= 0.01 * (bboxhi[0]-bboxlo[0]);
  bboxlo[1] -= 0.01 * (bboxhi[1]-bboxlo[1]);
  bboxlo[2] -= 0.01 * (bboxhi[2]-bboxlo[2]);
  bboxhi[0] += 0.01 * (bboxhi[0]-bboxlo[0]);
  bboxhi[1] += 0.01 * (bboxhi[1]-bboxlo[1]);
  bboxhi[2] += 0.01 * (bboxhi[2]-bboxlo[2]);

  printf("BBOX x %g %g y %g %g z %g %g\n",
	 bboxlo[0],bboxhi[0],
	 bboxlo[1],bboxhi[1],
	 bboxlo[2],bboxhi[2]);

  int nhits = 0;
  
  if (dim == 2) {
    for (int m = 0; m < nparticle; m++) {
      int side,minsurf;
      double param;
      double x[3],xnew[3];
      double xc[3],minxc[3];

      // set x,xnew randomly for each of N particles

      x[0] = bboxlo[0];
      xnew[0] = bboxhi[0];
      double rn = random->uniform();
      // entire face
      //xnew[1] = x[1] = bboxlo[1] + rn * (bboxhi[1]-bboxlo[1]);
      // lower half
      xnew[1] = x[1] = bboxlo[1] + 0.5 * rn * (bboxhi[1]-bboxlo[1]);
      // lower 10%
      //xnew[1] = x[1] = bboxlo[1] + 0.1 * rn * (bboxhi[1]-bboxlo[1]);
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
        MathExtra::axpy3(pmass,vpre,pforce);
        MathExtra::axpy3(-pmass,vpost,pforce);
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
    for (int m = 0; m < nparticle; m++) {
      int side,minsurf;
      double param;
      double x[3],xnew[3];
      double xc[3],minxc[3];

      // set x,xnew randomly for each of N particles
      
      x[0] = bboxlo[0];
      xnew[0] = bboxhi[0];
      double rn = random->uniform();
      xnew[1] = x[1] = bboxlo[1] + rn * (bboxhi[1]-bboxlo[1]);
      rn = random->uniform();
      xnew[2] = x[2] = bboxlo[2] + rn * (bboxhi[2]-bboxlo[2]);

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
        MathExtra::axpy3(pmass,vpre,pforce);
        MathExtra::axpy3(-pmass,vpost,pforce);
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

  printf("F/T %ld hits %d\n",update->ntimestep,nhits);
  printf("FCM %g %g %g\n",fcm[0],fcm[1],fcm[2]);
  printf("TQ %g %g %g\n",torque[0],torque[1],torque[2]);
    
  // reset xcm/quat to new xcm/quat from start_of_step()

  xcm[0] = xcmnew[0];
  xcm[1] = xcmnew[1];
  xcm[2] = xcmnew[2];

  quat[0] = quatnew[0];
  quat[1] = quatnew[1];
  quat[2] = quatnew[2];
  quat[3] = quatnew[3];

  // enforce2d on all body properties

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

  printf("XCM %g %g %g\n",xcm[0],xcm[1],xcm[2]);
  printf("VCM %g %g %g\n",vcm[0],vcm[1],vcm[2]);
  printf("ANG %g %g %g\n",angmom[0],angmom[1],angmom[2]);
  printf("OMG %g %g %g\n",omega[0],omega[1],omega[2]);

  // update body line and tri positions and orientations in Surf class
  // set via Line/Tri end/corner points and norm
  // matvec() converts displace vector from body frame to space frame
  
  //Surf::Line *lines = surf->lines;
  //Surf::Tri *tris = surf->tris;

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

  /*
  printf("LINE1 pt1 %g %g pt2 %g %g NORM %g %g %g\n",
	 lines[0].p1[0],lines[0].p1[1],
	 lines[0].p2[0],lines[0].p2[1],
	 lines[0].norm[0],lines[0].norm[1],lines[0].norm[2]);
  printf("LINE2 pt1 %g %g pt2 %g %g NORM %g %g %g\n",
	 lines[1].p1[0],lines[1].p1[1],
	 lines[1].p2[0],lines[1].p2[1],
	 lines[1].norm[0],lines[1].norm[1],lines[1].norm[2]);
  printf("LINE3 pt1 %g %g pt2 %g %g NORM %g %g %g\n",
	 lines[2].p1[0],lines[2].p1[1],
	 lines[2].p2[0],lines[2].p2[1],
	 lines[2].norm[0],lines[2].norm[1],lines[2].norm[2]);
  printf("LINE4 pt1 %g %g pt2 %g %g NORM %g %g %g\n",
	 lines[3].p1[0],lines[3].p1[1],
	 lines[3].p2[0],lines[3].p2[1],
	 lines[3].norm[0],lines[3].norm[1],lines[3].norm[2]);
  */
  
  printf("END %ld\n",update->ntimestep);
}

/* ----------------------------------------------------------------------
   one-time initialization of rigid body attributes from file
------------------------------------------------------------------------- */

void FixRigid::read_infile(char *filename)
{
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

  tensor[0][0] = ixx;
  tensor[1][1] = iyy;
  tensor[2][2] = izz;
  tensor[1][2] = tensor[2][1] = iyz;
  tensor[0][2] = tensor[2][0] = ixz;
  tensor[0][1] = tensor[1][0] = ixy;

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

  // store displacement of each end/corner point in each line/tri in body
  // delta = vector from COM to end/corner point
  // displace = delta rotated to be in basis of principal axes

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

  // zero body force/torque in case accessed via compute_vector() on step 0

  fcm[0] = fcm[1] = fcm[2] = 0.0;
  torque[0] = torque[1] = torque[2] = 0.0;
}

/* ----------------------------------------------------------------------
   return properties of the single rigid body
------------------------------------------------------------------------- */

double FixRigid::compute_vector(int index)
{
  if (index == 0) return xcm[0];
  if (index == 1) return xcm[1];
  if (index == 2) return xcm[2];
  
  if (index == 3) return vcm[0];
  if (index == 4) return vcm[1];
  if (index == 5) return vcm[2];
  
  if (index == 6) return fcm[0];
  if (index == 7) return fcm[1];
  if (index == 8) return fcm[2];
  
  if (index == 9) return torque[0];
  if (index == 10) return torque[1];
  if (index == 11) return torque[2];

  return 0.0;
}
