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

#ifdef COMMAND_CLASS

CommandStyle(create_particles,CreateParticles)

#else

#ifndef SPARTA_CREATE_PARTICLES_H
#define SPARTA_CREATE_PARTICLES_H

#include "pointers.h"

namespace SPARTA_NS {

class CreateParticles : protected Pointers {

 public:
  CreateParticles(class SPARTA *);
  void command(int, char **);
  int evib(int);
  double erot(int);

 protected:
  int imix,single,cutflag,mspecies,twopass;
  bigint np;
  double xp,yp,zp,vx,vy,vz;
  class Region *region;

  int speciesflag,densflag,velflag,tempflag,normflag;
  char *sstr,*sxstr,*systr,*szstr;
  char *dstr,*dxstr,*dystr,*dzstr;
  char *tstr,*txstr,*tystr,*tzstr;
  char *vxstr,*vystr,*vzstr,*vstrx,*vstry,*vstrz;
  int svar,sxvar,syvar,szvar;
  int dvar,dxvar,dyvar,dzvar;
  int tvar,txvar,tyvar,tzvar;
  int vxvar,vyvar,vzvar,vvarx,vvary,vvarz;
  char *sxstr_copy,*systr_copy,*szstr_copy;
  char *dxstr_copy,*dystr_copy,*dzstr_copy;
  char *txstr_copy,*tystr_copy,*tzstr_copy;
  char *vstrx_copy,*vstry_copy,*vstrz_copy;

  virtual void create_single();
  virtual void create_local();
  virtual void create_local_twopass();
  int species_variable(double *);
  double density_variable(double *, double *);
  double temperature_variable(double *);
  void velocity_variable(double *, double *, double *);
  int outside_region(int, double *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot create particles before simulation box is defined

Self-explanatory.

E: Cannot create particles before grid is defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Create_particles mixture ID does not exist

Self-explanatory.

E: Create_particles species ID does not exist

Self-explanatory.

E: Create_particles global option not yet implemented

Self-explanatory.

E: Created incorrect # of particles: %ld versus %ld

The create_particles command did not function
properly.

E: Create_particles single requires z = 0 for 2d simulation

Self-explanatory.

E: Could not create a single particle

The specified position was either not inside the simulation domain or
not inside a grid cell with no intersections with any defined surface
elements.

*/
