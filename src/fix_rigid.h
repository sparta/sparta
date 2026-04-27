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

#ifdef FIX_CLASS

FixStyle(rigid,FixRigid)

#else

#ifndef SPARTA_FIX_RIGID_H
#define SPARTA_FIX_RIGID_H

#include "fix.h"

namespace SPARTA_NS {

class FixRigid : public Fix {
 public:
  FixRigid(class SPARTA *, int, char **);
  virtual ~FixRigid() {}
  int setmask();
  void init();
  virtual void start_of_step();
  virtual void end_of_step();
  double compute_vector(int);

 protected:
  int igroup,groupbit;
  char *csurfID;
  int icompute;
  class Compute *csurf;
  char *infile;

  int nsurf;
  int *slist;
  
  int dim;
  int massflag,comflag,vcomflag,moiflag,angmomflag;

  int nparticleflag,pmassflag;
  int nparticle_user;
  double pmass_user;
  
  double massbody;
  double xcm[3],vcm[3];
  double moi[6];      // space frame
  double inertia[3];  // body frame
  double angmom[3];
  double quat[4];
  double ex_space[3],ey_space[3],ez_space[3];
  double fcm[3];
  double torque[3];
  double omega[3];

  double xcmnew[3];
  double quatnew[4];

  double ***displace;
  
  void read_infile(char *);
  void setup_body();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
