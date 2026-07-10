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
#include "my_page.h"

namespace SPARTA_NS {

class FixRigid : public Fix {
 public:
  double omega[3];        // omega in space frame
  double ***displace;     // displacement in body frame of line/tri points from COM

  double xcm[3],vcm[3];   // COM and velocity of COM
  double quat[4];         // quaternion for orientation of body

  double xcmnew[3];       // new COM at end of timestep
  double quatnew[4];      // new quaternion at end of timestep

  int *irigid;            // per-surf flags, indexed by local surf index:
                          // -1 = static surf,
                          // else index into body's list of surfs

  FixRigid(class SPARTA *, int, char **);
  virtual ~FixRigid();
  int setmask();
  void init();
  void setup();
  virtual void start_of_step();
  virtual void end_of_step();
  void grid_changed();
  double compute_scalar();
  double compute_vector(int);

 protected:
  int igroup,groupbit;
  char *csurfID;
  class ComputeSurf *csurf;
  char *infile;
  char *outfile;
  int outevery;

  int dim;
  int massflag,comflag,vcomflag,moiflag,angmomflag;

  int pseudoflag;
  int nparticle_user;
  double pmass_user,frac_user;
  class RanKnuth *random;   // RNG for pseudo particles

  int nsurf;     // # of surfs which comprise surface of rigid body
  int *slist;    // list of surf indices for rigid body surfs

  double rmaxbody;        // max distance of any body corner pt from COM
  double mincellsize;     // smallest edge length of any grid cell
  int warnrotate;         // 1 after warning about rotation rate
  int warntranslate;      // 1 after warning about translation rate
  int warnexit;           // 1 after warning that body exited the box

  double massbody;        // total mass of rigid body enclosed by surfs
  double moi[6];          // 6 MOI in space frame
  double inertia[3];      // 3 diagonalized MOI in body frame
  double angmom[3];       // angular momentum in space frame
  double ex_space[3],ey_space[3],ez_space[3];  // prinicpal axes of body
  double fcm[3];          // force on COM in space frame
  double torque[3];       // torque on body in space frame
  double fpush[3];        // push-off force on COM from static surfs

  int pushflag;           // 1 if push-off forces are enabled
  int pushboundflag;      // 1 to also push off non-periodic boundaries
  int pushstyle;          // LINEAR or HERTZ force law
  double kpush;           // spring constant for push-off force
  double pushcutoff;      // distance below which push-off is applied
  double gammapush;       // dashpot damping coefficient, 0 = elastic

  // remap of body surfs to grid cells

  int remapmode;          // OVERLAY or CUTCELL

  int nmodified;          // # of cells with body surfs overlaid this step
  int maxmodified;        // allocated size of overlay restore lists
  int *modified;          // indices of cells whose csurfs were overridden
  int *nsurf_saved;       // saved static nsurf of each modified cell
  surfint **csurfs_saved; // saved static csurfs ptr of each modified cell
  MyPage<surfint> *cpage; // storage for merged static+body surf lists

  double **elemlo;        // per-element swept bounding boxes for this step
  double **elemhi;
  surfint *tmplist;       // work list of body elements overlapping one cell

  // incremental cutcell remap data

  double pbodylo[3];      // bbox around body at end of previous step
  double pbodyhi[3];
  int pbodyflag;          // 1 if pbodylo/pbodyhi are set

  int noldinside;         // cells interior to the body before it moved
  int maxoldinside;
  int *oldinside;

  int nreg;               // registry of cells whose csurfs lists
  int maxreg;             //   are allocated by this fix
  int *regcell;
  surfint **reglist;

  surfint *newlist;       // work bufs for re-cutting one cell
  int *newmap;
  class Cut2d *cut2d;
  class Cut3d *cut3d;
  double bbodylo[3];      // bounding box around entire body,
  double bbodyhi[3];      //   swept over the current step

  bigint ndeleted;        // running count of deleted particles

  void read_infile(char *);
  void write_outfile();
  void setup_body();
  void check_watertight();

  void push_off();              // spring forces from close static surfs
  void grid_rebuild();          // full re-map of all surfs to grid cells
  void record_oldinside();      // cells interior to body, pre-move
  int incremental_recut();      // re-cut cells whose overlap changed
  void registry_replace(int, surfint *);
  void registry_remove(int);
  void free_registry();
  void overlay_assign();        // per-step overlay of body surfs into cells
  void overlay_restore();       // undo overlay of previous step
  void body_bbox(int);          // swept bbox of body elements over a step
  int inside_body(double *);    // 1 if point is inside rigid body, else 0
  bigint remove_inside_particles(int);
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
