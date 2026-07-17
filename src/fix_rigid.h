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
#include <map>

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
                          // non-distributed surfs only, else NULL
  int nsurfall;           // length of irigid = surf->nlocal when allocated
                          // Update::init() clamps its scan to this length

  // body_elem() = element index of a global surf ID, -1 if not in body
  // used by Update::build_rigidmap() for distributed surfs

  int body_elem(surfint id)
  {
    std::map<surfint,int>::iterator it = idmap.find(id);
    if (it == idmap.end()) return -1;
    return it->second;
  }

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
  int *slist;    // list of local surf indices for rigid body surfs
                 //   non-distributed surfs only, else NULL

  // replicated body geometry, the authoritative source for all body
  //   computations (bbox, inside tests, watertight, contacts)
  // for distributed surfs it is gathered from the owned copies at
  //   setup and regenerated from the body pose each step

  double ***bodypt;       // bodypt[i][j] = corner pt j of element i
  double **bodynorm;      // outward normal of element i
  surfint *sids;          // global surf ID of each element
  int *bodymask;          // per-element group mask
  int *bodytype;          // per-element surf type
  int *bodytrans;         // per-element transparent flag
  int *bodyisc;           // per-element collision model index
  int *bodyisr;           // per-element reaction model index
  std::map<surfint,int> idmap;   // global surf ID -> element index

  // distributed surfs: where this proc stores copies of body elements

  int *lblist;            // local surf index of each element, -1 if none
  int nolist;             // # of body elements this proc owns
  int *olist_own;         // owned-array index of each
  int *olist_elem;        // element index of each

  double rmaxbody;        // max distance of any body corner pt from COM
  double mincellsize;     // smallest edge length of any grid cell
  int warnrotate;         // 1 after warning about rotation rate
  int warntranslate;      // 1 after warning about translation rate
  int warnexit;           // 1 after warning that body exited the box

  double tqpush[3];       // torque from push-off contacts this step
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

  // bins over static surfs for push-off candidate pruning
  // built once per run in setup(), static surfs never move during a run

  int pushnbin[3];        // # of bins in each dim
  double pushbinlo[3];    // bin grid origin
  double pushbininv[3];   // inverse bin edge lengths
  int *pushbinstart;      // CSR offsets into pushbinlist per bin
  int *pushbinlist;       // static surf indices, binned by surf bbox
  int *pushstamp;         // per-surf visit stamp to dedup multi-bin surfs
  int pushstampcur;

  // work buffers for the fused force/torque Allreduce over all bodies

  double *ftbuf_mine;
  double *ftbuf_all;

  // remap of body surfs to grid cells

  int remapmode;          // CUTCELL or INCREMENTAL

  // swept collision assignment: each step every body's surfs are added
  //   to the collision lists (csurfs) of all cells they sweep through,
  //   so particles anywhere in a body's swept path are tested against
  //   the moving surfs and reflected rather than overtaken and deleted
  // this augments collision lists only; cut-cell volumes are unaffected
  // the single pass over grid cells for ALL bodies is performed by the
  //   last-defined rigid fix, which owns the save/restore bookkeeping

  double **elemlo;        // per-element bounding boxes: swept boxes for
  double **elemhi;        //   this step, or current boxes after commit

  int nmodified;          // # of cells augmented this step
  int maxmodified;        // allocated size of restore lists
  int *modified;          // indices of cells whose csurfs were augmented
  int *nsurf_saved;       // saved nsurf of each modified cell
  surfint **csurfs_saved; // saved csurfs ptr of each modified cell
  MyPage<surfint> *cpage; // storage for merged csurfs lists

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

  int nrcand;             // work list of cells overlapping the
  int maxrcand;           //   incremental re-cut region this step
  int *rcand;

  surfint *newlist;       // work bufs for re-cutting one cell
  int *newmap;
  class Cut2d *cut2d;
  class Cut3d *cut3d;
  double bbodylo[3];      // bounding box around entire body
  double bbodyhi[3];

  bigint ndeleted;        // per-proc count of deleted particles
  bigint ndeleted_all;    // cached global sum for compute_scalar()
  bigint ndelvalid;       // timestep the cached sum is valid for

  void read_infile(char *);
  void write_outfile();
  void setup_body();
  void check_watertight();

  void push_off();              // spring forces on this body, with
                                //   equal-opposite reactions on others
  void push_contact(double *, double *, double *, double *,
                    class FixRigid *);  // corner contacts vs one source elem
  void push_bins();             // bin static surfs for candidate pruning
  void gather_body();           // build replicated body element table
  void ensure_local_copies();   // distributed: local copies of body surfs
  void update_surf_copies();    // write bodypt/bodynorm into Surf storage
  void grid_rebuild();          // full re-map of all surfs to grid cells
  void record_oldinside();      // cells interior to body, pre-move
  int incremental_recut();      // re-cut cells whose overlap changed
  void registry_replace(int, surfint *);
  void registry_remove(int);
  void free_registry();
  void swept_assign_all();      // add all bodies' surfs to swept cells
  void swept_restore();         // undo swept_assign_all
  void body_bbox(int);          // bbox of body, current or swept over step
  int inside_body(double *);    // 1 if point is inside rigid body, else 0
  int inside_any_body(double *); // 1 if inside any rigid body
  bigint remove_inside_particles(int);  // per-body, used at setup
  void remove_inside_all(int);  // fused pass over all bodies, per step
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
