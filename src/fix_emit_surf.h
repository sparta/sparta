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

#ifdef FIX_CLASS

FixStyle(emit/surf,FixEmitSurf)

#else

#ifndef SPARTA_FIX_EMIT_SURF_H
#define SPARTA_FIX_EMIT_SURF_H

#include "fix_emit.h"
#include "surf.h"
#include "grid.h"

namespace SPARTA_NS {

class FixEmitSurf : public FixEmit {
 public:
  FixEmitSurf(class SPARTA *, int, char **);
  ~FixEmitSurf();
  void init();

  void grid_changed();

 private:
  int imix,groupbit,normalflag,subsonic,subsonic_style,subsonic_warning;
  int npertask,nthresh;
  double psubsonic,tsubsonic,nsubsonic;
  double tprefactor,soundspeed_mixture;

  int npmode,np;    // npmode = FLOW,CONSTANT,VARIABLE
  int npvar;
  char *npstr;

  // copies of data from other classes

  int dimension,nspecies;
  double fnum,dt;
  double nrho,temp_thermal,temp_rot,temp_vib;
  double *fraction,*cummulative;

  class Cut2d *cut2d;
  class Cut3d *cut3d;

  // one insertion task for a cell and a surf

  struct Task {
    double area;                // area of overlap of surf with cell
    double ntarget;             // # of mols to insert for all species
    double tan1[3],tan2[3];     // 2 normalized tangent vectors to surf normal
    double nrho;                // from mixture or adjacent subsonic cell
    double temp_thermal;        // from mixture or adjacent subsonic cell
    double temp_rot;            // from mixture or subsonic temp_thermal
    double temp_vib;            // from mixture or subsonic temp_thermal
    double vstream[3];          // from mixture or adjacent subsonic cell
    double *ntargetsp;          // # of mols to insert for each species,
                                //   only defined for PERSPECIES
    double *vscale;             // vscale for each species,
                                //   only defined for subsonic_style PONLY

    double *path;               // path of points for overlap of surf with cell
    double *fracarea;           // fractional area for each sub tri in path

    int icell;                  // associated cell index, unsplit or split cell
    surfint isurf;              // surf index, sometimes a surf ID
    int pcell;                  // associated cell index for particles
                                // unsplit or sub cell (not split cell)
    int npoint;                 // # of points in path
  };

                         // ntask = # of tasks is stored by parent class
  Task *tasks;           // list of particle insertion tasks
  int ntaskmax;          // max # of tasks allocated

  double magvstream;       // magnitude of mixture vstream
  double norm_vstream[3];  // direction of mixture vstream

  // active grid cells assigned to tasks, used by subsonic sorting

  int maxactive;
  int *activecell;

  // private methods

  void create_task(int);
  void perform_task();
  void grow_task();

  void subsonic_inflow();
  void subsonic_sort();
  void subsonic_grid();

  int option(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Fix emit/face mixture ID does not exist

Self-explanatory.

E: Cannot use fix emit/face in z dimension for 2d simulation

Self-explanatory.

E: Cannot use fix emit/face in y dimension for axisymmetric

This is because the y dimension boundaries cannot be
emit/face boundaries for an axisymmetric model.

E: Cannot use fix emit/face n > 0 with perspecies yes

This is because the perspecies option calculates the
number of particles to insert itself.

E: Cannot use fix emit/face on periodic boundary

Self-explanatory.

E: Fix emit/face used on outflow boundary

Self-explanatory.

*/
