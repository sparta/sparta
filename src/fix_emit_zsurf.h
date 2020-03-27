/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

#ifdef USE_ZSURF
FixStyle(emit/zsurf,FixEmitZSurf)
#endif

#else

#ifndef SPARTA_FIX_EMIT_ZSURF_H
#define SPARTA_FIX_EMIT_ZSURF_H




#include "fix_emit.h"
#include "surf.h"
#include "grid.h"

#ifdef USE_ZSURF
#include "zuzax/zeroD/SurfPropagationSparta.h"
#endif

namespace SPARTA_NS {

#ifdef USE_ZSURF

class FixEmitZSurf : public FixEmit {
 public:
  FixEmitZSurf(class SPARTA *, int, char **);
  ~FixEmitZSurf();
  void init();
  void setup();
  void post_compress_grid();

 private:
  int imix,groupbit,np,normalflag;
  int npertask,nthresh;
  double tprefactor,soundspeed_mixture;

  // copies of data from other classes

  int dimension,nspecies;
  double fnum,dt;
  double nrho,temp_thermal,temp_rot,temp_vib;
  double *fraction,*cummulative;

  Surf::Line *lines {nullptr};
  Surf::Tri *tris {nullptr};

  Zuzax::SurfPropagationSparta* net {nullptr};

  double* vscale {nullptr};    // vscale for each species, 
                               // evaluated at surface temperature
                               // which is obtained from surface state.

  class Cut2d *cut2d {nullptr};
  class Cut3d *cut3d {nullptr};

  //! The Task structure holds one insertion task for a single cell and a surf
  //! combination
  /*!
   *  The particle insertion task will consist of inserting all particles needed
   *  for implementation of reactions initiated by the surface by its itself,
   *  e.g. adsorbate dissociation, sublimation events.
   */
  struct Task {
    double area;                // area of overlap of surf with cell
    double tan1[3],tan2[3];     // 2 normalized tangent vectors to surf normal
    double nrho;                // from mixture or adjacent subsonic cell
    double temp_thermal;        // from mixture or adjacent subsonic cell
    double temp_rot;            // from mixture or subsonic temp_thermal
    double temp_vib;            // from mixture or subsonic temp_thermal
    double vstream[3];          // from mixture or adjacent subsonic cell
    double *path;               // path of points for overlap of surf with cell
    double *fracarea;           // fractional area for each sub tri in path

    int icell;                  // associated cell index, unsplit or split cell
    int isurf;                  // surf index
    int pcell;                  // associated cell index for particles
                                // unsplit or sub cell (not split cell)
    int npoint;                 // # of points in path
  };

  int isrZuzax;          // int representing the Zuzax reaction mechanism 
                         // assigned to this surface

  Task *tasks;           // List of particle insertion tasks.
                         
  int ntaskmax;          // max # of tasks allocated

  double magvstream;       // magnitude of mixture vstream
  double norm_vstream[3];  // direction of mixture vstream

  // active grid cells assigned to tasks, used by subsonic sorting

  int maxactive;
  int *activecell;

  //! Create a list of tasks to perform each time step
  /*!
   *  (virtual from FixEmit)
   */
  virtual void create_task(int) override;

  //! Perform the list of tasks to perform each time step
  /*!
   *  (virtual from FixEmit)
   */
  virtual void perform_task() override;

  //! Carry out final steps for Surface here
  /*!
   *
   */
  virtual void end_of_step() override;

  virtual int setmask() override;  

  virtual int pack_task(int, char *, int);
  virtual int unpack_task(char *, int);
  virtual void copy_task(int, int, int, int);
  void grow_task();

  int option(int, char **);
};
#endif

}

#endif
#endif


