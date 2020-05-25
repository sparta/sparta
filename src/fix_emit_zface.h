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

FixStyle(emit/zface,FixEmitZFace)

#else

#ifndef SPARTA_FIX_EMIT_ZFACE_H
#define SPARTA_FIX_EMIT_ZFACE_H

#include "fix_emit.h"
#include "surf.h"
#include "grid.h"

#ifdef USE_ZSURF
#include "zuzax/zeroD/SurfPropagationSparta.h"
#endif

namespace SPARTA_NS {

class FixEmitZFace : public FixEmit {
 public:
  FixEmitZFace(class SPARTA *, int, char **);
  virtual ~FixEmitZFace();
  virtual void init();
  virtual void post_compress_grid();



  //! one insertion task for each cell/face that shares a global boundary
  //! with the specified face of the domain

  struct Task {
    double lo[3];               // lower-left corner of face
    double hi[3];               // upper-right corner of face
    double normal[3];           // inward normal from external boundary
    double area;                // area of cell face
    //double ntarget;             // # of mols to insert for all species
    double nrho;                // from mixture or adjacent subsonic cell
    double temp_thermal;        // from mixture or adjacent subsonic cell
    double temp_rot;            // from mixture or subsonic temp_thermal
    double temp_vib;            // from mixture or subsonic temp_thermal
    double vstream[3];          // from mixture or adjacent subsonic cell

    int icell;                  // associated cell index, unsplit or split cell
    int iface;                  // which face of unsplit or split cell
    int pcell;                  // associated cell index for particles
                                // unsplit or sub cell (not split cell)
    int ndim;                   // dim (0,1,2) normal to face
    int pdim,qdim;              // 2 dims (0,1,2) parallel to face
  };

 protected:
  int isrZuzax;       // reaction int for the Zuzax reaction that's occuring on the surface
  int imix,np;

#ifdef USE_ZSURF
  mutable Zuzax::SurfPropagationSparta* net {nullptr};
#endif

  double *vscale;     //vscale for each species, calculated at the surface temperature

  int faces[6];       // One of these 6 faces will be true.
  
  int iFaceReact;     // The one face that is reacting in this single fix
                      // -> singlingly this down to one face for the initial installation.
  double areaLocal;   // local sum of areas

  SurfState* ssFaceReact {nullptr};
  int npertask,nthresh;
  int twopass;        // True if we are doing a two pass algorithm
  double tprefactor,soundspeed_mixture;

  // copies of data from other classes

  int dimension;
  int nspecies;
  double fnum,dt;
  double *fraction,*cummulative;

  Surf::Line *lines;     // not sure we need this for this bc
  Surf::Tri *tris;       // not sure we need this for this bc

                         // ntask = # of tasks is stored by parent class
  Task *tasks;           // list of particle insertion tasks
  int ntaskmax;          // max # of tasks allocated

  // active grid cells assigned to tasks, used by subsonic sorting

  int maxactive;
  int *activecell;

  // protected methods

  //! The routine, create_task() is run during the sparta initialization routines to
  //! create a placeholder task for each cell/face combination that needs a fix process
  //! to be run.
  /*!
   *  @return              Returns the # of tasks for this cell created
   */
  virtual void create_task(int icell) override;

  //! Perform_task() gets called during the modify->start_of_step() process
  /*!
   *  This loops over the previously created tasks, which are cell/face specific.
   */
  virtual void perform_task() override;

  void perform_task_onepass();

  virtual void perform_task_twopass();

  int split(int, int);


  virtual int setmask() override;

  //! Carry out final steps for Surface here
  /*!
   *
   */
  virtual void end_of_step() override;

  void getGlobalReactionEventsFromLocalEvents();



  virtual int pack_task(int, char *, int);
  virtual int unpack_task(char *, int);
  virtual void copy_task(int, int, int, int);
  void grow_task();
  virtual void realloc_nspecies();

  int option(int, char **);
};

}

#endif
#endif

