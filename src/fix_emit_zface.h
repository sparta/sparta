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

namespace SPARTA_NS {

class FixEmitZFace : public FixEmit {
 public:
  FixEmitZFace(class SPARTA *, int, char **);
  virtual ~FixEmitZFace();
  virtual void init();
  virtual void post_compress_grid();

  // one insertion task for a cell and a face

  struct Task {
    double lo[3];               // lower-left corner of face
    double hi[3];               // upper-right corner of face
    double normal[3];           // inward normal from external boundary
    double area;                // area of cell face
    double ntarget;             // # of mols to insert for all species
    double nrho;                // from mixture or adjacent subsonic cell
    double temp_thermal;        // from mixture or adjacent subsonic cell
    double temp_rot;            // from mixture or subsonic temp_thermal
    double temp_vib;            // from mixture or subsonic temp_thermal
    double vstream[3];          // from mixture or adjacent subsonic cell
    double *ntargetsp;          // # of mols to insert for each species,
                                //   only defined for PERSPECIES
    double *vscale;             // vscale for each species,
                                //   only defined for subsonic_style PONLY

    int icell;                  // associated cell index, unsplit or split cell
    int iface;                  // which face of unsplit or split cell
    int pcell;                  // associated cell index for particles
                                // unsplit or sub cell (not split cell)
    int ndim;                   // dim (0,1,2) normal to face
    int pdim,qdim;              // 2 dims (0,1,2) parallel to face
  };

 protected:
  int isrZuzax;
  int imix,np,subsonic,subsonic_style,subsonic_warning;
  int faces[6];
  int npertask,nthresh,twopass;
  double psubsonic,tsubsonic,nsubsonic;
  double tprefactor,soundspeed_mixture;

  // copies of data from other classes

  int dimension,nspecies;
  double fnum,dt;
  double *fraction,*cummulative;

  Surf::Line *lines;
  Surf::Tri *tris;

                         // ntask = # of tasks is stored by parent class
  Task *tasks;           // list of particle insertion tasks
  int ntaskmax;          // max # of tasks allocated

  // active grid cells assigned to tasks, used by subsonic sorting

  int maxactive;
  int *activecell;

  // protected methods

  int create_task(int);
  virtual void perform_task();
  void perform_task_onepass();
  virtual void perform_task_twopass();

  int split(int, int);

  void subsonic_inflow();
  void subsonic_sort();
  void subsonic_grid();

  virtual int pack_task(int, char *, int);
  virtual int unpack_task(char *, int);
  virtual void copy_task(int, int, int, int);
  virtual void grow_task();
  virtual void realloc_nspecies();

  int option(int, char **);
};

}

#endif
#endif

