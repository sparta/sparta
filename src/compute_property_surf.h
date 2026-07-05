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

#ifdef COMPUTE_CLASS

ComputeStyle(property/surf,ComputePropertySurf)

#else

#ifndef SPARTA_COMPUTE_PROPERTY_SURF_H
#define SPARTA_COMPUTE_PROPERTY_SURF_H

#include "compute.h"

namespace SPARTA_NS {

class ComputePropertySurf : public Compute {
 public:
  ComputePropertySurf(class SPARTA *, int, char **);
  ~ComputePropertySurf();
  void init();
  void compute_per_surf();
  bigint memory_usage();

 protected:
  int groupbit,nvalues;
  int distributed;

  int dimension;

  int firstflag;             // 1 until cglobal is setup
  int nsown;                 // # of surf elements owned by this proc
  int nchoose;               // # of surf elements output by this proc
  int *cglobal;              // indices of global elements for nchoose

  typedef void (ComputePropertySurf::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  double *buf;

  void pack_id(int);

  void pack_v1x(int);
  void pack_v1y(int);
  void pack_v1z(int);
  void pack_v2x(int);
  void pack_v2y(int);
  void pack_v2z(int);
  void pack_v3x(int);
  void pack_v3y(int);
  void pack_v3z(int);
  void pack_xc(int);
  void pack_yc(int);
  void pack_zc(int);

  void pack_area(int);
  void pack_normx(int);
  void pack_normy(int);
  void pack_normz(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Invalid compute property/grid field for 2d simulation

Fields that reference z-dimension properties cannot be used
in a 2d simulation.

E: Invalid keyword in compute property/grid command

Self-explanatory.

*/
