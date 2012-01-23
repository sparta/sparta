/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(property/cell,ComputePropertyCell)

#else

#ifndef DSMC_COMPUTE_PROPERTY_CELL_H
#define DSMC_COMPUTE_PROPERTY_CELL_H

#include "compute.h"

namespace DSMC_NS {

class ComputePropertyCell : public Compute {
 public:
  ComputePropertyCell(class DSMC *, int, char **);
  ~ComputePropertyCell();
  void init();
  void compute_per_cell();
  double memory_usage();

 private:
  int nvalues;
  int nglocal;

  void pack_standard();
  void pack_count();

  typedef void (ComputePropertyCell::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_u(int);
  void pack_v(int);
  void pack_w(int);
  void pack_usq(int);
  void pack_vsq(int);
  void pack_wsq(int);
};

}

#endif
#endif
