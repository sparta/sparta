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

#ifndef SPARTA_MARCHING_CUBES_H
#define SPARTA_MARCHING_CUBES_H

#include "pointers.h"
#include "surf.h"

namespace SPARTA_NS {

class MarchingCubes : protected Pointers {
 public:
  MarchingCubes(class SPARTA *, int, double);
  ~MarchingCubes() {}
  void invoke(double **, int *, int **);
  void cleanup();

 private:
  int me,ggroup;
  double thresh;

  double *lo,*hi;
  int v000,v001,v010,v011,v100,v101,v110,v111;
  double v000iso,v001iso,v010iso,v011iso,v100iso,v101iso,v110iso,v111iso;
  int bit0,bit1,bit2,bit3,bit4,bit5,bit6,bit7;
  double pt[36][3];

  int config;     // configuration of the active cube
  int subconfig;  // subconfiguration of the active cube

  // message datums for cleanup()

  struct SendDatum {
    int sendcell,sendface;
    int othercell,otherface;
    int inwardnorm;            // for sending cell
    Surf::Tri tri1,tri2;
  };

  double interpolate(double, double, double, double);
  int add_triangle(int *, int);
  bool test_face(int);
  bool test_interior(int, int);
  bool modified_test_interior(int, int);
  int interior_ambiguity(int, int);
  int interior_ambiguity_verification(int);
  bool interior_test_case13();
  void print_cube();
};

}

#endif

/* ERROR/WARNING messages:

*/
