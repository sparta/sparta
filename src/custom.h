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

#ifdef COMMAND_CLASS

CommandStyle(custom,Custom)

#else

#ifndef SPARTA_CUSTOM_H
#define SPARTA_CUSTOM_H

#include "pointers.h"

namespace SPARTA_NS {

class Custom : protected Pointers {

 public:
  Custom(class SPARTA *);
  virtual ~Custom();
  void command(int, char **);

  bigint process_actions(int, char **, int);
  bigint process_actions();

 private:
  int mode;
  
  struct Action {
    int action;
    int cindex,ctype,csize,ccol;
    int vindex,vstyle;
    int groupbit;
    class Mixture *mixture;
    class Region *region;
    char *fname;
    int colcount;
    int *cindex_file,*ctype_file,*csize_file,*ccol_file;
  };

  int naction;
  Action *actions;

  bigint action_set(int, int, int, int, int, int,
                    int, class Mixture *, class Region *);
  bigint set_particle(class Mixture *, class Region *,
                      int, int, int, int, double, double *);
  bigint set_grid(int, class Region *, int, int, int, int, double, double *);
  bigint set_surf(int, class Region *, int, int, int, int, double, double *);
  bigint read_file(int, int, int *, int *, int *, int *, char *);
  int attribute_bracket(char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
