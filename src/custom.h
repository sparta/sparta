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
  int mode;            // accessed by fix custom
  
  Custom(class SPARTA *);
  virtual ~Custom();
  void command(int, char **);

  bigint process_actions(int, char **, int);
  bigint process_actions();

 private:
  struct Action {
    int action;
    int cindex,ctype,csize,ccol;
    int vindex,vstyle;
    int groupbit;
    class Mixture *mixture;
    class Region *region;
    int numfile,filestyle;
    char *fname;
    int colcount;
    int *cindex_file,*ctype_file,*csize_file,*ccol_file;
  };

  int naction;
  Action *actions;

  int ncoarse;
  double **xyz_coarse;
  double **values_coarse;

  bigint action_set(int, int, int, int, int, int,
                    int, class Mixture *, class Region *);
  bigint set_particle(class Mixture *, class Region *,
                      int, int, int, int, double, double *);
  bigint set_grid(int, class Region *, int, int, int, int, double, double *);
  bigint set_surf(int, class Region *, int, int, int, int, double, double *);
  bigint read_file(int, int, int *, int *, int *, int *, char *);
  void read_coarse_files(char *, int, int);
  bigint coarse_tree_neighbor_assign(int, int, int *, int *, int *, int *);
  int attribute_bracket(char *);
};

// K-D tree class
  
class KDTree : protected Pointers {
 public:
  KDTree(class SPARTA *, int, int, double **);
  virtual ~KDTree();
  void create_tree(int, int, int *);
  int find_nearest(double *, int, double &);
  void find_within_cutoff(double *, int, double, int &, int *, double *);
  void stats_tree();
  void stats_search();
  void stats_neighbor();
  
 private:
  int dim;
  double **points;
  int npoints;

  int nsearch,nneigh;
  double avedist;
  int count_node,count_leaf;
  
  struct Node {
    int which;       // STUB, BRANCH, LEAF
    int iparent;     // index of parent node
    int ipoint;      // index of LEAF point
    int left,right;  // indices of 2 child nodes of this branch node
    int splitdim;    // dimension (0,1,2) of branch node splitting
    double split;    // coordinate of branch node splitting
  };

  Node *tree;
  int ntree;
  int maxtree;

  int depthwalk(int, int, int);
  int walk_to_leaf(int, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
