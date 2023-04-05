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

#ifdef COMMAND_CLASS

CommandStyle(move_surf,MoveSurf)

#else

#ifndef SPARTA_MOVE_SURF_H
#define SPARTA_MOVE_SURF_H

#include "stdio.h"
#include "pointers.h"
#include "surf.h"
#include "hash3.h"

namespace SPARTA_NS {

class MoveSurf : protected Pointers {
 public:
  int mode;                 // 0 = move_surf command, 1 = fix move/surf command
  int groupbit;             // FixMoveSurf sets surf group

  // hash for surface element points

#include "hash_options.h"

#ifdef SPARTA_MAP
  typedef std::map<OnePoint3d,int> MyHash;
#elif SPARTA_UNORDERED_MAP
  typedef std::unordered_map<OnePoint3d,int,OnePoint3dHash> MyHash;
#else
  typedef std::tr1::unordered_map<OnePoint3d,int,OnePoint3dHash> MyHash;
#endif

  MyHash *hash;

  MoveSurf(class SPARTA *);
  ~MoveSurf();
  void command(int, char **);
  void process_args(int, char **);
  void move_lines(double, Surf::Line *);
  void move_tris(double, Surf::Tri *);
  bigint remove_particles();

 private:
  int me,nprocs;
  int dim,action,connectflag;
  char *file,*entry;
  double theta;
  double delta[3],rvec[3],origin[3];
  FILE *fp;

  int *pselect;                    // 1 if point is moved, else 0

  int nread;
  int *readindex;
  double **oldcoord,**newcoord;

  void readfile();
  void update_points(double);
  void translate_2d(double, Surf::Line *);
  void translate_3d(double, Surf::Tri *);
  void rotate_2d(double, Surf::Line *);
  void rotate_3d(double, Surf::Tri *);
  void connect_2d_pre();
  void connect_2d_post();
  void connect_3d_pre();
  void connect_3d_post();
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
