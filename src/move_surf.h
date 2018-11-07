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

#ifdef COMMAND_CLASS

CommandStyle(move_surf,MoveSurf)

#else

#ifndef SPARTA_MOVE_SURF_H
#define SPARTA_MOVE_SURF_H

#include "stdio.h"
#include "pointers.h"
#include "surf.h"
#include "hash.h"

#ifdef SPARTA_MAP
#include <map>
#elif SPARTA_UNORDERED_MAP
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

namespace SPARTA_NS {

class MoveSurf : protected Pointers {
 public:
  int mode;                 // 0 = move_surf command, 1 = fix move/surf command
  int groupbit;             // FixMoveSurf sets surf group

#ifdef SPARTA_MAP
  struct Point3d {
    double pt[3];
    
    bool operator <(const Point3d& other) const {
      if (pt[0] < other.pt[0]) return 1;
      else if (pt[0] > other.pt[0]) return 0;
      if (pt[1] < other.pt[1]) return 1;
      else if (pt[1] > other.pt[1]) return 0;
      if (pt[2] < other.pt[2]) return 1;
      else if (pt[2] > other.pt[2]) return 0;
      return 0;
    }
  };

#else
  struct Point3d {
    double pt[3];

    bool operator ==(const Point3d &other) const { 
      if (pt[0] != other.pt[0]) return 0;
      if (pt[1] != other.pt[1]) return 0;
      if (pt[2] != other.pt[2]) return 0;
      return 1;
    }
  };

  struct Hasher3d {
    uint32_t operator ()(const Point3d& one) const {
      return hashlittle(one.pt,3*sizeof(double),0);
    }
  };
#endif

  // hash for surface element points

#ifdef SPARTA_MAP
  std::map<Point3d,int> *hash;
#elif SPARTA_UNORDERED_MAP
  std::unordered_map<Point3d,int,Hasher3d> *hash;
#else
  std::tr1::unordered_map<Point3d,int,Hasher3d> *hash;
#endif

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
