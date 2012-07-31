/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(read_surf,ReadSurf)

#else

#ifndef SPARTA_READ_SURF_H
#define SPARTA_READ_SURF_H

#include "pointers.h"
#include "surf.h"

namespace SPARTA_NS {

class ReadSurf : protected Pointers {
 public:
  ReadSurf(class SPARTA *);
  ~ReadSurf();
  void command(int, char **);

 private:
  int me;
  char *line,*keyword,*buffer;
  FILE *fp;
  int compressed;

  int dimension;
  int isc;
  double origin[3];

  Surf::Point *pts;
  Surf::Line *lines;
  Surf::Tri *tris;
  int npoint_old,nline_old,ntri_old;
  int npoint_new,nline_new,ntri_new;
  int maxpoint,maxline,maxtri;

  int **edge;
  int nedge,maxedge;

  void header();
  void read_points();
  void read_lines();
  void read_tris();

  void translate(double, double, double);
  void scale(double, double, double);
  void rotate(double, double, double, double);
  void invert();
  void clip2d();
  void clip3d();

  void check_point_inside();
  void check_watertight_2d();
  void check_watertight_3d();

  int find_edge(int, int);
  void add_edge(int, int, int);

  double shortest_line();
  void smallest_tri(double &, double &);

  void open(char *);
  void parse_keyword(int);
  int count_words(char *);
};

}

#endif
#endif
