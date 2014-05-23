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
  void check_neighbor_norm_2d();
  void check_neighbor_norm_3d();
  void check_point_near_surf_2d();
  void check_point_near_surf_3d();

  void point_line_compare(int, Surf::Line *, double, int &, int &);
  void point_tri_compare(int, Surf::Tri *, double, int &, int &, int, int, int);

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

/* ERROR/WARNING messages:

E: Cannot read_surf before grid is defined

UNDOCUMENTED

E: Cannot read_surf before grid ghost cells are defined

UNDOCUMENTED

E: Cannot read_surf after particles are defined

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Invalid reuse of surface ID in read_surf command

UNDOCUMENTED

E: Surf file cannot parse Points section

UNDOCUMENTED

E: Surf file cannot parse Lines section

UNDOCUMENTED

E: Surf file cannot parse Triangles section

UNDOCUMENTED

E: Invalid read_surf command

UNDOCUMENTED

E: Invalid read_surf geometry transformation for 2d simulation

UNDOCUMENTED

E: Unexpected end of data file

SPARTA hit the end of the data file while attempting to read a
section.  Something is wrong with the format of the data file.

E: Surf file cannot contain lines for 3d simulation

UNDOCUMENTED

E: Surf file cannot contain triangles for 2d simulation

UNDOCUMENTED

E: Surf files does not contain points

UNDOCUMENTED

E: Surf files does not contain lines

UNDOCUMENTED

E: Surf files does not contain triangles

UNDOCUMENTED

E: Unexpected end of surf file

UNDOCUMENTED

E: Incorrect point format in surf file

UNDOCUMENTED

E: Incorrect line format in surf file

UNDOCUMENTED

E: Invalid point index in line

UNDOCUMENTED

E: Incorrect triangle format in surf file

UNDOCUMENTED

E: Invalid point index in triangle

UNDOCUMENTED

E: %d read_surf points are not inside simulation box

UNDOCUMENTED

E: Surface check failed with %d duplicate points

UNDOCUMENTED

E: Surface check failed with %d unmatched points

UNDOCUMENTED

E: Surface check failed with %d duplicate edges

UNDOCUMENTED

E: Surface check failed with %d unmatched edges

UNDOCUMENTED

E: Surface check failed with %d infinitely thin line pairs

UNDOCUMENTED

W: Surface check found %d nearly infinitely thin line pairs

UNDOCUMENTED

E: Surface check failed with %d infinitely thin triangle pairs

UNDOCUMENTED

W: Surface check found %d nearly infinitely thin triangle pairs

UNDOCUMENTED

E: Surface check failed with %d points on lines

UNDOCUMENTED

W: Surface check found %d points nearly on lines

UNDOCUMENTED

E: Surface check failed with %d points on triangles

UNDOCUMENTED

W: Surface check found %d points nearly on triangles

UNDOCUMENTED

E: Cannot open gzipped file

SPARTA was compiled without support for reading and writing gzipped
files through a pipeline to the gzip program with -DSPARTA_GZIP.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: %d read_surf point pairs are too close

UNDOCUMENTED

*/
