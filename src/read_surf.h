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

#include "stdio.h"
#include "pointers.h"
#include "surf.h"

namespace SPARTA_NS {

class ReadSurf : protected Pointers {
 public:
  ReadSurf(class SPARTA *);
  virtual ~ReadSurf();
  virtual void command(int, char **);

 protected:
  int me;
  char *line,*keyword,*buffer;
  FILE *fp;
  int compressed;

  int dim;
  double origin[3];

  struct Point {
    double x[3];
  };

  struct Line {
    int type,mask;          // type and mask of the line
    int p1,p2;              // indices of points in line segment
  };

  struct Tri {
    int type,mask;          // type and mask of the triangle
    int p1,p2,p3;           // indices of points in triangle
  };

  Point *pts;
  Line *lines;
  Tri *tris;
  int npoint,nline,ntri;
  int maxpoint,maxline,maxtri;
  int nline_old,nline_new;
  int ntri_old,ntri_new;

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

  void push_points_to_boundary(double);
  void check_neighbor_norm_2d();
  void check_neighbor_norm_3d();
  void check_point_near_surf_2d();
  void check_point_near_surf_3d();

  void point_line_compare(double *, double *, double *, double, int &, int &);
  void point_tri_compare(double *, double *, double *, double *, double *,
                         double, int &, int &, int, int, int);

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

Self-explanatory.

E: Cannot read_surf before grid ghost cells are defined

This needs to be documented if keep this restriction.

E: Cannot read_surf after particles are defined

This is because the newly read surface objects may enclose particles.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Invalid reuse of surface ID in read_surf command

Surface IDs must be unique.

E: Read_surf did not find points section of surf file

Expected Parents section but did not find keyword.

E: Read_surf did not find lines section of surf file

Expected Lines section but did not find keyword.

E: Read_surf did not find triangles section of surf file

Expected Triangles section but did not find keyword.

E: Invalid read_surf command

Self-explanatory.

E: Invalid read_surf geometry transformation for 2d simulation

Cannot perform a transformation that changes z cooridinates of points
for a 2d simulation.

E: Unexpected end of data file

SPARTA hit the end of the data file while attempting to read a
section.  Something is wrong with the format of the data file.

E: Surf file cannot contain lines for 3d simulation

Self-explanatory.

E: Surf file cannot contain triangles for 2d simulation

Self-explanatory.

E: Surf file does not contain points

Self-explanatory.

E: Surf file does not contain lines

Required for a 2d simulation.

E: Surf file does not contain triangles

Required for a 3d simulation.

E: Unexpected end of surf file

Self-explanatory.

E: Incorrect point format in surf file

Self-explanatory.

E: Incorrect line format in surf file

Self-explanatory.

E: Invalid point index in line

Self-explanatory.

E: Incorrect triangle format in surf file

Self-explanatory.

E: Invalid point index in triangle

Self-explanatory.

E: %d read_surf points are not inside simulation box

If clipping was not performed, all points in surf file
must be inside (or on surface of) simulation box.

E: Surface check failed with %d duplicate points

One or more points appeared in more than 2 lines.

E: Surface check failed with %d unmatched points

One or more points did not appear in a line, or appeared only once and
point is not on surface of simulation box.

E: Surface check failed with %d duplicate edges

One or more edges appeared in more than 2 triangles.

E: Surface check failed with %d unmatched edges

One or more edges did not appear in a triangle, or appeared only once
and edge is not on surface of simulation box.

E: Surface check failed with %d infinitely thin line pairs

Two adjacent lines have normals in opposite directions
indicating the lines overlay each other.

W: Surface check found %d nearly infinitely thin line pairs

Two adjacent lines have normals in nearly opposite directions
indicating the lines nearly overlay each other.

E: Surface check failed with %d infinitely thin triangle pairs

Two adjacent triangles have normals in opposite directions indicating
the triangles overlay each other.

W: Surface check found %d nearly infinitely thin triangle pairs

Two adjacent triangles have normals in nearly opposite directions
indicating the triangles nearly overlay each other.

E: Surface check failed with %d points on lines

One or more points are on a line they are not an end point of, which
indicates an ill-formed surface.

W: Surface check found %d points nearly on lines

One or more points are nearly on a line they are not an end point of,
which indicates an ill-formed surface.

E: Surface check failed with %d points on triangles

One or more points are on a triangle they are not an end point of,
which indicates an ill-formed surface.

W: Surface check found %d points nearly on triangles

One or more points are nearly on a triangle they are not an end point
of, which indicates an ill-formed surface.

E: Cannot open gzipped file

SPARTA was compiled without support for reading and writing gzipped
files through a pipeline to the gzip program with -DSPARTA_GZIP.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

E: %d read_surf point pairs are too close

A pair of points is very close together, relative to grid size,
inidicating the grid is too large, or an ill-formed surface.

*/
