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

CommandStyle(read_surf,ReadSurf)

#else

#ifndef SPARTA_READ_SURF_H
#define SPARTA_READ_SURF_H

#include "mpi.h"
#include "stdio.h"
#include "pointers.h"
#include "hash3.h"
#include "surf.h"

namespace SPARTA_NS {

class ReadSurf : protected Pointers {
 public:
  ReadSurf(class SPARTA *);
  virtual ~ReadSurf();
  virtual void command(int, char **);

 protected:
  int me,nprocs;
  char *line,*keyword,*buffer;
  FILE *fp;
  int compressed;
  int distributed;

  int typeflag;
  int ncustom,nvalues_custom;
  char **name_custom;
  int *type_custom,*size_custom,*index_custom;
  int nclocal,ncmax;
  double **cvalues;
  
  int grouparg,typeadd,transparent_flag;
  int partflag,filearg;

  int multiproc;            // 1 if multiple files to read from
  int nfiles;               // # of proc files along with base file
  bigint nsurf_basefile;    // surface count in base file
  int me_file,nprocs_file;  // info for cluster of procs that read a file

  int dim;
  double origin[3];

  struct Point {
    double x[3];            // point coords
  };

  Point *pts;               // storage for points read from input file

  bigint nsurf_total_old;   // # of total system surfs before read
  int nsurf_old;            // # of surfs on this proc before read
  int nsurf_new;            // # of surfs on this proc after read
  int npoint_file;          // # of points in one file
  int nsurf_file;           // # of surfs in one file

  int filereader;
  MPI_Comm filecomm;

#ifdef SPARTA_MAP
  typedef std::map<bigint,int> MyHash;
  typedef std::map<bigint,int>::iterator MyIterator;
#elif defined SPARTA_UNORDERED_MAP
  typedef std::unordered_map<bigint,int> MyHash;
  typedef std::unordered_map<bigint,int>::iterator MyIterator;
#else
  typedef std::tr1::unordered_map<bigint,int> MyHash;
  typedef std::tr1::unordered_map<bigint,int>::iterator MyIterator;
#endif

  void read_single(char *);
  void read_multiple(char *);

  void surf_counts();
  void header();
  void base(char *file);

  void read_file(char *, int);
  void read_points();
  void read_lines(int);
  void read_tris(int);
  void add_custom(double *);
  void create_custom();
  
  void process_args(int, int, char **);

  void translate(double, double, double);
  void scale(double, double, double);
  void rotate(double, double, double, double);
  void invert();
  void clip2d();
  void clip3d();
  void transparent();

  void check_bounds();
  void push_points_to_boundary(double);
  void check_neighbor_norm_2d();
  void check_neighbor_norm_3d();

  void open(char *);
  void file_search(char *, char *);
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
