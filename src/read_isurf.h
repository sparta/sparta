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

CommandStyle(read_isurf,ReadISurf)

#else

#ifndef SPARTA_READ_ISURF_H
#define SPARTA_READ_ISURF_H

#include "stdio.h"
#include "pointers.h"
#include "surf.h"

#ifdef SPARTA_MAP
#include <map>
#elif defined SPARTA_UNORDERED_MAP
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

namespace SPARTA_NS {

class ReadISurf : protected Pointers {
 public:
  ReadISurf(class SPARTA *);
  virtual ~ReadISurf();
  virtual void command(int, char **);

 protected:
  int me;
  int dimension;
  int count,iggroup,sgrouparg,filearg,file_pflag;
  int nx,ny,nz;
  double thresh;
  double corner[3],xyzsize[3];
  char *typefile;

  int **cvalues;
  int *svalues;

  // extra data for 3d marching cubes

  double *lo,*hi;
  int v000,v001,v010,v011,v100,v101,v110,v111;
  double v000iso,v001iso,v010iso,v011iso,v100iso,v101iso,v110iso,v111iso;
  int bit0,bit1,bit2,bit3,bit4,bit5,bit6,bit7;
  double pt[36][3];
    
  int config;     // configuration of the active cube
  int subconfig;  // subconfiguration of the active cube
    
  // message datums for cleanup_MC()

  struct SendDatum {
    int sendcell,sendface;
    int othercell,otherface;
    int inwardnorm;            // for sending cell
    Surf::Tri tri1,tri2;
  };

  // hash for assigning grid corner points to grid cells

#ifdef SPARTA_MAP
  typedef std::map<bigint,int> MyHash;
#elif defined SPARTA_UNORDERED_MAP
  typedef std::unordered_map<bigint,int> MyHash;
#else
  typedef std::tr1::unordered_map<bigint,int> MyHash;
#endif

  MyHash *hash;

  void process_args(int, int, char **);

  void read_corners(char *);
  void read_types(char *);

  void create_hash(int, int);
  void destroy_hash();

  void assign_corners(int, bigint, uint8_t *);
  void assign_types(int, bigint, int *);

  void marching_cubes(int);
  void marching_squares(int);
  double interpolate(int, int, double, double);

  // extra functions for 3d marching cubes

  int add_triangle(int *, int);
  bool test_face(int);
  bool test_interior(int, int);
  bool modified_test_interior(int, int);
  int interior_ambiguity(int, int);
  int interior_ambiguity_verification(int);
  bool interior_test_case13();
  void cleanup_MC();
  void print_cube();
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
