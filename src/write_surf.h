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

CommandStyle(write_surf,WriteSurf)

#else

#ifndef SPARTA_WRITE_SURF_H
#define SPARTA_WRITE_SURF_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class WriteSurf : protected Pointers {
 public:
  int statflag;

  WriteSurf(class SPARTA *);
  void command(int, char **);

 private:
  int me,nprocs;
  int dim;
  FILE *fp;

  int pointflag;             // 1/0 to include/exclude Points section in file
  int typeflag;              // 1/0 to include/exclude element types

  int ncustom;               // number of custom per-surf attributes to output
  int *index_custom,*type_custom,*size_custom;  // flags for custom attributes
  int nvalues_custom;        // number of custom values per surf
  
  int multiproc;             // 0 = proc 0 writes for all
                             // else # of procs writing files
  int filewriter;            // 1 if this proc writes to file, else 0
  int icluster;              // which cluster I am in
  int nclusterprocs;         // # of procs in my cluster that write to one file
  int fileproc;              // ID of proc in my cluster who writes to file

  struct SurfIDType {
    surfint id;
    int type;
  };

  void write_file(char *);
  void write_file_all_points(char *);
  void write_file_all_nopoints(char *);
  void write_file_distributed_points(char *);
  void write_file_distributed_nopoints(char *);

  void write_base(char *);
  void open(char *);
  void write_header(int);

  void pack_custom(int, double **);
  void write_custom_all(int);
  void write_custom_distributed(int, double **);

  // union data struct for packing 32-bit and 64-bit ints into double bufs
  // this avoids aliasing issues by having 2 pointers (double,int)
  //   to same buf memory
  // constructor for 32-bit int prevents compiler
  //   from possibly calling the double constructor when passed an int
  // copy to a double *buf:
  //   buf[m++] = ubuf(foo).d, where foo is a 32-bit or 64-bit int
  // copy from a double *buf:
  //   foo = (int) ubuf(buf[m++]).i;, where (int) or (tagint) match foo
  //   the cast prevents compiler warnings about possible truncation

  union ubuf {
    double d;
    int64_t i;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };
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
