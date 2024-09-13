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

CommandStyle(read_restart,ReadRestart)

#else

#ifndef SPARTA_READ_RESTART_H
#define SPARTA_READ_RESTART_H

#include "stdio.h"
#include "pointers.h"
#include "surf.h"

namespace SPARTA_NS {

class ReadRestart : protected Pointers {
 public:
  ReadRestart(class SPARTA *);
  void command(int, char **);

 private:
  int me,nprocs,nprocs_file,multiproc_file;
  int mem_limit_flag,mem_limit_file;
  int nfix_restart_global,nfix_restart_peratom;
  FILE *fp;

  int multiproc;             // 0 = proc 0 writes for all
                             // else # of procs writing files
  int filereader;            // 1 if this proc reads file, else 0
  int procmatch_check;
  int procmatch;

  bigint nparticle_file;
  bigint nunsplit_file;
  int nsplit_file,nsub_file;
  int npoint_file,nsurf_file;

  // locally stored surfs

  int nsurf,maxsurf;
  Surf::Line *lines;
  Surf::Tri *tris;
  int ncustom_surf;
  int nvalues_custom_surf;
  double **cvalues;

  // local methods

  void file_search(char *, char *);
  void header(int);
  void box_params();
  void particle_params();
  void grid_params();
  int surf_params();

  void read_grid_particles(char *);
  void read_gp_single_file_same_procs();
  void read_gp_single_file_diff_procs();
  void read_gp_multi_file_less_procs(char *);
  void read_gp_multi_file_more_procs(char *);
  void read_gp_multi_file_less_procs_memlimit(char *);
  void read_gp_multi_file_more_procs_memlimit(char *);

  void read_surfs(char *);
  void read_surfs_single_file();
  void read_surfs_multi_file_less_procs(char *);
  void read_surfs_multi_file_more_procs(char *);
  void unpack_surfs(int, char *);
  void add_line(surfint, int, int, int, double *, double *);
  void add_tri(surfint, int, int, int, double *, double *, double *);
  void add_custom(surfint, double *);

  void create_child_cells(int);
  void assign_particles(int);

  void magic_string();
  void endian();
  int version_numeric();

  int read_int();
  bigint read_bigint();
  double read_double();
  char *read_string();
  void read_int_vec(int, int *);
  void read_double_vec(int, double *);
  void read_char_vec(bigint, char *);

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

E: Cannot read_restart after simulation box is defined

The read_restart command cannot be used after a read_data,
read_restart, or create_box command.

E: Cannot (yet) use global mem/limit without % in restart file name

This feature is not yet implemented.

E: Cannot open restart file %s

The specified file cannot be opened.  Check that the path and name are
correct.  If the file is a compressed file, also check that the gzip
executable can be found and run.

E: Invalid flag in peratom section of restart file

The format of this section of the file is not correct.

E: Did not assign all restart unsplit grid cells correctly

One or more unsplit grid cells in the restart file were not assigned
to a processor.  Please report the issue to the SPARTA developers.

E: Did not assign all restart split grid cells correctly

One or more split grid cells in the restart file were not assigned
to a processor.  Please report the issue to the SPARTA developers.

E: Did not assign all restart sub grid cells correctly

One or more sub grid cells in the restart file were not assigned to a
processor.  Please report the issue to the SPARTA developers.

E: Did not assign all restart particles correctly

One or more particles in the restart file were not assigned to a
processor.  Please report the issue to the SPARTA developers.

E: Cannot open dir to search for restart file

Using a "*" in the name of the restart file will open the current
directory to search for matching file names.

E: Found no restart file matching pattern

When using a "*" in the restart file name, no matching file was found.

E: Restart file incompatible with current version

This is probably because you are trying to read a file created with a
version of SPARTA that is too old compared to the current version.

E: Smallint setting in spatype.h is not compatible

Smallint size stored in restart file is not consistent with SPARTA
version you are running.

E: Cellint setting in spatype.h is not compatible

Cellint size stored in restart file is not consistent with SPARTA
version you are running.

E: Bigint setting in spatype.h is not compatible

Bigint size stored in restart file is not consistent with SPARTA
version you are running.

W: Restart file used different # of processors

The restart file was written out by a SPARTA simulation running on a
different number of processors.  This means you will likely want to
re-balance the grid cells and particles across processors.  This can
be done using the balance or fix balance commands.

E: Invalid flag in header section of restart file

Unrecognized entry in restart file.

E: Invalid flag in particle section of restart file

Unrecognized entry in restart file.

E: Invalid flag in grid section of restart file

Unrecognized entry in restart file.

E: Invalid flag in surf section of restart file

Unrecognized entry in restart file.

E: Restart file is not a multi-proc file

The file is inconsistent with the filename specified for it.

E: Restart file is a multi-proc file

The file is inconsistent with the filename specified for it.

E: Invalid flag in layout section of restart file

Unrecognized entry in restart file.

E: Invalid SPARTA restart file

The file does not appear to be a SPARTA restart file since it does not
have the expected magic string at the beginning.

E: Restart file byte ordering is swapped

The file was written on a machine with different byte-ordering than
the machine you are reading it on.

E: Restart file byte ordering is not recognized

The file does not appear to be a SPARTA restart file since it doesn't
contain a recognized byte-ordering flag at the beginning.

*/
