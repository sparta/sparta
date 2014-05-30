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

CommandStyle(read_restart,ReadRestart)

#else

#ifndef SPARTA_READ_RESTART_H
#define SPARTA_READ_RESTART_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class ReadRestart : protected Pointers {
 public:
  ReadRestart(class SPARTA *);
  void command(int, char **);

 private:
  int me,nprocs,nprocs_file,multiproc_file;
  FILE *fp;
  int nfix_restart_global,nfix_restart_peratom;

  int multiproc;             // 0 = proc 0 writes for all
                             // else # of procs writing files

  bigint nparticle_file;
  int nunsplit_file,nsplit_file,nsub_file;
  int npoint_file,nsurf_file;

  void file_search(char *, char *);
  void header(int);
  void box_params();
  void particle_params();
  void grid_params();
  int surf_params();
  void file_layout();

  void magic_string();
  void endian();
  int version_numeric();

  int read_int();
  bigint read_bigint();
  double read_double();
  char *read_string();
  void read_int_vec(int, int *);
  void read_double_vec(int, double *);
  void read_char_vec(int, char *);
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

E: Cannot open restart file %s

Self-explanatory.

E: Invalid flag in peratom section of restart file

The format of this section of the file is not correct.

E: Did not assign all restart atoms correctly

UNDOCUMENTED

E: Cannot open dir to search for restart file

Using a "*" in the name of the restart file will open the current
directory to search for matching file names.

E: Found no restart file matching pattern

When using a "*" in the restart file name, no matching file was found.

E: Restart file incompatible with current version

This is probably because you are trying to read a file created with a
version of SPARTA that is too old compared to the current version.
Use your older version of SPARTA and convert the restart file
to a data file.

E: Smallint setting in spatype.h is not compatible

UNDOCUMENTED

E: Cellint setting in spatype.h is not compatible

UNDOCUMENTED

E: Bigint setting in spatype.h is not compatible

UNDOCUMENTED

W: Restart file used different # of processors

The restart file was written out by a SPARTA simulation running on a
different number of processors.  Due to round-off, the trajectories of
your restarted simulation may diverge a little more quickly than if
you ran on the same # of processors.

E: Invalid flag in header section of restart file

Unrecognized entry in restart file.

E: Invalid flag in particle section of restart file

UNDOCUMENTED

E: Invalid flag in grid section of restart file

UNDOCUMENTED

E: Invalid flag in surf section of restart file

UNDOCUMENTED

E: Restart file is not a multi-proc file

The file is inconsistent with the filename you specified for it.

E: Restart file is a multi-proc file

The file is inconsistent with the filename you specified for it.

E: Invalid flag in layout section of restart file

UNDOCUMENTED

E: Invalid SPARTA restart file

UNDOCUMENTED

E: Restart file byte ordering is swapped

The file was written on a machine with different byte-ordering than
the machine you are reading it on.  Convert it to a text data file
instead, on the machine you wrote it on.

E: Restart file byte ordering is not recognized

The file does not appear to be a SPARTA restart file since it doesn't
contain a recognized byte-orderomg flag at the beginning.

*/
