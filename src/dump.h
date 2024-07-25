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

#ifndef SPARTA_DUMP_H
#define SPARTA_DUMP_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class Dump : protected Pointers {
 public:
  char *id;                  // user-defined name of Dump
  char *style;               // style of Dump

  int first_flag;            // 0 if no initial dump, 1 if yes initial dump
  int clearstep;             // 1 if dump invokes computes, 0 if not

  int comm_forward;          // size of forward communication (0 if none)
  int comm_reverse;          // size of reverse communication (0 if none)

  Dump(class SPARTA *, int, char **);
  virtual ~Dump();
  void init();
  virtual void write();
  virtual void reset_grid_count() {}
  void modify_params(int, char **);
  virtual bigint memory_usage();

 protected:
  int me,nprocs;             // proc info

  char *filename;            // user-specified file
  int compressed;            // 1 if dump file is written compressed, 0 no
  int binary;                // 1 if dump file is written binary, 0 no
  int multifile;             // 0 = one big file, 1 = one file per timestep
  int multiproc;             // 0 = proc 0 writes for all, 1 = one file/proc
                             // else # of procs writing files
  int nclusterprocs;         // # of procs in my cluster that write to one file
  int filewriter;            // 1 if this proc writes a file, else 0
  int fileproc;              // ID of proc in my cluster who writes to file
  char *multiname;           // filename with % converted to cluster ID
  MPI_Comm clustercomm;      // MPI communicator within my cluster of procs

  int header_flag;           // 0 = item, 2 = xyz
  int flush_flag;            // 0 if no flush, 1 if flush every dump
  int append_flag;           // 1 if open file in append mode, 0 if not
  int buffer_allow;          // 1 if style allows for buffer_flag, 0 if not
  int buffer_flag;           // 1 if buffer output as one big string, 0 if not
  int padflag;               // timestep padding in filename
  int singlefile_opened;     // 1 = one big file, already opened, else 0

  char boundstr[9];          // encoding of boundary flags

  char *format;              // format string for the file write
  char *format_default;      // default format string

  char *format_line_user;    // user-specified format strings
  char *format_float_user;
  char *format_int_user;
  char *format_bigint_user;
  char **format_column_user;

  FILE *fp;                  // file to write dump to
  int size_one;              // # of quantities for one entity
  int nme;                   // # of entities in this dump from me
  int nsme;                  // # of chars in string output from me

  double boxxlo,boxxhi;      // local copies of domain values
  double boxylo,boxyhi;
  double boxzlo,boxzhi;

  bigint ntotal;             // # of per-atom lines in snapshot

  int maxbuf;                // size of buf
  double *buf;               // memory for dumped quantities
  int maxsbuf;               // size of sbuf
  char *sbuf;                // memory for atom quantities in string format

  int *vtype;                // type of each field (INT, DOUBLE, etc)
  char **vformat;            // format string for each field

  int convert_string(int, double *);

  virtual void init_style() = 0;
  virtual void openfile();
  virtual int modify_param(int, char **) {return 0;}
  virtual void write_header(bigint) = 0;
  virtual int count() = 0;
  virtual void pack() = 0;
  virtual void write_data(int, double *) = 0;
};

}

#endif

/* ERROR/WARNING messages:

E: Too much per-proc info for dump

Number of local atoms times number of columns must fit in a 32-bit
integer for dump.

E: Too much buffered per-proc info for dump

Number of dumped values per processor cannot exceed a small integer
(~2 billion values).

E: Cannot open gzipped file

SPARTA was compiled without support for reading and writing gzipped
files through a pipeline to the gzip program with -DSPARTA_GZIP.

E: Cannot open dump file

The output file for the dump command cannot be opened.  Check that the
path and name are correct.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Dump_modify buffer yes not allowed for this style

Not all dump styles allow dump_modify buffer yes.  See the dump_modify
doc page.

E: Cannot use dump_modify fileper without % in dump file name

Self-explanatory.

E: Cannot use dump_modify nfile without % in dump file name

Self-explanatory.

*/
