/* ----------------------------------------------------------------------
   DSMC - Sandia parallel DSMC code
   www.sandia.gov/~sjplimp/dsmc.html
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2011) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level DSMC directory.
------------------------------------------------------------------------- */

#ifndef DSMC_DUMP_H
#define DSMC_DUMP_H

#include "stdio.h"
#include "pointers.h"

namespace DSMC_NS {

class Dump : protected Pointers {
 public:
  char *id;                  // user-defined name of Dump
  char *style;               // style of Dump

  int first_flag;            // 0 if no initial dump, 1 if yes initial dump

  int comm_forward;          // size of forward communication (0 if none)
  int comm_reverse;          // size of reverse communication (0 if none)

  Dump(class DSMC *, int, char **);
  virtual ~Dump();
  void init();
  virtual void write();
  void modify_params(int, char **);
  virtual bigint memory_usage();

 protected:
  int me,nprocs;             // proc info

  char *filename;            // user-specified file
  int compressed;            // 1 if dump file is written compressed, 0 no
  int binary;                // 1 if dump file is written binary, 0 no
  int multifile;             // 0 = one big file, 1 = one file per timestep
  int multiproc;             // 0 = proc 0 writes for all, 1 = one file/proc

  int header_flag;           // 0 = item, 2 = xyz
  int flush_flag;            // 0 if no flush, 1 if flush every dump
  int append_flag;           // 1 if open file in append mode, 0 if not
  int padflag;               // timestep padding in filename
  int singlefile_opened;     // 1 = one big file, already opened, else 0

  char boundstr[9];          // encoding of boundary flags
  char *format_default;      // default format string
  char *format_user;         // format string set by user
  char *format;              // format string for the file write
  FILE *fp;                  // file to write dump to
  int size_one;              // # of quantities for one particle
  int nme;                   // # of particles in this dump from me

  double boxxlo,boxxhi;      // local copies of domain values
  double boxylo,boxyhi;
  double boxzlo,boxzhi;

  bigint ntotal;             // # of per-atom lines in snapshot

  int maxbuf;                // size of buf
  double *buf;               // memory for atom quantities

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
