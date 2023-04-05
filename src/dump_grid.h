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

#ifdef DUMP_CLASS

DumpStyle(grid,DumpGrid)

#else

#ifndef SPARTA_DUMP_GRID_H
#define SPARTA_DUMP_GRID_H

#include "dump.h"

namespace SPARTA_NS {

class DumpGrid : public Dump {
 public:
  DumpGrid(class SPARTA *, int, char **);
  ~DumpGrid();
  void reset_grid_count();
  bigint memory_usage();

 private:
  int nevery;                // dump frequency to check Fix against
  int groupbit;              // mask for grid group

  char *columns;             // column labels

  int nfield;                // # of keywords listed by user
  int ioptional;             // index of start of optional args

  int *field2index;          // which compute,fix,variable calcs this field
  int *argindex;             // index into compute,fix scalar_atom,vector_atom
                             // 0 for scalar_atom, 1-N for vector_atom values
  int ncompute;              // # of Compute objects used by dump
  char **id_compute;         // their IDs
  class Compute **compute;   // list of ptrs to the Compute objects

  int nfix;                  // # of Fix objects used by dump
  char **id_fix;             // their IDs
  class Fix **fix;           // list of ptrs to the Fix objects

  int nvariable;             // # of Variables used by dump
  char **id_variable;        // their names
  int *variable;             // list of indices for the Variables
  double **vbuf;             // local storage for variable evaluation

  int *cpart;                // indices into grid->cells for cells with particles
  int ncpart;                // # of owned grid cells with particles
  int maxgrid;               // max length of per-grid variable vectors

  int dimension;

  void init_style();
  void write_header(bigint);
  int count();
  void pack();
  void write_data(int, double *);

  int parse_fields(int, char **);
  int add_compute(char *);
  int add_fix(char *);
  int add_variable(char *);
  void allocate_values(int);

  typedef void (DumpGrid::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(bigint);
  void header_item(bigint);

  typedef void (DumpGrid::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_string(int, double *);
  void write_text(int, double *);

  // customize by adding a method prototype

  typedef void (DumpGrid::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_compute(int);
  void pack_fix(int);
  void pack_variable(int);

  void pack_id(int);
  void pack_split(int);
  void pack_proc(int);

  void pack_xlo(int);
  void pack_ylo(int);
  void pack_zlo(int);
  void pack_xhi(int);
  void pack_yhi(int);
  void pack_zhi(int);
  void pack_xc(int);
  void pack_yc(int);
  void pack_zc(int);
  void pack_vol(int);
  void pack_svol(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: No dump grid attributes specified

Self-explanatory.

E: Invalid attribute in dump grid command

Self-explanatory.

E: Could not find dump grid compute ID

Self-explanatory.

E: Could not find dump grid fix ID

Self-explanatory.

E: Dump grid and fix not computed at compatible times

Fixes generate values on specific timesteps.  The dump grid output
does not match these timesteps.

E: Could not find dump grid variable name

Self-explanatory.

E: Invalid dump grid field for 2d simulation

Self-explanatory.

E: Dump grid compute does not compute per-grid info

Self-explanatory.

E: Dump grid compute does not calculate per-grid array

Self-explanatory.

E: Dump grid compute vector is accessed out-of-range

Self-explanatory.

E: Dump grid fix does not compute per-grid info

Self-explanatory.

E: Dump grid fix does not compute per-grid array

Self-explanatory.

E: Dump grid fix vector is accessed out-of-range

Self-explanatory.

E: Dump grid variable is not grid-style variable

Self-explanatory.

*/
