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

DumpStyle(surf,DumpSurf)

#else

#ifndef SPARTA_DUMP_SURF_H
#define SPARTA_DUMP_SURF_H

#include "dump.h"

namespace SPARTA_NS {

class DumpSurf : public Dump {
 public:
  DumpSurf(class SPARTA *, int, char **);
  ~DumpSurf();

 private:
  int nevery;                // dump frequency to check Fix against
  int groupbit;              // mask for surface group

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

  int ncustom;               // # of surf Custom attributes used by dump
  char **id_custom;          // their IDs
  int *custom;               // list of indices for the Custom attributes

  int dimension;
  int nown;                  // # of surf elements owned by this proc
  int nchoose;               // # of surf elements output by this proc
  int *cglobal;              // indices of global elements for nchoose
  int *clocal;               // indices of local owned elements for nchoose
  double *buflocal;          // buffer for per-surf element values

  int distributed,implicit;  // Surf settings

  int firstflag;

  // private methods

  void init_style();
  void write_header(bigint);
  int count();
  void pack();
  void write_data(int, double *);

  int parse_fields(int, char **);
  int add_custom(char *);
  int add_compute(char *);
  int add_fix(char *);
  int add_variable(char *);

  typedef void (DumpSurf::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(bigint);
  void header_item(bigint);

  typedef void (DumpSurf::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_string(int, double *);
  void write_text(int, double *);

  // customize by adding a method prototype

  typedef void (DumpSurf::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_compute(int);
  void pack_fix(int);
  void pack_variable(int);
  void pack_custom(int);

  void pack_id(int);
  void pack_type(int);
  void pack_v1x(int);
  void pack_v1y(int);
  void pack_v1z(int);
  void pack_v2x(int);
  void pack_v2y(int);
  void pack_v2z(int);
  void pack_v3x(int);
  void pack_v3y(int);
  void pack_v3z(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: No dump surf attributes specified

Self-explanatory.

E: Invalid attribute in dump surf command

Self-explanatory.

E: Could not find dump surf compute ID

Self-explanatory.

E: Could not find dump surf fix ID

Self-explanatory.

E: Dump surf and fix not computed at compatible times

Fixes generate values on specific timesteps.  The dump surf output
does not match these timesteps.

E: Could not find dump surf variable name

Self-explanatory.

E: Invalid dump surf field for 2d simulation

Self-explanatory.

E: Dump surf compute does not compute per-surf info

Self-explanatory.

E: Dump surf compute does not calculate per-surf array

Self-explanatory.

E: Dump surf compute vector is accessed out-of-range

Self-explanatory.

E: Dump surf fix does not compute per-surf info

Self-explanatory.

E: Dump surf fix does not compute per-surf array

Self-explanatory.

E: Dump surf fix vector is accessed out-of-range

Self-explanatory.

E: Dump surf variable is not surf-style variable

Self-explanatory.

*/
