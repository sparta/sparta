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

DumpStyle(particle,DumpParticle)

#else

#ifndef SPARTA_DUMP_PARTICLE_H
#define SPARTA_DUMP_PARTICLE_H

#include "dump.h"

namespace SPARTA_NS {

class DumpParticle : public Dump {
 public:
  DumpParticle(class SPARTA *, int, char **);
  virtual ~DumpParticle();
  bigint memory_usage();

 protected:
  int imix;                  // index of mixture to be dumped
  int nevery;                // dump frequency to check Fix against
  int iregion;               // -1 if no region, else which region
  int nthresh;               // # of defined threshholds
  int *thresh_array;         // array to threshhhold on for each nthresh
  int *thresh_op;            // threshhold operation for each nthresh
  double *thresh_value;      // threshhold value for each nthresh

  char *columns;             // column labels

  int nchoose;               // # of selected atoms
  int maxlocal;              // size of atom selection and variable arrays
  int *choose;               // local indices of selected atoms
  double *dchoose;           // value for each atom to threshhold against
  int *clist;                // compressed list of indices of selected atoms

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

  int ncustom;               // # of particle Custom attributes used by dump
  char **id_custom;          // their IDs
  int *custom;               // list of indices for the Custom attributes

  int ntypes;                // # of particle types
  char **typenames;             // array of element names for each type

  virtual void init_style();
  void write_header(bigint);
  int count();
  void pack();
  void write_data(int, double *);

  int parse_fields(int, char **);
  int add_custom(char *);
  int add_compute(char *);
  int add_fix(char *);
  int add_variable(char *);

  virtual int modify_param(int, char **);

  typedef void (DumpParticle::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(bigint);
  void header_item(bigint);

  typedef void (DumpParticle::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_string(int, double *);
  void write_text(int, double *);

  // customize by adding a method prototype

  typedef void (DumpParticle::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_compute(int);
  void pack_fix(int);
  void pack_variable(int);
  void pack_custom(int);

  void pack_id(int);
  void pack_type(int);
  void pack_proc(int);

  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_xs(int);
  void pack_ys(int);
  void pack_zs(int);

  void pack_vx(int);
  void pack_vy(int);
  void pack_vz(int);

  void pack_ke(int);
  void pack_erot(int);
  void pack_evib(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: No dump particle attributes specified

Self-explanatory.

E: Invalid attribute in dump particle command

Self-explanatory.

E: Could not find dump particle compute ID

Self-explanatory.

E: Could not find dump particle fix ID

Self-explanatory.

E: Dump particle and fix not computed at compatible times

Fixes generate values on specific timesteps.  The dump particle output
does not match these timesteps.

E: Could not find dump particle variable name

Self-explanatory.

E: Region ID for dump custom does not exist

Self-explanatory.

E: Dump particle compute does not compute per-particle info

Self-explanatory.

E: Dump particle compute does not calculate per-particle vector

Self-explanatory.

E: Dump particle compute does not calculate per-particle array

Self-explanatory.

E: Dump particle compute vector is accessed out-of-range

Self-explanatory.

E: Dump particle fix does not compute per-particle info

Self-explanatory.

E: Dump particle fix does not compute per-particle vector

Self-explanatory.

E: Dump particle fix does not compute per-particle array

Self-explanatory.

E: Dump particle fix vector is accessed out-of-range

Self-explanatory.

E: Dump particle variable is not particle-style variable

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Dump_modify region ID does not exist

Self-explanatory.

E: Invalid attribute in dump modify command

Self-explanatory.

E: Could not find dump modify compute ID

Self-explanatory.

E: Dump modify compute ID does not compute per-particle info

Self-explanatory.

E: Dump modify compute ID does not compute per-particle vector

Self-explanatory.

E: Dump modify compute ID does not compute per-particle array

Self-explanatory.

E: Dump modify compute ID vector is not large enough

Self-explanatory.

E: Could not find dump modify fix ID

Self-explanatory.

E: Dump modify fix ID does not compute per-particle info

Self-explanatory.

E: Dump modify fix ID does not compute per-particle vector

Self-explanatory.

E: Dump modify fix ID does not compute per-particle array

Self-explanatory.

E: Dump modify fix ID vector is not large enough

Self-explanatory.

E: Could not find dump modify variable name

Self-explanatory.

E: Dump modify variable is not particle-style variable

Self-explanatory.

E: Invalid dump_modify threshhold operator

Operator keyword used for threshold specification in not recognized.

*/
