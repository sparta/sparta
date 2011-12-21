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

#ifdef DUMP_CLASS

DumpStyle(particle,DumpParticle)

#else

#ifndef DSMC_DUMP_PARTICLE_H
#define DSMC_DUMP_PARTICLE_H

#include "dump.h"

namespace DSMC_NS {

class DumpParticle : public Dump {
 public:
  DumpParticle(class DSMC *, int, char **);
  virtual ~DumpParticle();

 protected:
  int nevery;                // dump frequency to check Fix against
  int nthresh;               // # of defined threshholds
  int *thresh_array;         // array to threshhhold on for each nthresh
  int *thresh_op;            // threshhold operation for each nthresh
  double *thresh_value;      // threshhold value for each nthresh

  int *vtype;                // type of each vector (INT, DOUBLE)
  char **vformat;            // format string for each vector element

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

  int ntypes;                // # of atom types
  char **typenames;	     // array of element names for each type

  // private methods

  virtual void init_style();
  virtual void write_header(bigint);
  int count();
  void pack();
  virtual void write_data(int, double *);
  bigint memory_usage();

  int parse_fields(int, char **);
  virtual int modify_param(int, char **);

  typedef void (DumpParticle::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(bigint);
  void header_item(bigint);

  typedef void (DumpParticle::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_text(int, double *);

  // customize by adding a method prototype

  typedef void (DumpParticle::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_id(int);
  void pack_type(int);

  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_xs(int);
  void pack_ys(int);
  void pack_zs(int);

  void pack_vx(int);
  void pack_vy(int);
  void pack_vz(int);
};

}

#endif
#endif
