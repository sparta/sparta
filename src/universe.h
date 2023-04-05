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

#ifndef SPARTA_UNIVERSE_H
#define SPARTA_UNIVERSE_H

#include "stdio.h"
#include "pointers.h"

namespace SPARTA_NS {

class Universe : protected Pointers {
 public:
  char *version;          // SPARTA version string = date

  MPI_Comm uworld;        // communicator for entire universe
  int me,nprocs;          // my place in universe

  FILE *uscreen;          // universe screen output
  FILE *ulogfile;         // universe logfile

  int existflag;          // 1 if universe exists due to -partition flag
  int nworlds;            // # of worlds in universe
  int iworld;             // which world I am in
  int *procs_per_world;   // # of procs in each world
  int *root_proc;         // root proc in each world

  MPI_Comm uorig;         // original communicator passed to SPARTA instance
  int *uni2orig;          // proc I in universe uworld is
                          // proc uni2orig[I] in original communicator

  Universe(class SPARTA *, MPI_Comm);
  ~Universe();
  void add_world(char *);
  int consistent();
};

}

#endif
