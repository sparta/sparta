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

#ifndef SPARTA_MPIIO_H
#define SPARTA_MPIIO_H

// true interface to MPIIO package
// used when MPIIO package is installed

#ifdef SPARTA_MPIIO

#if defined(MPI_STUBS)
#error "The MPIIO package cannot be compiled in serial with MPI STUBS"
#endif

#include "restart_mpiio.h"

#else

// dummy interface to MPIIO package
// needed for compiling when MPIIO package is not installed

namespace SPARTA_NS {

class RestartMPIIO {
 public:
  int mpiio_exists;

  RestartMPIIO(class SPARTA *) {mpiio_exists = 0;}
  ~RestartMPIIO() {}
  void openForRead(char *) {}
  void openForWrite(char *) {}
  void write(MPI_Offset,int,char *) {}
  void read(MPI_Offset,long,char *) {}
  void close() {}
};

}

#endif
#endif
