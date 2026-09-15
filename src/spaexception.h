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

#ifndef SPARTA_SPAEXCEPTION_H
#define SPARTA_SPAEXCEPTION_H

#include "mpi.h"

#include <exception>
#include <string>

namespace SPARTA_NS {

// thrown by Error::all() and Error::universe_all(),
// i.e. in situations where all MPI ranks fail collectively.
// callers embedding SPARTA as a library can catch it and recover,
// since every rank unwinds to the library interface

class SpartaException : public std::exception {
 public:
  explicit SpartaException(const char *msg) : message(msg) {}
  explicit SpartaException(const std::string &msg) : message(msg) {}
  const char *what() const noexcept override { return message.c_str(); }

 protected:
  std::string message;
};

// thrown by Error::one() and Error::universe_one(),
// i.e. when only one (or a subset of) MPI rank(s) hit the error.
// in parallel the catch site must call MPI_Abort() on the stored
// communicator since the other ranks cannot be unwound;
// on a single rank it can be handled like SpartaException

class SpartaAbortException : public SpartaException {
 public:
  SpartaAbortException(const char *msg, MPI_Comm comm) : SpartaException(msg), universe(comm) {}
  SpartaAbortException(const std::string &msg, MPI_Comm comm) :
      SpartaException(msg), universe(comm) {}
  MPI_Comm get_universe() const { return universe; }

 protected:
  MPI_Comm universe;
};

}    // namespace SPARTA_NS

#endif
