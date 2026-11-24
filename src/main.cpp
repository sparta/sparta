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

#include "mpi.h"
#include "sparta.h"
#include "input.h"
#include "library.h"
#include "spaexception.h"
#include "utils.h"

#include "stdio.h"
#include "stdlib.h"

using namespace SPARTA_NS;

// for convenience

static void finalize()
{
  sparta_kokkos_finalize();
}

/* ----------------------------------------------------------------------
   main program to drive SPARTA
------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);

  // warn if SPARTA is run with superuser or administrator privileges.
  // this is done before the SPARTA instance is created, so it is only
  // printed by the executable and not when SPARTA is used as a library,
  // and only by the first MPI rank

  int me = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  if ((me == 0) && utils::is_superuser())
    fputs("\nWARNING: SPARTA is run with superuser or administrator privileges.\n"
          "WARNING: This is STRONGLY discouraged, because typos or mistakes in an\n"
          "WARNING: input script or errors in SPARTA itself can damage the entire\n"
          "WARNING: system to the point of requiring a re-installation.\n"
          "WARNING: See the section \"Security considerations\" of the SPARTA manual.\n\n",stderr);

  try {
    SPARTA *sparta = new SPARTA(argc,argv,MPI_COMM_WORLD);
    sparta->input->file();
    delete sparta;
  } catch (SpartaAbortException &e) {
    finalize();
    MPI_Abort(e.get_universe(),1);
  } catch (SpartaException &) {
    // error message was already printed by the Error class
    finalize();
    MPI_Finalize();
    exit(1);
  }

  finalize();

  MPI_Finalize();
}
