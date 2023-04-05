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

#include "mpi.h"
#include "sparta.h"
#include "input.h"

using namespace SPARTA_NS;

/* ----------------------------------------------------------------------
   main program to drive SPARTA
------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);

  SPARTA *sparta = new SPARTA(argc,argv,MPI_COMM_WORLD);
  sparta->input->file();
  delete sparta;

  MPI_Finalize();
}
