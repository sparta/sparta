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

#include "mpi.h"
#include "timer.h"
#include "memory.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

Timer::Timer(SPARTA *sparta) : Pointers(sparta)
{
  memory->create(array,TIME_N,"array");
}

/* ---------------------------------------------------------------------- */

Timer::~Timer()
{
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void Timer::init()
{
  for (int i = 0; i < TIME_N; i++) array[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void Timer::stamp()
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  previous_time = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::stamp(int which)
{
  // uncomment if want synchronized timing
  // MPI_Barrier(world);
  double current_time = MPI_Wtime();
  array[which] += current_time - previous_time;
  previous_time = current_time;
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_start(int which)
{
  MPI_Barrier(world);
  array[which] = MPI_Wtime();
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_stop(int which)
{
  MPI_Barrier(world);
  double current_time = MPI_Wtime();
  array[which] = current_time - array[which];
}

/* ---------------------------------------------------------------------- */

double Timer::elapsed(int which)
{
  double current_time = MPI_Wtime();
  return (current_time - array[which]);
}
