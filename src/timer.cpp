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

/* ----------------------------------------------------------------------
   Timeout code adapted from LAMMPS (https://www.lammps.org), October 2024
   Ported to SPARTA by: Stan Moore (SNL)
   Original Author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "timer.h"
#include "comm.h"
#include "error.h"
#include "memory.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

Timer::Timer(SPARTA *sparta) : Pointers(sparta)
{
  memory->create(array,TIME_N,"array");
  // zero the timers at construction: init() only runs at the start of a run,
  // but commands such as "balance_grid rcb time" read the array before any run
  // and would otherwise use uninitialized values (garbage cell weights, and a
  // missed "no time history" warning) -- benign in a fresh process where the
  // heap is zeroed, but not when SPARTA is embedded in a long-lived process.
  for (int i = 0; i < TIME_N; i++) array[i] = 0.0;
  _timeout = -1.0;
  _s_timeout = -1.0;
  _checkfreq = 10;
  _nextcheck = -1;
  // init_timeout() sets this for real at the start of every run, but zero it
  // here as well: it is read by print_timeout() and get_timeout_remain(),
  // which the library interface may call before any run has started
  timeout_start = 0.0;
  last_cpu_secs = -1.0;
  last_cpu_wall = -1.0;
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

  // start the stamp reference at the same instant, so the first stamp()
  // inside the timed region charges from here rather than from whatever
  // ran before it.  otherwise the pre-loop setup is silently dropped.

  previous_time = array[which];
}

/* ---------------------------------------------------------------------- */

void Timer::barrier_stop(int which)
{
  // time spent waiting in the closing barrier is this rank's accumulated
  // load imbalance, not loop work.  charge it to TIME_SYNC so that
  // TIME_LOOP minus the named sections closes to ~zero and a nonzero
  // residual means something is genuinely uninstrumented.

  double pre_barrier = MPI_Wtime();
  MPI_Barrier(world);
  double current_time = MPI_Wtime();
  array[which] = current_time - array[which];
  if (which != TIME_SYNC) array[TIME_SYNC] += current_time - pre_barrier;
}

/* ---------------------------------------------------------------------- */

double Timer::elapsed(int which)
{
  double current_time = MPI_Wtime();
  return (current_time - array[which]);
}

/* ---------------------------------------------------------------------- */

void Timer::init_timeout()
{
  // open the timeout window here.  _check_timeout() and get_timeout_remain()
  //   measure from this instant, so it must be set before any run can call
  //   them -- leaving it unset made those reads use an uninitialized member

  timeout_start = MPI_Wtime();

  // snapshot the limit so reset_timeout() can restore it after a
  //   force_timeout().  callers must have already restored _timeout (see
  //   Run::command), otherwise an expired 0.0 would be latched here and every
  //   later run would exit immediately

  _s_timeout = _timeout;
  if (_timeout < 0)
    _nextcheck = -1;
  else
    _nextcheck = _checkfreq;
}

/* ---------------------------------------------------------------------- */

void Timer::print_timeout(FILE *fp)
{
  if (!fp) return;

  // format timeout setting
  if (_timeout > 0) {
    // time since init_timeout()
    const double d = MPI_Wtime() - timeout_start;
    // remaining timeout in seconds
    int s = _timeout - d;
    // remaining 1/100ths of seconds
    const int hs = 100 * ((_timeout - d) - s);
    // breaking s down into second/minutes/hours
    const int seconds = s % 60;
    s = (s - seconds) / 60;
    const int minutes = s % 60;
    const int hours = (s - minutes) / 60;
    fprintf(fp, "  Walltime left : %d:%02d:%02d.%02d\n", hours, minutes, seconds, hs);
  }
}

/* ---------------------------------------------------------------------- */

bool Timer::_check_timeout()
{
  double walltime =  MPI_Wtime() - timeout_start;
  // broadcast time to ensure all ranks act the same.
  MPI_Bcast(&walltime, 1, MPI_DOUBLE, 0, world);

  if (walltime < _timeout) {
    _nextcheck += _checkfreq;
    return false;
  } else {
    if (comm->me == 0) error->warning(FLERR, "Wall time limit reached");
    _timeout = 0.0;
    return true;
  }
}

/* ---------------------------------------------------------------------- */
double Timer::get_timeout_remain()
{
  double remain = _timeout + timeout_start - MPI_Wtime();
  // never report a negative remaining time.
  if (remain < 0.0) remain = 0.0;
  return (_timeout < 0.0) ? 0.0 : remain;
}

/* ----------------------------------------------------------------------
   return CPU utilization in percent since the previous call
   first call returns 0.0 and initializes the reference point
------------------------------------------------------------------------- */

#if defined(_WIN32)
#include <ctime>
#else
#include <sys/resource.h>
#include <sys/time.h>
#endif

double Timer::cpu_usage()
{
  double cpu_secs;

#if defined(_WIN32)
  cpu_secs = (double) clock() / CLOCKS_PER_SEC;
#else
  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
  cpu_secs = (double) ru.ru_utime.tv_sec + 1.0e-6*ru.ru_utime.tv_usec +
    (double) ru.ru_stime.tv_sec + 1.0e-6*ru.ru_stime.tv_usec;
#endif

  double wall = MPI_Wtime();

  double percent = 0.0;
  if (last_cpu_wall >= 0.0 && wall > last_cpu_wall)
    percent = 100.0*(cpu_secs-last_cpu_secs)/(wall-last_cpu_wall);

  last_cpu_secs = cpu_secs;
  last_cpu_wall = wall;
  return percent;
}
