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

#ifndef SPARTA_TIMER_H
#define SPARTA_TIMER_H

#include "pointers.h"

namespace SPARTA_NS {

// TIME_SYNC = time waiting in blocking global collectives inside the run
//   loop: the per-step MPI_Allreduce in move() and the closing barrier.
//   On a fast rank this is idle time, so it measures load imbalance.
//   It was previously discarded by an argument-less Timer::stamp() and so
//   silently landed in the "Other" residual printed by Finish::end(), which
//   is what made it invisible.
// Per-step bookkeeping and post-run work are still left to that residual;
//   measuring them as a section of their own showed them to be ~0, so they
//   do not need one.

enum{TIME_LOOP,TIME_MOVE,TIME_COLLIDE,TIME_SORT,TIME_COMM,TIME_MODIFY,
     TIME_OUTPUT,TIME_SYNC,TIME_N};

class Timer : protected Pointers {
 public:
  double *array;

  Timer(class SPARTA *);
  ~Timer();
  void init();
  void stamp();
  void stamp(int);
  void barrier_start(int);
  void barrier_stop(int);
  double elapsed(int);

  // initialize timeout timer
  // must be called after reset_timeout(), so that the setting saved for a
  //   later restore is the caller's real limit and never an expired 0.0
  void init_timeout();

  // trigger enforced timeout
  // persist = 1 also expires the saved setting, so the timeout survives the
  //   reset_timeout() that Run::command() performs at the start of every run.
  //   "fix halt ... error soft" needs this: it is documented to skip
  //   subsequent run commands, whereas an interrupt from the library
  //   interface must not outlive the run it interrupted
  void force_timeout(int persist = 0)
  {
    _timeout = 0.0;
    if (persist) _s_timeout = 0.0;
  }

  // restore original timeout setting after enforce timeout
  void reset_timeout() { _timeout = _s_timeout; }

  // get remaining time in seconds. 0.0 if inactive, negative if expired
  double get_timeout_remain();

  // CPU utilization in percent since the previous call
  // used by the library interface, e.g. for a GUI progress display
  double cpu_usage();

  // print timeout message
  void print_timeout(FILE *);

  // check for timeout. inline wrapper around internal
  // function to reduce overhead in case there is no check.
  bool check_timeout(int step)
  {
    if (_timeout == 0.0) return true;
    if (_nextcheck != step)
      return false;
    else
      return _check_timeout();
  }

 private:
  double previous_time;
  double timeout_start;  // wall time the timeout window opened, set by
                         //   init_timeout().  only meaningful when
                         //   _timeout >= 0.0
  double last_cpu_secs;     // process CPU time at last cpu_usage() call
  double last_cpu_wall;     // wall time at last cpu_usage() call
  double _timeout;      // max allowed wall time in seconds. infinity if negative
  double _s_timeout;    // copy of timeout for restoring after a forced timeout
  int _checkfreq;       // frequency of timeout checking
  int _nextcheck;       // loop number of next timeout check

  // check for timeout
  bool _check_timeout();
};

}

#endif
