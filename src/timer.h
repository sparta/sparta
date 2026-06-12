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

enum{TIME_LOOP,TIME_MOVE,TIME_COLLIDE,TIME_SORT,TIME_COMM,TIME_MODIFY,TIME_OUTPUT,TIME_N};

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
  void init_timeout();

  // trigger enforced timeout
  void force_timeout() { _timeout = 0.0; }

  // restore original timeout setting after enforce timeout
  void reset_timeout() { _timeout = _s_timeout; }

  // get remaining time in seconds. 0.0 if inactive, negative if expired
  double get_timeout_remain();

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
  double timeout_start;
  double _timeout;      // max allowed wall time in seconds. infinity if negative
  double _s_timeout;    // copy of timeout for restoring after a forced timeout
  int _checkfreq;       // frequency of timeout checking
  int _nextcheck;       // loop number of next timeout check

  // check for timeout
  bool _check_timeout();
};

}

#endif
