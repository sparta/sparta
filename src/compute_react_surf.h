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

#ifdef COMPUTE_CLASS

ComputeStyle(react/surf,ComputeReactSurf)

#else

#ifndef SPARTA_COMPUTE_REACT_SURF_H
#define SPARTA_COMPUTE_REACT_SURF_H

#include "compute.h"
#include "surf.h"
#include "hash3.h"

namespace SPARTA_NS {

class ComputeReactSurf : public Compute {
 public:
  ComputeReactSurf(class SPARTA *, int, char **);
  ~ComputeReactSurf();
  virtual void init();
  void compute_per_surf();
  virtual void clear();
  virtual void surf_tally(int, int, int, Particle::OnePart *,
                          Particle::OnePart *, Particle::OnePart *);
  virtual int tallyinfo(surfint *&);
  virtual void post_process_surf();
  bigint memory_usage();

 protected:
  int groupbit;
  int isr;                 // index of surface reaction model
  int ntotal,rpflag;
  int maxsurf,combined;

  int **reaction2col;      // 1 if ireaction triggers tally for icol

  int ntally;              // # of surfs I have tallied for
  int maxtally;            // # of tallies currently allocated
  surfint *tally2surf;     // tally2surf[I] = surf ID of Ith tally

  // hash for surf IDs

#ifdef SPARTA_MAP
  typedef std::map<surfint,int> MyHash;
#elif defined SPARTA_UNORDERED_MAP
  typedef std::unordered_map<surfint,int> MyHash;
#else
  typedef std::tr1::unordered_map<surfint,int> MyHash;
#endif

  MyHash *hash;

  int dim;                 // local copies
  Surf::Line *lines;
  Surf::Tri *tris;

  void grow_tally();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

*/
