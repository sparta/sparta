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

#ifdef COMPUTE_CLASS

ComputeStyle(surf,ComputeSurf)

#else

#ifndef SPARTA_COMPUTE_SURF_H
#define SPARTA_COMPUTE_SURF_H

#include "compute.h"
#include "surf.h"

namespace SPARTA_NS {

class ComputeSurf : public Compute {
 public:
  ComputeSurf(class SPARTA *, int, char **);
  ComputeSurf(class SPARTA* sparta) : Compute(sparta) {}
  ~ComputeSurf();
  virtual void init();
  void compute_per_surf();
  virtual void clear();
  virtual void surf_tally(int, Particle::OnePart *, 
                          Particle::OnePart *, Particle::OnePart *);
  virtual int tallyinfo(int *&);
  virtual void tallysum(int);
  bigint memory_usage();

 protected:
  int groupbit,imix,nvalue,ngroup,ntotal;
  int *which;
  bigint last_tallysum;    // last timestep tallysum was called

  int nsurf;               // # of lines/tris I own
                           // surf->nlocal+nghost for explicit all or distributed

  int ntally;              // # of surfs I have tallied for
  int maxtally;            // # of tallies currently allocated
  int *surf2tally;         // surf2tally[I] = tally index of Ith surf
  int *tally2surf;         // tally2surf[I] = surf index of Ith tally
  double **array;          // tally values, maxtally in length

  int dimension;           // local copies
  Surf::Line *lines;
  Surf::Tri *tris;

  int weightflag;          // 1 if cell weighting is enabled
  double weight;           // particle weight, based on initial cell
  double nfactor;          // dt/fnum for normalization
  double nfactor_inverse;  // fnum/dt for normalization
  double *normflux;        // normalization factor for each surf element
  double nfactor_previous; // nfactor from previous run
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

E: Compute surf mixture ID does not exist

Self-explanatory.

E: Number of groups in compute surf mixture has changed

This mixture property cannot be changed after this compute command is
issued.

*/
