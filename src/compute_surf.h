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
  ~ComputeSurf();
  void init();
  void compute_per_surf();
  void clear();
  void surf_tally(int, Particle::OnePart *, 
                  Particle::OnePart *, Particle::OnePart *);
  double *normptr(int);
  int surfinfo(int *&);
  bigint memory_usage();

 private:
  int groupbit,imix;
  int nvalue,ngroup,ntotal;
  int *which;

  int nsurf;               // # of global surfs, lines or triangles
  int nlocal;              // # of local surfs I have tallied for
  int maxlocal;            // # of local surfs currently allocated
  double **array;          // tally values for local surfs
  int *glob2loc;           // glob2loc[I] = local index of Ith global surf
  int *loc2glob;           // loc2glob[I] = global index of Ith local surf

  int dimension;           // local copies
  Surf::Line *lines;
  Surf::Tri *tris;

  int weightflag;          // 1 if cell weighting is enabled
  double weight;           // particle weight, based on initial cell
  double *normflux;        // normalization factor for all surf elements

  void grow();
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
