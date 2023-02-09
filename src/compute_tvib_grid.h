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

ComputeStyle(tvib/grid,ComputeTvibGrid)

#else

#ifndef SPARTA_COMPUTE_TVIB_GRID_H
#define SPARTA_COMPUTE_TVIB_GRID_H

#include "compute.h"

namespace SPARTA_NS {

class ComputeTvibGrid : public Compute {
 public:
  ComputeTvibGrid(class SPARTA *, int, char **);
  virtual ~ComputeTvibGrid();
  virtual void init();
  virtual void compute_per_grid();
  virtual int query_tally_grid(int, double **&, int *&);
  virtual void post_process_grid(int, int, double **, int *, double *, int);
  virtual void reallocate();
  bigint memory_usage();

 protected:
  int groupbit,imix,ngroup,mixspecies,nspecies;
  int modeflag;              // 1 when tallying stats for each vib mode
  int maxmode;               // max vib mode for any species

  int ntally;                // total # of columns in tally array
  int nglocal;               // # of owned grid cells
  double **tally;            // array of tally quantities, ncells by ntally

  int *nmap;               // nmap[i] = # of tally values for Ith output column
                           // size = # of outputs (Ngroup or Ngroup*Nmode)

  int **map;               // map[i][j] = tally columns for Ith output column
                           // size = # of outputs (Ngroup or Ngroup*Nmode)
                           //   by # of tally columns (2*nmax or 2*nmax*Nmode)
                           // nmax = max # of species in any group

  // modeflag = 0, tallying per grid cell and per species

  double *tspecies;          // per-species vibrational temps
                             // size = Nspecies

  int *s2t;                  // s2t[i] = first tally column for species I
                             // size = Nspecies
  int *t2s;                  // t2s[i] = species index for Ith tally column
                             // size = Ntally = 2*Nspecies

  // modeflag = 1 or 2, tallying per grid cell, per species, per mode

  int index_vibmode;         // index into extra particle values for vibmodes

  double **tspecies_mode;    // per-species per-mode vibrational temps
                             // size = Nspecies by Nmode
  int **s2t_mode;            // s2tmode[i][j] =
                             //   first tally column for species I, mode J
                             // length = Nspecies by Nmode
  int *t2s_mode;             // t2s_mode[i] = species index for Ith tally column
                             // size = Ntally = 2*Nspecies*Nmode
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Compute grid mixture ID does not exist

Self-explanatory.

E: Number of groups in compute tvib/grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

E: Number of species in compute tvib/grid mixture has changed

This mixture property cannot be changed after this compute command is
issued.

*/
