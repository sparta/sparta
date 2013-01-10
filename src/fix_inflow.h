/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@sandia.gov, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2012) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(inflow,FixInflow)

#else

#ifndef SPARTA_FIX_INFLOW_H
#define SPARTA_FIX_INFLOW_H

#include "stdio.h"
#include "fix.h"
#include "create_molecules.h"

namespace SPARTA_NS {

class FixInflow : public Fix {
 public:
  FixInflow(class SPARTA *, int, char **);
  ~FixInflow();
  int setmask();
  void init();
  void start_of_step();
  double compute_vector(int);

 private:
  int nevery,imix;
  int faces[6];
  int np,perspecies;
  int npercell,nthresh;
  int nsingle,ntotal;

  struct CellFace {
    double lo[3];               // lower-left corner of face
    double hi[3];               // upper-right corner of face
    double normal[3];           // inward normal from external boundary
    double ntarget;             // # of mols to insert for all species
    double *ntargetsp;          // # of mols to insert for each species
    int icell;                  // associated cell index, may be split cell
    int ndim;                   // dim (0,1,2) normal to face
    int pdim1,pdim2;            // 2 dims (0,1,2) parallel to face
  };

  CellFace *cellface;           // cell/face pairs to insert particles on
  int ncf;                      // # of cell/face pairs
  int ncfmax;                   // max # of cell/face pairs for allocation

  class RanPark *random;
  class CreateMolecules;

  double mol_inflow(int, double);
  int split(int, int);
};

}

#endif
#endif
