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

#ifdef FIX_CLASS

FixStyle(inflow/file,FixInflowFile)

#else

#ifndef SPARTA_FIX_INFLOW_FILE_H
#define SPARTA_FIX_INFLOW_FILE_H

#include "fix.h"

namespace SPARTA_NS {

class FixInflowFile : public Fix {
 public:
  FixInflowFile(class SPARTA *, int, char **);
  ~FixInflowFile();
  int setmask();
  void init();
  void start_of_step();
  int pack_grid_one(int, char *, int);
  int unpack_grid_one(int, char *);
  void compress_grid();
  void post_compress_grid();
  double compute_vector(int);

 private:
  int nevery,imix,face;
  int perspecies;
  int npercell,nthresh;
  int nsingle,ntotal;
  double frac_user;

  double normal[3];             // inward normal from external boundary face
  int ndim;                     // dim (0,1,2) normal to face
  int pdim,qdim;                // 2 dims (0,1,2) parallel to face

  struct Mesh {
    int ni,nj;                  // size of 2d grid, Nj = 1 for 1d grid
    int nvalues;                // # of values per grid point
    double lo[2],hi[2];         // bounds of mesh, in 2d/3d box coords
    int *which;                 // style of each value for which >= 0
                                // if negative, which = -(speciesIndex + 1)
    double *imesh;              // coordinate for each Ni
    double *jmesh;              // coordinate for each Nj
    double **values;            // 2d N,M values, where N = Ni*Nj
                                // I,J stored with J varying fastest in 1d N
  };

  Mesh mesh;

  struct CellFace {
    double lo[3];               // lower-left corner of overlap of cell/file
    double hi[3];               // upper-right corner of overlap of cell/file
    double ntarget;             // # of mols to insert for all species
    double *ntargetsp;          // # of mols to insert for each species
    cellint id;                 // ID of cell with insertion face
    int pcell;                  // associated cell index for particles
                                // unsplit or sub cell (not split cell)
    int icell;                  // associated cell index, unsplit or split cell

    // interpolated file values or defaults from mixture params

    double nrho;
    double vstream[3];
    double temp_thermal;
    double *fraction;
    double *cummulative;
    double *vscale;
  };

  CellFace *cellface;           // cell/face pairs to insert particles on
  int ncf;                      // # of cell/face pairs
  int ncfmax;                   // max # of cell/face pairs allocated

  int *c2f;                     // cellface index of cell I
                                // only for unsplit and split cells, not sub
                                // -1 if no insertions in that cell
  int nglocal;                  // # of owned grid cells
  int nglocalmax;               // max size of per-cell vectors/arrays

  class RanPark *random;

  void read_file(char *, char *);
  void bcast_mesh();
  void check_mesh_values();
  void interpolate();
  double linear_interpolation(double, int, int, int);
  double bilinear_interpolation(double, double, int, int, int, int, int);

  double mol_inflow(double, double, double);
  int split(int, int);
  void grow_percell(int);
  void grow_cellface(int);
  void print_face(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Fix inflow mixture ID does not exist

Self-explanatory.

E: Cannot use fix inflow in z dimension for 2d simulation

Self-explanatory.

E: Cannot use fix inflow in y dimension for axisymmetric

This is because the y dimension boundaries cannot be
inflow boundaries for an axisymmetric model.

E: Cannot use fix inflow n > 0 with perspecies yes

This is because the perspecies option calculates the
number of particles to insert itself.

E: Cannot use fix inflow on periodic boundary

Self-explanatory.

E: Fix inflow used on outflow boundary

Self-explanatory.

*/
