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

#ifdef FIX_CLASS

FixStyle(emit/face/file,FixEmitFaceFile)

#else

#ifndef SPARTA_FIX_EMIT_FACE_FILE_H
#define SPARTA_FIX_EMIT_FACE_FILE_H

#include "fix_emit.h"
#include "surf.h"
#include "grid.h"

namespace SPARTA_NS {

class FixEmitFaceFile : public FixEmit {
 public:
  FixEmitFaceFile(class SPARTA *, int, char **);
  ~FixEmitFaceFile();
  void init();

 private:
  int imix,iface,subsonic,subsonic_style,subsonic_warning;
  int npertask,nthresh;
  double frac_user;
  double tprefactor,soundspeed_mixture;

  double normal[3];             // inward normal from external boundary face
  int ndim;                     // dim (0,1,2) normal to face
  int pdim,qdim;                // 2 dims (0,1,2) parallel to face

  // copies of data from other classes

  int dimension,nspecies;
  double fnum,dt;
  double nrho_mix,temp_thermal_mix,temp_rot_mix,temp_vib_mix;
  double *vstream_mix,*vscale_mix,*fraction_mix;
  double *cummulative_mix,*fraction_user_mix;
  int *fraction_flag_mix,*species2species_mix;

  Surf::Line *lines;
  Surf::Tri *tris;

  // mesh data from input file

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

  // one insertion task for a cell and a face

  struct Task {
    double lo[3];               // lower-left corner of overlap of cell/file
    double hi[3];               // upper-right corner of overlap of cell/file
    double area;                // area of face
    double ntarget;             // # of mols to insert for all species
    double *ntargetsp;          // # of mols to insert for each species,
                                //   only defined for PERSPECIES

    int icell;                  // associated cell index, unsplit or split cell
    int pcell;                  // associated cell index for particles
                                // unsplit or sub cell (not split cell)

    // interpolated file values or defaults from mixture params

    double nrho;
    double temp_thermal,temp_rot,temp_vib;
    double press;
    double vstream[3];
    double *fraction;
    double *cummulative;
    double *vscale;
  };

                         // ntask = # of tasks is stored by parent class
  Task *tasks;           // list of particle insertion tasks
  int ntaskmax;          // max # of tasks allocated

  // active grid cells assigned to tasks, used by subsonic sorting

  int maxactive;
  int *activecell;

  // per-species vectors for species fractions on mesh

  int *fflag;
  double *fuser;

  // private methods

  void read_file(char *, char *);
  void bcast_mesh();
  void check_mesh_values();
  int interpolate(int);
  double linear_interpolation(double, int, int, int);
  double bilinear_interpolation(double, double, int, int, int, int, int);

  int split(int);

  void subsonic_inflow();
  void subsonic_sort();
  void subsonic_grid();

  void create_task(int);
  void perform_task();
  void grow_task();

  int option(int, char **);
  void print_task(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Fix emit/face/file mixture ID does not exist

Self-explanatory.

E: Cannot use emit/face/file in z dimension for 2d simulation

Self-explanatory.

E: Cannot use emit/face/file in y dimension for axisymmetric

This is because the y dimension boundaries cannot emit/face/file
boundaries for an axisymmetric model.

E: Cannot use emit/face/file n > 0 with perspecies yes

This is because the perspecies option calculates the
number of particles to insert itself.

E: Cannot use emit/face/file on periodic boundary

Self-explanatory.

E: emit/face/file used on outflow boundary

Self-explanatory.

*/
