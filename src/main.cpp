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

/* --------------------------------
// Virgile - Modif Start - 07/11/24
-----------------------------------
The species weighting scheme was implemented as described in: "Conservative 
Species Weighting Scheme for the Direct Simulation Monte Carlo Method", 
Iain D. Boyd, Journal of Thermophysics and Heat Transfer, 1996. The
particles' momentum, kinetic energy and internal energies are conserved 
throughout the collisions using a splitting-merging approach.
Chemical reactions occurring when two differentially weighted particles
collide are implemented. As reactant are deleted, it is not possible
to use the splitting-merging approach. Instead, each reactant and product 
is associated with a probability to remain or be created according to their weights. 
-----------------------------------
To use the Species Weighting Scheme, use one of the following keywords 
when calling the 'species' command in SPARTA input files: SWS, SWSmax.
The keyword is read along with the species list to be used in the 
simulation. When using one of these keyword, the weight (0<wi<1) are
read for every species in the gas.specie file (9th column). Usinga a low weight
results in an increase in the number of numerical particles of the related specie.
-----------------------------------
The modifications are the following:
0) Warning message about the current SPARTA version in sparta.cpp
1) Particles emition/creation:
- fix_emit_face: fnum weighted with species wi
   - fnum weighted with species wi in 'perspecies yes' mode
   - accumulated weighted fraction in 'perspecies no' mode
- create_particles : particles are created according to the accumulative
weighted fraction instead of the classical accumulative fraction to 
account for species weight
- fix_emit_surf: 
   - fnum weighted with species wi in 'perspecies yes' mode
   - accumulated weighted fraction in 'perspecies no' mode
   - compute local quantity accounting for the species wi
2) Particles collision:
- grid: creation of a new cinfo variable named count_wi to store the sum 
of the wi per cell
- particle: 
   - compute the count_wi variable for every particle in cells
   - add the SWS method keywords: SWS and SWSmax
- collide: 
   - add the count_wi and count_wi_group variables to the nattempt function, as
   well as the wi maximum if option SWSmax is selected
   - initialize the lost energy due to different weighted particle collision
   to zero for each cells before performing the collisions
   - add the count of every particle that are created after chemical reaction occur
   to the perform_collide function
- collide_vss: 
   - modify the random process to select particle pairs accounting for species weights
   using two possible options: 1-no modification (SWS) ; 2-using wi maximum (SWSmax)
   - use the count_wi and count_wi_group variables instead of np to compute nattempt, as
   well as the weight maximum if option SWSmax is selected
   - modified two body scattering function to consider wi in the post collision 
   velocities (splitting - merging approach)
   - modified two body scattering function to track lost energy when
   different weighted particles collide and add to the next nontrace species collision
   - modified the energy exchange functions to consider wi in the final 
   internal energies (splitting - merging approach). If energies are descrete
   levels, prevent any exchange from the highest weighted particle when colliding
   with a trace specie.
   - modified the perform_collision function to create the post chemical reaction
   particles according to a probability based on their weight. As a consequence: 
      - the major reactant might remain or be consumed after the reaction
      - several particles of the minor product will be produced to reach mass conservation
      - major product might be created or not after the reaction
   Thus, each reaction does not explicitly conserve the mass. However, with a large number
   of collision across the flow, the global mass is conserved. The particles weights
   are always the same (depending on the species) and are never modified. Furthermore,
   while using this approach the splitting-merging approach is not performed
   as the probabilistic approach ensure the momentum and energy conservation.
3) Flow quantities computation:
- compute_grid: computation of the cell quantities accounting for the 
species wi. Add a keyword to output the wi sum in each cells: NUMWI (nwi in input file)
- compute_thermal_grid: computation of the cells P and T accounting for the 
species wi
- compute_sirf: computation of the surface quantities accounting for
the species wi.
-----------------------------------
// Virgile - Modif End - 07/11/24
----------------------------------- */

#include "mpi.h"
#include "sparta.h"
#include "input.h"

using namespace SPARTA_NS;

/* ----------------------------------------------------------------------
   main program to drive SPARTA
------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);

  SPARTA *sparta = new SPARTA(argc,argv,MPI_COMM_WORLD);
  sparta->input->file();
  delete sparta;

  MPI_Finalize();
}
