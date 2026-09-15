SPARTA (24 Sep 2025)
Running on 4 MPI task(s)
################################################################################
# 2d axisymmetric flow through a duct
#
# Exercises the optimized move on an axisymmetric model.  The fast path folds
# the end-of-step position back into the axisymmetric plane once for the whole
# timestep, where the normal move does it at every cell crossing; both give the
# same point, since each remap is a rotation applied to the position and the
# velocity together and the radial coordinate is unchanged by it.
#
# The axis itself needs no boundary condition, and an outflow face at the
# outer radius is handled on the fast path.  A reflecting outer radial face
# would not be, since that face is a cylinder -- see the global doc page.
#
# Run the same script with "-var opt no" to compare against the normal move.
# The optimized move performs one remap per timestep where the normal move
# performs one per cell crossing, so the two are identical while a particle
# stays within a cell over a timestep, as it does at the timestep used here,
# and differ in the last few digits once particles cross cells.
#
# Note:
#  - The "comm/sort" option to the "global" command is used to match MPI runs.
#  - The "twopass" option is used to match Kokkos runs.
################################################################################

variable            opt index yes

seed	    	    12345
dimension   	    2
global              gridcut 0.0 comm/sort yes optmove ${opt}
global              gridcut 0.0 comm/sort yes optmove yes

boundary	    o ao p

create_box          -0.25 0.25 0.0 0.25 -0.5 0.5
Created orthogonal box = (-0.25 0 -0.5) to (0.25 0.25 0.5)
create_grid 	    20 10 1
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (../grid.cpp:483)
Created 200 child grid cells
  CPU time = 0.000919695 secs
  create/ghost percent = 90.885 9.11498
balance_grid        rcb cell
Balance grid migrated 120 cells
  CPU time = 0.000305147 secs
  reassign/sort/migrate/ghost percent = 66.4008 0.575133 14.5795 18.4446

global		    nrho 1.e20 fnum 1.e17 weight cell radius

species		    air.species N2
mixture		    air N2 vstream 3472.0 0.0 0.0 temp 300.0

fix                 in emit/face air xlo twopass
collide		    vss air air.vss

fix                 gcheck grid/check 1 error

stats		    100
stats_style	    step cpu np nattempt ncoll nbound nexit

timestep 	    1.e-6
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0 0 0
  modify    (ave,min,max) = 0 0 0
  total     (ave,min,max) = 1.51379 1.51379 1.51379
Step CPU Np Natt Ncoll Nbound Nexit 
       0            0        0        0        0        0        0 
     100  0.029197611    21360      529      315        0       20 
     200  0.059235925    30241      826      514        0      230 
     300  0.088868638    30286      847      500        0      230 
     400   0.12015281    30291      858      525        0      244 
     500   0.15141737    30238      866      540        0      213 
Loop time of 0.151476 on 4 procs for 500 steps with 30238 particles
Performance: 3300.854 timesteps/s, 99.811 Mparticle-step/s

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.034843   | 0.044605   | 0.051494   |   3.4 | 29.45
Coll    | 0.017529   | 0.021643   | 0.024508   |   2.0 | 14.29
Sort    | 0.0037611  | 0.0046732  | 0.0052989  |   0.9 |  3.09
Comm    | 0.023607   | 0.025035   | 0.026027   |   0.6 | 16.53
Modify  | 0.012604   | 0.020111   | 0.036191   |   6.7 | 13.28
Output  | 0.00014066 | 0.00023063 | 0.000335   |   0.0 |  0.15
MPI Sync| 0.0082152  | 0.035136   | 0.05879    |  10.3 | 23.20
Other   |            | 4.197e-05  |            |       |  0.03

Particle moves    = 13061213 (13.1M)
Cells touched (std move) = 0 (0K)
Particle comms    = 270435 (0.27M)
Boundary collides = 0 (0K)
Boundary exits    = 78364 (78.4K)
SurfColl checks   = 0 (0K)
SurfColl occurs   = 0 (0K)
Surf reactions    = 0 (0K)
Collide attempts  = 352741 (0.353M)
Collide occurs    = 213525 (0.214M)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.15566e+07
Particle-moves/step: 26122.4
Cell-touches/particle/step (std move): 0
Particle comm iterations/step: 1
Particle fraction communicated: 0.0207052
Particle fraction colliding with boundary: 0
Particle fraction exiting boundary: 0.00599975
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0.0270068
Collisions/particle/step: 0.016348
Reactions/particle/step: 0

Particles: 7559.5 ave 7776 max 7332 min
Histogram: 1 0 0 1 0 0 0 1 0 1
Cells:      50 ave 50 max 50 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 15 ave 20 max 10 min
Histogram: 2 0 0 0 0 0 0 0 0 2
EmptyCell: 15 ave 20 max 10 min
Histogram: 2 0 0 0 0 0 0 0 0 2
