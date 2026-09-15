SPARTA (24 Sep 2025)
Running on 1 MPI task(s)
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
Created 200 child grid cells
  CPU time = 0.0025411 secs
  create/ghost percent = 86.8522 13.1478
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.00068121 secs
  reassign/sort/migrate/ghost percent = 93.3596 0.222105 2.50716 3.91113

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
     100  0.059838157    21190      527      346        0       20 
     200   0.16533822    30189      809      487        0      220 
     300   0.28340559    30380      836      478        0      202 
     400   0.39768704    30154      848      497        0      228 
     500   0.51027475    30335      866      488        0      238 
Loop time of 0.510284 on 1 procs for 500 steps with 30335 particles
Performance: 979.847 timesteps/s, 29.724 Mparticle-step/s

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.20154    | 0.20154    | 0.20154    |   0.0 | 39.50
Coll    | 0.11443    | 0.11443    | 0.11443    |   0.0 | 22.43
Sort    | 0.028602   | 0.028602   | 0.028602   |   0.0 |  5.61
Comm    | 0.069411   | 0.069411   | 0.069411   |   0.0 | 13.60
Modify  | 0.095511   | 0.095511   | 0.095511   |   0.0 | 18.72
Output  | 7.8399e-05 | 7.8399e-05 | 7.8399e-05 |   0.0 |  0.02
MPI Sync| 0.0006117  | 0.0006117  | 0.0006117  |   0.0 |  0.12
Other   |            | 9.113e-05  |            |       |  0.02

Particle moves    = 13046351 (13M)
Cells touched (std move) = 0 (0K)
Particle comms    = 0 (0K)
Boundary collides = 0 (0K)
Boundary exits    = 78279 (78.3K)
SurfColl checks   = 0 (0K)
SurfColl occurs   = 0 (0K)
Surf reactions    = 0 (0K)
Collide attempts  = 350395 (0.35M)
Collide occurs    = 212336 (0.212M)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.55669e+07
Particle-moves/step: 26092.7
Cell-touches/particle/step (std move): 0
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0
Particle fraction exiting boundary: 0.00600007
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0.0268577
Collisions/particle/step: 0.0162755
Reactions/particle/step: 0

Particles: 30335 ave 30335 max 30335 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      200 ave 200 max 200 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
