SPARTA (18 Jul 2022)
Running on 4 MPI task(s)
################################################################################
# beam of particles striking the surface at an inclined angle - free molecular flow (no collisions)
#
# Note:
#  - The "comm/sort” option to the “global” command is used to match MPI runs.
#  - The “twopass” option is used to match Kokkos runs.
# The "comm/sort" and "twopass" options should not be used for production runs.
################################################################################

seed	    	    123456
dimension   	    3
global              gridcut 0.0 comm/sort yes

boundary	    oo oo so

create_box          -11 11 -11 11 0 10
Created orthogonal box = (-11 -11 0) to (11 11 10)
create_grid 	    1 1 1
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/me/sparta_master/src/grid.cpp:465)
Created 1 child grid cells
  CPU time = 0.00144802 secs
  create/ghost percent = 89.4133 10.5867
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000520223 secs
  reassign/sort/migrate/ghost percent = 80.7315 0.544574 7.61423 11.1097

global		    nrho 1e10 fnum 1e6

species		    air.species N O
mixture		    air N O vstream 0 1000 -1000

mixture             air N frac 0.8
mixture             air O frac 0.2

surf_collide	    1 td 1000.0 #barrier 1000
#surf_collide	    1 td 500.0 bond 500 0 0
#surf_collide	    1 td 300.0 initenergy 0.01 0 0
#surf_collide	    1 td 1000.0 barrier 1000 bond 500 0 0 initenergy 0.01 0 0

bound_modify        zlo collide 1

region              circle cylinder z 0 -10 1 -INF INF
fix                 in emit/face/file air zhi data.beam beam_area nevery 100 region circle

#dump                2 image all 10 image.*.ppm type type pdiam 0.2 surf proc 0.01 size 512 512 zoom 1.75 gline no 0.005
#dump_modify	    2 pad 4

timestep            0.0001

stats		    100
stats_style	    step cpu np nattempt ncoll nscoll nscheck
run 		    1000
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 0.435669 0.0762939 1.51379
  surf      (ave,min,max) = 0 0 0
  total     (ave,min,max) = 0.435669 0.0762939 1.51379
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
     100  0.004208003     3118        0        0        0        0 
     200  0.016700852     6225        0        0        0        0 
     300  0.036303551     7461        0        0        0        0 
     400  0.057615133     7628        0        0        0        0 
     500  0.078811457     7682        0        0        0        0 
     600   0.10026978     7701        0        0        0        0 
     700   0.12302557     7763        0        0        0        0 
     800   0.14525355     7772        0        0        0        0 
     900   0.16690667     7847        0        0        0        0 
    1000   0.18864097     7765        0        0        0        0 
Loop time of 0.188738 on 4 procs for 1000 steps with 7765 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.00024008 | 0.03957    | 0.15744    |  34.2 | 20.97
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0044474  | 0.004919   | 0.005633   |   0.6 |  2.61
Modify  | 0.00015485 | 0.0052713  | 0.020528   |  12.1 |  2.79
Output  | 0.00020677 | 0.00041307 | 0.0010213  |   0.0 |  0.22
Other   |            | 0.1386     |            |       | 73.42

Particle moves    = 5391255 (5.39M)
Cells touched     = 5391255 (5.39M)
Particle comms    = 0 (0K)
Boundary collides = 28030 (28K)
Boundary exits    = 23405 (23.4K)
SurfColl checks   = 0 (0K)
SurfColl occurs   = 0 (0K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 7.14121e+06
Particle-moves/step: 5391.26
Cell-touches/particle/step: 1
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00519916
Particle fraction exiting boundary: 0.00434129
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 1941.25 ave 7765 max 0 min
Histogram: 3 0 0 0 0 0 0 0 0 1
Cells:      0.25 ave 1 max 0 min
Histogram: 3 0 0 0 0 0 0 0 0 1
GhostCell: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0

