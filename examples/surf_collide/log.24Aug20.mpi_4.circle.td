SPARTA (6 Jul 2020)
################################################################################
# 2d flow around a circle
#
# Note:
#  - The "comm/sort” option to the “global” command is used to match MPI runs.
#  - The “twopass” option is used to match Kokkos runs.
# The "comm/sort" and "twopass" options should not be used for production runs.
################################################################################

seed	    	    12345
dimension   	    2
global              gridcut 0.0 comm/sort yes

boundary	    o r p

create_box  	    0 10 0 10 -0.5 0.5
Created orthogonal box = (0 0 -0.5) to (10 10 0.5)
create_grid 	    20 20 1
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/Users/eharvey/dev/SPARTA.base/sparta/src/grid.cpp:415)
Created 400 child grid cells
  parent cells = 1
  CPU time = 0.000919 secs
  create/ghost percent = 95.9739 4.02612
balance_grid        rcb cell
Balance grid migrated 280 cells
  CPU time = 0.001072 secs
  reassign/sort/migrate/ghost percent = 51.6791 0.373134 38.9925 8.95522

global		    nrho 1.0 fnum 0.001

species		    air.species N O
mixture		    air N O vstream 100.0 0 0

read_surf           data.circle
  50 points
  50 lines
  2 8 xlo xhi
  2.00592 7.99408 ylo yhi
  0 0 zlo zhi
  0.376743 min line length
  48 = cells with surfs
  104 = total surfs in all grid cells
  3 = max surfs in one grid cell
  0.753486 = min surf-size/cell-size ratio
  0 0 = number of pushed cells
  48 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  264 88 48 = cells outside/inside/overlapping surfs
  48 = surf cells with 1,2,etc splits
  71.8 71.8 = cell-wise and global flow volume
  CPU time = 0.000459 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 44.2266 9.36819 0.43573 30.9368 15.0327 10.0218 0.217865
  surf2grid time = 0.000142 secs
  map/rvous1/rvous2/split percent = 17.6056 52.1127 0.704225 16.9014

surf_collide	    1 td 1000.0 #barrier 1000
#surf_collide	    1 td 500.0 bond 500 0 0
#surf_collide	    1 td 300.0 initenergy 0.01 0 0
#surf_collide	    1 td 1000.0 barrier 1000 bond 500 0 0 initenergy 0.01 0 0

surf_modify         all collide 1

collide             vss air air.vss

fix		    in emit/face air xlo twopass # subsonic 0.1 NULL

timestep 	    0.0001

#dump                2 image all 50 image.*.ppm type type pdiam 0.1 #                    surf proc 0.01 size 512 512 zoom 1.75 #                    gline yes 0.005
#dump_modify	    2 pad 4

stats		    100
stats_style	    step cpu np nattempt ncoll nscoll nscheck
run 		    1000
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 1.51903 1.51903 1.51903
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
     100     0.015629    18821        0        0      131     3964 
     200     0.042225    28388        0        0      161     5480 
     300     0.077294    32168        0        0      180     6421 
     400     0.116272    33644        0        0      208     6677 
     500     0.152136    34631        0        0      215     7019 
     600     0.188613    35115        0        0      190     6972 
     700     0.223317    35271        0        0      213     6917 
     800     0.257314    35496        0        0      226     7088 
     900      0.29324    35565        0        0      197     6971 
    1000     0.331049    35708        0        0      194     6943 
Loop time of 0.331064 on 4 procs for 1000 steps with 35708 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.050217   | 0.13874    | 0.22871    |  23.7 | 41.91
Coll    | 0.006363   | 0.019026   | 0.031866   |   9.2 |  5.75
Sort    | 0.006613   | 0.017079   | 0.027621   |   8.0 |  5.16
Comm    | 0.013565   | 0.014316   | 0.014948   |   0.4 |  4.32
Modify  | 0.000153   | 0.011072   | 0.022175   |  10.4 |  3.34
Output  | 8.5e-05    | 0.0003305  | 0.000549   |   0.0 |  0.10
Other   |            | 0.1305     |            |       | 39.42

Particle moves    = 31073377 (31.1M)
Cells touched     = 35957079 (36M)
Particle comms    = 144128 (0.144M)
Boundary collides = 172051 (0.172M)
Boundary exits    = 174877 (0.175M)
SurfColl checks   = 6056814 (6.06M)
SurfColl occurs   = 173137 (0.173M)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 2.34648e+07
Particle-moves/step: 31073.4
Cell-touches/particle/step: 1.15717
Particle comm iterations/step: 2.13
Particle fraction communicated: 0.00463831
Particle fraction colliding with boundary: 0.00553693
Particle fraction exiting boundary: 0.00562787
Surface-checks/particle/step: 0.19492
Surface-collisions/particle/step: 0.00557188
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 8927 ave 14061 max 3743 min
Histogram: 2 0 0 0 0 0 0 0 0 2
Cells:      100 ave 100 max 100 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 21 ave 21 max 21 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 21 ave 21 max 21 min
Histogram: 4 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0

