SPARTA (18 Jul 2022)
Running on 4 MPI task(s)
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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/ascldap/users/stamoor/sparta_master/src/grid.cpp:465)
Created 400 child grid cells
  CPU time = 0.00112939 secs
  create/ghost percent = 84.0828 15.9172
balance_grid        rcb cell
Balance grid migrated 280 cells
  CPU time = 0.00112247 secs
  reassign/sort/migrate/ghost percent = 58.9635 0.637213 14.8683 25.531

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
  48 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  264 88 48 = cells outside/inside/overlapping surfs
  48 = surf cells with 1,2,etc splits
  71.8 71.8 = cell-wise and global flow volume
  CPU time = 0.00152159 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 20.3385 19.2259 1.7236 42.7139 15.9981 18.4738 0.125353
  surf2grid time = 0.000649929 secs
  map/comm1/comm2/comm3/comm4/split percent = 38.4813 7.85033 8.43727 3.70506 10.3081 23.8078

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
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 1.51894 1.51894 1.51894
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
     100  0.065774679    18935        0        0      112     3801 
     200   0.18979454    28432        0        0      161     5873 
     300   0.33184338    32176        0        0      183     6306 
     400   0.48120832    33776        0        0      187     6687 
     500   0.63368487    34642        0        0      215     6888 
     600   0.78819609    35095        0        0      140     6581 
     700   0.95245671    35382        0        0      176     6848 
     800      1.10883    35523        0        0      206     7186 
     900     1.265111    35651        0        0      186     6793 
    1000    1.4215178    35732        0        0      171     6942 
Loop time of 1.4216 on 4 procs for 1000 steps with 35732 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.22118    | 0.62435    | 1.0346     |  51.0 | 43.92
Coll    | 0.022651   | 0.05895    | 0.096079   |  14.9 |  4.15
Sort    | 0.034537   | 0.086333   | 0.13875    |  17.6 |  6.07
Comm    | 0.027684   | 0.030005   | 0.032305   |   1.0 |  2.11
Modify  | 0.00075531 | 0.052347   | 0.10396    |  22.5 |  3.68
Output  | 0.00031352 | 0.0013548  | 0.0021253  |   2.2 |  0.10
Other   |            | 0.5683     |            |       | 39.97

Particle moves    = 31095407 (31.1M)
Cells touched     = 35983645 (36M)
Particle comms    = 143676 (0.144M)
Boundary collides = 172693 (0.173M)
Boundary exits    = 175102 (0.175M)
SurfColl checks   = 6041072 (6.04M)
SurfColl occurs   = 173166 (0.173M)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 5.46838e+06
Particle-moves/step: 31095.4
Cell-touches/particle/step: 1.1572
Particle comm iterations/step: 2.145
Particle fraction communicated: 0.00462049
Particle fraction colliding with boundary: 0.00555365
Particle fraction exiting boundary: 0.00563112
Surface-checks/particle/step: 0.194275
Surface-collisions/particle/step: 0.00556886
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 8933 ave 14091 max 3766 min
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

