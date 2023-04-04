SPARTA (18 Jul 2022)
Running on 1 MPI task(s)
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
Created 400 child grid cells
  CPU time = 0.00139856 secs
  create/ghost percent = 58.0975 41.9025
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000625849 secs
  reassign/sort/migrate/ghost percent = 25.2952 0.952381 12.7619 60.9905

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
  CPU time = 0.00167084 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 11.3584 7.7911 3.496 53.8385 23.516 7.09189 0.0142694
  surf2grid time = 0.000899553 secs
  map/comm1/comm2/comm3/comm4/split percent = 45.375 5.56586 12.1919 3.207 6.91757 21.7069

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
     100    0.1165359    18931        0        0      119     3949 
     200   0.37285781    28398        0        0      163     5567 
     300   0.68632507    32180        0        0      191     6532 
     400    1.0278928    33853        0        0      188     6727 
     500    1.3759494    34643        0        0      181     6877 
     600    1.7299268    35281        0        0      177     6866 
     700    2.0884409    35551        0        0      179     6622 
     800    2.4468369    35458        0        0      189     7107 
     900    2.8056369    35446        0        0      197     6725 
    1000    3.1652975    35643        0        0      188     6802 
Loop time of 3.16532 on 1 procs for 1000 steps with 35643 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 2.3737     | 2.3737     | 2.3737     |   0.0 | 74.99
Coll    | 0.25328    | 0.25328    | 0.25328    |   0.0 |  8.00
Sort    | 0.33178    | 0.33178    | 0.33178    |   0.0 | 10.48
Comm    | 0.0081172  | 0.0081172  | 0.0081172  |   0.0 |  0.26
Modify  | 0.19581    | 0.19581    | 0.19581    |   0.0 |  6.19
Output  | 0.00021124 | 0.00021124 | 0.00021124 |   0.0 |  0.01
Other   |            | 0.002424   |            |       |  0.08

Particle moves    = 31107764 (31.1M)
Cells touched     = 35994224 (36M)
Particle comms    = 0 (0K)
Boundary collides = 172472 (0.172M)
Boundary exits    = 175080 (0.175M)
SurfColl checks   = 6067907 (6.07M)
SurfColl occurs   = 173578 (0.174M)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 9.82769e+06
Particle-moves/step: 31107.8
Cell-touches/particle/step: 1.15708
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00554434
Particle fraction exiting boundary: 0.00562818
Surface-checks/particle/step: 0.195061
Surface-collisions/particle/step: 0.00557989
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 35643 ave 35643 max 35643 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      400 ave 400 max 400 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0

