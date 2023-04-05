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
  CPU time = 0.00109506 secs
  create/ghost percent = 75.3538 24.6462
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.00026226 secs
  reassign/sort/migrate/ghost percent = 46.8182 0.545455 15.1818 37.4545

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
  CPU time = 0.00087142 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 22.4624 7.08618 2.02462 54.829 13.5978 8.29001 0.300958
  surf2grid time = 0.000477791 secs
  map/comm1/comm2/comm3/comm4/split percent = 42.7146 9.08184 16.2176 3.34331 12.1756 12.5749

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
     100  0.064087629    18931        0        0      119     3949 
     200   0.19977212    28398        0        0      163     5567 
     300    0.3678782    32180        0        0      191     6532 
     400   0.55167007    33853        0        0      188     6727 
     500   0.74070859    34643        0        0      181     6877 
     600    0.9306488    35281        0        0      177     6866 
     700    1.1230724    35551        0        0      179     6622 
     800    1.3171136    35458        0        0      189     7107 
     900    1.5100048    35446        0        0      197     6725 
    1000    1.7047114    35643        0        0      188     6802 
Loop time of 1.70474 on 1 procs for 1000 steps with 35643 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.2494     | 1.2494     | 1.2494     |   0.0 | 73.29
Coll    | 0.145      | 0.145      | 0.145      |   0.0 |  8.51
Sort    | 0.17905    | 0.17905    | 0.17905    |   0.0 | 10.50
Comm    | 0.0067759  | 0.0067759  | 0.0067759  |   0.0 |  0.40
Modify  | 0.12198    | 0.12198    | 0.12198    |   0.0 |  7.16
Output  | 0.00025153 | 0.00025153 | 0.00025153 |   0.0 |  0.01
Other   |            | 0.00224    |            |       |  0.13

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

Particle-moves/CPUsec/proc: 1.82478e+07
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

