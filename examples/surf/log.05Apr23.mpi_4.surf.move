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
create_grid 	    10 10 1
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/me/sparta_master/src/grid.cpp:465)
Created 100 child grid cells
  CPU time = 0.00152085 secs
  create/ghost percent = 85.7314 14.2686
balance_grid        rcb cell
Balance grid migrated 74 cells
  CPU time = 0.000733468 secs
  reassign/sort/migrate/ghost percent = 66.9934 0.545218 11.9299 20.5315

global		    nrho 1.0 fnum 0.001

species		    air.species N O
mixture		    air N O vstream 100.0 0 0

read_surf           data.circle origin 5 5 0 trans 0.0 2.0 0.0                     scale 0.33 0.33 1 group 1
  50 points
  50 lines
  4.01 5.99 xlo xhi
  6.01195 7.98805 ylo yhi
  0 0 zlo zhi
  0.124325 min line length
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  CPU time = 0.00127457 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 28.6595 22.5617 0.488558 37.9244 10.3659 21.9275 0.286215
  surf2grid time = 0.000483371 secs
  map/comm1/comm2/comm3/comm4/split percent = 29.8042 11.3881 7.51286 6.69444 16.0934 22.931
surf_collide	    1 diffuse 300.0 0.0
surf_modify         all collide 1

collide             vss air air.vss

fix		    in emit/face air xlo twopass

timestep 	    0.0001

#dump                2 image all 50 image.*.ppm type type pdiam 0.1 #                    surf proc 0.01 size 512 512 zoom 1.75
#dump_modify	    2 pad 4

stats		    100
stats_style	    step cpu np nattempt ncoll nscoll nscheck
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 1.51894 1.51894 1.51894
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
     100   0.03967462    20901        0        0       14     3122 
     200   0.13038062    35945        0        0       53     6650 
     300   0.25435578    43816        0        0       47     7812 
     400   0.39795286    47999        0        0       65     8414 
     500   0.54902848    50612        0        0       62     9296 
Loop time of 0.549147 on 4 procs for 500 steps with 50612 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.10209    | 0.17649    | 0.29004    |  18.1 | 32.14
Coll    | 0.014345   | 0.025585   | 0.035715   |   6.1 |  4.66
Sort    | 0.071802   | 0.10636    | 0.13092    |   7.8 | 19.37
Comm    | 0.04156    | 0.043517   | 0.04516    |   0.7 |  7.92
Modify  | 0.00076625 | 0.017454   | 0.035296   |  12.5 |  3.18
Output  | 0.00017346 | 0.00065979 | 0.0010097  |   0.0 |  0.12
Other   |            | 0.1791     |            |       | 32.61

Particle moves    = 17632049 (17.6M)
Cells touched     = 18843631 (18.8M)
Particle comms    = 122668 (0.123M)
Boundary collides = 59486 (59.5K)
Boundary exits    = 54770 (54.8K)
SurfColl checks   = 3034038 (3.03M)
SurfColl occurs   = 21569 (21.6K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 8.02701e+06
Particle-moves/step: 35264.1
Cell-touches/particle/step: 1.06871
Particle comm iterations/step: 2.366
Particle fraction communicated: 0.0069571
Particle fraction colliding with boundary: 0.00337374
Particle fraction exiting boundary: 0.00310628
Surface-checks/particle/step: 0.172075
Surface-collisions/particle/step: 0.00122328
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 12653 ave 16390 max 7983 min
Histogram: 1 0 0 1 0 0 0 0 0 2
Cells:      25 ave 25 max 25 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 11 ave 11 max 11 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 11 ave 11 max 11 min
Histogram: 4 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0

move_surf           all trans -1 0 0
Moving surfs ...
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1902 deleted particles
  CPU time = 0.00157156 secs
  sort/surf2grid/ghost/inout/particle percent = 22.0881 38.3849 6.08528 19.2827 14.159
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.10938 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 3.62832 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48710        0        0        0        0 
     600   0.14635559    50401        0        0       71     9800 
     700   0.29703804    51556        0        0       68     8722 
     800   0.45098572    52357        0        0       65     9576 
     900   0.60859773    52699        0        0       53     9072 
    1000   0.77155598    53033        0        0       65     9282 
Loop time of 0.771648 on 4 procs for 500 steps with 53033 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.17388    | 0.27175    | 0.42131    |  18.2 | 35.22
Coll    | 0.038089   | 0.046033   | 0.054735   |   3.7 |  5.97
Sort    | 0.13267    | 0.16383    | 0.18833    |   5.7 | 21.23
Comm    | 0.053151   | 0.055295   | 0.057185   |   0.8 |  7.17
Modify  | 0.00084997 | 0.018674   | 0.037007   |  13.0 |  2.42
Output  | 0.00019773 | 0.00068807 | 0.0010805  |   0.0 |  0.09
Other   |            | 0.2154     |            |       | 27.91

Particle moves    = 25898484 (25.9M)
Cells touched     = 27528828 (27.5M)
Particle comms    = 170023 (0.17M)
Boundary collides = 85994 (86K)
Boundary exits    = 101007 (0.101M)
SurfColl checks   = 4695502 (4.7M)
SurfColl occurs   = 30862 (30.9K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 8.39064e+06
Particle-moves/step: 51797
Cell-touches/particle/step: 1.06295
Particle comm iterations/step: 2.438
Particle fraction communicated: 0.00656498
Particle fraction colliding with boundary: 0.00332043
Particle fraction exiting boundary: 0.00390011
Surface-checks/particle/step: 0.181304
Surface-collisions/particle/step: 0.00119165
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 13258.2 ave 16254 max 9595 min
Histogram: 1 0 1 0 0 0 0 0 0 2
Cells:      25 ave 25 max 25 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 11 ave 11 max 11 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 11 ave 11 max 11 min
Histogram: 4 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0

move_surf           all trans 0 -1 0
Moving surfs ...
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1444 deleted particles
  CPU time = 0.00191881 secs
  sort/surf2grid/ghost/inout/particle percent = 20.4445 25.113 31.4869 10.6025 12.3532
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51589        0        0        0        0 
    1100   0.14812139    52493        0        0       48     8820 
    1200    0.3049943    53385        0        0       60     9618 
    1300   0.46692381    53890        0        0       70     9856 
    1400   0.62636709    54266        0        0       63    10094 
    1500   0.78629387    54428        0        0       64    10164 
Loop time of 0.786477 on 4 procs for 500 steps with 54428 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.18847    | 0.28145    | 0.42855    |  17.8 | 35.79
Coll    | 0.039646   | 0.048892   | 0.059623   |   3.9 |  6.22
Sort    | 0.14227    | 0.17192    | 0.19923    |   5.6 | 21.86
Comm    | 0.055905   | 0.058277   | 0.061235   |   0.9 |  7.41
Modify  | 0.00088006 | 0.017985   | 0.036202   |  12.7 |  2.29
Output  | 0.00019609 | 0.00073037 | 0.0012025  |   0.0 |  0.09
Other   |            | 0.2072     |            |       | 26.35

Particle moves    = 26791563 (26.8M)
Cells touched     = 28436601 (28.4M)
Particle comms    = 171396 (0.171M)
Boundary collides = 87221 (87.2K)
Boundary exits    = 102500 (0.102M)
SurfColl checks   = 4773216 (4.77M)
SurfColl occurs   = 30765 (30.8K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 8.51632e+06
Particle-moves/step: 53583.1
Cell-touches/particle/step: 1.0614
Particle comm iterations/step: 2.73
Particle fraction communicated: 0.00639739
Particle fraction colliding with boundary: 0.00325554
Particle fraction exiting boundary: 0.00382583
Surface-checks/particle/step: 0.178161
Surface-collisions/particle/step: 0.00114831
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 13607 ave 17039 max 10130 min
Histogram: 1 1 0 0 0 0 0 0 1 1
Cells:      25 ave 25 max 25 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 11 ave 11 max 11 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 11 ave 11 max 11 min
Histogram: 4 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0

move_surf           all trans 1 0 0
Moving surfs ...
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  696 deleted particles
  CPU time = 0.00139582 secs
  sort/surf2grid/ghost/inout/particle percent = 26.9879 30.4252 6.89279 20.008 15.6862
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53732        0        0        0        0 
    1600    0.1569103    54687        0        0       63     9464 
    1700   0.32225287    55802        0        0       72     9870 
    1800   0.48670232    56474        0        0       59     9408 
    1900   0.65486919    56965        0        0       57     9730 
    2000   0.82186047    57244        0        0       58     9618 
Loop time of 0.822064 on 4 procs for 500 steps with 57244 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.21723    | 0.29486    | 0.42531    |  15.7 | 35.87
Coll    | 0.041736   | 0.05279    | 0.065854   |   4.6 |  6.42
Sort    | 0.1439     | 0.18474    | 0.21579    |   7.2 | 22.47
Comm    | 0.0569     | 0.05952    | 0.061436   |   0.7 |  7.24
Modify  | 0.0010154  | 0.018875   | 0.036817   |  13.0 |  2.30
Output  | 0.00021718 | 0.00080426 | 0.0010727  |   0.0 |  0.10
Other   |            | 0.2105     |            |       | 25.60

Particle moves    = 28039815 (28M)
Cells touched     = 29761957 (29.8M)
Particle comms    = 176360 (0.176M)
Boundary collides = 92145 (92.1K)
Boundary exits    = 101797 (0.102M)
SurfColl checks   = 4657814 (4.66M)
SurfColl occurs   = 29999 (30K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 8.52726e+06
Particle-moves/step: 56079.6
Cell-touches/particle/step: 1.06142
Particle comm iterations/step: 2.868
Particle fraction communicated: 0.00628963
Particle fraction colliding with boundary: 0.00328622
Particle fraction exiting boundary: 0.00363044
Surface-checks/particle/step: 0.166114
Surface-collisions/particle/step: 0.00106987
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 14311 ave 17948 max 9821 min
Histogram: 1 0 1 0 0 0 0 0 0 2
Cells:      25 ave 25 max 25 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 11 ave 11 max 11 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 11 ave 11 max 11 min
Histogram: 4 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0

move_surf           all trans 0 1 0
Moving surfs ...
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1463 deleted particles
  CPU time = 0.00143013 secs
  sort/surf2grid/ghost/inout/particle percent = 27.618 29.5918 6.61777 18.4164 17.756
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55781        0        0        0        0 
    2100    0.1623245    56408        0        0       59     9408 
    2200   0.33468855    57315        0        0       61     9212 
    2300   0.50914592    57640        0        0       75     9954 
    2400    0.6797455    57767        0        0       70    10290 
    2500   0.85132494    57770        0        0       69    10584 
Loop time of 0.851555 on 4 procs for 500 steps with 57770 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.21715    | 0.30033    | 0.44422    |  16.7 | 35.27
Coll    | 0.041188   | 0.055373   | 0.068456   |   4.7 |  6.50
Sort    | 0.1426     | 0.19052    | 0.22526    |   7.9 | 22.37
Comm    | 0.059122   | 0.062697   | 0.064904   |   0.9 |  7.36
Modify  | 0.00097075 | 0.017768   | 0.03608    |  12.6 |  2.09
Output  | 0.00017055 | 0.00086507 | 0.0012997  |   0.0 |  0.10
Other   |            | 0.224      |            |       | 26.30

Particle moves    = 28683260 (28.7M)
Cells touched     = 30432587 (30.4M)
Particle comms    = 178617 (0.179M)
Boundary collides = 94726 (94.7K)
Boundary exits    = 103376 (0.103M)
SurfColl checks   = 4808622 (4.81M)
SurfColl occurs   = 30848 (30.8K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 8.42085e+06
Particle-moves/step: 57366.5
Cell-touches/particle/step: 1.06099
Particle comm iterations/step: 2.52
Particle fraction communicated: 0.00622722
Particle fraction colliding with boundary: 0.00330248
Particle fraction exiting boundary: 0.00360405
Surface-checks/particle/step: 0.167646
Surface-collisions/particle/step: 0.00107547
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 14442.5 ave 18396 max 9423 min
Histogram: 1 0 0 1 0 0 0 0 0 2
Cells:      25 ave 25 max 25 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 11 ave 11 max 11 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 11 ave 11 max 11 min
Histogram: 4 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0
