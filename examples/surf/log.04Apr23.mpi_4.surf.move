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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/ascldap/users/stamoor/sparta_master/src/grid.cpp:465)
Created 100 child grid cells
  CPU time = 0.00112033 secs
  create/ghost percent = 85.1245 14.8755
balance_grid        rcb cell
Balance grid migrated 74 cells
  CPU time = 0.000850439 secs
  reassign/sort/migrate/ghost percent = 66.4144 0.6448 14.2977 18.6431

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
  CPU time = 0.00126672 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 22.1532 20.6475 1.05402 42.6501 13.4952 16.262 0.376435
  surf2grid time = 0.000540257 secs
  map/comm1/comm2/comm3/comm4/split percent = 37.9965 10.6355 6.31068 4.98676 12.489 16.0194
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
     100  0.059676886    20901        0        0       14     3122 
     200   0.19077659    35945        0        0       53     6650 
     300   0.35146976    43816        0        0       47     7812 
     400   0.52632523    47999        0        0       65     8414 
     500   0.71021914    50612        0        0       62     9296 
Loop time of 0.710327 on 4 procs for 500 steps with 50612 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.1625     | 0.29243    | 0.51248    |  26.3 | 41.17
Coll    | 0.014779   | 0.03162    | 0.046373   |   7.6 |  4.45
Sort    | 0.025594   | 0.049407   | 0.069613   |   8.6 |  6.96
Comm    | 0.017982   | 0.0191     | 0.020195   |   0.6 |  2.69
Modify  | 0.00036526 | 0.027157   | 0.054059   |  16.3 |  3.82
Output  | 0.00015068 | 0.00057125 | 0.0009706  |   0.0 |  0.08
Other   |            | 0.29       |            |       | 40.83

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

Particle-moves/CPUsec/proc: 6.20561e+06
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
  CPU time = 0.00142741 secs
  sort/surf2grid/ghost/inout/particle percent = 14.5148 50.5762 7.06531 13.4625 14.3812
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.10938 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 3.62832 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48710        0        0        0        0 
     600   0.19598222    50401        0        0       71     9800 
     700   0.38956022    51556        0        0       68     8722 
     800   0.58523059    52357        0        0       65     9576 
     900   0.78132558    52699        0        0       53     9072 
    1000    0.9796176    53033        0        0       65     9282 
Loop time of 0.979702 on 4 procs for 500 steps with 53033 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.24716    | 0.42972    | 0.74308    |  29.4 | 43.86
Coll    | 0.031273   | 0.047664   | 0.060627   |   5.7 |  4.87
Sort    | 0.051474   | 0.072601   | 0.089586   |   6.0 |  7.41
Comm    | 0.018194   | 0.019542   | 0.020745   |   0.7 |  1.99
Modify  | 0.00036764 | 0.029217   | 0.064157   |  17.0 |  2.98
Output  | 0.00018907 | 0.00049591 | 0.00081539 |   0.0 |  0.05
Other   |            | 0.3805     |            |       | 38.83

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

Particle-moves/CPUsec/proc: 6.60877e+06
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
  CPU time = 0.00125527 secs
  sort/surf2grid/ghost/inout/particle percent = 16.3913 47.1035 7.97721 14.9858 13.5423
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51589        0        0        0        0 
    1100   0.18582368    52493        0        0       48     8820 
    1200   0.39226222    53385        0        0       60     9618 
    1300   0.59105468    53890        0        0       70     9856 
    1400   0.79144287    54266        0        0       63    10094 
    1500   0.99147177    54428        0        0       64    10164 
Loop time of 0.991544 on 4 procs for 500 steps with 54428 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.26489    | 0.44206    | 0.75433    |  29.1 | 44.58
Coll    | 0.034233   | 0.051995   | 0.072471   |   6.8 |  5.24
Sort    | 0.054975   | 0.074766   | 0.094041   |   6.1 |  7.54
Comm    | 0.018974   | 0.020499   | 0.021774   |   0.7 |  2.07
Modify  | 0.00036478 | 0.026386   | 0.052845   |  16.0 |  2.66
Output  | 0.00023055 | 0.00052136 | 0.00082231 |   0.0 |  0.05
Other   |            | 0.3753     |            |       | 37.85

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

Particle-moves/CPUsec/proc: 6.75501e+06
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
  CPU time = 0.0011065 secs
  sort/surf2grid/ghost/inout/particle percent = 19.6725 46.1108 8.5111 15.32 10.3857
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53732        0        0        0        0 
    1600   0.17933369    54687        0        0       63     9464 
    1700   0.38165236    55802        0        0       72     9870 
    1800   0.57777071    56474        0        0       59     9408 
    1900   0.77784991    56965        0        0       57     9730 
    2000    0.9770658    57244        0        0       58     9618 
Loop time of 0.977141 on 4 procs for 500 steps with 57244 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.30989    | 0.45935    | 0.71887    |  24.1 | 47.01
Coll    | 0.033301   | 0.052389   | 0.067323   |   6.5 |  5.36
Sort    | 0.054009   | 0.079325   | 0.099539   |   7.1 |  8.12
Comm    | 0.019865   | 0.021295   | 0.022505   |   0.7 |  2.18
Modify  | 0.00037098 | 0.026681   | 0.053818   |  16.1 |  2.73
Output  | 0.00015259 | 0.00055063 | 0.00094366 |   0.0 |  0.06
Other   |            | 0.3376     |            |       | 34.54

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

Particle-moves/CPUsec/proc: 7.17394e+06
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
  CPU time = 0.00119996 secs
  sort/surf2grid/ghost/inout/particle percent = 19.9483 42.9764 8.48401 14.8222 13.7691
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55781        0        0        0        0 
    2100   0.19332981    56408        0        0       59     9408 
    2200   0.40608144    57315        0        0       61     9212 
    2300   0.60979819    57640        0        0       75     9954 
    2400   0.81438184    57767        0        0       70    10290 
    2500    1.0184391    57770        0        0       69    10584 
Loop time of 1.01853 on 4 procs for 500 steps with 57770 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.3242     | 0.47058    | 0.75269    |  25.5 | 46.20
Coll    | 0.031925   | 0.053697   | 0.071348   |   7.1 |  5.27
Sort    | 0.05289    | 0.081497   | 0.10402    |   7.5 |  8.00
Comm    | 0.019227   | 0.020598   | 0.02177    |   0.7 |  2.02
Modify  | 0.00035834 | 0.026739   | 0.054458   |  16.1 |  2.63
Output  | 0.00014424 | 0.00060529 | 0.001055   |   0.0 |  0.06
Other   |            | 0.3648     |            |       | 35.82

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

Particle-moves/CPUsec/proc: 7.04035e+06
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
