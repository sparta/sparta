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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/runner/work/sparta/sparta/src/grid.cpp:465)
Created 100 child grid cells
  CPU time = 0.00302214 secs
  create/ghost percent = 92.075 7.92495
balance_grid        rcb cell
Balance grid migrated 74 cells
  CPU time = 0.00140622 secs
  reassign/sort/migrate/ghost percent = 59.579 1.03113 14.6565 24.7333

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
  CPU time = 0.00213523 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 23.5809 11.7741 3.04886 47.1759 14.4202 17.6892 1.00224
  surf2grid time = 0.00100731 secs
  map/comm1/comm2/comm3/comm4/split percent = 30.2591 15.5564 9.01418 6.85993 15.894 12.1116
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
     100  0.056261918    20901        0        0       14     3122 
     200   0.13155518    35945        0        0       53     6650 
     300   0.21714487    43816        0        0       47     7812 
     400   0.31088587    47999        0        0       65     8414 
     500   0.41077684    50612        0        0       62     9296 
Loop time of 0.410658 on 4 procs for 500 steps with 50612 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.036625   | 0.0662     | 0.10938    |  11.5 | 16.12
Coll    | 0.0046729  | 0.012774   | 0.019959   |   5.9 |  3.11
Sort    | 0.0079836  | 0.017169   | 0.024932   |   5.9 |  4.18
Comm    | 0.076314   | 0.14228    | 0.2228     |  13.9 | 34.65
Modify  | 0.0001177  | 0.0077252  | 0.015643   |   8.6 |  1.88
Output  | 0.00055061 | 0.0013546  | 0.0019861  |   1.5 |  0.33
Other   |            | 0.1632     |            |       | 39.73

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

Particle-moves/CPUsec/proc: 1.0734e+07
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
  CPU time = 0.00123521 secs
  sort/surf2grid/ghost/inout/particle percent = 13.7306 48.7451 11.9819 18.2643 7.27817
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.10938 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 3.62832 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48710        0        0        0        0 
     600  0.095425118    50401        0        0       71     9800 
     700   0.19356867    51556        0        0       68     8722 
     800   0.29262913    52357        0        0       65     9576 
     900    0.3920153    52699        0        0       53     9072 
    1000   0.49157887    53033        0        0       65     9282 
Loop time of 0.49174 on 4 procs for 500 steps with 53033 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.058747   | 0.096508   | 0.15325    |  11.7 | 19.63
Coll    | 0.011323   | 0.020476   | 0.028031   |   4.9 |  4.16
Sort    | 0.017049   | 0.025196   | 0.032414   |   4.2 |  5.12
Comm    | 0.078992   | 0.14545    | 0.25931    |  17.9 | 29.58
Modify  | 0.0001396  | 0.0071301  | 0.015033   |   8.3 |  1.45
Output  | 0.0001814  | 0.0014978  | 0.0025967  |   2.6 |  0.30
Other   |            | 0.1955     |            |       | 39.75

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

Particle-moves/CPUsec/proc: 1.31668e+07
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
  CPU time = 0.00104631 secs
  sort/surf2grid/ghost/inout/particle percent = 15.3111 41.8618 13.1129 21.4948 8.21943
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51589        0        0        0        0 
    1100  0.095405018    52493        0        0       48     8820 
    1200    0.1963878    53385        0        0       60     9618 
    1300   0.30058343    53890        0        0       70     9856 
    1400   0.40325384    54266        0        0       63    10094 
    1500   0.50547924    54428        0        0       64    10164 
Loop time of 0.505637 on 4 procs for 500 steps with 54428 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.064759   | 0.099815   | 0.15598    |  11.5 | 19.74
Coll    | 0.013017   | 0.021793   | 0.030032   |   4.8 |  4.31
Sort    | 0.018279   | 0.026195   | 0.034395   |   4.3 |  5.18
Comm    | 0.07611    | 0.1485     | 0.27023    |  18.9 | 29.37
Modify  | 0.0001463  | 0.0066505  | 0.013389   |   8.0 |  1.32
Output  | 0.0001775  | 0.0015355  | 0.0026421  |   2.5 |  0.30
Other   |            | 0.2011     |            |       | 39.78

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

Particle-moves/CPUsec/proc: 1.32464e+07
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
  CPU time = 0.00101021 secs
  sort/surf2grid/ghost/inout/particle percent = 15.7395 42.6648 13.5716 21.4908 6.53337
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53732        0        0        0        0 
    1600  0.095627217    54687        0        0       63     9464 
    1700   0.19774522    55802        0        0       72     9870 
    1800   0.30127003    56474        0        0       59     9408 
    1900   0.40617057    56965        0        0       57     9730 
    2000   0.51140261    57244        0        0       58     9618 
Loop time of 0.511583 on 4 procs for 500 steps with 57244 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.073804   | 0.10278    | 0.15239    |  10.1 | 20.09
Coll    | 0.012273   | 0.023171   | 0.031841   |   5.7 |  4.53
Sort    | 0.017974   | 0.027636   | 0.03596    |   4.7 |  5.40
Comm    | 0.073079   | 0.14819    | 0.27162    |  19.3 | 28.97
Modify  | 0.0001375  | 0.0066236  | 0.013247   |   8.0 |  1.29
Output  | 0.0001784  | 0.0015402  | 0.0027073  |   2.4 |  0.30
Other   |            | 0.2016     |            |       | 39.41

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

Particle-moves/CPUsec/proc: 1.37025e+07
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
  CPU time = 0.00102301 secs
  sort/surf2grid/ghost/inout/particle percent = 16.6862 42.4926 13.3725 19.5015 7.94721
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55781        0        0        0        0 
    2100    0.1018877    56408        0        0       59     9408 
    2200   0.20873046    57315        0        0       61     9212 
    2300   0.31561882    57640        0        0       75     9954 
    2400   0.42345039    57767        0        0       70    10290 
    2500   0.53056116    57770        0        0       69    10584 
Loop time of 0.530741 on 4 procs for 500 steps with 57770 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.072226   | 0.10518    | 0.15926    |  10.7 | 19.82
Coll    | 0.011727   | 0.024017   | 0.033396   |   5.9 |  4.53
Sort    | 0.017584   | 0.028297   | 0.036377   |   4.9 |  5.33
Comm    | 0.084063   | 0.15602    | 0.27726    |  18.4 | 29.40
Modify  | 0.0001393  | 0.0066679  | 0.013359   |   8.0 |  1.26
Output  | 0.0001814  | 0.0016744  | 0.0029329  |   2.6 |  0.32
Other   |            | 0.2089     |            |       | 39.36

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

Particle-moves/CPUsec/proc: 1.3511e+07
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
