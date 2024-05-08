SPARTA (13 Apr 2023)
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
  CPU time = 0.0018235 secs
  create/ghost percent = 92.5473 7.4527
balance_grid        rcb cell
Balance grid migrated 74 cells
  CPU time = 0.000755401 secs
  reassign/sort/migrate/ghost percent = 64.178 0.741328 13.794 21.2867

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
  CPU time = 0.001195 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 25.6318 13.9665 0.677824 45.9833 13.7406 16.5439 0.60251
  surf2grid time = 0.0005495 secs
  map/comm1/comm2/comm3/comm4/split percent = 35.3958 12.7753 8.20746 6.49682 17.2884 11.9745
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
     100  0.031113814    20901        0        0       14     3122 
     200  0.095664544    35945        0        0       53     6650 
     300   0.17620748    43816        0        0       47     7812 
     400   0.26416232    47999        0        0       65     8414 
     500   0.35944437    50612        0        0       62     9296 
Loop time of 0.359438 on 4 procs for 500 steps with 50612 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.036924   | 0.066977   | 0.1124     |  11.8 | 18.63
Coll    | 0.0052798  | 0.012141   | 0.018469   |   5.3 |  3.38
Sort    | 0.0070249  | 0.012867   | 0.017763   |   4.0 |  3.58
Comm    | 0.054628   | 0.14607    | 0.21322    |  17.9 | 40.64
Modify  | 0.0001002  | 0.0077599  | 0.016277   |   8.7 |  2.16
Output  | 0.0003949  | 0.00079    | 0.0011804  |   0.0 |  0.22
Other   |            | 0.1128     |            |       | 31.39

Particle moves    = 17577331 (17.6M)
Cells touched     = 18843630 (18.8M)
Particle comms    = 122668 (0.123M)
Boundary collides = 59486 (59.5K)
Boundary exits    = 54770 (54.8K)
SurfColl checks   = 3034052 (3.03M)
SurfColl occurs   = 21569 (21.6K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 1.22256e+07
Particle-moves/step: 35154.7
Cell-touches/particle/step: 1.07204
Particle comm iterations/step: 2.366
Particle fraction communicated: 0.00697876
Particle fraction colliding with boundary: 0.00338425
Particle fraction exiting boundary: 0.00311595
Surface-checks/particle/step: 0.172612
Surface-collisions/particle/step: 0.00122709
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
  CPU time = 0.0011118 secs
  sort/surf2grid/ghost/inout/particle percent = 10.2716 51.2052 12.9969 17.881 7.64525
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.10938 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 3.62832 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48710        0        0        0        0 
     600  0.092619843    50401        0        0       71     9800 
     700   0.18743329    51556        0        0       68     8722 
     800   0.28363393    52357        0        0       65     9576 
     900   0.38006238    52699        0        0       53     9072 
    1000   0.48047502    53033        0        0       65     9282 
Loop time of 0.480654 on 4 procs for 500 steps with 53033 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.059341   | 0.097475   | 0.15709    |  12.1 | 20.28
Coll    | 0.010743   | 0.018703   | 0.025999   |   4.9 |  3.89
Sort    | 0.014044   | 0.020256   | 0.024883   |   3.3 |  4.21
Comm    | 0.077301   | 0.19188    | 0.28843    |  19.7 | 39.92
Modify  | 0.0001224  | 0.007542   | 0.01596    |   8.6 |  1.57
Output  | 0.0005725  | 0.00085765 | 0.0011009  |   0.0 |  0.18
Other   |            | 0.1439     |            |       | 29.95

Particle moves    = 25797573 (25.8M)
Cells touched     = 27528826 (27.5M)
Particle comms    = 170023 (0.17M)
Boundary collides = 85994 (86K)
Boundary exits    = 101007 (0.101M)
SurfColl checks   = 4695488 (4.7M)
SurfColl occurs   = 30862 (30.9K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 1.34179e+07
Particle-moves/step: 51595.1
Cell-touches/particle/step: 1.06711
Particle comm iterations/step: 2.438
Particle fraction communicated: 0.00659066
Particle fraction colliding with boundary: 0.00333341
Particle fraction exiting boundary: 0.00391537
Surface-checks/particle/step: 0.182013
Surface-collisions/particle/step: 0.00119631
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
  CPU time = 0.0010067 secs
  sort/surf2grid/ghost/inout/particle percent = 11.7016 41.6709 14.0757 19.2013 13.3505
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51589        0        0        0        0 
    1100  0.092795943    52493        0        0       48     8820 
    1200   0.19256869    53385        0        0       60     9618 
    1300   0.29043583    53890        0        0       70     9856 
    1400   0.38892708    54266        0        0       63    10094 
    1500   0.48756733    54428        0        0       64    10164 
Loop time of 0.48769 on 4 procs for 500 steps with 54428 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.063278   | 0.10023    | 0.1605     |  12.1 | 20.55
Coll    | 0.011632   | 0.019801   | 0.028397   |   5.2 |  4.06
Sort    | 0.015096   | 0.020612   | 0.025701   |   3.1 |  4.23
Comm    | 0.081891   | 0.19007    | 0.27677    |  18.9 | 38.97
Modify  | 0.0001209  | 0.0070283  | 0.013993   |   8.2 |  1.44
Output  | 0.0008413  | 0.0009596  | 0.0010794  |   0.0 |  0.20
Other   |            | 0.149      |            |       | 30.55

Particle moves    = 26689096 (26.7M)
Cells touched     = 28436600 (28.4M)
Particle comms    = 171396 (0.171M)
Boundary collides = 87221 (87.2K)
Boundary exits    = 102500 (0.102M)
SurfColl checks   = 4773188 (4.77M)
SurfColl occurs   = 30765 (30.8K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 1.36814e+07
Particle-moves/step: 53378.2
Cell-touches/particle/step: 1.06548
Particle comm iterations/step: 2.73
Particle fraction communicated: 0.00642195
Particle fraction colliding with boundary: 0.00326804
Particle fraction exiting boundary: 0.00384052
Surface-checks/particle/step: 0.178844
Surface-collisions/particle/step: 0.00115272
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
  CPU time = 0.0009901 secs
  sort/surf2grid/ghost/inout/particle percent = 12.7361 39.9657 16.6448 21.21 9.44349
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53732        0        0        0        0 
    1600  0.093160943    54687        0        0       63     9464 
    1700   0.19192919    55802        0        0       72     9870 
    1800   0.29175204    56474        0        0       59     9408 
    1900   0.39309898    56965        0        0       57     9730 
    2000   0.49538403    57244        0        0       58     9618 
Loop time of 0.495467 on 4 procs for 500 steps with 57244 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.073806   | 0.10427    | 0.15653    |  10.5 | 21.05
Coll    | 0.011368   | 0.021747   | 0.031031   |   6.0 |  4.39
Sort    | 0.014747   | 0.021973   | 0.027683   |   3.8 |  4.43
Comm    | 0.09321    | 0.18659    | 0.27909    |  18.2 | 37.66
Modify  | 0.0001251  | 0.006955   | 0.014077   |   8.2 |  1.40
Output  | 0.0009881  | 0.0013232  | 0.0018069  |   0.9 |  0.27
Other   |            | 0.1526     |            |       | 30.80

Particle moves    = 27938061 (27.9M)
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

Particle-moves/CPUsec/proc: 1.40968e+07
Particle-moves/step: 55876.1
Cell-touches/particle/step: 1.06528
Particle comm iterations/step: 2.868
Particle fraction communicated: 0.00631254
Particle fraction colliding with boundary: 0.00329819
Particle fraction exiting boundary: 0.00364367
Surface-checks/particle/step: 0.166719
Surface-collisions/particle/step: 0.00107377
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
  CPU time = 0.000978401 secs
  sort/surf2grid/ghost/inout/particle percent = 13.5221 40.8115 14.1149 20.7379 10.8137
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55781        0        0        0        0 
    2100  0.099304046    56408        0        0       59     9408 
    2200   0.20348939    57315        0        0       61     9212 
    2300   0.30780544    57640        0        0       75     9954 
    2400   0.41224629    57767        0        0       70    10290 
    2500   0.51638774    57770        0        0       69    10584 
Loop time of 0.516447 on 4 procs for 500 steps with 57770 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.073004   | 0.10675    | 0.16496    |  11.3 | 20.67
Coll    | 0.010959   | 0.022504   | 0.032678   |   6.2 |  4.36
Sort    | 0.014423   | 0.022633   | 0.029193   |   4.0 |  4.38
Comm    | 0.087548   | 0.2004     | 0.30397    |  19.5 | 38.80
Modify  | 0.0001212  | 0.0069157  | 0.013847   |   8.2 |  1.34
Output  | 0.0007697  | 0.0012834  | 0.001788   |   1.1 |  0.25
Other   |            | 0.156      |            |       | 30.20

Particle moves    = 28579962 (28.6M)
Cells touched     = 30432584 (30.4M)
Particle comms    = 178617 (0.179M)
Boundary collides = 94726 (94.7K)
Boundary exits    = 103376 (0.103M)
SurfColl checks   = 4808608 (4.81M)
SurfColl occurs   = 30848 (30.8K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 1.38349e+07
Particle-moves/step: 57159.9
Cell-touches/particle/step: 1.06482
Particle comm iterations/step: 2.52
Particle fraction communicated: 0.00624973
Particle fraction colliding with boundary: 0.00331442
Particle fraction exiting boundary: 0.00361708
Surface-checks/particle/step: 0.168251
Surface-collisions/particle/step: 0.00107936
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
