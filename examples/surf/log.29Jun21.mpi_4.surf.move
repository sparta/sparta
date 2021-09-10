SPARTA (26 Feb 2021)
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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/ascldap/users/stamoor/sparta_stanmoore1/src/grid.cpp:410)
Created 100 child grid cells
  CPU time = 0.00108981 secs
  create/ghost percent = 84.2485 15.7515
balance_grid        rcb cell
Balance grid migrated 74 cells
  CPU time = 0.000790358 secs
  reassign/sort/migrate/ghost percent = 69.7436 0.784314 12.006 17.4661

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
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  CPU time = 0.00104165 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 25.4063 17.7157 0.686656 46.6011 9.5903 16.8917 0.320439
  surf2grid time = 0.00048542 secs
  map/comm1/comm2/comm3/comm4/split percent = 26.6699 10.609 5.64833 4.56778 10.0688 32.9568
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
     100  0.033916712    20901        0        0       14     3122 
     200    0.1001997    35945        0        0       53     6650 
     300   0.19177508    43816        0        0       47     7812 
     400   0.28949904    47999        0        0       65     8414 
     500   0.38408399    50612        0        0       62     9296 
Loop time of 0.384179 on 4 procs for 500 steps with 50612 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.079534   | 0.14961    | 0.23665    |  16.8 | 38.94
Coll    | 0.006043   | 0.015653   | 0.024828   |   6.4 |  4.07
Sort    | 0.018365   | 0.032112   | 0.043346   |   6.0 |  8.36
Comm    | 0.016084   | 0.017218   | 0.018188   |   0.6 |  4.48
Modify  | 0.00022531 | 0.017527   | 0.035843   |  13.1 |  4.56
Output  | 0.00013566 | 0.0027124  | 0.0095782  |   7.6 |  0.71
Other   |            | 0.1494     |            |       | 38.88

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

Particle-moves/CPUsec/proc: 1.14738e+07
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
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1902 deleted particles
  CPU time = 0.00100994 secs
  sort/surf2grid/ghost/inout/particle percent = 14.542 47.4032 6.20869 19.0982 12.7479
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.10938 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 3.62832 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48710        0        0        0        0 
     600  0.092022419    50401        0        0       71     9800 
     700   0.18600965    51556        0        0       68     8722 
     800   0.28100276    52357        0        0       65     9576 
     900   0.37642574    52699        0        0       53     9072 
    1000   0.47471309    53033        0        0       65     9282 
Loop time of 0.474794 on 4 procs for 500 steps with 53033 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.13591    | 0.21565    | 0.32555    |  15.8 | 45.42
Coll    | 0.014035   | 0.024503   | 0.032903   |   5.1 |  5.16
Sort    | 0.036485   | 0.048648   | 0.058772   |   4.2 | 10.25
Comm    | 0.016048   | 0.017452   | 0.018577   |   0.7 |  3.68
Modify  | 0.00024796 | 0.017366   | 0.035022   |  13.0 |  3.66
Output  | 0.00017667 | 0.00036263 | 0.00053191 |   0.0 |  0.08
Other   |            | 0.1508     |            |       | 31.76

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

Particle-moves/CPUsec/proc: 1.36367e+07
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
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1444 deleted particles
  CPU time = 0.000891924 secs
  sort/surf2grid/ghost/inout/particle percent = 15.8514 45.576 7.35098 19.0591 12.1625
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51589        0        0        0        0 
    1100  0.091908693    52493        0        0       48     8820 
    1200   0.18849063    53385        0        0       60     9618 
    1300   0.28580928    53890        0        0       70     9856 
    1400   0.40406227    54266        0        0       63    10094 
    1500   0.50223327    54428        0        0       64    10164 
Loop time of 0.502306 on 4 procs for 500 steps with 54428 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.14545    | 0.22444    | 0.33103    |  15.8 | 44.68
Coll    | 0.016014   | 0.025905   | 0.034829   |   5.1 |  5.16
Sort    | 0.040361   | 0.051863   | 0.063148   |   4.2 | 10.33
Comm    | 0.016409   | 0.017987   | 0.019258   |   0.8 |  3.58
Modify  | 0.00024557 | 0.017172   | 0.034889   |  12.9 |  3.42
Output  | 0.00019526 | 0.00039089 | 0.00055981 |   0.0 |  0.08
Other   |            | 0.1645     |            |       | 32.76

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

Particle-moves/CPUsec/proc: 1.33343e+07
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
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  696 deleted particles
  CPU time = 0.00083828 secs
  sort/surf2grid/ghost/inout/particle percent = 19.2264 43.4585 8.24801 19.8805 9.18658
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53732        0        0        0        0 
    1600  0.093466043    54687        0        0       63     9464 
    1700   0.19234276    55802        0        0       72     9870 
    1800   0.29269004    56474        0        0       59     9408 
    1900   0.39409947    56965        0        0       57     9730 
    2000   0.49596977    57244        0        0       58     9618 
Loop time of 0.496048 on 4 procs for 500 steps with 57244 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.16372    | 0.23139    | 0.33236    |  14.4 | 46.65
Coll    | 0.015436   | 0.027496   | 0.037112   |   5.8 |  5.54
Sort    | 0.03859    | 0.053943   | 0.065894   |   5.1 | 10.87
Comm    | 0.016884   | 0.018542   | 0.019869   |   0.8 |  3.74
Modify  | 0.00024867 | 0.017401   | 0.034729   |  13.0 |  3.51
Output  | 0.00013852 | 0.00041157 | 0.00064969 |   0.0 |  0.08
Other   |            | 0.1469     |            |       | 29.61

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

Particle-moves/CPUsec/proc: 1.41316e+07
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
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1463 deleted particles
  CPU time = 0.000971317 secs
  sort/surf2grid/ghost/inout/particle percent = 17.2067 35.9352 17.8694 17.673 11.3157
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05019 3.20644 4.89394
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55781        0        0        0        0 
    2100   0.10916519    56408        0        0       59     9408 
    2200   0.21731043    57315        0        0       61     9212 
    2300   0.32172465    57640        0        0       75     9954 
    2400   0.44650126    57767        0        0       70    10290 
    2500   0.55125737    57770        0        0       69    10584 
Loop time of 0.551374 on 4 procs for 500 steps with 57770 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.16321    | 0.24333    | 0.35092    |  15.5 | 44.13
Coll    | 0.015712   | 0.029068   | 0.039716   |   5.8 |  5.27
Sort    | 0.039325   | 0.057848   | 0.073704   |   5.8 | 10.49
Comm    | 0.018356   | 0.02199    | 0.023971   |   1.5 |  3.99
Modify  | 0.00026298 | 0.017499   | 0.035097   |  13.0 |  3.17
Output  | 0.00014925 | 0.0005247  | 0.00074387 |   0.0 |  0.10
Other   |            | 0.1811     |            |       | 32.85

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

Particle-moves/CPUsec/proc: 1.30053e+07
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
