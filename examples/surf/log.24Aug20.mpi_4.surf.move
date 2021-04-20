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
create_grid 	    10 10 1
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/Users/eharvey/dev/SPARTA.base/sparta/src/grid.cpp:415)
Created 100 child grid cells
  parent cells = 1
  CPU time = 0.000949 secs
  create/ghost percent = 96.1012 3.89884
balance_grid        rcb cell
Balance grid migrated 74 cells
  CPU time = 0.00067 secs
  reassign/sort/migrate/ghost percent = 83.2836 0.447761 7.01493 9.25373

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
  4 = cells with surfs
  56 = total surfs in all grid cells
  14 = max surfs in one grid cell
  0.124325 = min surf-size/cell-size ratio
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  CPU time = 0.000406 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 40.6404 11.8227 0.246305 34.4828 12.8079 11.8227 0.246305
  surf2grid time = 0.00014 secs
  map/rvous1/rvous2/split percent = 9.28571 64.2857 0.714286 15
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
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 1.51903 1.51903 1.51903
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
     100     0.016158    20845        0        0       30     2898 
     200     0.045336    35962        0        0       58     6916 
     300     0.081447    43638        0        0       61     7966 
     400     0.123138    47916        0        0       65     8680 
     500     0.168424    50599        0        0       60     9016 
Loop time of 0.168491 on 4 procs for 500 steps with 50599 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.035862   | 0.066212   | 0.10808    |  11.3 | 39.30
Coll    | 0.004316   | 0.010903   | 0.016622   |   5.1 |  6.47
Sort    | 0.005858   | 0.01105    | 0.015432   |   3.9 |  6.56
Comm    | 0.010892   | 0.011363   | 0.011855   |   0.3 |  6.74
Modify  | 7.3e-05    | 0.006606   | 0.013995   |   8.1 |  3.92
Output  | 0.00012    | 0.000246   | 0.000373   |   0.0 |  0.15
Other   |            | 0.06211    |            |       | 36.86

Particle moves    = 17599838 (17.6M)
Cells touched     = 18809534 (18.8M)
Particle comms    = 122801 (0.123M)
Boundary collides = 59274 (59.3K)
Boundary exits    = 54713 (54.7K)
SurfColl checks   = 3023328 (3.02M)
SurfColl occurs   = 21446 (21.4K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 2.61139e+07
Particle-moves/step: 35199.7
Cell-touches/particle/step: 1.06873
Particle comm iterations/step: 2.39
Particle fraction communicated: 0.00697739
Particle fraction colliding with boundary: 0.00336787
Particle fraction exiting boundary: 0.00310872
Surface-checks/particle/step: 0.171782
Surface-collisions/particle/step: 0.00121853
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 12649.8 ave 16269 max 7947 min
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
  4 = cells with surfs
  56 = total surfs in all grid cells
  14 = max surfs in one grid cell
  0.124325 = min surf-size/cell-size ratio
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1933 deleted particles
  CPU time = 0.000599 secs
  sort/surf2grid/ghost/inout/particle percent = 12.3539 41.0684 11.1853 14.0234 21.3689
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.10938 1.6875 3.375
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 3.6284 3.20653 4.89403
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48666        0        0        0        0 
     600     0.041403    50442        0        0       52     9254 
     700     0.079737    51613        0        0       55     9576 
     800     0.117274    52348        0        0       58    10024 
     900      0.15633    52886        0        0       59     9478 
    1000     0.194653    53190        0        0       59     9100 
Loop time of 0.194667 on 4 procs for 500 steps with 53190 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.054677   | 0.087711   | 0.13424    |  10.4 | 45.06
Coll    | 0.009579   | 0.01603    | 0.021603   |   4.0 |  8.23
Sort    | 0.011319   | 0.01514    | 0.018087   |   2.3 |  7.78
Comm    | 0.008235   | 0.0085155  | 0.008696   |   0.2 |  4.37
Modify  | 8.6e-05    | 0.0056547  | 0.011962   |   7.4 |  2.90
Output  | 5.5e-05    | 0.00015925 | 0.00026    |   0.0 |  0.08
Other   |            | 0.06146    |            |       | 31.57

Particle moves    = 25930347 (25.9M)
Cells touched     = 27560834 (27.6M)
Particle comms    = 170043 (0.17M)
Boundary collides = 85804 (85.8K)
Boundary exits    = 100823 (0.101M)
SurfColl checks   = 4705442 (4.71M)
SurfColl occurs   = 30529 (30.5K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 3.33009e+07
Particle-moves/step: 51860.7
Cell-touches/particle/step: 1.06288
Particle comm iterations/step: 2.466
Particle fraction communicated: 0.00655768
Particle fraction colliding with boundary: 0.00330902
Particle fraction exiting boundary: 0.00388822
Surface-checks/particle/step: 0.181465
Surface-collisions/particle/step: 0.00117735
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 13297.5 ave 16390 max 9693 min
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
  4 = cells with surfs
  56 = total surfs in all grid cells
  14 = max surfs in one grid cell
  0.124325 = min surf-size/cell-size ratio
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1464 deleted particles
  CPU time = 0.000287 secs
  sort/surf2grid/ghost/inout/particle percent = 12.5436 37.2822 15.331 16.0279 18.8153
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05028 3.20653 4.89403
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51726        0        0        0        0 
    1100     0.034864    52248        0        0       69    10262 
    1200     0.070937    52963        0        0       50     8890 
    1300     0.108583    53485        0        0       60     9870 
    1400     0.147393    53822        0        0       57     9450 
    1500     0.187917    54063        0        0       71     9548 
Loop time of 0.187955 on 4 procs for 500 steps with 54063 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.056579   | 0.086894   | 0.13035    |   9.8 | 46.23
Coll    | 0.010418   | 0.016218   | 0.021841   |   3.8 |  8.63
Sort    | 0.011284   | 0.014577   | 0.017914   |   2.4 |  7.76
Comm    | 0.008266   | 0.008578   | 0.008974   |   0.3 |  4.56
Modify  | 7.4e-05    | 0.0051273  | 0.010187   |   7.0 |  2.73
Output  | 6.2e-05    | 0.00013575 | 0.000209   |   0.0 |  0.07
Other   |            | 0.05643    |            |       | 30.02

Particle moves    = 26653037 (26.7M)
Cells touched     = 28299095 (28.3M)
Particle comms    = 171936 (0.172M)
Boundary collides = 86961 (87K)
Boundary exits    = 102968 (0.103M)
SurfColl checks   = 4736116 (4.74M)
SurfColl occurs   = 30432 (30.4K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 3.54514e+07
Particle-moves/step: 53306.1
Cell-touches/particle/step: 1.06176
Particle comm iterations/step: 2.756
Particle fraction communicated: 0.0064509
Particle fraction colliding with boundary: 0.00326271
Particle fraction exiting boundary: 0.00386327
Surface-checks/particle/step: 0.177695
Surface-collisions/particle/step: 0.00114178
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 13515.8 ave 16689 max 10176 min
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
  4 = cells with surfs
  56 = total surfs in all grid cells
  14 = max surfs in one grid cell
  0.124325 = min surf-size/cell-size ratio
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  702 deleted particles
  CPU time = 0.000455 secs
  sort/surf2grid/ghost/inout/particle percent = 12.7473 30.7692 14.7253 24.6154 17.1429
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05028 3.20653 4.89403
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53361        0        0        0        0 
    1600     0.036661    54322        0        0       45     9296 
    1700     0.074739    55673        0        0       53     9590 
    1800     0.112421    56285        0        0       52     9072 
    1900     0.152251    56792        0        0       64     9688 
    2000     0.191984    57285        0        0       63     9800 
Loop time of 0.192007 on 4 procs for 500 steps with 57285 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.064731   | 0.090089   | 0.13064    |   9.0 | 46.92
Coll    | 0.009901   | 0.017059   | 0.022826   |   4.4 |  8.88
Sort    | 0.010273   | 0.014657   | 0.018165   |   2.9 |  7.63
Comm    | 0.008374   | 0.0085037  | 0.008711   |   0.1 |  4.43
Modify  | 7.5e-05    | 0.0051053  | 0.010288   |   7.0 |  2.66
Output  | 4.5e-05    | 0.000144   | 0.000248   |   0.0 |  0.07
Other   |            | 0.05645    |            |       | 29.40

Particle moves    = 27939324 (27.9M)
Cells touched     = 29656592 (29.7M)
Particle comms    = 176024 (0.176M)
Boundary collides = 92087 (92.1K)
Boundary exits    = 101446 (0.101M)
SurfColl checks   = 4633258 (4.63M)
SurfColl occurs   = 29822 (29.8K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 3.6378e+07
Particle-moves/step: 55878.6
Cell-touches/particle/step: 1.06146
Particle comm iterations/step: 2.852
Particle fraction communicated: 0.00630022
Particle fraction colliding with boundary: 0.00329596
Particle fraction exiting boundary: 0.00363094
Surface-checks/particle/step: 0.165833
Surface-collisions/particle/step: 0.00106738
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 14321.2 ave 17849 max 9819 min
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
  4 = cells with surfs
  56 = total surfs in all grid cells
  14 = max surfs in one grid cell
  0.124325 = min surf-size/cell-size ratio
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1531 deleted particles
  CPU time = 0.00072 secs
  sort/surf2grid/ghost/inout/particle percent = 17.3611 32.2222 14.4444 18.3333 17.6389
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 2.53125 1.6875 3.375
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 4.05028 3.20653 4.89403
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55754        0        0        0        0 
    2100     0.041862    56560        0        0       61    10052 
    2200     0.085473    57024        0        0       57     9730 
    2300     0.125459    57626        0        0       67     9646 
    2400     0.165155    57768        0        0       67     9422 
    2500     0.203416    57972        0        0       54    10206 
Loop time of 0.203431 on 4 procs for 500 steps with 57972 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.062268   | 0.092447   | 0.13861    |   9.9 | 45.44
Coll    | 0.009312   | 0.017733   | 0.024243   |   4.7 |  8.72
Sort    | 0.010021   | 0.015437   | 0.019989   |   3.3 |  7.59
Comm    | 0.008255   | 0.0084203  | 0.008618   |   0.1 |  4.14
Modify  | 8e-05      | 0.0050843  | 0.010138   |   7.0 |  2.50
Output  | 3.8e-05    | 0.0001465  | 0.000267   |   0.0 |  0.07
Other   |            | 0.06416    |            |       | 31.54

Particle moves    = 28684528 (28.7M)
Cells touched     = 30434939 (30.4M)
Particle comms    = 178896 (0.179M)
Boundary collides = 94449 (94.4K)
Boundary exits    = 103092 (0.103M)
SurfColl checks   = 4782652 (4.78M)
SurfColl occurs   = 30787 (30.8K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 3.5251e+07
Particle-moves/step: 57369.1
Cell-touches/particle/step: 1.06102
Particle comm iterations/step: 2.446
Particle fraction communicated: 0.00623667
Particle fraction colliding with boundary: 0.00329268
Particle fraction exiting boundary: 0.00359399
Surface-checks/particle/step: 0.166733
Surface-collisions/particle/step: 0.0010733
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 14493 ave 18307 max 9593 min
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
