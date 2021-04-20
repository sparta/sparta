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
Created 100 child grid cells
  parent cells = 1
  CPU time = 0.000817 secs
  create/ghost percent = 91.9217 8.07834
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000173 secs
  reassign/sort/migrate/ghost percent = 79.7688 0.578035 8.67052 10.9827

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
  CPU time = 0.000263 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 44.1065 14.4487 0.760456 28.8973 11.7871 7.22433 0
  surf2grid time = 7.6e-05 secs
  map/rvous1/rvous2/split percent = 25 40.7895 0 21.0526
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
     100     0.022854    20898        0        0       37     3332 
     200     0.081191    35935        0        0       44     6622 
     300     0.153248    43670        0        0       61     8092 
     400       0.2716    47842        0        0       64     8456 
     500     0.392426    50468        0        0       63     8960 
Loop time of 0.392436 on 1 procs for 500 steps with 50468 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.26243    | 0.26243    | 0.26243    |   0.0 | 66.87
Coll    | 0.058244   | 0.058244   | 0.058244   |   0.0 | 14.84
Sort    | 0.046634   | 0.046634   | 0.046634   |   0.0 | 11.88
Comm    | 0.00111    | 0.00111    | 0.00111    |   0.0 |  0.28
Modify  | 0.02379    | 0.02379    | 0.02379    |   0.0 |  6.06
Output  | 4.9e-05    | 4.9e-05    | 4.9e-05    |   0.0 |  0.01
Other   |            | 0.000178   |            |       |  0.05

Particle moves    = 17610723 (17.6M)
Cells touched     = 18820431 (18.8M)
Particle comms    = 0 (0K)
Boundary collides = 59114 (59.1K)
Boundary exits    = 54845 (54.8K)
SurfColl checks   = 2998660 (3M)
SurfColl occurs   = 21600 (21.6K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 4.48754e+07
Particle-moves/step: 35221.4
Cell-touches/particle/step: 1.06869
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.0033567
Particle fraction exiting boundary: 0.0031143
Surface-checks/particle/step: 0.170275
Surface-collisions/particle/step: 0.00122653
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 50468 ave 50468 max 50468 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      100 ave 100 max 100 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0

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
  2031 deleted particles
  CPU time = 0.000391 secs
  sort/surf2grid/ghost/inout/particle percent = 37.5959 25.5754 4.85934 3.58056 28.3887
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26903 8.26903 8.26903
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48437        0        0        0        0 
     600     0.094281    50401        0        0       61     9534 
     700     0.194283    51491        0        0       66     9870 
     800     0.299254    52100        0        0       52     9590 
     900     0.405242    52633        0        0       63    10276 
    1000     0.504421    53228        0        0       55     9926 
Loop time of 0.504431 on 1 procs for 500 steps with 53228 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.34166    | 0.34166    | 0.34166    |   0.0 | 67.73
Coll    | 0.089031   | 0.089031   | 0.089031   |   0.0 | 17.65
Sort    | 0.052135   | 0.052135   | 0.052135   |   0.0 | 10.34
Comm    | 0.001742   | 0.001742   | 0.001742   |   0.0 |  0.35
Modify  | 0.019627   | 0.019627   | 0.019627   |   0.0 |  3.89
Output  | 5.7e-05    | 5.7e-05    | 5.7e-05    |   0.0 |  0.01
Other   |            | 0.00018    |            |       |  0.04

Particle moves    = 25855761 (25.9M)
Cells touched     = 27483462 (27.5M)
Particle comms    = 0 (0K)
Boundary collides = 85464 (85.5K)
Boundary exits    = 100485 (0.1M)
SurfColl checks   = 4755268 (4.76M)
SurfColl occurs   = 30843 (30.8K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 5.12573e+07
Particle-moves/step: 51711.5
Cell-touches/particle/step: 1.06295
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00330541
Particle fraction exiting boundary: 0.00388637
Surface-checks/particle/step: 0.183915
Surface-collisions/particle/step: 0.00119289
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 53228 ave 53228 max 53228 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      100 ave 100 max 100 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0

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
  1448 deleted particles
  CPU time = 0.000378 secs
  sort/surf2grid/ghost/inout/particle percent = 36.7725 24.6032 5.02646 3.7037 29.8942
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26903 8.26903 8.26903
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51780        0        0        0        0 
    1100     0.104299    52330        0        0       47     9002 
    1200     0.215713    53105        0        0       67    10192 
    1300     0.318203    53850        0        0       56     9310 
    1400     0.415243    53970        0        0       65     9814 
    1500     0.513817    54135        0        0       64    10094 
Loop time of 0.513825 on 1 procs for 500 steps with 54135 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.34716    | 0.34716    | 0.34716    |   0.0 | 67.56
Coll    | 0.092286   | 0.092286   | 0.092286   |   0.0 | 17.96
Sort    | 0.052887   | 0.052887   | 0.052887   |   0.0 | 10.29
Comm    | 0.001633   | 0.001633   | 0.001633   |   0.0 |  0.32
Modify  | 0.019619   | 0.019619   | 0.019619   |   0.0 |  3.82
Output  | 6.9e-05    | 6.9e-05    | 6.9e-05    |   0.0 |  0.01
Other   |            | 0.000172   |            |       |  0.03

Particle moves    = 26719208 (26.7M)
Cells touched     = 28364396 (28.4M)
Particle comms    = 0 (0K)
Boundary collides = 87100 (87.1K)
Boundary exits    = 103032 (0.103M)
SurfColl checks   = 4746966 (4.75M)
SurfColl occurs   = 30703 (30.7K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 5.20006e+07
Particle-moves/step: 53438.4
Cell-touches/particle/step: 1.06157
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00325983
Particle fraction exiting boundary: 0.0038561
Surface-checks/particle/step: 0.177661
Surface-collisions/particle/step: 0.0011491
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 54135 ave 54135 max 54135 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      100 ave 100 max 100 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0

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
  700 deleted particles
  CPU time = 0.000353 secs
  sort/surf2grid/ghost/inout/particle percent = 37.9603 24.3626 5.66572 3.96601 28.0453
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26903 8.26903 8.26903
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53435        0        0        0        0 
    1600     0.098288    54435        0        0       56     9408 
    1700     0.217675    55805        0        0       59     9534 
    1800     0.355773    56500        0        0       73     9772 
    1900     0.468715    56995        0        0       51     9492 
    2000     0.609474    57205        0        0       65     9940 
Loop time of 0.609485 on 1 procs for 500 steps with 57205 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.40035    | 0.40035    | 0.40035    |   0.0 | 65.69
Coll    | 0.10776    | 0.10776    | 0.10776    |   0.0 | 17.68
Sort    | 0.078401   | 0.078401   | 0.078401   |   0.0 | 12.86
Comm    | 0.001951   | 0.001951   | 0.001951   |   0.0 |  0.32
Modify  | 0.020677   | 0.020677   | 0.020677   |   0.0 |  3.39
Output  | 9.5e-05    | 9.5e-05    | 9.5e-05    |   0.0 |  0.02
Other   |            | 0.000252   |            |       |  0.04

Particle moves    = 27997286 (28M)
Cells touched     = 29716467 (29.7M)
Particle comms    = 0 (0K)
Boundary collides = 92246 (92.2K)
Boundary exits    = 101524 (0.102M)
SurfColl checks   = 4645466 (4.65M)
SurfColl occurs   = 30158 (30.2K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 4.5936e+07
Particle-moves/step: 55994.6
Cell-touches/particle/step: 1.06141
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00329482
Particle fraction exiting boundary: 0.00362621
Surface-checks/particle/step: 0.165926
Surface-collisions/particle/step: 0.00107718
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 57205 ave 57205 max 57205 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      100 ave 100 max 100 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0

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
  1435 deleted particles
  CPU time = 0.000503 secs
  sort/surf2grid/ghost/inout/particle percent = 44.1352 21.67 4.17495 2.98211 27.0378
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26903 8.26903 8.26903
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55770        0        0        0        0 
    2100     0.119662    56429        0        0       65     9842 
    2200     0.230169    57127        0        0       57     9828 
    2300     0.344181    57616        0        0       59     9268 
    2400     0.449076    57978        0        0       61     9688 
    2500      0.55796    57908        0        0       74     9786 
Loop time of 0.557969 on 1 procs for 500 steps with 57908 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.37213    | 0.37213    | 0.37213    |   0.0 | 66.69
Coll    | 0.10509    | 0.10509    | 0.10509    |   0.0 | 18.83
Sort    | 0.059489   | 0.059489   | 0.059489   |   0.0 | 10.66
Comm    | 0.001761   | 0.001761   | 0.001761   |   0.0 |  0.32
Modify  | 0.01929    | 0.01929    | 0.01929    |   0.0 |  3.46
Output  | 5.4e-05    | 5.4e-05    | 5.4e-05    |   0.0 |  0.01
Other   |            | 0.000155   |            |       |  0.03

Particle moves    = 28696403 (28.7M)
Cells touched     = 30448988 (30.4M)
Particle comms    = 0 (0K)
Boundary collides = 95250 (95.2K)
Boundary exits    = 103205 (0.103M)
SurfColl checks   = 4779250 (4.78M)
SurfColl occurs   = 30730 (30.7K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 5.14301e+07
Particle-moves/step: 57392.8
Cell-touches/particle/step: 1.06107
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00331923
Particle fraction exiting boundary: 0.00359644
Surface-checks/particle/step: 0.166545
Surface-collisions/particle/step: 0.00107087
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 57908 ave 57908 max 57908 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      100 ave 100 max 100 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
