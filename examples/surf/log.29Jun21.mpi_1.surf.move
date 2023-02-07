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
Created 100 child grid cells
  CPU time = 0.000945091 secs
  create/ghost percent = 84.561 15.439
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.00015974 secs
  reassign/sort/migrate/ghost percent = 57.4627 1.64179 20 20.8955

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
  CPU time = 0.000664949 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 30.4769 11.151 1.0398 46.5041 10.8283 5.88024 0.0358551
  surf2grid time = 0.000309229 secs
  map/comm1/comm2/comm3/comm4/split percent = 41.6345 12.4133 8.40401 3.93215 16.4225 13.1843
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
     100  0.055136919    20918        0        0       31     3080 
     200   0.19117332    36005        0        0       47     6622 
     300   0.37282896    43617        0        0       62     7700 
     400   0.58409691    48013        0        0       71     8610 
     500   0.81000447    50729        0        0       55     8456 
Loop time of 0.810027 on 1 procs for 500 steps with 50729 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.55683    | 0.55683    | 0.55683    |   0.0 | 68.74
Coll    | 0.084132   | 0.084132   | 0.084132   |   0.0 | 10.39
Sort    | 0.098967   | 0.098967   | 0.098967   |   0.0 | 12.22
Comm    | 0.0025785  | 0.0025785  | 0.0025785  |   0.0 |  0.32
Modify  | 0.066316   | 0.066316   | 0.066316   |   0.0 |  8.19
Output  | 0.00014734 | 0.00014734 | 0.00014734 |   0.0 |  0.02
Other   |            | 0.001052   |            |       |  0.13

Particle moves    = 17640116 (17.6M)
Cells touched     = 18850705 (18.9M)
Particle comms    = 0 (0K)
Boundary collides = 59466 (59.5K)
Boundary exits    = 54593 (54.6K)
SurfColl checks   = 3017210 (3.02M)
SurfColl occurs   = 21530 (21.5K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.17772e+07
Particle-moves/step: 35280.2
Cell-touches/particle/step: 1.06863
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00337107
Particle fraction exiting boundary: 0.00309482
Surface-checks/particle/step: 0.171043
Surface-collisions/particle/step: 0.00122051
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 50729 ave 50729 max 50729 min
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
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1881 deleted particles
  CPU time = 0.000876427 secs
  sort/surf2grid/ghost/inout/particle percent = 37.7312 25.0544 5.4951 4.27095 27.4483
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48848        0        0        0        0 
     600   0.22455573    50538        0        0       67     9912 
     700   0.45578265    51573        0        0       70     9898 
     800   0.69130349    52461        0        0       70    10136 
     900   0.93599987    53020        0        0       60     9492 
    1000     1.176517    53326        0        0       63     9590 
Loop time of 1.17653 on 1 procs for 500 steps with 53326 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.82094    | 0.82094    | 0.82094    |   0.0 | 69.78
Coll    | 0.14332    | 0.14332    | 0.14332    |   0.0 | 12.18
Sort    | 0.14707    | 0.14707    | 0.14707    |   0.0 | 12.50
Comm    | 0.0041859  | 0.0041859  | 0.0041859  |   0.0 |  0.36
Modify  | 0.059793   | 0.059793   | 0.059793   |   0.0 |  5.08
Output  | 0.00011611 | 0.00011611 | 0.00011611 |   0.0 |  0.01
Other   |            | 0.00111    |            |       |  0.09

Particle moves    = 25976617 (26M)
Cells touched     = 27611294 (27.6M)
Particle comms    = 0 (0K)
Boundary collides = 86022 (86K)
Boundary exits    = 100890 (0.101M)
SurfColl checks   = 4739308 (4.74M)
SurfColl occurs   = 30852 (30.9K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.20789e+07
Particle-moves/step: 51953.2
Cell-touches/particle/step: 1.06293
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00331152
Particle fraction exiting boundary: 0.00388388
Surface-checks/particle/step: 0.182445
Surface-collisions/particle/step: 0.00118768
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 53326 ave 53326 max 53326 min
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
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1405 deleted particles
  CPU time = 0.000961304 secs
  sort/surf2grid/ghost/inout/particle percent = 37.252 26.1409 4.21627 4.43948 27.9514
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51921        0        0        0        0 
    1100   0.23346949    52243        0        0       58     9786 
    1200    0.4715116    53130        0        0       65     9632 
    1300   0.71859479    53567        0        0       76    10150 
    1400   0.96152353    53940        0        0       78    10052 
    1500    1.2068017    54231        0        0       69     9758 
Loop time of 1.20682 on 1 procs for 500 steps with 54231 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.83963    | 0.83963    | 0.83963    |   0.0 | 69.57
Coll    | 0.14967    | 0.14967    | 0.14967    |   0.0 | 12.40
Sort    | 0.15186    | 0.15186    | 0.15186    |   0.0 | 12.58
Comm    | 0.0043628  | 0.0043628  | 0.0043628  |   0.0 |  0.36
Modify  | 0.060035   | 0.060035   | 0.060035   |   0.0 |  4.97
Output  | 0.00012517 | 0.00012517 | 0.00012517 |   0.0 |  0.01
Other   |            | 0.001137   |            |       |  0.09

Particle moves    = 26705888 (26.7M)
Cells touched     = 28349568 (28.3M)
Particle comms    = 0 (0K)
Boundary collides = 87110 (87.1K)
Boundary exits    = 102973 (0.103M)
SurfColl checks   = 4705162 (4.71M)
SurfColl occurs   = 30430 (30.4K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.21291e+07
Particle-moves/step: 53411.8
Cell-touches/particle/step: 1.06155
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00326183
Particle fraction exiting boundary: 0.00385582
Surface-checks/particle/step: 0.176184
Surface-collisions/particle/step: 0.00113945
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 54231 ave 54231 max 54231 min
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
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  786 deleted particles
  CPU time = 0.000905037 secs
  sort/surf2grid/ghost/inout/particle percent = 41.5964 25.5005 4.42571 4.03056 24.4468
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53445        0        0        0        0 
    1600   0.24068642    54383        0        0       59     9142 
    1700   0.49522042    55542        0        0       54     8806 
    1800   0.74906659    56407        0        0       57     9436 
    1900    1.0060201    56697        0        0       59    10150 
    2000    1.2658842    57299        0        0       65    10472 
Loop time of 1.2659 on 1 procs for 500 steps with 57299 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.87365    | 0.87365    | 0.87365    |   0.0 | 69.01
Coll    | 0.16406    | 0.16406    | 0.16406    |   0.0 | 12.96
Sort    | 0.1616     | 0.1616     | 0.1616     |   0.0 | 12.77
Comm    | 0.0044456  | 0.0044456  | 0.0044456  |   0.0 |  0.35
Modify  | 0.060802   | 0.060802   | 0.060802   |   0.0 |  4.80
Output  | 0.00013161 | 0.00013161 | 0.00013161 |   0.0 |  0.01
Other   |            | 0.001221   |            |       |  0.10

Particle moves    = 27939283 (27.9M)
Cells touched     = 29655039 (29.7M)
Particle comms    = 0 (0K)
Boundary collides = 91806 (91.8K)
Boundary exits    = 101497 (0.101M)
SurfColl checks   = 4651612 (4.65M)
SurfColl occurs   = 29739 (29.7K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.20706e+07
Particle-moves/step: 55878.6
Cell-touches/particle/step: 1.06141
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00328591
Particle fraction exiting boundary: 0.00363277
Surface-checks/particle/step: 0.16649
Surface-collisions/particle/step: 0.00106442
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 57299 ave 57299 max 57299 min
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
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1505 deleted particles
  CPU time = 0.00101471 secs
  sort/surf2grid/ghost/inout/particle percent = 42.7397 22.11 3.97086 4.3938 26.7857
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55794        0        0        0        0 
    2100   0.25925422    56563        0        0       60     9464 
    2200   0.51816368    57262        0        0       52     9464 
    2300   0.78029251    57706        0        0       62     9394 
    2400    1.0440998    57892        0        0       69    10556 
    2500     1.317132    58073        0        0       64    10150 
Loop time of 1.31715 on 1 procs for 500 steps with 58073 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.90008    | 0.90008    | 0.90008    |   0.0 | 68.33
Coll    | 0.17833    | 0.17833    | 0.17833    |   0.0 | 13.54
Sort    | 0.1708     | 0.1708     | 0.1708     |   0.0 | 12.97
Comm    | 0.0045698  | 0.0045698  | 0.0045698  |   0.0 |  0.35
Modify  | 0.061935   | 0.061935   | 0.061935   |   0.0 |  4.70
Output  | 0.00013137 | 0.00013137 | 0.00013137 |   0.0 |  0.01
Other   |            | 0.001315   |            |       |  0.10

Particle moves    = 28728272 (28.7M)
Cells touched     = 30479602 (30.5M)
Particle comms    = 0 (0K)
Boundary collides = 94930 (94.9K)
Boundary exits    = 103154 (0.103M)
SurfColl checks   = 4825366 (4.83M)
SurfColl occurs   = 30946 (30.9K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.18109e+07
Particle-moves/step: 57456.5
Cell-touches/particle/step: 1.06096
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00330441
Particle fraction exiting boundary: 0.00359068
Surface-checks/particle/step: 0.167966
Surface-collisions/particle/step: 0.0010772
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 58073 ave 58073 max 58073 min
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
