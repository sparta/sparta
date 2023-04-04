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
create_grid 	    10 10 1
Created 100 child grid cells
  CPU time = 0.00102854 secs
  create/ghost percent = 74.8957 25.1043
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000256538 secs
  reassign/sort/migrate/ghost percent = 40.5204 1.39405 21.0037 37.0818

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
  CPU time = 0.0010879 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 17.8611 13.5218 1.73132 47.5345 19.3513 5.06246 0.0219154
  surf2grid time = 0.00051713 secs
  map/comm1/comm2/comm3/comm4/split percent = 45.0899 11.4338 8.71369 3.82665 12.4481 13.9696
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
     100    0.1016295    20918        0        0       31     3080 
     200   0.36440849    36005        0        0       47     6622 
     300   0.71349144    43617        0        0       62     7700 
     400    1.1126316    48013        0        0       71     8610 
     500    1.5412369    50729        0        0       55     8456 
Loop time of 1.54126 on 1 procs for 500 steps with 50729 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.1071     | 1.1071     | 1.1071     |   0.0 | 71.83
Coll    | 0.14197    | 0.14197    | 0.14197    |   0.0 |  9.21
Sort    | 0.18532    | 0.18532    | 0.18532    |   0.0 | 12.02
Comm    | 0.0030382  | 0.0030382  | 0.0030382  |   0.0 |  0.20
Modify  | 0.10249    | 0.10249    | 0.10249    |   0.0 |  6.65
Output  | 0.00012159 | 0.00012159 | 0.00012159 |   0.0 |  0.01
Other   |            | 0.001254   |            |       |  0.08

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

Particle-moves/CPUsec/proc: 1.14453e+07
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
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1881 deleted particles
  CPU time = 0.00173211 secs
  sort/surf2grid/ghost/inout/particle percent = 35.9532 32.7323 6.0702 2.76669 22.4776
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48848        0        0        0        0 
     600   0.42798162    50538        0        0       67     9912 
     700   0.87143731    51573        0        0       70     9898 
     800    1.3203011    52461        0        0       70    10136 
     900    1.7763095    53020        0        0       60     9492 
    1000    2.2343752    53326        0        0       63     9590 
Loop time of 2.23447 on 1 procs for 500 steps with 53326 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.6381     | 1.6381     | 1.6381     |   0.0 | 73.31
Coll    | 0.21922    | 0.21922    | 0.21922    |   0.0 |  9.81
Sort    | 0.27447    | 0.27447    | 0.27447    |   0.0 | 12.28
Comm    | 0.0050812  | 0.0050812  | 0.0050812  |   0.0 |  0.23
Modify  | 0.095977   | 0.095977   | 0.095977   |   0.0 |  4.30
Output  | 0.00013494 | 0.00013494 | 0.00013494 |   0.0 |  0.01
Other   |            | 0.001495   |            |       |  0.07

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

Particle-moves/CPUsec/proc: 1.16254e+07
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
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1405 deleted particles
  CPU time = 0.00169945 secs
  sort/surf2grid/ghost/inout/particle percent = 36.2654 31.748 5.79405 2.74972 23.4428
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51921        0        0        0        0 
    1100   0.44155359    52243        0        0       58     9786 
    1200   0.89519191    53130        0        0       65     9632 
    1300    1.3543477    53567        0        0       76    10150 
    1400    1.8172803    53940        0        0       78    10052 
    1500     2.283216    54231        0        0       69     9758 
Loop time of 2.28323 on 1 procs for 500 steps with 54231 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.6695     | 1.6695     | 1.6695     |   0.0 | 73.12
Coll    | 0.22884    | 0.22884    | 0.22884    |   0.0 | 10.02
Sort    | 0.28211    | 0.28211    | 0.28211    |   0.0 | 12.36
Comm    | 0.0052974  | 0.0052974  | 0.0052974  |   0.0 |  0.23
Modify  | 0.095896   | 0.095896   | 0.095896   |   0.0 |  4.20
Output  | 0.00012183 | 0.00012183 | 0.00012183 |   0.0 |  0.01
Other   |            | 0.001432   |            |       |  0.06

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

Particle-moves/CPUsec/proc: 1.16965e+07
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
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  786 deleted particles
  CPU time = 0.00167298 secs
  sort/surf2grid/ghost/inout/particle percent = 38.1075 32.8345 5.94271 2.82172 20.2936
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53445        0        0        0        0 
    1600   0.45157981    54383        0        0       59     9142 
    1700   0.92165875    55542        0        0       54     8806 
    1800    1.4024944    56407        0        0       57     9436 
    1900    1.8867753    56697        0        0       59    10150 
    2000    2.3780968    57299        0        0       65    10472 
Loop time of 2.37812 on 1 procs for 500 steps with 57299 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.729      | 1.729      | 1.729      |   0.0 | 72.70
Coll    | 0.24723    | 0.24723    | 0.24723    |   0.0 | 10.40
Sort    | 0.2983     | 0.2983     | 0.2983     |   0.0 | 12.54
Comm    | 0.005403   | 0.005403   | 0.005403   |   0.0 |  0.23
Modify  | 0.096481   | 0.096481   | 0.096481   |   0.0 |  4.06
Output  | 0.00013351 | 0.00013351 | 0.00013351 |   0.0 |  0.01
Other   |            | 0.001559   |            |       |  0.07

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

Particle-moves/CPUsec/proc: 1.17485e+07
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
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  96 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  96.929 96.929 = cell-wise and global flow volume
  1505 deleted particles
  CPU time = 0.00182962 secs
  sort/surf2grid/ghost/inout/particle percent = 36.995 31.1181 5.69455 2.54105 23.6513
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55794        0        0        0        0 
    2100   0.47516513    56563        0        0       60     9464 
    2200   0.96407986    57262        0        0       52     9464 
    2300    1.4550986    57706        0        0       62     9394 
    2400    1.9514256    57892        0        0       69    10556 
    2500    2.4481733    58073        0        0       64    10150 
Loop time of 2.44819 on 1 procs for 500 steps with 58073 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.7799     | 1.7799     | 1.7799     |   0.0 | 72.70
Coll    | 0.25712    | 0.25712    | 0.25712    |   0.0 | 10.50
Sort    | 0.3068     | 0.3068     | 0.3068     |   0.0 | 12.53
Comm    | 0.0055079  | 0.0055079  | 0.0055079  |   0.0 |  0.22
Modify  | 0.097141   | 0.097141   | 0.097141   |   0.0 |  3.97
Output  | 0.00013542 | 0.00013542 | 0.00013542 |   0.0 |  0.01
Other   |            | 0.001629   |            |       |  0.07

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

Particle-moves/CPUsec/proc: 1.17345e+07
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
