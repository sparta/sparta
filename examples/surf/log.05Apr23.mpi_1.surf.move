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
  CPU time = 0.00119569 secs
  create/ghost percent = 92.5562 7.44383
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000215491 secs
  reassign/sort/migrate/ghost percent = 67.4724 0.572646 17.8467 14.1082

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
  CPU time = 0.00105937 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 25.8091 19.9273 0.895058 45.6403 7.72825 13.639 0.0814633
  surf2grid time = 0.000483501 secs
  map/comm1/comm2/comm3/comm4/split percent = 33.292 10.9934 8.97578 5.05956 22.4335 16.6641
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
     100  0.089142939    20918        0        0       31     3080 
     200    0.3263613    36005        0        0       47     6622 
     300   0.65665248    43617        0        0       62     7700 
     400    1.0272245    48013        0        0       71     8610 
     500    1.4403602    50729        0        0       55     8456 
Loop time of 1.44051 on 1 procs for 500 steps with 50729 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.78093    | 0.78093    | 0.78093    |   0.0 | 54.21
Coll    | 0.22399    | 0.22399    | 0.22399    |   0.0 | 15.55
Sort    | 0.35009    | 0.35009    | 0.35009    |   0.0 | 24.30
Comm    | 0.0068739  | 0.0068739  | 0.0068739  |   0.0 |  0.48
Modify  | 0.075853   | 0.075853   | 0.075853   |   0.0 |  5.27
Output  | 0.0011473  | 0.0011473  | 0.0011473  |   0.0 |  0.08
Other   |            | 0.001632   |            |       |  0.11

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

Particle-moves/CPUsec/proc: 1.22458e+07
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
  CPU time = 0.00322556 secs
  sort/surf2grid/ghost/inout/particle percent = 38.4667 13.8461 1.24344 4.82682 41.6169
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48848        0        0        0        0 
     600   0.38923094    50538        0        0       67     9912 
     700   0.81263896    51573        0        0       70     9898 
     800    1.1849679    52461        0        0       70    10136 
     900    1.5539786    53020        0        0       60     9492 
    1000    1.9257403    53326        0        0       63     9590 
Loop time of 1.92588 on 1 procs for 500 steps with 53326 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.105      | 1.105      | 1.105      |   0.0 | 57.38
Coll    | 0.25223    | 0.25223    | 0.25223    |   0.0 | 13.10
Sort    | 0.49702    | 0.49702    | 0.49702    |   0.0 | 25.81
Comm    | 0.0096155  | 0.0096155  | 0.0096155  |   0.0 |  0.50
Modify  | 0.058997   | 0.058997   | 0.058997   |   0.0 |  3.06
Output  | 0.0011872  | 0.0011872  | 0.0011872  |   0.0 |  0.06
Other   |            | 0.001802   |            |       |  0.09

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

Particle-moves/CPUsec/proc: 1.34882e+07
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
  CPU time = 0.00272866 secs
  sort/surf2grid/ghost/inout/particle percent = 35.9629 15.5947 1.32585 5.16096 41.9555
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51921        0        0        0        0 
    1100   0.34552243    52243        0        0       58     9786 
    1200    0.6970671    53130        0        0       65     9632 
    1300    1.0563758    53567        0        0       76    10150 
    1400    1.4171591    53940        0        0       78    10052 
    1500    1.7786954    54231        0        0       69     9758 
Loop time of 1.77883 on 1 procs for 500 steps with 54231 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.0339     | 1.0339     | 1.0339     |   0.0 | 58.12
Coll    | 0.22157    | 0.22157    | 0.22157    |   0.0 | 12.46
Sort    | 0.45482    | 0.45482    | 0.45482    |   0.0 | 25.57
Comm    | 0.0088362  | 0.0088362  | 0.0088362  |   0.0 |  0.50
Modify  | 0.057125   | 0.057125   | 0.057125   |   0.0 |  3.21
Output  | 0.00089134 | 0.00089134 | 0.00089134 |   0.0 |  0.05
Other   |            | 0.00164    |            |       |  0.09

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

Particle-moves/CPUsec/proc: 1.50132e+07
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
  CPU time = 0.00270623 secs
  sort/surf2grid/ghost/inout/particle percent = 36.9604 15.6127 1.37741 5.17011 40.8794
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53445        0        0        0        0 
    1600   0.43106116    54383        0        0       59     9142 
    1700   0.92481661    55542        0        0       54     8806 
    1800    1.4078466    56407        0        0       57     9436 
    1900    1.8561204    56697        0        0       59    10150 
    2000    2.2780054    57299        0        0       65    10472 
Loop time of 2.27825 on 1 procs for 500 steps with 57299 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.2642     | 1.2642     | 1.2642     |   0.0 | 55.49
Coll    | 0.34184    | 0.34184    | 0.34184    |   0.0 | 15.00
Sort    | 0.59577    | 0.59577    | 0.59577    |   0.0 | 26.15
Comm    | 0.010727   | 0.010727   | 0.010727   |   0.0 |  0.47
Modify  | 0.062598   | 0.062598   | 0.062598   |   0.0 |  2.75
Output  | 0.0010772  | 0.0010772  | 0.0010772  |   0.0 |  0.05
Other   |            | 0.00201    |            |       |  0.09

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

Particle-moves/CPUsec/proc: 1.22635e+07
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
  CPU time = 0.00283091 secs
  sort/surf2grid/ghost/inout/particle percent = 50.1968 16.6012 1.50867 5.55879 26.1345
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55794        0        0        0        0 
    2100    0.4276616    56563        0        0       60     9464 
    2200   0.84402427    57262        0        0       52     9464 
    2300    1.2753542    57706        0        0       62     9394 
    2400    1.6980971    57892        0        0       69    10556 
    2500    2.1130264    58073        0        0       64    10150 
Loop time of 2.11307 on 1 procs for 500 steps with 58073 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.2245     | 1.2245     | 1.2245     |   0.0 | 57.95
Coll    | 0.27928    | 0.27928    | 0.27928    |   0.0 | 13.22
Sort    | 0.53397    | 0.53397    | 0.53397    |   0.0 | 25.27
Comm    | 0.01045    | 0.01045    | 0.01045    |   0.0 |  0.49
Modify  | 0.062509   | 0.062509   | 0.062509   |   0.0 |  2.96
Output  | 0.0005416  | 0.0005416  | 0.0005416  |   0.0 |  0.03
Other   |            | 0.001806   |            |       |  0.09

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

Particle-moves/CPUsec/proc: 1.35955e+07
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
