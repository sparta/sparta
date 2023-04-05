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
  CPU time = 0.000921512 secs
  create/ghost percent = 95.8655 4.13451
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 9.7301e-05 secs
  reassign/sort/migrate/ghost percent = 70.7094 0.308322 14.3883 14.5939

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
  CPU time = 0.000503007 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 28.2704 16.1232 0.497011 47.9721 7.13728 11.6303 0.0198804
  surf2grid time = 0.000241303 secs
  map/comm1/comm2/comm3/comm4/split percent = 41.9394 9.82168 7.91577 3.35678 21.4258 13.6758
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
     100  0.051362301    20918        0        0       31     3080 
     200   0.19918602    36005        0        0       47     6622 
     300     0.402934    43617        0        0       62     7700 
     400   0.62777477    48013        0        0       71     8610 
     500   0.86722625    50729        0        0       55     8456 
Loop time of 0.877272 on 1 procs for 500 steps with 50729 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.52041    | 0.52041    | 0.52041    |   0.0 | 59.32
Coll    | 0.16653    | 0.16653    | 0.16653    |   0.0 | 18.98
Sort    | 0.10545    | 0.10545    | 0.10545    |   0.0 | 12.02
Comm    | 0.0012263  | 0.0012263  | 0.0012263  |   0.0 |  0.14
Modify  | 0.053088   | 0.053088   | 0.053088   |   0.0 |  6.05
Output  | 0.030241   | 0.030241   | 0.030241   |   0.0 |  3.45
Other   |            | 0.0003224  |            |       |  0.04

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

Particle-moves/CPUsec/proc: 2.01079e+07
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
  CPU time = 0.00113531 secs
  sort/surf2grid/ghost/inout/particle percent = 46.349 16.8414 1.55023 5.30258 29.9568
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48848        0        0        0        0 
     600   0.24747768    50538        0        0       67     9912 
     700   0.50039943    51573        0        0       70     9898 
     800   0.76691408    52461        0        0       70    10136 
     900    1.0431758    53020        0        0       60     9492 
    1000    1.3028676    53326        0        0       63     9590 
Loop time of 1.3029 on 1 procs for 500 steps with 53326 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.84714    | 0.84714    | 0.84714    |   0.0 | 65.02
Coll    | 0.22632    | 0.22632    | 0.22632    |   0.0 | 17.37
Sort    | 0.18999    | 0.18999    | 0.18999    |   0.0 | 14.58
Comm    | 0.0020714  | 0.0020714  | 0.0020714  |   0.0 |  0.16
Modify  | 0.027363   | 0.027363   | 0.027363   |   0.0 |  2.10
Output  | 0.0095791  | 0.0095791  | 0.0095791  |   0.0 |  0.74
Other   |            | 0.0004352  |            |       |  0.03

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

Particle-moves/CPUsec/proc: 1.99375e+07
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
  CPU time = 0.00075871 secs
  sort/surf2grid/ghost/inout/particle percent = 46.4873 20.5221 2.28019 3.70365 27.0068
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51921        0        0        0        0 
    1100   0.25197994    52243        0        0       58     9786 
    1200   0.53110406    53130        0        0       65     9632 
    1300   0.80857564    53567        0        0       76    10150 
    1400    1.0958821    53940        0        0       78    10052 
    1500    1.3788617    54231        0        0       69     9758 
Loop time of 1.3789 on 1 procs for 500 steps with 54231 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.76471    | 0.76471    | 0.76471    |   0.0 | 55.46
Coll    | 0.28652    | 0.28652    | 0.28652    |   0.0 | 20.78
Sort    | 0.24339    | 0.24339    | 0.24339    |   0.0 | 17.65
Comm    | 0.0022864  | 0.0022864  | 0.0022864  |   0.0 |  0.17
Modify  | 0.08124    | 0.08124    | 0.08124    |   0.0 |  5.89
Output  | 0.0002673  | 0.0002673  | 0.0002673  |   0.0 |  0.02
Other   |            | 0.000479   |            |       |  0.03

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

Particle-moves/CPUsec/proc: 1.93676e+07
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
  CPU time = 0.000654509 secs
  sort/surf2grid/ghost/inout/particle percent = 41.7875 23.5601 2.7196 4.4308 27.502
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53445        0        0        0        0 
    1600   0.28016613    54383        0        0       59     9142 
    1700   0.55180494    55542        0        0       54     8806 
    1800   0.82357105    56407        0        0       57     9436 
    1900    1.1111329    56697        0        0       59    10150 
    2000    1.3982248    57299        0        0       65    10472 
Loop time of 1.39826 on 1 procs for 500 steps with 57299 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.84473    | 0.84473    | 0.84473    |   0.0 | 60.41
Coll    | 0.26042    | 0.26042    | 0.26042    |   0.0 | 18.62
Sort    | 0.22786    | 0.22786    | 0.22786    |   0.0 | 16.30
Comm    | 0.0021237  | 0.0021237  | 0.0021237  |   0.0 |  0.15
Modify  | 0.056212   | 0.056212   | 0.056212   |   0.0 |  4.02
Output  | 0.0064636  | 0.0064636  | 0.0064636  |   0.0 |  0.46
Other   |            | 0.0004492  |            |       |  0.03

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

Particle-moves/CPUsec/proc: 1.99815e+07
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
  CPU time = 0.000653109 secs
  sort/surf2grid/ghost/inout/particle percent = 40.0091 23.5645 2.58778 4.54748 29.2911
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55794        0        0        0        0 
    2100   0.29774377    56563        0        0       60     9464 
    2200   0.59165591    57262        0        0       52     9464 
    2300   0.88772828    57706        0        0       62     9394 
    2400    1.1741324    57892        0        0       69    10556 
    2500     1.466268    58073        0        0       64    10150 
Loop time of 1.47068 on 1 procs for 500 steps with 58073 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.88561    | 0.88561    | 0.88561    |   0.0 | 60.22
Coll    | 0.28105    | 0.28105    | 0.28105    |   0.0 | 19.11
Sort    | 0.22399    | 0.22399    | 0.22399    |   0.0 | 15.23
Comm    | 0.002352   | 0.002352   | 0.002352   |   0.0 |  0.16
Modify  | 0.06758    | 0.06758    | 0.06758    |   0.0 |  4.60
Output  | 0.0094276  | 0.0094276  | 0.0094276  |   0.0 |  0.64
Other   |            | 0.0006691  |            |       |  0.05

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

Particle-moves/CPUsec/proc: 1.9534e+07
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
