SPARTA (13 Apr 2023)
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
  CPU time = 0.000864802 secs
  create/ghost percent = 96.2419 3.75809
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 9.29e-05 secs
  reassign/sort/migrate/ghost percent = 71.4747 0.430571 15.5005 12.5942

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
  CPU time = 0.000479501 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 27.3203 16.3295 0.45881 49.5098 6.38163 12.1587 0.04171
  surf2grid time = 0.0002374 secs
  map/comm1/comm2/comm3/comm4/split percent = 40.396 9.64617 7.49789 3.58045 22.1567 14.8273
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
     100  0.059666635    20918        0        0       31     3080 
     200   0.20024515    36005        0        0       47     6622 
     300   0.39352849    43617        0        0       62     7700 
     400   0.62601031    48013        0        0       71     8610 
     500   0.87762548    50729        0        0       55     8456 
Loop time of 0.877691 on 1 procs for 500 steps with 50729 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.56859    | 0.56859    | 0.56859    |   0.0 | 64.78
Coll    | 0.16899    | 0.16899    | 0.16899    |   0.0 | 19.25
Sort    | 0.071578   | 0.071578   | 0.071578   |   0.0 |  8.16
Comm    | 0.0010602  | 0.0010602  | 0.0010602  |   0.0 |  0.12
Modify  | 0.054952   | 0.054952   | 0.054952   |   0.0 |  6.26
Output  | 0.012299   | 0.012299   | 0.012299   |   0.0 |  1.40
Other   |            | 0.0002139  |            |       |  0.02

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

Particle-moves/CPUsec/proc: 2.00983e+07
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
  CPU time = 0.000564202 secs
  sort/surf2grid/ghost/inout/particle percent = 30.2199 28.1637 2.53455 10.2446 28.8374
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
     500            0    48848        0        0        0        0 
     600   0.24661095    50538        0        0       67     9912 
     700   0.50434184    51573        0        0       70     9898 
     800    0.7547388    52461        0        0       70    10136 
     900    1.0253664    53020        0        0       60     9492 
    1000    1.2950168    53326        0        0       63     9590 
Loop time of 1.2951 on 1 procs for 500 steps with 53326 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.76477    | 0.76477    | 0.76477    |   0.0 | 59.05
Coll    | 0.2787     | 0.2787     | 0.2787     |   0.0 | 21.52
Sort    | 0.17256    | 0.17256    | 0.17256    |   0.0 | 13.32
Comm    | 0.0016901  | 0.0016901  | 0.0016901  |   0.0 |  0.13
Modify  | 0.050947   | 0.050947   | 0.050947   |   0.0 |  3.93
Output  | 0.026165   | 0.026165   | 0.026165   |   0.0 |  2.02
Other   |            | 0.0002682  |            |       |  0.02

Particle moves    = 25976616 (26M)
Cells touched     = 27611293 (27.6M)
Particle comms    = 0 (0K)
Boundary collides = 86022 (86K)
Boundary exits    = 100890 (0.101M)
SurfColl checks   = 4739322 (4.74M)
SurfColl occurs   = 30852 (30.9K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.00576e+07
Particle-moves/step: 51953.2
Cell-touches/particle/step: 1.06293
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00331152
Particle fraction exiting boundary: 0.00388388
Surface-checks/particle/step: 0.182446
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
  CPU time = 0.000646901 secs
  sort/surf2grid/ghost/inout/particle percent = 34.6885 25.2436 2.27237 8.85761 28.938
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1000            0    51921        0        0        0        0 
    1100   0.24851836    52243        0        0       58     9786 
    1200   0.51211825    53130        0        0       65     9632 
    1300   0.76506822    53567        0        0       76    10150 
    1400    1.0331651    53940        0        0       78    10052 
    1500    1.3146307    54231        0        0       69     9758 
Loop time of 1.31471 on 1 procs for 500 steps with 54231 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.80401    | 0.80401    | 0.80401    |   0.0 | 61.15
Coll    | 0.2349     | 0.2349     | 0.2349     |   0.0 | 17.87
Sort    | 0.19597    | 0.19597    | 0.19597    |   0.0 | 14.91
Comm    | 0.0097084  | 0.0097084  | 0.0097084  |   0.0 |  0.74
Modify  | 0.044805   | 0.044805   | 0.044805   |   0.0 |  3.41
Output  | 0.025049   | 0.025049   | 0.025049   |   0.0 |  1.91
Other   |            | 0.0002641  |            |       |  0.02

Particle moves    = 26705887 (26.7M)
Cells touched     = 28349567 (28.3M)
Particle comms    = 0 (0K)
Boundary collides = 87110 (87.1K)
Boundary exits    = 102973 (0.103M)
SurfColl checks   = 4705148 (4.71M)
SurfColl occurs   = 30430 (30.4K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.03131e+07
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
  CPU time = 0.000909002 secs
  sort/surf2grid/ghost/inout/particle percent = 42.3763 18.1298 1.67216 6.36962 31.4522
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    1500            0    53445        0        0        0        0 
    1600   0.21277618    54383        0        0       59     9142 
    1700    0.3980817    55542        0        0       54     8806 
    1800   0.63294823    56407        0        0       57     9436 
    1900   0.92622379    56697        0        0       59    10150 
    2000    1.2221025    57299        0        0       65    10472 
Loop time of 1.22219 on 1 procs for 500 steps with 57299 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.7613     | 0.7613     | 0.7613     |   0.0 | 62.29
Coll    | 0.25365    | 0.25365    | 0.25365    |   0.0 | 20.75
Sort    | 0.1355     | 0.1355     | 0.1355     |   0.0 | 11.09
Comm    | 0.0020268  | 0.0020268  | 0.0020268  |   0.0 |  0.17
Modify  | 0.055244   | 0.055244   | 0.055244   |   0.0 |  4.52
Output  | 0.014157   | 0.014157   | 0.014157   |   0.0 |  1.16
Other   |            | 0.000316   |            |       |  0.03

Particle moves    = 27939282 (27.9M)
Cells touched     = 29655039 (29.7M)
Particle comms    = 0 (0K)
Boundary collides = 91806 (91.8K)
Boundary exits    = 101497 (0.101M)
SurfColl checks   = 4651626 (4.65M)
SurfColl occurs   = 29739 (29.7K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.28601e+07
Particle-moves/step: 55878.6
Cell-touches/particle/step: 1.06141
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00328591
Particle fraction exiting boundary: 0.00363277
Surface-checks/particle/step: 0.166491
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
  CPU time = 0.000919502 secs
  sort/surf2grid/ghost/inout/particle percent = 42.8929 18.2491 1.66394 6.34039 30.8537
run 		    500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 6.75 6.75 6.75
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 8.26894 8.26894 8.26894
Step CPU Np Natt Ncoll Nscoll Nscheck 
    2000            0    55794        0        0        0        0 
    2100   0.27788373    56563        0        0       60     9464 
    2200   0.57017018    57262        0        0       52     9464 
    2300   0.85653753    57706        0        0       62     9394 
    2400    1.1439505    57892        0        0       69    10556 
    2500    1.4388064    58073        0        0       64    10150 
Loop time of 1.43883 on 1 procs for 500 steps with 58073 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.90061    | 0.90061    | 0.90061    |   0.0 | 62.59
Coll    | 0.29293    | 0.29293    | 0.29293    |   0.0 | 20.36
Sort    | 0.18686    | 0.18686    | 0.18686    |   0.0 | 12.99
Comm    | 0.0018506  | 0.0018506  | 0.0018506  |   0.0 |  0.13
Modify  | 0.046429   | 0.046429   | 0.046429   |   0.0 |  3.23
Output  | 0.0098612  | 0.0098612  | 0.0098612  |   0.0 |  0.69
Other   |            | 0.0002932  |            |       |  0.02

Particle moves    = 28728268 (28.7M)
Cells touched     = 30479598 (30.5M)
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

Particle-moves/CPUsec/proc: 1.99664e+07
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
