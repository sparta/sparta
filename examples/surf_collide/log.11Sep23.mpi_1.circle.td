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
create_grid 	    20 20 1
Created 400 child grid cells
  CPU time = 0.000914102 secs
  create/ghost percent = 91.1169 8.88314
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.0001339 secs
  reassign/sort/migrate/ghost percent = 60.717 0.672143 12.696 25.9149

global		    nrho 1.0 fnum 0.001

species		    air.species N O
mixture		    air N O vstream 100.0 0 0

read_surf           data.circle
  50 points
  50 lines
  2 8 xlo xhi
  2.00592 7.99408 ylo yhi
  0 0 zlo zhi
  0.376743 min line length
  48 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  264 88 48 = cells outside/inside/overlapping surfs
  48 = surf cells with 1,2,etc splits
  71.8 71.8 = cell-wise and global flow volume
  CPU time = 0.000555902 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 22.738 14.1572 1.15128 54.2002 7.75335 13.7254 0.0179888
  surf2grid time = 0.0003013 secs
  map/comm1/comm2/comm3/comm4/split percent = 41.5865 7.10256 14.1719 3.21938 13.4418 17.3913

surf_collide	    1 td 1000.0 #barrier 1000
#surf_collide	    1 td 500.0 bond 500 0 0
#surf_collide	    1 td 300.0 initenergy 0.01 0 0
#surf_collide	    1 td 1000.0 barrier 1000 bond 500 0 0 initenergy 0.01 0 0

surf_modify         all collide 1

collide             vss air air.vss

fix		    in emit/face air xlo twopass # subsonic 0.1 NULL

timestep 	    0.0001

#dump                2 image all 50 image.*.ppm type type pdiam 0.1 #                    surf proc 0.01 size 512 512 zoom 1.75 #                    gline yes 0.005
#dump_modify	    2 pad 4

stats		    100
stats_style	    step cpu np nattempt ncoll nscoll nscheck
run 		    1000
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 1.51894 1.51894 1.51894
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
     100  0.064713446    18931        0        0      119     3949 
     200   0.20427766    28398        0        0      163     5567 
     300   0.38907657    32180        0        0      190     6542 
     400   0.59710814    33853        0        0      185     6629 
     500   0.82319695    34533        0        0      195     6819 
     600   0.99326423    35159        0        0      181     6994 
     700    1.2077339    35410        0        0      184     6587 
     800    1.4139359    35434        0        0      187     7032 
     900    1.6279959    35545        0        0      193     6819 
    1000     1.842262    35581        0        0      174     6748 
Loop time of 1.84234 on 1 procs for 1000 steps with 35581 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.2612     | 1.2612     | 1.2612     |   0.0 | 68.46
Coll    | 0.22901    | 0.22901    | 0.22901    |   0.0 | 12.43
Sort    | 0.21891    | 0.21891    | 0.21891    |   0.0 | 11.88
Comm    | 0.003193   | 0.003193   | 0.003193   |   0.0 |  0.17
Modify  | 0.095522   | 0.095522   | 0.095522   |   0.0 |  5.18
Output  | 0.033992   | 0.033992   | 0.033992   |   0.0 |  1.85
Other   |            | 0.0004775  |            |       |  0.03

Particle moves    = 31094570 (31.1M)
Cells touched     = 35976871 (36M)
Particle comms    = 0 (0K)
Boundary collides = 171944 (0.172M)
Boundary exits    = 175142 (0.175M)
SurfColl checks   = 6061194 (6.06M)
SurfColl occurs   = 173528 (0.174M)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 1.68777e+07
Particle-moves/step: 31094.6
Cell-touches/particle/step: 1.15701
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00552971
Particle fraction exiting boundary: 0.00563256
Surface-checks/particle/step: 0.194928
Surface-collisions/particle/step: 0.00558065
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 35581 ave 35581 max 35581 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      400 ave 400 max 400 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0

