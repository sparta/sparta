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
create_grid 	    20 20 1
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/runner/work/sparta/sparta/src/grid.cpp:465)
Created 400 child grid cells
  CPU time = 0.0018946 secs
  create/ghost percent = 92.3995 7.60054
balance_grid        rcb cell
Balance grid migrated 280 cells
  CPU time = 0.000937101 secs
  reassign/sort/migrate/ghost percent = 56.2694 0.757656 21.097 21.876

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
  CPU time = 0.0011936 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 21.875 12.1649 1.04725 50.1173 14.7956 15.3736 0.569705
  surf2grid time = 0.000598201 secs
  map/comm1/comm2/comm3/comm4/split percent = 34.9882 12.3871 9.59544 5.34937 16.182 13.3567

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
     100  0.035419521    18935        0        0      112     3801 
     200   0.10314016    28432        0        0      161     5873 
     300    0.1803037    32176        0        0      183     6306 
     400   0.26178225    33776        0        0      187     6687 
     500    0.3463658    34642        0        0      215     6888 
     600   0.43138905    35095        0        0      140     6581 
     700    0.5170689    35410        0        0      189     6885 
     800   0.60353015    35496        0        0      202     7032 
     900    0.6893619    35577        0        0      174     6782 
    1000   0.77537265    35734        0        0      187     7107 
Loop time of 0.775369 on 4 procs for 1000 steps with 35734 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.056838   | 0.15629    | 0.25618    |  25.1 | 20.16
Coll    | 0.0085185  | 0.02285    | 0.037444   |   9.5 |  2.95
Sort    | 0.0089709  | 0.020959   | 0.033291   |   8.2 |  2.70
Comm    | 0.066721   | 0.2627     | 0.47451    |  35.4 | 33.88
Modify  | 0.0001989  | 0.01436    | 0.028529   |  11.8 |  1.85
Output  | 0.0014907  | 0.0017735  | 0.0022708  |   0.7 |  0.23
Other   |            | 0.2964     |            |       | 38.23

Particle moves    = 30914149 (30.9M)
Cells touched     = 35975211 (36M)
Particle comms    = 143491 (0.143M)
Boundary collides = 172623 (0.173M)
Boundary exits    = 175100 (0.175M)
SurfColl checks   = 6038112 (6.04M)
SurfColl occurs   = 173223 (0.173M)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 9.96756e+06
Particle-moves/step: 30914.1
Cell-touches/particle/step: 1.16371
Particle comm iterations/step: 2.154
Particle fraction communicated: 0.0046416
Particle fraction colliding with boundary: 0.00558395
Particle fraction exiting boundary: 0.00566407
Surface-checks/particle/step: 0.195319
Surface-collisions/particle/step: 0.00560336
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 8933.5 ave 14085 max 3775 min
Histogram: 2 0 0 0 0 0 0 0 0 2
Cells:      100 ave 100 max 100 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 21 ave 21 max 21 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 21 ave 21 max 21 min
Histogram: 4 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0

