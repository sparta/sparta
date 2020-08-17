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
create_grid 	    20 20 1
Created 400 child grid cells
  parent cells = 1
  CPU time = 0.000954 secs
  create/ghost percent = 85.8491 14.1509
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.00028 secs
  reassign/sort/migrate/ghost percent = 71.7857 0.714286 6.42857 21.0714

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
  48 = cells with surfs
  104 = total surfs in all grid cells
  3 = max surfs in one grid cell
  0.753486 = min surf-size/cell-size ratio
  0 0 = number of pushed cells
  48 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  264 88 48 = cells outside/inside/overlapping surfs
  48 = surf cells with 1,2,etc splits
  71.8 71.8 = cell-wise and global flow volume
  CPU time = 0.000396 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 31.3131 8.58586 1.26263 42.4242 16.4141 8.83838 0
  surf2grid time = 0.000168 secs
  map/rvous1/rvous2/split percent = 30.9524 44.0476 0.595238 17.8571

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
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 1.51903 1.51903 1.51903
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
     100     0.026976    18917        0        0      114     3863 
     200     0.083163    28423        0        0      169     5690 
     300     0.151907    32220        0        0      170     6437 
     400     0.226616    33622        0        0      195     6564 
     500      0.30047    34565        0        0      191     6887 
     600     0.376667    35086        0        0      188     6666 
     700     0.454522    35363        0        0      209     6961 
     800     0.532314    35418        0        0      192     6473 
     900     0.607611    35419        0        0      206     6679 
    1000     0.705723    35756        0        0      211     7236 
Loop time of 0.705737 on 1 procs for 1000 steps with 35756 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.51485    | 0.51485    | 0.51485    |   0.0 | 72.95
Coll    | 0.086506   | 0.086506   | 0.086506   |   0.0 | 12.26
Sort    | 0.060632   | 0.060632   | 0.060632   |   0.0 |  8.59
Comm    | 0.002698   | 0.002698   | 0.002698   |   0.0 |  0.38
Modify  | 0.040752   | 0.040752   | 0.040752   |   0.0 |  5.77
Output  | 0.000103   | 0.000103   | 0.000103   |   0.0 |  0.01
Other   |            | 0.000199   |            |       |  0.03

Particle moves    = 31064477 (31.1M)
Cells touched     = 35942727 (35.9M)
Particle comms    = 0 (0K)
Boundary collides = 171851 (0.172M)
Boundary exits    = 175003 (0.175M)
SurfColl checks   = 6029960 (6.03M)
SurfColl occurs   = 172928 (0.173M)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 4.40171e+07
Particle-moves/step: 31064.5
Cell-touches/particle/step: 1.15704
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00553207
Particle fraction exiting boundary: 0.00563354
Surface-checks/particle/step: 0.194111
Surface-collisions/particle/step: 0.00556674
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 35756 ave 35756 max 35756 min
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

