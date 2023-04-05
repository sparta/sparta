SPARTA (18 Jul 2022)
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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/ascldap/users/stamoor/sparta_master/src/grid.cpp:465)
Created 400 child grid cells
  CPU time = 0.00112152 secs
  create/ghost percent = 83.5672 16.4328
balance_grid        rcb cell
Balance grid migrated 280 cells
  CPU time = 0.000954151 secs
  reassign/sort/migrate/ghost percent = 66.092 0.274863 14.1429 19.4903

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
  CPU time = 0.00110817 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 26.8072 20.9983 1.31239 40.7272 10.1549 24.5912 1.07573
  surf2grid time = 0.000451326 secs
  map/comm1/comm2/comm3/comm4/split percent = 33.9672 10.2483 9.77285 5.70523 10.6709 21.7644

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
     100  0.037733316    18935        0        0      112     3801 
     200     0.104105    28432        0        0      161     5873 
     300    0.1808064    32176        0        0      183     6306 
     400   0.26159501    33776        0        0      187     6687 
     500    0.3440609    34642        0        0      215     6888 
     600   0.42778516    35095        0        0      140     6581 
     700   0.51186824    35382        0        0      176     6848 
     800   0.59609151    35523        0        0      206     7186 
     900   0.68048215    35651        0        0      186     6793 
    1000   0.76523352    35732        0        0      171     6942 
Loop time of 0.765309 on 4 procs for 1000 steps with 35732 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.11831    | 0.32418    | 0.5321     |  36.1 | 42.36
Coll    | 0.0097067  | 0.028425   | 0.047562   |  11.1 |  3.71
Sort    | 0.023839   | 0.052228   | 0.080976   |  12.4 |  6.82
Comm    | 0.023714   | 0.025929   | 0.028095   |   1.0 |  3.39
Modify  | 0.00051403 | 0.032204   | 0.063956   |  17.6 |  4.21
Output  | 0.00027633 | 0.00089735 | 0.0012615  |   0.0 |  0.12
Other   |            | 0.3014     |            |       | 39.39

Particle moves    = 31095407 (31.1M)
Cells touched     = 35983645 (36M)
Particle comms    = 143676 (0.144M)
Boundary collides = 172693 (0.173M)
Boundary exits    = 175102 (0.175M)
SurfColl checks   = 6041072 (6.04M)
SurfColl occurs   = 173166 (0.173M)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 1.01578e+07
Particle-moves/step: 31095.4
Cell-touches/particle/step: 1.1572
Particle comm iterations/step: 2.145
Particle fraction communicated: 0.00462049
Particle fraction colliding with boundary: 0.00555365
Particle fraction exiting boundary: 0.00563112
Surface-checks/particle/step: 0.194275
Surface-collisions/particle/step: 0.00556886
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 8933 ave 14091 max 3766 min
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

