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

boundary	    	o r p

create_box  	    0 10 0 10 -0.5 0.5
Created orthogonal box = (0 0 -0.5) to (10 10 0.5)
create_grid 	    20 20 1
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/me/sparta_master/src/grid.cpp:465)
Created 400 child grid cells
  CPU time = 0.00202165 secs
  create/ghost percent = 89.219 10.781
balance_grid        rcb cell
Balance grid migrated 280 cells
  CPU time = 0.000883571 secs
  reassign/sort/migrate/ghost percent = 56.602 0.479079 20.9132 22.0057

global		    	nrho 1.0 fnum 0.001

species		    	air.species O CO CO2 C
mixture		    	air O vstream 100.0 0 0

mixture             air O   frac 1.0
mixture             air CO  frac 0.0
mixture             air CO2 frac 0.0
mixture             air C   frac 0.0

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
  CPU time = 0.00128403 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 26.6579 20.815 0.60279 41.6573 10.2671 20.7898 0.296333
  surf2grid time = 0.000534892 secs
  map/comm1/comm2/comm3/comm4/split percent = 30.1715 10.665 10.5328 6.5344 15.2109 21.0263

surf_collide        1 cll 300.0 0.5 0.5 0.5 0.5

################################### SURF REACT ADSORB ######################################

#surf_react          adsorb_test_gs1 adsorb gs sample-GS_1.surf nsync 1 surf 1000 6.022e9 O CO
#surf_modify         all collide 1 react adsorb_test_gs1

surf_react          adsorb_test_gs2 adsorb gs sample-GS_2.surf nsync 1 surf 1000 6.022e9 O CO
surf_modify         all collide 1 react adsorb_test_gs2

############################################################################################

#collide            vss air air.vss

fix		    		in emit/face air xlo nevery 100 n 10000 perspecies no twopass # subsonic 0.1 NULL

timestep 	    	0.0001

#dump                2 image all 50 image.*.ppm type type pdiam 0.1 surf proc 0.01 size 512 512 zoom 1.75 gline yes 0.005
#dump_modify	     2 pad 4

stats		    	10
stats_style	    	step cpu np nattempt ncoll nscoll nscheck
run 		    	500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 1.51894 1.51894 1.51894
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
      10  0.000437681        0        0        0        0        0 
      20  0.000812159        0        0        0        0        0 
      30  0.001169112        0        0        0        0        0 
      40  0.001525587        0        0        0        0        0 
      50  0.001877406        0        0        0        0        0 
      60  0.002226673        0        0        0        0        0 
      70  0.002584403        0        0        0        0        0 
      80  0.002934794        0        0        0        0        0 
      90   0.00328802        0        0        0        0        0 
     100  0.008337219    10000        0        0        0        0 
     110  0.010774031    10000        0        0        0        4 
     120  0.013485808     9884        0        0       40      687 
     130  0.018150696     9087        0        0      109     1623 
     140  0.022852725     8103        0        0      105     1935 
     150  0.025877095     7239        0        0       88     1895 
     160  0.028482689     6496        0        0       75     1742 
     170  0.030729936     5887        0        0       63     1455 
     180  0.032701844     5388        0        0       45     1150 
     190  0.034372843     4965        0        0       29      949 
     200  0.039006325    14554        0        0       32      917 
     210  0.042252435    14122        0        0       24      776 
     220  0.045640908    13622        0        0       53     1297 
     230   0.04937174    12532        0        0      117     2167 
     240  0.053175161    11257        0        0      111     2542 
     250  0.056553338    10057        0        0       92     2401 
     260  0.059526617     9047        0        0       69     2111 
     270  0.062124857     8231        0        0       61     1942 
     280  0.064406126     7457        0        0       63     1553 
     290  0.066352784     6825        0        0       45     1309 
     300  0.071396219    16246        0        0       26     1039 
     310  0.074839751    15765        0        0       27      966 
     320  0.078358027    15128        0        0       48     1460 
     330  0.082200251    13929        0        0      104     2360 
     340  0.086075279    12579        0        0      104     2706 
     350  0.089577397    11329        0        0      105     2589 
     360   0.09292693    10292        0        0       64     2251 
     370  0.095706485     9319        0        0       65     1916 
     380  0.098078968     8476        0        0       61     1679 
     390    0.1002586     7716        0        0       46     1418 
     400   0.10532429    17098        0        0       40     1207 
     410   0.10893717    16466        0        0       39     1130 
     420   0.11263957    15738        0        0       64     1665 
     430   0.11664843    14466        0        0      100     2585 
     440   0.12064163    13036        0        0      116     2907 
     450    0.1242782    11721        0        0      111     2638 
     460   0.12748639    10671        0        0       71     2245 
     470   0.13031494     9761        0        0       58     1949 
     480   0.13274615     8991        0        0       57     1752 
     490   0.13492521     8246        0        0       42     1565 
     500   0.14019701    17569        0        0       44     1437 
Loop time of 0.140335 on 4 procs for 500 steps with 17569 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.014692   | 0.053065   | 0.091729   |  16.5 | 37.81
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.012794   | 0.013186   | 0.013913   |   0.4 |  9.40
Modify  | 0.0001182  | 0.0082251  | 0.016483   |   8.9 |  5.86
Output  | 0.0018749  | 0.0025116  | 0.0043898  |   2.2 |  1.79
Other   |            | 0.06335    |            |       | 45.14

Particle moves    = 4222016 (4.22M)
Cells touched     = 4774461 (4.77M)
Particle comms    = 17445 (17.4K)
Boundary collides = 16286 (16.3K)
Boundary exits    = 7138 (7.14K)
SurfColl checks   = 675752 (0.676M)
SurfColl occurs   = 26709 (26.7K)
Surf reactions    = 26643 (26.6K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 7.52134e+06
Particle-moves/step: 8444.03
Cell-touches/particle/step: 1.13085
Particle comm iterations/step: 1.804
Particle fraction communicated: 0.00413191
Particle fraction colliding with boundary: 0.0038574
Particle fraction exiting boundary: 0.00169066
Surface-checks/particle/step: 0.160054
Surface-collisions/particle/step: 0.00632612
Surf-reactions/particle/step: 0.00631049
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 26643
    reaction O(g) --> O(s): 19687
    reaction CO2(g) --> C(b) + 2O(g): 1
    reaction O(g) + O(s) --> CO2(g): 98
    reaction O(g) --> CO(s): 5955
    reaction O(g) --> CO(g): 554
    reaction O(g) + O(s) --> O(g) + O(g): 348

Particles: 4392.25 ave 7359 max 1437 min
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

