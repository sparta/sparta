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

boundary	    	o r p

create_box  	    0 10 0 10 -0.5 0.5
Created orthogonal box = (0 0 -0.5) to (10 10 0.5)
create_grid 	    20 20 1
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/ascldap/users/stamoor/sparta_stanmoore1/src/grid.cpp:410)
Created 400 child grid cells
  CPU time = 0.00120592 secs
  create/ghost percent = 86.0024 13.9976
balance_grid        rcb cell
Balance grid migrated 280 cells
  CPU time = 0.00140738 secs
  reassign/sort/migrate/ghost percent = 56.0901 0.440454 11.9939 31.4755

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
  0 0 = number of pushed cells
  48 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  264 88 48 = cells outside/inside/overlapping surfs
  48 = surf cells with 1,2,etc splits
  71.8 71.8 = cell-wise and global flow volume
  CPU time = 0.00218296 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 11.457 14.0236 1.39799 57.82 15.3014 9.13062 0.436872
  surf2grid time = 0.00126219 secs
  map/comm1/comm2/comm3/comm4/split percent = 56.2146 4.19343 3.85342 6.64904 5.61012 18.5115

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
      10 0.00033760071        0        0        0        0        0 
      20 0.00064492226        0        0        0        0        0 
      30  0.001032114        0        0        0        0        0 
      40 0.0013394356        0        0        0        0        0 
      50 0.0016596317        0        0        0        0        0 
      60 0.0019943714        0        0        0        0        0 
      70 0.0022888184        0        0        0        0        0 
      80 0.0025963783        0        0        0        0        0 
      90 0.0029227734        0        0        0        0        0 
     100 0.0091843605    10000        0        0        0        0 
     110  0.015108824    10000        0        0        0        4 
     120  0.019160509     9884        0        0       40      687 
     130  0.024505377     9087        0        0      109     1623 
     140  0.029738426     8103        0        0      105     1935 
     150  0.033743382     7239        0        0       88     1895 
     160  0.037379503     6496        0        0       75     1742 
     170  0.040663242     5887        0        0       63     1455 
     180  0.046078444     5388        0        0       45     1150 
     190  0.048611164     4965        0        0       29      949 
     200  0.055605888    14554        0        0       32      917 
     210  0.060357332    14122        0        0       24      776 
     220  0.065296412    13622        0        0       53     1297 
     230  0.070945024    12532        0        0      117     2167 
     240  0.078140497    11257        0        0      111     2542 
     250  0.085474014    10057        0        0       92     2401 
     260  0.091901779     9047        0        0       69     2111 
     270  0.095951557     8231        0        0       61     1942 
     280  0.099069357     7457        0        0       63     1553 
     290   0.10185289     6825        0        0       45     1309 
     300   0.10909581    16246        0        0       26     1039 
     310   0.11674571    15765        0        0       27      966 
     320   0.12177563    15128        0        0       48     1460 
     330   0.12952662    13929        0        0      104     2360 
     340   0.13491392    12579        0        0      104     2706 
     350   0.13980436    11329        0        0      105     2589 
     360    0.1447165    10292        0        0       64     2251 
     370   0.15026784     9319        0        0       65     1916 
     380   0.15415788     8476        0        0       61     1679 
     390   0.15885019     7716        0        0       46     1418 
     400    0.1663301    17098        0        0       40     1207 
     410   0.17147231    16466        0        0       39     1130 
     420   0.17712569    15738        0        0       64     1665 
     430   0.18263698    14466        0        0      100     2585 
     440   0.18809819    13036        0        0      116     2907 
     450   0.19306874    11721        0        0      111     2638 
     460   0.19989157    10671        0        0       71     2245 
     470    0.2046392     9761        0        0       58     1949 
     480   0.20989323     8991        0        0       57     1752 
     490   0.21304965     8246        0        0       42     1565 
     500   0.22082376    17569        0        0       44     1437 
Loop time of 0.220904 on 4 procs for 500 steps with 17569 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.026866   | 0.082042   | 0.13789    |  19.1 | 37.14
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.010593   | 0.030515   | 0.0378     |   6.6 | 13.81
Modify  | 0.00014973 | 0.012222   | 0.024393   |  10.9 |  5.53
Output  | 0.0023644  | 0.0026904  | 0.0036471  |   1.1 |  1.22
Other   |            | 0.09344    |            |       | 42.30

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

Particle-moves/CPUsec/proc: 4.7781e+06
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

