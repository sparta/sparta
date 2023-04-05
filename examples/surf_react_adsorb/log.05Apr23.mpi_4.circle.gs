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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/runner/work/sparta/sparta/src/grid.cpp:465)
Created 400 child grid cells
  CPU time = 0.00236343 secs
  create/ghost percent = 93.704 6.29598
balance_grid        rcb cell
Balance grid migrated 280 cells
  CPU time = 0.000936112 secs
  reassign/sort/migrate/ghost percent = 53.2743 0.630266 20.8204 25.2751

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
  CPU time = 0.00125772 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 20.2353 11.0996 1.07337 51.5862 16.0054 12.8886 1.62994
  surf2grid time = 0.000648808 secs
  map/comm1/comm2/comm3/comm4/split percent = 33.9396 11.4673 12.762 5.51781 13.7793 10.7738

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
      10  0.000561007        0        0        0        0        0 
      20  0.001084214        0        0        0        0        0 
      30   0.00160782        0        0        0        0        0 
      40  0.002128327        0        0        0        0        0 
      50  0.002696634        0        0        0        0        0 
      60  0.003238841        0        0        0        0        0 
      70  0.003761048        0        0        0        0        0 
      80  0.004291355        0        0        0        0        0 
      90  0.004813061        0        0        0        0        0 
     100  0.008925114    10000        0        0        0        0 
     110  0.011423046    10000        0        0        0        4 
     120  0.014336983     9884        0        0       40      687 
     130  0.018123531     9087        0        0      109     1623 
     140  0.021691377     8103        0        0      105     1935 
     150  0.024422312     7239        0        0       88     1895 
     160  0.026831142     6496        0        0       75     1742 
     170  0.028904569     5887        0        0       63     1455 
     180  0.030888894     5388        0        0       45     1150 
     190  0.032548415     4965        0        0       29      949 
     200  0.036653868    14554        0        0       32      917 
     210  0.039782808    14122        0        0       24      776 
     220  0.042932048    13622        0        0       53     1297 
     230  0.046387892    12532        0        0      117     2167 
     240  0.049825936    11257        0        0      111     2542 
     250  0.052932276    10057        0        0       92     2401 
     260  0.055905814     9047        0        0       69     2111 
     270  0.058374345     8231        0        0       61     1942 
     280  0.060588374     7457        0        0       63     1553 
     290  0.062488698     6825        0        0       45     1309 
     300  0.066687951    16246        0        0       26     1039 
     310  0.069895492    15765        0        0       27      966 
     320  0.073275036    15128        0        0       48     1460 
     330  0.076977783    13929        0        0      104     2360 
     340  0.081387839    12579        0        0      104     2706 
     350   0.08461558    11329        0        0      105     2589 
     360  0.087448617    10292        0        0       64     2251 
     370  0.090016849     9319        0        0       65     1916 
     380   0.09241298     8476        0        0       61     1679 
     390  0.094406605     7716        0        0       46     1418 
     400  0.098832462    17098        0        0       40     1207 
     410   0.10229561    16466        0        0       39     1130 
     420   0.10550945    15738        0        0       64     1665 
     430   0.10912509    14466        0        0      100     2585 
     440   0.11268134    13036        0        0      116     2907 
     450   0.11606068    11721        0        0      111     2638 
     460   0.11885502    10671        0        0       71     2245 
     470   0.12149185     9761        0        0       58     1949 
     480   0.12386448     8991        0        0       57     1752 
     490   0.12593281     8246        0        0       42     1565 
     500   0.13025526    17569        0        0       44     1437 
Loop time of 0.130325 on 4 procs for 500 steps with 17569 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0063249  | 0.021047   | 0.035521   |   9.8 | 16.15
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.030195   | 0.034876   | 0.040385   |   2.0 | 26.76
Modify  | 5.4303e-05 | 0.0038625  | 0.0085171  |   6.2 |  2.96
Output  | 0.004229   | 0.0045267  | 0.0049212  |   0.4 |  3.47
Other   |            | 0.06601    |            |       | 50.65

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

Particle-moves/CPUsec/proc: 8.09904e+06
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

