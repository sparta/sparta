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

boundary	    	o r p

create_box  	    0 10 0 10 -0.5 0.5
Created orthogonal box = (0 0 -0.5) to (10 10 0.5)
create_grid 	    20 20 1
Created 400 child grid cells
  CPU time = 0.00101861 secs
  create/ghost percent = 89.9568 10.0432
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000150502 secs
  reassign/sort/migrate/ghost percent = 56.7441 0.531554 11.2962 31.4282

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
  CPU time = 0.000605508 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 23.1544 13.4599 1.15605 53.2948 8.93481 20.1487 0.0165151
  surf2grid time = 0.000322704 secs
  map/comm1/comm2/comm3/comm4/split percent = 43.9418 7.25123 14.6887 3.1298 12.9224 14.9673

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
      10     1.62e-05        0        0        0        0        0 
      20   4.4401e-05        0        0        0        0        0 
      30   7.2401e-05        0        0        0        0        0 
      40   9.9901e-05        0        0        0        0        0 
      50  0.000127202        0        0        0        0        0 
      60  0.000154602        0        0        0        0        0 
      70  0.000182703        0        0        0        0        0 
      80  0.000209803        0        0        0        0        0 
      90  0.000237103        0        0        0        0        0 
     100  0.003683249    10000        0        0        0        0 
     110  0.008782617    10000        0        0        0        6 
     120  0.018672049     9877        0        0       36      651 
     130  0.020652876     9145        0        0       88     1558 
     140  0.022676003     8137        0        0       95     1991 
     150  0.024477327     7226        0        0       77     1949 
     160  0.026093848     6468        0        0       81     1668 
     170   0.03445966     5893        0        0       56     1439 
     180  0.035700676     5410        0        0       43     1138 
     190  0.036799391     5016        0        0       35     1026 
     200  0.040254337    14590        0        0       26      891 
     210  0.050371772    14169        0        0       16      716 
     220  0.052801105    13675        0        0       60     1354 
     230   0.05544994    12583        0        0      119     2190 
     240  0.062348532    11304        0        0      107     2520 
     250  0.064724964    10092        0        0      105     2370 
     260  0.075141203     9153        0        0       57     1994 
     270  0.077043628     8297        0        0       67     1856 
     280   0.07872195     7583        0        0       53     1534 
     290  0.080236171     6910        0        0       49     1410 
     300  0.092837039    16318        0        0       36     1137 
     310  0.095504574    15790        0        0       28     1006 
     320  0.098221611    15154        0        0       58     1653 
     330   0.10398459    13973        0        0      112     2475 
     340    0.1126781    12601        0        0      124     2751 
     350   0.11531464    11323        0        0      109     2578 
     360   0.11935199    10254        0        0       75     2230 
     370   0.13109535     9350        0        0       56     1907 
     380   0.13296137     8556        0        0       51     1808 
     390    0.1346558     7787        0        0       54     1512 
     400   0.13964746    17151        0        0       38     1268 
     410   0.15179402    16541        0        0       31     1126 
     420   0.15467876    15842        0        0       55     1669 
     430     0.157708    14627        0        0      117     2409 
     440   0.16798834    13245        0        0      128     2837 
     450   0.17076138    11926        0        0       98     2655 
     460   0.17475703    10792        0        0       87     2344 
     470   0.18318264     9880        0        0       73     2042 
     480   0.18513177     9042        0        0       58     1813 
     490   0.18687109     8264        0        0       48     1470 
     500   0.19241567    17600        0        0       33     1339 
Loop time of 0.19243 on 1 procs for 500 steps with 17600 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.12833    | 0.12833    | 0.12833    |   0.0 | 66.69
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00044771 | 0.00044771 | 0.00044771 |   0.0 |  0.23
Modify  | 0.013581   | 0.013581   | 0.013581   |   0.0 |  7.06
Output  | 0.049111   | 0.049111   | 0.049111   |   0.0 | 25.52
Other   |            | 0.0009586  |            |       |  0.50

Particle moves    = 4244208 (4.24M)
Cells touched     = 4797333 (4.8M)
Particle comms    = 0 (0K)
Boundary collides = 16239 (16.2K)
Boundary exits    = 7197 (7.2K)
SurfColl checks   = 680716 (0.681M)
SurfColl occurs   = 26659 (26.7K)
Surf reactions    = 26588 (26.6K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.20559e+07
Particle-moves/step: 8488.42
Cell-touches/particle/step: 1.13032
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00382616
Particle fraction exiting boundary: 0.00169572
Surface-checks/particle/step: 0.160387
Surface-collisions/particle/step: 0.00628127
Surf-reactions/particle/step: 0.00626454
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 26588
    reaction O(g) --> O(s): 19515
    reaction CO2(g) --> C(b) + 2O(g): 4
    reaction O(g) + O(s) --> CO2(g): 88
    reaction O(g) --> CO(s): 6046
    reaction O(g) --> CO(g): 581
    reaction O(g) + O(s) --> O(g) + O(g): 354

Particles: 17600 ave 17600 max 17600 min
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

