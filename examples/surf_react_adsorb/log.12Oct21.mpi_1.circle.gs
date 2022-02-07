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
Created 400 child grid cells
  CPU time = 0.00146699 secs
  create/ghost percent = 60.2633 39.7367
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000772238 secs
  reassign/sort/migrate/ghost percent = 35.0108 0.802717 13.8314 50.355

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
  CPU time = 0.0017581 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 11.0523 7.39083 3.48522 55.099 22.9726 6.78058 0.0135612
  surf2grid time = 0.000968695 secs
  map/comm1/comm2/comm3/comm4/split percent = 41.8902 5.21782 13.2907 3.29806 7.11297 23.2587

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
      10 0.00017666817        0        0        0        0        0 
      20 0.00027704239        0        0        0        0        0 
      30 0.0003824234        0        0        0        0        0 
      40 0.00047850609        0        0        0        0        0 
      50 0.00057411194        0        0        0        0        0 
      60 0.00067019463        0        0        0        0        0 
      70 0.00076627731        0        0        0        0        0 
      80 0.00093650818        0        0        0        0        0 
      90 0.0010697842        0        0        0        0        0 
     100  0.012127161    10000        0        0        0        0 
     110  0.017892599    10000        0        0        0        6 
     120  0.024326086     9877        0        0       36      651 
     130  0.031729937     9145        0        0       88     1558 
     140  0.039266109     8137        0        0       95     1991 
     150  0.046339512     7226        0        0       77     1949 
     160  0.052739143     6468        0        0       81     1668 
     170  0.058418036     5893        0        0       56     1439 
     180  0.063284874     5410        0        0       43     1138 
     190  0.067784309     5016        0        0       35     1026 
     200  0.082171202    14590        0        0       26      891 
     210  0.091641188    14169        0        0       16      716 
     220   0.10151839    13675        0        0       60     1354 
     230   0.11205006    12583        0        0      119     2190 
     240   0.12239194    11304        0        0      107     2520 
     250   0.13219309    10092        0        0      105     2370 
     260    0.1404717     9153        0        0       57     1994 
     270   0.14818478     8297        0        0       67     1856 
     280   0.15510631     7583        0        0       53     1534 
     290   0.16131616     6910        0        0       49     1410 
     300   0.17656922    16318        0        0       36     1137 
     310    0.1874423    15790        0        0       28     1006 
     320   0.19873166    15154        0        0       58     1653 
     330   0.21043777    13973        0        0      112     2475 
     340   0.22175431    12601        0        0      124     2751 
     350    0.2322731    11323        0        0      109     2578 
     360   0.24192572    10254        0        0       75     2230 
     370   0.25014925     9350        0        0       56     1907 
     380   0.25782728     8556        0        0       51     1808 
     390   0.26458049     7787        0        0       54     1512 
     400   0.28229785    17151        0        0       38     1268 
     410   0.29386258    16541        0        0       31     1126 
     420   0.30556989    15842        0        0       55     1669 
     430   0.31810427    14627        0        0      117     2409 
     440   0.33031917    13245        0        0      128     2837 
     450    0.3417604    11926        0        0       98     2655 
     460    0.3517704    10792        0        0       87     2344 
     470   0.36085153     9880        0        0       73     2042 
     480    0.3688674     9042        0        0       58     1813 
     490   0.37616944     8264        0        0       48     1470 
     500   0.39296222    17600        0        0       33     1339 
Loop time of 0.392977 on 1 procs for 500 steps with 17600 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.33089    | 0.33089    | 0.33089    |   0.0 | 84.20
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0019362  | 0.0019362  | 0.0019362  |   0.0 |  0.49
Modify  | 0.050269   | 0.050269   | 0.050269   |   0.0 | 12.79
Output  | 0.0017917  | 0.0017917  | 0.0017917  |   0.0 |  0.46
Other   |            | 0.00809    |            |       |  2.06

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

Particle-moves/CPUsec/proc: 1.08002e+07
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

