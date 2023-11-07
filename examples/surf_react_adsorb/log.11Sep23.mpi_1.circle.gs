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

boundary	    	o r p

create_box  	    0 10 0 10 -0.5 0.5
Created orthogonal box = (0 0 -0.5) to (10 10 0.5)
create_grid 	    20 20 1
Created 400 child grid cells
  CPU time = 0.000898502 secs
  create/ghost percent = 91.4191 8.58095
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.0001409 secs
  reassign/sort/migrate/ghost percent = 61.7459 0.567779 13.0589 24.6274

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
  CPU time = 0.000486301 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 22.4554 7.62902 1.29549 59.9217 8.69832 9.08902 0.0205634
  surf2grid time = 0.0002914 secs
  map/comm1/comm2/comm3/comm4/split percent = 44.1318 8.06452 14.9279 3.32876 13.7268 12.5257

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
      10     1.69e-05        0        0        0        0        0 
      20     4.41e-05        0        0        0        0        0 
      30     7.04e-05        0        0        0        0        0 
      40     9.69e-05        0        0        0        0        0 
      50    0.0001234        0        0        0        0        0 
      60  0.000149701        0        0        0        0        0 
      70  0.000175601        0        0        0        0        0 
      80  0.000202101        0        0        0        0        0 
      90  0.000228201        0        0        0        0        0 
     100  0.015105535    10000        0        0        0        0 
     110  0.016563338    10000        0        0        0        6 
     120  0.018249142     9877        0        0       36      651 
     130  0.020231347     9145        0        0       88     1558 
     140  0.027175663     8137        0        0       95     1991 
     150  0.029040167     7226        0        0       77     1949 
     160  0.030678871     6468        0        0       81     1668 
     170  0.038259188     5893        0        0       56     1439 
     180  0.039531191     5410        0        0       43     1138 
     190  0.040652394     5016        0        0       35     1026 
     200  0.044336202    14590        0        0       26      891 
     210  0.052572621    14169        0        0       16      716 
     220  0.059547537    13675        0        0       60     1354 
     230  0.067862057    12583        0        0      119     2190 
     240   0.07376037    11304        0        0      107     2520 
     250  0.079434783    10092        0        0      105     2370 
     260  0.084047894     9153        0        0       57     1994 
     270  0.089284806     8297        0        0       67     1856 
     280  0.093053115     7583        0        0       53     1534 
     290  0.096698623     6910        0        0       49     1410 
     300   0.10546264    16318        0        0       36     1137 
     310   0.11120336    15790        0        0       28     1006 
     320   0.11744797    15154        0        0       58     1653 
     330   0.12406609    13973        0        0      112     2475 
     340    0.1306086    12601        0        0      124     2751 
     350   0.13635471    11323        0        0      109     2578 
     360   0.14635084    10254        0        0       75     2230 
     370   0.14850994     9350        0        0       56     1907 
     380   0.15039245     8556        0        0       51     1808 
     390   0.15209825     7787        0        0       54     1512 
     400   0.16514148    17151        0        0       38     1268 
     410   0.16795249    16541        0        0       31     1126 
     420   0.17083289    15842        0        0       55     1669 
     430   0.18192552    14627        0        0      117     2409 
     440   0.18493073    13245        0        0      128     2837 
     450   0.18772573    11926        0        0       98     2655 
     460   0.19820696    10792        0        0       87     2344 
     470   0.20040526     9880        0        0       73     2042 
     480   0.20235577     9042        0        0       58     1813 
     490   0.20410447     8264        0        0       48     1470 
     500    0.2163328    17600        0        0       33     1339 
Loop time of 0.216349 on 1 procs for 500 steps with 17600 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.16697    | 0.16697    | 0.16697    |   0.0 | 77.18
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0007462  | 0.0007462  | 0.0007462  |   0.0 |  0.34
Modify  | 0.028653   | 0.028653   | 0.028653   |   0.0 | 13.24
Output  | 0.018319   | 0.018319   | 0.018319   |   0.0 |  8.47
Other   |            | 0.001657   |            |       |  0.77

Particle moves    = 4244566 (4.24M)
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

Particle-moves/CPUsec/proc: 1.96191e+07
Particle-moves/step: 8489.13
Cell-touches/particle/step: 1.13023
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00382583
Particle fraction exiting boundary: 0.00169558
Surface-checks/particle/step: 0.160374
Surface-collisions/particle/step: 0.00628074
Surf-reactions/particle/step: 0.00626401
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

