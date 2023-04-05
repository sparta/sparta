SPARTA (18 Jul 2022)
Running on 4 MPI task(s)
################################################################################
# beam of particles striking the surface at an inclined angle
# free molecular flow (no collisions)
#
# Note:
#  - The "comm/sort” option to the “global” command is used to match MPI runs.
# The "comm/sort" option should not be used for production runs.
################################################################################

seed	    	    123456
dimension   	    3
global              gridcut 0.0 comm/sort yes

boundary	    	oo oo oo


create_box          -11 11 -11 11 0 10
Created orthogonal box = (-11 -11 0) to (11 11 10)
create_grid 	    2 2 2
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/me/sparta_master/src/grid.cpp:465)
Created 8 child grid cells
  CPU time = 0.00148976 secs
  create/ghost percent = 85.3918 14.6082
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.00066688 secs
  reassign/sort/migrate/ghost percent = 71.4439 0.400222 10.2402 17.9157

global		    	nrho 1e10 fnum 1e6

species		    	air.species O CO CO2 O2 C
mixture		    	air O O2 vstream 0 1000 -1000

mixture             air O   frac 1.0
mixture             air CO  frac 0.0
mixture             air CO2 frac 0.0
mixture             air C   frac 0.0
mixture 			air O2 	frac 0.0


surf_collide        1 cll 300.0 0.5 0.5 0.5 0.5

read_surf			base_plate.surf
  8 points
  12 triangles
  -11 11 xlo xhi
  -11 11 ylo yhi
  0 1 zlo zhi
  1 min triangle edge length
  11 min triangle area
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  4 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  4356 4356 = cell-wise and global flow volume
  CPU time = 0.000965509 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 25.6877 10.872 0.514754 52.7093 10.2163 11.5207 3.93337
  surf2grid time = 0.000508913 secs
  map/comm1/comm2/comm3/comm4/split percent = 25.6409 11.3096 6.89548 4.87942 16.8532 29.1582

##################################### SURF REACT ADSORB ######################################
##################################### SURF OPTION ############################################

#surf_react        	 adsorb_test_gs1 adsorb gs sample-GS_1.surf nsync 1 surf 1000 6.022e18 O CO
#surf_modify 		 all collide 1 react adsorb_test_gs1

surf_react        	adsorb_test_gs2 adsorb gs sample-GS_2.surf nsync 1 surf 1000 6.022e18 O CO
surf_modify 		all collide 1 react adsorb_test_gs2

########################## BEAM ############################################################
# Beam at multiple points so that different processors handle the surface collisions

region              circle2 cylinder z  6 -10 1 -INF INF
region              circle3 cylinder z -6 -10 1 -INF INF

fix                 in2 emit/face/file air zhi data.beam beam_area_2 nevery 100 region circle2
fix                 in3 emit/face/file air zhi data.beam beam_area_3 nevery 100 region circle3

################################################################################################

#dump                2 image all 10 image.*.ppm type type pdiam 0.2 surf proc 0.01 size 512 512 zoom 1.75 gline no 0.005
#dump_modify	     	2 pad 4

timestep            0.0001

stats		    	10
stats_style	    	step cpu np nattempt ncoll nscoll nscheck
run 		    	1000
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00151062 0.00151062 0.00151062
  total     (ave,min,max) = 1.5153 1.5153 1.5153
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
      10  0.000341479        0        0        0        0        0 
      20  0.000662129        0        0        0        0        0 
      30  0.000972388        0        0        0        0        0 
      40  0.001316458        0        0        0        0        0 
      50  0.001637703        0        0        0        0        0 
      60  0.001961207        0        0        0        0        0 
      70  0.002275085        0        0        0        0        0 
      80  0.002593855        0        0        0        0        0 
      90  0.002900401        0        0        0        0        0 
     100  0.006927798     6273        0        0        0        0 
     110  0.008355641     6273        0        0        0        0 
     120   0.00978672     6273        0        0        0        0 
     130  0.011194566     6273        0        0        0        0 
     140  0.012615818     6273        0        0        0        0 
     150  0.014723191     6273        0        0        0    49952 
     160  0.021543373     6273        0        0        0    50184 
     170  0.029760042     6273        0        0        0    50184 
     180  0.036582136     6273        0        0        0    50184 
     190  0.048390892      182        0        0     6178    50888 
     200  0.052475801     6450        0        0        0     1040 
     210  0.054102952     6450        0        0        0     1016 
     220  0.055758969     6450        0        0        0     1016 
     230  0.057470506     6450        0        0        0     1016 
     240  0.059054726     6450        0        0        0     1000 
     250   0.06174992     6450        0        0        0    51336 
     260  0.068708491     6450        0        0        0    51432 
     270  0.075651479     6449        0        0        0    51336 
     280  0.082739353     6448        0        0        0    51248 
     290  0.092473294      310        0        0     6233    51928 
     300  0.095113823     6572        0        0        0     1656 
     310  0.096884822     6569        0        0        0     1608 
     320  0.098646977     6564        0        0        0     1544 
     330   0.10040193     6556        0        0        0     1480 
     340   0.10214755     6551        0        0        0     1424 
     350   0.10453266     6545        0        0        0    51608 
     360    0.1116349     6537        0        0        0    51688 
     370    0.1187141     6529        0        0        0    51536 
     380   0.12580833     6520        0        0        0    51352 
     390   0.13557687      414        0        0     6211    52152 
     400   0.13830059     6616        0        0        0     1976 
     410   0.14011377     6604        0        0        0     1880 
     420   0.14189727     6598        0        0        0     1792 
     430   0.14363028     6587        0        0        0     1776 
     440   0.14541736     6572        0        0        0     1736 
     450   0.14779231     6563        0        0        0    51400 
     460   0.15483083     6555        0        0        0    51496 
     470   0.16203063     6541        0        0        0    51328 
     480   0.16907914     6533        0        0        0    51224 
     490   0.17966168      454        0        0     6141    51800 
     500   0.18412388     6662        0        0        0     2064 
     510   0.18690286     6643        0        0        0     1968 
     520   0.18872594     6632        0        0        0     1888 
     530   0.19047684     6622        0        0        0     1824 
     540   0.19262116     6607        0        0        0     1784 
     550   0.19502261     6588        0        0        0    51632 
     560    0.2021513     6575        0        0        0    51736 
     570   0.20943939     6555        0        0        0    51552 
     580   0.21653163     6546        0        0        0    51448 
     590   0.22633838      467        0        0     6171    52136 
     600   0.22909097     6684        0        0        0     2208 
     610   0.23092988     6675        0        0        0     2032 
     620   0.23273026     6659        0        0        0     1920 
     630   0.23449806     6645        0        0        0     1840 
     640   0.23619966     6628        0        0        0     1720 
     650   0.23862513     6615        0        0        0    51688 
     660   0.24578575     6602        0        0        0    51704 
     670   0.25295247     6587        0        0        0    51536 
     680     0.260044     6572        0        0        0    51312 
     690   0.27002842      480        0        0     6196    52064 
     700   0.27276521     6714        0        0        0     2064 
     710   0.27463503     6700        0        0        0     1976 
     720   0.27645981     6686        0        0        0     1888 
     730   0.27822409     6673        0        0        0     1840 
     740   0.27997002     6650        0        0        0     1696 
     750   0.28272675     6630        0        0        0    51664 
     760   0.28981203     6620        0        0        0    51744 
     770   0.29686011     6601        0        0        0    51600 
     780   0.30391683     6583        0        0        0    51440 
     790   0.31360983      454        0        0     6205    51992 
     800   0.31636363     6684        0        0        0     1968 
     810     0.318145     6661        0        0        0     1888 
     820   0.31993942     6640        0        0        0     1768 
     830   0.32168133     6624        0        0        0     1696 
     840   0.32338965     6615        0        0        0     1600 
     850   0.32574385     6606        0        0        0    51504 
     860   0.33282431     6590        0        0        0    51600 
     870    0.3399268     6573        0        0        0    51496 
     880   0.34698251     6558        0        0        0    51352 
     890    0.3567383      478        0        0     6183    52152 
     900   0.35942225     6581        0        0        0     2160 
     910   0.36122194     6562        0        0        0     2088 
     920   0.36299396     6550        0        0        0     2032 
     930    0.3648249     6535        0        0        0     1904 
     940   0.36655688     6510        0        0        0     1736 
     950   0.36893736     6496        0        0        0    50656 
     960   0.37588605     6485        0        0        0    50760 
     970   0.38293268     6472        0        0        0    50600 
     980   0.38987768     6459        0        0        0    50480 
     990   0.39939853      476        0        0     6060    51096 
    1000   0.40211044     6645        0        0        0     2040 
Loop time of 0.402219 on 4 procs for 1000 steps with 6645 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0054091  | 0.16748    | 0.32967    |  39.6 | 41.64
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.018017   | 0.019276   | 0.020765   |   0.8 |  4.79
Modify  | 0.00026536 | 0.010453   | 0.020816   |  10.0 |  2.60
Output  | 0.0030888  | 0.0045382  | 0.0088104  |   3.7 |  1.13
Other   |            | 0.2005     |            |       | 49.84

Particle moves    = 5388927 (5.39M)
Cells touched     = 5447071 (5.45M)
Particle comms    = 725 (0.725K)
Boundary collides = 0 (0K)
Boundary exits    = 941 (0.941K)
SurfColl checks   = 19760352 (19.8M)
SurfColl occurs   = 56465 (56.5K)
Surf reactions    = 56465 (56.5K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 3.3495e+06
Particle-moves/step: 5388.93
Cell-touches/particle/step: 1.01079
Particle comm iterations/step: 1.456
Particle fraction communicated: 0.000134535
Particle fraction colliding with boundary: 0
Particle fraction exiting boundary: 0.000174617
Surface-checks/particle/step: 3.66684
Surface-collisions/particle/step: 0.010478
Surf-reactions/particle/step: 0.010478
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 56465
    reaction O(g) --> O(s): 41989
    reaction O(g) + O(s) --> CO2(g): 9
    reaction O(g) --> CO(s): 13164
    reaction O(g) --> CO(g): 1275
    reaction O(g) + O(s) --> O(g) + O(g): 28

Particles: 1661.25 ave 3277 max 54 min
Histogram: 2 0 0 0 0 0 0 0 0 2
Cells:      2 ave 2 max 2 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 6 ave 6 max 6 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 6 ave 6 max 6 min
Histogram: 4 0 0 0 0 0 0 0 0 0
Surfs:    12 ave 12 max 12 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0
