SPARTA (26 Feb 2021)
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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/ascldap/users/stamoor/sparta_stanmoore1/src/grid.cpp:410)
Created 8 child grid cells
  CPU time = 0.0014627 secs
  create/ghost percent = 70.8883 29.1117
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000996113 secs
  reassign/sort/migrate/ghost percent = 76.0172 0.909526 12.4222 10.651

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
  0 0 = number of pushed cells
  4 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  4 0 4 = cells outside/inside/overlapping surfs
  4 = surf cells with 1,2,etc splits
  4356 4356 = cell-wise and global flow volume
  CPU time = 0.0169191 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 1.51908 1.07519 0.0591849 96.525 0.821543 5.37033 0.217011
  surf2grid time = 0.0163312 secs
  map/comm1/comm2/comm3/comm4/split percent = 59.5915 0.309498 0.181027 1.73436 0.382493 37.3339

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
      10 0.00034308434        0        0        0        0        0 
      20 0.00062680244        0        0        0        0        0 
      30 0.00090026855        0        0        0        0        0 
      40 0.0011868477        0        0        0        0        0 
      50 0.0014705658        0        0        0        0        0 
      60 0.0017488003        0        0        0        0        0 
      70 0.0020489693        0        0        0        0        0 
      80 0.0023267269        0        0        0        0        0 
      90 0.0025911331        0        0        0        0        0 
     100 0.0072956085     6273        0        0        0        0 
     110  0.012465239     6273        0        0        0        0 
     120  0.015183449     6273        0        0        0        0 
     130  0.017612696     6273        0        0        0        0 
     140  0.019683838     6273        0        0        0        0 
     150  0.023407698     6273        0        0        0    49952 
     160  0.041374445     6273        0        0        0    50184 
     170  0.059551477     6273        0        0        0    50184 
     180  0.077636003     6273        0        0        0    50184 
     190    0.0995996      182        0        0     6178    50888 
     200    0.1045773     6450        0        0        0     1040 
     210   0.10700679     6450        0        0        0     1016 
     220   0.10990953     6450        0        0        0     1016 
     230   0.11256742     6450        0        0        0     1016 
     240   0.11499882     6450        0        0        0     1000 
     250    0.1193769     6450        0        0        0    51336 
     260   0.13797593     6450        0        0        0    51432 
     270   0.15430689     6449        0        0        0    51336 
     280   0.17359471     6448        0        0        0    51248 
     290    0.1973052      310        0        0     6233    51928 
     300   0.20286012     6572        0        0        0     1656 
     310   0.20617008     6569        0        0        0     1608 
     320   0.20937347     6564        0        0        0     1544 
     330   0.21465635     6556        0        0        0     1480 
     340    0.2179985     6551        0        0        0     1424 
     350    0.2220819     6545        0        0        0    51608 
     360   0.24757457     6537        0        0        0    51688 
     370   0.26873255     6529        0        0        0    51536 
     380   0.29332256     6520        0        0        0    51352 
     390   0.31388998      414        0        0     6211    52152 
     400    0.3254354     6616        0        0        0     1976 
     410   0.32910109     6604        0        0        0     1880 
     420   0.33461761     6598        0        0        0     1792 
     430   0.33916521     6587        0        0        0     1776 
     440   0.34593177     6572        0        0        0     1736 
     450   0.35433626     6563        0        0        0    51400 
     460   0.37758136     6555        0        0        0    51496 
     470     0.397048     6541        0        0        0    51328 
     480   0.41634297     6533        0        0        0    51224 
     490   0.43569064      454        0        0     6141    51800 
     500   0.44011903     6662        0        0        0     2064 
     510   0.44466233     6643        0        0        0     1968 
     520   0.44774842     6632        0        0        0     1888 
     530   0.45043564     6622        0        0        0     1824 
     540   0.45304036     6607        0        0        0     1784 
     550   0.46043015     6588        0        0        0    51632 
     560   0.48095202     6575        0        0        0    51736 
     570   0.49815392     6555        0        0        0    51552 
     580   0.51825047     6546        0        0        0    51448 
     590   0.53772807      467        0        0     6171    52136 
     600   0.54224563     6684        0        0        0     2208 
     610   0.54538941     6675        0        0        0     2032 
     620   0.55046296     6659        0        0        0     1920 
     630    0.5533638     6645        0        0        0     1840 
     640   0.55626106     6628        0        0        0     1720 
     650   0.56083584     6615        0        0        0    51688 
     660   0.57751584     6602        0        0        0    51704 
     670   0.59509873     6587        0        0        0    51536 
     680   0.61271739     6572        0        0        0    51312 
     690   0.63596964      480        0        0     6196    52064 
     700    0.6427145     6714        0        0        0     2064 
     710   0.64554167     6700        0        0        0     1976 
     720   0.64834833     6686        0        0        0     1888 
     730   0.65098262     6673        0        0        0     1840 
     740   0.65358257     6650        0        0        0     1696 
     750   0.65794182     6630        0        0        0    51664 
     760   0.67435908     6620        0        0        0    51744 
     770   0.69088531     6601        0        0        0    51600 
     780   0.70920563     6583        0        0        0    51440 
     790   0.73030782      454        0        0     6205    51992 
     800   0.73490787     6684        0        0        0     1968 
     810    0.7380352     6661        0        0        0     1888 
     820   0.74080729     6640        0        0        0     1768 
     830   0.74343634     6624        0        0        0     1696 
     840   0.74603105     6615        0        0        0     1600 
     850   0.75409889     6606        0        0        0    51504 
     860   0.77747846     6590        0        0        0    51600 
     870   0.80081034     6573        0        0        0    51496 
     880   0.81937647     6558        0        0        0    51352 
     890   0.84119201      478        0        0     6183    52152 
     900   0.84570146     6581        0        0        0     2160 
     910   0.84850454     6562        0        0        0     2088 
     920    0.8511498     6550        0        0        0     2032 
     930   0.85370612     6535        0        0        0     1904 
     940   0.85638809     6510        0        0        0     1736 
     950   0.86043596     6496        0        0        0    50656 
     960   0.87678623     6485        0        0        0    50760 
     970   0.89350128     6472        0        0        0    50600 
     980    0.9093461     6459        0        0        0    50480 
     990   0.93043208      476        0        0     6060    51096 
    1000   0.93480206     6645        0        0        0     2040 
Loop time of 0.934875 on 4 procs for 1000 steps with 6645 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.012518   | 0.37286    | 0.77766    |  59.2 | 39.88
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.020862   | 0.063227   | 0.07884    |   9.7 |  6.76
Modify  | 0.00025105 | 0.017837   | 0.038248   |  13.2 |  1.91
Output  | 0.0043612  | 0.0055885  | 0.0092356  |   2.8 |  0.60
Other   |            | 0.4754     |            |       | 50.85

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

Particle-moves/CPUsec/proc: 1.44108e+06
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
