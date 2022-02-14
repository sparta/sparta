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
Created 8 child grid cells
  CPU time = 0.00101829 secs
  create/ghost percent = 89.0658 10.9342
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.00014782 secs
  reassign/sort/migrate/ghost percent = 57.4194 1.6129 28.871 12.0968

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
  CPU time = 0.00138736 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 17.7866 6.99433 0.463997 61.9522 12.8029 2.57776 0.0171851
  surf2grid time = 0.000859499 secs
  map/comm1/comm2/comm3/comm4/split percent = 23.9112 3.93897 1.58114 1.96949 4.77115 60.6103

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
      10 0.00011825562        0        0        0        0        0 
      20 0.00018000603        0        0        0        0        0 
      30 0.00023674965        0        0        0        0        0 
      40 0.00029301643        0        0        0        0        0 
      50 0.00037145615        0        0        0        0        0 
      60 0.00044155121        0        0        0        0        0 
      70 0.00050091743        0        0        0        0        0 
      80 0.00055861473        0        0        0        0        0 
      90 0.0006377697        0        0        0        0        0 
     100 0.0084490776     6218        0        0        0        0 
     110  0.011951685     6218        0        0        0        0 
     120  0.015459299     6218        0        0        0        0 
     130  0.018955946     6218        0        0        0        0 
     140  0.022462368     6218        0        0        0        0 
     150  0.028987885     6218        0        0        0    49504 
     160  0.059415579     6218        0        0        0    49744 
     170  0.089929342     6218        0        0        0    49744 
     180   0.12062478     6218        0        0        0    49744 
     190   0.15625453      189        0        0     6134    50584 
     200   0.16393614     6432        0        0        0     1088 
     210   0.16814089     6432        0        0        0     1072 
     220   0.17233753     6432        0        0        0     1072 
     230   0.17666149     6432        0        0        0     1064 
     240   0.18086314     6432        0        0        0     1056 
     250   0.18867087     6432        0        0        0    51168 
     260     0.220577     6432        0        0        0    51320 
     270   0.25264335     6430        0        0        0    51224 
     280   0.28359222     6425        0        0        0    51112 
     290   0.32187915      326        0        0     6211    51920 
     300   0.33041406     6607        0        0        0     1736 
     310   0.33560777     6604        0        0        0     1640 
     320   0.34031773     6601        0        0        0     1608 
     330   0.34535766     6597        0        0        0     1584 
     340   0.35049915     6589        0        0        0     1568 
     350   0.35838056     6588        0        0        0    51824 
     360   0.39139032     6579        0        0        0    51944 
     370   0.42327476     6573        0        0        0    51816 
     380   0.45438004     6563        0        0        0    51680 
     390    0.4904778      423        0        0     6240    52480 
     400   0.49861932     6636        0        0        0     2104 
     410   0.50347662     6622        0        0        0     1984 
     420   0.50827265     6607        0        0        0     1920 
     430   0.51302838     6592        0        0        0     1840 
     440   0.51778817     6580        0        0        0     1768 
     450   0.52548099     6563        0        0        0    51496 
     460   0.55649161     6548        0        0        0    51536 
     470   0.58742476     6535        0        0        0    51392 
     480   0.61827326     6526        0        0        0    51232 
     490   0.65475273      447        0        0     6155    51728 
     500    0.6630106     6716        0        0        0     1944 
     510   0.66783381     6705        0        0        0     1824 
     520   0.67259645     6690        0        0        0     1744 
     530   0.67730021     6679        0        0        0     1672 
     540   0.68195915     6665        0        0        0     1560 
     550   0.68961954     6652        0        0        0    51832 
     560    0.7209208     6638        0        0        0    51992 
     570   0.75214839     6625        0        0        0    51864 
     580   0.78330469     6609        0        0        0    51704 
     590   0.81935334      469        0        0     6218    52312 
     600   0.82736349     6764        0        0        0     1928 
     610   0.83213305     6747        0        0        0     1840 
     620   0.83675528     6739        0        0        0     1744 
     630    0.8413229     6722        0        0        0     1704 
     640    0.8458817     6706        0        0        0     1600 
     650   0.85339427     6696        0        0        0    52296 
     660    0.8835783     6689        0        0        0    52336 
     670   0.91416383     6671        0        0        0    52176 
     680   0.94439125     6657        0        0        0    52016 
     690   0.97975016      472        0        0     6271    52688 
     700   0.98751116     6701        0        0        0     1920 
     710   0.99223042     6685        0        0        0     1872 
     720   0.99695182     6667        0        0        0     1760 
     730     1.001585     6649        0        0        0     1648 
     740    1.0060828     6630        0        0        0     1552 
     750    1.0135169     6622        0        0        0    51616 
     760    1.0435734     6607        0        0        0    51752 
     770     1.073312     6593        0        0        0    51584 
     780    1.1032517     6582        0        0        0    51496 
     790    1.1378229      450        0        0     6205    52120 
     800    1.1456504     6628        0        0        0     1928 
     810    1.1502678     6609        0        0        0     1832 
     820    1.1547573     6598        0        0        0     1744 
     830    1.1591773     6579        0        0        0     1640 
     840    1.1635509     6568        0        0        0     1528 
     850    1.1707211     6549        0        0        0    51080 
     860    1.2006075     6530        0        0        0    51160 
     870    1.2301877     6511        0        0        0    51016 
     880    1.2596879     6499        0        0        0    50904 
     890    1.2940338      449        0        0     6140    51632 
     900    1.3016789     6634        0        0        0     1920 
     910    1.3062468     6621        0        0        0     1824 
     920    1.3107178     6603        0        0        0     1720 
     930     1.315134     6587        0        0        0     1640 
     940    1.3196354     6571        0        0        0     1544 
     950    1.3271105     6556        0        0        0    51280 
     960     1.357326     6541        0        0        0    51384 
     970    1.3875732     6529        0        0        0    51264 
     980    1.4174125     6519        0        0        0    51136 
     990    1.4520173      470        0        0     6165    52048 
    1000    1.4600065     6662        0        0        0     2280 
Loop time of 1.46002 on 1 procs for 1000 steps with 6662 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.3833     | 1.3833     | 1.3833     |   0.0 | 94.75
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0020394  | 0.0020394  | 0.0020394  |   0.0 |  0.14
Modify  | 0.063816   | 0.063816   | 0.063816   |   0.0 |  4.37
Output  | 0.0024953  | 0.0024953  | 0.0024953  |   0.0 |  0.17
Other   |            | 0.008328   |            |       |  0.57

Particle moves    = 5397837 (5.4M)
Cells touched     = 5456087 (5.46M)
Particle comms    = 0 (0K)
Boundary collides = 0 (0K)
Boundary exits    = 938 (0.938K)
SurfColl checks   = 19767352 (19.8M)
SurfColl occurs   = 56634 (56.6K)
Surf reactions    = 56634 (56.6K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 3.6971e+06
Particle-moves/step: 5397.84
Cell-touches/particle/step: 1.01079
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0
Particle fraction exiting boundary: 0.000173773
Surface-checks/particle/step: 3.66209
Surface-collisions/particle/step: 0.010492
Surf-reactions/particle/step: 0.010492
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 56634
    reaction O(g) --> O(s): 42304
    reaction O(g) + O(s) --> CO2(g): 8
    reaction O(g) --> CO(s): 13008
    reaction O(g) --> CO(g): 1292
    reaction O(g) + O(s) --> O(g) + O(g): 22

Particles: 6662 ave 6662 max 6662 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      8 ave 8 max 8 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Surfs:    12 ave 12 max 12 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
