SPARTA (18 Jul 2022)
Running on 1 MPI task(s)
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
  CPU time = 0.00125111 secs
  create/ghost percent = 97.6776 2.32241
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000134763 secs
  reassign/sort/migrate/ghost percent = 76.0513 0.437064 16.4556 7.05609

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
  CPU time = 0.00100466 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 17.1788 25.8917 0.19121 48.118 8.62026 9.92677 0.0217985
  surf2grid time = 0.000483421 secs
  map/comm1/comm2/comm3/comm4/split percent = 26.8966 6.26824 3.44793 3.75884 14.1343 44.364

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
      10   3.8211e-05        0        0        0        0        0 
      20  0.000108621        0        0        0        0        0 
      30  0.000175586        0        0        0        0        0 
      40  0.000270296        0        0        0        0        0 
      50  0.000340499        0        0        0        0        0 
      60   0.00040642        0        0        0        0        0 
      70  0.000472314        0        0        0        0        0 
      80  0.000537823        0        0        0        0        0 
      90  0.000628574        0        0        0        0        0 
     100  0.005473618     6218        0        0        0        0 
     110  0.007725902     6218        0        0        0        0 
     120   0.00998318     6218        0        0        0        0 
     130  0.012121103     6218        0        0        0        0 
     140  0.014265231     6218        0        0        0        0 
     150  0.017580437     6218        0        0        0    49504 
     160  0.029184001     6218        0        0        0    49744 
     170  0.040822483     6218        0        0        0    49744 
     180   0.05243086     6218        0        0        0    49744 
     190  0.068557845      189        0        0     6134    50584 
     200  0.072513032     6432        0        0        0     1088 
     210  0.075006637     6432        0        0        0     1072 
     220  0.077506148     6432        0        0        0     1072 
     230  0.079998859     6432        0        0        0     1064 
     240   0.08249419     6432        0        0        0     1056 
     250  0.086125712     6432        0        0        0    51168 
     260  0.098206494     6432        0        0        0    51320 
     270   0.11041618     6430        0        0        0    51224 
     280   0.12253045     6425        0        0        0    51112 
     290   0.13931996      326        0        0     6211    51920 
     300   0.14365981     6607        0        0        0     1736 
     310   0.14630496     6604        0        0        0     1640 
     320   0.14900123     6601        0        0        0     1608 
     330   0.15162798     6597        0        0        0     1584 
     340   0.15431283     6589        0        0        0     1568 
     350   0.15816154     6588        0        0        0    51824 
     360   0.17048968     6579        0        0        0    51944 
     370   0.18285611     6573        0        0        0    51816 
     380   0.19506155     6563        0        0        0    51680 
     390   0.21205155      423        0        0     6240    52480 
     400   0.21640077     6636        0        0        0     2104 
     410   0.21907456     6622        0        0        0     1984 
     420   0.22175276     6607        0        0        0     1920 
     430   0.22439271     6592        0        0        0     1840 
     440   0.22672515     6580        0        0        0     1768 
     450   0.22917565     6563        0        0        0    51496 
     460   0.23594408     6548        0        0        0    51536 
     470   0.24274761     6535        0        0        0    51392 
     480   0.24956924     6526        0        0        0    51232 
     490   0.25942388      447        0        0     6155    51728 
     500   0.26216235     6716        0        0        0     1944 
     510   0.26379208     6705        0        0        0     1824 
     520   0.26541163     6690        0        0        0     1744 
     530    0.2669506     6679        0        0        0     1672 
     540   0.26864756     6665        0        0        0     1560 
     550   0.27209086     6652        0        0        0    51832 
     560   0.28382484     6638        0        0        0    51992 
     570   0.29546466     6625        0        0        0    51864 
     580   0.30765943     6609        0        0        0    51704 
     590   0.32461174      469        0        0     6218    52312 
     600   0.32887312     6764        0        0        0     1928 
     610   0.33170604     6747        0        0        0     1840 
     620   0.33437656     6739        0        0        0     1744 
     630   0.33704272     6722        0        0        0     1704 
     640   0.33964537     6706        0        0        0     1600 
     650   0.34343814     6696        0        0        0    52296 
     660    0.3556407     6689        0        0        0    52336 
     670   0.36781593     6671        0        0        0    52176 
     680   0.37981482     6657        0        0        0    52016 
     690   0.39578681      472        0        0     6271    52688 
     700   0.39847457     6701        0        0        0     1920 
     710    0.4000669     6685        0        0        0     1872 
     720   0.40168857     6667        0        0        0     1760 
     730   0.40323773     6649        0        0        0     1648 
     740   0.40480689     6630        0        0        0     1552 
     750    0.4069434     6622        0        0        0    51616 
     760   0.41484315     6607        0        0        0    51752 
     770   0.42669739     6593        0        0        0    51584 
     780   0.43854717     6582        0        0        0    51496 
     790   0.45514465      450        0        0     6205    52120 
     800   0.45913912     6628        0        0        0     1928 
     810   0.46168116     6609        0        0        0     1832 
     820   0.46414786     6598        0        0        0     1744 
     830   0.46661522     6579        0        0        0     1640 
     840   0.46906492     6568        0        0        0     1528 
     850   0.47264948     6549        0        0        0    51080 
     860   0.48446489     6530        0        0        0    51160 
     870   0.49616135     6511        0        0        0    51016 
     880   0.50780709     6499        0        0        0    50904 
     890   0.52420566      449        0        0     6140    51632 
     900   0.52893405     6634        0        0        0     1920 
     910   0.53144974     6621        0        0        0     1824 
     920   0.53407677     6603        0        0        0     1720 
     930   0.53654579     6587        0        0        0     1640 
     940   0.53896354     6571        0        0        0     1544 
     950   0.54254519     6556        0        0        0    51280 
     960   0.55426425     6541        0        0        0    51384 
     970   0.56602272     6529        0        0        0    51264 
     980   0.57812761     6519        0        0        0    51136 
     990   0.59484375      470        0        0     6165    52048 
    1000   0.59922549     6662        0        0        0     2280 
Loop time of 0.599243 on 1 procs for 1000 steps with 6662 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.55422    | 0.55422    | 0.55422    |   0.0 | 92.49
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0016485  | 0.0016485  | 0.0016485  |   0.0 |  0.28
Modify  | 0.03263    | 0.03263    | 0.03263    |   0.0 |  5.45
Output  | 0.0073114  | 0.0073114  | 0.0073114  |   0.0 |  1.22
Other   |            | 0.003438   |            |       |  0.57

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

Particle-moves/CPUsec/proc: 9.00776e+06
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
