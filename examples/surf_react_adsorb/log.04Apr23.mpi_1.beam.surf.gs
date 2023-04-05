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
  CPU time = 0.000994682 secs
  create/ghost percent = 88.6625 11.3375
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000126123 secs
  reassign/sort/migrate/ghost percent = 64.0832 0.189036 27.9773 7.75047

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
  CPU time = 0.000626087 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 24.6002 8.07312 1.02818 56.1691 10.1295 5.40746 0.380807
  surf2grid time = 0.000351667 secs
  map/comm1/comm2/comm3/comm4/split percent = 31.7966 8.0678 3.38983 3.38983 9.49153 40.1356

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
      10 0.00011396408        0        0        0        0        0 
      20 0.00016021729        0        0        0        0        0 
      30 0.00020909309        0        0        0        0        0 
      40 0.00025510788        0        0        0        0        0 
      50 0.00030088425        0        0        0        0        0 
      60 0.00034666061        0        0        0        0        0 
      70 0.00039219856        0        0        0        0        0 
      80 0.00043797493        0        0        0        0        0 
      90 0.0005671978        0        0        0        0        0 
     100 0.0056979656     6218        0        0        0        0 
     110 0.0076072216     6218        0        0        0        0 
     120  0.009447813     6218        0        0        0        0 
     130  0.011292934     6218        0        0        0        0 
     140  0.013139725     6218        0        0        0        0 
     150  0.015848398     6218        0        0        0    49504 
     160  0.024423599     6218        0        0        0    49744 
     170  0.032932997     6218        0        0        0    49744 
     180  0.041404009     6218        0        0        0    49744 
     190  0.052113533      189        0        0     6134    50584 
     200  0.055992365     6432        0        0        0     1088 
     210  0.058073759     6432        0        0        0     1072 
     220  0.060152531     6432        0        0        0     1072 
     230  0.062231064     6432        0        0        0     1064 
     240  0.064319372     6432        0        0        0     1056 
     250  0.067189932     6432        0        0        0    51168 
     260   0.07594943     6432        0        0        0    51320 
     270  0.084700108     6430        0        0        0    51224 
     280  0.093493223     6425        0        0        0    51112 
     290   0.10427952      326        0        0     6211    51920 
     300   0.10832667     6607        0        0        0     1736 
     310   0.11054707     6604        0        0        0     1640 
     320   0.11275601     6601        0        0        0     1608 
     330   0.11495733     6597        0        0        0     1584 
     340   0.11716199     6589        0        0        0     1568 
     350   0.12015271     6588        0        0        0    51824 
     360   0.12905097     6579        0        0        0    51944 
     370   0.13790822     6573        0        0        0    51816 
     380   0.14715481     6563        0        0        0    51680 
     390   0.15820765      423        0        0     6240    52480 
     400   0.16231132     6636        0        0        0     2104 
     410    0.1646111     6622        0        0        0     1984 
     420   0.16687489     6607        0        0        0     1920 
     430   0.16912293     6592        0        0        0     1840 
     440   0.17135978     6580        0        0        0     1768 
     450   0.17436981     6563        0        0        0    51496 
     460   0.18319321     6548        0        0        0    51536 
     470   0.19204926     6535        0        0        0    51392 
     480   0.20113802     6526        0        0        0    51232 
     490   0.21183324      447        0        0     6155    51728 
     500   0.21599269     6716        0        0        0     1944 
     510   0.21827674     6705        0        0        0     1824 
     520   0.22055101     6690        0        0        0     1744 
     530   0.22280622     6679        0        0        0     1672 
     540   0.22504544     6665        0        0        0     1560 
     550   0.22806883     6652        0        0        0    51832 
     560   0.23699975     6638        0        0        0    51992 
     570   0.24596834     6625        0        0        0    51864 
     580   0.25497127     6609        0        0        0    51704 
     590   0.26580763      469        0        0     6218    52312 
     600   0.26995063     6764        0        0        0     1928 
     610   0.27224731     6747        0        0        0     1840 
     620    0.2745347     6739        0        0        0     1744 
     630   0.27680397     6722        0        0        0     1704 
     640   0.27905941     6706        0        0        0     1600 
     650   0.28210807     6696        0        0        0    52296 
     660   0.29109764     6689        0        0        0    52336 
     670    0.3001771     6671        0        0        0    52176 
     680   0.30915523     6657        0        0        0    52016 
     690   0.32007694      472        0        0     6271    52688 
     700   0.32426262     6701        0        0        0     1920 
     710   0.32655263     6685        0        0        0     1872 
     720   0.32882023     6667        0        0        0     1760 
     730    0.3310709     6649        0        0        0     1648 
     740    0.3333025     6630        0        0        0     1552 
     750   0.33631182     6622        0        0        0    51616 
     760   0.34526157     6607        0        0        0    51752 
     770   0.35433626     6593        0        0        0    51584 
     780   0.36319518     6582        0        0        0    51496 
     790   0.37394357      450        0        0     6205    52120 
     800   0.37795472     6628        0        0        0     1928 
     810   0.38021779     6609        0        0        0     1832 
     820   0.38246799     6598        0        0        0     1744 
     830   0.38467741     6579        0        0        0     1640 
     840   0.38687181     6568        0        0        0     1528 
     850   0.38985872     6549        0        0        0    51080 
     860   0.39884114     6530        0        0        0    51160 
     870   0.40766573     6511        0        0        0    51016 
     880   0.41640115     6499        0        0        0    50904 
     890   0.42706394      449        0        0     6140    51632 
     900   0.43111539     6634        0        0        0     1920 
     910   0.43337584     6621        0        0        0     1824 
     920   0.43562341     6603        0        0        0     1720 
     930   0.43782711     6587        0        0        0     1640 
     940   0.44003868     6571        0        0        0     1544 
     950   0.44305992     6556        0        0        0    51280 
     960    0.4519887     6541        0        0        0    51384 
     970    0.4607892     6529        0        0        0    51264 
     980   0.46958184     6519        0        0        0    51136 
     990   0.48028708      470        0        0     6165    52048 
    1000   0.48438287     6662        0        0        0     2280 
Loop time of 0.484394 on 1 procs for 1000 steps with 6662 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.43944    | 0.43944    | 0.43944    |   0.0 | 90.72
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0011408  | 0.0011408  | 0.0011408  |   0.0 |  0.24
Modify  | 0.035561   | 0.035561   | 0.035561   |   0.0 |  7.34
Output  | 0.0019951  | 0.0019951  | 0.0019951  |   0.0 |  0.41
Other   |            | 0.006257   |            |       |  1.29

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

Particle-moves/CPUsec/proc: 1.11435e+07
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
