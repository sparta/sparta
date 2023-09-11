SPARTA (13 Apr 2023)
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
  CPU time = 0.000805202 secs
  create/ghost percent = 97.69 2.30998
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.0001155 secs
  reassign/sort/migrate/ghost percent = 83.2035 0.4329 11.9481 4.41558

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
  CPU time = 0.000810501 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 15.7557 10.2776 0.111042 68.3282 5.52745 6.16915 0.0246761
  surf2grid time = 0.000553801 secs
  map/comm1/comm2/comm3/comm4/split percent = 14.3012 2.61827 1.76959 1.37233 5.03791 74.2326

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
      10     1.35e-05        0        0        0        0        0 
      20     3.69e-05        0        0        0        0        0 
      30     5.95e-05        0        0        0        0        0 
      40     8.19e-05        0        0        0        0        0 
      50    0.0001042        0        0        0        0        0 
      60    0.0001268        0        0        0        0        0 
      70    0.0001492        0        0        0        0        0 
      80    0.0001717        0        0        0        0        0 
      90    0.0002046        0        0        0        0        0 
     100  0.002523905     6218        0        0        0        0 
     110  0.005990912     6218        0        0        0        0 
     120  0.006741314     6218        0        0        0        0 
     130  0.010045721     6218        0        0        0        0 
     140  0.010795622     6218        0        0        0        0 
     150  0.011928624     6218        0        0        0    49504 
     160  0.017191236     6218        0        0        0    49744 
     170  0.029204661     6218        0        0        0    49744 
     180  0.033228969     6218        0        0        0    49744 
     190  0.047002598      189        0        0     6134    50584 
     200  0.054774514     6432        0        0        0     1088 
     210  0.055652416     6432        0        0        0     1072 
     220  0.056486018     6432        0        0        0     1072 
     230   0.05733512     6432        0        0        0     1064 
     240  0.058170921     6432        0        0        0     1056 
     250  0.059418424     6432        0        0        0    51168 
     260  0.069668845     6432        0        0        0    51320 
     270  0.073837054     6430        0        0        0    51224 
     280  0.085346578     6425        0        0        0    51112 
     290   0.09121159      326        0        0     6211    51920 
     300   0.10286171     6607        0        0        0     1736 
     310   0.10380942     6604        0        0        0     1640 
     320   0.10471062     6601        0        0        0     1608 
     330   0.10560922     6597        0        0        0     1584 
     340   0.10650222     6589        0        0        0     1568 
     350   0.11454564     6588        0        0        0    51824 
     360   0.11880265     6579        0        0        0    51944 
     370   0.12299696     6573        0        0        0    51816 
     380   0.13342078     6563        0        0        0    51680 
     390   0.14714221      423        0        0     6240    52480 
     400   0.15486322     6636        0        0        0     2104 
     410   0.15583613     6622        0        0        0     1984 
     420   0.15676283     6607        0        0        0     1920 
     430   0.15770493     6592        0        0        0     1840 
     440   0.15862113     6580        0        0        0     1768 
     450   0.15993013     6563        0        0        0    51496 
     460   0.17346536     6548        0        0        0    51536 
     470   0.17778257     6535        0        0        0    51392 
     480   0.18192478     6526        0        0        0    51232 
     490   0.19508461      447        0        0     6155    51728 
     500   0.20115492     6716        0        0        0     1944 
     510   0.20213502     6705        0        0        0     1824 
     520   0.21017584     6690        0        0        0     1744 
     530   0.21113124     6679        0        0        0     1672 
     540   0.21204864     6665        0        0        0     1560 
     550   0.21336145     6652        0        0        0    51832 
     560   0.22558837     6638        0        0        0    51992 
     570   0.23039828     6625        0        0        0    51864 
     580   0.23459479     6609        0        0        0    51704 
     590   0.24720172      469        0        0     6218    52312 
     600   0.25888184     6764        0        0        0     1928 
     610   0.25986814     6747        0        0        0     1840 
     620   0.26080475     6739        0        0        0     1744 
     630   0.26174495     6722        0        0        0     1704 
     640   0.26267275     6706        0        0        0     1600 
     650   0.26398915     6696        0        0        0    52296 
     660   0.26827436     6689        0        0        0    52336 
     670   0.28147199     6671        0        0        0    52176 
     680    0.2857703     6657        0        0        0    52016 
     690   0.29955713      472        0        0     6271    52688 
     700   0.30692404     6701        0        0        0     1920 
     710   0.30790284     6685        0        0        0     1872 
     720   0.30883435     6667        0        0        0     1760 
     730   0.30976765     6649        0        0        0     1648 
     740   0.31068545     6630        0        0        0     1552 
     750   0.31856367     6622        0        0        0    51616 
     760   0.32278258     6607        0        0        0    51752 
     770   0.32697198     6593        0        0        0    51584 
     780   0.33746901     6582        0        0        0    51496 
     790   0.35107853      450        0        0     6205    52120 
     800   0.35882325     6628        0        0        0     1928 
     810   0.35977655     6609        0        0        0     1832 
     820   0.36069375     6598        0        0        0     1744 
     830   0.36162096     6579        0        0        0     1640 
     840   0.36253066     6568        0        0        0     1528 
     850   0.36382956     6549        0        0        0    51080 
     860   0.37737839     6530        0        0        0    51160 
     870    0.3815309     6511        0        0        0    51016 
     880   0.38565481     6499        0        0        0    50904 
     890   0.39902733      449        0        0     6140    51632 
     900   0.40682285     6634        0        0        0     1920 
     910   0.40778605     6621        0        0        0     1824 
     920   0.40870406     6603        0        0        0     1720 
     930   0.40962256     6587        0        0        0     1640 
     940   0.41815147     6571        0        0        0     1544 
     950   0.41949478     6556        0        0        0    51280 
     960   0.42370879     6541        0        0        0    51384 
     970    0.4278607     6529        0        0        0    51264 
     980   0.44137262     6519        0        0        0    51136 
     990   0.44724114      470        0        0     6165    52048 
    1000   0.45886616     6662        0        0        0     2280 
Loop time of 0.458925 on 1 procs for 1000 steps with 6662 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.21823    | 0.21823    | 0.21823    |   0.0 | 47.55
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0003319  | 0.0003319  | 0.0003319  |   0.0 |  0.07
Modify  | 0.014124   | 0.014124   | 0.014124   |   0.0 |  3.08
Output  | 0.22549    | 0.22549    | 0.22549    |   0.0 | 49.13
Other   |            | 0.000748   |            |       |  0.16

Particle moves    = 5397859 (5.4M)
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

Particle-moves/CPUsec/proc: 1.1762e+07
Particle-moves/step: 5397.86
Cell-touches/particle/step: 1.01079
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0
Particle fraction exiting boundary: 0.000173773
Surface-checks/particle/step: 3.66207
Surface-collisions/particle/step: 0.0104919
Surf-reactions/particle/step: 0.0104919
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
