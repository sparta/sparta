SPARTA (24 Sep 2025)
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
  CPU time = 0.00116704 secs
  create/ghost percent = 97.3621 2.63788
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 9.2887e-05 secs
  reassign/sort/migrate/ghost percent = 82.1568 0.573815 12.0071 5.26231

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
  CPU time = 0.000625412 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 8.40614 28.1674 0.302041 56.8513 6.27314 7.26097 0.022705
  surf2grid time = 0.000355555 secs
  map/comm1/comm2/comm3/comm4/split percent = 29.2371 6.48845 2.71463 2.8575 20.2843 37.5377

##################################### SURF REACT ADSORB ######################################
##################################### SURF OPTION ############################################

#surf_react        	 adsorb_test_gs1 adsorb gs sample-GS_1.surf nsync 1 surf 1000 6.022e18 O CO
#surf_modify 		 all collide 1 react adsorb_test_gs1

surf_react        	adsorb_test_gs2 adsorb gs sample-GS_2.surf nsync 1 surf 1000 6.022e18 O CO
surf_modify 		all collide 1 react adsorb_test_gs2

########################## BEAM ############################################################
# Beam at multiple points so that different processors handle the surface collisions

region              circle2 cylinder z  6 -10 1 INF INF
region              circle3 cylinder z -6 -10 1 INF INF

fix                 in2 emit/face/file air zhi data.beam beam_area_2 nevery 100 region circle2 twopass
fix                 in3 emit/face/file air zhi data.beam beam_area_3 nevery 100 region circle3 twopass

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
  modify    (ave,min,max) = 0 0 0
  total     (ave,min,max) = 1.5153 1.5153 1.5153
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
      10   1.7374e-05        0        0        0        0        0 
      20   4.8229e-05        0        0        0        0        0 
      30   7.8491e-05        0        0        0        0        0 
      40  0.000107757        0        0        0        0        0 
      50  0.000139655        0        0        0        0        0 
      60  0.000167327        0        0        0        0        0 
      70   0.00019857        0        0        0        0        0 
      80  0.000256784        0        0        0        0        0 
      90  0.000290148        0        0        0        0        0 
     100   0.00268705     6216        0        0        0        0 
     110  0.003382528     6216        0        0        0        0 
     120  0.004088701     6216        0        0        0        0 
     130  0.004918386     6216        0        0        0        0 
     140  0.005739727     6216        0        0        0        0 
     150  0.006893626     6216        0        0        0    49488 
     160   0.01084866     6216        0        0        0    49728 
     170  0.014819581     6216        0        0        0    49728 
     180  0.018812655     6216        0        0        0    49728 
     190  0.024388528      188        0        0     6133    50568 
     200  0.025728792     6432        0        0        0     1080 
     210  0.026627987     6432        0        0        0     1072 
     220   0.02749882     6432        0        0        0     1080 
     230  0.028265944     6432        0        0        0     1064 
     240  0.029052449     6432        0        0        0     1040 
     250  0.030226546     6432        0        0        0    51184 
     260  0.034200726     6432        0        0        0    51312 
     270  0.038236526     6427        0        0        0    51200 
     280  0.042222287     6423        0        0        0    51088 
     290  0.049496015      324        0        0     6212    51904 
     300  0.051482382     6606        0        0        0     1744 
     310  0.052668206     6604        0        0        0     1656 
     320  0.053815506     6599        0        0        0     1608 
     330  0.055019291     6595        0        0        0     1584 
     340  0.056172131     6587        0        0        0     1568 
     350  0.057879135     6584        0        0        0    51816 
     360  0.063494161     6575        0        0        0    51928 
     370  0.069046138     6569        0        0        0    51800 
     380   0.07467771     6560        0        0        0    51648 
     390  0.082242796      420        0        0     6239    52464 
     400  0.084199384     6633        0        0        0     2096 
     410   0.08541291     6618        0        0        0     1992 
     420  0.086632961     6604        0        0        0     1928 
     430  0.087786242     6589        0        0        0     1840 
     440  0.088933717     6576        0        0        0     1776 
     450  0.090660985     6560        0        0        0    51488 
     460  0.096172041     6546        0        0        0    51544 
     470   0.10197149     6535        0        0        0    51400 
     480    0.1072203     6525        0        0        0    51224 
     490   0.11327255      445        0        0     6154    51728 
     500   0.11482124     6716        0        0        0     1920 
     510   0.11582284     6700        0        0        0     1808 
     520   0.11674046     6692        0        0        0     1720 
     530   0.11763488     6682        0        0        0     1688 
     540   0.11855342     6671        0        0        0     1576 
     550   0.11991309     6657        0        0        0    51856 
     560   0.12428431     6643        0        0        0    52000 
     570   0.12863988     6628        0        0        0    51864 
     580   0.13285456     6608        0        0        0    51712 
     590   0.13854661      472        0        0     6218    52312 
     600   0.14007937     6769        0        0        0     1936 
     610   0.14105059     6753        0        0        0     1856 
     620   0.14192519     6745        0        0        0     1776 
     630   0.14282802     6727        0        0        0     1680 
     640   0.14368836     6712        0        0        0     1576 
     650   0.14498568     6697        0        0        0    52304 
     660   0.14904443     6691        0        0        0    52376 
     670   0.15324774     6676        0        0        0    52208 
     680   0.15739822     6659        0        0        0    52024 
     690   0.16325702      470        0        0     6275    52688 
     700   0.16477754     6701        0        0        0     1928 
     710   0.16572805     6681        0        0        0     1832 
     720    0.1666681     6668        0        0        0     1736 
     730   0.16806141     6654        0        0        0     1648 
     740   0.16898307     6638        0        0        0     1568 
     750   0.17031971     6623        0        0        0    51624 
     760    0.1744397     6605        0        0        0    51752 
     770   0.17854852     6593        0        0        0    51592 
     780   0.18258663     6579        0        0        0    51480 
     790   0.18816435      449        0        0     6204    52096 
     800   0.18967692     6626        0        0        0     1920 
     810   0.19063454     6610        0        0        0     1856 
     820   0.19152028     6599        0        0        0     1728 
     830   0.19234537     6580        0        0        0     1648 
     840   0.19359969     6564        0        0        0     1520 
     850    0.1949392     6546        0        0        0    51048 
     860   0.19907621     6535        0        0        0    51160 
     870    0.2030909     6512        0        0        0    51008 
     880   0.20710713     6498        0        0        0    50896 
     890   0.21307935      449        0        0     6141    51624 
     900   0.21459233     6628        0        0        0     1920 
     910   0.21546683     6615        0        0        0     1800 
     920   0.21629751     6603        0        0        0     1728 
     930   0.21714136     6590        0        0        0     1656 
     940   0.21798448     6577        0        0        0     1584 
     950   0.21928043     6557        0        0        0    51304 
     960    0.2232973     6542        0        0        0    51360 
     970   0.22727855     6529        0        0        0    51216 
     980   0.23150473     6519        0        0        0    51144 
     990   0.23750588      471        0        0     6165    52048 
    1000   0.23982323     6669        0        0        0     2312 
Loop time of 0.239878 on 1 procs for 1000 steps with 6669 particles
Performance: 4168.791 timesteps/s, 27.802 Mparticle-step/s

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.2199     | 0.2199     | 0.2199     |   0.0 | 91.67
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00066918 | 0.00066918 | 0.00066918 |   0.0 |  0.28
Modify  | 0.014353   | 0.014353   | 0.014353   |   0.0 |  5.98
Output  | 0.0033003  | 0.0033003  | 0.0033003  |   0.0 |  1.38
MPI Sync| 0.00028946 | 0.00028946 | 0.00028946 |   0.0 |  0.12
Other   |            | 0.001366   |            |       |  0.57

Particle moves    = 5397929 (5.4M)
Cells touched     = 5456156 (5.46M)
Particle comms    = 0 (0K)
Boundary collides = 0 (0K)
Boundary exits    = 934 (0.934K)
SurfColl checks   = 19766792 (19.8M)
SurfColl occurs   = 56634 (56.6K)
Surf reactions    = 56634 (56.6K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.25028e+07
Particle-moves/step: 5397.93
Cell-touches/particle/step: 1.01079
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0
Particle fraction exiting boundary: 0.000173029
Surface-checks/particle/step: 3.66192
Surface-collisions/particle/step: 0.0104918
Surf-reactions/particle/step: 0.0104918
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 56634
    reaction O(g) --> O(s): 42304
    reaction O(g) + O(s) --> CO2(g): 8
    reaction O(g) --> CO(s): 13008
    reaction O(g) --> CO(g): 1291
    reaction O(g) + O(s) --> O(g) + O(g): 23

Particles: 6669 ave 6669 max 6669 min
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
