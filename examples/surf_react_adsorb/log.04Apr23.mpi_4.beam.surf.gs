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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/ascldap/users/stamoor/sparta_master/src/grid.cpp:465)
Created 8 child grid cells
  CPU time = 0.00110793 secs
  create/ghost percent = 84.1403 15.8597
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000782251 secs
  reassign/sort/migrate/ghost percent = 72.356 0.365742 11.399 15.8793

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
  CPU time = 0.00106382 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 23.4872 24.5406 0.672344 43.3438 7.95607 25.6163 0.627521
  surf2grid time = 0.000461102 secs
  map/comm1/comm2/comm3/comm4/split percent = 26.0083 11.3237 5.11892 5.53257 8.53154 34.0228

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
      10 0.0002822876        0        0        0        0        0 
      20 0.00054860115        0        0        0        0        0 
      30 0.00082111359        0        0        0        0        0 
      40 0.0011005402        0        0        0        0        0 
      50 0.0013580322        0        0        0        0        0 
      60  0.001635313        0        0        0        0        0 
      70 0.0019330978        0        0        0        0        0 
      80 0.0021941662        0        0        0        0        0 
      90 0.0024549961        0        0        0        0        0 
     100 0.0056619644     6273        0        0        0        0 
     110 0.0069377422     6273        0        0        0        0 
     120  0.008212328     6273        0        0        0        0 
     130 0.0094738007     6273        0        0        0        0 
     140  0.010744572     6273        0        0        0        0 
     150  0.021093369     6273        0        0        0    49952 
     160   0.02594614     6273        0        0        0    50184 
     170  0.030725241     6273        0        0        0    50184 
     180  0.035509825     6273        0        0        0    50184 
     190  0.042361498      182        0        0     6178    50888 
     200  0.045572996     6450        0        0        0     1040 
     210  0.047000647     6450        0        0        0     1016 
     220  0.048449039     6450        0        0        0     1016 
     230  0.049862862     6450        0        0        0     1016 
     240  0.051260471     6450        0        0        0     1000 
     250  0.053056717     6450        0        0        0    51336 
     260  0.058002234     6450        0        0        0    51432 
     270  0.062905788     6449        0        0        0    51336 
     280  0.067798376     6448        0        0        0    51248 
     290  0.073769093      310        0        0     6233    51928 
     300  0.076169729     6572        0        0        0     1656 
     310  0.077634096     6569        0        0        0     1608 
     320  0.079109907     6564        0        0        0     1544 
     330   0.08058238     6556        0        0        0     1480 
     340  0.082043409     6551        0        0        0     1424 
     350  0.083934546     6545        0        0        0    51608 
     360  0.088943005     6537        0        0        0    51688 
     370  0.093931437     6529        0        0        0    51536 
     380  0.098940372     6520        0        0        0    51352 
     390   0.10496497      414        0        0     6211    52152 
     400   0.10744643     6616        0        0        0     1976 
     410    0.1089685     6604        0        0        0     1880 
     420   0.11054158     6598        0        0        0     1792 
     430   0.11201143     6587        0        0        0     1776 
     440   0.11346936     6572        0        0        0     1736 
     450    0.1153214     6563        0        0        0    51400 
     460    0.1202352     6555        0        0        0    51496 
     470   0.12516212     6541        0        0        0    51328 
     480   0.13002586     6533        0        0        0    51224 
     490   0.13596177      454        0        0     6141    51800 
     500   0.13843727     6662        0        0        0     2064 
     510   0.13997173     6643        0        0        0     1968 
     520   0.14146161     6632        0        0        0     1888 
     530   0.14293885     6622        0        0        0     1824 
     540    0.1444056     6607        0        0        0     1784 
     550   0.14626026     6588        0        0        0    51632 
     560   0.15127397     6575        0        0        0    51736 
     570   0.15627718     6555        0        0        0    51552 
     580   0.16126943     6546        0        0        0    51448 
     590   0.16726685      467        0        0     6171    52136 
     600    0.1697371     6684        0        0        0     2208 
     610   0.17124104     6675        0        0        0     2032 
     620   0.17272258     6659        0        0        0     1920 
     630   0.17420006     6645        0        0        0     1840 
     640   0.17563152     6628        0        0        0     1720 
     650   0.17752481     6615        0        0        0    51688 
     660   0.18253756     6602        0        0        0    51704 
     670    0.1875627     6587        0        0        0    51536 
     680   0.19253731     6572        0        0        0    51312 
     690   0.19856811      480        0        0     6196    52064 
     700   0.20105505     6714        0        0        0     2064 
     710   0.20258403     6700        0        0        0     1976 
     720   0.20409846     6686        0        0        0     1888 
     730   0.20560002     6673        0        0        0     1840 
     740   0.20707583     6650        0        0        0     1696 
     750    0.2089889     6630        0        0        0    51664 
     760   0.21393752     6620        0        0        0    51744 
     770   0.21887112     6601        0        0        0    51600 
     780   0.22380948     6583        0        0        0    51440 
     790    0.2297368      454        0        0     6205    51992 
     800    0.2322166     6684        0        0        0     1968 
     810   0.23370147     6661        0        0        0     1888 
     820   0.23522162     6640        0        0        0     1768 
     830   0.23669124     6624        0        0        0     1696 
     840   0.23813868     6615        0        0        0     1600 
     850   0.24014616     6606        0        0        0    51504 
     860    0.2451129     6590        0        0        0    51600 
     870   0.25008583     6573        0        0        0    51496 
     880   0.25502205     6558        0        0        0    51352 
     890   0.26108289      478        0        0     6183    52152 
     900   0.26349163     6581        0        0        0     2160 
     910   0.26498008     6562        0        0        0     2088 
     920   0.26645398     6550        0        0        0     2032 
     930   0.26795244     6535        0        0        0     1904 
     940   0.26940846     6510        0        0        0     1736 
     950   0.27126193     6496        0        0        0    50656 
     960   0.27614808     6485        0        0        0    50760 
     970   0.28098989     6472        0        0        0    50600 
     980   0.28581715     6459        0        0        0    50480 
     990   0.29165602      476        0        0     6060    51096 
    1000   0.29410839     6645        0        0        0     2040 
Loop time of 0.294178 on 4 procs for 1000 steps with 6645 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0043628  | 0.11534    | 0.22713    |  32.7 | 39.21
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.011389   | 0.018519   | 0.021544   |   3.0 |  6.30
Modify  | 0.00018859 | 0.0097092  | 0.019351   |   9.7 |  3.30
Output  | 0.0030689  | 0.0042235  | 0.0076494  |   3.0 |  1.44
Other   |            | 0.1464     |            |       | 49.76

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

Particle-moves/CPUsec/proc: 4.57965e+06
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
