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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/runner/work/sparta/sparta/src/grid.cpp:465)
Created 8 child grid cells
  CPU time = 0.00294154 secs
  create/ghost percent = 94.1832 5.81675
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.00112051 secs
  reassign/sort/migrate/ghost percent = 66.6043 1.22265 9.17445 22.9986

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
  CPU time = 0.00200213 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 21.5824 9.04049 2.24267 50.477 16.6575 13.5557 0.963975
  surf2grid time = 0.00101061 secs
  map/comm1/comm2/comm3/comm4/split percent = 28.6364 12.3392 6.3032 6.1053 18.1774 18.9788

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
      10  0.000869011        0        0        0        0        0 
      20  0.001756923        0        0        0        0        0 
      30  0.002661934        0        0        0        0        0 
      40  0.003519045        0        0        0        0        0 
      50  0.004389356        0        0        0        0        0 
      60  0.005226967        0        0        0        0        0 
      70  0.006155378        0        0        0        0        0 
      80  0.007171291        0        0        0        0        0 
      90  0.008189604        0        0        0        0        0 
     100  0.014445484     6273        0        0        0        0 
     110   0.01654981     6273        0        0        0        0 
     120  0.018612736     6273        0        0        0        0 
     130  0.020676563     6273        0        0        0        0 
     140  0.022761189     6273        0        0        0        0 
     150  0.025490124     6273        0        0        0    49952 
     160  0.033133321     6273        0        0        0    50184 
     170   0.04099822     6273        0        0        0    50184 
     180  0.048781319     6273        0        0        0    50184 
     190  0.062440592      182        0        0     6178    50888 
     200  0.067423455     6450        0        0        0     1040 
     210  0.068386368     6450        0        0        0     1016 
     220  0.069455081     6450        0        0        0     1016 
     230  0.070440394     6450        0        0        0     1016 
     240  0.071770011     6450        0        0        0     1000 
     250  0.073539133     6450        0        0        0    51336 
     260  0.078655898     6450        0        0        0    51432 
     270  0.083741962     6449        0        0        0    51336 
     280  0.088932428     6448        0        0        0    51248 
     290  0.095926717      310        0        0     6233    51928 
     300  0.098057144     6572        0        0        0     1656 
     310  0.099530063     6569        0        0        0     1608 
     320   0.10103398     6564        0        0        0     1544 
     330    0.1025292     6556        0        0        0     1480 
     340   0.10394662     6551        0        0        0     1424 
     350   0.10583994     6545        0        0        0    51608 
     360   0.11104761     6537        0        0        0    51688 
     370   0.11622917     6529        0        0        0    51536 
     380   0.12145414     6520        0        0        0    51352 
     390   0.12848093      414        0        0     6211    52152 
     400   0.13068926     6616        0        0        0     1976 
     410   0.13222148     6604        0        0        0     1880 
     420    0.1337306     6598        0        0        0     1792 
     430   0.13516432     6587        0        0        0     1776 
     440   0.13659643     6572        0        0        0     1736 
     450   0.13843466     6563        0        0        0    51400 
     460   0.14358982     6555        0        0        0    51496 
     470   0.14874539     6541        0        0        0    51328 
     480   0.15380715     6533        0        0        0    51224 
     490   0.16072824      454        0        0     6141    51800 
     500   0.16296917     6662        0        0        0     2064 
     510   0.16449009     6643        0        0        0     1968 
     520   0.16601481     6632        0        0        0     1888 
     530   0.16747712     6622        0        0        0     1824 
     540   0.16894964     6607        0        0        0     1784 
     550   0.17081487     6588        0        0        0    51632 
     560   0.17600783     6575        0        0        0    51736 
     570    0.1811988     6555        0        0        0    51552 
     580   0.18693607     6546        0        0        0    51448 
     590   0.19427766      467        0        0     6171    52136 
     600   0.19656519     6684        0        0        0     2208 
     610   0.19808331     6675        0        0        0     2032 
     620   0.19964053     6659        0        0        0     1920 
     630   0.20111885     6645        0        0        0     1840 
     640   0.20250967     6628        0        0        0     1720 
     650   0.20444009     6615        0        0        0    51688 
     660   0.20974386     6602        0        0        0    51704 
     670   0.21508083     6587        0        0        0    51536 
     680   0.22025619     6572        0        0        0    51312 
     690   0.22730228      480        0        0     6196    52064 
     700   0.22953551     6714        0        0        0     2064 
     710   0.23112773     6700        0        0        0     1976 
     720   0.23264265     6686        0        0        0     1888 
     730   0.23410057     6673        0        0        0     1840 
     740   0.23556489     6650        0        0        0     1696 
     750   0.23751451     6630        0        0        0    51664 
     760   0.24270558     6620        0        0        0    51744 
     770   0.24794544     6601        0        0        0    51600 
     780   0.25313581     6583        0        0        0    51440 
     790    0.2600714      454        0        0     6205    51992 
     800   0.26233853     6684        0        0        0     1968 
     810   0.26383515     6661        0        0        0     1888 
     820   0.26537347     6640        0        0        0     1768 
     830   0.26718599     6624        0        0        0     1696 
     840   0.26860561     6615        0        0        0     1600 
     850   0.27044413     6606        0        0        0    51504 
     860    0.2756445     6590        0        0        0    51600 
     870   0.28081746     6573        0        0        0    51496 
     880   0.28596643     6558        0        0        0    51352 
     890   0.29294742      478        0        0     6183    52152 
     900   0.29507744     6581        0        0        0     2160 
     910   0.29660126     6562        0        0        0     2088 
     920   0.29806708     6550        0        0        0     2032 
     930    0.2995115     6535        0        0        0     1904 
     940   0.30095472     6510        0        0        0     1736 
     950   0.30286464     6496        0        0        0    50656 
     960   0.30801831     6485        0        0        0    50760 
     970   0.31308387     6472        0        0        0    50600 
     980   0.31814114     6459        0        0        0    50480 
     990   0.32492942      476        0        0     6060    51096 
    1000   0.32712275     6645        0        0        0     2040 
Loop time of 0.327205 on 4 procs for 1000 steps with 6645 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0022356  | 0.058869   | 0.11614    |  23.3 | 17.99
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.052349   | 0.10873    | 0.15898    |  14.9 | 33.23
Modify  | 0.0001021  | 0.0038641  | 0.0076832  |   6.0 |  1.18
Output  | 0.005908   | 0.0077952  | 0.0096602  |   1.5 |  2.38
Other   |            | 0.1479     |            |       | 45.22

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

Particle-moves/CPUsec/proc: 4.1174e+06
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
