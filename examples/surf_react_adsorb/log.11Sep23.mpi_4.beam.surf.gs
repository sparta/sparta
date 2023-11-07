SPARTA (13 Apr 2023)
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
  CPU time = 0.0018162 secs
  create/ghost percent = 92.3302 7.66985
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000645301 secs
  reassign/sort/migrate/ghost percent = 66.6202 0.821322 11.359 21.1994

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
  CPU time = 0.0010869 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 21.0875 9.10855 0.60723 55.5249 13.6719 11.7306 1.00285
  surf2grid time = 0.000603501 secs
  map/comm1/comm2/comm3/comm4/split percent = 32.1292 10.3562 5.23611 5.36867 16.5534 23.8442

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
      10  0.000435501        0        0        0        0        0 
      20  0.000858401        0        0        0        0        0 
      30  0.001286002        0        0        0        0        0 
      40  0.001713803        0        0        0        0        0 
      50  0.002131903        0        0        0        0        0 
      60  0.002576704        0        0        0        0        0 
      70  0.003013904        0        0        0        0        0 
      80  0.003432605        0        0        0        0        0 
      90  0.003847506        0        0        0        0        0 
     100  0.007605011     6273        0        0        0        0 
     110  0.008798413     6273        0        0        0        0 
     120  0.009956914     6273        0        0        0        0 
     130  0.011152016     6273        0        0        0        0 
     140  0.012304518     6273        0        0        0        0 
     150   0.01385412     6273        0        0        0    49952 
     160  0.018253626     6273        0        0        0    50184 
     170  0.022685933     6273        0        0        0    50184 
     180  0.027090439     6273        0        0        0    50184 
     190  0.034282649      182        0        0     6178    50888 
     200  0.037477854     6450        0        0        0     1040 
     210  0.038834056     6450        0        0        0     1016 
     220  0.040243558     6450        0        0        0     1016 
     230   0.04156996     6450        0        0        0     1016 
     240  0.042885661     6450        0        0        0     1000 
     250  0.044558064     6450        0        0        0    51336 
     260   0.04910687     6450        0        0        0    51432 
     270  0.053632877     6449        0        0        0    51336 
     280  0.058159183     6448        0        0        0    51248 
     290  0.064553192      310        0        0     6233    51928 
     300  0.066657695     6572        0        0        0     1656 
     310  0.068065997     6569        0        0        0     1608 
     320  0.069483699     6564        0        0        0     1544 
     330  0.070940802     6556        0        0        0     1480 
     340  0.072354904     6551        0        0        0     1424 
     350  0.074111006     6545        0        0        0    51608 
     360  0.078767413     6537        0        0        0    51688 
     370  0.083373019     6529        0        0        0    51536 
     380  0.088165126     6520        0        0        0    51352 
     390  0.094565135      414        0        0     6211    52152 
     400  0.096708938     6616        0        0        0     1976 
     410  0.098188441     6604        0        0        0     1880 
     420  0.099665543     6598        0        0        0     1792 
     430   0.10103894     6587        0        0        0     1776 
     440   0.10243875     6572        0        0        0     1736 
     450   0.10419845     6563        0        0        0    51400 
     460   0.10881336     6555        0        0        0    51496 
     470   0.11345416     6541        0        0        0    51328 
     480   0.11797257     6533        0        0        0    51224 
     490   0.12433728      454        0        0     6141    51800 
     500   0.12654568     6662        0        0        0     2064 
     510   0.12802228     6643        0        0        0     1968 
     520   0.12950619     6632        0        0        0     1888 
     530   0.13124809     6622        0        0        0     1824 
     540   0.13341009     6607        0        0        0     1784 
     550   0.13531119     6588        0        0        0    51632 
     560    0.1399347     6575        0        0        0    51736 
     570   0.14461811     6555        0        0        0    51552 
     580   0.14921931     6546        0        0        0    51448 
     590   0.15568462      467        0        0     6171    52136 
     600   0.15781183     6684        0        0        0     2208 
     610   0.15926943     6675        0        0        0     2032 
     620   0.16075813     6659        0        0        0     1920 
     630   0.16219033     6645        0        0        0     1840 
     640   0.16355783     6628        0        0        0     1720 
     650   0.16539634     6615        0        0        0    51688 
     660   0.17005184     6602        0        0        0    51704 
     670   0.17481635     6587        0        0        0    51536 
     680   0.17944476     6572        0        0        0    51312 
     690   0.18587147      480        0        0     6196    52064 
     700   0.18800467     6714        0        0        0     2064 
     710   0.18952947     6700        0        0        0     1976 
     720   0.19102727     6686        0        0        0     1888 
     730   0.19244518     6673        0        0        0     1840 
     740   0.19391038     6650        0        0        0     1696 
     750   0.19578418     6630        0        0        0    51664 
     760   0.20043169     6620        0        0        0    51744 
     770   0.20508859     6601        0        0        0    51600 
     780    0.2097361     6583        0        0        0    51440 
     790   0.21616761      454        0        0     6205    51992 
     800   0.21832481     6684        0        0        0     1968 
     810   0.21978062     6661        0        0        0     1888 
     820   0.22128052     6640        0        0        0     1768 
     830   0.22271612     6624        0        0        0     1696 
     840   0.22410482     6615        0        0        0     1600 
     850   0.22586512     6606        0        0        0    51504 
     860   0.23049193     6590        0        0        0    51600 
     870   0.23514114     6573        0        0        0    51496 
     880   0.23971404     6558        0        0        0    51352 
     890   0.24607915      478        0        0     6183    52152 
     900   0.24818235     6581        0        0        0     2160 
     910   0.24966116     6562        0        0        0     2088 
     920   0.25110526     6550        0        0        0     2032 
     930   0.25250516     6535        0        0        0     1904 
     940   0.25391626     6510        0        0        0     1736 
     950   0.25575137     6496        0        0        0    50656 
     960   0.26031517     6485        0        0        0    50760 
     970   0.26484298     6472        0        0        0    50600 
     980   0.26935568     6459        0        0        0    50480 
     990   0.27567349      476        0        0     6060    51096 
    1000    0.2777582     6645        0        0        0     2040 
Loop time of 0.277808 on 4 procs for 1000 steps with 6645 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0020757  | 0.052934   | 0.10429    |  22.1 | 19.05
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.020101   | 0.022731   | 0.026772   |   1.6 |  8.18
Modify  | 8.69e-05   | 0.0038357  | 0.0075932  |   6.1 |  1.38
Output  | 0.0049898  | 0.0077064  | 0.015178   |   4.9 |  2.77
Other   |            | 0.1906     |            |       | 68.61

Particle moves    = 5369855 (5.37M)
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

Particle-moves/CPUsec/proc: 4.83235e+06
Particle-moves/step: 5369.85
Cell-touches/particle/step: 1.01438
Particle comm iterations/step: 1.456
Particle fraction communicated: 0.000135013
Particle fraction colliding with boundary: 0
Particle fraction exiting boundary: 0.000175238
Surface-checks/particle/step: 3.67987
Surface-collisions/particle/step: 0.0105152
Surf-reactions/particle/step: 0.0105152
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
