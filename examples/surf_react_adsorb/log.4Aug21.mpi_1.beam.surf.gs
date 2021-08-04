SPARTA (26 Feb 2021)
# beam of particles striking the surface at an inclined angle
# free molecular flow (no collisions)

seed	    	    123456
dimension   	    3
global              gridcut 0.0 comm/sort yes

boundary	    	oo oo oo


create_box          -11 11 -11 11 0 10
Created orthogonal box = (-11 -11 0) to (11 11 10)
create_grid 	    2 2 2
Created 8 child grid cells
  CPU time = 0.00171248 secs
  create/ghost percent = 93.9067 6.09325
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000304516 secs
  reassign/sort/migrate/ghost percent = 79.8851 0.229545 14.5641 5.32123

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
  CPU time = 0.00246001 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 67.154 5.38867 0.210121 24.1071 3.14006 4.04294 0.0369104
  surf2grid time = 0.000593038 secs
  map/comm1/comm2/comm3/comm4/split percent = 29.3015 6.52454 3.46251 2.69696 12.743 42.0798

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
      10   8.7583e-05        0        0        0        0        0 
      20  0.000161966        0        0        0        0        0 
      30    0.0002316        0        0        0        0        0 
      40  0.000300395        0        0        0        0        0 
      50  0.000369959        0        0        0        0        0 
      60  0.000440012        0        0        0        0        0 
      70  0.000510623        0        0        0        0        0 
      80  0.000580536        0        0        0        0        0 
      90   0.00064996        0        0        0        0        0 
     100  0.006587393     6253        0        0        0        0 
     110  0.008345834     6253        0        0        0        0 
     120  0.010115169     6253        0        0        0        0 
     130  0.011879965     6253        0        0        0        0 
     140  0.013527704     6253        0        0        0        0 
     150  0.015966696     6253        0        0        0    49816 
     160  0.023985733     6253        0        0        0    50024 
     170  0.031560708     6253        0        0        0    50024 
     180  0.038815243     6253        0        0        0    50024 
     190  0.047521677      187        0        0     6152    50712 
     200  0.050999026     6454        0        0        0     1072 
     210  0.052525797     6454        0        0        0     1072 
     220  0.054070098     6454        0        0        0     1088 
     230  0.055555102     6454        0        0        0     1080 
     240  0.057028233     6454        0        0        0     1056 
     250  0.059125762     6454        0        0        0    51216 
     260  0.065692901     6453        0        0        0    51464 
     270  0.072018522     6453        0        0        0    51344 
     280  0.078137548     6453        0        0        0    51256 
     290  0.085588621      345        0        0     6208    52000 
     300  0.088569665     6621        0        0        0     1808 
     310  0.089987969     6616        0        0        0     1760 
     320  0.091355847     6608        0        0        0     1672 
     330  0.092713178     6601        0        0        0     1592 
     340  0.094075817     6596        0        0        0     1552 
     350  0.095966401     6588        0        0        0    51792 
     360   0.10182805     6579        0        0        0    51904 
     370   0.10764787     6574        0        0        0    51808 
     380   0.11346832     6566        0        0        0    51720 
     390   0.12066649      451        0        0     6229    52552 
     400   0.12367987     6642        0        0        0     2152 
     410   0.12509615     6633        0        0        0     2096 
     420   0.12649036     6623        0        0        0     1968 
     430   0.12788624     6603        0        0        0     1912 
     440    0.1292602     6592        0        0        0     1792 
     450   0.13117893     6579        0        0        0    51504 
     460   0.13698631     6567        0        0        0    51488 
     470   0.14276129     6556        0        0        0    51280 
     480   0.14853194     6543        0        0        0    51128 
     490   0.15564099      464        0        0     6152    51728 
     500   0.15862958     6671        0        0        0     2040 
     510   0.16005172     6651        0        0        0     1912 
     520   0.16143902     6631        0        0        0     1856 
     530   0.16282736     6612        0        0        0     1736 
     540   0.16419656     6593        0        0        0     1672 
     550   0.16610132     6577        0        0        0    51368 
     560   0.17190529     6563        0        0        0    51472 
     570   0.17768494     6548        0        0        0    51368 
     580   0.18345517     6539        0        0        0    51208 
     590   0.19056458      476        0        0     6158    51968 
     600   0.19356622     6639        0        0        0     2096 
     610    0.1949744     6619        0        0        0     1984 
     620    0.1963596     6606        0        0        0     1880 
     630   0.19773146     6587        0        0        0     1816 
     640   0.19911002     6576        0        0        0     1752 
     650   0.20101625     6569        0        0        0    51336 
     660   0.20681344     6557        0        0        0    51376 
     670   0.21258499     6544        0        0        0    51240 
     680     0.218346     6532        0        0        0    51136 
     690   0.22543221      463        0        0     6142    51672 
     700   0.22846669     6771        0        0        0     2056 
     710   0.22989554     6751        0        0        0     1928 
     720   0.23129506     6734        0        0        0     1832 
     730   0.23268926     6716        0        0        0     1736 
     740   0.23409507     6701        0        0        0     1624 
     750   0.23602665     6687        0        0        0    52192 
     760   0.24191924     6671        0        0        0    52296 
     770   0.24779465     6662        0        0        0    52136 
     780   0.25365826     6648        0        0        0    52016 
     790    0.2608982      481        0        0     6262    52712 
     800   0.26388665     6665        0        0        0     2096 
     810   0.26528952     6650        0        0        0     2008 
     820   0.26668596     6637        0        0        0     1952 
     830   0.26809267     6628        0        0        0     1864 
     840   0.26947165     6613        0        0        0     1760 
     850   0.27139632     6596        0        0        0    51440 
     860    0.2772218     6580        0        0        0    51520 
     870   0.28301528     6568        0        0        0    51384 
     880   0.28880116     6551        0        0        0    51248 
     890   0.29592529      500        0        0     6148    51968 
     900   0.29893574     6704        0        0        0     2264 
     910   0.30035956     6688        0        0        0     2160 
     920   0.30178611     6674        0        0        0     2096 
     930    0.3031903     6661        0        0        0     1984 
     940   0.30458403     6648        0        0        0     1880 
     950   0.30652252     6632        0        0        0    51624 
     960   0.31236224     6621        0        0        0    51744 
     970   0.31818569     6607        0        0        0    51632 
     980   0.32399273     6587        0        0        0    51448 
     990   0.33114837      537        0        0     6163    52336 
    1000   0.33417935     6694        0        0        0     2440 
Loop time of 0.334191 on 1 procs for 1000 steps with 6694 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.29908    | 0.29908    | 0.29908    |   0.0 | 89.49
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0012268  | 0.0012268  | 0.0012268  |   0.0 |  0.37
Modify  | 0.029755   | 0.029755   | 0.029755   |   0.0 |  8.90
Output  | 0.0014479  | 0.0014479  | 0.0014479  |   0.0 |  0.43
Other   |            | 0.002682   |            |       |  0.80

Particle moves    = 5405885 (5.41M)
Cells touched     = 5464228 (5.46M)
Particle comms    = 0 (0K)
Boundary collides = 0 (0K)
Boundary exits    = 966 (0.966K)
SurfColl checks   = 19832568 (19.8M)
SurfColl occurs   = 56577 (56.6K)
Surf reactions    = 56577 (56.6K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 1.6176e+07
Particle-moves/step: 5405.89
Cell-touches/particle/step: 1.01079
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0
Particle fraction exiting boundary: 0.000178694
Surface-checks/particle/step: 3.6687
Surface-collisions/particle/step: 0.0104658
Surf-reactions/particle/step: 0.0104658
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 56577
    reaction O(g) --> O(s): 42256
    reaction O(g) + O(s) --> CO2(g): 6
    reaction O(g) --> CO(s): 12923
    reaction O(g) --> CO(g): 1362
    reaction O(g) + O(s) --> O(g) + O(g): 30

Particles: 6694 ave 6694 max 6694 min
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
