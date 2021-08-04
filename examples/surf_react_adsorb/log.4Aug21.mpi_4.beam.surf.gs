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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (../grid.cpp:410)
Created 8 child grid cells
  CPU time = 0.00188423 secs
  create/ghost percent = 94.818 5.18202
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000951193 secs
  reassign/sort/migrate/ghost percent = 74.5136 0.198172 11.8511 13.4371

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
  CPU time = 0.00119264 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 23.4774 9.18835 0.919472 52.2195 14.1953 11.3844 0.386537
  surf2grid time = 0.000622791 secs
  map/comm1/comm2/comm3/comm4/split percent = 31.4793 11.2707 5.82025 5.99977 15.3191 22.6421

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
      10    0.0005387        0        0        0        0        0 
      20  0.000975639        0        0        0        0        0 
      30  0.001383174        0        0        0        0        0 
      40  0.001797762        0        0        0        0        0 
      50   0.00222192        0        0        0        0        0 
      60  0.002619327        0        0        0        0        0 
      70  0.003037827        0        0        0        0        0 
      80  0.003442358        0        0        0        0        0 
      90  0.003873709        0        0        0        0        0 
     100  0.008415957     6308        0        0        0        0 
     110  0.009797104     6308        0        0        0        0 
     120  0.011168963     6308        0        0        0        0 
     130  0.012235467     6308        0        0        0        0 
     140  0.013188127     6308        0        0        0        0 
     150  0.014459971     6308        0        0        0    50248 
     160  0.017777658     6308        0        0        0    50464 
     170  0.021083473     6308        0        0        0    50464 
     180  0.024304777     6308        0        0        0    50464 
     190  0.029709447      192        0        0     6203    51168 
     200   0.03280685     6481        0        0        0     1128 
     210  0.033904225     6481        0        0        0     1136 
     220   0.03494391     6481        0        0        0     1136 
     230   0.03598094     6481        0        0        0     1128 
     240  0.037046816     6481        0        0        0     1120 
     250  0.038307834     6481        0        0        0    51512 
     260  0.041508883     6480        0        0        0    51656 
     270  0.044876439     6480        0        0        0    51544 
     280  0.048128614     6479        0        0        0    51464 
     290  0.052224075      332        0        0     6238    52128 
     300  0.054075756     6572        0        0        0     1680 
     310  0.055238574     6568        0        0        0     1624 
     320  0.056357251     6561        0        0        0     1496 
     330  0.057420193     6554        0        0        0     1448 
     340  0.058461764     6545        0        0        0     1376 
     350  0.059713283     6538        0        0        0    51352 
     360  0.063008551     6527        0        0        0    51560 
     370   0.06624983     6520        0        0        0    51400 
     380  0.069563677     6514        0        0        0    51240 
     390   0.07359593      437        0        0     6191    52096 
     400  0.075421699     6629        0        0        0     2040 
     410  0.076479612     6622        0        0        0     1952 
     420  0.077594937     6608        0        0        0     1920 
     430  0.078678972     6599        0        0        0     1840 
     440  0.079702174     6589        0        0        0     1768 
     450  0.081010754     6576        0        0        0    51488 
     460  0.084225843     6560        0        0        0    51512 
     470  0.087515943     6553        0        0        0    51328 
     480  0.090683398     6543        0        0        0    51248 
     490  0.094541881      465        0        0     6145    51808 
     500  0.096363459     6665        0        0        0     2112 
     510   0.09744037     6654        0        0        0     2048 
     520  0.098440105     6633        0        0        0     1920 
     530  0.099469802     6617        0        0        0     1824 
     540   0.10044558     6602        0        0        0     1752 
     550   0.10170485     6585        0        0        0    51448 
     560   0.10490786     6569        0        0        0    51584 
     570   0.10809221     6553        0        0        0    51400 
     580   0.11125199     6546        0        0        0    51320 
     590   0.11515985      467        0        0     6146    51816 
     600   0.11696879     6672        0        0        0     1968 
     610    0.1180552     6661        0        0        0     1888 
     620   0.11907267     6641        0        0        0     1792 
     630   0.12006577     6623        0        0        0     1664 
     640   0.12106257     6612        0        0        0     1576 
     650   0.12230976     6602        0        0        0    51456 
     660   0.12552897     6588        0        0        0    51520 
     670    0.1287018     6573        0        0        0    51304 
     680   0.13186416     6559        0        0        0    51176 
     690   0.13568863      499        0        0     6174    52064 
     700   0.13746781     6781        0        0        0     2184 
     710   0.13846468     6763        0        0        0     2120 
     720   0.13947755     6746        0        0        0     2024 
     730   0.14048371     6733        0        0        0     1984 
     740    0.1414782     6714        0        0        0     1832 
     750   0.14273293     6706        0        0        0    52272 
     760    0.1459768     6697        0        0        0    52288 
     770   0.14920488     6676        0        0        0    52168 
     780   0.15247193     6662        0        0        0    51960 
     790   0.15631868      488        0        0     6247    52616 
     800   0.15809486     6706        0        0        0     2144 
     810   0.15911184     6685        0        0        0     2016 
     820   0.16014608     6671        0        0        0     1888 
     830   0.16113981     6651        0        0        0     1808 
     840   0.16209484     6626        0        0        0     1720 
     850   0.16331654     6611        0        0        0    51544 
     860   0.16647471     6593        0        0        0    51664 
     870   0.16967233     6581        0        0        0    51512 
     880   0.17283364     6570        0        0        0    51344 
     890   0.17662731      480        0        0     6157    51880 
     900    0.1784667     6730        0        0        0     2104 
     910   0.17957707     6714        0        0        0     1960 
     920   0.18067046     6699        0        0        0     1808 
     930   0.18185989     6676        0        0        0     1728 
     940     0.182925     6657        0        0        0     1640 
     950    0.1842238     6639        0        0        0    51704 
     960   0.18740809     6622        0        0        0    51840 
     970   0.19059999     6615        0        0        0    51728 
     980   0.19384427     6607        0        0        0    51616 
     990   0.19767956      496        0        0     6203    52304 
    1000   0.19949688     6638        0        0        0     2072 
Loop time of 0.199514 on 4 procs for 1000 steps with 6638 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0068313  | 0.072084   | 0.14563    |  24.4 | 36.13
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.012713   | 0.013577   | 0.014337   |   0.5 |  6.81
Modify  | 0.0017952  | 0.0087667  | 0.017206   |   7.5 |  4.39
Output  | 0.0025734  | 0.002861   | 0.0036126  |   0.8 |  1.43
Other   |            | 0.1022     |            |       | 51.24

Particle moves    = 5416065 (5.42M)
Cells touched     = 5474472 (5.47M)
Particle comms    = 724 (0.724K)
Boundary collides = 0 (0K)
Boundary exits    = 954 (0.954K)
SurfColl checks   = 19841512 (19.8M)
SurfColl occurs   = 56689 (56.7K)
Surf reactions    = 56689 (56.7K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 6.78657e+06
Particle-moves/step: 5416.06
Cell-touches/particle/step: 1.01078
Particle comm iterations/step: 1.461
Particle fraction communicated: 0.000133676
Particle fraction colliding with boundary: 0
Particle fraction exiting boundary: 0.000176143
Surface-checks/particle/step: 3.66346
Surface-collisions/particle/step: 0.0104668
Surf-reactions/particle/step: 0.0104668
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 56689
    reaction O(g) --> O(s): 42500
    reaction O(g) + O(s) --> CO2(g): 8
    reaction O(g) --> CO(s): 12845
    reaction O(g) --> CO(g): 1307
    reaction O(g) + O(s) --> O(g) + O(g): 29

Particles: 1659.5 ave 3273 max 46 min
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
