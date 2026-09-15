SPARTA (24 Sep 2025)
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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/user/sparta/src/grid.cpp:486)
Created 8 child grid cells
  CPU time = 0.0017971 secs
  create/ghost percent = 91.0973 8.90274
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000372161 secs
  reassign/sort/migrate/ghost percent = 67.7825 0.613982 11.2508 20.3528

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
  CPU time = 0.000799601 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 5.75475 33.7525 0.343296 51.606 8.54351 7.19696 0.125438
  surf2grid time = 0.000412642 secs
  map/comm1/comm2/comm3/comm4/split percent = 36.8421 13.0658 5.47133 7.02691 18.3539 14.8288

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
      10  0.000137864        0        0        0        0        0 
      20  0.000274756        0        0        0        0        0 
      30  0.000412627        0        0        0        0        0 
      40  0.000541888        0        0        0        0        0 
      50  0.000829103        0        0        0        0        0 
      60  0.001045218        0        0        0        0        0 
      70  0.001246048        0        0        0        0        0 
      80  0.001497987        0        0        0        0        0 
      90  0.001760178        0        0        0        0        0 
     100  0.004343575     6270        0        0        0        0 
     110  0.005414411     6270        0        0        0        0 
     120   0.00641353     6270        0        0        0        0 
     130  0.007366954     6270        0        0        0        0 
     140  0.008321646     6270        0        0        0        0 
     150  0.009790516     6270        0        0        0    49928 
     160  0.014513643     6270        0        0        0    50160 
     170  0.019233217     6270        0        0        0    50160 
     180  0.024028423     6270        0        0        0    50160 
     190  0.031465869      181        0        0     6176    50864 
     200  0.034209028     6450        0        0        0     1024 
     210  0.035292864     6450        0        0        0     1008 
     220  0.036410767     6450        0        0        0     1024 
     230  0.037523471     6450        0        0        0     1016 
     240  0.038564999     6450        0        0        0     1000 
     250  0.040117406     6450        0        0        0    51336 
     260  0.045069072     6450        0        0        0    51432 
     270  0.049947977     6449        0        0        0    51336 
     280  0.054817433     6448        0        0        0    51248 
     290  0.060995606      312        0        0     6233    51936 
     300  0.062533784     6573        0        0        0     1656 
     310  0.063665993     6568        0        0        0     1584 
     320  0.064876543     6562        0        0        0     1528 
     330  0.065604482     6555        0        0        0     1472 
     340  0.066275919     6550        0        0        0     1416 
     350  0.067256016     6543        0        0        0    51584 
     360  0.069855727     6535        0        0        0    51672 
     370  0.072387326     6527        0        0        0    51520 
     380  0.074906841     6520        0        0        0    51368 
     390  0.078373806      412        0        0     6211    52160 
     400  0.079427881     6614        0        0        0     1984 
     410  0.080208767     6602        0        0        0     1880 
     420  0.081010259     6599        0        0        0     1824 
     430  0.081688563     6587        0        0        0     1784 
     440  0.082320199     6572        0        0        0     1728 
     450  0.083252341     6559        0        0        0    51392 
     460  0.085860588     6555        0        0        0    51496 
     470  0.088400632     6540        0        0        0    51328 
     480  0.090946662     6530        0        0        0    51224 
     490   0.09435115      453        0        0     6141    51792 
     500  0.095384456     6664        0        0        0     2048 
     510  0.096041204     6648        0        0        0     2008 
     520   0.09672998     6636        0        0        0     1936 
     530  0.097408892     6624        0        0        0     1856 
     540  0.098049483     6606        0        0        0     1760 
     550  0.098958611     6586        0        0        0    51616 
     560   0.10152543     6576        0        0        0    51728 
     570   0.10416765     6560        0        0        0    51568 
     580   0.10671101     6546        0        0        0    51400 
     590   0.11012873      462        0        0     6171    52104 
     600   0.11118094     6681        0        0        0     2168 
     610   0.11190403     6670        0        0        0     2048 
     620   0.11252463     6655        0        0        0     1928 
     630   0.11326655     6643        0        0        0     1856 
     640   0.11393542     6629        0        0        0     1720 
     650   0.11481945     6614        0        0        0    51672 
     660   0.11741351     6597        0        0        0    51648 
     670   0.11996215     6585        0        0        0    51544 
     680   0.12257015     6573        0        0        0    51344 
     690   0.12609812      484        0        0     6196    52088 
     700   0.12715104     6715        0        0        0     2088 
     710   0.12791471     6706        0        0        0     1968 
     720   0.12867701     6689        0        0        0     1904 
     730   0.12938005     6675        0        0        0     1864 
     740   0.13022623     6654        0        0        0     1704 
     750   0.13115827     6629        0        0        0    51680 
     760   0.13386975     6616        0        0        0    51744 
     770   0.13637484     6599        0        0        0    51600 
     780   0.13894225     6586        0        0        0    51432 
     790   0.14325682      458        0        0     6205    52016 
     800   0.14433656     6689        0        0        0     1984 
     810   0.14508294     6666        0        0        0     1920 
     820   0.14585499     6648        0        0        0     1784 
     830   0.14651944     6628        0        0        0     1720 
     840   0.14715817     6613        0        0        0     1600 
     850   0.14808249     6603        0        0        0    51504 
     860   0.15069034     6590        0        0        0    51624 
     870   0.15334745     6576        0        0        0    51488 
     880   0.15606628     6565        0        0        0    51368 
     890   0.15959123      481        0        0     6183    52152 
     900   0.16067269     6582        0        0        0     2144 
     910   0.16133542     6561        0        0        0     2088 
     920   0.16218244     6551        0        0        0     2016 
     930   0.16290252     6536        0        0        0     1896 
     940   0.16372776     6513        0        0        0     1768 
     950   0.16460089     6504        0        0        0    50680 
     960   0.16734574     6492        0        0        0    50792 
     970   0.16988921     6475        0        0        0    50608 
     980   0.17232908     6461        0        0        0    50472 
     990   0.17583875      477        0        0     6060    51112 
    1000   0.17691915     6645        0        0        0     2040 
Loop time of 0.176973 on 4 procs for 1000 steps with 6645 particles
Performance: 5650.588 timesteps/s, 37.548 Mparticle-step/s

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0023472  | 0.065535   | 0.13964    |  24.8 | 37.03
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.008541   | 0.0089149  | 0.0093901  |   0.4 |  5.04
Modify  | 0.00014889 | 0.004435   | 0.0093131  |   6.5 |  2.51
Output  | 0.0017069  | 0.0027988  | 0.0060288  |   3.5 |  1.58
MPI Sync| 0.0096508  | 0.087962   | 0.15681    |  23.3 | 49.70
Other   |            | 0.007327   |            |       |  4.14

Particle moves    = 5351732 (5.35M)
Cells touched     = 5447142 (5.45M)
Particle comms    = 723 (0.723K)
Boundary collides = 0 (0K)
Boundary exits    = 941 (0.941K)
SurfColl checks   = 19762256 (19.8M)
SurfColl occurs   = 56462 (56.5K)
Surf reactions    = 56462 (56.5K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 7.56011e+06
Particle-moves/step: 5351.73
Cell-touches/particle/step: 1.01783
Particle comm iterations/step: 1.455
Particle fraction communicated: 0.000135096
Particle fraction colliding with boundary: 0
Particle fraction exiting boundary: 0.000175831
Surface-checks/particle/step: 3.69268
Surface-collisions/particle/step: 0.0105502
Surf-reactions/particle/step: 0.0105502
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 56462
    reaction O(g) --> O(s): 41986
    reaction O(g) + O(s) --> CO2(g): 9
    reaction O(g) --> CO(s): 13164
    reaction O(g) --> CO(g): 1275
    reaction O(g) + O(s) --> O(g) + O(g): 28

Particles: 1661.25 ave 3281 max 51 min
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
