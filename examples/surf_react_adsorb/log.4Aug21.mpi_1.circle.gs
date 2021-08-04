SPARTA (26 Feb 2021)
# 2d flow around a circle

seed	    	    12345
dimension   	    2
global              gridcut 0.0 comm/sort yes

boundary	    	o r p

create_box  	    0 10 0 10 -0.5 0.5
Created orthogonal box = (0 0 -0.5) to (10 10 0.5)
create_grid 	    20 20 1
Created 400 child grid cells
  CPU time = 0.00187005 secs
  create/ghost percent = 86.1923 13.8077
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.00044979 secs
  reassign/sort/migrate/ghost percent = 61.91 0.714333 12.0496 25.326

global		    	nrho 1.0 fnum 0.001

species		    	air.species O CO CO2 C
mixture		    	air O vstream 100.0 0 0

mixture             air O   frac 1.0
mixture             air CO  frac 0.0
mixture             air CO2 frac 0.0
mixture             air C   frac 0.0

read_surf           data.circle
  50 points
  50 lines
  2 8 xlo xhi
  2.00592 7.99408 ylo yhi
  0 0 zlo zhi
  0.376743 min line length
  0 0 = number of pushed cells
  48 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  264 88 48 = cells outside/inside/overlapping surfs
  48 = surf cells with 1,2,etc splits
  71.8 71.8 = cell-wise and global flow volume
  CPU time = 0.00123008 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 22.2235 11.4922 1.92476 53.8667 10.4928 11.1742 0.0965793
  surf2grid time = 0.000662602 secs
  map/comm1/comm2/comm3/comm4/split percent = 41.7202 6.91486 14.5673 3.83684 12.0164 15.9797

surf_collide        1 cll 300.0 0.5 0.5 0.5 0.5

################################### SURF REACT ADSORB ######################################

#surf_react          adsorb_test_gs1 adsorb gs sample-GS_1.surf nsync 1 surf 1000 6.022e9 O CO
#surf_modify         all collide 1 react adsorb_test_gs1

surf_react          adsorb_test_gs2 adsorb gs sample-GS_2.surf nsync 1 surf 1000 6.022e9 O CO
surf_modify         all collide 1 react adsorb_test_gs2

############################################################################################

#collide            vss air air.vss

fix		    		in emit/face air xlo nevery 100 n 10000 perspecies no # subsonic 0.1 NULL

timestep 	    	0.0001

#dump                2 image all 50 image.*.ppm type type pdiam 0.1 surf proc 0.01 size 512 512 zoom 1.75 gline yes 0.005
#dump_modify	     2 pad 4

stats		    	10
stats_style	    	step cpu np nattempt ncoll nscoll nscheck
run 		    	500
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00514984 0.00514984 0.00514984
  total     (ave,min,max) = 1.51894 1.51894 1.51894
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
      10   9.8968e-05        0        0        0        0        0 
      20  0.000179916        0        0        0        0        0 
      30  0.000271551        0        0        0        0        0 
      40  0.000348657        0        0        0        0        0 
      50  0.000424367        0        0        0        0        0 
      60  0.000514535        0        0        0        0        0 
      70  0.000594296        0        0        0        0        0 
      80  0.000673428        0        0        0        0        0 
      90  0.000750325        0        0        0        0        0 
     100  0.009425958    10000        0        0        0        0 
     110  0.012356994    10000        0        0        0        8 
     120  0.015390002     9892        0        0       32      672 
     130   0.01894774     9088        0        0      121     1780 
     140  0.022305308     8110        0        0       98     1996 
     150  0.025246332     7223        0        0       84     1945 
     160  0.027886682     6470        0        0       60     1669 
     170  0.030297247     5839        0        0       50     1439 
     180  0.032373753     5335        0        0       38     1134 
     190  0.034233605     4909        0        0       33     1020 
     200  0.041233143    14507        0        0       22      830 
     210  0.044935876    14099        0        0       27      790 
     220   0.04875385    13588        0        0       61     1378 
     230  0.052764381    12488        0        0      110     2209 
     240  0.056519287    11202        0        0      110     2577 
     250  0.059915967    10021        0        0       83     2364 
     260  0.062976353     9031        0        0       86     2106 
     270  0.065557126     8239        0        0       58     1809 
     280  0.067862718     7544        0        0       47     1551 
     290  0.069947046     6906        0        0       43     1267 
     300  0.076389096    16312        0        0       31     1134 
     310  0.079938383    15762        0        0       33      970 
     320  0.083581819    15139        0        0       58     1541 
     330  0.087418232    13951        0        0      120     2552 
     340  0.091163499    12512        0        0      119     2668 
     350  0.094572752    11196        0        0      120     2489 
     360   0.09757405    10114        0        0       85     2198 
     370   0.10027691     9206        0        0       59     1934 
     380   0.10274831     8427        0        0       61     1680 
     390   0.10491365     7726        0        0       51     1437 
     400   0.11388577    17076        0        0       49     1331 
     410    0.1176309    16485        0        0       35     1166 
     420   0.12142345    15824        0        0       65     1682 
     430   0.12541359    14598        0        0      121     2538 
     440   0.12934401    13195        0        0      134     2889 
     450   0.13293339    11863        0        0       94     2653 
     460   0.13615001    10730        0        0       91     2330 
     470   0.13901023     9779        0        0       70     2101 
     480   0.14158848     8945        0        0       58     1841 
     490   0.14391992     8181        0        0       39     1510 
     500   0.15044578    17496        0        0       37     1282 
Loop time of 0.150459 on 1 procs for 500 steps with 17496 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.11608    | 0.11608    | 0.11608    |   0.0 | 77.15
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00176    | 0.00176    | 0.00176    |   0.0 |  1.17
Modify  | 0.029208   | 0.029208   | 0.029208   |   0.0 | 19.41
Output  | 0.00097662 | 0.00097662 | 0.00097662 |   0.0 |  0.65
Other   |            | 0.002437   |            |       |  1.62

Particle moves    = 4219340 (4.22M)
Cells touched     = 4773077 (4.77M)
Particle comms    = 0 (0K)
Boundary collides = 16423 (16.4K)
Boundary exits    = 7176 (7.18K)
SurfColl checks   = 681975 (0.682M)
SurfColl occurs   = 26738 (26.7K)
Surf reactions    = 26665 (26.7K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 2.80432e+07
Particle-moves/step: 8438.68
Cell-touches/particle/step: 1.13124
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00389231
Particle fraction exiting boundary: 0.00170074
Surface-checks/particle/step: 0.161631
Surface-collisions/particle/step: 0.00633701
Surf-reactions/particle/step: 0.00631971
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 26665
    reaction O(g) --> O(s): 19493
    reaction CO2(g) --> C(b) + 2O(g): 3
    reaction O(g) + O(s) --> CO2(g): 110
    reaction O(g) --> CO(s): 6152
    reaction O(g) --> CO(g): 593
    reaction O(g) + O(s) --> O(g) + O(g): 314

Particles: 17496 ave 17496 max 17496 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      400 ave 400 max 400 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0

