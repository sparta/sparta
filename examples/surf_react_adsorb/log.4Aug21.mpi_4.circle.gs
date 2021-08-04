SPARTA (26 Feb 2021)
# 2d flow around a circle

seed	    	    12345
dimension   	    2
global              gridcut 0.0 comm/sort yes

boundary	    	o r p

create_box  	    0 10 0 10 -0.5 0.5
Created orthogonal box = (0 0 -0.5) to (10 10 0.5)
create_grid 	    20 20 1
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (../grid.cpp:410)
Created 400 child grid cells
  CPU time = 0.00170697 secs
  create/ghost percent = 93.4616 6.53844
balance_grid        rcb cell
Balance grid migrated 280 cells
  CPU time = 0.00118719 secs
  reassign/sort/migrate/ghost percent = 59.1776 0.411812 25.5265 14.8841

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
  CPU time = 0.00122358 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 30.3099 10.7255 0.964709 44.8028 13.1971 14.0933 0.348158
  surf2grid time = 0.000548199 secs
  map/comm1/comm2/comm3/comm4/split percent = 35.3549 12.2563 8.76543 4.54835 16.002 15.1866

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
      10  0.000508947        0        0        0        0        0 
      20  0.000925142        0        0        0        0        0 
      30  0.001356982        0        0        0        0        0 
      40  0.001755717        0        0        0        0        0 
      50  0.002146349        0        0        0        0        0 
      60  0.002550811        0        0        0        0        0 
      70  0.002938649        0        0        0        0        0 
      80  0.003348489        0        0        0        0        0 
      90  0.003737166        0        0        0        0        0 
     100  0.009024499    10000        0        0        0        0 
     110   0.01120584    10000        0        0        0        6 
     120  0.013557737     9903        0        0       32      683 
     130  0.017422157     9169        0        0       99     1622 
     140  0.021213171     8207        0        0      112     2015 
     150  0.023410925     7311        0        0       79     1862 
     160  0.025347535     6595        0        0       60     1621 
     170  0.027067422     5992        0        0       65     1470 
     180  0.028591329     5511        0        0       40     1243 
     190  0.029950685     5058        0        0       39     1065 
     200  0.033934467    14639        0        0       34      952 
     210  0.036404608    14201        0        0       31      831 
     220  0.038927272    13673        0        0       61     1382 
     230  0.041555888    12626        0        0      104     2336 
     240  0.044019604    11291        0        0      127     2593 
     250  0.046320585    10073        0        0       86     2387 
     260  0.048378094     9110        0        0       77     2000 
     270  0.050192548     8285        0        0       64     1766 
     280  0.051838541     7581        0        0       49     1485 
     290  0.053288344     6934        0        0       37     1234 
     300  0.057028094    16355        0        0       40     1171 
     310  0.059434609    15817        0        0       35      996 
     320  0.061831066    15196        0        0       75     1641 
     330  0.064281442    13971        0        0      107     2434 
     340   0.06663264    12593        0        0      110     2719 
     350  0.068790024    11263        0        0      102     2513 
     360  0.070723142    10182        0        0       98     2251 
     370  0.072481652     9315        0        0       56     2023 
     380  0.074024766     8491        0        0       61     1714 
     390  0.075391736     7798        0        0       43     1419 
     400  0.078827947    17138        0        0       36     1285 
     410  0.081045536    16539        0        0       39     1170 
     420  0.083294834    15863        0        0       62     1570 
     430  0.085650991    14641        0        0      122     2685 
     440  0.087969294    13217        0        0      116     2843 
     450  0.090071223    11921        0        0       99     2809 
     460  0.091989953    10780        0        0       92     2418 
     470  0.093693147     9819        0        0       62     1947 
     480   0.09523668     9002        0        0       46     1696 
     490  0.096620691     8250        0        0       43     1544 
     500    0.1000965    17512        0        0       43     1352 
Loop time of 0.10012 on 4 procs for 500 steps with 17512 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.011912   | 0.032569   | 0.057006   |  10.8 | 32.53
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.011736   | 0.012385   | 0.012647   |   0.3 | 12.37
Modify  | 0.00082763 | 0.0072985  | 0.014752   |   7.6 |  7.29
Output  | 0.0018659  | 0.0020219  | 0.0024452  |   0.5 |  2.02
Other   |            | 0.04585    |            |       | 45.79

Particle moves    = 4247271 (4.25M)
Cells touched     = 4801983 (4.8M)
Particle comms    = 17795 (17.8K)
Boundary collides = 16318 (16.3K)
Boundary exits    = 7264 (7.26K)
SurfColl checks   = 683350 (0.683M)
SurfColl occurs   = 26641 (26.6K)
Surf reactions    = 26565 (26.6K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 1.06055e+07
Particle-moves/step: 8494.54
Cell-touches/particle/step: 1.1306
Particle comm iterations/step: 1.804
Particle fraction communicated: 0.00418975
Particle fraction colliding with boundary: 0.003842
Particle fraction exiting boundary: 0.00171027
Surface-checks/particle/step: 0.160892
Surface-collisions/particle/step: 0.0062725
Surf-reactions/particle/step: 0.0062546
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 26565
    reaction O(g) --> O(s): 19542
    reaction CO2(g) --> C(b) + 2O(g): 1
    reaction O(g) + O(s) --> CO2(g): 94
    reaction O(g) --> CO(s): 6018
    reaction O(g) --> CO(g): 575
    reaction O(g) + O(s) --> O(g) + O(g): 335

Particles: 4378 ave 7311 max 1437 min
Histogram: 2 0 0 0 0 0 0 0 0 2
Cells:      100 ave 100 max 100 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 21 ave 21 max 21 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 21 ave 21 max 21 min
Histogram: 4 0 0 0 0 0 0 0 0 0
Surfs:    50 ave 50 max 50 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0

