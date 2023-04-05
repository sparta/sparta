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

boundary	    	oo oo so


create_box          -11 11 -11 11 0 10
Created orthogonal box = (-11 -11 0) to (11 11 10)
create_grid 	    2 2 2
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/runner/work/sparta/sparta/src/grid.cpp:465)
Created 8 child grid cells
  CPU time = 0.00200533 secs
  create/ghost percent = 93.3925 6.6075
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000639808 secs
  reassign/sort/migrate/ghost percent = 68.0837 0.828374 10.0813 21.0066

global		    	nrho 1e10 fnum 1e6

species		    	air.species O CO CO2 O2 C
mixture		    	air O O2 vstream 0 1000 -1000

mixture             air O   frac 1.0
mixture             air CO  frac 0.0
mixture             air CO2 frac 0.0
mixture             air C   frac 0.0
mixture 			air O2 frac 0.0


surf_collide        1 cll 300.0 0.5 0.5 0.5 0.5

bound_modify 		zlo collide 1

##################################### SURF REACT ADSORB ######################################
##################################### FACE/BOUNDARY OPTION ###################################

#surf_react        	adsorb_test_gs1 adsorb gs sample-GS_1.surf nsync 1 face 1000 6.022e18 O CO
#bound_modify        zlo react adsorb_test_gs1


surf_react        	adsorb_test_gs2 adsorb gs sample-GS_2.surf nsync 1 face 1000 6.022e18 O CO
bound_modify        zlo react adsorb_test_gs2

########################## BEAM ############################################################
# Beam at multiple points so that different processors handle the surface collisions

region              circle1 cylinder z  0 -10 1 -INF INF

fix                 in1 emit/face/file air zhi data.beam beam_area_1 nevery 100 region circle1

################################################################################################

#dump                2 image all 10 image.*.ppm type type pdiam 0.2 surf proc 0.01 size 512 512 zoom 1.75 gline no 0.005
#dump_modify	     2 pad 4

timestep            0.0001

stats		    	10
stats_style	    	step cpu np nattempt ncoll nscoll nscheck
run 		    	1000
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0 0 0
  total     (ave,min,max) = 1.51379 1.51379 1.51379
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
      10  0.000261404        0        0        0        0        0 
      20  0.000550207        0        0        0        0        0 
      30  0.000831011        0        0        0        0        0 
      40  0.001136515        0        0        0        0        0 
      50  0.001435619        0        0        0        0        0 
      60  0.001729822        0        0        0        0        0 
      70  0.002048926        0        0        0        0        0 
      80   0.00233163        0        0        0        0        0 
      90  0.002615533        0        0        0        0        0 
     100  0.005770173     3140        0        0        0        0 
     110  0.006355181     3140        0        0        0        0 
     120  0.006962389     3140        0        0        0        0 
     130  0.007542096     3140        0        0        0        0 
     140  0.008099003     3140        0        0        0        0 
     150  0.008689611     3140        0        0        0        0 
     160  0.009274618     3140        0        0        0        0 
     170  0.009856125     3140        0        0        0        0 
     180  0.010430833     3140        0        0        0        0 
     190  0.013107167     3140        0        0        0        0 
     200  0.015177293     3213        0        0        0        0 
     210  0.015899002     3187        0        0        0        0 
     220  0.016618911     3187        0        0        0        0 
     230   0.01736142     3187        0        0        0        0 
     240  0.017968328     3187        0        0        0        0 
     250  0.018600136     3187        0        0        0        0 
     260  0.019197344     3187        0        0        0        0 
     270  0.019750351     3187        0        0        0        0 
     280  0.020305258     3187        0        0        0        0 
     290  0.020903065     3187        0        0        0        0 
     300  0.022887591     3305        0        0        0        0 
     310     0.023647     3275        0        0        0        0 
     320  0.024372509     3273        0        0        0        0 
     330  0.025094319     3273        0        0        0        0 
     340  0.025807928     3272        0        0        0        0 
     350  0.026437736     3271        0        0        0        0 
     360  0.027054243     3269        0        0        0        0 
     370  0.027648051     3268        0        0        0        0 
     380  0.028293059     3266        0        0        0        0 
     390  0.028920867     3264        0        0        0        0 
     400  0.030915092     3355        0        0        0        0 
     410  0.031694402     3326        0        0        0        0 
     420  0.032393811     3323        0        0        0        0 
     430   0.03307322     3319        0        0        0        0 
     440  0.033704828     3315        0        0        0        0 
     450  0.034348936     3312        0        0        0        0 
     460  0.034958044     3302        0        0        0        0 
     470  0.035572051     3295        0        0        0        0 
     480  0.036156459     3291        0        0        0        0 
     490  0.036788467     3284        0        0        0        0 
     500  0.038824393     3341        0        0        0        0 
     510  0.039582302     3307        0        0        0        0 
     520  0.040349412     3301        0        0        0        0 
     530  0.041039621     3294        0        0        0        0 
     540  0.041660629     3290        0        0        0        0 
     550  0.042309237     3284        0        0        0        0 
     560  0.042931245     3278        0        0        0        0 
     570  0.043504352     3272        0        0        0        0 
     580  0.044082559     3272        0        0        0        0 
     590  0.044686667     3264        0        0        0        0 
     600  0.046637992     3404        0        0        0        0 
     610  0.047403102     3367        0        0        0        0 
     620  0.048151111     3361        0        0        0        0 
     630   0.04886972     3356        0        0        0        0 
     640  0.049530328     3347        0        0        0        0 
     650  0.050156236     3341        0        0        0        0 
     660  0.050809645     3334        0        0        0        0 
     670  0.051479653     3325        0        0        0        0 
     680   0.05204316     3317        0        0        0        0 
     690  0.052686569     3313        0        0        0        0 
     700  0.054663494     3406        0        0        0        0 
     710  0.055442603     3369        0        0        0        0 
     720  0.056190113     3365        0        0        0        0 
     730  0.056869122     3358        0        0        0        0 
     740   0.05756603     3351        0        0        0        0 
     750  0.058223339     3346        0        0        0        0 
     760  0.058810146     3340        0        0        0        0 
     770  0.059413954     3328        0        0        0        0 
     780  0.060036762     3319        0        0        0        0 
     790   0.06067187     3312        0        0        0        0 
     800  0.062656295     3444        0        0        0        0 
     810  0.063458005     3390        0        0        0        0 
     820  0.064189314     3387        0        0        0        0 
     830  0.064837723     3381        0        0        0        0 
     840  0.065479531     3374        0        0        0        0 
     850   0.06618404     3365        0        0        0        0 
     860  0.066849448     3359        0        0        0        0 
     870  0.067467756     3354        0        0        0        0 
     880  0.068115364     3346        0        0        0        0 
     890  0.068753272     3343        0        0        0        0 
     900  0.070736397     3358        0        0        0        0 
     910  0.071513107     3317        0        0        0        0 
     920  0.072213116     3308        0        0        0        0 
     930  0.072930125     3298        0        0        0        0 
     940  0.073634234     3291        0        0        0        0 
     950  0.074288042     3281        0        0        0        0 
     960  0.074929351     3272        0        0        0        0 
     970  0.075563659     3261        0        0        0        0 
     980  0.076127466     3258        0        0        0        0 
     990  0.076746774     3249        0        0        0        0 
    1000  0.078708999     3393        0        0        0        0 
Loop time of 0.0787771 on 4 procs for 1000 steps with 3393 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0025732  | 0.0088131  | 0.015049   |   6.6 | 11.19
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.020354   | 0.022242   | 0.025124   |   1.3 | 28.23
Modify  | 7.58e-05   | 0.0022626  | 0.0044531  |   4.6 |  2.87
Output  | 0.0025521  | 0.0043178  | 0.006716   |   2.8 |  5.48
Other   |            | 0.04114    |            |       | 52.23

Particle moves    = 2983638 (2.98M)
Cells touched     = 3027028 (3.03M)
Particle comms    = 14713 (14.7K)
Boundary collides = 628 (0.628K)
Boundary exits    = 401 (0.401K)
SurfColl checks   = 0 (0K)
SurfColl occurs   = 0 (0K)
Surf reactions    = 28231 (28.2K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 9.4686e+06
Particle-moves/step: 2983.64
Cell-touches/particle/step: 1.01454
Particle comm iterations/step: 1.42
Particle fraction communicated: 0.00493123
Particle fraction colliding with boundary: 0.000210481
Particle fraction exiting boundary: 0.0001344
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0.00946194
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 28231
    reaction O(g) --> O(s): 20944
    reaction O(g) + O(s) --> CO2(g): 4
    reaction O(g) --> CO(s): 6659
    reaction O(g) --> CO(g): 616
    reaction O(g) + O(s) --> O(g) + O(g): 8

Particles: 848.25 ave 1642 max 62 min
Histogram: 2 0 0 0 0 0 0 0 0 2
Cells:      2 ave 2 max 2 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 6 ave 6 max 6 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 6 ave 6 max 6 min
Histogram: 4 0 0 0 0 0 0 0 0 0
