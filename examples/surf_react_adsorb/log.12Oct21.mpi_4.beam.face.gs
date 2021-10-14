SPARTA (26 Feb 2021)
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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/ascldap/users/stamoor/sparta_stanmoore1/src/grid.cpp:410)
Created 8 child grid cells
  CPU time = 0.00110841 secs
  create/ghost percent = 85.4592 14.5408
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000848293 secs
  reassign/sort/migrate/ghost percent = 71.1917 1.29286 13.7156 13.7999

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
      10 0.00016283989        0        0        0        0        0 
      20 0.00034284592        0        0        0        0        0 
      30 0.00052642822        0        0        0        0        0 
      40  0.000695467        0        0        0        0        0 
      50 0.00087213516        0        0        0        0        0 
      60 0.0010392666        0        0        0        0        0 
      70 0.0012197495        0        0        0        0        0 
      80 0.0013895035        0        0        0        0        0 
      90 0.0015630722        0        0        0        0        0 
     100 0.0043594837     3140        0        0        0        0 
     110 0.0053339005     3140        0        0        0        0 
     120 0.0064661503     3140        0        0        0        0 
     130 0.0078406334     3140        0        0        0        0 
     140 0.0088276863     3140        0        0        0        0 
     150 0.0098791122     3140        0        0        0        0 
     160  0.010854483     3140        0        0        0        0 
     170  0.011837006     3140        0        0        0        0 
     180  0.012796402     3140        0        0        0        0 
     190  0.014796257     3140        0        0        0        0 
     200  0.018142462     3213        0        0        0        0 
     210   0.01931715     3187        0        0        0        0 
     220  0.020368814     3187        0        0        0        0 
     230  0.021421432     3187        0        0        0        0 
     240  0.022422075     3187        0        0        0        0 
     250  0.023502111     3187        0        0        0        0 
     260  0.024490118     3187        0        0        0        0 
     270   0.02804184     3187        0        0        0        0 
     280  0.029119492     3187        0        0        0        0 
     290  0.030567646     3187        0        0        0        0 
     300  0.034149647     3305        0        0        0        0 
     310  0.035222054     3275        0        0        0        0 
     320  0.036554337     3273        0        0        0        0 
     330  0.037599564     3273        0        0        0        0 
     340  0.038662195     3272        0        0        0        0 
     350  0.039768934     3271        0        0        0        0 
     360  0.041181564     3269        0        0        0        0 
     370  0.044560671     3268        0        0        0        0 
     380  0.046038628     3266        0        0        0        0 
     390  0.047752142     3264        0        0        0        0 
     400  0.051796198     3355        0        0        0        0 
     410  0.052966118     3326        0        0        0        0 
     420  0.054037809     3323        0        0        0        0 
     430  0.055087328     3319        0        0        0        0 
     440  0.056118488     3315        0        0        0        0 
     450  0.057226896     3312        0        0        0        0 
     460   0.05867815     3302        0        0        0        0 
     470  0.059844017     3295        0        0        0        0 
     480  0.061305285     3291        0        0        0        0 
     490  0.062709808     3284        0        0        0        0 
     500  0.067785263     3341        0        0        0        0 
     510  0.068876982     3307        0        0        0        0 
     520  0.070078611     3301        0        0        0        0 
     530  0.071543932     3294        0        0        0        0 
     540  0.073207855     3290        0        0        0        0 
     550  0.074692249     3284        0        0        0        0 
     560  0.075720549     3278        0        0        0        0 
     570  0.076719522     3272        0        0        0        0 
     580  0.077718973     3272        0        0        0        0 
     590  0.078721285     3264        0        0        0        0 
     600  0.081975937     3404        0        0        0        0 
     610   0.08306098     3367        0        0        0        0 
     620  0.084128618     3361        0        0        0        0 
     630    0.0873456     3356        0        0        0        0 
     640  0.088941336     3347        0        0        0        0 
     650  0.090639114     3341        0        0        0        0 
     660  0.091988802     3334        0        0        0        0 
     670   0.09500289     3325        0        0        0        0 
     680  0.096183777     3317        0        0        0        0 
     690  0.097597837     3313        0        0        0        0 
     700   0.10284209     3406        0        0        0        0 
     710   0.10393906     3369        0        0        0        0 
     720   0.10500932     3365        0        0        0        0 
     730   0.10605907     3358        0        0        0        0 
     740   0.10712552     3351        0        0        0        0 
     750   0.10842299     3346        0        0        0        0 
     760   0.10958076     3340        0        0        0        0 
     770   0.11059642     3328        0        0        0        0 
     780   0.11179328     3319        0        0        0        0 
     790   0.11281705     3312        0        0        0        0 
     800   0.11602092     3444        0        0        0        0 
     810   0.11719346     3390        0        0        0        0 
     820   0.11828327     3387        0        0        0        0 
     830   0.11932135     3381        0        0        0        0 
     840   0.12042832     3374        0        0        0        0 
     850   0.12210536     3365        0        0        0        0 
     860   0.12314963     3359        0        0        0        0 
     870   0.12435007     3354        0        0        0        0 
     880   0.12823272     3346        0        0        0        0 
     890   0.13008952     3343        0        0        0        0 
     900   0.13347554     3358        0        0        0        0 
     910   0.13455915     3317        0        0        0        0 
     920    0.1355927     3308        0        0        0        0 
     930   0.13665152     3298        0        0        0        0 
     940   0.13769317     3291        0        0        0        0 
     950   0.13879037     3281        0        0        0        0 
     960   0.13980508     3272        0        0        0        0 
     970   0.14083505     3261        0        0        0        0 
     980   0.14181066     3258        0        0        0        0 
     990   0.14284563     3249        0        0        0        0 
    1000   0.14602947     3393        0        0        0        0 
Loop time of 0.146074 on 4 procs for 1000 steps with 3393 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0074108  | 0.041757   | 0.07597    |  16.7 | 28.59
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.019573   | 0.03029    | 0.036      |   3.8 | 20.74
Modify  | 0.00021911 | 0.0090774  | 0.016927   |   8.2 |  6.21
Output  | 0.0029387  | 0.003676   | 0.0054295  |   1.7 |  2.52
Other   |            | 0.06127    |            |       | 41.95

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

Particle-moves/CPUsec/proc: 5.10639e+06
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
