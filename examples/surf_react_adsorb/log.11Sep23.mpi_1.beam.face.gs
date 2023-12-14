SPARTA (13 Apr 2023)
Running on 1 MPI task(s)
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
Created 8 child grid cells
  CPU time = 0.000824202 secs
  create/ghost percent = 98.1558 1.84421
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 8.06e-05 secs
  reassign/sort/migrate/ghost percent = 78.0397 0.248139 16.3772 5.33499

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
      10      6.5e-06        0        0        0        0        0 
      20     2.59e-05        0        0        0        0        0 
      30     4.45e-05        0        0        0        0        0 
      40      6.3e-05        0        0        0        0        0 
      50     8.15e-05        0        0        0        0        0 
      60     9.99e-05        0        0        0        0        0 
      70    0.0001184        0        0        0        0        0 
      80    0.0001371        0        0        0        0        0 
      90    0.0001556        0        0        0        0        0 
     100  0.001853604     3150        0        0        0        0 
     110  0.006447714     3150        0        0        0        0 
     120   0.01446033     3150        0        0        0        0 
     130  0.014801531     3150        0        0        0        0 
     140  0.015110132     3150        0        0        0        0 
     150  0.015438932     3150        0        0        0        0 
     160  0.015746633     3150        0        0        0        0 
     170  0.016049034     3150        0        0        0        0 
     180  0.016351534     3150        0        0        0        0 
     190  0.016659335     3150        0        0        0        0 
     200  0.018635039     3230        0        0        0        0 
     210  0.026457155     3204        0        0        0        0 
     220  0.026798556     3204        0        0        0        0 
     230  0.027109757     3204        0        0        0        0 
     240  0.027421957     3204        0        0        0        0 
     250  0.027755358     3204        0        0        0        0 
     260  0.028069859     3204        0        0        0        0 
     270  0.028379959     3204        0        0        0        0 
     280   0.02869016     3204        0        0        0        0 
     290  0.028999361     3204        0        0        0        0 
     300  0.030921165     3299        0        0        0        0 
     310  0.031252265     3272        0        0        0        0 
     320  0.031568866     3272        0        0        0        0 
     330  0.031886267     3271        0        0        0        0 
     340  0.032204767     3270        0        0        0        0 
     350  0.032551368     3268        0        0        0        0 
     360  0.032867769     3267        0        0        0        0 
     370   0.03318667     3266        0        0        0        0 
     380   0.03350697     3264        0        0        0        0 
     390  0.042504789     3257        0        0        0        0 
     400  0.044470893     3379        0        0        0        0 
     410  0.044809294     3346        0        0        0        0 
     420  0.045133295     3345        0        0        0        0 
     430  0.045458895     3343        0        0        0        0 
     440  0.045782896     3339        0        0        0        0 
     450  0.046141797     3330        0        0        0        0 
     460  0.046465297     3322        0        0        0        0 
     470  0.046785298     3315        0        0        0        0 
     480  0.047103499     3312        0        0        0        0 
     490  0.047427399     3306        0        0        0        0 
     500  0.049336903     3384        0        0        0        0 
     510  0.049671804     3343        0        0        0        0 
     520  0.050009405     3337        0        0        0        0 
     530  0.050341405     3328        0        0        0        0 
     540  0.050666406     3318        0        0        0        0 
     550  0.051012807     3312        0        0        0        0 
     560  0.058534723     3302        0        0        0        0 
     570  0.058889423     3292        0        0        0        0 
     580  0.059207024     3285        0        0        0        0 
     590  0.059525025     3278        0        0        0        0 
     600  0.061415129     3356        0        0        0        0 
     610  0.061746629     3317        0        0        0        0 
     620   0.06206613     3311        0        0        0        0 
     630  0.062394031     3308        0        0        0        0 
     640  0.062716831     3303        0        0        0        0 
     650  0.063065232     3297        0        0        0        0 
     660  0.063383533     3288        0        0        0        0 
     670  0.063701833     3280        0        0        0        0 
     680  0.064023734     3274        0        0        0        0 
     690  0.064341035     3268        0        0        0        0 
     700  0.066226739     3374        0        0        0        0 
     710  0.074477956     3340        0        0        0        0 
     720  0.074831657     3335        0        0        0        0 
     730  0.075152857     3326        0        0        0        0 
     740  0.075472158     3320        0        0        0        0 
     750  0.075819459     3314        0        0        0        0 
     760  0.076137559     3310        0        0        0        0 
     770   0.07645746     3303        0        0        0        0 
     780  0.076776961     3295        0        0        0        0 
     790  0.077094961     3286        0        0        0        0 
     800  0.078995565     3406        0        0        0        0 
     810  0.079349466     3368        0        0        0        0 
     820  0.079675167     3356        0        0        0        0 
     830  0.079996467     3351        0        0        0        0 
     840  0.080318768     3329        0        0        0        0 
     850  0.080671669     3321        0        0        0        0 
     860   0.08099157     3319        0        0        0        0 
     870   0.08131047     3306        0        0        0        0 
     880  0.090460689     3299        0        0        0        0 
     890   0.09081429     3294        0        0        0        0 
     900  0.092723294     3385        0        0        0        0 
     910  0.093055795     3350        0        0        0        0 
     920  0.093378695     3341        0        0        0        0 
     930  0.093704496     3333        0        0        0        0 
     940  0.094029297     3326        0        0        0        0 
     950  0.094386598     3318        0        0        0        0 
     960  0.094708598     3310        0        0        0        0 
     970  0.095028599     3299        0        0        0        0 
     980    0.0953473     3293        0        0        0        0 
     990     0.095664     3286        0        0        0        0 
    1000  0.097552104     3388        0        0        0        0 
Loop time of 0.0975686 on 1 procs for 1000 steps with 3388 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.034708   | 0.034708   | 0.034708   |   0.0 | 35.57
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0003552  | 0.0003552  | 0.0003552  |   0.0 |  0.36
Modify  | 0.0084601  | 0.0084601  | 0.0084601  |   0.0 |  8.67
Output  | 0.053829   | 0.053829   | 0.053829   |   0.0 | 55.17
Other   |            | 0.0002159  |            |       |  0.22

Particle moves    = 2982819 (2.98M)
Cells touched     = 3026072 (3.03M)
Particle comms    = 0 (0K)
Boundary collides = 645 (0.645K)
Boundary exits    = 432 (0.432K)
SurfColl checks   = 0 (0K)
SurfColl occurs   = 0 (0K)
Surf reactions    = 28160 (28.2K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 3.05715e+07
Particle-moves/step: 2982.82
Cell-touches/particle/step: 1.0145
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.000216238
Particle fraction exiting boundary: 0.000144829
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0.00944073
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 28160
    reaction O(g) --> O(s): 20986
    reaction O(g) + O(s) --> CO2(g): 2
    reaction O(g) --> CO(s): 6529
    reaction O(g) --> CO(g): 638
    reaction O(g) + O(s) --> O(g) + O(g): 5

Particles: 3388 ave 3388 max 3388 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      8 ave 8 max 8 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
