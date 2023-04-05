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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/me/sparta_master/src/grid.cpp:465)
Created 8 child grid cells
  CPU time = 0.0015072 secs
  create/ghost percent = 85.1177 14.8823
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000620257 secs
  reassign/sort/migrate/ghost percent = 70.9572 0.447879 11.1723 17.4226

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
      10  0.000202313        0        0        0        0        0 
      20  0.000444654        0        0        0        0        0 
      30  0.000670498        0        0        0        0        0 
      40  0.000890838        0        0        0        0        0 
      50  0.001117971        0        0        0        0        0 
      60  0.001334449        0        0        0        0        0 
      70  0.001546081        0        0        0        0        0 
      80  0.001763748        0        0        0        0        0 
      90  0.001977982        0        0        0        0        0 
     100  0.004785998     3140        0        0        0        0 
     110  0.005459476     3140        0        0        0        0 
     120  0.006103613     3140        0        0        0        0 
     130   0.00678383     3140        0        0        0        0 
     140  0.007391889     3140        0        0        0        0 
     150  0.008046971     3140        0        0        0        0 
     160   0.00869452     3140        0        0        0        0 
     170  0.009620828     3140        0        0        0        0 
     180  0.010225948     3140        0        0        0        0 
     190  0.012952101     3140        0        0        0        0 
     200   0.01553588     3213        0        0        0        0 
     210  0.016238289     3187        0        0        0        0 
     220  0.016961931     3187        0        0        0        0 
     230  0.017660298     3187        0        0        0        0 
     240  0.018281154     3187        0        0        0        0 
     250   0.01895076     3187        0        0        0        0 
     260  0.019571232     3187        0        0        0        0 
     270  0.020167229     3187        0        0        0        0 
     280  0.020794408     3187        0        0        0        0 
     290  0.021421152     3187        0        0        0        0 
     300  0.023872978     3305        0        0        0        0 
     310  0.024609012     3275        0        0        0        0 
     320  0.025300672     3273        0        0        0        0 
     330  0.025993492     3273        0        0        0        0 
     340  0.026701645     3272        0        0        0        0 
     350  0.027371995     3271        0        0        0        0 
     360  0.028010366     3269        0        0        0        0 
     370  0.028654565     3268        0        0        0        0 
     380  0.029300833     3266        0        0        0        0 
     390  0.029937729     3264        0        0        0        0 
     400  0.032390774     3355        0        0        0        0 
     410   0.03311553     3326        0        0        0        0 
     420  0.033807506     3323        0        0        0        0 
     430  0.034477462     3319        0        0        0        0 
     440  0.035117768     3315        0        0        0        0 
     450   0.03579372     3312        0        0        0        0 
     460  0.036454624     3302        0        0        0        0 
     470  0.037099505     3295        0        0        0        0 
     480  0.037721095     3291        0        0        0        0 
     490  0.038365695     3284        0        0        0        0 
     500  0.040857614     3341        0        0        0        0 
     510  0.041612652     3307        0        0        0        0 
     520  0.042330878     3301        0        0        0        0 
     530  0.043004675     3294        0        0        0        0 
     540  0.043651717     3290        0        0        0        0 
     550  0.044351474     3284        0        0        0        0 
     560   0.04499986     3278        0        0        0        0 
     570  0.045620679     3272        0        0        0        0 
     580  0.046245815     3272        0        0        0        0 
     590  0.046885275     3264        0        0        0        0 
     600  0.049274807     3404        0        0        0        0 
     610  0.050003278     3367        0        0        0        0 
     620  0.050720769     3361        0        0        0        0 
     630  0.051412476     3356        0        0        0        0 
     640  0.052079654     3347        0        0        0        0 
     650  0.052779177     3341        0        0        0        0 
     660   0.05344105     3334        0        0        0        0 
     670  0.054117876     3325        0        0        0        0 
     680  0.054734788     3317        0        0        0        0 
     690  0.055396352     3313        0        0        0        0 
     700  0.058090442     3406        0        0        0        0 
     710  0.058827685     3369        0        0        0        0 
     720  0.059546566     3365        0        0        0        0 
     730  0.060231644     3358        0        0        0        0 
     740  0.060920091     3351        0        0        0        0 
     750  0.061623212     3346        0        0        0        0 
     760  0.062254736     3340        0        0        0        0 
     770   0.06289458     3328        0        0        0        0 
     780  0.063544845     3319        0        0        0        0 
     790  0.064219018     3312        0        0        0        0 
     800  0.066643604     3444        0        0        0        0 
     810  0.067397328     3390        0        0        0        0 
     820   0.06811536     3387        0        0        0        0 
     830   0.06878663     3381        0        0        0        0 
     840  0.069438472     3374        0        0        0        0 
     850  0.070156799     3365        0        0        0        0 
     860  0.072639152     3359        0        0        0        0 
     870  0.073302526     3354        0        0        0        0 
     880  0.073972621     3346        0        0        0        0 
     890  0.074629408     3343        0        0        0        0 
     900  0.077102705     3358        0        0        0        0 
     910  0.077840627     3317        0        0        0        0 
     920  0.078537934     3308        0        0        0        0 
     930  0.079235078     3298        0        0        0        0 
     940  0.079922797     3291        0        0        0        0 
     950  0.080615606     3281        0        0        0        0 
     960  0.081265662     3272        0        0        0        0 
     970  0.081922096     3261        0        0        0        0 
     980  0.082535428     3258        0        0        0        0 
     990  0.083184113     3249        0        0        0        0 
    1000  0.085607631     3393        0        0        0        0 
Loop time of 0.0856534 on 4 procs for 1000 steps with 3393 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0067121  | 0.02533    | 0.043885   |  11.6 | 29.57
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.014715   | 0.016016   | 0.017793   |   0.9 | 18.70
Modify  | 0.00021503 | 0.0054653  | 0.010751   |   7.1 |  6.38
Output  | 0.0018653  | 0.0023662  | 0.0037196  |   1.6 |  2.76
Other   |            | 0.03648    |            |       | 42.58

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

Particle-moves/CPUsec/proc: 8.70846e+06
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
