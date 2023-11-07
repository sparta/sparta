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

boundary	    	oo oo so


create_box          -11 11 -11 11 0 10
Created orthogonal box = (-11 -11 0) to (11 11 10)
create_grid 	    2 2 2
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/runner/work/sparta/sparta/src/grid.cpp:465)
Created 8 child grid cells
  CPU time = 0.0027604 secs
  create/ghost percent = 91.1426 8.85741
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.0011735 secs
  reassign/sort/migrate/ghost percent = 67.2092 0.945888 11.0865 20.7584

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
      10    0.0005699        0        0        0        0        0 
      20  0.001200801        0        0        0        0        0 
      30  0.001798901        0        0        0        0        0 
      40  0.002397101        0        0        0        0        0 
      50  0.002973802        0        0        0        0        0 
      60  0.003560302        0        0        0        0        0 
      70  0.004167903        0        0        0        0        0 
      80  0.004787103        0        0        0        0        0 
      90  0.005393503        0        0        0        0        0 
     100  0.009793606     3140        0        0        0        0 
     110  0.010726607     3140        0        0        0        0 
     120  0.011731008     3140        0        0        0        0 
     130  0.012701708     3140        0        0        0        0 
     140  0.013630109     3140        0        0        0        0 
     150  0.014586109     3140        0        0        0        0 
     160   0.01551461     3140        0        0        0        0 
     170  0.016513911     3140        0        0        0        0 
     180  0.017446911     3140        0        0        0        0 
     190  0.020320413     3140        0        0        0        0 
     200  0.022716915     3213        0        0        0        0 
     210  0.023733015     3187        0        0        0        0 
     220  0.024804616     3187        0        0        0        0 
     230  0.026050517     3187        0        0        0        0 
     240  0.027060718     3187        0        0        0        0 
     250  0.028090918     3187        0        0        0        0 
     260  0.029110219     3187        0        0        0        0 
     270   0.03003672     3187        0        0        0        0 
     280   0.03096572     3187        0        0        0        0 
     290  0.031978721     3187        0        0        0        0 
     300  0.035134123     3305        0        0        0        0 
     310  0.036462524     3275        0        0        0        0 
     320  0.037776325     3273        0        0        0        0 
     330  0.039049226     3273        0        0        0        0 
     340  0.040316226     3272        0        0        0        0 
     350  0.041384227     3271        0        0        0        0 
     360  0.042444828     3269        0        0        0        0 
     370  0.043467028     3268        0        0        0        0 
     380  0.044567029     3266        0        0        0        0 
     390   0.04569223     3264        0        0        0        0 
     400  0.049934633     3355        0        0        0        0 
     410  0.051337634     3326        0        0        0        0 
     420  0.052550534     3323        0        0        0        0 
     430  0.053693235     3319        0        0        0        0 
     440  0.054753536     3315        0        0        0        0 
     450  0.055836537     3312        0        0        0        0 
     460  0.056894937     3302        0        0        0        0 
     470  0.057926838     3295        0        0        0        0 
     480  0.058926539     3291        0        0        0        0 
     490  0.060040239     3284        0        0        0        0 
     500  0.063323841     3341        0        0        0        0 
     510  0.064603842     3307        0        0        0        0 
     520  0.065919743     3301        0        0        0        0 
     530  0.067072444     3294        0        0        0        0 
     540  0.068129545     3290        0        0        0        0 
     550  0.069242345     3284        0        0        0        0 
     560  0.070309146     3278        0        0        0        0 
     570  0.071280247     3272        0        0        0        0 
     580  0.072254947     3272        0        0        0        0 
     590  0.073272648     3264        0        0        0        0 
     600   0.07645225     3404        0        0        0        0 
     610  0.077793451     3367        0        0        0        0 
     620  0.079095352     3361        0        0        0        0 
     630  0.080311153     3356        0        0        0        0 
     640  0.081415553     3347        0        0        0        0 
     650  0.082480154     3341        0        0        0        0 
     660  0.083584455     3334        0        0        0        0 
     670  0.084708155     3325        0        0        0        0 
     680  0.085637556     3317        0        0        0        0 
     690  0.086752457     3313        0        0        0        0 
     700  0.089928759     3406        0        0        0        0 
     710   0.09124986     3369        0        0        0        0 
     720  0.092554861     3365        0        0        0        0 
     730  0.093762861     3358        0        0        0        0 
     740  0.094941262     3351        0        0        0        0 
     750  0.096043363     3346        0        0        0        0 
     760  0.096833263     3340        0        0        0        0 
     770  0.097343864     3328        0        0        0        0 
     780  0.097874964     3319        0        0        0        0 
     790  0.098392364     3312        0        0        0        0 
     800  0.099933065     3444        0        0        0        0 
     810   0.10062387     3390        0        0        0        0 
     820   0.10125467     3387        0        0        0        0 
     830   0.10180437     3381        0        0        0        0 
     840   0.10233067     3374        0        0        0        0 
     850   0.10291807     3365        0        0        0        0 
     860   0.10350487     3359        0        0        0        0 
     870   0.10415827     3354        0        0        0        0 
     880   0.10485187     3346        0        0        0        0 
     890   0.10551647     3343        0        0        0        0 
     900   0.10748787     3358        0        0        0        0 
     910   0.10825567     3317        0        0        0        0 
     920   0.10898567     3308        0        0        0        0 
     930   0.10973117     3298        0        0        0        0 
     940   0.11046017     3291        0        0        0        0 
     950   0.11111947     3281        0        0        0        0 
     960   0.11175667     3272        0        0        0        0 
     970   0.11242357     3261        0        0        0        0 
     980   0.11302647     3258        0        0        0        0 
     990   0.11367737     3249        0        0        0        0 
    1000   0.11566528     3393        0        0        0        0 
Loop time of 0.115715 on 4 procs for 1000 steps with 3393 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0026749  | 0.0093169  | 0.016714   |   6.9 |  8.05
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.034137   | 0.053343   | 0.064234   |   5.2 | 46.10
Modify  | 9.97e-05   | 0.0022262  | 0.0043976  |   4.5 |  1.92
Output  | 0.0052052  | 0.006056   | 0.0072656  |   1.0 |  5.23
Other   |            | 0.04477    |            |       | 38.69

Particle moves    = 2956516 (2.96M)
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

Particle-moves/CPUsec/proc: 6.38751e+06
Particle-moves/step: 2956.52
Cell-touches/particle/step: 1.02385
Particle comm iterations/step: 1.42
Particle fraction communicated: 0.00497647
Particle fraction colliding with boundary: 0.000212412
Particle fraction exiting boundary: 0.000135633
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0.00954874
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
