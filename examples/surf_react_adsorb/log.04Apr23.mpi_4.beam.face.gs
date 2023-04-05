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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/ascldap/users/stamoor/sparta_master/src/grid.cpp:465)
Created 8 child grid cells
  CPU time = 0.00111604 secs
  create/ghost percent = 83.401 16.599
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000789165 secs
  reassign/sort/migrate/ghost percent = 72.7795 0.302115 12.0544 14.864

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
      10 0.0001282692        0        0        0        0        0 
      20 0.00030255318        0        0        0        0        0 
      30 0.00048017502        0        0        0        0        0 
      40 0.00065922737        0        0        0        0        0 
      50 0.00083518028        0        0        0        0        0 
      60 0.0010139942        0        0        0        0        0 
      70 0.0011985302        0        0        0        0        0 
      80 0.0013833046        0        0        0        0        0 
      90 0.0015618801        0        0        0        0        0 
     100 0.0036549568     3140        0        0        0        0 
     110 0.0042345524     3140        0        0        0        0 
     120 0.0048434734     3140        0        0        0        0 
     130 0.0054168701     3140        0        0        0        0 
     140 0.0059745312     3140        0        0        0        0 
     150 0.0065860748     3140        0        0        0        0 
     160 0.0071220398     3140        0        0        0        0 
     170 0.0076985359     3140        0        0        0        0 
     180 0.0082728863     3140        0        0        0        0 
     190 0.0097703934     3140        0        0        0        0 
     200    0.0117836     3213        0        0        0        0 
     210  0.012414217     3187        0        0        0        0 
     220  0.013043404     3187        0        0        0        0 
     230  0.013691187     3187        0        0        0        0 
     240  0.014280558     3187        0        0        0        0 
     250  0.014907837     3187        0        0        0        0 
     260  0.015479326     3187        0        0        0        0 
     270  0.016048908     3187        0        0        0        0 
     280  0.016606092     3187        0        0        0        0 
     290   0.01718235     3187        0        0        0        0 
     300   0.01904583     3305        0        0        0        0 
     310  0.019711256     3275        0        0        0        0 
     320  0.020339251     3273        0        0        0        0 
     330  0.020995617     3273        0        0        0        0 
     340  0.021647453     3272        0        0        0        0 
     350  0.022267342     3271        0        0        0        0 
     360  0.022870302     3269        0        0        0        0 
     370  0.023459911     3268        0        0        0        0 
     380  0.024042606     3266        0        0        0        0 
     390   0.02462244     3264        0        0        0        0 
     400  0.026512384     3355        0        0        0        0 
     410  0.027171373     3326        0        0        0        0 
     420  0.027871609     3323        0        0        0        0 
     430  0.028462172     3319        0        0        0        0 
     440  0.029047012     3315        0        0        0        0 
     450  0.029654026     3312        0        0        0        0 
     460  0.030227423     3302        0        0        0        0 
     470  0.030820608     3295        0        0        0        0 
     480   0.03137517     3291        0        0        0        0 
     490   0.03196454     3284        0        0        0        0 
     500  0.033887386     3341        0        0        0        0 
     510  0.034527063     3307        0        0        0        0 
     520  0.035184383     3301        0        0        0        0 
     530  0.035790443     3294        0        0        0        0 
     540  0.036394835     3290        0        0        0        0 
     550  0.037017584     3284        0        0        0        0 
     560  0.037597179     3278        0        0        0        0 
     570  0.038153648     3272        0        0        0        0 
     580  0.038712502     3272        0        0        0        0 
     590  0.039270639     3264        0        0        0        0 
     600   0.04112792     3404        0        0        0        0 
     610  0.041782856     3367        0        0        0        0 
     620  0.042433023     3361        0        0        0        0 
     630  0.043076992     3356        0        0        0        0 
     640  0.043691158     3347        0        0        0        0 
     650  0.044295549     3341        0        0        0        0 
     660  0.044891119     3334        0        0        0        0 
     670   0.04547739     3325        0        0        0        0 
     680   0.04604435     3317        0        0        0        0 
     690  0.046635628     3313        0        0        0        0 
     700  0.048480749     3406        0        0        0        0 
     710  0.049134731     3369        0        0        0        0 
     720  0.049795866     3365        0        0        0        0 
     730  0.050405741     3358        0        0        0        0 
     740  0.051025152     3351        0        0        0        0 
     750  0.051649809     3346        0        0        0        0 
     760  0.052210569     3340        0        0        0        0 
     770  0.052781343     3328        0        0        0        0 
     780  0.053360939     3319        0        0        0        0 
     790  0.053952217     3312        0        0        0        0 
     800  0.055808067     3444        0        0        0        0 
     810  0.056463957     3390        0        0        0        0 
     820  0.057102442     3387        0        0        0        0 
     830  0.057708979     3381        0        0        0        0 
     840  0.058309555     3374        0        0        0        0 
     850  0.058954716     3365        0        0        0        0 
     860  0.059599161     3359        0        0        0        0 
     870  0.060183048     3354        0        0        0        0 
     880  0.060804129     3346        0        0        0        0 
     890  0.061392546     3343        0        0        0        0 
     900  0.063265324     3358        0        0        0        0 
     910  0.063924551     3317        0        0        0        0 
     920  0.064545393     3308        0        0        0        0 
     930  0.065173388     3298        0        0        0        0 
     940   0.06581521     3291        0        0        0        0 
     950  0.066431999     3281        0        0        0        0 
     960   0.06701088     3272        0        0        0        0 
     970  0.067613363     3261        0        0        0        0 
     980  0.068157673     3258        0        0        0        0 
     990  0.068739414     3249        0        0        0        0 
    1000  0.070586443     3393        0        0        0        0 
Loop time of 0.0706483 on 4 procs for 1000 steps with 3393 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0037436  | 0.019215   | 0.035021   |  11.2 | 27.20
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.010218   | 0.010693   | 0.011371   |   0.5 | 15.14
Modify  | 0.00014281 | 0.0051678  | 0.010235   |   7.0 |  7.31
Output  | 0.0019493  | 0.0030513  | 0.0063317  |   3.4 |  4.32
Other   |            | 0.03252    |            |       | 46.03

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

Particle-moves/CPUsec/proc: 1.05581e+07
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
