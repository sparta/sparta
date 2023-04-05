SPARTA (18 Jul 2022)
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
  CPU time = 0.000903112 secs
  create/ghost percent = 98.2615 1.73854
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 8.6602e-05 secs
  reassign/sort/migrate/ghost percent = 78.4058 0.230942 16.398 4.96524

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
      20     2.68e-05        0        0        0        0        0 
      30     4.64e-05        0        0        0        0        0 
      40   6.5801e-05        0        0        0        0        0 
      50   8.5201e-05        0        0        0        0        0 
      60  0.000104501        0        0        0        0        0 
      70  0.000124201        0        0        0        0        0 
      80  0.000143702        0        0        0        0        0 
      90  0.000163002        0        0        0        0        0 
     100  0.001888525     3150        0        0        0        0 
     110  0.010411339     3150        0        0        0        0 
     120  0.010767044     3150        0        0        0        0 
     130  0.011069848     3150        0        0        0        0 
     140  0.011371452     3150        0        0        0        0 
     150  0.011704157     3150        0        0        0        0 
     160  0.012008961     3150        0        0        0        0 
     170  0.012310765     3150        0        0        0        0 
     180  0.012611969     3150        0        0        0        0 
     190  0.012914573     3150        0        0        0        0 
     200  0.014782998     3230        0        0        0        0 
     210  0.015119602     3204        0        0        0        0 
     220  0.015428806     3204        0        0        0        0 
     230  0.015735511     3204        0        0        0        0 
     240  0.016041515     3204        0        0        0        0 
     250  0.016373719     3204        0        0        0        0 
     260  0.016680423     3204        0        0        0        0 
     270  0.016986927     3204        0        0        0        0 
     280  0.017293531     3204        0        0        0        0 
     290  0.017600836     3204        0        0        0        0 
     300   0.01939176     3299        0        0        0        0 
     310  0.021241384     3272        0        0        0        0 
     320  0.023052809     3272        0        0        0        0 
     330  0.023394113     3271        0        0        0        0 
     340  0.023708317     3270        0        0        0        0 
     350  0.024050622     3268        0        0        0        0 
     360  0.024363726     3267        0        0        0        0 
     370   0.03438846     3266        0        0        0        0 
     380  0.034749865     3264        0        0        0        0 
     390   0.03506487     3257        0        0        0        0 
     400  0.036868294     3379        0        0        0        0 
     410  0.037207098     3346        0        0        0        0 
     420  0.037529703     3345        0        0        0        0 
     430  0.037849307     3343        0        0        0        0 
     440  0.038190611     3339        0        0        0        0 
     450  0.038547516     3330        0        0        0        0 
     460   0.03886642     3322        0        0        0        0 
     470  0.039184425     3315        0        0        0        0 
     480  0.039501129     3312        0        0        0        0 
     490  0.039818633     3306        0        0        0        0 
     500  0.041606457     3384        0        0        0        0 
     510  0.041946062     3343        0        0        0        0 
     520  0.042331467     3337        0        0        0        0 
     530  0.042661071     3328        0        0        0        0 
     540  0.042979176     3318        0        0        0        0 
     550   0.04332418     3312        0        0        0        0 
     560  0.043640984     3302        0        0        0        0 
     570  0.043956589     3292        0        0        0        0 
     580  0.050391975     3285        0        0        0        0 
     590   0.05075498     3278        0        0        0        0 
     600  0.052557704     3356        0        0        0        0 
     610  0.052893508     3317        0        0        0        0 
     620  0.053212713     3311        0        0        0        0 
     630  0.053529317     3308        0        0        0        0 
     640  0.053845121     3303        0        0        0        0 
     650  0.054200826     3297        0        0        0        0 
     660   0.05451933     3288        0        0        0        0 
     670  0.054833834     3280        0        0        0        0 
     680  0.055146639     3274        0        0        0        0 
     690  0.055459643     3268        0        0        0        0 
     700  0.057225666     3374        0        0        0        0 
     710  0.066413389     3340        0        0        0        0 
     720  0.066778594     3335        0        0        0        0 
     730  0.067099199     3326        0        0        0        0 
     740  0.067417703     3320        0        0        0        0 
     750  0.067763008     3314        0        0        0        0 
     760  0.068079812     3310        0        0        0        0 
     770  0.068397216     3303        0        0        0        0 
     780   0.06871262     3295        0        0        0        0 
     790  0.069028224     3286        0        0        0        0 
     800  0.070823249     3406        0        0        0        0 
     810  0.071188553     3368        0        0        0        0 
     820  0.071518458     3356        0        0        0        0 
     830  0.071840662     3351        0        0        0        0 
     840  0.072161266     3329        0        0        0        0 
     850  0.072509271     3321        0        0        0        0 
     860  0.072826975     3319        0        0        0        0 
     870   0.07314438     3306        0        0        0        0 
     880  0.073460284     3299        0        0        0        0 
     890  0.073775688     3294        0        0        0        0 
     900  0.075569612     3385        0        0        0        0 
     910  0.082400004     3350        0        0        0        0 
     920  0.082755008     3341        0        0        0        0 
     930  0.083076613     3333        0        0        0        0 
     940  0.083394617     3326        0        0        0        0 
     950  0.083739322     3318        0        0        0        0 
     960  0.084057126     3310        0        0        0        0 
     970   0.08437303     3299        0        0        0        0 
     980  0.084688834     3293        0        0        0        0 
     990  0.085003438     3286        0        0        0        0 
    1000  0.086797063     3388        0        0        0        0 
Loop time of 0.0868234 on 1 procs for 1000 steps with 3388 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.03428    | 0.03428    | 0.03428    |   0.0 | 39.48
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00033571 | 0.00033571 | 0.00033571 |   0.0 |  0.39
Modify  | 0.0076422  | 0.0076422  | 0.0076422  |   0.0 |  8.80
Output  | 0.044309   | 0.044309   | 0.044309   |   0.0 | 51.03
Other   |            | 0.0002562  |            |       |  0.30

Particle moves    = 2982814 (2.98M)
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

Particle-moves/CPUsec/proc: 3.4355e+07
Particle-moves/step: 2982.81
Cell-touches/particle/step: 1.0145
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.000216239
Particle fraction exiting boundary: 0.00014483
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0.00944075
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
