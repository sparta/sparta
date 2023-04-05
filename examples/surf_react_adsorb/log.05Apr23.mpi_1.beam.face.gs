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
  CPU time = 0.00111209 secs
  create/ghost percent = 96.9057 3.09426
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000151024 secs
  reassign/sort/migrate/ghost percent = 77.8393 0.282737 15.6028 6.27516

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
      10   3.8839e-05        0        0        0        0        0 
      20  0.000107345        0        0        0        0        0 
      30   0.00016431        0        0        0        0        0 
      40  0.000252239        0        0        0        0        0 
      50  0.000301884        0        0        0        0        0 
      60  0.000346114        0        0        0        0        0 
      70  0.000390992        0        0        0        0        0 
      80   0.00043737        0        0        0        0        0 
      90  0.000483916        0        0        0        0        0 
     100  0.003448498     3150        0        0        0        0 
     110  0.004439965     3150        0        0        0        0 
     120  0.005256996     3150        0        0        0        0 
     130   0.00606752     3150        0        0        0        0 
     140  0.006867166     3150        0        0        0        0 
     150  0.007728679     3150        0        0        0        0 
     160  0.008514597     3150        0        0        0        0 
     170  0.009298165     3150        0        0        0        0 
     180  0.010081933     3150        0        0        0        0 
     190  0.010879436     3150        0        0        0        0 
     200  0.015455889     3230        0        0        0        0 
     210  0.016392455     3204        0        0        0        0 
     220  0.017252964     3204        0        0        0        0 
     230  0.018095971     3204        0        0        0        0 
     240  0.018924988     3204        0        0        0        0 
     250  0.019830296     3204        0        0        0        0 
     260  0.020672081     3204        0        0        0        0 
     270  0.021490226     3204        0        0        0        0 
     280  0.022309443     3204        0        0        0        0 
     290  0.023128053     3204        0        0        0        0 
     300  0.027771952     3299        0        0        0        0 
     310  0.028676774     3272        0        0        0        0 
     320  0.029522629     3272        0        0        0        0 
     330  0.030365827     3271        0        0        0        0 
     340  0.031323036     3270        0        0        0        0 
     350  0.032259993     3268        0        0        0        0 
     360  0.033098086     3267        0        0        0        0 
     370  0.033936247     3266        0        0        0        0 
     380  0.034773942     3264        0        0        0        0 
     390  0.035611784     3257        0        0        0        0 
     400  0.040293775     3379        0        0        0        0 
     410  0.041181807     3346        0        0        0        0 
     420  0.042041499     3345        0        0        0        0 
     430  0.042871802     3343        0        0        0        0 
     440  0.043708256     3339        0        0        0        0 
     450  0.044611346     3330        0        0        0        0 
     460  0.045449432     3322        0        0        0        0 
     470  0.046301107     3315        0        0        0        0 
     480   0.04715034     3312        0        0        0        0 
     490   0.04801453     3306        0        0        0        0 
     500  0.052689235     3384        0        0        0        0 
     510  0.053591623     3343        0        0        0        0 
     520  0.054451855     3337        0        0        0        0 
     530  0.055315957     3328        0        0        0        0 
     540  0.056179103     3318        0        0        0        0 
     550  0.057114823     3312        0        0        0        0 
     560  0.057960986     3302        0        0        0        0 
     570  0.058808824     3292        0        0        0        0 
     580  0.059650492     3285        0        0        0        0 
     590  0.060498144     3278        0        0        0        0 
     600  0.065148157     3356        0        0        0        0 
     610  0.066046869     3317        0        0        0        0 
     620   0.06690502     3311        0        0        0        0 
     630  0.067758664     3308        0        0        0        0 
     640  0.068622848     3303        0        0        0        0 
     650   0.06955769     3297        0        0        0        0 
     660  0.070414229     3288        0        0        0        0 
     670  0.071273223     3280        0        0        0        0 
     680  0.072124325     3274        0        0        0        0 
     690  0.072964947     3268        0        0        0        0 
     700   0.07760966     3374        0        0        0        0 
     710   0.07849682     3340        0        0        0        0 
     720  0.079355788     3335        0        0        0        0 
     730  0.080199575     3326        0        0        0        0 
     740  0.081025438     3320        0        0        0        0 
     750  0.081923855     3314        0        0        0        0 
     760  0.082737423     3310        0        0        0        0 
     770  0.083593764     3303        0        0        0        0 
     780  0.084455422     3295        0        0        0        0 
     790  0.085320495     3286        0        0        0        0 
     800  0.089976226     3406        0        0        0        0 
     810    0.0909168     3368        0        0        0        0 
     820  0.091786929     3356        0        0        0        0 
     830  0.092660241     3351        0        0        0        0 
     840  0.093521201     3329        0        0        0        0 
     850  0.094465507     3321        0        0        0        0 
     860  0.095314636     3319        0        0        0        0 
     870  0.096173918     3306        0        0        0        0 
     880  0.097019464     3299        0        0        0        0 
     890  0.097862008     3294        0        0        0        0 
     900   0.10252256     3385        0        0        0        0 
     910   0.10341403     3350        0        0        0        0 
     920   0.10429286     3341        0        0        0        0 
     930   0.10515587     3333        0        0        0        0 
     940   0.10601495     3326        0        0        0        0 
     950   0.10695523     3318        0        0        0        0 
     960    0.1078044     3310        0        0        0        0 
     970   0.10865975     3299        0        0        0        0 
     980   0.10950666     3293        0        0        0        0 
     990   0.11035072     3286        0        0        0        0 
    1000   0.11499003     3388        0        0        0        0 
Loop time of 0.115 on 1 procs for 1000 steps with 3388 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.093949   | 0.093949   | 0.093949   |   0.0 | 81.69
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0010182  | 0.0010182  | 0.0010182  |   0.0 |  0.89
Modify  | 0.017725   | 0.017725   | 0.017725   |   0.0 | 15.41
Output  | 0.0014011  | 0.0014011  | 0.0014011  |   0.0 |  1.22
Other   |            | 0.0009068  |            |       |  0.79

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

Particle-moves/CPUsec/proc: 2.59374e+07
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
