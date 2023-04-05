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
  CPU time = 0.000915051 secs
  create/ghost percent = 87.4935 12.5065
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000128031 secs
  reassign/sort/migrate/ghost percent = 64.9907 0.18622 27.3743 7.44879

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
      10 1.8358231e-05        0        0        0        0        0 
      20 4.0054321e-05        0        0        0        0        0 
      30 5.6743622e-05        0        0        0        0        0 
      40 7.4148178e-05        0        0        0        0        0 
      50 9.059906e-05        0        0        0        0        0 
      60 0.00010824203        0        0        0        0        0 
      70 0.00012564659        0        0        0        0        0 
      80 0.00014233589        0        0        0        0        0 
      90 0.00015950203        0        0        0        0        0 
     100 0.0028719902     3150        0        0        0        0 
     110 0.0035755634     3150        0        0        0        0 
     120 0.0042796135     3150        0        0        0        0 
     130 0.0049843788     3150        0        0        0        0 
     140  0.005689621     3150        0        0        0        0 
     150 0.0064325333     3150        0        0        0        0 
     160 0.0071387291     3150        0        0        0        0 
     170 0.0078425407     3150        0        0        0        0 
     180 0.0085363388     3150        0        0        0        0 
     190 0.0092360973     3150        0        0        0        0 
     200  0.012667656     3230        0        0        0        0 
     210  0.013384819     3204        0        0        0        0 
     220  0.014101505     3204        0        0        0        0 
     230  0.014817476     3204        0        0        0        0 
     240  0.015523672     3204        0        0        0        0 
     250  0.016319036     3204        0        0        0        0 
     260   0.01706171     3204        0        0        0        0 
     270  0.017777681     3204        0        0        0        0 
     280  0.018507481     3204        0        0        0        0 
     290  0.019219875     3204        0        0        0        0 
     300  0.022522926     3299        0        0        0        0 
     310  0.023283005     3272        0        0        0        0 
     320  0.024036169     3272        0        0        0        0 
     330  0.024790049     3271        0        0        0        0 
     340  0.025533438     3270        0        0        0        0 
     350  0.026309013     3268        0        0        0        0 
     360  0.027038097     3267        0        0        0        0 
     370   0.02776742     3266        0        0        0        0 
     380  0.028486729     3264        0        0        0        0 
     390  0.029237509     3257        0        0        0        0 
     400  0.032561302     3379        0        0        0        0 
     410  0.033308268     3346        0        0        0        0 
     420  0.034054756     3345        0        0        0        0 
     430  0.034800053     3343        0        0        0        0 
     440  0.035536289     3339        0        0        0        0 
     450  0.036326885     3330        0        0        0        0 
     460  0.037068605     3322        0        0        0        0 
     470  0.037808418     3315        0        0        0        0 
     480  0.038537979     3312        0        0        0        0 
     490  0.039273024     3306        0        0        0        0 
     500  0.042605877     3384        0        0        0        0 
     510  0.043354511     3343        0        0        0        0 
     520  0.044100046     3337        0        0        0        0 
     530  0.044843912     3328        0        0        0        0 
     540  0.045583963     3318        0        0        0        0 
     550  0.046362638     3312        0        0        0        0 
     560  0.047100306     3302        0        0        0        0 
     570  0.047835827     3292        0        0        0        0 
     580  0.048568487     3285        0        0        0        0 
     590  0.049290895     3278        0        0        0        0 
     600   0.05257988     3356        0        0        0        0 
     610  0.053321838     3317        0        0        0        0 
     620  0.054061413     3311        0        0        0        0 
     630  0.054807425     3308        0        0        0        0 
     640  0.055564404     3303        0        0        0        0 
     650  0.056362867     3297        0        0        0        0 
     660  0.057118416     3288        0        0        0        0 
     670  0.057872772     3280        0        0        0        0 
     680  0.058624268     3274        0        0        0        0 
     690  0.059361696     3268        0        0        0        0 
     700   0.06280756     3374        0        0        0        0 
     710  0.063597679     3340        0        0        0        0 
     720  0.064332485     3335        0        0        0        0 
     730  0.065111637     3326        0        0        0        0 
     740   0.06585288     3320        0        0        0        0 
     750  0.066693068     3314        0        0        0        0 
     760  0.067425966     3310        0        0        0        0 
     770  0.068161726     3303        0        0        0        0 
     780  0.068920612     3295        0        0        0        0 
     790   0.06965661     3286        0        0        0        0 
     800  0.072989225     3406        0        0        0        0 
     810  0.073752165     3368        0        0        0        0 
     820  0.074515343     3356        0        0        0        0 
     830  0.075281858     3351        0        0        0        0 
     840  0.076028824     3329        0        0        0        0 
     850   0.07682085     3321        0        0        0        0 
     860  0.077560425     3319        0        0        0        0 
     870  0.078289509     3306        0        0        0        0 
     880  0.079025745     3299        0        0        0        0 
     890  0.079761267     3294        0        0        0        0 
     900  0.083076239     3385        0        0        0        0 
     910  0.083833933     3350        0        0        0        0 
     920  0.084578753     3341        0        0        0        0 
     930  0.085312366     3333        0        0        0        0 
     940  0.086053371     3326        0        0        0        0 
     950  0.086843967     3318        0        0        0        0 
     960   0.08758235     3310        0        0        0        0 
     970  0.088308811     3299        0        0        0        0 
     980  0.089043856     3293        0        0        0        0 
     990  0.089777946     3286        0        0        0        0 
    1000  0.093062639     3388        0        0        0        0 
Loop time of 0.0930734 on 1 procs for 1000 steps with 3388 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.071496   | 0.071496   | 0.071496   |   0.0 | 76.82
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00089931 | 0.00089931 | 0.00089931 |   0.0 |  0.97
Modify  | 0.018113   | 0.018113   | 0.018113   |   0.0 | 19.46
Output  | 0.0013649  | 0.0013649  | 0.0013649  |   0.0 |  1.47
Other   |            | 0.0012     |            |       |  1.29

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

Particle-moves/CPUsec/proc: 3.2048e+07
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
