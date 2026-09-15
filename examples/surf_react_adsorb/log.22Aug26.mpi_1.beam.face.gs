SPARTA (24 Sep 2025)
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
  CPU time = 0.00108553 secs
  create/ghost percent = 97.5187 2.48128
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000116641 secs
  reassign/sort/migrate/ghost percent = 85.2522 0.53326 9.90732 4.30723

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

region              circle1 cylinder z  0 -10 1 INF INF

fix                 in1 emit/face/file air zhi data.beam beam_area_1 nevery 100 region circle1 twopass

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
  modify    (ave,min,max) = 0 0 0
  total     (ave,min,max) = 1.51379 1.51379 1.51379
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
      10    9.906e-06        0        0        0        0        0 
      20   3.8526e-05        0        0        0        0        0 
      30   6.7503e-05        0        0        0        0        0 
      40    7.657e-05        0        0        0        0        0 
      50   8.6626e-05        0        0        0        0        0 
      60  0.000119889        0        0        0        0        0 
      70  0.000132249        0        0        0        0        0 
      80  0.000168524        0        0        0        0        0 
      90  0.000201818        0        0        0        0        0 
     100  0.002038751     3149        0        0        0        0 
     110  0.002333005     3149        0        0        0        0 
     120  0.002623935     3149        0        0        0        0 
     130  0.002903586     3149        0        0        0        0 
     140  0.003169759     3149        0        0        0        0 
     150   0.00350742     3149        0        0        0        0 
     160  0.003790605     3149        0        0        0        0 
     170  0.004052949     3149        0        0        0        0 
     180  0.004317764     3149        0        0        0        0 
     190  0.004609908     3149        0        0        0        0 
     200  0.006249508     3230        0        0        0        0 
     210  0.006536667     3204        0        0        0        0 
     220  0.006823635     3204        0        0        0        0 
     230  0.007178355     3204        0        0        0        0 
     240  0.007496845     3204        0        0        0        0 
     250  0.007806934     3204        0        0        0        0 
     260  0.008079409     3204        0        0        0        0 
     270   0.00840513     3204        0        0        0        0 
     280  0.008715336     3204        0        0        0        0 
     290  0.008986759     3204        0        0        0        0 
     300  0.010594032     3301        0        0        0        0 
     310  0.010883363     3274        0        0        0        0 
     320  0.011181155     3274        0        0        0        0 
     330  0.011708801     3273        0        0        0        0 
     340  0.012036167     3272        0        0        0        0 
     350  0.012378939     3270        0        0        0        0 
     360  0.012682255     3268        0        0        0        0 
     370  0.012968119     3266        0        0        0        0 
     380  0.013242897     3266        0        0        0        0 
     390  0.013530864     3261        0        0        0        0 
     400  0.015158344     3379        0        0        0        0 
     410  0.015538132     3345        0        0        0        0 
     420  0.015847083     3342        0        0        0        0 
     430  0.016139419     3341        0        0        0        0 
     440  0.016428981     3336        0        0        0        0 
     450  0.016769512     3327        0        0        0        0 
     460  0.017056166     3321        0        0        0        0 
     470  0.017341897     3315        0        0        0        0 
     480  0.017640525     3309        0        0        0        0 
     490  0.017917437     3305        0        0        0        0 
     500  0.019554039     3385        0        0        0        0 
     510  0.019879931     3344        0        0        0        0 
     520  0.020162399     3335        0        0        0        0 
     530  0.020447837     3327        0        0        0        0 
     540  0.020768419     3319        0        0        0        0 
     550  0.021156716     3313        0        0        0        0 
     560  0.021447301     3302        0        0        0        0 
     570  0.021749149     3290        0        0        0        0 
     580  0.022036556     3283        0        0        0        0 
     590  0.022418605     3276        0        0        0        0 
     600  0.024177568     3356        0        0        0        0 
     610  0.024489141     3317        0        0        0        0 
     620  0.024829422     3313        0        0        0        0 
     630  0.025138416     3311        0        0        0        0 
     640  0.025430094     3303        0        0        0        0 
     650  0.025770695     3296        0        0        0        0 
     660  0.026066276     3290        0        0        0        0 
     670  0.026348038     3283        0        0        0        0 
     680  0.026645485     3277        0        0        0        0 
     690  0.026920619     3270        0        0        0        0 
     700  0.028634495     3375        0        0        0        0 
     710  0.028942977     3337        0        0        0        0 
     720  0.029225319     3335        0        0        0        0 
     730  0.029521893     3326        0        0        0        0 
     740  0.029830816     3319        0        0        0        0 
     750   0.03014121     3314        0        0        0        0 
     760  0.030425419     3310        0        0        0        0 
     770  0.030723146     3304        0        0        0        0 
     780  0.030999713     3295        0        0        0        0 
     790  0.031276533     3286        0        0        0        0 
     800  0.033056974     3410        0        0        0        0 
     810  0.033367827     3372        0        0        0        0 
     820  0.033675261     3359        0        0        0        0 
     830  0.033955735     3353        0        0        0        0 
     840  0.034235655     3330        0        0        0        0 
     850  0.034548338     3326        0        0        0        0 
     860  0.034892404     3323        0        0        0        0 
     870  0.035181232     3312        0        0        0        0 
     880  0.035512887     3301        0        0        0        0 
     890  0.035813044     3296        0        0        0        0 
     900  0.037409455     3388        0        0        0        0 
     910  0.037739991     3354        0        0        0        0 
     920  0.038030187     3344        0        0        0        0 
     930  0.038315552     3337        0        0        0        0 
     940  0.038620252     3329        0        0        0        0 
     950  0.038931641     3320        0        0        0        0 
     960  0.039278835     3314        0        0        0        0 
     970  0.039699567     3305        0        0        0        0 
     980   0.03999796     3296        0        0        0        0 
     990  0.040378393     3288        0        0        0        0 
    1000  0.042009308     3387        0        0        0        0 
Loop time of 0.0420447 on 1 procs for 1000 steps with 3387 particles
Performance: 23784.187 timesteps/s, 80.557 Mparticle-step/s

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.031121   | 0.031121   | 0.031121   |   0.0 | 74.02
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.00040146 | 0.00040146 | 0.00040146 |   0.0 |  0.95
Modify  | 0.007106   | 0.007106   | 0.007106   |   0.0 | 16.90
Output  | 0.0028999  | 0.0028999  | 0.0028999  |   0.0 |  6.90
MPI Sync| 0.00016559 | 0.00016559 | 0.00016559 |   0.0 |  0.39
Other   |            | 0.0003508  |            |       |  0.83

Particle moves    = 2983471 (2.98M)
Cells touched     = 3026752 (3.03M)
Particle comms    = 0 (0K)
Boundary collides = 646 (0.646K)
Boundary exits    = 434 (0.434K)
SurfColl checks   = 0 (0K)
SurfColl occurs   = 0 (0K)
Surf reactions    = 28166 (28.2K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 7.09594e+07
Particle-moves/step: 2983.47
Cell-touches/particle/step: 1.01451
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.000216526
Particle fraction exiting boundary: 0.000145468
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0.00944068
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 28166
    reaction O(g) --> O(s): 20991
    reaction O(g) + O(s) --> CO2(g): 2
    reaction O(g) --> CO(s): 6529
    reaction O(g) --> CO(g): 639
    reaction O(g) + O(s) --> O(g) + O(g): 5

Particles: 3387 ave 3387 max 3387 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      8 ave 8 max 8 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
