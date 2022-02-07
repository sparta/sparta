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
Created 8 child grid cells
  CPU time = 0.000919342 secs
  create/ghost percent = 87.3185 12.6815
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000153303 secs
  reassign/sort/migrate/ghost percent = 57.3872 1.39969 29.0824 12.1306

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
      10 2.7894974e-05        0        0        0        0        0 
      20 8.3684921e-05        0        0        0        0        0 
      30 0.00013661385        0        0        0        0        0 
      40 0.00018930435        0        0        0        0        0 
      50 0.00025105476        0        0        0        0        0 
      60 0.00030565262        0        0        0        0        0 
      70 0.00038599968        0        0        0        0        0 
      80 0.00043964386        0        0        0        0        0 
      90 0.00049233437        0        0        0        0        0 
     100 0.0048794746     3150        0        0        0        0 
     110 0.0065524578     3150        0        0        0        0 
     120 0.0082764626     3150        0        0        0        0 
     130 0.0099341869     3150        0        0        0        0 
     140  0.011614323     3150        0        0        0        0 
     150  0.013439655     3150        0        0        0        0 
     160  0.015128613     3150        0        0        0        0 
     170  0.016793013     3150        0        0        0        0 
     180  0.018508196     3150        0        0        0        0 
     190  0.020168066     3150        0        0        0        0 
     200  0.026780128     3230        0        0        0        0 
     210  0.028511763     3204        0        0        0        0 
     220  0.030200005     3204        0        0        0        0 
     230  0.031893015     3204        0        0        0        0 
     240  0.033742666     3204        0        0        0        0 
     250  0.035599947     3204        0        0        0        0 
     260   0.03729701     3204        0        0        0        0 
     270  0.038981199     3204        0        0        0        0 
     280  0.040676117     3204        0        0        0        0 
     290  0.042370796     3204        0        0        0        0 
     300  0.048898458     3299        0        0        0        0 
     310  0.050650835     3272        0        0        0        0 
     320  0.052418232     3272        0        0        0        0 
     330  0.054188251     3271        0        0        0        0 
     340  0.055963755     3270        0        0        0        0 
     350  0.057893991     3268        0        0        0        0 
     360  0.059642076     3267        0        0        0        0 
     370  0.061374664     3266        0        0        0        0 
     380  0.063093424     3264        0        0        0        0 
     390  0.064819574     3257        0        0        0        0 
     400  0.071420193     3379        0        0        0        0 
     410  0.073209047     3346        0        0        0        0 
     420  0.074982405     3345        0        0        0        0 
     430  0.076752663     3343        0        0        0        0 
     440   0.07852149     3339        0        0        0        0 
     450  0.080447435     3330        0        0        0        0 
     460  0.082199812     3322        0        0        0        0 
     470  0.083984852     3315        0        0        0        0 
     480  0.085782051     3312        0        0        0        0 
     490  0.087537527     3306        0        0        0        0 
     500  0.094262123     3384        0        0        0        0 
     510  0.096057892     3343        0        0        0        0 
     520  0.097849846     3337        0        0        0        0 
     530  0.099646807     3328        0        0        0        0 
     540    0.1014142     3318        0        0        0        0 
     550   0.10333514     3312        0        0        0        0 
     560   0.10507703     3302        0        0        0        0 
     570   0.10682034     3292        0        0        0        0 
     580   0.10856271     3285        0        0        0        0 
     590   0.11030316     3278        0        0        0        0 
     600   0.11686397     3356        0        0        0        0 
     610    0.1186583     3317        0        0        0        0 
     620   0.12041903     3311        0        0        0        0 
     630   0.12216473     3308        0        0        0        0 
     640   0.12391186     3303        0        0        0        0 
     650    0.1258564     3297        0        0        0        0 
     660   0.12760067     3288        0        0        0        0 
     670   0.12941933     3280        0        0        0        0 
     680   0.13114786     3274        0        0        0        0 
     690   0.13302112     3268        0        0        0        0 
     700   0.13963485     3374        0        0        0        0 
     710   0.14144897     3340        0        0        0        0 
     720   0.14325047     3335        0        0        0        0 
     730   0.14500999     3326        0        0        0        0 
     740   0.14684176     3320        0        0        0        0 
     750   0.14876175     3314        0        0        0        0 
     760   0.15051675     3310        0        0        0        0 
     770   0.15226603     3303        0        0        0        0 
     780   0.15400076     3295        0        0        0        0 
     790   0.15581799     3286        0        0        0        0 
     800   0.16245794     3406        0        0        0        0 
     810   0.16428876     3368        0        0        0        0 
     820   0.16606355     3356        0        0        0        0 
     830    0.1678412     3351        0        0        0        0 
     840   0.16961145     3329        0        0        0        0 
     850   0.17156696     3321        0        0        0        0 
     860   0.17332363     3319        0        0        0        0 
     870   0.17506552     3306        0        0        0        0 
     880    0.1768105     3299        0        0        0        0 
     890   0.17855597     3294        0        0        0        0 
     900   0.18525147     3385        0        0        0        0 
     910   0.18707705     3350        0        0        0        0 
     920   0.18888879     3341        0        0        0        0 
     930   0.19066572     3333        0        0        0        0 
     940   0.19242978     3326        0        0        0        0 
     950   0.19434714     3318        0        0        0        0 
     960   0.19613814     3310        0        0        0        0 
     970   0.19788337     3299        0        0        0        0 
     980   0.19965029     3293        0        0        0        0 
     990   0.20140982     3286        0        0        0        0 
    1000   0.20792127     3388        0        0        0        0 
Loop time of 0.207971 on 1 procs for 1000 steps with 3388 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.16714    | 0.16714    | 0.16714    |   0.0 | 80.37
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0015428  | 0.0015428  | 0.0015428  |   0.0 |  0.74
Modify  | 0.032544   | 0.032544   | 0.032544   |   0.0 | 15.65
Output  | 0.0050285  | 0.0050285  | 0.0050285  |   0.0 |  2.42
Other   |            | 0.001711   |            |       |  0.82

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

Particle-moves/CPUsec/proc: 1.43424e+07
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
