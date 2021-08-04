SPARTA (26 Feb 2021)
# beam of particles striking the surface at an inclined angle
# free molecular flow (no collisions)

seed	    	    123456
dimension   	    3
global              gridcut 0.0 comm/sort yes

boundary	    	oo oo so


create_box          -11 11 -11 11 0 10
Created orthogonal box = (-11 -11 0) to (11 11 10)
create_grid 	    2 2 2
Created 8 child grid cells
  CPU time = 0.00171709 secs
  create/ghost percent = 95.1597 4.84033
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.00032065 secs
  reassign/sort/migrate/ghost percent = 79.8952 0.239825 14.746 5.11898

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
      10   5.6294e-05        0        0        0        0        0 
      20  0.000115451        0        0        0        0        0 
      30  0.000169999        0        0        0        0        0 
      40   0.00022294        0        0        0        0        0 
      50  0.000275601        0        0        0        0        0 
      60  0.000328193        0        0        0        0        0 
      70  0.000380576        0        0        0        0        0 
      80  0.000432958        0        0        0        0        0 
      90   0.00048548        0        0        0        0        0 
     100  0.004863875     3133        0        0        0        0 
     110  0.005811157     3133        0        0        0        0 
     120   0.00674901     3133        0        0        0        0 
     130  0.007692591     3133        0        0        0        0 
     140  0.008625904     3133        0        0        0        0 
     150  0.009655462     3133        0        0        0        0 
     160  0.010592058     3133        0        0        0        0 
     170  0.011515873     3133        0        0        0        0 
     180  0.012388423     3133        0        0        0        0 
     190  0.013256852     3133        0        0        0        0 
     200  0.017363278     3255        0        0        0        0 
     210  0.018266559     3228        0        0        0        0 
     220  0.019156151     3228        0        0        0        0 
     230  0.020047768     3228        0        0        0        0 
     240  0.020942039     3228        0        0        0        0 
     250   0.02185901     3228        0        0        0        0 
     260   0.02272681     3228        0        0        0        0 
     270  0.023570725     3228        0        0        0        0 
     280  0.024408702     3228        0        0        0        0 
     290  0.025247728     3228        0        0        0        0 
     300  0.029083581     3275        0        0        0        0 
     310  0.029939718     3240        0        0        0        0 
     320  0.030758978     3240        0        0        0        0 
     330  0.031596537     3240        0        0        0        0 
     340  0.032432628     3238        0        0        0        0 
     350  0.033283946     3235        0        0        0        0 
     360   0.03408826     3231        0        0        0        0 
     370  0.034872318     3228        0        0        0        0 
     380  0.035661615     3226        0        0        0        0 
     390  0.036448328     3222        0        0        0        0 
     400  0.040049858     3360        0        0        0        0 
     410  0.040857105     3327        0        0        0        0 
     420  0.041667984     3325        0        0        0        0 
     430  0.042487104     3319        0        0        0        0 
     440  0.043248533     3316        0        0        0        0 
     450  0.044067164     3312        0        0        0        0 
     460  0.044819095     3308        0        0        0        0 
     470  0.045577592     3303        0        0        0        0 
     480  0.046336996     3296        0        0        0        0 
     490  0.047094025     3290        0        0        0        0 
     500  0.050449009     3384        0        0        0        0 
     510  0.051231741     3364        0        0        0        0 
     520  0.052007698     3355        0        0        0        0 
     530  0.052756626     3347        0        0        0        0 
     540   0.05349829     3341        0        0        0        0 
     550  0.054275364     3333        0        0        0        0 
     560  0.054996215     3326        0        0        0        0 
     570   0.05570617     3314        0        0        0        0 
     580  0.056433586     3307        0        0        0        0 
     590  0.057167917     3299        0        0        0        0 
     600  0.060350528     3391        0        0        0        0 
     610  0.061114472     3355        0        0        0        0 
     620   0.06184084     3351        0        0        0        0 
     630  0.062587393     3344        0        0        0        0 
     640  0.063307476     3334        0        0        0        0 
     650  0.064064226     3327        0        0        0        0 
     660  0.064749457     3319        0        0        0        0 
     670  0.065445793     3311        0        0        0        0 
     680  0.066143176     3306        0        0        0        0 
     690  0.066824845     3295        0        0        0        0 
     700  0.069843185     3435        0        0        0        0 
     710  0.070556353     3393        0        0        0        0 
     720  0.071262467     3389        0        0        0        0 
     730  0.071963901     3382        0        0        0        0 
     740  0.072660586     3375        0        0        0        0 
     750  0.073386885     3366        0        0        0        0 
     760  0.074059404     3355        0        0        0        0 
     770  0.074718514     3345        0        0        0        0 
     780  0.075387261     3336        0        0        0        0 
     790  0.076056359     3330        0        0        0        0 
     800  0.078946048     3433        0        0        0        0 
     810  0.079624294     3397        0        0        0        0 
     820  0.080304775     3384        0        0        0        0 
     830  0.080979041     3377        0        0        0        0 
     840  0.081644297     3368        0        0        0        0 
     850   0.08236843     3361        0        0        0        0 
     860  0.083050378     3350        0        0        0        0 
     870  0.083678128     3342        0        0        0        0 
     880  0.084311885     3334        0        0        0        0 
     890  0.084948155     3333        0        0        0        0 
     900  0.087660932     3359        0        0        0        0 
     910  0.088303139     3322        0        0        0        0 
     920  0.088925511     3315        0        0        0        0 
     930  0.089537057     3308        0        0        0        0 
     940  0.090159639     3299        0        0        0        0 
     950  0.090810996     3290        0        0        0        0 
     960  0.091427082     3282        0        0        0        0 
     970  0.092041911     3272        0        0        0        0 
     980  0.092659883     3261        0        0        0        0 
     990  0.093278064     3258        0        0        0        0 
    1000  0.095868475     3429        0        0        0        0 
Loop time of 0.0958799 on 1 procs for 1000 steps with 3429 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.07112    | 0.07112    | 0.07112    |   0.0 | 74.18
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0015554  | 0.0015554  | 0.0015554  |   0.0 |  1.62
Modify  | 0.020047   | 0.020047   | 0.020047   |   0.0 | 20.91
Output  | 0.0014248  | 0.0014248  | 0.0014248  |   0.0 |  1.49
Other   |            | 0.001733   |            |       |  1.81

Particle moves    = 2989180 (2.99M)
Cells touched     = 3032441 (3.03M)
Particle comms    = 0 (0K)
Boundary collides = 685 (0.685K)
Boundary exits    = 448 (0.448K)
SurfColl checks   = 0 (0K)
SurfColl occurs   = 0 (0K)
Surf reactions    = 28166 (28.2K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 3.11763e+07
Particle-moves/step: 2989.18
Cell-touches/particle/step: 1.01447
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00022916
Particle fraction exiting boundary: 0.000149874
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0.00942265
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 28166
    reaction O(g) --> O(s): 21181
    reaction O(g) + O(s) --> CO2(g): 1
    reaction O(g) --> CO(s): 6300
    reaction O(g) --> CO(g): 676
    reaction O(g) + O(s) --> O(g) + O(g): 8

Particles: 3429 ave 3429 max 3429 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      8 ave 8 max 8 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
