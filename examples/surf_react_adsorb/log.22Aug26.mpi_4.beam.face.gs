SPARTA (24 Sep 2025)
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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/home/user/sparta/src/grid.cpp:486)
Created 8 child grid cells
  CPU time = 0.00145738 secs
  create/ghost percent = 90.0735 9.92648
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000362139 secs
  reassign/sort/migrate/ghost percent = 66.2276 2.19142 10.3466 21.2344

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
      10  0.000176694        0        0        0        0        0 
      20  0.000327012        0        0        0        0        0 
      30  0.000498606        0        0        0        0        0 
      40    0.0006214        0        0        0        0        0 
      50  0.000778711        0        0        0        0        0 
      60  0.000896651        0        0        0        0        0 
      70  0.001015864        0        0        0        0        0 
      80  0.001130213        0        0        0        0        0 
      90   0.00127682        0        0        0        0        0 
     100  0.003177437     3141        0        0        0        0 
     110  0.003584792     3141        0        0        0        0 
     120  0.004002611     3141        0        0        0        0 
     130  0.004297821     3141        0        0        0        0 
     140  0.004595496     3141        0        0        0        0 
     150  0.004862795     3141        0        0        0        0 
     160  0.005099541     3141        0        0        0        0 
     170  0.005382073     3141        0        0        0        0 
     180  0.005710775     3141        0        0        0        0 
     190  0.007271319     3141        0        0        0        0 
     200  0.008598137     3213        0        0        0        0 
     210  0.009115735     3187        0        0        0        0 
     220   0.00957046     3187        0        0        0        0 
     230  0.009949352     3187        0        0        0        0 
     240  0.010311606     3187        0        0        0        0 
     250  0.010698143     3187        0        0        0        0 
     260  0.011037813     3187        0        0        0        0 
     270  0.011314185     3187        0        0        0        0 
     280  0.011651977     3187        0        0        0        0 
     290   0.01194862     3187        0        0        0        0 
     300  0.013111442     3304        0        0        0        0 
     310  0.013497513     3274        0        0        0        0 
     320  0.013819228     3272        0        0        0        0 
     330   0.01411953     3272        0        0        0        0 
     340  0.014496838     3271        0        0        0        0 
     350  0.014842899     3270        0        0        0        0 
     360  0.015149289     3268        0        0        0        0 
     370   0.01540583     3267        0        0        0        0 
     380  0.015766973     3265        0        0        0        0 
     390  0.016042216     3263        0        0        0        0 
     400  0.017178609     3354        0        0        0        0 
     410  0.017594432     3325        0        0        0        0 
     420  0.017899836     3322        0        0        0        0 
     430  0.018217843     3318        0        0        0        0 
     440  0.018719146     3314        0        0        0        0 
     450  0.019245263     3311        0        0        0        0 
     460  0.019777146     3302        0        0        0        0 
     470  0.020173204     3294        0        0        0        0 
     480  0.020422914     3290        0        0        0        0 
     490   0.02076718     3283        0        0        0        0 
     500  0.021808763     3340        0        0        0        0 
     510   0.02212597     3306        0        0        0        0 
     520  0.022573524     3299        0        0        0        0 
     530  0.022913572     3292        0        0        0        0 
     540  0.023256225     3288        0        0        0        0 
     550  0.023640712     3282        0        0        0        0 
     560   0.02400112     3276        0        0        0        0 
     570  0.024359219     3270        0        0        0        0 
     580  0.024739879     3270        0        0        0        0 
     590  0.025046327     3261        0        0        0        0 
     600  0.026078083     3403        0        0        0        0 
     610  0.026525886     3366        0        0        0        0 
     620  0.026857782     3361        0        0        0        0 
     630  0.027173369     3356        0        0        0        0 
     640  0.027504164     3348        0        0        0        0 
     650  0.027817657     3342        0        0        0        0 
     660  0.028094057     3333        0        0        0        0 
     670  0.028360598     3324        0        0        0        0 
     680   0.02873432     3319        0        0        0        0 
     690  0.029012862     3311        0        0        0        0 
     700  0.030041332     3403        0        0        0        0 
     710  0.030505113     3369        0        0        0        0 
     720  0.030855768     3365        0        0        0        0 
     730  0.031172213     3357        0        0        0        0 
     740    0.0314825     3351        0        0        0        0 
     750  0.031842863     3346        0        0        0        0 
     760  0.032101801     3340        0        0        0        0 
     770  0.032364309     3328        0        0        0        0 
     780  0.032692538     3318        0        0        0        0 
     790  0.032993145     3313        0        0        0        0 
     800  0.034327281     3446        0        0        0        0 
     810  0.034754402     3393        0        0        0        0 
     820  0.035089584     3387        0        0        0        0 
     830  0.035422395     3382        0        0        0        0 
     840  0.035802768     3376        0        0        0        0 
     850  0.036081358     3369        0        0        0        0 
     860  0.036361994     3360        0        0        0        0 
     870  0.036680547     3354        0        0        0        0 
     880  0.037051293     3347        0        0        0        0 
     890   0.03738311     3341        0        0        0        0 
     900  0.038529525     3359        0        0        0        0 
     910   0.03890522     3316        0        0        0        0 
     920  0.039186576     3309        0        0        0        0 
     930  0.039509704     3299        0        0        0        0 
     940   0.03982448     3291        0        0        0        0 
     950  0.040138196     3281        0        0        0        0 
     960  0.040409795     3272        0        0        0        0 
     970  0.040746389     3263        0        0        0        0 
     980  0.041015896     3257        0        0        0        0 
     990  0.041275833     3251        0        0        0        0 
    1000  0.042277372     3393        0        0        0        0 
Loop time of 0.0422895 on 4 procs for 1000 steps with 3393 particles
Performance: 23646.559 timesteps/s, 80.233 Mparticle-step/s

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0026514  | 0.009688   | 0.016873   |   7.1 | 22.91
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.010018   | 0.010323   | 0.010666   |   0.2 | 24.41
Modify  | 7.3928e-05 | 0.0025649  | 0.0050964  |   4.9 |  6.07
Output  | 0.0011582  | 0.0014406  | 0.0022417  |   1.2 |  3.41
MPI Sync| 0.0041678  | 0.01458    | 0.024291   |   8.0 | 34.48
Other   |            | 0.003693   |            |       |  8.73

Particle moves    = 2956450 (2.96M)
Cells touched     = 3026859 (3.03M)
Particle comms    = 14716 (14.7K)
Boundary collides = 629 (0.629K)
Boundary exits    = 402 (0.402K)
SurfColl checks   = 0 (0K)
SurfColl occurs   = 0 (0K)
Surf reactions    = 28229 (28.2K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 1.74775e+07
Particle-moves/step: 2956.45
Cell-touches/particle/step: 1.02382
Particle comm iterations/step: 1.431
Particle fraction communicated: 0.00497759
Particle fraction colliding with boundary: 0.000212755
Particle fraction exiting boundary: 0.000135974
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0.00954828
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 28229
    reaction O(g) --> O(s): 20941
    reaction O(g) + O(s) --> CO2(g): 4
    reaction O(g) --> CO(s): 6659
    reaction O(g) --> CO(g): 617
    reaction O(g) + O(s) --> O(g) + O(g): 8

Particles: 848.25 ave 1643 max 59 min
Histogram: 2 0 0 0 0 0 0 0 0 2
Cells:      2 ave 2 max 2 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 6 ave 6 max 6 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 6 ave 6 max 6 min
Histogram: 4 0 0 0 0 0 0 0 0 0
