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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (../grid.cpp:410)
Created 8 child grid cells
  CPU time = 0.00170452 secs
  create/ghost percent = 94.0012 5.99875
balance_grid        rcb cell
Balance grid migrated 4 cells
  CPU time = 0.000946374 secs
  reassign/sort/migrate/ghost percent = 70.3764 0.265646 14.3395 15.0185

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
      10  0.000298719        0        0        0        0        0 
      20  0.000545474        0        0        0        0        0 
      30  0.000798376        0        0        0        0        0 
      40  0.001041081        0        0        0        0        0 
      50  0.001280503        0        0        0        0        0 
      60  0.001519018        0        0        0        0        0 
      70  0.001768148        0        0        0        0        0 
      80   0.00201204        0        0        0        0        0 
      90   0.00226571        0        0        0        0        0 
     100  0.005514184     3112        0        0        0        0 
     110  0.006164214     3112        0        0        0        0 
     120  0.006814523     3112        0        0        0        0 
     130  0.007439689     3112        0        0        0        0 
     140  0.008055007     3112        0        0        0        0 
     150  0.008700357     3112        0        0        0        0 
     160  0.009312811     3112        0        0        0        0 
     170  0.009929596     3112        0        0        0        0 
     180  0.010539466     3112        0        0        0        0 
     190  0.012737569     3112        0        0        0        0 
     200  0.015136959     3273        0        0        0        0 
     210  0.015968092     3242        0        0        0        0 
     220  0.016739859     3242        0        0        0        0 
     230  0.017421667     3242        0        0        0        0 
     240  0.018136442     3242        0        0        0        0 
     250  0.018784865     3242        0        0        0        0 
     260  0.019362747     3242        0        0        0        0 
     270  0.019960465     3242        0        0        0        0 
     280  0.020558392     3242        0        0        0        0 
     290  0.021168052     3242        0        0        0        0 
     300  0.023286813     3357        0        0        0        0 
     310   0.02400955     3329        0        0        0        0 
     320  0.024668938     3329        0        0        0        0 
     330  0.025338175     3326        0        0        0        0 
     340  0.026005526     3324        0        0        0        0 
     350  0.026641727     3321        0        0        0        0 
     360  0.027243286     3318        0        0        0        0 
     370  0.027834019     3315        0        0        0        0 
     380  0.028421819     3311        0        0        0        0 
     390  0.029017092     3307        0        0        0        0 
     400  0.031044288     3415        0        0        0        0 
     410  0.031710382     3381        0        0        0        0 
     420  0.032409302     3378        0        0        0        0 
     430  0.033045852     3372        0        0        0        0 
     440  0.033673393     3367        0        0        0        0 
     450  0.034340115     3364        0        0        0        0 
     460  0.034932455     3357        0        0        0        0 
     470  0.035526889     3352        0        0        0        0 
     480  0.036134664     3345        0        0        0        0 
     490  0.036697879     3336        0        0        0        0 
     500  0.038580081     3392        0        0        0        0 
     510  0.039278233     3351        0        0        0        0 
     520  0.039934549     3347        0        0        0        0 
     530  0.040550006     3341        0        0        0        0 
     540  0.041141019     3331        0        0        0        0 
     550  0.041709612     3330        0        0        0        0 
     560  0.042346651     3320        0        0        0        0 
     570  0.042903859     3316        0        0        0        0 
     580  0.043432502     3305        0        0        0        0 
     590  0.044018346     3300        0        0        0        0 
     600  0.045932746     3406        0        0        0        0 
     610   0.04663425     3361        0        0        0        0 
     620  0.047273804     3356        0        0        0        0 
     630  0.047923903     3349        0        0        0        0 
     640  0.048483975     3339        0        0        0        0 
     650  0.049090633     3331        0        0        0        0 
     660  0.049631498     3320        0        0        0        0 
     670  0.050180534     3312        0        0        0        0 
     680  0.050701843     3304        0        0        0        0 
     690  0.051259052     3299        0        0        0        0 
     700  0.053095367     3402        0        0        0        0 
     710  0.053775918     3368        0        0        0        0 
     720  0.054395567     3360        0        0        0        0 
     730  0.055048251     3356        0        0        0        0 
     740  0.055634723     3348        0        0        0        0 
     750  0.056252555     3343        0        0        0        0 
     760  0.056780011     3334        0        0        0        0 
     770  0.057325066     3332        0        0        0        0 
     780  0.057829404     3323        0        0        0        0 
     790  0.058376275     3315        0        0        0        0 
     800   0.06013269     3421        0        0        0        0 
     810  0.060794104     3397        0        0        0        0 
     820  0.061431562     3389        0        0        0        0 
     830  0.062078519     3385        0        0        0        0 
     840  0.062673024     3379        0        0        0        0 
     850  0.063266621     3367        0        0        0        0 
     860  0.063868808     3361        0        0        0        0 
     870  0.064385088     3353        0        0        0        0 
     880  0.064893686     3346        0        0        0        0 
     890  0.065444608     3342        0        0        0        0 
     900  0.067225119     3388        0        0        0        0 
     910   0.06787466     3351        0        0        0        0 
     920  0.068430192     3343        0        0        0        0 
     930  0.069047256     3338        0        0        0        0 
     940  0.069560603     3331        0        0        0        0 
     950  0.070088407     3325        0        0        0        0 
     960   0.07055028     3318        0        0        0        0 
     970  0.071096733     3315        0        0        0        0 
     980  0.071588498     3309        0        0        0        0 
     990  0.072116442     3299        0        0        0        0 
    1000  0.073785483     3381        0        0        0        0 
Loop time of 0.0738098 on 4 procs for 1000 steps with 3381 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.0084472  | 0.020039   | 0.03455    |   8.3 | 27.15
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.013568   | 0.014815   | 0.015852   |   0.7 | 20.07
Modify  | 0.0018823  | 0.0068144  | 0.012506   |   6.0 |  9.23
Output  | 0.0016836  | 0.0019717  | 0.0027673  |   1.0 |  2.67
Other   |            | 0.03017    |            |       | 40.88

Particle moves    = 3004552 (3M)
Cells touched     = 3048162 (3.05M)
Particle comms    = 14743 (14.7K)
Boundary collides = 644 (0.644K)
Boundary exits    = 431 (0.431K)
SurfColl checks   = 0 (0K)
SurfColl occurs   = 0 (0K)
Surf reactions    = 28387 (28.4K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0
Axisymm bad moves = 0

Particle-moves/CPUsec/proc: 1.01767e+07
Particle-moves/step: 3004.55
Cell-touches/particle/step: 1.01451
Particle comm iterations/step: 1.419
Particle fraction communicated: 0.00490689
Particle fraction colliding with boundary: 0.000214341
Particle fraction exiting boundary: 0.000143449
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0.009448
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Surface reaction tallies:
  id adsorb_test_gs2 style adsorb #-of-reactions 9
    reaction all: 28387
    reaction O(g) --> O(s): 21263
    reaction O(g) + O(s) --> CO2(g): 1
    reaction O(g) --> CO(s): 6480
    reaction O(g) --> CO(g): 637
    reaction O(g) + O(s) --> O(g) + O(g): 6

Particles: 845.25 ave 1642 max 61 min
Histogram: 2 0 0 0 0 0 0 0 0 2
Cells:      2 ave 2 max 2 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 6 ave 6 max 6 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 6 ave 6 max 6 min
Histogram: 4 0 0 0 0 0 0 0 0 0
