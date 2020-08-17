SPARTA (6 Jul 2020)
################################################################################
# beam of particles striking the surface at an inclined angle - free molecular flow (no collisions)
#
# Note:
#  - The "comm/sort” option to the “global” command is used to match MPI runs.
#  - The “twopass” option is used to match Kokkos runs.
# The "comm/sort" and "twopass" options should not be used for production runs.
################################################################################

seed	    	    123456
dimension   	    3
global              gridcut 0.0 comm/sort yes

boundary	    oo oo so

create_box          -11 11 -11 11 0 10
Created orthogonal box = (-11 -11 0) to (11 11 10)
create_grid 	    1 1 1
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/Users/eharvey/dev/SPARTA.base/sparta/src/grid.cpp:415)
Created 1 child grid cells
  parent cells = 1
  CPU time = 0.000901 secs
  create/ghost percent = 95.5605 4.43951
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000616 secs
  reassign/sort/migrate/ghost percent = 89.6104 0.487013 3.57143 6.33117

global		    nrho 1e10 fnum 1e6

species		    air.species N O
mixture		    air N O vstream 0 1000 -1000

mixture             air N frac 0.8
mixture             air O frac 0.2

surf_collide	    1 td 1000.0 #barrier 1000
#surf_collide	    1 td 500.0 bond 500 0 0
#surf_collide	    1 td 300.0 initenergy 0.01 0 0
#surf_collide	    1 td 1000.0 barrier 1000 bond 500 0 0 initenergy 0.01 0 0

bound_modify        zlo collide 1

region              circle cylinder z 0 -10 1 -INF INF
fix                 in emit/face/file air zhi data.beam beam_area nevery 100 region circle

#dump                2 image all 10 image.*.ppm type type pdiam 0.2 surf proc 0.01 size 512 512 zoom 1.75 gline no 0.005
#dump_modify	    2 pad 4

timestep            0.0001

stats		    100
stats_style	    step cpu np nattempt ncoll nscoll nscheck
run 		    1000
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 0.435753 0.0763779 1.51388
  surf      (ave,min,max) = 0 0 0
  total     (ave,min,max) = 0.435753 0.0763779 1.51388
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
     100     0.001362     3138        0        0        0        0 
     200     0.004644     6292        0        0        0        0 
     300     0.009998     7629        0        0        0        0 
     400     0.015496     7828        0        0        0        0 
     500     0.022354     7864        0        0        0        0 
     600     0.029714     7884        0        0        0        0 
     700     0.035314     7810        0        0        0        0 
     800     0.040831     7826        0        0        0        0 
     900     0.047015     7763        0        0        0        0 
    1000     0.053292     7773        0        0        0        0 
Loop time of 0.053302 on 4 procs for 1000 steps with 7773 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 6.9e-05    | 0.010938   | 0.043489   |  18.0 | 20.52
Coll    | 0          | 0          | 0          |   0.0 |  0.00
Sort    | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.001159   | 0.0012812  | 0.0015     |   0.4 |  2.40
Modify  | 3.6e-05    | 0.00181    | 0.007104   |   7.2 |  3.40
Output  | 4.6e-05    | 5.625e-05  | 8.2e-05    |   0.0 |  0.11
Other   |            | 0.03922    |            |       | 73.57

Particle moves    = 5453001 (5.45M)
Cells touched     = 5453001 (5.45M)
Particle comms    = 0 (0K)
Boundary collides = 28255 (28.3K)
Boundary exits    = 23659 (23.7K)
SurfColl checks   = 0 (0K)
SurfColl occurs   = 0 (0K)
Surf reactions    = 0 (0K)
Collide attempts  = 0 (0K)
Collide occurs    = 0 (0K)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 2.5576e+07
Particle-moves/step: 5453
Cell-touches/particle/step: 1
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00518155
Particle fraction exiting boundary: 0.00433871
Surface-checks/particle/step: 0
Surface-collisions/particle/step: 0
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0
Collisions/particle/step: 0
Reactions/particle/step: 0

Particles: 1943.25 ave 7773 max 0 min
Histogram: 3 0 0 0 0 0 0 0 0 1
Cells:      0.25 ave 1 max 0 min
Histogram: 3 0 0 0 0 0 0 0 0 1
GhostCell: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0

