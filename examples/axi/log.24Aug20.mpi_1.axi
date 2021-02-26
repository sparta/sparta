SPARTA (6 Jul 2020)
################################################################################
# 2d axisymmetric flow around a circle with specular reflections
#
# Note:
#  - The "comm/sort” option to the “global” command is used to match MPI runs.
#  - The “twopass” option is used to match Kokkos runs.
# The "comm/sort" and "twopass" options should not be used for production runs.
################################################################################

seed	    	    12345
dimension   	    2
global              gridcut 0.0 comm/sort yes

boundary	    o ar p

create_box          -0.25 0.25 0.0 0.25 -0.5 0.5
Created orthogonal box = (-0.25 0 -0.5) to (0.25 0.25 0.5)

create_grid 	    20 10 1
Created 200 child grid cells
  parent cells = 1
  CPU time = 0.001502 secs
  create/ghost percent = 91.012 8.98802
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000357 secs
  reassign/sort/migrate/ghost percent = 77.8711 0.840336 8.12325 13.1653

global		    nrho 1.e20 fnum 1.e17 weight cell radius

species		    air.species N2
mixture		    air N2 vstream 3472.0 0.0 0.0 temp 300.0

fix                 in emit/face air xlo twopass
collide		    vss air air.vss

read_surf           data.circle origin 5 5 0                     trans -5 -5 0 scale 0.05 0.05 1 clip
  50 points
  50 lines
  clipped to 25 lines
  -0.15 0.15 xlo xhi
  0 0.149704 ylo yhi
  0 0 zlo zhi
  0.0188372 min line length
  24 = cells with surfs
  48 = total surfs in all grid cells
  3 = max surfs in one grid cell
  0.753486 = min surf-size/cell-size ratio
  0 0 = number of pushed cells
  24 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  132 44 24 = cells outside/inside/overlapping surfs
  24 = surf cells with 1,2,etc splits
  0.0840933 0.0840933 = cell-wise and global flow volume
  CPU time = 0.000528 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 37.8788 10.0379 1.89394 41.0985 9.09091 8.90152 0
  surf2grid time = 0.000217 secs
  map/rvous1/rvous2/split percent = 19.8157 33.1797 0 40.553

surf_collide	    1 specular
surf_modify         all collide 1

timestep 	    1e-6

#dump                2 image all 50 image.*.ppm type type pdiam 0.002 #	            size 512 512 particle yes #                    gline yes 0.01 #                    surf proc 0.02 zoom 3.5
#dump_modify	    2 pad 4

stats		    100
stats_style	    step cpu np nattempt ncoll nscoll nscheck

run 		    1000
Memory usage per proc in Mbytes:
  particles (ave,min,max) = 0 0 0
  grid      (ave,min,max) = 1.51388 1.51388 1.51388
  surf      (ave,min,max) = 0.00257492 0.00257492 0.00257492
  total     (ave,min,max) = 1.51645 1.51645 1.51645
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
     100     0.094122    18343     1517      797       82     4885 
     200     0.309926    27124     2797     1543      125     5921 
     300      0.51919    30869     3454     1837      116     6177 
     400     0.762453    32366     3677     1906      102     6249 
     500     1.021153    33246     3844     1991      106     6600 
     600     1.278854    33697     3945     2025      101     6365 
     700     1.524732    33891     3977     2114      116     6593 
     800      1.73976    33811     3950     2046       98     6435 
     900     1.986161    33858     4058     2122      102     6896 
    1000     2.225281    33600     3955     2035      117     6644 
Loop time of 2.22529 on 1 procs for 1000 steps with 33600 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 1.3446     | 1.3446     | 1.3446     |   0.0 | 60.42
Coll    | 0.63109    | 0.63109    | 0.63109    |   0.0 | 28.36
Sort    | 0.072937   | 0.072937   | 0.072937   |   0.0 |  3.28
Comm    | 0.12673    | 0.12673    | 0.12673    |   0.0 |  5.69
Modify  | 0.049107   | 0.049107   | 0.049107   |   0.0 |  2.21
Output  | 0.000146   | 0.000146   | 0.000146   |   0.0 |  0.01
Other   |            | 0.000693   |            |       |  0.03

Particle moves    = 29814692 (29.8M)
Cells touched     = 33581089 (33.6M)
Particle comms    = 0 (0K)
Boundary collides = 81123 (81.1K)
Boundary exits    = 131872 (0.132M)
SurfColl checks   = 5951496 (5.95M)
SurfColl occurs   = 100932 (0.101M)
Surf reactions    = 0 (0K)
Collide attempts  = 3326303 (3.33M)
Collide occurs    = 1741265 (1.74M)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 1.33981e+07
Particle-moves/step: 29814.7
Cell-touches/particle/step: 1.12633
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.00272091
Particle fraction exiting boundary: 0.00442305
Surface-checks/particle/step: 0.199616
Surface-collisions/particle/step: 0.00338531
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0.111566
Collisions/particle/step: 0.0584029
Reactions/particle/step: 0

Particles: 33600 ave 33600 max 33600 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Cells:      200 ave 200 max 200 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
EmptyCell: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Surfs:    25 ave 25 max 25 min
Histogram: 1 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
