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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (/Users/eharvey/dev/SPARTA.base/sparta/src/grid.cpp:415)
Created 200 child grid cells
  parent cells = 1
  CPU time = 0.000969 secs
  create/ghost percent = 96.0784 3.92157
balance_grid        rcb cell
Balance grid migrated 120 cells
  CPU time = 0.000684 secs
  reassign/sort/migrate/ghost percent = 82.0175 0.584795 6.72515 10.6725

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
  CPU time = 0.000476 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 41.5966 9.45378 0.420168 30.4622 18.0672 9.03361 0.210084
  surf2grid time = 0.000145 secs
  map/rvous1/rvous2/split percent = 11.0345 63.4483 0.689655 15.1724

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
     100     0.049865    18354     1533      776       94     4808 
     200      0.12481    27288     2808     1503       89     5719 
     300     0.202188    31059     3445     1867      107     6258 
     400     0.279464    32619     3776     1983       94     6480 
     500     0.364844    33668     3888     2045      101     6508 
     600      0.44358    33716     3913     1944       95     6641 
     700     0.526573    33751     3956     2015      105     6424 
     800     0.615888    33785     3951     2042      118     6580 
     900     0.708015    33676     4010     2025      116     6538 
    1000     0.795262    33672     3998     2056      113     6407 
Loop time of 0.795277 on 4 procs for 1000 steps with 33672 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.22369    | 0.30876    | 0.46738    |  17.9 | 38.82
Coll    | 0.090394   | 0.13881    | 0.19846    |  10.8 | 17.45
Sort    | 0.010892   | 0.014426   | 0.020138   |   2.9 |  1.81
Comm    | 0.040734   | 0.045951   | 0.056481   |   2.9 |  5.78
Modify  | 0.000148   | 0.010617   | 0.041957   |  17.6 |  1.33
Output  | 0.000137   | 0.00078375 | 0.001263   |   0.0 |  0.10
Other   |            | 0.2759     |            |       | 34.70

Particle moves    = 29910872 (29.9M)
Cells touched     = 33685601 (33.7M)
Particle comms    = 403445 (0.403M)
Boundary collides = 80871 (80.9K)
Boundary exits    = 133392 (0.133M)
SurfColl checks   = 5965475 (5.97M)
SurfColl occurs   = 101044 (0.101M)
Surf reactions    = 0 (0K)
Collide attempts  = 3336412 (3.34M)
Collide occurs    = 1752308 (1.75M)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 9.40266e+06
Particle-moves/step: 29910.9
Cell-touches/particle/step: 1.1262
Particle comm iterations/step: 2.224
Particle fraction communicated: 0.0134882
Particle fraction colliding with boundary: 0.00270373
Particle fraction exiting boundary: 0.00445965
Surface-checks/particle/step: 0.199442
Surface-collisions/particle/step: 0.00337817
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0.111545
Collisions/particle/step: 0.0585843
Reactions/particle/step: 0

Particles: 8418 ave 11132 max 6501 min
Histogram: 1 0 1 1 0 0 0 0 0 1
Cells:      50 ave 50 max 50 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostCell: 15 ave 20 max 10 min
Histogram: 2 0 0 0 0 0 0 0 0 2
EmptyCell: 15 ave 20 max 10 min
Histogram: 2 0 0 0 0 0 0 0 0 2
Surfs:    25 ave 25 max 25 min
Histogram: 4 0 0 0 0 0 0 0 0 0
GhostSurf: 0 ave 0 max 0 min
Histogram: 4 0 0 0 0 0 0 0 0 0
