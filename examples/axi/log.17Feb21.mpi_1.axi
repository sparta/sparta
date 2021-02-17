SPARTA (20 Nov 2020)
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
  CPU time = 0.00100708 secs
  create/ghost percent = 86.6951 13.3049
balance_grid        rcb cell
Balance grid migrated 0 cells
  CPU time = 0.000123978 secs
  reassign/sort/migrate/ghost percent = 44.2308 0 26.7308 29.0385

global		    nrho 1.e20 fnum 2.5e15 weight cell radius

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
  0 0 = number of pushed cells
  24 0 = cells overlapping surfs, overlap cells with unmarked corner pts
  132 44 24 = cells outside/inside/overlapping surfs
  24 = surf cells with 1,2,etc splits
  0.0840933 0.0840933 = cell-wise and global flow volume
  CPU time = 0.000784159 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 38.644 12.6178 1.27698 42.0797 5.38157 9.0301 0.152022
  surf2grid time = 0.000329971 secs
  map/comm1/comm2/comm3/comm4/split percent = 38.5116 8.52601 9.9711 3.61272 11.4884 22.9769

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
  grid      (ave,min,max) = 1.51379 1.51379 1.51379
  surf      (ave,min,max) = 0.00257492 0.00257492 0.00257492
  total     (ave,min,max) = 1.51637 1.51637 1.51637
Step CPU Np Natt Ncoll Nscoll Nscheck 
       0            0        0        0        0        0        0 
     100   0.17789698    18448     1543      806       91     4772 
     200   0.60940099    27357     2861     1585      105     6005 
     300    1.1561711    31252     3518     1872      110     6275 
     400     1.764668    32885     3773     2007      117     6324 
     500     2.388835    33571     3893     2042       93     6340 
     600    3.0254631    33713     3944     2051      108     6647 
     700    3.6611409    33808     3974     2076       78     6552 
     800    4.3015249    33947     3997     2035      106     6441 
     900    4.9431751    34166     4083     2106      113     6501 
    1000     5.587605    33849     4078     2063       98     6482 
Loop time of 5.58763 on 1 procs for 1000 steps with 33849 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 3.6133     | 3.6133     | 3.6133     |   0.0 | 64.67
Coll    | 1.5171     | 1.5171     | 1.5171     |   0.0 | 27.15
Sort    | 0.14623    | 0.14623    | 0.14623    |   0.0 |  2.62
Comm    | 0.1807     | 0.1807     | 0.1807     |   0.0 |  3.23
Modify  | 0.12788    | 0.12788    | 0.12788    |   0.0 |  2.29
Output  | 0.00045085 | 0.00045085 | 0.00045085 |   0.0 |  0.01
Other   |            | 0.001997   |            |       |  0.04

Particle moves    = 30042153 (30M)
Cells touched     = 33831033 (33.8M)
Particle comms    = 0 (0K)
Boundary collides = 81925 (81.9K)
Boundary exits    = 133585 (0.134M)
SurfColl checks   = 5967352 (5.97M)
SurfColl occurs   = 100914 (0.101M)
Surf reactions    = 0 (0K)
Collide attempts  = 3365880 (3.37M)
Collide occurs    = 1761561 (1.76M)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 5.37654e+06
Particle-moves/step: 30042.2
Cell-touches/particle/step: 1.12612
Particle comm iterations/step: 1
Particle fraction communicated: 0
Particle fraction colliding with boundary: 0.002727
Particle fraction exiting boundary: 0.00444659
Surface-checks/particle/step: 0.198633
Surface-collisions/particle/step: 0.00335908
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0.112039
Collisions/particle/step: 0.0586363
Reactions/particle/step: 0

Particles: 33849 ave 33849 max 33849 min
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
