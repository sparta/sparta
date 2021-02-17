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
WARNING: Could not acquire nearby ghost cells b/c grid partition is not clumped (../grid.cpp:410)
Created 200 child grid cells
  CPU time = 0.00109696 secs
  create/ghost percent = 81.8518 18.1482
balance_grid        rcb cell
Balance grid migrated 120 cells
  CPU time = 0.000472069 secs
  reassign/sort/migrate/ghost percent = 52.9798 0.40404 20.7576 25.8586

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
  CPU time = 0.00135303 secs
  read/check/sort/surf2grid/ghost/inout/particle percent = 21.8678 24.6167 0.740088 45.3744 7.40088 20.8458 1.9207
  surf2grid time = 0.000613928 secs
  map/comm1/comm2/comm3/comm4/split percent = 18.7573 8.46602 5.35922 5.20388 6.99029 45.4369

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
     100   0.14855003    18434     1529      823       95     4797 
     200     0.366045    27306     2857     1487      106     5789 
     300   0.59696698    31002     3516     1833      117     6275 
     400   0.83822703    32478     3698     1989       94     6097 
     500     1.070297    33507     3897     2034       95     6527 
     600    1.3129559    33571     3922     2011       96     6319 
     700     1.546097    33735     3944     2017      119     6464 
     800    1.7792599    33995     4005     2135      114     6379 
     900      2.01068    34129     4024     2097      100     6621 
    1000    2.2533979    34105     4072     2102      123     6899 
Loop time of 2.25353 on 4 procs for 1000 steps with 34105 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 0.66999    | 0.93647    | 1.3935     |  30.4 | 41.56
Coll    | 0.23539    | 0.36857    | 0.52473    |  17.7 | 16.36
Sort    | 0.027198   | 0.034981   | 0.047678   |   4.1 |  1.55
Comm    | 0.081746   | 0.093183   | 0.1187     |   4.9 |  4.14
Modify  | 0.00049782 | 0.032667   | 0.12911    |  30.8 |  1.45
Output  | 0.00087714 | 0.002323   | 0.0036592  |   2.2 |  0.10
Other   |            | 0.7853     |            |       | 34.85

Particle moves    = 29904535 (29.9M)
Cells touched     = 33680543 (33.7M)
Particle comms    = 403718 (0.404M)
Boundary collides = 81541 (81.5K)
Boundary exits    = 133149 (0.133M)
SurfColl checks   = 5958956 (5.96M)
SurfColl occurs   = 101047 (0.101M)
Surf reactions    = 0 (0K)
Collide attempts  = 3347338 (3.35M)
Collide occurs    = 1755152 (1.76M)
Reactions         = 0 (0K)
Particles stuck   = 0

Particle-moves/CPUsec/proc: 3.31753e+06
Particle-moves/step: 29904.5
Cell-touches/particle/step: 1.12627
Particle comm iterations/step: 2.233
Particle fraction communicated: 0.0135002
Particle fraction colliding with boundary: 0.00272671
Particle fraction exiting boundary: 0.00445247
Surface-checks/particle/step: 0.199266
Surface-collisions/particle/step: 0.00337899
Surf-reactions/particle/step: 0
Collision-attempts/particle/step: 0.111934
Collisions/particle/step: 0.0586918
Reactions/particle/step: 0

Particles: 8526.25 ave 11115 max 6760 min
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
