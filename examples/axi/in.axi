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

create_grid 	    20 10 1 
balance_grid        rcb cell

global		    nrho 1.e20 fnum 1.e17 weight cell radius

species		    air.species N2
mixture		    air N2 vstream 3472.0 0.0 0.0 temp 300.0

fix                 in emit/face air xlo twopass
collide		    vss air air.vss

read_surf           data.circle origin 5 5 0 &
                    trans -5 -5 0 scale 0.05 0.05 1 clip

surf_collide	    1 specular
surf_modify         all collide 1

timestep 	    1e-6

#dump                2 image all 50 image.*.ppm type type pdiam 0.002 &
#	            size 512 512 particle yes &
#                    gline yes 0.01 &
#                    surf proc 0.02 zoom 3.5
#dump_modify	    2 pad 4 

stats		    100
stats_style	    step cpu np nattempt ncoll nscoll nscheck

run 		    1000
