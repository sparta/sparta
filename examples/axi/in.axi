# 2d axisymmetric flow around a circle with specular reflections

seed	    	    12345
dimension   	    2
global              gridcut 0.0

boundary	    o ar p

create_box          -0.25 0.5 0.0 0.5 -0.5 0.5

create_grid 	    30 20 1 
balance_grid        rcb cell

global		    nrho 1.e20 fnum 1.e17 weight cell radius

species		    air.species N2
mixture		    air N2 vstream 3472.0 0.0 0.0 temp 300.0

fix                 in emit/face air xlo
collide		    vss air air.vss

read_surf           data.circle origin 5 5 0 &
                    trans -5 -5 0 scale 1.666666e-2 1.666666e-2 1 clip

surf_collide	    1 specular
surf_modify         all collide 1

create_particles    air n 100000

timestep 	    1e-6

#dump                2 image all 100 tmp.*.ppm type type pdiam 0.001 &
#                    surf proc 0.0 &
#		    size 512 512 axes yes 0.9 0.02 particle yes &
#                    gline yes 0.005 &
#                    surf proc 0.005 zoom 4.0
#dump_modify	    2 pad 4 

stats		    100
stats_style	    step cpu np nattempt ncoll nscoll nscheck

run 		    1000

