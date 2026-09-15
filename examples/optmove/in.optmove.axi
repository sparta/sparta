################################################################################
# 2d axisymmetric flow through a duct
#
# Exercises the optimized move on an axisymmetric model.  The fast path folds
# the end-of-step position back into the axisymmetric plane once for the whole
# timestep, where the normal move does it at every cell crossing; both give the
# same point, since each remap is a rotation applied to the position and the
# velocity together and the radial coordinate is unchanged by it.
#
# The axis itself needs no boundary condition, and an outflow face at the
# outer radius is handled on the fast path.  A reflecting outer radial face
# would not be, since that face is a cylinder -- see the global doc page.
#
# Run the same script with "-var opt no" to compare against the normal move.
# The optimized move performs one remap per timestep where the normal move
# performs one per cell crossing, so the two are identical while a particle
# stays within a cell over a timestep, as it does at the timestep used here,
# and differ in the last few digits once particles cross cells.
#
# Note:
#  - The "comm/sort" option to the "global" command is used to match MPI runs.
#  - The "twopass" option is used to match Kokkos runs.
################################################################################

variable            opt index yes

seed	    	    12345
dimension   	    2
global              gridcut 0.0 comm/sort yes optmove ${opt}

boundary	    o ao p

create_box          -0.25 0.25 0.0 0.25 -0.5 0.5
create_grid 	    20 10 1
balance_grid        rcb cell

global		    nrho 1.e20 fnum 1.e17 weight cell radius

species		    air.species N2
mixture		    air N2 vstream 3472.0 0.0 0.0 temp 300.0

fix                 in emit/face air xlo twopass
collide		    vss air air.vss

fix                 gcheck grid/check 1 error

stats		    100
stats_style	    step cpu np nattempt ncoll nbound nexit

timestep 	    1.e-6
run 		    500
