.. index:: compute gas/collision/grid

compute gas/collision/grid command
==================================

**Syntax:**


.. parsed-literal::

   compute ID gas/collision/grid group-ID mix-ID

* ID is documented in :doc:`compute <compute>` command
* gas/collision/grid = style name of this compute command
* group-ID = group ID for which grid cells to perform calculation on
* mix-ID = mixture ID to perform calculation on

**Examples:**


.. parsed-literal::

   compute 1 gas/collision/grid all all
   compute 2 gas/collision/grid subset mymixture

**Description:**

Count the number of gas-phase collisions bewteen pairs of particles
which occur in each grid cell during the current timestep.  Only gas
collisions which do not result in chemical reactions are counted by
this command.  See the related :doc:`compute gas/reaction/grid <compute_gas_reaction_grid>` command to count
collisions which induce reactions.

Only collisions within grid cells in the grid group specified by
*group-ID* and pairs of particles with both species in the mixture
specified by *mix-ID* are included.  See the :doc:`group grid <group>`
command for info on how grid cells can be assigned to grid groups.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump grid <dump>` command or used as inputs to the :doc:`compute reduce <compute_reduce>` command.  The values can also be time
averaged by the :doc:`fix ave/grid <fix_ave_grid>` command.


----------


**Output info:**

This compute calculates a per-grid vector with the count of collisions
for each grid cell.

The vector can be accessed by any command that uses per-grid values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

----------


Styles with a *kk* suffix are functionally the same as the
corresponding style without the suffix.  They have been optimized to
run faster, depending on your available hardware, as discussed in the
:doc:`Accelerating SPARTA <Section_accelerate>` section of the manual.
The accelerated styles take the same arguments and should produce the
same results, except for different random number, round-off and
precision issues.

These accelerated styles are part of the KOKKOS package. They are only
enabled if SPARTA was built with that package.  See the `Making SPARTA <Section_start.html#start_3>`_ section for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the `-suffix command-line switch <Section_start.html#start_7>`_ when you invoke SPARTA, or you can
use the :doc:`suffix <suffix>` command in your input script.

See the :doc:`Accelerating SPARTA <Section_accelerate>` section of the
manual for more instructions on how to use the accelerated styles
effectively.


----------


**Restrictions:** none

**Related commands:**

:doc:`compute gas/reaction/grid <compute_gas_reaction_grid>`, :doc:`compute gas/collision/tally <compute_gas_collision_tally>`, :doc:`dump grid <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
