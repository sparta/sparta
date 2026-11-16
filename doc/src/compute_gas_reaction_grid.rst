.. index:: compute gas/reaction/grid

compute gas/reaction/grid command
=================================

**Syntax:**


.. parsed-literal::

   compute ID gas/reaction/grid group-ID mix-ID mode value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* gas/collision/grid = style name of this compute command
* group-ID = group ID for which grid cells to perform calculation on
* mix-ID = mixture ID to perform calculation on
* mode = *all* or *every* or *select*
  
  .. parsed-literal::
  
       all = single count of all reactions in each grid cell
         values = none
       every = count of each defined reaction in each grid cell
         values = none
       select = count of selected reactions in each grid cell
         values = one or more numeric indices from 1 to M where M = # of defined reactions



**Examples:**


.. parsed-literal::

   compute 1 gas/reaction/grid all all all
   compute 2 gas/reaction/grid subset mymixture all
   compute 4 gas/reaction/grid all all every
   compute 4 gas/reaction/grid all all select 5 7 14

**Description:**

Count the number of gas-phase reactions bewteen pairs of particles
which occur in each grid cell during the current timestep.  Only gas
collisions which result in chemical reactions are counted by this
command, as defined by the :doc:`react <react>` command and the list of
reactions it reads from a file.  See the related :doc:`compute gas/collision/grid <compute_gas_reaction_grid>` command to count
collisions which do not induce reactions.

Only reactions within grid cells in the grid group specified by
*group-ID* and pairs of reactant particles with both species in the
mixture specified by *mix-ID* are included.  See the :doc:`group grid <group>` command for info on how grid cells can be assigned to
grid groups.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump grid <dump>` command or used as inputs to the :doc:`compute reduce <compute_reduce>` command.  The values can also be time
averaged by the :doc:`fix ave/grid <fix_ave_grid>` command.


----------


The *mode* argument determines how the reactions in each grid cell are
counted.

If *mode* is specified as *all* then a single count of all reactions
is produced for each grid cell.  The compute calculates this count as
a per-grid vector.

If *mode* is specified as *every* then an individual count for each
reaction is produced for each grid cell.  The number of possible
reactions M are those listed in the file read by the
:doc:`react <react>` commmand.  The compute calculates these counts as a
per-grid array, even if M = 1.

If *mode* is specified as *select* then a count for each reaction in
the specified list of numeric indices is produced for each grid cell.
Each index must be a number from 1 to M, where M is the number of
reactions listed in the file read by the :doc:`react <react>` commmand.
The compute calculates these counts as a per-grid array, even if only
a single index is specified.


----------


**Output info:**

If no numeric keywords are used, then this compute calculates a
per-grid vector with the count of all reactions for each grid cell.

If one or more numeric keywords are used, then this compute calculates
a per-grid array.  The first column wlll be the count of all reactions.
The remaining columns (one per numeric keyword) will be the counts of only
the specific reactions listed.

The vector or array can be accessed by any command that uses per-grid
values from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_ for an overview of SPARTA output
options.

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

:doc:`compute gas/collision/grid <compute_gas_collision_grid>`, :doc:`compute gas/reaction/tally <compute_gas_reaction_tally>`, :doc:`dump grid <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
