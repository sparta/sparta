.. index:: compute telec/grid

compute telec/grid command
==========================

compute telec/grid/kk command
=============================

**Syntax:**


.. parsed-literal::

   compute ID telec/grid group-ID mix-ID

* ID is documented in :doc:`compute <compute>` command
* telec/grid = style name of this compute command
* group-ID = group ID for which grid cells to perform calculation on
* mix-ID = mixture ID to perform calculation on

**Examples:**


.. parsed-literal::

   compute 1 telec/grid all species
   compute 1 telec/grid subset all

**Description:**

Define a computation that calculates the electronic temperature for
each grid cell in a grid cell group, based on the particles in the
cell.  How the electronic temperature is computed is explained below.
The temperature is calculated separately for each group of species in
the specified mixture, as described in the Output section below.  See
the mixture command for how a set of species can be partitioned into
groups.

Only grid cells in the grid group specified by *group-ID* are included
in the calculations.  See the :doc:`group grid <group>` command for info
on how grid cells can be assigned to grid groups.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump grid <dump>` command.

The values over many sampling timesteps can be averaged by the :doc:`fix ave/grid <fix_ave_grid>` command.  It does its averaging as if the
particles in the cell at each sampling timestep were combined together
into one large set to compute the formulas below.  Note that this is a
different normalization than taking the values produced by the
formulas below for a single timestep, summing them over the sampling
timesteps, and then dividing by the number of sampling steps.


----------


The electronic temperature in a grid cell for a group of particles
comprised of different species is defined as a weighted average as
follows:


.. parsed-literal::

   T_group = (T1\*N1 + T2\*N2 + ...) / (N1 + N2 + ...)

The electronic temperature Tsp for particles of a single species is
defined by calculating the average electronic energy per particle.
Then an electronic temperature is solved for where the average
per-particle electronic energy at that equilibrium temperature matches
the value calculated from the particles in the grid cell.


----------


**Output info:**

This compute calculates a per-grid array. The number of columns is
equal to the number of groups in the specified mixture.

This compute performs calculations for all flavors of child grid cells
in the simulation, which includes unsplit, cut, split, and sub cells.
See `Section 6.8 <Section_howto.html#howto_8>`_ of the manual gives
details of how SPARTA defines child, unsplit, split, and sub cells.
Note that cells inside closed surfaces contain no particles.  These
could be unsplit or cut cells (if they have zero flow volume).  Both
of these kinds of cells will compute a zero result for all their
values.  Likewise, split cells store no particles and will produce a
zero result.  This is because their sub-cells actually contain the
particles that are geometrically inside the split cell.

Grid cells not in the specified *group-ID* will output zeroes for all
their values.

The array can be accessed by any command that uses per-grid values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The per-grid array values will be in temperature :doc:`units <units>`.


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

:doc:`compute grid <compute_grid>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
