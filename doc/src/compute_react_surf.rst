.. index:: compute react/surf

compute react/surf command
==========================

compute react/surf/kk command
=============================

**Syntax:**


.. parsed-literal::

   compute ID react/surf group-ID reaction-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* react/surf = style name of this compute command
* group-ID = group ID for which surface elements to perform calculation on
* reaction-ID = surface reaction ID which defines surface reactions
* zero or more values can be appended
* value = *r:s1/s2/s3 ...* or *p:s1/s2/s3 ...*
  
  .. parsed-literal::
  
       r: or p: = list of reactant species or product species
       s1,s2,s3 = one or more species IDs, separated by "/" character



**Examples:**


.. parsed-literal::

   surf_react air prob air.surf
   compute 1 react/surf all air
   compute 2 react/surf all air r:N/O/N2/O2 p:N/O/NO

These commands will dump time averages for each surface element to a
dump file every 1000 steps:


.. parsed-literal::

   compute 2 react/surf all air r:N/O/N2/O2 p:N/O/NO
   fix 1 ave/surf all 10 100 1000 c_2[\*]
   dump 1 surf all 1000 tmp.surf id f_1[\*]

**Description:**

Define a computation that tallies counts of reactions for each
explicit surface element in a surface element group, based on the
particles that collide with that element.  Only surface elements in
the surface group specified by *group-ID* are included in the
tallying.  See the :doc:`group surf <group>` command for info on how
surface elements can be assigned to surface groups.  Likewise only
surface elements assigned to the surface reaction model specified by
*reaction-ID* are included in the tallying.

Explicit surface elements are triangles for 3d simulations and line
segments for 2d simulations.  Unlike implicit surface elements, each
explicit triangle or line segment may span multiple grid cells.  See
the :doc:`read\_surf <read_surf>` command for details.

This command can only be used for simulations with explicit surface
elements.  See the similar :doc:`compute react/isurf/grid <compute_react_isurf_grid>` command for use with
simulations with implicit surface elements.

Note that when a particle collides with a surface element, it can
bounce off (possibly as a different species), be captured by the
surface (vanish), or a 2nd particle can also be emitted.

The doc page for the :doc:`surf\_react <surf_react>` command explains the
different reactions that can occur for each specified style.

If no values are specified each reaction specified by the
:doc:`surf\_react <surf_react>` style is tallied individually for each
surface element.

If M values are specified, then M tallies are made for each surface
element, one per value.  If the value starts with "r:" then any
reaction which occurs with one (or more) of the listed species as a
reactant is counted as part of that tally.  If the value starts with
"p:" then any reaction which occurs with one (or more) of the listed
species as a product is counted as part of that tally.  Note that
these rules mean that a single reaction may be tallied multiple times
depending on which values it matches.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump surf <dump>` command.

The values over many sampling timesteps can be averaged by the :doc:`fix ave/surf <fix_ave_surf>` command.


----------


**Output info:**

This compute calculates a per-surf array, with the number of columns
either equal to the number of reactions defined by the
:doc:`surf\_react <surf_react>` style (if no values are specified) or equal to
M = the # of values specified.

Surface elements not in the specified *group-ID* or not assigned to
the specified *reaction-ID* will output zeroes for all their values.

The array can be accessed by any command that uses per-surf values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The per-surf array values are counts of the number of reactions that
occurred.


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

:doc:`fix ave/surf <fix_ave_surf>`, :doc:`dump surf <dump>`, :doc:`compute react/isurf/grid <compute_react_isurf_grid>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
