.. index:: compute react/boundary

compute react/boundary command
==============================

compute react/boundary/kk command
=================================

**Syntax:**


.. parsed-literal::

   compute ID react/boundary reaction-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* react/boundary = style name of this compute command
* reaction-ID = surface reaction ID which defines surface reactions
* zero or more values can be appended
* value = *r:s1/s2/s3 ...* or *p:s1/s2/s3 ...*
  
  .. parsed-literal::
  
       r: or p: = list of reactant species or product species
       s1,s2,s3 = one or more species IDs, separated by "/" character



**Examples:**


.. parsed-literal::

   surf_react air prob air.surf
   compute 1 react/boundary air
   compute 2 react/boundary air r:N/O/N2/O2 p:N/O/NO

These commands will time average the reaction tallies for each face
and output the results as part of statistical output:


.. parsed-literal::

   compute 2 react/boundary air r:N/O/N2/O2 p:N/O/NO

   fix 1 ave/time all 10 100 1000 c_2[\*]
   stats_style step np f_1[1][\*] f_1[2][\*] f_1[3][\*] f_1[4][\*]

**Description:**

Define a computation that tallies counts of reactions for each
boundary (i.e. face) of the simulation box, based on the particles
that collide with the boundary.  Only faces assigned to the surface
reaction model specified by *reaction-ID* are included in the
tallying.

Note that when a particle collides with a face, it can bounce off
(possibly as a different species), be captured by the surface
(vanish), or a 2nd particle can also be emitted.

The doc page for the :doc:`surf\_react <surf_react>` command explains the
different reactions that can occur for each specified style.

If no values are specified each reaction specified by the
:doc:`surf\_react <surf_react>` style is tallied individually for each
boundary.

If M values are specified, then M tallies are made for each face, one
per value.  If the value starts with "r:" then any reaction which
occurs with one (or more) of the listed species as a reactant is
counted as part of that tally.  If the value starts with "p:" then any
reaction which occurs with one (or more) of the listed species as a
product is counted as part of that tally.  Note that these rules mean
that a single reaction may be tallied multiple times depending on
which values it matches.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`stats\_style <stats_style>` command.  The values over many sampling
timesteps can be averaged by the :doc:`fix ave/time <fix_ave_time>`
command.


----------


**Output info:**

This compute calculates a global array, with the number of columns
either equal to the number of reactions defined by the
:doc:`surf\_react <surf_react>` style (if no values are specified) or equal to
M = the # of values specified.  The number of rows is 4 for a 2d
simulation for the 4 faces (xlo, xhi, ylo, yhi), and it is 6 for a 3d
simulation (xlo, xhi, ylo, yhi, zlo, zhi).

The array can be accessed by any command that uses global array values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The array values are counts of the number of reactions that occurred
on each face.

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

:doc:`fix ave/time <fix_ave_time>`, :doc:`compute react/surf <compute_react_surf>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
