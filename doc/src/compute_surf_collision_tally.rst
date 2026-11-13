.. index:: compute surf/collision/tally

compute surf/collision/tally command
====================================

**Syntax:**


.. parsed-literal::

   compute ID surf/collision/tally group-ID mix-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* surf/collision/tally = style name of this compute command
* group-ID = group ID for which surface elements to perform calculation on
* mix-ID = mixture ID to perform calculation on
* one or more values can be appended
* value = *id/surf* or *id* or *type* or *time* or *xc* or *yc* or *zc* or *vx/pre* or *vy/pre* or *vz/pre* or *vx/post* or *vy/post* or *vz/post*
  
  .. parsed-literal::
  
       id/surf = surface element ID
       id = particle ID
       type = particle species
       time = time of collision (within timestep, see below)
       xc\ , yc\ , zc = coordinates of collision point
       vx/pre\ , vy/pre\ , vz/pre = velocity components before collision
       vx/post\ , vy/post\ , vz/post = velocity components after collision



**Examples:**


.. parsed-literal::

   compute 1 surf/collision/tally all all id id/surf time xc yc zc

This command will dump the tallies in the previous command to a dump
file every 10 steps:


.. parsed-literal::

   dump 1 tally all 10 tmp.tally c_1[\*]

**Description:**

Tally various values for each collision of a particle with a surface
element during the current timestep.  Only surface collisions which do
not result in chemical reactions are tallied by this command.  See the
related :doc:`compute surf/reaction/tally <compute_surf_reaction_tally>`
command to tally collisions which induce reactions.  It provides
different tally value options than this compute.

Only surface elements in the surface group specified by *group-ID* and
particles with species in the mixture specified by *mix-ID* are
included.  See the :doc:`group grid <group>` command for info on how
surface elements can be assigned to surface groups.

This compute can be used for simulations with either explicit or
implicit surface elements.  See `Section 6.13 <Section_howto.html#howto_13>`_ for
more discussion of explicit and implicit surfaces.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump tally <dump>` command.  Currently this is the only option for
accessing the tallied values.


----------


The *id/surf* value is the ID of the surface element.  For explicit
surfaces, surface element IDs are unique.  For implicit surfaces, the
surface ID of all the surface elements within a grid cell are the grid
cell ID.  See `Section 6.13 <Section_howto.html#howto_13>`_ for more discussion of
explicit and implicit surfaces.

The *id* value is the ID of the particle.  Note that particle IDs are
generated randomly.  Thus multiple particles in the system can
potentially have the same ID.  See `Section 6.19 <Section_howto.html#howto_19>`_ for
more details on particle IDs.

The *type* value is an integer index for the particle species.  It is
a value from 1 to Nspecies. The value corresponds to the order in
which species were defined via the :doc:`species <species>` command.
See `Section 6.19 <Section_howto.html#howto_19>`_ for more details on particle
types.

The *time* value is the point in time within the timestep when the
collision occured.  Thus if the timestep size is DT, as set by the
:doc:`timestep <timestep>` command, the *time* values will all be within
the range of 0.0 to DT, for all collisions tallied by this compute on
each timestep it is invoked.

The *xc*\ , *yc*\ , and *zc* are the coordinates of the collision point
between the particle and the surface element.

The *vx/pre*\ , *vy/pre*\ , *vz/pre* values are the velocity components of
the incident particle before it collides with the surface element.

The *vx/post*\ , *vy/post*\ , *vz/post* values are the velocity components
of the particle after it collides with the surface element.


----------


**Output info:**

This compute calculates a per-tally array, with the number of columns
equal to the number of values.

The array can be accessed by any command that uses per-tally values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The per-tally array values will be in the :doc:`units <units>`
appropriate to the individual values as described above.  *Time* is in
time units; see above for details.  *Xc*\ , *yc*\ , and *zc* are in
distance units.  *Vx/pre*\ , *vy/pre*\ , *vz/pre* are in velocity units, as
are *vx/post*\ , *vy/post*\ , *vz/post*\ .

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

:doc:`compute surf/reaction/tally <compute_surf_reaction_tally>`,
:doc:`compute gas/collision/tally <compute_gas_collision_tally>`,
:doc:`dump tally <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
