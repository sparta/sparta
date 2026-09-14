.. index:: compute surf/reaction/tally

compute surf/reaction/tally command
===================================

**Syntax:**


.. parsed-literal::

   compute ID surf/reaction/tally group-ID mix-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* surf/reaction/tally = style name of this compute command
* group-ID = group ID for which surface elements to perform calculation on
* mix-ID = mixture ID to perform calculation on
* one or more values can be appended
* value = *reaction* or *id/surf* or *id/pre* or *id1/post* or *id2/post* or *type/pre* or *type1/post* or *type2/post* or *time* or *xc* or *yc* or *zc* or *vx/pre* or *vy/pre* or *vz/pre* or *vx1/post* or *vy1/post* or *vz1/post* or *vx2/post* or *vy2/post* or *vz2/post*
  
  .. parsed-literal::
  
       reaction = which reaction occurred (1 to N)
       id/surf = surface element ID
       id/pre = particle ID before reaction
       id1/post = particle ID of first particle after reaction
       id2/post = particle ID of second particle after reaction
       type/pre = particle species before reaction
       type1/post = particle species of first particle after reaction
       type2/post = particle species of second particle after reaction
       time = time of collision (within timestep, see below)
       xc\ , yc\ , zc = coordinates of collision point
       vx/pre\ , vy/pre\ , vz/pre = velocity components before reaction
       vx1/post\ , vy1/post\ , vz1/post = velocity components of first particle after reaction
       vx2/post\ , vy2/post\ , vz2/post = velocity components of second particle after reaction



**Examples:**


.. parsed-literal::

   compute 1 surf/reaction/tally all all reaction id/pre id1/post id2/post type/pre type1/post type2/post id/surf time xc yc zc

This command will dump the tallies in the previous command to a dump
file every 10 steps:


.. parsed-literal::

   dump 1 tally all 10 tmp.tally c_1[\*]

**Description:**

Tally various values for each reaction of a particle with a surface
element during the current timestep.  Only surface collisions which
result in a chemical reaction are tallied by this command, as enabled
by the :doc:`surf\_react <surf_react>` command.  See the related :doc:`compute surf/collision/tally <compute_surf_collision_tally>` command to
tally collisions which do not induce reactions.  It provides different
tally value options than this compute.

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

Note that if a gas/surface reaction takes place, then a collision of a
single particle with a surface can result in 0, 1, or 2 particles
being produced.  The species of the product particles may be different
than the species of the reactant.  All of these attributes can be
tallied using the values explained next.


----------


The *reaction* value is the index of which reaction occurred.  It is a
number from 1 to N, where N is the number of reactions defined by the
:doc:`surf\_react <surf_react>` command.

The *id/surf* value is the ID of the surface element.  For explicit
surfaces, surface element IDs are unique.  For implicit surfaces, the
surface ID of all the surface elements within a grid cell are the grid
cell ID.  See `Section 6.13 <Section_howto.html#howto_13>`_ for more discussion of
explicit and implicit surfaces.

The *id/pre* value is the ID of the incident particle.  Note that
particle IDs are generated randomly.  Thus multiple particles in the
system can potentially have the same ID.  See `Section 6.19 <Section_howto.html#howto_19>`_ for more details on particle IDs.

The *id1/post* and *id2/post* values are the post-collision IDs of the
0, 1, or 2 resulting particles.  If the incident particle vanished in
the reaction, then *id1/post* will be 0.  Otherwise it will be the
same as *id/pre*\ .  If *id2/post* is 0, no 2nd particle was created by
the reaction, otherwise it is the random ID assigned to the newly
created particle.

The *type/pre* value is the integer index for the incident particle
species.  It is a value from 1 to Nspecies. The value corresponds to
the order in which species were defined via the :doc:`species <species>`
command.  See `Section 6.19 <Section_howto.html#howto_19>`_ for more details on
particle types.

The *type1/post* and *type2/post* values are the post-collision
integer indices for the particle species of the 0, 1, or 2 resulting
particles.  If the incident particle vanished in the reaction, then
*type1/post* will be 0.  Otherwise *type1/post* is the new species
index of the incident particle.  If *type2/post* is 0, no 2nd particle
was created by the reaction, otherwise it is the species index of the
newly created particle.

The *time* value is the point in time within the timestep when the
collision occured.  Thus if the timestep size is DT, as set by the
:doc:`timestep <timestep>` command, the *time* values will all be within
the range of 0.0 to DT, for all collisions tallied by this compute on
each timestep it is invoked.

The *xc*\ , *yc*\ , and *zc* are the coordinates of the collision point
between the particle and the surface element.

The *vx/pre*\ , *vy/pre*\ , *vz/pre* values are the velocity components of
the incident particle before it collides with the surface element.

The *vx1/post*\ , *vy1/post*\ , *vz1/post* values are the velocity
components of the incident particle after it collides with the surface
element and has reacted.  These values will be zero if the paticle
vanished in the reaction.

If a 2nd particle was created by the reaction, the *vx2/post*\ ,
*vy2/post*\ , *vz2/post* values are its velocity components.  If no 2nd
particle was created, these values will be zero.


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
distance units.  *Vx/pre*\ , *vy/pre*\ , *vz/pre* are in velocity units,
as are *vx1/post*\ , *vy1/post*\ , *vz1/post* and *vx2/post*\ , *vy2/post*\ ,
*vz2/post*\ .

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

:doc:`compute surf/collision/tally <compute_surf_collision_tally>`,
:doc:`compute gas/reaction/tally <compute_gas_reaction_tally>`,
:doc:`dump tally <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
