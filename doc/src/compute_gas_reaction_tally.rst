.. index:: compute reaction/tally

compute reaction/tally command
==============================

**Syntax:**


.. parsed-literal::

   compute ID reaction/tally group-ID mix-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* reaction/tally = style name of this compute command
* group-ID = group ID for which surface elements to perform calculation on
* mix-ID = mixture ID to perform calculation on
* one or more values can be appended
* value = *id/pre* or *id/post* or *id/post2* or *type/pre* or *type/post* or *type/post2* or *id/surf* or *time* or *xc* or *yc* or *zc* or *vx/pre* or *vy/pre* or *vz/pre* or *vx/post* or *vy/post* or *vz/post* or *vx/post2* or *vy/post2* or *vz/post2*
  
  .. parsed-literal::
  
       reaction = which reaction occurred (1 to N)
       id/cell = grid cell ID
       id1/pre = particle ID of first particle before reaction
       id2/pre = particle ID of second particle before reaction
       id1/post = particle ID of first particle after reaction
       id2/post = particle ID of second particle after reaction
       id3/post = particle ID of third particle after reaction
       type1/pre = particle species of first particle before reaction
       type2/pre = particle species of second particle before reaction
       type1/post = particle species of first particle before reaction
       type2/post = particle species of second particle before reaction
       type3/post = particle species of third particle before reaction
       vx1/pre\ , vy1/pre\ , vz1/pre = velocity components of first particle before reaction
       vx2/pre\ , vy2/pre\ , vz2/pre = velocity components of second particle before reaction
       vx1/post\ , vy1/post\ , vz1/post = velocity components of first particle after reaction
       vx2/post\ , vy2/post\ , vz2/post = velocity components of second particle after reaction
       vx3/post\ , vy3/post\ , vz3/post = velocity components of thrid particle after reaction



**Examples:**


.. parsed-literal::

   compute 1 reaction/tally all all id/pre id/post id/post2 type/pre type/post type/post2 id/surf time xc yc zc

This command will dump the tallies in the previous command to a dump
file every 10 steps:


.. parsed-literal::

   dump 1 tally all 10 tmp.tally c_1[\*]

**Description:**

Tally various values for each gas-phase reaction of two particles
during the current timestep.  Only gas collisions which result in a
chemical reactions are tallied by this command, as enabled by the
:doc:`react <react>` command.  See the related :doc:`compute gas/collide/tally <compute_gas_collision_tally>` command to tally
collisions which do not induce reactions.  It provides different tally
value options than this compute.

Only collisions within grid cells in the grid group specified by
*group-ID* and pairs of particles with both species in the mixture
specified by *mix-ID* are included.  See the :doc:`group grid <group>`
command for info on how grid cells can be assigned to grid groups.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump tally <dump>` command.  Currently this is the only option for
accessing the tallied values.

Note that if a gas-phase reaction takes place then a collision of two
particles can result in 1, 2, or 3 particles being produced.  The
species of the 1, 2, or product particles may be different than the
species of the reactants.  All of these attributes can be tallied
using the values explained next.


----------


The *reaction* value is the index of which reaction occurred.  It is a
number from 1 to N, where N is the number of reactions defined by the
:doc:`react <react>` command.

The *id1/pre* and *id2/pre* values are the IDs of the reacting
particles before the reaction.  Note that particle IDs are generated
randomly.  Thus multiple particles in the system can potentially have
the same ID.  See `Section 6.19 <Section_howto.html#howto_19>`_ for more details on
particle IDs.

The *id1/post*\ , *id2/post*\ , and *id3/post* values are the
post-reaction IDs of the 1, 2, or 3 resulting particles.  If one of
the particles vanished in the reaction, then *id2/post* will be 0.  If
*id3/post* is 0, no 3rd particle was created by the reaction,
otherwise it is the random ID assigned to the newly created particle.

The *type1/pre* and *type2/pre* value are the integer indices for the
species of the reacting particles before the reaction.  They are
values from 1 to Nspecies.  The values correspond to the order in
which species were defined via the :doc:`species <species>` command.
See `Section 6.19 <Section_howto.html#howto_19>`_ for more details on particle
types.

The *type1/post*\ , *type/post2*\ , and *type/post3* values are the
post-collision integer indices for the particle species of the 0, 1,
or 2 resulting particles.  If one of the particles vanished in the
reaction, then *type2/post* will be 0.  If *type/post3* is 0, no 3rd
particle was created by the reaction, otherwise it is the species
index of the newly created particle.

The values starting with *v* for particles 1 and 2 with a *pre* suffix
are the velocity components of the two reacting particles before the
reaction.

The values starting with *v* for particles 1 and 2 and 3 with a *post*
suffix are the velocity components of the 1, 2, or 3 resulting
particles after the reaction.  Velocity components will be zero if the
2nd or 3rd particle does not exist due to the reaction.


----------


**Output info:**

This compute calculates a per-tally array, with the number of columns
equal to the number of values.

The array can be accessed by any command that uses per-tally values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The per-tally array values will be in the :doc:`units <units>`
appropriate to the individual values as described above.  All the
velocity components are in velocity units.

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

:doc:`compute gas/collision/tally <compute_gas_collision_tally>`, :doc:`dump tally <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
