.. index:: compute gas/collision/tally

compute gas/collision/tally command
===================================

**Syntax:**


.. parsed-literal::

   compute ID gas/collision/tally group-ID mix-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* gas/collision/tally = style name of this compute command
* group-ID = group ID for which grid cells to perform calculation on
* mix-ID = mixture ID to perform calculation on
* one or more values can be appended
* value = *id/cell* or *id1* or *id2* or *type1* or *type2* or *vx1/pre* or *vy1/pre* or *vz1/pre* or *vx2/pre* or *vy2/pre* or *vz2/pre* or *vx1/post* or *vy1/post* or *vz1/post* or *vx2/post* or *vy2/post* or *vz2/post*
  
  .. parsed-literal::
  
       id/cell = grid cell ID
       id1\ , id2 = particle IDs of two particles
       type1\ , type2 = particle species of two particles
       vx1/pre\ , vy1/pre\ , vz1/pre = velocity components of first particle before collision
       vx2/pre\ , vy2/pre\ , vz2/pre = velocity components of second particle before collision
       vx1/post\ , vy1/post\ , vz1/post = velocity components of first particle after collision
       vx2/post\ , vy2/post\ , vz2/post = velocity components of second particle after collision



**Examples:**


.. parsed-literal::

   compute 1 gas/collision/tally all all id id/cell time xc yc zc

This command will dump the tallies in the previous command to a dump
file every 10 steps:


.. parsed-literal::

   dump 1 tally all 10 tmp.tally c_1[\*]

**Description:**

Tally various values for each gas-phase collision of two particles
during the current timestep.  Only gas collisions which do not result
in chemical reactions are tallied by this command.  See the related
:doc:`compute gas/reaction/tally <compute_gas_reaction_tally>` command
to tally collisions which induce reactions.  It provides different
tally value options than this compute.

Only collisions within grid cells in the grid group specified by
*group-ID* and pairs of particles with both species in the mixture
specified by *mix-ID* are included.  See the :doc:`group grid <group>`
command for info on how grid cells can be assigned to grid groups.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump tally <dump>` command.  Currently this is the only option for
accessing the tallied values.


----------


The *id/cell* value is the ID of the grid cell in which the collision
takes place.

The *id1* and *id2* values are the IDs of the colliding particles.
Note that particle IDs are generated randomly.  Thus multiple
particles in the system can potentially have the same ID.  See
`Section 6.19 <Section_howto.html#howto_19>`_ for more details on particle IDs.

The *type1* and *type2* values are the integer index for the particle
species of the colliding particles.  They are values from 1 to
Nspecies. The values corresponds to the order in which species were
defined via the :doc:`species <species>` command.  See `Section 6.19 <Section_howto.html#howto_19>`_ for more details on particle types.

The *vx1/pre*\ , *vy1/pre*\ , *vz1/pre*\ , *vx2/pre*\ , *vy2/pre*\ , *vz2/pre*
values are the velocity components of the two particles before the
collision.

The *vx1/post*\ , *vy1/post*\ , *vz1/post*\ , *vx2/post*\ , *vy2/post*\ ,
*vz2/post* values are the velocity components of the two particles
after the collision.


----------


**Output info:**

This compute calculates a per-tally array, with the number of columns
equal to the number of values.

The array can be accessed by any command that uses per-tally values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The per-tally array values will be in the :doc:`units <units>`
appropriate to the individual values as described above.  All the
velocity keywords are in velocity units.

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

:doc:`compute gas/reaction/tally <compute_gas_reaction_tally>`,
:doc:`compute surf/collision/tally <compute_surf_collision_tally>`,
:doc:`dump tally <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
