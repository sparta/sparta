.. index:: read\_particles

read\_particles command
=======================

**Syntax:**


.. parsed-literal::

   read_particles file Nstep

* file = dump file to read snapshot from
* Nstep = timestep to read

**Examples:**


.. parsed-literal::

   read_particles dump.sphere 10500

**Description:**

Read a snapshot of particles from a previously created dump file and
add them to the simulation domain.  This is a means of reading in
particles from a previous SPARTA simulation or created as output by
another code.  The :doc:`create\_particles <create_particles>`, :doc:`fix emit/face <fix_emit_face>`, and :doc:`read\_restart <read_restart>`
commands are alternate ways to generate particles for a simulation.

The dump file must be in the SPARTA format created by the :doc:`dump particles <dump>` command which is described on its doc page.

Currently, each line of particle data in the file must have 8 fields
in the following order.  At some point we may generalize this format.


.. parsed-literal::

   id, type, x, y, z, vx, vy, vz

The *id* is any positive integer, which can simply be set to values
from 1 to Nparticles if desired.  The type is the species ID from 1 to
Nspecies.  The value corresponds to the order in which species are
defined in the current input script via the :doc:`species <species>`
command.  The x,y,z values are the particle coordinates which must be
inside (or on the surface of) the simulation box.  If a particle is
outside the box it will be skipped when the file is read.  For 2d or
axisymmetric simulations z = 0.0 should be used, though SPARTA does
not check for this.  The vx,vy,vz values are the particle velocity.
The rotational and vibrational energies for the new particles are set
to 0.0.

When the reading of particles is complete, the number of particles
read is printed to the screen.  If the number is smaller than the
particles in the file, it is because some were outside the simulation
box.

A check is made for any particle inside a surface object which
triggers an error.  However the check is only for grid cells entirely
inside a surface object.  Particles in grid cells which are cut by
surfaces are not checked.  It is your responsibility to insure
particles close to surfaces are actually outside the surface object.
If this is not the case, errors may be triggered once particles begin
to move.


----------


**Restrictions:** none

**Related commands:**

:doc:`create\_particles <create_particles>`, :doc:`fix emit/face <fix_emit_face>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
