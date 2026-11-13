.. index:: scale\_particles

scale\_particles command
========================

**Syntax:**


.. parsed-literal::

   scale_particles mix-ID factor

* mix-ID = ID of mixture to use when scaling particles
* factor = scale factor

**Examples:**


.. parsed-literal::

   scale_particles air 0.5
   scale_particles air 4.0

**Description:**

Scale the number of particles in the simulation by cloning or deleting
individual particles.  This can be useful between runs, or after
reading a restart file, to increase or decrease the particle count
before a new :doc:`run <run>` command is issued, as if the :doc:`global fnum <global>` value had been changed.  For example, an initial
coarse simulation can be performed, followed by a simulation at
higher resolution.

Only particles of species in the specified mixture are considered for
cloning/deleting.  See the :doc:`mixture <mixture>` command for how it
defines a collection of species.

The specified *factor* can be any value >= 0.0.

If *factor* < 1.0, then for each particle, a random number *R* is
generated.  If R > factor, the particle is deleted.

If *factor* > 1.0, then for each particle additional particles may be
created, by cloning all attributes of the original particle, except
for a new random particle ID assigned to each new particle.  E.g. if
*factor* = 3.4, then two extra particles are created, and a 3rd is
created with probability 0.4.


----------


**Restrictions:** none

**Related commands:**

:doc:`create\_particles <create_particles>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
