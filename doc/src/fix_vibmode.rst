.. index:: fix vibmode

fix vibmode command
===================

**Syntax:**


.. parsed-literal::

   fix ID vibmode

* ID is documented in :doc:`fix <fix>` command
* vibmode = style name of this fix command

**Examples:**


.. parsed-literal::

   fix 1 vibmode

**Description:**

Enable multiple vibrational energy levels, defined on a per-species
basis, to be used in a simulation.  This fix is meant to be used with
the :doc:`collide\_modify vibrate discrete <collide_modify>` setting
which means that the vibrational energy of each (non-monatomic)
particle is discretized across one or more energy modes, each with its
own characteristic vibrational temperature.  This fix allocates
per-particle storage for the mode indices and also has code to
populate the multiple levels appropriately when particles are created.
Collisions between pairs of particles will then transfer energy
between the different modes of the two particles.

An overview of how to run simulations with multiple vibrational energy
modes is given in the `Section 6.12 <Section_howto.html#howto_12>`_.
This includes use of the :doc:`species <species>` command with its
*vibfile* option, and the use of the :doc:`collide\_modify vibrate discrete <collide_modify>` command.  The section also lists all the
commands that can be used in an input script to invoke various options
associated with the vibrational energy modes.  All of them depend on
this fix vibmode command being defined.

Internally, this fix defines a custom particle attribute named
"vibmode".  It is an integer array with N values per particle.  N is
the maximum number of energy modes for any species defined in the
simulation.  The number of energy modes is half the vibrational
degrees of freedom defined for each species.  See the "species"
command for how the degrees of freedom and associated vibrational
temperatures and other properties are defined for each mode for each
species.

Each of the N values is an integer count for the


----------


**Restart, output info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

However, the values of the custom particle attribute defined by this
fix are written to the restart file.  Namely the integer values stored
in "vibmode" for each particle.  As explained on the
:doc:`read\_restart <read_restart>` doc page these values can be
re-assigned to particles when a restart file is read, if a new fix
vibmode command is specified in the restart script before the first
:doc:`run <run>` command is used.

No global or per-particle or per-grid quantities are stored by this
fix for access by various output commands.

However, the custom particle attributes defined by this fix can be
accessed by the :doc:`dump particle <dump>` command, as p\_vibmode.  That
means those per-particle values can be written to particle dump files.

**Restrictions:**

This fix is required if "collide\_modify vibrate discrete" is used and
there is one or more species defined which haave multiple vibrational
energy modes (2 or more).  In this scenario, if it is not defined, an
error will occur when a "create\_particles" or :doc:`run <run>` command
is issued.  Conversely, if no species has multiple vibrational modes,
this fix cannot be used.

Defining this fix after particles have been created will not populate
the vibrational energy modes of particles that already exist.  An
exception is if the :doc:`read\_restart <read_restart>` command is used
to read in particles from a previous simulation where this fix was
used.  In that case, defining this fix after reading the restart file
will enable the particles to keep their previous vibrational energy
mode values.

**Related commands:**

:doc:`collide\_modify vibrate discrete <collide_modify>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
