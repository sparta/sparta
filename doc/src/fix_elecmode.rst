.. index:: fix elecmode

fix elecmode command
====================

fix elecmode command/kk
=======================

**Syntax:**


.. parsed-literal::

   fix ID elecmode

* ID is documented in :doc:`fix <fix>` command
* elecmode = style name of this fix command

**Examples:**


.. parsed-literal::

   fix 1 elecmode

**Description:**

Enable electronic energy states, defined on a per-species
basis, to be used in a simulation.  This fix is meant to be used with
the :doc:`collide\_modify electronic discrete <collide_modify>` setting.
This fix allocates per-particle storage for the electronic state.
Collisions between pairs of particles will then transfer energy
between the different electronic states of the two particles.

Internally, this fix defines custom particle attributes named
"elecstate" and "eelec". They are integer and double typed arrays
with a single value for each particle. See the "species"
command for how the electronic states are defined for each
species.


----------


**Restart, output info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

However, the values of the custom particle attribute defined by this
fix are written to the restart file. As explained on the
:doc:`read\_restart <read_restart>` doc page these values can be
re-assigned to particles when a restart file is read, if a new fix
elecmode command is specified in the restart script before the first
:doc:`run <run>` command is used.

The per-species electronic state tables read via the *elecfile* option
of the :doc:`species <species>` command are also stored in the restart
file, so a restarted run can continue performing electronic relaxation
(and reactions) without re-issuing the species command; only the
:doc:`collide <collide>`, :doc:`collide\_modify <collide_modify>` and fix
elecmode commands need to be repeated in the restart script.

No global or per-particle or per-grid quantities are stored by this
fix for access by various output commands.

However, the custom particle attributes defined by this fix can be
accessed by the :doc:`dump particle <dump>` command, as p\_elecstate and p\_eelec.
That means those per-particle values can be written to particle dump files.


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


**Restrictions:**

This fix is required if "collide\_modify electronic discrete" is used.
In this scenario, if it is not defined, an
error will occur when a "create\_particles" or :doc:`run <run>` command
is issued.

Defining this fix after particles have been created will not populate
the electronic energy modes of particles that already exist.  An
exception is if the :doc:`read\_restart <read_restart>` command is used
to read in particles from a previous simulation where this fix was
used.  In that case, defining this fix after reading the restart file
will enable the particles to keep their previous electronic energy
mode values.

**Related commands:**

:doc:`collide\_modify electronic discrete <collide_modify>`,
:doc:`react\_modify partial\_energy <react_modify>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
