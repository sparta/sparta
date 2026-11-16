.. index:: units

units command
=============

**Syntax:**


.. parsed-literal::

   units style

* style = *cgs* or *si*

**Examples:**


.. parsed-literal::

   units cgs

**Description:**

This command sets the style of units used for a simulation.  It
determines the units of all quantities specified in the input script
and various input files read by SPARTA, as well as the units of all
quantities output to the screen, log file, dump files, and other
output files.  Typically, this command is used at the very beginning
of an input script.

IMPORTANT NOTE: Internally, this command simply sets the numeric
values of conversion factors used by SPARTA, e.g. the Boltzmann
constant used to convert temperature to energy.  It is up to you to
insure that all input values used in the input script and other input
files (surface data, species files, reaction files) contain numeric
values consistent with the chosen units.

For style *cgs*\ , these are the units:

* mass = grams
* distance = centimeters
* area = cm\^2
* volume = cm\^3
* time = seconds
* energy = ergs
* velocity = centimeters/second
* acceleration = centimeters/second\^2
* pressure = barye (dyne/cm\^2 = 0.1 pascals)
* magnetic moment = ??
* temperature = degrees K

For style *si*\ , these are the units:

* mass = kilograms
* distance = meters
* area = m\^2
* volume = m\^3
* time = seconds
* energy = Joules
* velocity = meters/second
* acceleration = meters/second\^2
* pressure = pascals (newton/meter\^2)
* magnetic moment = ??
* temperature = degrees K

The units command also sets a default timestep size; see the
:doc:`timestep <timestep>` command to change this value.

* For style *cgs* this is dt = 1.0 sec.
* For style *si* this is dt = 1.0 sec.

**Restrictions:**

This command must be used before the simulation box is defined by a
:doc:`create\_box <create_box>` command.

**Related commands:** none

**Default:**


.. parsed-literal::

   units si


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
