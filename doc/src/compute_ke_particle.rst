.. index:: compute ke/particle

compute ke/particle command
===========================

compute ke/particle/kk command
==============================

**Syntax:**


.. parsed-literal::

   compute ID ke/particle

* ID is documented in :doc:`compute <compute>` command
* ke/particle = style name of this compute command

**Examples:**


.. parsed-literal::

   compute 1 ke/particle

**Description:**

Define a computation that calculates the per-atom translational
kinetic energy for each particle.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump particle <dump>` command.

The kinetic energy is


.. parsed-literal::

   Vsq = Vx\*Vx + Vy\*Vy + Vz\*Vz
   KE = 1/2 m Vsq

where m is the mass and (Vx,Vy,Vz) are the velocity components of the
particle.

**Output info:**

This compute calculates a per-particle vector, which can be accessed
by any command that uses per-particle values from a compute as input.

The vector can be accessed by any command that uses per-particle
values from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_ for an overview of SPARTA output
options.

The per-particle vector values will be in energy :doc:`units <units>`.


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

:doc:`dump particle <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
