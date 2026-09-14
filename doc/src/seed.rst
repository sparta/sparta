.. index:: seed

seed command
============

**Syntax:**


.. parsed-literal::

   seed Nvalue

* Nvalue = seed for a random number generator (positive integer)

**Examples:**


.. parsed-literal::

   seed 5838959

**Description:**

This command sets the random number seed for a master random number
generator.  This generator is used by SPARTA to initialize auxiliary
random number generators, which in turn are used for all operations in
the code requiring random numbers.  This means you can effectively run
a statistically-independent simulation by simply changing this single
seed.

The various random number generators used in SPARTA are portable,
which means they produce the same random number streams on any
machine.

This command is required to perform a SPARTA simulation.

**Restrictions:** none

**Related commands:** none

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
