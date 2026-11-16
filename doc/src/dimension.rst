.. index:: dimension

dimension command
=================

**Syntax:**


.. parsed-literal::

   dimension N

* N = 2 or 3

**Examples:**


.. parsed-literal::

   dimension 2
   dimension 3

**Description:**

Set the dimensionality of the simulation.  By default SPARTA runs 3d
simulations, but 2d simulations can also be run.

2d axi-symmetric models can be run by setting the dimension to 2, and
defining the lower boundary in the y-dimension to axi-symmetric via
the :doc:`boundary <boundary>` command.

**Restrictions:**

This command must be used before the simulation box is defined by a
:doc:`create\_box <create_box>` command.

**Related commands:** none

**Default:**


.. parsed-literal::

   dimension 3


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
