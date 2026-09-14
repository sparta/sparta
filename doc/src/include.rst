.. index:: include

include command
===============

**Syntax:**


.. parsed-literal::

   include file

* file = filename of new input script to switch to

**Examples:**


.. parsed-literal::

   include newfile
   include in.run2

**Description:**

This command opens a new input script file and begins reading SPARTA
commands from that file.  When the new file is finished, the original
file is returned to.  Include files can be nested as deeply as
desired.  If input script A includes script B, and B includes A, then
SPARTA could run for a long time.

If the filename is a variable (see the :doc:`variable <variable>`
command), different processor partitions can run different input
scripts.

**Restrictions:** none

**Related commands:**

:doc:`variable <variable>`, :doc:`jump <jump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
