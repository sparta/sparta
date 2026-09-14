.. index:: quit

quit command
============

**Syntax:**


.. parsed-literal::

   quit

**Examples:**


.. parsed-literal::

   quit
   if "$n > 10000" then quit

**Description:**

This command causes SPARTA to exit, after shutting down all
output cleanly.

It can be used as a debug statement in an input script, to terminate
the script at some intermediate point.

It can also be used as an invoked command inside the
"then" or "else" portion of an :doc:`if <if>` command.

**Restrictions:** none

**Related commands:**

:doc:`if <if>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
