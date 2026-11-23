.. index:: stats

stats command
=============

**Syntax:**


.. parsed-literal::

   stats N

* N = output statistics every N timesteps

**Examples:**


.. parsed-literal::

   stats 100

**Description:**

Compute and print statistical info (e.g. particle count, temperature)
on timesteps that are a multiple of N and at the beginning and end of
a simulation run.  A value of 0 will only print statistics at the
beginning and end.

The content and format of what is printed is controlled by the
:doc:`stats\_style <stats_style>` and :doc:`stats\_modify <stats_modify>`
commands.

The timesteps on which statistical output is written can also be
controlled by a :doc:`variable <variable>`.  See the :doc:`stats\_modify every <stats_modify>` command.

**Restrictions:** none

**Related commands:**

:doc:`stats\_style <stats_style>`, :doc:`stats\_modify <stats_modify>`

**Default:**


.. parsed-literal::

   stats 0


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
