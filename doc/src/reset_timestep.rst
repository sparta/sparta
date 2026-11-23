.. index:: reset\_timestep

reset\_timestep command
=======================

**Syntax:**


.. parsed-literal::

   reset_timestep N

* N = timestep number

**Examples:**


.. parsed-literal::

   reset_timestep 0
   reset_timestep 4000000

**Description:**

Set the timestep counter to the specified value.  This command
normally comes after the timestep has been set by reading a restart
file via the :doc:`read\_restart <read_restart>` command, or a previous
simulation advanced the timestep.

The :doc:`create\_box <create_box>` command sets the timestep to 0; the
:doc:`read\_restart <read_restart>` command sets the timestep to the
value it had when the restart file was written.

**Restrictions:** none

This command cannot be used when any fixes are defined that keep track
of elapsed time to perform certain kinds of time-dependent operations.
Examples are the :doc:`fix ave/time <fix_ave_time>`, :doc:`fix ave/grid <fix_ave_grid>`, and :doc:`fix ave/surf <fix_ave_surf>`
commands.  Thus these fixes should be specified after the timestep has
been reset.

Resetting the timestep clears flags for :doc:`computes <compute>` that
may have calculated some quantity from a previous run.  This means
these quantity cannot be accessed by a variable in between runs until
a new run is performed.  See the :doc:`variable <variable>` command for
more details.

**Related commands:** none

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
