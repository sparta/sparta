.. index:: fix

fix command
===========

**Syntax:**


.. parsed-literal::

   fix ID style args

* ID = user-assigned name for the fix
* style = one of a long list of possible style names (see below)
* args = arguments used by a particular style

**Examples:**


.. parsed-literal::

   fix 1 grid/check 100 warn
   fix 1 ave/time all 100 5 1000 c_myTemp c_thermo_temp file temp.profile

**Description:**

Set a fix that will be applied to the system.  In SPARTA, a "fix" is
an operation that is applied to the system during timestepping.
Examples include adding particles via inlet boundary conditions or
computing diagnostics.  Code for new fixes can be added to SPARTA; see
:doc:`Section 10 <Section_modify>` of the manual for details.

Fixes perform their operations at different stages of the timestep.
If 2 or more fixes operate at the same stage of the timestep, they are
invoked in the order they were specified in the input script.

The ID for a fix is used to identify the fix in other commands.  Each
fix ID must be unique; see an exception below.  The ID can only
contain alphanumeric characters and underscores.  You can specify
multiple fixes of the same style so long as they have different IDs.
A fix can be deleted with the :doc:`unfix <unfix>` command, after which
its ID can be re-used.

IMPORTANT NOTE: The :doc:`unfix <unfix>` command is the only way to turn
off a fix; simply specifying a new fix with the same style and a
different ID will not turn off the first one.

If you specify a new fix with the same ID and style as an existing
fix, the old fix is deleted and the new one is created (presumably
with new settings).  This is the same as if an "unfix" command were
first performed on the old fix, except that the new fix is kept in the
same order relative to the existing fixes as the old one originally
was.

Some fixes store an internal "state" which is written to binary
restart files via the :doc:`restart <restart>` or
:doc:`write\_restart <write_restart>` commands.  This allows the fix to
continue on with its calculations in a restarted simulation.  See the
:doc:`read\_restart <read_restart>` command for info on how to re-specify
a fix in an input script that reads a restart file.  See the doc pages
for individual fixes for info on which ones can be restarted.


----------


Each fix style has its own doc page which describes its arguments and
what it does, as listed below.  Here is an alphabetic list of fix
styles available in SPARTA:

* :doc:`ablate <fix_ablate>` - alter implicit surfaces within each grid cell
* :doc:`adapt <fix_adapt>` - on-the-fly grid adaptation
* :doc:`ambipolar <fix_ambipolar>` - ambipolar approximation for ionized plasmas
* :doc:`ave/grid <fix_ave_grid>` - compute per grid cell time-averaged quantities
* :doc:`ave/histo <fix_ave_histo>` - compute/output time averaged histograms
* :doc:`ave/histo/weight <fix_ave_histo>` - compute/output weighted histograms
* :doc:`ave/surf <fix_ave_surf>` - compute per surface element time-averaged quantities
* :doc:`ave/time <fix_ave_time>` - compute/output global time-averaged quantities
* :doc:`balance <fix_balance>` - perform dynamic load-balancing
* :doc:`custom <fix_custom>` - periodically reset the values of one or more custom attributes
* :doc:`dt/reset <fix_dt_reset>` - adjust global timestep dynamically
* :doc:`emit/face <fix_emit_face>` - emit particles at global boundaries
* :doc:`emit/face/file <fix_emit_face_file>` - emit particles at global boundaries using a distribution defined in a file
* :doc:`emit/surf <fix_emit_surf>` - emit particles at surfaces
* :doc:`field/grid <fix_field_grid>` - apply an external field on a per grid cell basis
* :doc:`field/particle <fix_field_particle>` - apply an external field on a per particle basis
* :doc:`grid/check <fix_grid_check>` - check if particles are in the correct grid cell
* :doc:`halt <fix_halt>` - stop the simulation early if a condition is met
* :doc:`move/surf <fix_move_surf>` - move surfaces dynamically during a simulation
* :doc:`print <fix_print>` - print text and variables during a simulation
* :doc:`surf/temp <fix_surf_temp>` - compute per-surf temperatures dynamically
* :doc:`temp/global/rescale <fix_temp_global_rescale>` - rescale particle temperatures
* :doc:`temp/rescale <fix_temp_rescale>` - rescale particle temperatures within each grid cell
* :doc:`surf/temp <fix_surf_temp>` - compute per-surf temperatures dynamically
* :doc:`vibmode <fix_vibmode>` - discrete vibrational energy modes

There are also additional accelerated compute styles included in the
SPARTA distribution for faster performance on specific hardware.  The
list of these with links to the individual styles are given in the
pair section of `this page <Section_commands.html#cmd_5>`_.


----------


In addition to the operation they perform, some fixes also produce one
of four styles of quantities: global, per-particle, per-grid, or
per-surf.  These can be used by other commands or output as described
below.  A global quantity is one or more system-wide values, e.g. the
temperature of the system.  A per-particle quantity is one or more
values per particle, e.g. the kinetic energy of each particle.  A
per-grid quantity is one or more values per grid cell.  A per-surf
quantity is one or more values per surface element.

Global, per-particle, per-grid, and per-surf quantities each come in
two forms: a single scalar value or a vector of values.  Additionaly,
global quantities can also be a 2d array of values.  The doc page for
each fix describes the style and kind of values it produces, e.g. a
per-particle vector.  Some fixes can produce more than one form of a
single style, e.g. a global scalar and a global vector.

When a fix quantity is accessed, as in many of the output commands
discussed below, it can be referenced via the following bracket
notation, where ID is the ID of the fix:

+-------------+--------------------------------------------+
| f\_ID       | entire scalar, vector, or array            |
+-------------+--------------------------------------------+
| f\_ID[I]    | one element of vector, one column of array |
+-------------+--------------------------------------------+
| f\_ID[I][J] | one element of array                       |
+-------------+--------------------------------------------+

In other words, using one bracket reduces the dimension of the
quantity once (vector -> scalar, array -> vector).  Using two brackets
reduces the dimension twice (array -> scalar).  Thus a command that
uses scalar fix values as input can also process elements of a vector
or array.

Note that commands and :doc:`variables <variable>` which use fix
quantities typically do not allow for all kinds, e.g. a command may
require a vector of values, not a scalar.  This means there is no
ambiguity about referring to a fix quantity as f\_ID even if it
produces, for example, both a scalar and vector.  The doc pages for
various commands explain the details.


----------


Any values generated by a fix can be used in several ways:

* Global values can be output via the :doc:`stats\_style <stats_style>`
  command.  Or the values can be referenced in a :doc:`variable equal <variable>` or :doc:`variable atom <variable>` command.
* Per-particle values can be output via the :doc:`dump particle <dump>`
  command.  Or the per-particle values can be referenced in an
  :doc:`particle-style variable <variable>`.
* Per-grid values can be output via the :doc:`dump grid <dump>` command.
  Or the per-grid values can be referenced in a :doc:`grid-style variable <variable>`.


----------


**Restrictions:** none

**Related commands:**

:doc:`unfix <unfix>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
