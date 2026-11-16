.. index:: compute

compute command
===============

**Syntax:**


.. parsed-literal::

   compute ID style args

* ID = user-assigned name for the computation
* style = one of a list of possible style names (see below)
* args = arguments used by a particular style

**Examples:**


.. parsed-literal::

   compute 1 ke/particle 
   compute myGrid all n mass u usq temp

**Description:**

Define a computation that will be performed on a collection of
particles or grid cells or surface elements.  Quantities calculated by
a compute are instantaneous values, meaning they are calculated from
information about the current timestep.  Examples include calculation
of the system temperature or counting collisions of particles with
surface elements.  Code for new computes can be added to SPARTA; see
:doc:`Section 10 <Section_modify>` of the manual for details.

Note that defining a compute does not perform a computation.  Instead
computes are invoked by other SPARTA commands as needed, e.g. to
generate statistics or dump file output.  See `Section 6.4 <Section_howto.html#howto_4>`_ for a summary of various SPARTA output
options, many of which involve computes.

The ID for a compute is used to identify the compute in other
commands.  Each compute ID must be unique.  The ID can only contain
alphanumeric characters and underscores.  You can specify multiple
computees of the same style so long as they have different IDs.  A
compute can be deleted with the :doc:`uncompute <uncompute>` command,
after which its ID can be re-used.


----------


Each compute style has its own doc page which describes its arguments
and what it does.  Here is an alphabetic list of compute styles
available in SPARTA:

* :doc:`boundary <compute_boundary>` - various quantities on each global boundary
* :doc:`count <compute_count>` - particle counts for species and mixtures and mixture groups
* :doc:`distsurf/grid <compute_distsurf_grid>` - distance from grid cells to surface
* :doc:`dt/grid <compute_dt_grid>` - optimal timestep per grid cell
* :doc:`eflux/grid <compute_eflux_grid>` - energy flux density per grid cell
* :doc:`fft/grid <compute_fft_grid>` - FFTs across grid cells
* :doc:`gas/collision/grid <compute_gas_collision_grid>` - gas particle collisions per grid cell
* :doc:`gas/collision/tally <compute_gas_collision_tally>` - tallies for gas particle collisions
* :doc:`gas/reaction/grid <compute_gas_reaction_grid>` - gas particle reactions per grid cell
* :doc:`gas/reaction/tally <compute_gas_reaction_tally>` - tallies for gas particle reactions
* :doc:`grid <compute_grid>` - various per grid cell quantities
* :doc:`isurf/grid <compute_isurf_grid>` - various implicit surface element quantities
* :doc:`ke/particle <compute_ke_particle>` - temperature per particle
* :doc:`lambda/grid <compute_lambda_grid>` - mean-free path per grid cell
* :doc:`pflux/grid <compute_pflux_grid>` - momentum flux density per grid cell
* :doc:`property/grid <compute_property_grid>` - per grid cell properties
* :doc:`property/surf <compute_property_surf>` - per surface element properties
* :doc:`react/boundary <compute_react_boundary>` - reaction stats on global boundary
* :doc:`react/surf <compute_react_surf>` = reaction stats for explicit surfs
* :doc:`react/isurf/grid <compute_react_isurf_grid>` - reactions stats for implicit surfs
* :doc:`reduce <compute_reduce>` - reduce vectors to scalars
* :doc:`sonine/grid <compute_sonine_grid>` - Sonine moments per grid cell
* :doc:`surf <compute_surf>` - various explicit surface element quantities
* :doc:`surf/collision/tally <compute_surf_collision_tally>` - tallies for particle/surface collisions
* :doc:`surf/reaction/tally <compute_surf_reaction_tally>` - tallies for particle/surface reactions
* :doc:`thermal/grid <compute_thermal_grid>` - thermal temperature per grid cell
* :doc:`temp <compute_temp>` - temperature of particles
* :doc:`tvib/grid <compute_tvib_grid>` - vibrational temperature per grid cell

There are also additional accelerated compute styles included in the
SPARTA distribution for faster performance on specific hardware.  The
list of these with links to the individual styles are given in the
pair section of `this page <Section_commands.html#cmd_5>`_.


----------


Computes calculate one of five styles of quantities: global,
per-particle, per-grid, per-surf, or per-tally.  A global quantity is
one or more system-wide values, e.g. the temperature of the system.  A
per-particle quantity is one or more values per particle, e.g. the
kinetic energy of each particle.  A per-grid quantity is one or more
values per grid cell.  A per-surf quantity is one or more values per
surface element.  A per-tally quantity is one or more values per
event, e.g. a particle colliding or reacting with a surface element.

Global, per-particle, per-grid, per-surf, and per-tally quantities
each come in two forms: a single scalar value or a vector of values.
Additionaly, global quantities can also be a 2d array of values.  The
doc page for each compute describes the style and kind of values it
produces, e.g. a per-particle vector.  Some computes can produce more
than one form of a single style, e.g. a global scalar and a global
vector.

When a compute quantity is accessed, as in many of the output commands
discussed below, it can be referenced via the following bracket
notation, where ID is the ID of the compute:

+-------------+--------------------------------------------+
| c\_ID       | entire scalar, vector, or array            |
+-------------+--------------------------------------------+
| c\_ID[I]    | one element of vector, one column of array |
+-------------+--------------------------------------------+
| c\_ID[I][J] | one element of array                       |
+-------------+--------------------------------------------+

In other words, using one bracket reduces the dimension of the
quantity once (vector -> scalar, array -> vector).  Using two brackets
reduces the dimension twice (array -> scalar).  Thus a command that
uses scalar compute values as input can also process elements of a
vector or array.

Note that commands and :doc:`variables <variable>` which use compute
quantities typically do not allow for all kinds, e.g. a command may
require a vector of values, not a scalar.  This means there is no
ambiguity about referring to a compute quantity as f\_ID even if it
produces, for example, both a scalar and vector.  The doc pages for
various commands explain the details.


----------


The values generated by a compute can be used in several ways:

* Global values can be output via the :doc:`stats\_style <stats_style>`
  command.  Or the values can be referenced in a :doc:`variable equal <variable>` or :doc:`variable atom <variable>` command.
* Per-particle values can be output via the :doc:`dump particle <dump>`
  command.  Or the values can be referenced in a :doc:`particle-style variable <variable>`.
* Per-grid values can be output via the :doc:`dump grid <dump>` command.
  They can be time-averaged via the :doc:`fix ave/grid <fix_ave_grid>`
  command.
* Per-surf values can be output via the :doc:`dump surf <dump>` command.
  They can be time-averaged via the :doc:`fix ave/surf <fix_ave_surf>`
  command.
* Per-tally values can be output via the :doc:`dump tally <dump>`
  command.


----------


**Restrictions:** none

**Related commands:**

:doc:`uncompute <uncompute>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
