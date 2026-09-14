.. index:: fix field/grid

fix field/grid command
======================

**Syntax:**


.. parsed-literal::

   fix ID field/grid axvar ayvar azvar

* ID is documented in :doc:`fix <fix>` command
* field/grid = style name of this fix command
* axvar,ayvar,azvar = names of grid-style variables for acceleration components

**Examples:**


.. parsed-literal::

   fix 1 field/grid gradBx gradBy NULL

**Description:**

Specify the formulas used to calculate the acceleration effect of an
external field on particle motion.  The ID of this fix can be used by
the :doc:`global field grid <global>` command which applies the field
when particles are advected during a simulation run.  This is done by
invoking a method in this fix, which evaluates the specified
grid-style variables.

Each of the *axvar*\ , *ayvar*\ , and *azvar* arguments is the name of a
grid-style variable.  The variables should compute the x,y,z
components of acceleration applied at the center point of each grid
cell in the simulation.  Any of the three variables can be specified
as NULL, which means there is no acceleration in that dimension.

Each timestep when a particle is advected the acceleration vector (a)
for the grid cell it is in acts as a perturbation on straight-line
motion which affects both the end-of-timestep position (x) and
velocity (v) vectors of the particle:


.. parsed-literal::

   xnew = x + dt\*v + 0.5\*a\*dt\^2
   vnew = v + dt\*a

Note that the formulas encoded by the *axvar*\ , *ayvar*\ , and *azvar*
variables should produce values that are in units of acceleration
(distance/time\^2, see the :doc:`units <units>` command), not force.  And
they should not include the timestep (dt) value in the formulas above.
That is applied by SPARTA during advection.

See the :doc:`variable <variable>` doc page for a description of the
formula syntax allowed for grid-style variables.  They can include
references to the grid vectors xc, yc, and zc for the grid cell center
point.  Using these values in a formula can enable a
spatially-dependent field.  The formulas can also include the current
timestep and timestep size (dt) to enable a time-dependent field.

NOTE still need to figure this out: And they can include properties of
the particle, such as its mass or magnetic moment.


----------


Note that the :doc:`global field <global>` command provides three
alternatives for specifying an external field:


.. parsed-literal::

   global field constant ...     # field is constant in space and time
   global field particle ...     # field is applied on a per particle basis
   global field grid ...         # field is applied on a per grid cell basis

This fix is only used for per-grid fields.  It should only be used for
fields which vary spatially or in time; otherwise use the constant
option which will be much more efficient.  The use of per-grid
variables allows the field to vary spatially as a function of the grid
cell center point. It also allows the field to vary in time by having
the variables use the current timestep.

NOTE still need to figure out how to do this: The field can also
depend on particle attributes, such as its mass and magnetic moment
(for a B field).

Note that use of the :doc:`global field grid <global>` command with this
fix will evaluate the specified grid-style variables as often as
requested.  For a field that has no time-ependence, you can specificy
it only be evaluated once at the beginning of a run.  For a field that
is time-dependent you can choose how often to recompute the field,
depending on how fast it varies.

The :doc:`fix field/particle <fix_field_particle>` command is an
alternative which will typically run much slower, but be more
accurate.  When used with the :doc:`global field particle <global>`
command, the particle-style variables it uses are invoked every
timestep using current particle positions.  And the field calculation
is done for each grid particle, not for each grid cell.  The trade-off
is that the fields it calculates for each particle is more accurate,
but the simulation will typically run several times slower than it
would with this fix.

**Restart, output info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

This fix stores a per-grid cell array of values which can be accessed
by various output commands on any timestep, e.g. by the :doc:`dump grid <dump>` command.  The values are those produced by evaluating
the grid-style variables.  The number of rows in the array is the
number of grid cells this processor owns.  The number of columns in
the array is the number of non-NULL variables specified.

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

:doc:`fix field/particle <fix_field_particle>`, :doc:`global field <global>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
