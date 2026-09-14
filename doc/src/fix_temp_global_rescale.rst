.. index:: fix temp/global/rescale

fix temp/global/rescale command
===============================

fix temp/global/rescale/kk command
==================================

**Syntax:**


.. parsed-literal::

   fix ID temp/global/rescale N Tstart Tstop fraction

* ID is documented in :doc:`fix <fix>` command
* temp/global/rescale = style name of this fix command
* N = thermostat every N timesteps
* Tstart,Tstop = desired temperature at start/end of run (temperature units)
* fraction = rescale to target temperature by this fraction

**Examples:**


.. parsed-literal::

   fix 1 temp/global/rescale 100 300.0 300.0 0.5
   fix 5 temp/global/rescale 10 300.0 10.0 1.0

**Description:**

Reset the temperature of all the particles in the entire simulation by
explicitly rescaling their velocities.  This is a simple
thermostatting operation to keep the temperature of the gas near the
desired target temperature.  This can be useful if an external driving
force is adding energy to the system.  Or if you wish the heat or cool
the temperature of the system over time.

The rescaling is applied to only the translational degrees of freedom
for the particles.  Their rotational or vibrational degrees of freedom
are not altered.

Rescaling is performed every N timesteps. The target temperature is a
ramped value between the Tstart and Tstop temperatures at the
beginning and end of the run.

From the current global temperature and the current target
temperature, a velocity scale factor is calculated.  The amount of
rescaling that is applied is adjusted by the *fraction* parameter
which is a value from 0.0 to 1.0.  difference between the actual and
desired temperature.  If *fraction* = 1.0, the temperature is reset to
exactly the desired value.  If *fraction* = 0.5, the temperature is
reset to a value halfway between the current global and target
temperatures.

The rescaling factor is applied to each of the components of the
translational velocity for every particle in the simulation.

Note that this fix performs thermostatting using the same formula for
temperature as calculated by the :doc:`compute temp <compute_temp>`
command.  It does not currently subtract out a net streaming velocity
to measure a thermal temperature since it assumes the net center of
mass velocity for the entire system is zero.  An option for this may
be added in the future.  See the :doc:`fix temp/rescale <fix_temp_rescale>` doc page for a command that
thermostats the thermal temperature on a per-grid-cell basis.


----------


**Restart, output info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

This fix produces no output.

This fix can ramp its target temperature over multiple runs, using the
start and stop keywords of the run command. See the run command for
details of how to do this.


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

:doc:`fix temp/rescale <fix_temp_rescale>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
