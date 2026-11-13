.. index:: fix temp/rescale

fix temp/rescale command
========================

fix temp/rescale/kk command
===========================

**Syntax:**


.. parsed-literal::

   fix ID temp/rescale N Tstart Tstop keyword value ...

* ID is documented in :doc:`fix <fix>` command
* temp/rescale = style name of this fix command
* N = thermostat every N timesteps
* Tstart,Tstop = desired temperature at start/end of run (temperature units)
* zero or more keyword/args pairs may be appended
* keyword = *ave*
  
  .. parsed-literal::
  
       ave values = yes or no



**Examples:**


.. parsed-literal::

   fix 1 temp/rescale 100 300.0 300.0
   fix 5 temp/rescale 10 300.0 10.0 ave yes

**Description:**

Reset the thermal temperature of all the particles within each grid
cell by explicitly rescaling their thermal velocities.  This is a
simple thermostatting operation to keep the thermal temperature of the
gas near the desired target temperature. This can be useful if an
external driving force is adding energy to the system.  Or if you wish
the thermal temperature of the system to heat or cool over time.

The rescaling is applied to only the translational degrees of freedom
for the particles.  Their rotational or vibrational degrees of freedom
are not altered.

Rescaling is performed every N timesteps. The target temperature
(Ttarget) is a ramped value between the Tstart and Tstop temperatures
at the beginning and end of the run.

This fix performs thermostatting on a per grid cell basis.  For each
grid cell, the center-of-mass velocity and thermal temperature of the
particles in the grid cell is computed.  The thermal temperature is
defined as the kinetic temperature after any center-of-mass motion
(e.g. a flow velocity) is subtracted from the collection of particles
in the grid cell.  These are the same calculations as are performed by
the :doc:`compute thermal/grid <compute_thermal_grid>` command.  See its
doc page for the equations.  See the :doc:`fix temp/global/rescale <fix_temp_global_rescale>` doc page for a
command that thermostats the temperature of the global system.

How the rescaling of particle velocities is done depends on the value
of the *ave* keyword.

For *ave* with a value *no* (the default), the thermal temperature
(Tthermal) of each cell is used to compute a velocity scale factor for
that cell, which is


.. parsed-literal::

   vscale = sqrt(Ttarget/Tthermal)

The vscale factor is applied to each of the components of the thermal
velocity for each particle in the grid cell.  Only cells with 2 or
more particles have their particle velocities rescaled.

For *ave* with a value *yes*\ , the thermal temperatures of all the
cells are averaged.  The average thermal temperature is simply the sum
of cell thermal temperatures divided by the number of cells.  Cells
with less than 2 particles or whose thermal temperature = 0.0
contribute a thermal temperaure = Ttarget to the average.  The average
thermal temperature (Tthermal\_ave) for all cells is used to compute a
velocity scale factor for all cells, which is


.. parsed-literal::

   vscale = sqrt(Ttarget/Tthermal_ave)

This single vscale factor is applied to each of the components of the
thermal velocity for each particle in all the grid cells, including
the particles in single-particle cells.

After rescaling, for either *ave* = *no* or *yes*\ , if the thermal
temperature were re-computed for any grid cell with more than one
particle, it would be exactly the target temperature.


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

:doc:`fix temp/global/rescale <fix_temp_global_rescale>`

**Default:**

The default is ave = no.


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
