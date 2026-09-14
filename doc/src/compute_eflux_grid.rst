.. index:: compute eflux/grid

compute eflux/grid command
==========================

compute eflux/grid/kk command
=============================

**Syntax:**


.. parsed-literal::

   compute ID eflux/grid group-ID mix-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* eflux/grid = style name of this compute command
* group-ID = group ID for which grid cells to perform calculation on
* mix-ID = mixture ID to perform calculation on
* one or more values can be appended
* values = *heatx* or *heaty* or *heatz*
  
  .. parsed-literal::
  
       heatx\ ,\ heaty\ ,\ heatz = xyz components of energy flux density tensor



**Examples:**


.. parsed-literal::

   compute 1 eflux/grid all species heatx heaty heatz
   compute 1 eflux/grid subset species heaty

These commands will dump time averaged energy flux densities for
each species and each grid cell to a dump file every 1000 steps:


.. parsed-literal::

   compute 1 eflux/grid all species heatx heaty heatz
   fix 1 ave/grid 10 100 1000 c_1[\*]
   dump 1 grid all 1000 tmp.grid id f_1[\*]

**Description:**

Define a computation that calculates components of the energy flux
density vector for each grid cell in a grid cell group.  This is also
called the heat flux density vector, and is based on the thermal
velocity of the particles in each grid cell.  The values are tallied
separately for each group of species in the specified mixture, as
described in the Output section below.  See the mixture command for
how a set of species can be partitioned into groups.

Only grid cells in the grid group specified by *group-ID* are included
in the calculations.  See the :doc:`group grid <group>` command for info
on how grid cells can be assigned to grid groups.

The values listed above rely on first computing and subtracting the
center-of-mass (COM) velocity for all particles in the group and grid
cell from each particle to yield a thermal velocity.  This thermal
velocity is used to compute the components of the energy flux density
vector, as described below.  This is in contrast to some of the values
tallied by the :doc:`compute grid temp <compute_grid>` command which
simply uses the full velocity of each particle to compute a momentum
or kinetic energy density.  For non-streaming simulations, the two
results should be similar, but for streaming flows, they will be
different.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump grid <dump>` command.

The values over many sampling timesteps can be averaged by the :doc:`fix ave/grid <fix_ave_grid>` command.  It does its averaging as if the
particles in the cell at each sampling timestep were combined together
into one large set of particles to compute the formulas below.

Note that the center-of-mass (COM) velocity that is subtracted from
each particle to yield a thermal velocity for each particle, as
described below, is also computed over one large set of particles
(across all timesteps), in contrast to using a COM velocity computed
only for particles in the current timestep, which is what the :doc:`compute sonine/grid <compute_sonine_grid>` command does.

Note that this is a different form of averaging than taking the values
produced by the formulas below for a single timestep, summing those
values over the sampling timesteps, and then dividing by the number of
sampling steps.


----------


Calculation of the energy flux density is done by first calcuating the
center-of-mass (COM) velocity of particles for each group with a grid
cell.  This is done as follows:


.. parsed-literal::

   COMx = Sum_i (mass_i Vx_i) / Sum_i (mass_i)
   COMy = Sum_i (mass_i Vy_i) / Sum_i (mass_i)
   COMz = Sum_i (mass_i Vz_i) / Sum_i (mass_i)
   Cx = Vx - COMx
   Cy = Vy - COMy
   Cz = Vz - COMz
   Csq = Cx\*Cx + Cy\*Cy + Cz\*Cz

The COM velocity is (COMx,COMy,COMz).  The thermal velocity of each
particle is (Cx,Cy,Cz), i.e. its velocity minus the COM velocity of
particles in its group and cell.

The *heatx*\ , *heaty*\ , *heatz* values compute the components of the
energy flux density vector due to particles in the group as follows:


.. parsed-literal::

   heatx = 0.5 \* fnum/volume Sum_i (mass_i Cx Csq)
   heaty = 0.5 \* fnum/volume Sum_i (mass_i Cy Csq)
   heatz = 0.5 \* fnum/volume Sum_i (mass_i Cz Csq)

Note that if particle weighting is enabled via the :doc:`global weight <global>` command, then the volume used in the formula is
divided by the weight assigned to the grid cell.


----------


**Output info:**

This compute calculates a per-grid array, with the number of columns
equal to the number of values times the number of groups.  The
ordering of columns is first by values, then by groups.  I.e. if
*momxx* and *momxy* values were specified as keywords, then the first
two columns would be *momxx* and *momxy* for the first group, the 3rd
and 4th columns would be *momxx* and *momxy* for the second group, etc.

This compute performs calculations for all flavors of child grid cells
in the simulation, which includes unsplit, cut, split, and sub cells.
See `Section 6.8 <Section_howto.html#howto_8>`_ of the manual gives
details of how SPARTA defines child, unsplit, split, and sub cells.
Note that cells inside closed surfaces contain no particles.  These
could be unsplit or cut cells (if they have zero flow volume).  Both
of these kinds of cells will compute a zero result for all their
values.  Likewise, split cells store no particles and will produce a
zero result.  This is because their sub-cells actually contain the
particles that are geometrically inside the split cell.

Grid cells not in the specified *group-ID* will output zeroes for all
their values.

The array can be accessed by any command that uses per-grid values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The per-grid array values will be in the :doc:`units <units>` of
energy flux density = energy-velocity/volume units.


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

:doc:`compute grid <compute_grid>`, :doc:`compute thermal/grid <compute_thermal_grid>`, :doc:`compute pflux/grid <compute_pflux_grid>`, :doc:`fix ave/grid <fix_ave_grid>`,
:doc:`dump grid <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
