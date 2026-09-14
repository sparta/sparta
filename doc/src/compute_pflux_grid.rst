.. index:: compute pflux/grid

compute pflux/grid command
==========================

compute pflux/grid/kk command
=============================

**Syntax:**


.. parsed-literal::

   compute ID pflux/grid group-ID mix-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* pflux/grid = style name of this compute command
* group-ID = group ID for which grid cells to perform calculation on
* mix-ID = mixture ID to perform calculation on
* one or more values can be appended
* values = *momxx* or *momyy* or *momzz* or *momxy* or *momyz* or *momxz*
  
  .. parsed-literal::
  
       momxx\ ,\ momyy\ ,\ momzz = diagonal components of momentum flux density tensor
       momxy\ ,\ momyz\ ,\ momxz = off-diagonal components of momentum flux density tensor



**Examples:**


.. parsed-literal::

   compute 1 pflux/grid all species momxx momyy momzz
   compute 1 pflux/grid subset species momxx momxy

These commands will dump time averaged momentum flux densities for
each species and each grid cell to a dump file every 1000 steps:


.. parsed-literal::

   compute 1 pflux/grid all species momxx momyy momzz
   fix 1 ave/grid 10 100 1000 c_1[\*]
   dump 1 grid all 1000 tmp.grid id f_1[\*]

**Description:**

Define a computation that calculates components of the momemtum flux
density tensor for each grid cell in a grid cell group.  This is
equivalent to the kinetic energy density tensor, and is based on the
thermal velocity of the particles in each grid cell.  The values are
tallied separately for each group of species in the specified mixture,
as described in the Output section below.  See the mixture command for
how a set of species can be partitioned into groups.

Only grid cells in the grid group specified by *group-ID* are included
in the calculations.  See the :doc:`group grid <group>` command for info
on how grid cells can be assigned to grid groups.

The values listed above rely on first computing and subtracting the
center-of-mass (COM) velocity for all particles in the group and grid
cell from each particle to yield a thermal velocity.  This thermal
velocity is used to compute the components of the momentum flux
density tensor, as described below.  This is in contrast to some of
the values tallied by the :doc:`compute grid temp <compute_grid>`
command which simply uses the full velocity of each particle to
compute a momentum or kinetic energy density.  For non-streaming
simulations, the two results should be similar, but for streaming
flows, they will be different.

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


Calculation of the momentum flux density is done by first calcuating
the center-of-mass (COM) velocity of particles for each group within a
grid cell.  This is done as follows:


.. parsed-literal::

   COMx = Sum_i (mass_i Vx_i) / Sum_i (mass_i)
   COMy = Sum_i (mass_i Vy_i) / Sum_i (mass_i)
   COMz = Sum_i (mass_i Vz_i) / Sum_i (mass_i)
   Cx = Vx - COMx
   Cy = Vy - COMy
   Cz = Vz - COMz

The COM velocity is (COMx,COMy,COMz).  The thermal velocity of each
particle is (Cx,Cy,Cz), i.e. its velocity minus the COM velocity of
particles in its group and cell.

The *momxx*\ , *momyy*\ , *momzz* values compute the diagonal components
of the momentum flux density tensor due to particles in the group as
follows:


.. parsed-literal::

   momxx = fnum/volume Sum_i (mass_i Cx\^2)
   momyy = fnum/volume Sum_i (mass_i Cy\^2)
   momzz = fnum/volume Sum_i (mass_i Cz\^2)

The *momxy*\ , *momyz*\ , *momxz* values compute the off-diagonal
components of the momentum flux density tensor due to particles in the
group as follows:


.. parsed-literal::

   momxy = fnum/volume Sum_i (mass_i Cx Cy)
   momyz = fnum/volume Sum_i (mass_i Cy Cz)
   momxz = fnum/volume Sum_i (mass_i Cx Cz)

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
momentum flux density = energy density = energy/volume units.


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

:doc:`compute grid <compute_grid>`, :doc:`compute thermal/grid <compute_thermal_grid>`, :doc:`compute eflux/grid <compute_eflux_grid>`, :doc:`fix ave/grid <fix_ave_grid>`,
:doc:`dump grid <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
