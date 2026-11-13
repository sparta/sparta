.. index:: compute isurf/grid

compute isurf/grid command
==========================

compute isurf/grid/kk command
=============================

**Syntax:**


.. parsed-literal::

   compute ID isurf/grid group-ID mix-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* isurf/grid = style name of this compute command
* group-ID = group ID for which grid cells to perform calculation on
* mix-ID = mixture ID for particles to perform calculation on
* one or more values can be appended
* value = *n* or *nwt* or *mflux* or *fx* or *fy* or *fz* or *press* or *px* or *py* or *pz* or *shx* or *shy* or *shz* or *ke*
  
  .. parsed-literal::
  
       n = count of particles hitting surface elements in a grid cell
       nwt = weighted count of particles hitting surface elements in a grid cell
       mflux = flux of mass on surface elements in a grid cell
       fx,fy,fz = components of force on surface elements in a grid cell
       press = magnitude of normal pressure on surface elements in a grid cell
       px,py,pz = components of normal pressure on surface elements in a grid cell
       shx,shy,shz = components of shear stress on surface elements in a grid cell
       ke = flux of particle kinetic energy on surface elements in a grid cell
       erot = flux of particle rotational energy on surface elements in a grid cell
       evib = flux of particle vibrational energy on surface elements in a grid cell
       etot = flux of particle total energy on surface elements in a grid cell



**Examples:**


.. parsed-literal::

   compute 1 isurf/grid all all n press eng
   compute mine isurf/grid sphere species press shx shy shz

These commands will dump time averages for each species and each
grid cell to a dump file every 1000 steps:


.. parsed-literal::

   compute 1 isurfgrid all species n press shx shy shz
   fix 1 ave/grid all 10 100 1000 c_1[\*]
   dump 1 grid all 1000 tmp.grid id f_1[\*]

These commands will time-average the force surface elements in each
grid cell, then sum them across grid cells to compute drag (fx) and
lift (fy) on the set of implicit surfs:


.. parsed-literal::

   compute 1 isurf/grid all all fx fy
   fix 1 ave/grid all 10 100 1000 c_1[\*]
   compute 2 reduce sum f_1[1] f_1[2]
   stats 1000
   stats_style step cpu np c_2[1] c_2[2]

**Description:**

Define a computation that calculates one or more values for each grid
cell in a grid cell group, based on the particles that collide with
the implicit surfaces in that grid cell.  The values are summed for
each group of species in the specified mixture.  See the
:doc:`mixture <mixture>` command for how a set of species can be
partitioned into groups.  Only grid cells in the grid group specified
by *group-ID* are included in the calculations.  See the :doc:`group grid <group>` command for info on how grid cells can be assigned to
grid groups.

Implicit surface elements are triangles for 3d simulations and line
segments for 2d simulations.  Unlike explicit surface elements, each
triangle or line segment is wholly contained within a single grid
cell.  See the :doc:`read\_isurf <read_isurf>` command for details.

This command can only be used for simulations with implicit surface
elements.  See the similar :doc:`compute surf <compute_surf>` command
for use with simulations with explicit surface elements.

Note that when a particle collides with a surface element, it can
bounce off (possibly as a different species), be captured by the
surface (vanish), or a 2nd particle can also be emitted.  The formulas
below account for all the possible outcomes.  For example, the kinetic
energy flux *ke* onto a suface element for a single collision includes
a positive contribution from the incoming particle and negative
contributions from 0, 1, or 2 outgoing particles.  The exception is
the *n* and *nwt* values which simply tally counts of particles
colliding with the surface element.

Also note that all values for a collision are tallied based on the
species group of the incident particle.  Quantities associated with
outgoing particles are part of the same tally, even if they are in
different species groups.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump grid <dump>` command.

The values over many sampling timesteps can be averaged by the :doc:`fix ave/grid <fix_ave_grid>` command.  It does its averaging as if the
particles striking the surface elements within the grid cell at each
sampling timestep were combined together into one large set to compute
the formulas below.  The answer is then divided by the number of
sampling timesteps if it is not otherwise normalized by the number of
particles.  Note that in general this is a different normalization
than taking the values produced by the formulas below for a single
timestep, summing them over the sampling timesteps, and then dividing
by the number of sampling steps.  However for the current values
listed below, the two normalization methods are the same.

NOTE: If particle weighting is enabled via the :doc:`global weight <global>` command, then all of the values below are scaled
by the weight assigned to the grid cell in which the particle
collision with the surface element occurs.  The only exception is the
the *n* value, which is NOT scaled by the weight; it is a simple count
of particle collisions with surface elements in the grid cell.


----------


The meaning of all the value keywords and the formulas for calculating
these quantities is exactly the same as described by the :doc:`compute surf <compute_surf>` command.

The only difference is that the quantities are calculated on a per
grid cell basis, summing over all the surface elements in that grid
cell.


----------


**Output info:**

This compute calculates a per-grid array, with the number of columns
equal to the number of values times the number of groups.  The
ordering of columns is first by values, then by groups.  I.e. if the
*n* and *u* values were specified as keywords, then the first two
columns would be *n* and *u* for the first group, the 3rd and 4th
columns would be *n* and *u* for the second group, etc.

Grid cells not in the specified *group-ID* will output zeroes for all
their values.

The array can be accessed by any command that uses per-grid values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The per-grid array values will be in the :doc:`units <units>`
appropriate to the individual values as described above. *N* is
unitless.  *Press*\ , *px*\ , *py*\ , *pz*\ , *shx*\ , *shy*\ , *shz* are in in
pressure units.  *Ke*\ , *erot*\ , *evib*\ , and *etot* are in
energy/area-time units for 3d simulations and energy/length-time units
for 2d simulations.


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

:doc:`fix ave/grid <fix_ave_grid>`, :doc:`dump grid <dump>`, :doc:`compute surf <compute_surf>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
