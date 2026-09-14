.. index:: compute grid

compute grid command
====================

compute grid/kk command
=======================

**Syntax:**


.. parsed-literal::

   compute ID grid group-ID mix-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* grid = style name of this compute command
* group-ID = group ID for which grid cells to perform calculation on
* mix-ID = mixture ID to perform calculation on
* one or more values can be appended
* value = *n* or *nrho* or *nfrac* or *mass* or *massrho* or *massfrac* or *u* or *v* or *w* or *usq* or *vsq* or *wsq* of *ke* or *temp* or *erot* or *trot* or *evib* or *tvib* or *pxrho* or *pyrho* or *pzrho* or *kerho*
  
  .. parsed-literal::
  
       n = particle count
       nrho = number density
       nfrac = number fraction
       mass = mass
       massrho = mass density
       massfrac = mass fraction
       u = x component of velocity
       v = y component of velocity
       w = z component of velocity
       usq = x component of velocity squared
       vsq = y component of velocity squared
       wsq = z component of velocity squared
       ke = kinetic energy
       temp = temperature
       erot = rotational energy
       trot = rotational temperature
       evib = vibrational energy 
       tvib = vibrational temperature (classical definition)
       pxrho = x component of momentum density
       pyrho = y component of momentum density
       pzrho = z component of momentum density
       kerho = kinetic energy density



**Examples:**


.. parsed-literal::

   compute 1 grid all species n u v w usq vsq wsq
   compute 1 grid subset air n u v w

These commands will dump time averages for each species and each grid
cell to a dump file every 1000 steps:


.. parsed-literal::

   compute 1 grid all species n u v w usq vsq wsq
   fix 1 ave/grid 10 100 1000 c_1[\*]
   dump 1 grid all 1000 tmp.grid id f_1[\*]

**Description:**

Define a computation that calculates one or more values for each grid
cell in a grid cell group, based on the particles in the cell.  The
values are tallied separately for each group of species in the
specified mixture, as described in the Ouput section below.  See the
:doc:`mixture <mixture>` command for how a set of species can be
partitioned into groups.

Only grid cells in the grid group specified by *group-ID* are included
in the calculations.  See the :doc:`group grid <group>` command for info
on how grid cells can be assigned to grid groups.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump grid <dump>` command.

The values over many sampling timesteps can be averaged by the :doc:`fix ave/grid <fix_ave_grid>` command.  It does its averaging as if the
particles in the cell at each sampling timestep were combined together
into one large set of particles to compute the formulas below.

Note that for most of the values, this is a different form of
averaging than taking the values produced by the formulas below for a
single timestep, summing those values over the sampling timesteps, and
then dividing by the number of sampling steps.





The *n* value counts the number of particles in each group.  When
accumulated over multiple sampling steps, this value is normalized by
the number of sampling steps.

The *nrho* value computes the number density for the grid cell volume
due to particles in each group:


.. parsed-literal::

   Nrho = fnum/volume \* N

N is the number of particles (same as the *n* keyword), fnum is the
real/simulated particle ratio set by the :doc:`global fnum <global>`
command, and volume is the flow volume of the grid cell.  When
accumulated over multiple sampling steps, this value is normalized by
the number of sampling steps.  Note that if particle weighting is
enabled via the :doc:`global weight <global>` command, then the volume
used in the formula is divided by the weight assigned to the grid
cell.

The *nfrac* value computes the number fraction of particles in each
group:

Nfrac = Ngroup / Ntotal

Ngroup is the count of particles in the group and Ntotal is the total
number of particles in all groups in the mixture.  Note that this
total is not (necessarily) all particles in the cell.


----------


The *mass* value computes the average mass of particles in each group:


.. parsed-literal::

   Mass = Sum_i (mass_i) / N

where Sum\_i is a sum over particles in the group.

The *massrho* value computes the mass density for the grid cell volume
due to particles in each group:


.. parsed-literal::

   Massrho = fnum/volume \* Sum_i (mass_i)

where Sum\_i is a sum over particles in the group, fnum is the
real/simulated particle ratio set by the :doc:`global fnum <global>`
command, and volume is the flow volume of the grid cell.  When
accumulated over multiple sampling steps, this value is normalized by
the number of sampling steps.  Note that if particle weighting is
enabled via the :doc:`global weight <global>` command, then the volume
used in the formula is divided by the weight assigned to the grid
cell.

The *massfrac* value computes the mass fraction of particles in each
group:


.. parsed-literal::

   Massfrac = Sum_i (mass_i) / Masstotal

where Sum\_i is a sum over particles in the group and Masstotal is the
total mass of particles in all groups in the mixture.  Note that this
total is not (necessarily) the mass of all particles in the cell.


----------


The *u*\ , *v*\ , *w* values compute the components of the mass-weighted
average velocity of particles in each group:


.. parsed-literal::

   U = Sum_i (mass_i Vx_i) / Sum_i (mass_i)
   V = Sum_i (mass_i Vy_i) / Sum_i (mass_i)
   W = Sum_i (mass_i Vz_i) / Sum_i (mass_i)

This is the same as the center-of-mass velocity of particles in each
group.

The *usq*\ , *vsq*\ , *wsq* values compute the average mass-weighted
squared components of the velocity of particles in each group:


.. parsed-literal::

   Usq = Sum_i (mass_i Vx_i Vx_i) / Sum_i (mass_i)
   Vsq = Sum_i (mass_i Vy_i Vy_i) / Sum_i (mass_i)
   Wsq = Sum_i (mass_i Vz_i Vz_i) / Sum_i (mass_i)


----------


The *ke* value computes the average kinetic energy of particles in
each group:


.. parsed-literal::

   Vsq = Vx\*Vx + Vy\*Vy + Vz\*Vz
   KE = Sum_i (1/2 mass_i Vsq_i) / N

Note that this is different than the group's contribution to the
average kinetic energy of entire grid cells.  That can be calculated
by multiplying the *ke* quantity by the *n* quantity.

The *temp* value first computes the average kinetic energy of
particles in each group, as for the *ke* value.  This is then
converted to a temperature *T* by the following formula where *kB* is
the Boltzmann factor:


.. parsed-literal::

   Vsq = Vx\*Vx + Vy\*Vy + Vz\*Vz
   KE = Sum_i (1/2 mass_i Vsq_i) / N
   T = KE / (3/2 kB)

Note that this definition of temperature does not subtract out a net
streaming velocity for particles in the grid cell, so it is not a
thermal temperature when the particles have a non-zero streaming
velocity.  See the :doc:`compute thermal/grid <compute_thermal_grid>`
command to calculate thermal temperatures after subtracting out
streaming components of velocity.


----------


The *erot* value computes the average rotational energy of particles
in each group:


.. parsed-literal::

   Erot = Sum_i (erot_i) / N

Note that this is different than the group's contribution to the
average rotational energy of entire grid cells.  That can be
calculated by multiplying the *erot* quantity by the *n* quantity.

The *trot* value computes a rotational temperature by the following
formula where *kB* is the Boltzmann factor:


.. parsed-literal::

   Trot = (2/kB) Sum_i (erot_i) / Sum_i (dof_i)

Dof\_i is the number of rotational degrees of freedom for particle i.


----------


The *evib* value computes the average vibrational energy of particles
in each group:


.. parsed-literal::

   Evib = Sum_i (evib_i) / N

Note that this is different than the group's contribution to the
average vibrational energy of entire grid cells.  That can be
calculated by multiplying the *evib* quantity by the *n* quantity.

The *tvib* value computes a classical definition of vibrational
temperature, valid for continous distributions of vibrational energy,
by the following formula where *kB* is the Boltzmann factor:


.. parsed-literal::

   Tvib = (2/kB) Sum_i (evib_i) / Sum_i (dof_i)

Dof\_i is the number of vibrational degrees of freedom for particle i.


----------


The *pxrho*\ , *pyrho*\ , *pzrho* values compute components of momentum
density for the grid cell volume due to particles in each group:


.. parsed-literal::

   Pxrho = fnum/volume \* Sum_i (mass_i \* Vx_i)
   Pyrho = fnum/volume \* Sum_i (mass_i \* Vy_i)
   Pzrho = fnum/volume \* Sum_i (mass_i \* Vz_i)

where Sum\_i is a sum over particles in the group, fnum is the
real/simulated particle ratio set by the :doc:`global fnum <global>`
command, and volume is the flow volume of the grid cell.  When
accumulated over multiple sampling steps, this value is normalized by
the number of sampling steps.  Note that if particle weighting is
enabled via the :doc:`global weight <global>` command, then the volume
used in the formula is divided by the weight assigned to the grid
cell.

The *kerho* value computes the kinetic energy density for the grid
cell volume due to particles in each group:


.. parsed-literal::

   Vsq = Vx\*Vx + Vy\*Vy + Vz\*Vz
   KErho = fnum/volume \* Sum_i (mass_i \* Vsq_i)

where Sum\_i is a sum over particles in the group, fnum is the
real/simulated particle ratio set by the :doc:`global fnum <global>`
command, and volume is the flow volume of the grid cell.  When
accumulated over multiple sampling steps, this value is normalized by
the number of sampling steps.  Note that if particle weighting is
enabled via the :doc:`global weight <global>` command, then the volume
used in the formula is divided by the weight assigned to the grid
cell.





**Output info:**

This compute calculates a per-grid array, with the number of columns
equal to the number of values times the number of groups.  The
ordering of columns is first by values, then by groups.  I.e. if the
*n* and *u* values were specified as keywords, then the first two
columns would be *n* and *u* for the first group, the 3rd and 4th
columns would be *n* and *u* for the second group, etc.

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

The per-grid array values will be in the :doc:`units <units>`
appropriate to the individual values as described above.  *N* is
unitless.  *Nrho* is in 1/distance\^3 units for 3d simulations and
1/distance\^2 units for 2d simulations.  *Mass* is in mass units.
*Massrho* is in is in mass/distance\^3 units for 3d simulations and
mass/distance\^2 units for 2d simulations.  *U*\ , *v*\ , and *w* are in
velocity units.  *Usq*\ , *vsq*\ , and *wsq* are in velocity squared
units.  *Ke*\ , *erot*\ , and *evib* are in energy units.  *Temp* and
*trot* and *tvib* are in temperature units.  *Pxrho*\ , *pyrho*\ , *pzrho*
are in momentum/distance\^3 units for 3d simulations and
momentum/distance\^2 units for 2d simulations, where momentum is in
units of mass\*velocity.  *Kerho* is in units of energy/distance\^3
units for 3d simulations and energy/distance\^2 units for 2d
simulations.


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

:doc:`fix ave/grid <fix_ave_grid>`, :doc:`dump grid <dump>`, :doc:`compute thermal/grid <compute_thermal_grid>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
