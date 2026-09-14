.. index:: compute boundary

compute boundary command
========================

**Syntax:**


.. parsed-literal::

   compute ID boundary mix-ID value1 value2 ...

   compute ID boundary/kk mix-ID value1 value2 ...

* ID is documented in :doc:`compute <compute>` command
* boundary = style name of this compute command
* mix-ID = mixture ID to perform calculation on
* one or more values can be appended
* value = *n* or *nwt* or *nflux* or *mflux* or *press* or *shx* or *shy* or *shz* or *ke* or *erot* or *evib* or *etot*
  
  .. parsed-literal::
  
       n = count of particles hitting boundary
       nwt = weighted count of particles hitting boundary
       nflux = flux of particles on boundary
       mflux = flux of mass on boundary
       press = magnitude of normal pressure on boundary
       shx,shy,shz = components of shear stress on boundary
       ke = flux of particle kinetic energy on boundary 
       erot = flux of particle rotational energy on boundary 
       evib = flux of particle vibrational energy on boundary 
       etot = flux of particle total energy on boundary



**Examples:**


.. parsed-literal::

   compute 1 boundary all n press eng
   compute mine boundary species press shx shy shz

These commands will print values for the current timestep for 
the xlo and xhi boundaryies, as part of statistical output:


.. parsed-literal::

   compute 1 boundary all n press
   stats_style step np c_1[1][1] c_1[1][2] c_1[2][1] c_1[2][2]

These commands will dump time averages for each species and each
boundary to a file every 1000 steps:


.. parsed-literal::

   compute 1 boundary species n press shx shy shz
   fix 1 ave/time 10 100 1000 c_1[\*] mode vector file tmp.boundary

**Description:**

Define a computation that calculates one or more values for each
boundary (i.e. face) of the simulation box, based on the particles
that cross or collide with the boundary.  The values are summed for
each group of species in the specified mixture.  See the
:doc:`mixture <mixture>` command for how a set of species can be
partitioned into groups.

Note that depending on the settings for the :doc:`boundary <boundary>`
command, when a particle collides with a boundary, it can exit the
simulation box (outflow), re-enter from the other side (periodic),
reflect specularly from the boundary, or interact with it as if it
were a surface.  In the surface case, the incident particle may bounce
off (possibly as a different species), be captured by the boundary
(vanish), or a 2nd particle can also be emitted.  The formulas below
account for all these possible scenarios.  As an example, the pressure
exerted on an outflow boundary versus a specularly reflecting boundary
is different, since in the former case there is no net momentum flux
back into the simulation box by reflected particles.

Also note that all values for a boundary collision are tallied based
on the species group of the incident particle.  Quantities associated
with outgoing particles are part of the same tally, even if they are
in different species groups.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`stats\_style <stats_style>` command.

The values over many sampling timesteps can be averaged by the :doc:`fix ave/time <fix_ave_time>` command.  It does its averaging as if the
particles striking the boundary at each sampling timestep were
combined together into one large set to compute the formulas below.
The answer is then divided by the number of sampling timesteps if it
is not otherwise normalized by the number of particles.  Note that in
general this is a different normalization than taking the values
produced by the formulas below for a single timestep, summing them
over the sampling timesteps, and then dividing by the number of
sampling steps.  However for the current values listed below, the two
normalization methods are the same.

NOTE: If particle weighting is enabled via the :doc:`global weight <global>` command, then all of the values below are scaled
by the weight assigned to the grid cell in which the particle
collision with the boundary occurs.  The only exception is the the *n*
value, which is NOT scaled by the weight; it is a simple count of
particle crossings or collisions with the boundary.


----------


The *n* value counts the number of particles in the group crossing or
colliding with the boundary.

The *nwt* value counts the number of particles in the group crossing
or colliding with the boundary and weights the count by the weight
assigned to the grid cell in which the particle collision with the
boundary occurs.  The *nwt* quantity will only be different than *n*
if particle weighting is enabled via the :doc:`global weight <global>`
command.

The *nflux* value calculates the number flux imparted to the boundary by
particles in the group.  This is computed as


.. parsed-literal::

   Nflux = N / (A \* dt / fnum)

where N is the number of all contributing particles, normalized by
A = the area of the surface element, dt = the timestep, and fnum = the
real/simulated particle ratio set by the :doc:`global fnum <global>`
command.

The *mflux* value calculates the mass flux imparted to the boundary by
particles in the group.  This is computed as


.. parsed-literal::

   Mflux = Sum_i (mass_i) / (A \* dt / fnum)

where the sum is over all contributing particle masses, normalized by
the area of the surface element, dt and fnum as defined before.

The *press* value calculates the pressure *P* exerted on the boundary
in the normal direction by particles in the group, such that outward
pressure is positive.  This is computed as


.. parsed-literal::

   p_delta = mass \* (V_post - V_pre)
   P = Sum_i (p_delta_i dot N) / (A \* dt / fnum)

where A, dt, fnum are defined as before.  P\_delta is the change in
momentum of a particle, whose velocity changes from V\_pre to V\_post
when colliding with the boundary.  The pressure exerted on the
boundary is the sum over all contributing p\_delta dotted into the
normal N of the boundary which is directed into the box, normalized by
A = the area of the boundary face and dt = the timestep and fnum = the
real/simulated particle ratio set by the :doc:`global fnum <global>`
command.

The *shx*\ , *shy*\ , *shz* values calculate the shear pressure components
Sx, Sy, Sz extered on the boundary in the tangential direction to its
normal by particles in the group, with respect to the x, y, z
coordinate axes.  These are computed as


.. parsed-literal::

   p_delta = mass \* (V_post - V_pre)
   p_delta_t = p_delta - (p_delta dot N) N
   Sx = - Sum_i (p_delta_t_x) / (A \* dt / fnum)
   Sy = - Sum_i (p_delta_t_y) / (A \* dt / fnum)
   Sz = - Sum_i (p_delta_t_z) / (A \* dt / fnum)

where p\_delta, V\_pre, V\_post, N, A, dt, and fnum are defined as
before.  P\_delta\_t is the tangential component of the change in
momentum vector p\_delta of a particle.  P\_delta\_t_x (and y,z) are its
x, y, z components.

The *ke* value calculates the kinetic energy flux *Eflux* imparted to
the boundary by particles in the group, such that energy lost by a
particle is a positive flux.  This is computed as


.. parsed-literal::

   e_delta = 1/2 mass (V_post\^2 - V_pre\^2)
   Eflux = - Sum_i (e_delta) / (A \* dt / fnum)

where e\_delta is the kinetic energy change in a particle, whose
velocity changes from V\_pre to V\_post when colliding with the
boundary.  The energy flux imparted to the boundary is the sum over
all contributing e\_delta, normalized by A = the area of the boundary
face and dt = the timestep and fnum = the real/simulated particle
ratio set by the :doc:`global fnum <global>` command.

The *erot* value calculates the rotational energy flux *Eflux*
imparted to the boundary by particles in the group, such that energy
lost by a particle is a positive flux.  This is computed as


.. parsed-literal::

   e_delta = Erot_post - Erot_pre
   Eflux = - Sum_i (e_delta) / (A \* dt / fnum)

where e\_delta is the rotational energy change in a particle, whose
internal rotational energy changes from Erot\_pre to Erot\_post when
colliding with the boundary.  The flux equation is the same as for the
*ke* value.

The *evib* value calculates the vibrational energy flux *Eflux*
imparted to the boundary by particles in the group, such that energy
lost by a particle is a positive flux.  This is computed as


.. parsed-literal::

   e_delta = Evib_post - Evib_pre
   Eflux = - Sum_i (e_delta) / (A \* dt / fnum)

where e\_delta is the vibrational energy change in a particle, whose
internal vibrational energy changes from Evib\_pre to Evib\_post when
colliding with the boundary.  The flux equation is the same as for the
*ke* value.

The *etot* value calculates the total energy flux imparted to the
boundary by particles in the group, such that energy lost by a
particle is a positive flux.  This is simply the sum of kinetic,
rotational, and vibrational energies.  Thus the total energy flux is
the sum of what is computed by the *ke*\ , *erot*\ , and *evib* values.


----------


**Output info:**

This compute calculates a global array, with the number of columns
equal to the number of values times the number of groups.  The
ordering of columns is first by values, then by groups.  I.e. if the
*n* and *u* values were specified as keywords, then the first two
columns would be *n* and *u* for the first group, the 3rd and 4th
columns would be *n* and *u* for the second group, etc.  The number of
rows is 4 for a 2d simulation for the 4 faces (xlo, xhi, ylo, yhi),
and it is 6 for a 3d simulation (xlo, xhi, ylo, yhi, zlo, zhi).

The array can be accessed by any command that uses global array values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The array values will be in the :doc:`units <units>` appropriate to the
individual values as described above.  *N* is unitless. *Press*\ ,
*shx*\ , *shy*\ , *shz* are in pressure units.  *Ke*\ , *erot*\ , *evib*\ , and
*etot* are in energy/area-time units for 3d simulations and
energy/length-time units for 2d simulations.


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


**Restrictions:**

If specified with a *kk* suffix, this compute can be used no more than
twice in the same input script (active at the same time).

**Related commands:**

:doc:`fix ave/time <fix_ave_time>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
