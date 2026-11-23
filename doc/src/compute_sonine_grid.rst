.. index:: compute sonine/grid

compute sonine/grid command
===========================

compute sonine/grid/kk command
==============================

**Syntax:**


.. parsed-literal::

   compute ID sonine/grid group-ID mix-ID keyword values ...

* ID is documented in :doc:`compute <compute>` command
* sonine/grid = style name of this compute command
* group-ID = group ID for which grid cells to perform calculation on
* mix-ID = mixture ID to perform calculation on
* one or more keywords may be appended, multiple times
* keyword = *a* or *b*
* values = values for specific keyword
  
  .. parsed-literal::
  
       a args = dim order = sonine A moment
         dim = x or y or z
         order = number from 1 to 5
       b args = dim2 order = sonine B moment
         dim2 = xx or yy or zz or xy or yz or xz
         order = number from 1 to 5



**Examples:**


.. parsed-literal::

   compute 1 sonine/grid all air a x 5 b xy 5
   compute 1 sonine/grid subset air a x 5

These commands will dump time averaged sonine moments for each
species and each grid cell to a dump file every 1000 steps:


.. parsed-literal::

   compute 1 sonine/grid all species a x 5 b xy 5
   fix 1 ave/grid 10 100 1000 c_1[\*]
   dump 1 grid all 1000 tmp.grid id f_1[\*]

**Description:**

Define a computation that calculates the sonine moments of the
velocity distribution of the particles in each grid cell in a grid
cell group.  The values are tallied separately for each group of
species in the specified mixture, as described in the Output section
below.  See the mixture command for how a set of species can be
partitioned into groups.

Only grid cells in the grid group specified by *group-ID* are included
in the calculations.  See the :doc:`group grid <group>` command for info
on how grid cells can be assigned to grid groups.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump grid <dump>` command.

The values over many sampling timesteps can be averaged by the :doc:`fix ave/grid <fix_ave_grid>` command.  It does its averaging as if the
particles in the cell at each sampling timestep were combined together
into one large set of particles to compute the A,B formulas below.

Note however that the center-of-mass (COM) velocity that is subtracted
from each particle to yield a squared thermal velocity Csq for each
particle, as described below, is the COM velocity for only the
particles in the current timestep.  When time-averaging it is NOT the
COM velocity for all particles across all timesteps.

Note that this is a different form of averaging than taking the values
produced by the formulas below for a single timestep, summing those
values over the sampling timesteps, and then dividing by the number of
sampling steps.


----------


Calculation of both the A and B sonine moments is done by first
calcuating the center-of-mass (COM) velocity of particles for each
group within a grid cell.  This is done as follows:


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
particles in its group and cell.  This allows computation of Csq for
each particle which is used in the formulas below to calculate the
sonine moments.


----------


The *a* keyword calculates the average of one or more sonine A moments
for all particles in each group:


.. parsed-literal::

   A1 = Sum_i (mass_i \* Vdim \* pow(Csq,1)) / Sum_i (mass_i)
   A2 = Sum_i (mass_i \* Vdim \* pow(Csq,2)) / Sum_i (mass_i)
   A3 = Sum_i (mass_i \* Vdim \* pow(Csq,3)) / Sum_i (mass_i)
   A4 = Sum_i (mass_i \* Vdim \* pow(Csq,4)) / Sum_i (mass_i)
   A5 = Sum_i (mass_i \* Vdim \* pow(Csq,5)) / Sum_i (mass_i)

Vdim is Vx or Vy or Vz as specified by the *dim* value.  *Csq* is the
squared thermal velocity of the particle, as in the COM equations
above.  The number of moments computed is specified by the *order*
value, e.g. for order = 3, the first 3 moments are computed, which
leads to 3 columns of output as explained below.

The *b* keyword calculates the average of one or more sonine B moments
for all particles in each group:


.. parsed-literal::

   B1 = Sum_i (mass_i \* Vdim1 \* Vdim2 \* pow(Csq,1)) / Sum_i (mass_i)
   B2 = Sum_i (mass_i \* Vdim1 \* Vdim2 \* pow(Csq,2)) / Sum_i (mass_i)
   B3 = Sum_i (mass_i \* Vdim1 \* Vdim2 \* pow(Csq,3)) / Sum_i (mass_i)
   B4 = Sum_i (mass_i \* Vdim1 \* Vdim2 \* pow(Csq,4)) / Sum_i (mass_i)
   B5 = Sum_i (mass_i \* Vdim1 \* Vdim2 \* pow(Csq,5)) / Sum_i (mass_i)

Vdim is Vx or Vy or Vz as specified by the *dim* value.  *Csq* is the
squared thermal velocity of the particle, as in the COM equations
above.  The number of moments computed is specified by the *order*
value, e.g. for order = 2, the first 2 moments are computed, which
leads to 2 columns of output as explained below.


----------


**Output info:**

This compute calculates a per-grid array, with the number of columns
equal to the number of values times the number of groups.  The
ordering of columns is first by values, then by groups.  I.e. if the
*a z 3* and *b xy 2* moments were specified as keywords, then the 1st
thru 3rd columns would be the A1, A2, A3 moments of the first group,
the 4th and 5th columns would be the B1 and B2 moments of the first
group, the 6th thru 8th columns would be the A1, A2, A3 moments of the
2nd group, etc.

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

Grid cells not in the specified *group-ID* will have zeroes for all
their values.

The array can be accessed by any command that uses per-grid values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The per-grid array values will be in the :doc:`units <units>`
appropriate to the individual values as described above.  These are
units like velocity cubed or velocity to the 6th power.


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

:doc:`fix ave/grid <fix_ave_grid>`, :doc:`dump grid <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
