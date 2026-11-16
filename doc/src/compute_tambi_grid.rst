.. index:: compute tambi/grid

compute tambi/grid command
==========================

**Syntax:**


.. parsed-literal::

   compute ID tambi/grid group-ID

* ID is documented in :doc:`compute <compute>` command
* tambi/grid = style name of this compute command
* group-ID = group ID for which grid cells to perform calculation on

**Examples:**


.. parsed-literal::

   compute 1 tambi/grid all
   compute te tambi/grid all
   compute ter reduce ave c_te

These commands will dump the time-averaged ambipolar electron
temperature for each grid cell to a dump file every 1000 steps:


.. parsed-literal::

   compute te tambi/grid all
   fix 1 ave/grid 10 100 1000 c_te
   dump 1 grid all 1000 tmp.grid id f_1

**Description:**

Define a computation that calculates the ambipolar electron
temperature in each grid cell in a grid cell group.  It requires the
:doc:`fix ambipolar <fix_ambipolar>` command: every ambipolar ion
carries the velocity of its attached electron in the *velambi* custom
particle attribute, initialized when the ion is created and updated by
the electron's elastic collisions with heavy particles (see
:doc:`collide\_modify ambipolar yes <collide_modify>`).  The electron
velocity distribution in a cell therefore defines a kinetic electron
temperature


.. parsed-literal::

   T_e = m_e <\|v_e - <v_e>\|^2> / (3 kB)

where m\_e is the mass of the ambipolar electron species, the averages
run over the electrons attached to the ambipolar ions in the cell, and
the center-of-mass electron velocity is subtracted so that a streaming
electron gas is not misread as a hot one.

This is the collisional-exchange electron temperature of the ambipolar
approximation: electron energy changes only through elastic collisions
with heavy particles, so T\_e relaxes toward the heavy-particle
temperature at the mass-ratio-limited rate ~2 m\_e/M per electron
collision.  Field (Joule) heating, electron heat conduction, and
inelastic electron-energy losses are not represented.  Cells
containing fewer than 2 ambipolar ions report a temperature of 0.0.

Only grid cells in the grid group specified by *group-ID* are included
in the calculation.  See the :doc:`group grid <group>` command for info
on how grid cells can be assigned to grid groups.


----------


**Output info:**

This compute calculates a per-grid vector (the electron temperature in
each cell).

This compute performs calculations for all flavors of child grid cells
in the simulation, which includes unsplit, cut, split, and sub cells.

The vector can be accessed by any command that uses per-grid values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The per-grid vector values will be in temperature :doc:`units <units>`
(i.e. Kelvin).

**Restrictions:**

Requires the :doc:`fix ambipolar <fix_ambipolar>` command.

**Related commands:**

:doc:`fix ambipolar <fix_ambipolar>`, :doc:`compute thermal/grid <compute_thermal_grid>`, :doc:`compute grid <compute_grid>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
