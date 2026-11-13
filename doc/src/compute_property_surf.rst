.. index:: compute property/surf

compute property/surf command
=============================

compute property/surf/kk command
================================

**Syntax:**


.. parsed-literal::

   compute ID property/surf group-ID input1 input2 ...

* ID is documented in :doc:`compute <compute>` command
* property/surf = style name of this compute command
* group-ID = group ID for which surface elements to perform calculation on
* input = one or more surface element attributes
  
  .. parsed-literal::
  
       possible attributes = id, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v2y, v3z, xc, yc, zc, area, normx, normy, normz

  
  .. parsed-literal::
  
       id = surface element ID
       v1x,v1y,v1z = coords of first line end point or triangle corner point
       v3x,v2y,v2z = coords of second line end point or triangle corner point
       v3x,v3y,v3z = coords of third triangle corner point
       xc,yc,zc = coords of center of line segment or triangle
       area = length of line segment or area of triangle
       normx, normy, normz = unit normal vector for line segment or triangle



**Examples:**


.. parsed-literal::

   compute 1 property/surf all id xc yc zc

**Description:**

Define a computation that simply stores surface element attributes for
each explicit surface element in a surface group.  This is useful for
values which can be used by other `output commands <Section_howto.html#howto_4>`_ that take computes as inputs.
See for example, the :doc:`compute reduce <compute_reduce>`, :doc:`fix ave/surf <fix_ave_surf>`, :doc:`dump surf <dump>`, and :doc:`surf-style variable <variable>` commands.

Only surface elements in the surface group specified by *group-ID* are
included in the calculation.  See the :doc:`group surf <group>` command
for info on how surface elements can be assigned to surface groups.

This command can only be used for simulations with explicit surface
elements.  Explicit surface elements are triangles for 3d simulations
and line segments for 2d simulations.  Unlike implicit surface
elements, each explicit triangle or line segment may span multiple
grid cells.  See `Section 6.9 <Section_howto.html#howto_9>`_ of the
manual for details.


----------


*Id* is the surface element ID, as defined in the surface data file
read by the :doc:`read\_surf <read_surf>` comand.

The *v1x*\ , *v1y*\ , *v1z* attributes are the coordinates of the first
end point of a line segment (2d) or first corner point of a triangle
(3d).  Likewise, the *v2x*\ , *v2y*\ , *v2z* attributes are the
coordinates of the second end point of a line segment (2d) or second
corner point of a triangle (3d).  The *v3x*\ , *v3y*\ , *v23z* attributes
are the coordinates of the third corner point of a triangle (3d).

The *xc*\ , *yc*\ , *zc* attributes are the coordinates of the center
point of a line segment or tringle.

The *area* attribute is the length of a line segment (distance units
in 2d), or area of a triangle (area units in 3d).

The *normx*\ , *normy*\ , *normz* attributes are components of a unit
normal perpendicular to the line segment or face of the trangle. It
points into the flow volume of the simulation.


----------


**Output info:**

This compute calculates a per-surf vector or per-surf array depending
on the number of input values.  If a single input is specified, a
per-surf vector is produced.  If two or more inputs are specified, a
per-surf array is produced where the number of columns = the number of
inputs.

This compute performs calculations for each explicit surface element
in the simulation.

Surface elements not in the specified *group-ID* will output zeroes
for all their values.

The vector or array can be accessed by any command that uses per-surf
values from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_ for an overview of SPARTA output
options.

The vector or array values will be in whatever :doc:`units <units>` the
corresponding attribute is in, e.g. distance units for *v1x* or *xc*\ ,
length units for *area* in 2d, area units for *area* in 3d.


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

For 2d simulations, none of the attributes which refer to the 3rd
dimension may be used.  Likewise *v3x*\ , *v3y*\ , *v3z* may not be used
since they refer to triangles.

**Related commands:**

:doc:`dump surf <dump>`, :doc:`fix ave/surf <fix_ave_surf>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
