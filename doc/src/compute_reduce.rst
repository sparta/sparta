.. index:: compute reduce

compute reduce command
======================

compute reduce/kk command
=========================

**Syntax:**


.. parsed-literal::

   compute ID reduce mode input1 input2 ... keyword args ...

* ID is documented in :doc:`compute <compute>` command
* reduce = style name of this compute command
* mode = *sum* or *min* or *max* or *ave* or *sumsq* or *avesq* or *sum-area* or *ave-area*
* one or more inputs can be listed
* input = x, y, z, vx, vy, vz, ke, erot, evib, c\_ID, c\_ID[N], f\_ID, f\_ID[N], v\_name, p\_name, p\_name[N], g\_name, g\_name[N], s\_name, s\_name[N]
  
  .. parsed-literal::
  
       x,y,z,vx,vy,vz = particle position or velocity component
       ke,erot,evib = particle energy component
       c_ID = per-particle or per-grid vector calculated by a compute with ID
       c_ID[N] = Nth column of per-particle or per-grid array calculated by a compute with ID, N can include wildcard (see below)
       f_ID = per-particle or per-grid or per-surf vector calculated by a fix with ID
       f_ID[N] = Nth column of per-particle or per-grid or per-surf array calculated by a fix with ID, N can include wildcard (see below)
       v_name = per-particle or per-grid or per-surf vector calculated by a particle-style or grid-style or surf-style variable with name
       p_name = custom per-particle vector with name
       p_name[N] = Nth column of per-particle custom array with name, N can include wildcard (see below)
       g_name = custom per-grid vector with name
       g_name[N] = Nth column of per-grid custom array with name, N can include wildcard (see below)
       s_name = custom per-surf vector with name
       s_name[N] = Nth column of per-surf custom array with name, N can include wildcard (see below)

* zero or more keyword/args pairs may be appended
* keyword = *replace* or *subset*
  
  .. parsed-literal::
  
       replace args = vec1 vec2
         vec1 = reduced value from this input vector will be replaced
         vec2 = replace it with vec1[N] where N is index of max/min value from vec2
       subset arg = subsetID
         subsetID = mixture-ID or grid group-ID or surface group-ID



**Examples:**


.. parsed-literal::

   compute 1 reduce sum c_grid[\*]
   compute 2 reduce min f_ave v_myKE subset trace_species
   compute 3 reduce max c_mine[1] c_mine[2] c_temp replace 1 3 replace 2 3

These commands will include the average grid cell temperature, across
all grid cells, in the stats output:


.. parsed-literal::

   compute 1 temp
   compute	2 grid all all temp
   compute 3 reduce ave c_2[1]
   stats_style step c_temp c_3

**Description:**

Define a calculation that "reduces" one or more vector inputs into
scalar values, one per listed input.  The inputs can be per-particle
or per-grid or per-surf quantities; they cannot be global quantities.
Particle attributes are per-particle quantities.
:doc:`Computes <compute>` and :doc:`fixes <fix>` may generate any of the
three kinds of quantities.  :doc:`Particle-style, grid-style, and surf-style variables <variable>` generate per-particle, per-grid,
or per-surf quantities respectively.  Custom attributes can be
per-particle, per-grid, or per-surf quantities.  See the
:doc:`variable <variable>` command and its special functions which can
perform the same operations as the compute reduce command on global
vectors.

IMPORTANT NOTE: All inputs to a compute reduce command must be the
same type: per-particle, per-grid, or per-surf.  You can use the
command multiple times if you need to reduce values of different
types.

The reduction operation is specified by the *mode* setting.  The *sum*
option adds the values in the vector into a global total.  The *min*
or *max* operations find the minimum or maximum value across all
vector values.  The *ave* operation adds the vector values into a
global total, then divides by the number of values in the vector.  The
*sumsq* operation sums the square of the values in the vector into a
global total.  The *avesq* oepration does the same as *sumsq*\ , then
divdes the sum of squares by the number of values.  These two
operations can be useful for calculating the variance of some
quantity, e.g. variance = sumsq - ave\^2.

The *sum-area* or *ave-area* options can only be used for per-surf
inputs.  Both options multiply each per-surf value by the area of the
surface element (triangle in 3d, line segment in 2d) and sum the
resulting values over all surface elements.  That is the output for
the *sum-area* option.  For the *ave-area* option the summed value is
divided by the summed area of all elements.  Note that both of these
options are designed to work with flux values (e.g. mass per area per
time) produced by the :doc:`compute surf <compute_surf>` command with
its default *norm* = yes option.


----------


Each listed input vector is operated on independently.

Each listed input vector can be a particle attribute or can be the
result of a :doc:`compute <compute>` or :doc:`fix <fix>` or the evaluation
of a :doc:`variable <variable>`.  Or it can be a custom attribute of a
particle, grid cell, or surface element.

Note that for values from a compute or fix or custom attribute, the
bracketed index I can be specified using a wildcard asterisk with the
index to effectively specify multiple values.  This takes the form "\*"
or "\*n" or "n\*" or "m\*n".  If N = the size of the vector (for *mode* =
scalar) or the number of columns in the array (for *mode* = vector),
then an asterisk with no numeric values means all indices from 1 to N.
A leading asterisk means all indices from 1 to n (inclusive).  A
trailing asterisk means all indices from n to N (inclusive).  A middle
asterisk means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual columns of the array
had been listed one by one.  E.g. these 2 compute reduce commands are
equivalent, since the :doc:`compute grid <compute_grid>` command creates
a per-grid array with 3 columns:


.. parsed-literal::

   compute myGrid grid all all u v w
   compute 2 all reduce min c_myGrid[\*]
   compute 2 all reduce min c_myGrid[1] c_myGrid[2] c_myGrid[3]


----------


The particle attributes x,y,z,vx,vy,vz are position and velocity
components.  The ke,erot,evib attributes are for kinetic, rotational,
and vibrational energy of particles.

If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  Computes can generate
per-particle or per-grid quantities.  See the individual
:doc:`compute <compute>` doc page for details.  If no bracketed integer
is appended, the vector calculated by the compute is used.  If a
bracketed integer is appended, the Nth column of the array calculated
by the compute is used.  Users can also write code for their own
compute styles and :doc:`add them to SPARTA <Section_modify>`.  See the
discussion above for how N can be specified with a wildcard asterisk
to effectively specify multiple values.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  Fixes can generate
per-particle or per-grid or per-surf quantities.  See the individual
:doc:`fix <fix>` doc page for details.  Note that some fixes only
produce their values on certain timesteps, which must be compatible
with when this compute references the values, else an error results.
If no bracketed integer is appended, the vector calculated by the fix
is used.  If a bracketed integer is appended, the Nth column of the
array calculated by the fix is used.  Users can also write code for
their own fix style and :doc:`add them to SPARTA <Section_modify>`.  See
the discussion above for how N can be specified with a wildcard
asterisk to effectively specify multiple values.

If a value begins with "v\_", a variable name must follow which has
been previously defined in the input script.  It must be a
:doc:`particle-style or grid-style or surf-style variable <variable>`.
These styles define formulas which can reference stats keywords or
invoke other computes, fixes, or variables when they are evaluated.
Particle-style variables can also reference various per-particle
attributes (position, velocity, etc).  So these variables are a very
general means of creating per-particle or per-grid or per-surf
quantities to reduce.

If a value begins with "p\_" or "g\_" or "s\_", then a custom
per-particle, per-grid, or per-surf attribute with the specified name
is used.  Particles, grid cells, and surface elements can have custom
attributes which store either single or multiple values per particle,
per grid cell, or per surface element.  They can be defined and
initialized in data files, e.g. via the :doc:`read\_surf <read_surf>`
command.  Or they can be defined and used by specific commands,
e.g. :doc:`fix ambipolar <fix_ambipolar>` or :doc:`fix surf/temp <fix_surf_temp>` or :doc:`surf\_react adsorb <surf_react_adsorb>`.  The name of each attribute is set by
the user or defined by the command.  See `Section 6.17 <Section_howto.html#howto_17>`_ for more discussion of custom
attributes.

If no bracketed integer is appended, the custom attribute must be a
per-particle, per-grid, or per-surf vector (single value).  If a
bracketed integer is appended, the custom attribute must be a
per-particle, per-grid, or per-surf arayy (multiple values) and the
Nth column of the custom array is used.  See the discussion above for
how N can be specified with a wildcard asterisk to effectively specify
multiple values.


----------


If the *replace* keyword is used, two indices *vec1* and *vec2* are
specified, where each index ranges from 1 to the # of input values.
The replace keyword can only be used if the *mode* is *min* or *max*\ .
It works as follows.  A min/max is computed as usual on the *vec2*
input vector.  The index N of that value within *vec2* is also stored.
Then, instead of performing a min/max on the *vec1* input vector, the
stored index is used to select the Nth element of the *vec1* vector.

Here is an example which prints out both the grid cell ID and number
of particles for the grid cell with the maximum number of particles:


.. parsed-literal::

   compute 1 property/grid id
   compute	2 grid all n
   compute	3 reduce max c_1 c_2[1] replace 1 2
   stats_style step c_temp c_3[1] c_3[2]

The first two input values in the compute reduce command are vectors
with the ID and particle count of each grid cell.  Instead of taking
the max of the ID vector, which does not yield useful information in
this context, the *replace* keyword will extract the ID for the grid
cell which has the maximum number of particles.  This ID and the
cell's particle count will be printed with the statistical output.

Note that the *replace* keyword can be used multiple times with
different pairs of indices.


----------


The *subset* keyword allows selection of a subset of each input
vectors quantities to be used for the reduce operation.  This may
affect all of the reduction operations.  E.g. the ave and avesq
operations will become averages for only a subset of numerical values.

If inputs are per-particle values, then a mixture ID should be
specified.  Only particle species belonging to the mixture will be
included in the calculations.  See the :doc:`mixture <mixture>` command
for how a set of species is included in a mixture.

If inputs are per-grid values, then a grid group ID should be
specified.  Only grid cells in the grid group will be included in the
calculations.  See the :doc:`group grid <group>` command for info on how
grid cells can be assigned to grid groups.

If inputs are per-surf values, then a surface group ID should be
specified.  Only surface elements in the surface group will be
included in the calculations.  See the :doc:`group surf <group>` command
for info on how surface elements can be assigned to surface groups.

IMPORTANT NOTE: If computes or fixes are used as inputs to compute
reduce, they may define their own subsets of particles, grid cells, or
surface elements which contribute to their output.  Typically output
from those computes or fixes will be zero for grid cells or surface
elements not in the grid or surface group specified for those
commands.  Thus you may want to use an argument for the *subset*
keyword which is consistent with the inputs, but that is not required.


----------


If a single input is specified this compute produces a global scalar
value.  If multiple inputs are specified, this compute produces a
global vector of values, the length of which is equal to the number of
inputs specified.


----------


**Output info:**

This compute calculates a global scalar if a single input value is
specified or a global vector of length N where N is the number of
inputs, and which can be accessed by indices 1 to N.  These values can
be used by any command that uses global scalar or vector values from a
compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_ for an
overview of SPARTA output options.

The scalar or vector values will be in whatever :doc:`units <units>` the
quantities being reduced are in.

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

For the *kk* style, an input is reduced on the device when its source is
device-resident: the explicit per-particle attributes (*x*, *vx*, *ke*,
*erot*, *evib*, etc), per-particle and per-grid custom attributes, and
the per-particle or per-grid output of a compute or fix that itself has
a *kk* variant.  Every other input -- per-surf inputs, *v_name*
variables, and the output of a non-Kokkos compute or fix -- is reduced
on the host, after the data it reads is copied back.  This is a
performance distinction only; the result is the same either way.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the `-suffix command-line switch <Section_start.html#start_7>`_ when you invoke SPARTA, or you can
use the :doc:`suffix <suffix>` command in your input script.

See the :doc:`Accelerating SPARTA <Section_accelerate>` section of the
manual for more instructions on how to use the accelerated styles
effectively.


----------


**Restrictions:** none

**Related commands:**

:doc:`compute <compute>`, :doc:`fix <fix>`, :doc:`variable <variable>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
