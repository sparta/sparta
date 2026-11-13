.. index:: compute fft/grid

compute fft/grid command
========================

compute fft/grid/kk command
===========================

**Syntax:**


.. parsed-literal::

   compute ID fft/grid value1 value2 ... keyword args ...

* ID is documented in :doc:`compute <compute>` command
* fft/grid = style name of this compute command
* one or more values can be appended
* value = c\_ID, c\_ID[N], f\_ID, f\_ID[N], v\_name
  
  .. parsed-literal::
  
       c_ID = per-grid vector calculated by a compute with ID
       c_ID[I] = Ith column of per-grid array calculated by a compute with ID
       f_ID = per-grid vector calculated by a fix with ID
       f_ID[I] = Ith column of per-grid or array calculated by a fix with ID
       v_name = per-grid vector calculated by a grid-style variable with name

* zero or more keyword/arg pairs can be appended
* keyword = *sum* or *scale* or *conjugate* or *kmag*
  
  .. parsed-literal::
  
       sum = yes or no to sum all FFTs into a single output
       scale = sfactor = numeric value to scale results by
       conjugate = yes or no = perform complex conjugate multiply or not
       kx = yes or no = calculate x-component of wavelength or not
       kx = yes or no = calculate y-component of wavelength or not
       kx = yes or no = calculate z-component of wavelength or not
       kmag = yes or no = calculate wavelength magnitude or not



**Examples:**


.. parsed-literal::

   compute 1 fft/grid c_1

These commands will dump FFTs of instantaneous and time-averaged
velocity components in each grid cell to a dump file every 1000 steps:


.. parsed-literal::

   compute 1 grid all u v w
   fix 1 ave/grid 10 100 1000 c_1
   compute 2 fft/grid f_1[1] f_1[2] f_1[3]
   dump 1 grid all 1000 tmp.grid id c_2 f_1

**Description:**

Define a computation that performs forward FFTs on per-grid values.  
This can be useful, for example, in calculating the energy spectrum
of a turbulent flow.

The defined grid must be a regular one-level grid (not hierarchical)
with an even number of grid cells in each dimension.  Depending on the
:doc:`dimension <dimension>` of the simulation, either 2d or 3d FFTs
will be performed.  Because FFTs assume a periodic field, the
simulation domain should be periodic in all dimensions, as set by the
:doc:`boundary <boundary>` command, though SPARTA does not check for
that.

The results of this compute can be used by different commands in
different ways.  The values for a single timestep can be output by the
:doc:`dump grid <dump>` command.  The values over many sampling
timesteps can be averaged by the :doc:`fix ave/grid <fix_ave_grid>`
command.

A forward FFT is perfomed on each input value independently.

Each listed input can be the result of a :doc:`compute <compute>` or
:doc:`fix <fix>` or the evaluation of a :doc:`variable <variable>`, all of
which must generate per-grid quantities.

If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  The compute must generate a
per-grid vector or array.  See the individual :doc:`compute <compute>`
doc page for details.  If no bracketed integer is appended, the vector
calculated by the compute is used.  If a bracketed integer is
appended, the Ith column of the array calculated by the compute is
used.  Users can also write code for their own compute styles and :doc:`add them to SPARTA <Section_modify>`.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  The fix must generate a
per-grid vector or array.  See the individual :doc:`fix <fix>` doc page
for details.  Note that some fixes only produce their values on
certain timesteps, which must be compatible with when this compute
references the values, else an error results.  If no bracketed integer
is appended, the vector calculated by the fix is used.  If a bracketed
integer is appended, the Ith column of the array calculated by the fix
is used.  Users can also write code for their own fix style and :doc:`add them to SPARTA <Section_modify>`.

If a value begins with "v\_", a variable name must follow which has
been previously defined in the input script.  It must be a :doc:`grid-style variable <variable>`.  Such a variable defines a formula which can
reference stats keywords or invoke other computes, fixes, or variables
when they are evaluated.  So this is a very general means of creating
a per-grid input to perform an FFT on.


----------


If the *sum* keyword is set to *yes*\ , the results of all FFTs
will be summed together, grid value by grid value, to create
a single output.

The result of each FFT is scaled by the *sfactor* value of
the *scale* keyword, whose default is 1.0.  Note that forward FFTs do
not perform any scaling of their own; backward FFTs scale each output
value by N = # of points in the FFT grid.

If the *conjugate* keyword is set to *no*\ , the result of each FFT is 2
values for each grid point, the real and imaginary parts of a complex
number.  If the *conjugate* keyword is set to *yes*\ , the complex value
for each grid point is multiplied by its complex conjugate to yield a
single real-valued number for each grid point.  Note that this value
is effectively the squared length of the complex 2-vector with real
and imaginary components.

If one or more of the *kx*\ , *ky*\ , *kz*\ , or *kmag* keywords are set to
*yes*\ , then one or more extra columns of per-grid output is generated.
For *kx* the x-component of the K-space wavevector is generated.
Similarly for *ky* and *kz*\ .  For *kmag* the length of each K-space
wavevector is generated.  These values can be useful, for example, for
histogramming an energy spectrum computed from the FFT of a velocity
field, as a function of wavelength or a component of the wavelength.

Note that the wavevector for each grid cell is indexed as (Kx,Ky,Kz).
Those indices are the x,y,z components output by the *kx*\ , *ky*\ , *kz*
keywords.  The total wavelength, which is output by the *kmag*
keyword, is sqrt(Kx\^2 + Ky\^2 + Kz\^2) for 3d models and sqrt(Kx\^2 +
Ky\^2) for 2d models.  For all keywords, the Kx,Ky,Kz represent
distance from the origin in a periodic sense.  Thus for a grid that is
NxMxP, the Kx values associated with the x-dimension and used in those
formulas are not Kx = 0,1,2 ... N-2,N-1.  Rather they are Kx = 0,1,2,
... N/2-1, N/2, N/2-1, ... 2,1.  Similary for Ky in the y-dimension
with a max index of M/2, and Kz in the z-dimension with a max index of
P/2.





**Output info:**

The number of per-grid values ouptut by this compute depends on the
optional keyword settings.  The number of FFTs is equal to the number
of specified input values.

There are 2 columns of output per FFT if *sum* = no and *conjugate* =
no, with real and imaginary components for each FFT.  There is 1
column of output per FFT if *sum* = no and *conjugate* = yes.  There
are 2 columns of output if *sum* = yes and *conjugate* = no, with real
and imaginary components for the sum of all the FFTs.  There is one
column of output for *sum* = yes and *conjugate* = yes.  For all these
cases, there is one extra column of output for each of the *kx*\ , *ky*\ ,
*kz*\ , *kmag* keywords if they are set to *yes*\ .  The extra columns
come before the FFT columns, in the order *kx*\ , *ky*\ , *kz*\ , *kmag*\ .
Thus is only *ky* and *kmag* are set to yes, there will be 2 extra
columns, the first for *ky* and the 2nd for *kmag*\ .

If the total number of output columns = 1, then this compute produces
a per-grid vector as output.  Otherwise it produces a per-grid array.

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

The array can be accessed by any command that uses per-grid values
from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_
for an overview of SPARTA output options.

The per-grid vector or array values will be in the :doc:`units <units>`
appropriate to the FFT operations as described above.  The K-space
wavevector magnitudes are effectively unitless, e.g. sqrt(Kx\^2 + Ky\^2
+ Kz\^2) where Kx,Ky,Kz are integers.  The FFT values can be real or
imaginary or squared values in K-space resulting from FFTs of per-grid
quantities in whatever units the specified input values represent.


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

This style is part of the FFT package.  It is only enabled if SPARTA
was built with that package.  See the `Getting Started <Section_start.html#start_3>`_ section for more info.

**Related commands:**

:doc:`fix ave/grid <fix_ave_grid>`, :doc:`dump grid <dump>`, :doc:`compute grid <compute_grid>`

**Default:**

The option defaults are sum = no, scale = 1.0, conjugate = no, kmag =
no.


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
