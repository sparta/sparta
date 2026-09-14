.. index:: fix ave/time

fix ave/time command
====================

**Syntax:**


.. parsed-literal::

   fix ID ave/time Nevery Nrepeat Nfreq value1 value2 ... keyword args ...

* ID is documented in :doc:`fix <fix>` command
* ave/time = style name of this fix command
* Nevery = use input values every this many timesteps
* Nrepeat = # of times to use input values for calculating averages
* Nfreq = calculate averages every this many timesteps
* one or more input values can be listed
* value = c\_ID, c\_ID[N], f\_ID, f\_ID[N], v\_name
  
  .. parsed-literal::
  
       c_ID = global scalar or vector or array calculated by a compute with ID
       c_ID[I] = Ith component of global vector or Ith column of global array calculated by a compute with ID, I can include wildcard (see below)
       f_ID = global scalar or vector or array calculated by a fix with ID
       f_ID[I] = Ith component of global vector or Ith column of global array calculated by a fix with ID, I can include wildcard (see below)
       v_name = global value calculated by an equal-style variable with name

* zero or more keyword/arg pairs may be appended
* keyword = *mode* or *file* or *ave* or *start* or *off* or *title1* or *title2* or *title3*
  
  .. parsed-literal::
  
       mode arg = scalar or vector
         scalar = all input values are global scalars
         vector = all input values are global vectors or global arrays
       ave args = one or running or window M
         one = output a new average value every Nfreq steps
         running = output cummulative average of all previous Nfreq steps
         window M = output average of M most recent Nfreq steps
       start args = Nstart
         Nstart = start averaging on this timestep
       off arg = M = do not average this value
         M = value # from 1 to Nvalues
       file arg = filename
         filename = name of file to output time averages to
       title1 arg = string
         string = text to print as 1st line of output file
       title2 arg = string
         string = text to print as 2nd line of output file
       title3 arg = string
         string = text to print as 3rd line of output file, only for vector mode



**Examples:**


.. parsed-literal::

   fix 1 ave/time 100 5 1000 c_myTemp c_thermo_temp file temp.profile
   fix 1 ave/time 100 5 1000 c_myCount[2] c_myCount[3] ave window 20 &
                                 title1 "My output values"
   fix 1 ave/time 100 5 1000 c_myCount[\*] ave window 20
   fix 1 ave/time 1 100 1000 f_indent f_indent[1] file temp.indent off 1

**Description:**

Use one or more global values as inputs every few timesteps, and
average them over longer timescales.  The resulting averages can be
used by other output commands such as :doc:`stats\_style custom <stats_style>`, and can also be written to a file.  Note
that if no time averaging is done, this command can be used as a
convenient way to simply output one or more global values to a file.

Each listed value can be the result of a :doc:`compute <compute>` or
:doc:`fix <fix>` or the evaluation of an equal-style
:doc:`variable <variable>`.  In each case, the compute, fix, or variable
must produce a global quantity, not a per-grid or per-surf quantity.
If you wish to time-average those quantities, see the :doc:`fix ave/grid <fix_ave_grid>` and :doc:`fix ave/surf <fix_ave_surf>`
commands.

:doc:`Computes <compute>` that produce global quantities are those which
do not have the word *particle* or *grid* or *surf* in their style
name.  Only a few :doc:`fixes <fix>` produce global quantities.  See the
doc pages for individual fixes for info on which ones produce such
values.  :doc:`Variables <variable>` of style *equal* are the only ones
that can be used with this fix.  Variables of style *particle* cannot
be used, since they produce per-particle values.

The input values must either be all scalars or all vectors (or
arrays), depending on the setting of the *mode* keyword.  In both
cases, the averaging is performed independently on each input value.
I.e. each input scalar is averaged independently and each element of
each input vector (or array) is averaged independently.

If *mode* = scalar, then the input values must be scalars, or vectors
with a bracketed term appended, indicating the Ith value of the vector
is used.

If *mode* = vector, then the input values must be vectors, or arrays
with a bracketed term appended, indicating the Ith column of the array
is used.  All vectors must be the same length, which is the length of
the vector or number of rows in the array.

Note that for values from a compute or fix, the bracketed index I can
be specified using a wildcard asterisk with the index to effectively
specify multiple values.  This takes the form "\*" or "\*n" or "n\*" or
"m\*n".  If N = the size of the vector (for *mode* = scalar) or the
number of columns in the array (for *mode* = vector), then an asterisk
with no numeric values means all indices from 1 to N.  A leading
asterisk means all indices from 1 to n (inclusive).  A trailing
asterisk means all indices from n to N (inclusive).  A middle asterisk
means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual elements of the
vector or columns of the array had been listed one by one.  E.g. these
2 fix ave/time commands are equivalent, since the :doc:`compute count <compute_count>` command creates, in this case, a global
vector with 3 values.


.. parsed-literal::

   compute 1 count Ar He O
   fix 1 ave/time 100 1 100 c_1 file tmp.count
   fix 1 ave/time 100 1 100 c_1[1] c_1[2] c_1[3] file tmp.count


----------


The *Nevery*\ , *Nrepeat*\ , and *Nfreq* arguments specify on what
timesteps the input values will be used in order to contribute to the
average.  The final averaged quantities are generated on timesteps
that are a mlutiple of *Nfreq*\ .  The average is over *Nrepeat*
quantities, computed in the preceding portion of the simulation every
*Nevery* timesteps.  *Nfreq* must be a multiple of *Nevery* and
*Nevery* must be non-zero even if *Nrepeat* is 1.  Also, the timesteps
contributing to the average value cannot overlap, i.e. Nfreq >
(Nrepeat-1)\*Nevery is required.

For example, if Nevery=2, Nrepeat=6, and Nfreq=100, then values on
timesteps 90,92,94,96,98,100 will be used to compute the final average
on timestep 100.  Similarly for timesteps 190,192,194,196,198,200 on
timestep 200, etc.  If Nrepeat=1 and Nfreq = 100, then no time
averaging is done; values are simply generated on timesteps
100,200,etc.


----------


If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  If *mode* = scalar, then if
no bracketed term is appended, the global scalar calculated by the
compute is used.  If a bracketed term is appended, the Ith element of
the global vector calculated by the compute is used.  If *mode* =
vector, then if no bracketed term is appended, the global vector
calculated by the compute is used.  If a bracketed term is appended,
the Ith column of the global array calculated by the compute is used.
See the discussion above for how I can be specified with a wildcard
asterisk to effectively specify multiple values.

Note that there is a :doc:`compute reduce <compute_reduce>` command
which can sum per-particle or per-grid or per-surf quantities into a
global scalar or vector which can thus be accessed by fix ave/time.
Also Note that users can also write code for their own compute styles
and :doc:`add them to SPARTA <Section_modify>`; their output can then be
processed by this fix.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  If *mode* = scalar, then if
no bracketed term is appended, the global scalar calculated by the fix
is used.  If a bracketed term is appended, the Ith element of the
global vector calculated by the fix is used.  If *mode* = vector, then
if no bracketed term is appended, the global vector calculated by the
fix is used.  If a bracketed term is appended, the Ith column of the
global array calculated by the fix is used.  See the discussion above
for how I can be specified with a wildcard asterisk to effectively
specify multiple values.

Note that some fixes only produce their values on certain timesteps,
which must be compatible with *Nevery*\ , else an error will result.
Users can also write code for their own fix styles and :doc:`add them to SPARTA <Section_modify>`.

If a value begins with "v\_", a variable name must follow which has
been previously defined in the input script.  Variables can only be
used as input for *mode* = scalar.  Only equal-style variables can be
referenced.  See the :doc:`variable <variable>` command for details.
Note that variables of style *equal* define a formula which can
reference :doc:`stats\_style <stats_style>` keywords, or they can invoke
other computes, fixes, or variables when they are evaluated, so this
is a very general means of specifying quantities to time average.


----------


Additional optional keywords also affect the operation of this fix.

If the *mode* keyword is set to *scalar*\ , then all input values must
be global scalars, or elements of global vectors.  If the *mode*
keyword is set to *vector*\ , then all input values must be global
vectors, or columns of global arrays.  They can also be global arrays,
which are converted into a series of global vectors (one per column),
as explained above.

The *ave* keyword determines how the values produced every *Nfreq*
steps are averaged with values produced on previous steps that were
multiples of *Nfreq*\ , before they are accessed by another output
command or written to a file.

If the *ave* setting is *one*\ , then the values produced on timesteps
that are multiples of *Nfreq* are independent of each other; they are
output as-is without further averaging.

If the *ave* setting is *running*\ , then the values produced on
timesteps that are multiples of *Nfreq* are summed and averaged in a
cummulative sense before being output.  Each output value is thus the
average of the value produced on that timestep with all preceding
values.  This running average begins when the fix is defined; it can
only be restarted by deleting the fix via the :doc:`unfix <unfix>`
command, or by re-defining the fix by re-specifying it.

If the *ave* setting is *window*\ , then the values produced on
timesteps that are multiples of *Nfreq* are summed and averaged within
a moving "window" of time, so that the last M values are used to
produce the output.  E.g. if M = 3 and Nfreq = 1000, then the output
on step 10000 will be the average of the individual values on steps
8000,9000,10000.  Outputs on early steps will average over less than M
values if they are not available.

The *start* keyword specifies what timestep averaging will begin on.
The default is step 0.  Often input values can be 0.0 at time 0, so
setting *start* to a larger value can avoid including a 0.0 in a
running or windowed average.

The *off* keyword can be used to flag any of the input values.  If a
value is flagged, it will not be time averaged.  Instead the most
recent input value will always be stored and output.  This is useful
if one of more of the inputs produced by a compute or fix or variable
are effectively constant or are simply current values.  E.g. they are
being written to a file with other time-averaged values for purposes
of creating well-formatted output.

The *file* keyword allows a filename to be specified.  Every *Nfreq*
steps, one quantity or vector of quantities is written to the file for
each input value specified in the fix ave/time command.  For *mode* =
scalar, this means a single line is written each time output is
performed.  Thus the file ends up to be a series of lines, i.e. one
column of numbers for each input value.  For *mode* = vector, an array
of numbers is written each time output is performed.  The number of
rows is the length of the input vectors, and the number of columns is
the number of values.  Thus the file ends up to be a series of these
array sections.

The *title1* and *title2* and *title3* keywords allow specification of
the strings that will be printed as the first 2 or 3 lines of the
output file, assuming the *file* keyword was used.  SPARTA uses
default values for each of these, so they do not need to be specified.

By default, these header lines are as follows for *mode* = scalar:


.. parsed-literal::

   # Time-averaged data for fix ID
   # TimeStep value1 value2 ...

In the first line, ID is replaced with the fix-ID.  In the second line
the values are replaced with the appropriate fields from the fix
ave/time command.  There is no third line in the header of the file,
so the *title3* setting is ignored when *mode* = scalar.

By default, these header lines are as follows for *mode* = vector:


.. parsed-literal::

   # Time-averaged data for fix ID
   # TimeStep Number-of-rows
   # Row value1 value2 ...

In the first line, ID is replaced with the fix-ID.  The second line
describes the two values that are printed at the first of each section
of output.  In the third line the values are replaced with the
appropriate fields from the fix ave/time command.


----------


**Restart, output info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

This fix produces a global scalar or global vector or global array
which can be accessed by various output commands.  The values can only
be accessed on timesteps that are multiples of *Nfreq* since that is
when averaging is performed.

A scalar is produced if only a single input value is averaged and
*mode* = scalar.  A vector is produced if multiple input values are
averaged for *mode* = scalar, or a single input value for *mode* =
vector.  In the first case, the length of the vector is the number of
inputs.  In the second case, the length of the vector is the same as
the length of the input vector.  An array is produced if multiple
input values are averaged and *mode* = vector.  The global array has #
of rows = length of the input vectors and # of columns = number of
inputs.

**Restrictions:** none

**Related commands:**

:doc:`compute <compute>`, :doc:`fix ave/grid <fix_ave_grid>`, :doc:`fix ave/surf <fix_ave_surf>`, :doc:`variable <variable>`

**Default:**

The option defaults are mode = scalar, ave = one, start = 0, no file
output, title 1,2,3 = strings as described above, and no off settings
for any input values.


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
