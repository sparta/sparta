.. index:: fix halt

fix halt command
================

**Syntax:**


.. parsed-literal::

   fix ID halt N attribute operator avalue keyword value ...

* ID is documented in :doc:`fix <fix>` command
* halt = style name of this fix command
* N = check halt condition every N timesteps
* attribute = *tlimit* or v\_name
  
  .. parsed-literal::
  
       tlimit = elapsed CPU time (in seconds)
       v_name = name of :doc:`equal-style variable <variable>`

* operator = "<" or "<=" or ">" or ">=" or "==" or "!=" or "\|\^"
* avalue = numeric value to compare attribute to
* zero or more keyword/value pairs may be appended
* keyword = *error* or *message*
  
  .. parsed-literal::
  
       error value = hard or soft or continue
       message value = yes or no



**Examples:**


.. parsed-literal::

   fix 1 halt 10 tlimit > 3600.0
   fix 1 halt 10 v_myCheck != 0 error soft

**Description:**

Check a condition every N steps during a simulation run. N must be
>=1. If the condition is met, exit the run.

The specified *attribute* can be one of the options listed above,
namely *tlimit*\ , or an :doc:`equal-style variable <variable>` referenced
as *v\_name*, where "name" is the name of a variable that has been
defined previously in the input script.

The *tlimit* attribute queries the elapsed CPU time (in seconds) since
the current run began, and sets *attribute* to that value. The clock
starts at the beginning of the current run (not when the fix command
is specified), so that any setup time for the run is not included in
the elapsed time. The timer invocation and syncing across all
processors (via MPI\_Allreduce) is performed (typically) only a small
number of times and the elapsed times are used to predict when the
end-of-the-run will be.  This can be useful when performing benchmark
calculations for a desired length of time with minimal overhead.

Equal-style variables evaluate to a numeric value. See the
:doc:`variable <variable>` command for a description. They calculate
formulas which can involve mathematical operations, particle
properties, grid properties, surface properties, global values
calculated by a :doc:`compute <compute>` or :doc:`fix <fix>`, or
references to other :doc:`variables <variable>`. Thus they are a very
general means of computing some attribute of the current system.  For
example, the following two versions of a fix halt command will both
stop the run after an hour of walltime:


.. parsed-literal::

   fix 1 halt 10 tlimit > 3600.0

   variable cpu equal cpu
   fix 1 halt 10 v_cpu > 3600.0

The commands above apply only to the time spent in the current run
command. If multiple run commands are used in the same input script,
one can also stop the run after a predetermined amount of *total*
walltime:


.. parsed-literal::

   variable wall equal wall
   fix 1 halt 10 v_wall > 3600.0

Similarly one can stop the run after a predetermined amount of
simulation time, which is useful when using a :doc:`variable timestep <fix_dt_reset>`:


.. parsed-literal::

   variable time equal time
   fix 1 halt 10 v_time > 1.0e-3

One can also stop the run if a file exists. The example below
checks for a file named "EXIT" in the current directory:


.. parsed-literal::

   variable file equal is_file(EXIT)
   fix 1 halt 10 v_file > 0

The choice of operators listed above are the usual comparison
operators. The XOR operation (exclusive or) is also included as "\|\^".
In this context, XOR means that if either the attribute or avalue is
0.0 and the other is non-zero, then the result is "true". Otherwise it
is "false".

The specified *avalue* must be a numeric value.

The optional *error* keyword determines how the current run is halted.
If its value is *hard*\ , then SPARTA will stop with an error message.

If its value is *soft*\ , SPARTA will exit the current run, but continue
to execute subsequent commands in the input script. However,
additional :doc:`run <run>` commands will be skipped. For example, this
allows a script to output the current state of the system, e.g. via a
:doc:`write\_grid <write_grid>` or :doc:`write\_restart <write_restart>`
command.

If its value is *continue*\ , the behavior is the same as for soft,
except subsequent run commands are executed. This allows your script
to remedy the condition that triggered the halt, if necessary.  Note
that you may wish use the :doc:`unfix <unfix>` command on the fix halt
ID, so that the same condition is not immediately triggered in a
subsequent run.


----------


The optional *message* keyword determines whether a message is printed
to the screen and logfile when the halt condition is triggered. If
*message* is set to *yes*\ , a one line message with the values that
triggered the halt is printed. If *message* is set to *no*\ , no message
is printed; the run simply exits. The latter may be desirable for
post-processing tools that extract statistical information from log
files.


----------


**Restart, output info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

This fix produces no output.


----------


**Restrictions:** none

**Related commands:**

:doc:`run <run>`

**Default:**

The option defaults are error = soft and message = yes.


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
