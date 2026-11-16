.. index:: fix controller

fix controller command
======================

**Syntax:**


.. parsed-literal::

   fix ID controller Nevery alpha Kp Ki Kd pv setpoint cv

* ID is documented in :doc:`fix <fix>` command
* controller = style name of this fix command
* Nevery = invoke this fix every this many timesteps
* alpha = coupling constant for converting the control variable to a change (see below)
* Kp = proportional gain coefficient (see below)
* Ki = integral gain coefficient (see below)
* Kd = derivative gain coefficient (see below)
* pv = process variable of the feedback loop, can be a compute, fix, or variable
  
  .. parsed-literal::
  
       c_ID = global scalar calculated by a compute with ID
       c_ID[I] = Ith component of global vector calculated by a compute with ID
       f_ID = global scalar calculated by a fix with ID
       f_ID[I] = Ith component of global vector calculated by a fix with ID
       v_name = value calculated by an equal-style variable with name

* setpoint = desired value of the process variable (target of the controller)
* cv = name of the control variable, an internal-style :doc:`variable <variable>`


**Examples:**


.. parsed-literal::

   compute myTemp all temp
   variable drive internal 0.0
   fix pid controller 100 0.1 1.0 0.0 0.0 c_myTemp 300.0 drive

**Description:**

Apply a "proportional-integral-derivative" (PID) feedback controller to
a simulation.  Every *Nevery* timesteps the fix samples a *process
variable* *pv*\ , compares it to a desired *setpoint*\ , and updates a
*control variable* *cv* so as to drive the process variable toward the
setpoint.  The control variable is an :doc:`internal-style variable <variable>` whose value can be referenced by any other
command in the input script that evaluates a variable, so the feedback
loop can be used to steer essentially any property of a running
simulation.

A typical use is to maintain a target temperature, pressure, or flux by
adjusting an input such as a wall temperature, an inflow velocity, or an
applied field that is itself expressed through the control variable.
The process variable to be controlled is computed by another command and
passed to this fix as *pv*\ ; the value that this fix adjusts is the
internal-style variable named *cv*\ .

The *process variable* *pv* can be the output of a :doc:`compute <compute>`,
a :doc:`fix <fix>`, or an equal-style :doc:`variable <variable>`:

* if *pv* begins with "c\_", it is the global scalar (no bracket) or the
  Ith component of the global vector (with bracket) computed by the named
  compute
* if *pv* begins with "f\_", it is the global scalar or Ith vector
  component computed by the named fix
* if *pv* begins with "v\_", it is the value of the named equal-style
  variable


The *control variable* *cv* must be the name of an :doc:`internal-style variable <variable>`.  This fix overwrites the value of that variable
each time it is invoked.  Some other command in the input script should
reference this variable (e.g. as *v\_name*) so that the new value
produced by the controller takes effect.

The controller works as follows.  Each time the fix is invoked, the
current value of the process variable *pv* is compared to the
*setpoint*\ , and the difference is the current error:


.. parsed-literal::

   err = pv - setpoint

The control variable is then updated from its previous value using the
three PID terms:


.. parsed-literal::

   cv = cv - Kp\*alpha\*tau\*err - Ki\*alpha\*tau\*tau\*sumerr - Kd\*alpha\*deltaerr

where *tau* = *Nevery* \* dt is the elapsed time between invocations (dt
is the timestep size), *sumerr* is the running sum (time integral) of
the error, and *deltaerr* is the change in the error since the previous
invocation (a finite-difference approximation to the time derivative).

The *Kp*\ , *Ki*\ , and *Kd* coefficients are the dimensionless gains of
the proportional, integral, and derivative terms, respectively.  The
*alpha* coefficient is an overall coupling constant that carries the
units of the equation, namely the control-variable units divided by the
process-variable units divided by time units.  The advantage of this
convention is that the values of *Kp*\ , *Ki*\ , *Kd* are invariant under a
change of the timestep
size or of *Nevery*\ ; and if the :doc:`units <units>` style is changed,
only *alpha* needs to be adjusted, leaving the three gains unaltered.
Setting one or more of *Kp*\ , *Ki*\ , *Kd* to zero disables the
corresponding term, e.g. *Ki* = *Kd* = 0.0 yields a purely proportional
controller.

NOTE: The sign convention in the update formula above (the PID terms are
subtracted) is appropriate when an increase of the control variable
*decreases* the process variable.  If increasing the control variable
*increases* the process variable, use a negative *alpha* (or negative
gains) so that the feedback remains stabilizing.

Choosing good values for the four constants depends on the system and
typically requires some experimentation.  It is best to first choose a
value and sign for *alpha* consistent with the magnitudes and signs of
the process and control variables, then tune *Kp* (with *Ki* = *Kd* =
0.0) for a fast response that does not overshoot the setpoint.  For many
applications a purely proportional controller is sufficient.  A non-zero
*Ki* can remove a steady-state offset when the process variable plateaus
before reaching the setpoint, and a non-zero *Kd* can counteract lag in
the response of the process variable to a change in the control variable.

Because this fix updates the control variable but does not initialize
it, its initial value is whatever the user assigned to the internal-style
:doc:`variable <variable>` in the input script.  That value is used (by
every other command that references the variable) until this fix
performs its first update after *Nevery* timesteps.  On that first
update the derivative term is set to zero, because the previous error is
not yet defined.


----------


**Restart, output info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

This fix computes a global vector of length 3 which can be accessed by
various output commands.  The values can be accessed on any timestep,
though they are only updated on timesteps that are a multiple of
*Nevery*\ .  The 3 quantities are the most recent updates made to the
control variable by each of the three terms in the PID equation above:

* 1 = proportional term = -Kp\*alpha\*tau\*err
* 2 = integral term = -Ki\*alpha\*tau\*tau\*sumerr
* 3 = derivative term = -Kd\*alpha\*deltaerr

These values can be useful for monitoring and tuning the relative
contributions of the three terms.  The units of the vector values are
whatever units the control variable is in.  The vector values calculated
by this fix are "extensive".


----------


**Restrictions:** none

**Related commands:**

:doc:`compute <compute>`, :doc:`variable <variable>`, :doc:`fix adapt <fix_adapt>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
