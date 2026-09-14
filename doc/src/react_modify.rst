.. index:: react\_modify

react\_modify command
=====================

**Syntax:**


.. parsed-literal::

   react_modify keyword values ...

* one or more keyword/value pairs may be listed
* keywords = *recomb* or *rboost* or *compute\_chem\_rates* or *partial\_energy*
  
  .. parsed-literal::
  
       recomb value = yes or no = enable or disable defined recombination reactions
       rboost value = rfactor
         rfactor = boost probability of recombination reactions by this factor
       compute_chem_rates value = yes or no = enable or disable computation of Arrhenius rate for chemical
       reaction without performing the reaction
       partial_energy = yes or no = use partial energy or total energy for TCE chemistry



**Examples:**


.. parsed-literal::

   react_modify recomb no
   react_modify rboost 100.0

**Description:**

Set parameters that affect how reactions are performed.

The *recomb* keyword turns on or off recombination reactions.  It is
only relevant if recombination reactions were defined in the reaction
file read in by the :doc:`react <react>` command.  If the setting is
*no* then they will be disabled even if they were listed in the
reaction file.  This is useful to turn recombination reactions off, to
see if they affect simulation results.

The *rboost* keyword is a setting for recombination reactions.  It is
ignored if no recombination reactions exist, or the *recomb* keyword
is set to *no*\ .  The *rboost* setting does not affect the overalll
statistical results of recombination reactions, but tries to improve
their computational efficiency.  Recombination reactions typically
occur with very low probability, which means the code spends time
testing for reactions that rarely occur.  If the *rfactor* is set to N
> 1, then recombination reactions are skipped N-1 out of N times, when
one or more such reactions is defined for a pair of colliding
particles.  A random number us used to select on that probability.  To
compensate, when a recombination reaction is actually tested for
occurrence, its rate is boosted by a factor of N, making it N times
more likely to occur.

The smallest value *rboost* can be set to is 1.0, which effectively
applies no boost factor.

IMPORTANT NOTE: Setting *rboost* too large could meant the probability
of a recombination reaction becomes > 1.0, when it is does occur.
SPARTA does not check for this, so you should estimate the largest
boost factor that is safe to use for your model.

The *compute\_chem\_rates* keyword is a setting that allows the user to
only compute Arrhenius rates for chemical reactions without performing them.
Currently only the TCE reaction model supports this keyword; an error
will occur when using the QK or TCE/QK reaction model with this keyword.

The *partial\_energy* keyword is a setting that allows the user to
choose the amount of internal energy and internal degrees of freedom
used in the TCE model.

If the *partial\_energy* keyword is set to *yes*\ , the rDOF model of
Bird is used, and only the sum of the relative translational energy
between the partcles and a fraction of the rotational energy is
used. The participating internal degrees of freedom are either set to
1 (dissociation reactions), or 0 (recombination, exchange, ionization
reactions).

Conversely, if the *partial\_energy* keyword is set to *no*\ , then the
total energy model is used, i.e. the sum of the relative translational
energy between the partcles and the rotational and vibrational
energies. The participating internal degrees of freedom are computed
directly by the code and do not need to be inputted by the user. The
vibrational energy model used has an impact on the internal degrees of
freedom used in the TCE model in that case. This option is ignored for
the QK reaction model.


----------


**Restrictions:** none

**Related commands:**

:doc:`react <react>`

**Default:**

The option defaults are recomb = yes, rboost = 1000.0,
compute\_chem\_rates = no, partial\_energy = yes.


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
