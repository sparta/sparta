.. index:: react\_modify

react\_modify command
=====================

**Syntax:**


.. parsed-literal::

   react_modify keyword values ...

* one or more keyword/value pairs may be listed
* keywords = *recomb* or *rboost* or *compute\_chem\_rates* or *partial\_energy* or *reverse* or *keq\_file*
  
  .. parsed-literal::
  
       recomb value = yes or no = enable or disable defined recombination reactions
       rboost value = rfactor
         rfactor = boost probability of recombination reactions by this factor
       compute_chem_rates value = yes or no = enable or disable computation of Arrhenius rate for chemical
       reaction without performing the reaction
       partial_energy = yes or no = use partial energy or total energy for TCE chemistry
       reverse value = auto or no = auto-generate a detailed-balance reverse for every eligible forward reaction
       keq_file value = kfile or none
         kfile = file of equilibrium-constant curve fits for reverse reactions



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

Reverse (B-style) recombination reactions, derived by detailed balance
(see the :doc:`react <react>` command), have a probability that grows
without bound for soft collisions, so any boost moves part of their
rate into the cap at 1 and lowers the realized recombination rate
(by about 25% at 1e24 m^-3 and 6000 K with rboost = 1000).  Use
*rboost* = 1.0 with such reactions; SPARTA prints a warning when they
are active with *rboost* > 1.0.

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

With *partial\_energy* set to *no*\ , the TCE reaction energy is the
total collision energy, including the discrete vibrational and
electronic energy of the colliding pair, and the energy-dependent
factor of the TCE probability is the microcanonical
(density-of-states weighted) average of the standard factor over the
reactants' discrete energy ladders at fixed total collision energy:
the simple-harmonic levels of each vibrational mode (with :doc:`collide\_modify vibrate discrete <collide_modify>`) and the electronic states of
each reactant (with :doc:`collide\_modify electronic discrete <collide_modify>`), with translation, rotation, and smooth
vibration forming the continuous part.  The average is precomputed
once per reaction as a lookup table, so the runtime cost is one table
interpolation per potential reaction.

This construction keeps the equilibrium reaction rate on the input
Arrhenius rate for any mix of continuous and discrete internal modes
-- input rate sets (Park, Gupta, etc.) are equilibrium-calibrated, so
their coefficients already contain the equilibrium internal-state
contributions -- while out of equilibrium the reaction rate responds
to the actual internal-state populations in both directions: an
internally cold gas reacts more slowly than the thermal rate and an
internally hot gas reacts faster, with excited particles
preferentially selected, with no adjustable parameters.  Reactions
whose reactants carry no discrete ladders use the standard analytic
TCE factor, which is the same average with an empty ladder.  The
per-state DOF column of the electronic species file is read but not
used.


----------


The *reverse* keyword set to *auto* generates a B-style
(detailed-balance) reverse reaction for every eligible active forward
reaction in the reaction file read by the :doc:`react <react>` command,
as if a corresponding B line had been written by hand (see
:doc:`react <react>`): every exchange reaction gains the exchange
reverse with reactants and products swapped, every dissociation
reaction with an explicit third-body species gains the corresponding
recombination, and every 2-product (associative) ionization gains the
corresponding dissociative recombination as an exchange reverse.  A forward reaction is skipped if the file already
contains any reverse for it -- an explicit B line, an independently
fitted reverse rate, or a wildcard (atom/mol) recombination covering
the same product and third body (atom covers third bodies without
vibrational modes, mol those with) -- so existing reaction sets are
not double counted.  Generated recombinations obey the *recomb*
keyword like any other recombination.  Generated reactions appear in
the reaction tallies under their constructed formula, and are
appended after the reactions of the file, so :doc:`compute gas/reaction/grid <compute_gas_reaction_grid>` with *every* gains
a column for each of them.  Setting *reverse* back to *no* before a
later run makes the generated reactions inert.  *Reverse* *auto*
requires the :doc:`react tce <react>` style and *partial\_energy* *no*\ ;
it is an error otherwise.

The *keq\_file* keyword supplies equilibrium-constant curve fits that
replace the internally computed partition functions for the reverse
(B-style) reactions whose forward reaction they describe, so that a
published chemistry set's backward rates are reproduced exactly.  The
file contains 2-line entries (comments and blank lines are allowed):
the formula of the FORWARD reaction, then a fit line


.. parsed-literal::

   A + B --> C + D
   park c0 c1 c2 c3 c4

defining ln Keq = c0/Z + c1 + c2 ln(Z) + c3 Z + c4 Z^2 with Z = 10000
K / T (the form used by Park's chemistry reviews).  Keq = k\_f/k\_b must
be in SI units: dimensionless for an exchange pair, 1/m^3 for a
dissociation/recombination pair; published fits in molar or cgs units
must be converted.  The fit stays in SI under :doc:`units cgs <units>`:
SPARTA converts it to 1/cm^3 itself.  SPARTA warns if the fit and the
statistical-mechanics Keq differ by more than a factor of 10 at every
temperature of the fitting window (1000-60000 K), the signature of a fit
in the wrong units.  A reverse reaction matched by an entry keeps its
energy-resolved detailed-balance table (which by itself reproduces the
statistical-mechanics equilibrium constant Keq\_sm) and multiplies the
per-collision probability by the residual factor R(T) =
Keq\_sm(T)/Keq\_fit(T), evaluated at the representative cell
translational temperature.  Its thermal average therefore reproduces
the supplied Keq exactly, while the collision-energy selectivity stays
microscopically reversible and only the (usually small and smooth) Keq
discrepancy is applied thermally; ln R is fit once to the same Park
form so the tce and tce/kk styles evaluate it identically.  Entries
naming species that are not defined are skipped; entries matching no
reverse reaction produce a warning.  Setting *keq\_file* to *none*
restores the internal partition functions.


----------


**Restrictions:** none

**Related commands:**

:doc:`react <react>`

**Default:**

The option defaults are recomb = yes, rboost = 1000.0,
compute\_chem\_rates = no, partial\_energy = yes, reverse = no,
keq\_file = none.


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
