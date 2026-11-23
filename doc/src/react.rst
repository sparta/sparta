.. index:: react

react command
=============

**Syntax:**


.. parsed-literal::

   react style args

* style = *none* or *tce* or *qk* or *tce/qk*
* args = arguments for that style
  
  .. parsed-literal::
  
       none args = none
       tce args = infile
         infile = file with list of gas-phase chemistry reactions
       qk args = infile
         infile = file with list of gas-phase chemistry reactions
       tce/qk args = infile
         infile = file with list of gas-phase chemistry reactions
       tce/kk args = infile
         infile = file with list of gas-phase chemistry reactions
       qk/kk args = infile
         infile = file with list of gas-phase chemistry reactions
       tce/qk/kk args = infile
         infile = file with list of gas-phase chemistry reactions



**Examples:**


.. parsed-literal::

   react none
   react tce air.tce
   react qk air.tce

**Description:**

Define chemical reactions to perform in the gas phase when
particle-particle collisions occur.  See the
:doc:`surf\_react <surf_react>` command for specification of surface
chemistry reactions.

The *none* style means that no chemistry will be performed, which is
the default.

For other styles, a file is specified which contains a list of
chemical reactions, with their associated parameters.  The reactions
are read into SPARTA and stored in a list.  Each time a simulation is
run via the :doc:`run <run>` command, the list is scanned.  Only
reactions for which all the reactants and all the products are
currently defined as species-IDs will be active for the simulation.
Thus the file can contain more reactions than are used in a particular
simulation.  See the :doc:`species <species>` command for how species
IDs are defined.

The reaction models for the various styles are described below.  When
a pair of particles collide, the list of all reactions with those two
species as reactants is looped over.  A probability for each reaction
is calculated, using the formulas discussed below, and a random number
is used to decide which reaction (if any) takes place.

For the *tce* style, each reaction's probability is individually
capped at 1: the TCE derivation assumes the reaction cross-section is
a fraction of the VHS cross-section, so a computed probability above 1
is outside the model and is limited to the collision rate (following
Strand and Goldstein's total-cross-section correction as simplified by
Higdon).  When the summed probability of all of a pair's reactions
exceeds 1, the collision reacts with probability 1 and the reaction is
chosen proportionally to the individual probabilities, so the outcome
does not depend on the order in which reactions are listed in the
reaction file.  A warning is printed at the end of each run in which
either condition occurred, naming a reaction at which it did, since
the simulated rate then saturates at the collision rate and no longer
follows the input Arrhenius rate; reverse (B-style) reactions get a
warning of their own (see below).  The probability is per collision,
so the timestep, fnum and grid do not change it.

The format of the reaction file is the same for all three of the
currently defined styles, and is also described below.  The various
styles interpret and compute the specified reactions in different
ways.  The data directory in the SPARTA distribution contains reaction
files for these reaction models, all with the suffix ".tce".


----------


The *tce* style is Bird's Total Collision Energy (TCE) model.  When
this style is specified, all computed reactions will use the TCE
model.

Using kinetic theory, the TCE model allows for reaction probabilities
to be defined based on known, measured, reaction rates.  The model is
described in detail in `(Bird94) <#Bird94>`_; see chapter 6. The required
input parameters for each reaction (discussed below) are values that
permit its effective Arrhenius rate to be calculated, namely

.. math::

   K(T)  =  A T^b e^{-E_a/kT }

where K(T) is the forward reaction rate, T is the temperature of the
participating molecules which is a function of their velocities and
internal energy states, k the Boltzmann constant, and A,b,Ea are input
parameters as discussed below.

All 5 reactions coefficients read from the reaction file (described
below) are used to calculate terms in equation 6.10 of
`(Bird94) <#Bird94>`_ for the probability that a reaction takes place.

The C2, C3, C4 values are the Arrhenius activation energy Ea,
prefactor A, and exponent b, used in the rate formula above.

A negative activation energy is allowed, but the TCE model then
reproduces the Arrhenius rate only approximately.  The probability is
proportional to (Ec - Ea)\^(b-1+omega) (1 - Ea/Ec)\^(z+3/2-omega).  Its
thermal average over collision energies Ec > 0 misses the part of the
Arrhenius integral below Ec = 0, a fraction that grows with \|Ea\|/kT.
The probability also diverges as Ec approaches 0, so it is capped (see
below).

The TCE model can only reproduce the Arrhenius rate of a reaction if
its exponent b is greater than -(z+3/2), where z is the number of
internal degrees of freedom contributing to the collision energy:
coefficient C1 when using the default partial energy model (see the
:doc:`react\_modify <react_modify>` partial\_energy keyword), or the
average number of rotational (and, for smooth vibration, vibrational)
degrees of freedom of the reactants otherwise.  Below this bound the
gamma function in the TCE reaction probability is non-positive; SPARTA
clamps it so that the reaction still runs, but its rate is incorrect.
SPARTA prints one warning listing any such reactions at the start of
the first run.

Above the bound, the TCE reaction probability reproduces the Arrhenius
rate on average, but it is a per-collision probability and can exceed
1 for some collision energies, e.g. when the Arrhenius rate approaches
the collision rate, or near the activation energy when b < -(z+1/2).
It is then capped as described above, and the rates of the capped
reactions are under-predicted.

----------


The *qk* style is Bird's Quantum-Kinetic model (QK).  When this style
is specified, all computed reactions will use the QK model.

The QK model implemented is that of `(Bird09) <#Bird09>`_ as validated
`(Gallis09) <#Gallis09>`_ and modified `(Gallis10) <#Gallis10>`_.

The QK model depends solely on properties of the colliding molecules
and unlike the TCE model makes no use of measured reaction rates or 
adjustable parameters. The macroscopic properties used in the QK model 
are the available collision energy, activation energies, and quantized
vibrational energy levels.

According to the QK model dissociation reactions take place when the 
maximum obtainable vibrational energy after an inelastic energy exchange
is higher than the dissociation level `(Bird09) <#Bird09>`_.

.. math::

   int [E_c/(k \Theta_v)] > \Theta_d / \Theta_v

Exchange reactions take place when the vibrational energy after a trial 
energy exchange is above the activation energy of the exchange reaction
`(Gallis10) <#Gallis10>`_.

.. math::

   i_v > int[E_a/ (k \Theta_v)]

A new version of the QK model for exchange reactions has been proposed
by `(Bird11) <#Bird11>`_.  This will be implemented in future releases of
SPARTA.

For the QK model, SPARTA reads the same 5 coefficients per reaction
from the reaction file (described below) as for the TCE model.  Three
of the coefficients (C1,C2,C5) are used to calculate terms in equation
6.10 of `(Bird94) <#Bird94>`_ for the probability that a reaction takes
place.  The Arrhenius rate parameters C3 and C4 are ignored by the QK
model.


----------


The *tce/qk* style is a hybrid model which can be used to compute
reactions using both the TCE and QK models.  When this style is
specified, reactions from the input file that are flagged with an A =
Arrhenius style will be computed using the TCE model.  Reactions from
the input file that are flagged with a Q = Quantum style will be
computed using the QK model.


----------


The format of the input reaction file is as follows.  Comments or
blank lines are allowed in the file.  Comment lines start with a "#"
character.  All other entries must come in 2-line pairs with values
separated by whitespace in the following format


.. parsed-literal::

   R1 + R2 + ... --> P1 + P2 + ...
   type style C1 C2 ...

The first line is a text-based description of a single reaction.  R1,
R2, etc are one or more reactants, listed as :doc:`species <species>`
IDs.  P1, P2, etc are one or more products, also listed as
:doc:`species <species>` IDs.  The number of allowed reactants and
products depends on the reaction type, as discussed below.  In most
cases there is no restriction on the order or listed reactants or
products on what species are listed.  Exceptions are detailed below.
Note that individual reactants and products must be separated by
whitespace and a "+" sign.  The left-hand and right-hand sides of the
equation must be separated by whitespace and "-->".

The *type* of each reaction is a single character (upper or lower
case) with the following meaning.  The type determines how many
reactants and products can be specified in the first line.


.. parsed-literal::

   D = dissociation = 2 reactants and 3 products
   E = exchange = 2 reactants and 2 products
   I = ionization = 2 reactants and 2 or 3 products
   R = recombination = 2 reactants and 1 product (see below)

A dissociation reaction means that R1 dissociates into P1 and P2 when
it collides with R2.  R2 is preserved in the collision, so P3 = R2 is
required.

An exchange reaction is a collision between R1 and R2 that results in
new products P1 and P2.  There is no restriction on the species
involved in the reaction.

An ionization reaction with 2 products is typically a collision
between R1 and R2 that results in a positively charged ion and an
election.  See the discussion on ambipolar reactions below.  However,
SPARTA does not check for this, so there is no restriction on the
species involved in the reaction.

An ionization reaction with 3 products is typically a collision
between a neutral R1 and an electon R2 which ejects an electron from
the neutral species, resulting in an ion P1 and a new electron P2.
See the discussion on ambipolar reactions below.  Again, SPARTA does
not check for this, so there is no restriction on the species involved
in the reaction.  R2 is preserved in the collision, so P3 = R2 is
required.

A recombination reaction is a collision between R1 and R2 that results
in P1.  There is no restriction on the species involved in the
reaction.

Note that recombination reactions actually involve a 3rd particle
whose species is not altered by the reaction but whose velocity is, in
order to balance energy and momentum.  So conceptually it can be
thought of as both a reactant and a product.  There are 3 ways you can
specify recombination reactions, to include information about which
species of 3rd particles are eligible to participate:


.. parsed-literal::

   R1 + R2 -> P1
   R1 + R2 -> P1 + atom/mol
   R1 + R2 -> P1 + P2

In the first case, no info for a 3rd particle is listed.  This means
any species of 3rd particle can be used.  In the second case, a
non-species keyword is used, either "atom" or "mol".  This means the
3rd particle must be either an atomic species, or a molecular species.
This is based on the vibrational degrees of freedom listed in the
:doc:`species file <species>`.  A non-zero DOF is molecular; zero DOF is
atomic.  In the third case, a specific species P2 is listed.  This
means the 3rd particle must be that species.

Note that for the same R1 and R2, multiple recombination reactions can
be listed in the reaction file.  When two particles R1 and R2 are
selected for collision and a possible reaction, if any recombination
reaction is defined for R1 and R2, then a 3rd particle in the same
grid cell is randomly selected.  Its species P2 is used to match at
most one of the possibly multiple recombination reactions for R1 and
R2.  Only that recombination reaction is checked for a reaction as a
possible outcome of the collision.

This matching is done from most-specific to least-specific, i.e. the
reverse ordering of the 3 cases above.  If there is a defined reaction
that lists P2 (third case, most specific), it is used.  If not, and
there is a defined reaction for "atom" or "mol" that corresponds to P2
(second case, intermediate specificity), then it is used.  If not, and
there is a defined reaction with no P2 (first case, least specific),
then it is used.  If none of these matches occur, no recombination
reaction is possible for that collision between R1 and R2.  Note that
these matching rules means that for the same R1 and R2, you can list
two reactions, one with P2 = "atom", and one with P2 = "mol".  And/or
you can list multiple reactions of the third kind, each with a unique
P2.

IMPORTANT NOTE: If the ambipolar approximation is being used, via the
:doc:`fix ambipolar <fix_ambipolar>` and :doc:`collide\_modify ambipolar yes <collide_modify>` commands, then reactions which involve either
ambipolar ions or the ambipolar electron have more restricitve rules
about the ordering of reactants and products.  See the next section
for a discussion of these requirements.

The *style* of each reaction is a single character (upper or lower
case) with the following meaning:

* A = Arrhenius
* Q = Quantum
* B = Arrhenius Backward (see below)

The style determines how many reaction coefficients are listed as C1,
C2, etc, and how they are interpreted by SPARTA.

For both the A = Arrhenius style and Q = Quantum style, there are 5
coefficients:

* C1 = number of internal degrees of freedom (as defined by the TCE model)
* C2 = Arrhenius activation energy Ea
* C3 = Arrhenius prefactor A
* C4 = Arrhenius exponent b
* C5 = overall reaction energy (positive for exothermic)

The different reaction styles use these values in different ways, as
explained above.


----------


**Backward (detailed-balance) reactions**

The B = Arrhenius Backward style derives a reverse reaction rate from
its forward reaction using the principle of detailed balance, rather
than supplying an independent Arrhenius fit for the reverse rate.  See
`(Bird94) <#Bird94>`_ section 6.6 and Boyd & Schwartzentruber sections
7.5.2-7.5.3 for background.  A single Arrhenius form cannot exactly
reproduce k\_b(T) = k\_f(T)/K\_eq(T), because the equilibrium constant
K\_eq carries the non-Arrhenius temperature dependence of the
partition-function ratios, so independently-fitted forward and
backward rates do not satisfy detailed balance and the gas does not
relax to the correct chemical-equilibrium composition at elevated
temperatures.

Two kinds of pairs are supported.  A B-style exchange reaction pairs
with the forward A-style exchange reaction whose reactants/products
are its products/reactants:


.. parsed-literal::

   N2 + O --> NO + N
   E A 0.0 5.175e-19 1.069e-12 -1.0 -5.175e-19


.. parsed-literal::

   NO + N --> N2 + O
   E B 0.0 0.0 0.0 0.0 0.0

A B-style recombination reaction pairs with the forward dissociation
reaction with the same explicit third-body species (the atom/mol
wildcards cannot be used on a B recombination, since per-partner
forward dissociation rates differ):


.. parsed-literal::

   N2 + N --> N + N + N
   D A 1.0 1.561e-18 4.980e-8 -1.5 -1.561e-18


.. parsed-literal::

   N + N --> N2 + N
   R B 0.0 0.0 0.0 0.0 0.0

The 5 coefficients for a B reaction are placeholders; they are ignored
and overwritten internally.  The backward rate is constructed to
satisfy


.. parsed-literal::

   k_b(T) = k_f(T) / K_eq(T)

where K\_eq(T) is the equilibrium constant computed from the species
partition functions: translational, rotational (rigid rotor, with the
symmetry number from the :doc:`species <species>` *rotfile*\ ), vibrational
(harmonic oscillator, including degenerate modes), and electronic
(summed over the ladder from the *elecfile*\ , when defined).  These
include only the internal states the collision model carries: with
:doc:`collide\_modify <collide_modify>` *vibrate no* the vibrational
factor is 1, and with *electronic no* the electronic factor is the
ground-state degeneracy, since the detailed-balance tables can only
reproduce the equilibrium of the states the particles actually
occupy.  For the full equilibrium constants use *vibrate discrete*
(or *smooth*\ ) and *electronic discrete* with :doc:`fix elecmode <fix_elecmode>`.

Rotation is the exception: the tables always count the rotational
degrees of freedom, so reverse reactions need rotation modelled.  With
*rotate no* particles carry no rotational energy, the reverse rates of
rotating species are biased, and a recombination with a nonlinear
third body never occurs; SPARTA warns about this.  It also warns when
the *rotfile* gives no symmetry number for a rotating species used by
a reverse reaction, since the default sigma = 1 is wrong by a factor
of 2 in K\_eq for each homonuclear molecule.

Two properties of the data files are required for energy conservation
once reverse reactions are active, and SPARTA warns at initialization
when either fails.  First, every reaction must conserve mass: if the
product masses do not sum to the reactant masses (e.g. a species file
with m(O2) different from 2 m(O)), each reaction event changes the
centre-of-mass kinetic energy, and a gas that cycles through forward
and reverse reactions gains or loses energy steadily.  Second, the
reaction energies must be consistent around reaction cycles, such as
N2 + O -> NO + N, NO + M -> N + O + M, N + N + M -> N2 + M: if they do
not sum to zero around a cycle, every loop creates or destroys energy
and the reactions share no common equilibrium.  Both hold when all
reaction energies are computed from one set of species formation
energies, as in the *data/air.tce* file, whose header lists them.

This is enforced without evaluating any temperature.  At
initialization an exchange pair is given a detailed-balance table
relating the microcanonical densities of states of the two channels
at matched total collision energy (the same density-of-states
machinery that builds the microcanonical TCE energy factors), so each
collision's backward probability follows from the forward reaction by
energy-resolved microscopic reversibility.  The thermal average of
this probability reproduces k\_f(T)/K\_eq(T) at every temperature
simultaneously, and the per-collision probability remains meaningful
arbitrarily far from equilibrium, where no temperature is defined.

A dissociation/recombination pair is handled the same way, extended
to three bodies: the backward activation energy is seeded as Ea\_f +
dHrxn and the backward reaction energy as -dHrxn, and the
recombination probability is resolved in the total available energy
of all three particles - the collision energy of the recombining
pair, plus the third particle's internal energies, plus its
translational energy relative to the pair's center of mass (which is
exactly the relative translational energy of the forward
dissociation collision).  The probability is the ratio of the
forward channel's microcanonical reaction numerator at that total
energy to the density of states of the (pair energy, third-body
energy) decomposition, so the third particle's energy genuinely
counts toward the energy balance: an energetic third body enhances
recombination exactly as microscopic reversibility of the forward
dissociation demands, and the thermal average again reproduces
k\_f(T)/K\_eq(T) at every temperature with no cell temperature
evaluated.  Reverse reactions require
:doc:`react\_modify <react_modify>` *partial\_energy* *no*\ ; the default
is *yes*\ .  With *partial\_energy* *yes*\ , B lines in the reaction file
are switched off with a warning and the forward reactions run alone,
so a reaction file that carries reverses (such as data/air.tce) can
still be used with the default partial-energy model.

The cap on the forward probability described above is carried into the
reverse tables.  Where a forward reaction's probability is capped at
1, its rate falls below the Arrhenius fit, and a reverse reaction
built from the uncapped fit would outrun it and shift the equilibrium.
Each reverse table is therefore built from the forward probability
that is actually realized: the value capped at 1, divided by the
summed probability of the forward pair's reactions where that sum
exceeds 1.  A capped pair then reacts more slowly than its Arrhenius
fit in both directions, but relaxes to the correct equilibrium.  With
the *data/air.tce* set this matters for NO + O and NO + N
dissociation, whose rate fits exceed the VHS collision rate over much
of the thermal energy range.  In a test with those two reactions
alone, the equilibrium constant realized by the simulation was 0.57
of the partition-function value without the limit and 1.01 with it
(three seeds, 1.006-1.025), and 1.01 at a tenth of the density; an
adiabatic air reactor that should equilibrate at 6240 K settled 3.5%
too hot with 25% too much NO without the limit.  A test like this
depends on the post-reaction energy disposal as well: when a reaction
proceeds at the collision rate, its products react again before they
relax, so they must leave the reaction with the distribution that the
reverse samples (see below).  The limit works in one direction only: a
cap on the reverse reaction itself (three-body recombination in dense
gas, or a barrierless reverse at low collision energy) is not carried
back into its forward reaction.  Where it acts, the reverse rate falls
below its detailed-balance value and the equilibrium of the pair is
shifted toward the products of the forward reaction.  SPARTA prints a
separate warning at the end of each run in which a reverse reaction's
probability was capped.  This is mainly a low-temperature limit.  With
the *data/air.tce* set, NO + N -> N2 + O (the reverse of N2 + O -> NO +
N) reaches 0.74 of its detailed-balance rate at 1000 K and 0.90 at
2000 K.  The dissociative recombinations NO+ + e -> N + O, N2+ + e ->
N + N and O2+ + e -> O + O would need average probabilities of 11, 15
and 23 at 2000 K; they reach 0.54-0.69 of their detailed-balance rates
at 10000 K and 0.72-0.87 at 15000 K, so their ionization equilibria are
shifted toward the ions below about 15000 K.  (Rates of the reverse
reaction averaged over a thermal reverse pair, from
examples/tce\_air/tce\_ref.py.)

After a reaction the products' discrete levels (vibration, electronic
states, discrete rotation) are drawn jointly from the conditional, at
the collision's total energy, of the products' equilibrium collision
ensemble, and the remaining energy is split between translation and
the continuous modes; this is the distribution the reverse reaction
samples, so a reacting gas keeps equilibrium populations in every mode
even where relaxation is slow (vibration) or forbidden (spin-changing
electronic transitions).

For quantitatively correct equilibrium constants, define a *rotfile*
with the rotational symmetry numbers of homonuclear molecules and an
*elecfile* with at least the ground-state degeneracies of each species
(see the :doc:`species <species>` command): for N2 + O <-> NO + N the
ground-state electronic degeneracies alone change K\_eq by a factor of
~1.8 and the N2 symmetry number by a further factor of 2.  A
*rotfile* entry is required for every rotating species that takes part
in an active reverse reaction; it is an error if one is missing.  See
examples/reverse for a complete worked example and
examples/reverse/validate\_reverse.py for a quantitative validation
against independently computed equilibrium constants.  Instead of
writing B lines by hand, :doc:`react\_modify <react_modify>` *reverse*
*auto* generates the reverse of every eligible forward reaction
automatically, and :doc:`react\_modify <react_modify>` *keq\_file* can
supply equilibrium-constant curve fits that override the internal
partition functions when a published chemistry set's backward rates
must be reproduced exactly.

The B style requires the :doc:`react tce <react>` style (with or without
:doc:`react\_modify <react_modify>` *compute\_chem\_rates*\ ); with the qk
and tce/qk styles, B lines are switched off with a warning.  A reverse
reaction is also switched off if the exponent b of its forward
reaction is below the -(z+3/2) bound described above, since the
forward probability then does not reproduce its Arrhenius rate; the
warning about the forward reaction says so.  Reactions between charged and neutral
heavy species are fully supported: a charge-exchange reaction such as
A+ + B --> A + B+ is an exchange (E) reaction and pairs with a B-style
reverse like any other exchange, since no free electron participates
and all partition functions are those of heavy particles.  The reverse
of a 2-product ionization (I) reaction - associative ionization A + B
--> AB+ + e, whose reverse is dissociative recombination AB+ + e --> A
+ B - is also supported: write the reverse as an E-style B line (it is
a 2 -> 2 reaction).  Its Saha-type equilibrium constant is built from
the same partition functions, with the free electron contributing its
translational partition function and spin degeneracy 2; the reverse
probability is resolved in the electron-ion collision energy, so with
:doc:`fix ambipolar <fix_ambipolar>` the recombination rate follows the
ambipolar electron velocities (the electron energies maintained by the
elastic-collision exchange) rather than the heavy-particle cell
temperature.  What remains unsupported is a B-style reverse of a
3-product electron-impact ionization (A + e --> A+ + e + e): its
reverse is a 3-body process in the electron gas and must be supplied
explicitly.  The KOKKOS-accelerated
tce/kk style implements the identical path and reproduces the CPU code
bit-for-bit when built with -DSPARTA\_KOKKOS\_EXACT.


----------


If the ambipolar approximation is being used, via the :doc:`fix ambipolar <fix_ambipolar>` command, then reactions which involve
either ambipolar ions or the ambipolar electron have more restricitve
rules about the ordering of reactants and products, than those
described in the preceeding section.

Note that ambipolar collisions are turned on via the :doc:`collide\_modify ambipolar yes <collide_modify>` commands, which in turn requries
that the :doc:`fix ambipolar <fix_ambipolar>` is defined in your input
script.  This fix defines a particular species as an ambipolar
electron, written as "e" in the reactions that follow.  It also
defines a list of ambipolar ions, which are written as species with a
trailing "+" sign in the rules that follow.  Neutral species (without
"+") can be any non-ambipolar species.

These rules only apply to reactions that involve ambipolar species
(ions or electrons) as a reactant or product.  Note that every
ambipolar reaction written here conserves charge.  I.e. the net charge
of the reactants equals the net charge of the products.

Ambipolar dissociation reactions must list their reactants and
products in one of the following orders:


.. parsed-literal::

   AB + e -> A + e + B
   AB+ + e -> A+ + e + B
   AB + C+ -> A + C+ + B
   AB+ + C+ -> A+ + C+ + B

In the last two forms an ambipolar ion C+ acts as an inert third
body (collision partner).  As with the electron "e", the third body
must be written as the middle product, so that it maps to the
unaltered second colliding particle.

Ambipolar ionization reactions with 2 or 3 products must be in one of
the following orders:


.. parsed-literal::

   A + B -> AB+ + e
   A + e -> A+ + e + e

Ambipolar exchange reactions must be in one of the following orders:


.. parsed-literal::

   AB+ + e -> A + B
   AB+ + C -> A + BC+
   C + AB+ -> A + BC+

Ambipolar recombination reactions must be in the following order:


.. parsed-literal::

   A+ + e -> A
   A + B+ -> AB+
   A+ + B -> AB+
   A + B -> AB + C+

In the last form an ambipolar ion C+ acts as the inert third body,
written as the second product P2 (see the non-ambipolar discussion of
3rd particles above).  The two reactants and the recombined molecule
P1 must be neutral; the ion C+ is a spectator whose species selects the
recombination rate.

A third particle for recombination reactions can be specified in
the same way as described above for non-ambipolar recombination.


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

:doc:`collide <collide>`, :doc:`surf\_react <surf_react>`

**Default:**

style = none


----------


.. _Bird94:

.. raw:: html

   <span id="Bird94"></span>



**(Bird94)** G. A. Bird, Molecular Gas Dynamics and the Direct
Simulation of Gas Flows, Clarendon Press, Oxford (1994).

.. _Bird09:

.. raw:: html

   <span id="Bird09"></span>



**(Bird09)** G. A. Bird, Chemical Reactions in DSMC Rarefied Gas
Dynamics, Editor T Abe, AIP Conference Proceedings (2009).

.. _Bird11:

.. raw:: html

   <span id="Bird11"></span>



**(Bird11)** G. A. Bird, "The Q-K model for gas-phase chemical reaction
rates", Physics of Fluids, 23, 106101, (2011).

.. _Gallis09:

.. raw:: html

   <span id="Gallis09"></span>



**(Gallis09)** M. A. Gallis, R. B. Bond, and J. R. Torczynski, "A
Kinetic-Theory Approach for Computing Chemical-Reaction Rates in
Upper-Atmosphere Hypersonic Flows", J Chem Phys, 131, 124311, (2009).

.. _Gallis10:

.. raw:: html

   <span id="Gallis10"></span>



**(Gallis10)** M. A. Gallis, R. B. Bond, and J. R.Torczynski,
"Assessment of Collision-Energy-Based Models for Atmospheric-Species
Reactions in Hypersonic Flows", J Thermophysics and Heat Transfer,
(2010).


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
