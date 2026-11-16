.. index:: species

species command
===============

**Syntax:**


.. parsed-literal::

   species file ID1 ID2 ... keyword value ...

* file = filename with species info
* ID1, ID2, ... = one or more species names listed in file
* multi-species abbreviations can also be used (see below)
* zero or more keyword/value pairs may be appended
* keyword = *rotfile*
  
  .. parsed-literal::
  
       rotfile value = rfile = filename for extra rotational info

* keyword = *vibfile*
  
  .. parsed-literal::
  
       vibfile value = vfile = filename for extra vibrational info

* keyword = *elecfile*
  
  .. parsed-literal::
  
       elecfile value = efile = filename for extra electronic info



**Examples:**


.. parsed-literal::

   species air.species air
   species ar.species Ar
   species air.species air CO2 CO vibfile co2.species.vib
   species myfile H+ Cl- HCl

**Description:**

Define one or more particle species to use in the simulation.  This
command can be used as many times as desired to add species to the
list of species that the simulation recognizes.

The specified *file* is the name of a file containing definitions for
a list of species, not all of which need to specified in this command,
or used in a simulation.  Only those requested by ID will be extracted
from the file and they must be present in the file.  The format of the
species file is discussed below.  The data directory in the SPARTA
distribution contains several species files, all with the suffix
".species".

Each *ID* is a character string used to identify the species, such as
N or O2 or NO or D or Fe-.  The string can be any combination of
alphanumeric characters or "+", "-", or underscore.

Instead of specifying IDs for single species, one of several
pre-defined multi-species names can be used, each of which is expanded
into a list of several individual species IDs.  The list of currently
recognized abbreviations is as follows:

* air = N, O, NO

These abbreviations can be used in combination with single-species IDs
as in the 3rd example above.


----------


The format of a species file is as follows.  Comments or blank lines
are allowed in the file.  Comment lines start with a "#" character.
All other lines must have the following format with values separated
by whitespace:


.. parsed-literal::

   species-ID prop1 prop2 ... prop9 prop10

The species-ID is a string that will be matched against the requested
species-ID, as described above.  The properties are as follows:

* prop1 = molecular weight (atomic mass units, e.g. 16 for oxygen)
* prop2 = molecular mass (mass units)
* prop3 = rotational degrees of freedom (integer, unitless)
* prop4 = inverse rotational relaxtion number (unitless)
* prop5 = vibrational degrees of freedom (integer, unitless)
* prop6 = inverse vibrational relaxation number (unitless)
* prop7 = vibrational temperature (temperature units)
* prop8 = species weight (unitless)
* prop9 = multiple of electon charge (1 for a proton)

The allowed values for rotational degrees of freedom (rotdof = prop3)
are 0,2,3.  Typically, 0 = monatomic species, 2 = diatomic, 3 =
anything else.

The allowed values for vibrational degrees of freedom (vibdof = prop5)
are 0,2,4,6,8.  The associated number of vibrational modes is vibdof
divided by 2.  Typically, 0 modes = monatomic species, 1 mode =
diatomic, 2/3/4 modes = triatomic or larger molecules.

Note that all the listed rotational and vibrational values must be
specified for each species, but in cases where they are not used by
SPARTA, they can simply be specified as 0.0.  Whether or not the
values are used for a species depends on the value of rotdof and
vibdof.  Whether the values are used in a simulation also depends on
the settings specified for the *rotation* and *vibration* keywords of
the :doc:`collide\_modify <collide_modify>` command.

Specifically, if prop3 for rotdof = 0, then prop4 is ignored.  If
prop5 for vibdof = 0, then prop6 and prop7 are ignored.

If vibdof = 4,6,8, then information for 2,3,4 vibrational modes can be
specified for the species in a separate file using the optional
*vibfile* keyword, as discussed below.  If the :doc:`collide\_modify vibration <collide_modify>` command is used with a setting of
*discrete*\ , then this vibrational mode info must be specified for each
species with a vibdof = 4,6,8.  Note that the :doc:`fix vibmode <fix_vibmode>` command must also be used to allocate
per-particle storage for these additional modes.

NOTE: By default the maximum allowed number of vibrational modes is 4
(dof = 8). If you have a model with species which need more, you can
change the settings at the top of src/particle.h in the enum command
for MAXVIBMODE=4 to a larger value and re-compile the code.  The
format of the *vibfile*\ , as described next, is then enhanced
accordingly.


----------


The optional *rotfile* keyword can be used to specify the
characteristic rotational temperature(s) of species in the specified
*rfile*\ .  It is required for each linear (rotdof = 2) species when the
:doc:`collide\_modify rotate <collide_modify>` setting is *discrete*\ : the
rotational temperature Theta\_r defines the rigid-rotor level energies
E\_J = J(J+1) k Theta\_r used by the discrete rotational model.

The format of a species rotational file is as follows.  Comments or
blank lines are allowed in the file.  Comment lines start with a "#"
character.  All other lines must have the following format with values
separated by whitespace:


.. parsed-literal::

   species-ID N temp1 ... tempN sigma

The species-ID is a string that will be matched against the requested
species-ID, as described above.  N is the number of rotational
temperatures that follow: 1 for a linear molecule (rotdof = 2), or 3
for a nonlinear molecule (rotdof = 3).  The optional trailing value
*sigma* is the rotational symmetry number (e.g. 2 for homonuclear
diatomics, default 1), used in the rotational partition functions of
detailed-balance reverse reactions (see the :doc:`react <react>`
command).  Note that only the single temperature of linear molecules
is used by the discrete rotational model; nonlinear molecules keep the
continuous rotational model.


----------


The optional *vibfile* keyword can be used to specify additional
vibrational information in the specified *vfile*\ .  If this option is
used, then an entry must appear in *vfile* for every species in this
command with a vibdof value = 4,6,8.  Note that even if this
vibrational info is read, it is ignored by SPARTA unless the
:doc:`collide\_modify vibrate <collide_modify>` setting is specified as
*discrete*\ .

The format of a species vibrational file is as follows.  See
data/co2.species.vib for an example. Comments or blank lines are
allowed in the file.  Comment lines start with a "#" character.  All
other lines must have the following format with values separated by
whitespace:


.. parsed-literal::

   species-ID N temp1 relax1 degen1 temp2 relax2 degen2 ...

The species-ID is a string that will be matched against the requested
species-ID, as described above.  N is the total number of vibrational
oscillators, which must be 1,2,3,4 and must match the corresponding
vibdof value = 2,4,6,8 (divided by two) used in the species file.  A
vibration file entry is required for every species with vibdof > 2 and is
optional for a diatomic species (vibdof = 2, N = 1); when present it
overrides that species' vibrational temperature and relaxation from the
species file.

One or more distinct vibrational modes then follow, each listed as 3
values:

* tempI = vibrational temperature of mode I (temperature units)
* relaxI = inverse vibrational relaxation number of mode I (unitless)
* degenI = degeneracy of mode I (integer >= 1, unitless)

The degeneracy degenI is the number of independent oscillators of equal
frequency that make up mode I.  For example, the doubly-degenerate
bending mode of CO2 has degenI = 2.  The degeneracies of the listed
modes must sum to N, so the number of modes listed may be fewer than N.
Internally each mode is expanded into degenI independent oscillators, so
a degenerate mode may equivalently be listed once with degenI = g or as
g separate modes each with degeneracy 1.  For example, these two CO2
entries (vibdof = 8, so N = 4) are equivalent:


.. parsed-literal::

   CO2 4 1918.6 0.05 1 3382.0 0.05 1 959.0 0.05 2
   CO2 4 1918.6 0.05 1 3382.0 0.05 1 959.0 0.05 1 959.0 0.05 1

These quantities are used during collisions if vibrational energy is
modeled in discrete levels.

Note that the values for temp1 and relax1 override the same values
defined in the species file (prop7 and prop6) when they are listed for
the same species in the *vibfile*\ .


----------


The optional *elecfile* keyword can be used in conjuction with
"collide\_modify electronic discrete" and "fix elecmode" to model
electronic excitation. If this option is used, then an entry must
appear in *elecfile* for every species in this command.

The format of a species electronic file is as follows. Comments or
blank lines are allowed in the file. Comment lines start with a "#"
character. All other lines must have one of the following formats:


.. parsed-literal::

   species-ID N temp1 relax1 degen1 spin1 dof1 temp2 relax2 degen2 spin2 dof2 ... tempN relaxN degenN spinN dofN

The species-ID is a string that will be matched against the requested species-ID.
N is the number of electronic states. For each of the N states, 5 values are
listed:

* tempI = electronic temperature of state I (temperature units)
* relaxI = inverse electronic relaxation number of state I (unitless)
* degenI = total degeneracy of state I (integer, unitless)
* spinI = spin degeneracy of state I (e.g. singlet = 1, triplet = 3, etc.)
* dofI = effective DOF of state I for the TCE chemistry model (unitless)

The *dofI* value is read for file compatibility but not used: the TCE
reaction model couples electronic energy to the reaction probability
through the microcanonical ladder average (see
:doc:`react\_modify <react_modify>`), which requires no per-state DOF.

Note that the *relaxI* values enter a state-to-state transition twice:
a relaxation event is attempted with probability relax(current state),
and the destination state is then selected with a weight proportional
to degen \* relax(destination) times an available-energy factor.  The
effective i->f transition rate therefore scales as relaxI\*relaxF; this
product form is what preserves detailed balance for state-dependent
relaxation numbers, and should be accounted for when calibrating the
inputs against measured relaxation times or state-to-state rates.

The other accepted format is:


.. parsed-literal::

   species-ID1 species-ID2 conserve_spin relax1 ... relaxN

species-ID1 is a string of a species ID with electronic modes defined.
species-ID2 is another species ID with which species-ID1 will collide.
conserve\_spin is either a "T" or a "F". When "F", then any electronic state
for the relaxing particle is considered for the post-collision state. If "T"
then only states with the same spin degeneracy number will be considered for relaxation.
relaxI is the same as relaxI for the other format, i.e. inverse electronic
relaxation number for state I of species-ID1 when colliding with species-ID2.
These values override the default values set in the line with the first format.

The states of a species must be listed in ascending energy order, and the
first (ground) state must have a temperature of zero. A line in the second
format must come after the line in the first format that defines the states
of species-ID1.

**Restrictions:** none

**Related commands:** none

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
