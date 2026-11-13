.. index:: compute count

compute count command
=====================

compute count/kk command
========================

**Syntax:**


.. parsed-literal::

   compute ID count id1 id2 ...

* ID is documented in :doc:`compute <compute>` command
* count = style name of this compute command
* id1,id2,... = species ID or mixture ID or mixture/group
  
  .. parsed-literal::
  
       species ID = ID used with the :doc:`species <species>` command
       mixture ID = ID used with the :doc:`mixture <mixture>` command, expands to all groups in mixture
       mixture/group = ID of mixture followed by name of a group within mixture



**Examples:**


.. parsed-literal::

   compute 1 count species
   compute Ncounts count N N2 N+ air/O

**Description:**

Define a computation that counts the number of particles currently in
the simulation for various species or groups within mixtures.  Groups
are collections of one or more species within a mixture.  See the
"mixture" command for an explanation of how species are added to a
mixture and how groups of species within the mixture are defined.

Each of the listed *ids* (id1, id2, etc) can be in one of three
formats.  Any of the ids can be in any of the formats.

An *id* can be a species ID, in which case the count is for particles
of that species.

An *id* can be a mixture ID, in which case one count is performed for
each of the groups within the mixture.  In the first example above,
"species" is the name of a default mixture which assigns every species
defined for the simulation to its own group.  If there are 10 species
in the simulation, there will thus be 10 counts calculated, the same
as if the command had been specified with explicit names for all 10
species, e.g.


.. parsed-literal::

   compute 1 count O2 N2 O N NO O2+ N2+ O+ N+ NO+

An *id* can also be of the form mix-ID/name where mix-ID is a mixture
ID and name is the name of a group in that mixture.


----------


**Output info:**

If there is a single count accumulated, this compute calculates a
global scalar.  If there are multiple counts accumulated, it
calculates a global vector with a length = number of counts.  These
results can be used by any command that uses global scalar or vector
values from a compute as input.  See `Section 6.4 <Section_howto.html#howto_4>`_ for an overview of SPARTA output
options.

The values will all be unitless counts.


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

It is an error if a listed *id* is both a species ID and a mixture ID,
since this command cannot distinguish between them.

**Related commands:** none

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
