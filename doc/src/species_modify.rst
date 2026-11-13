.. index:: species\_modify

species\_modify command
=======================

**Syntax:**


.. parsed-literal::

   species_modify ID property value ...

* ID, property, value can be repeated one or more times
* ID = species ID
* property = *mu*
  
  .. parsed-literal::
  
       mu = magnetic moment

* value = value of property for that species
  
  .. parsed-literal::
  
       value for mu (magnetic moment units)



**Examples:**


.. parsed-literal::

   species_modify Fe mu 2.0 Cr mu 3.0

**Description:**

Set additional properties of one or more species used in a simulation.
This can be used as many times as desired for different species and
properties.  Currently it only supports setting of a single optional
property (the magnetic moment) which is not included in the species
files read in by the :doc:`species <species>` command.

Each *ID* is a character string used to identify a species, such as N
or O2 or NO or D or Fe-.  See the :doc:`species <species>` command for
how species are added to a simulation model by reading their
properties from a species file.

The only property currently recognized is *mu* or the scalar magnetic
moment of each particle of the species.  The *value* for the *mu*
property should be specified in the units described on the
:doc:`units <units>` doc page.


----------


**Restrictions:** none

**Related commands:** none

**Default:**

No magnetic moments are defined for any species (all 0.0).


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
