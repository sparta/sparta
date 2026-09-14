10. Modifying & extending SPARTA
================================

This section describes how to extend SPARTA by modifying its source code.

| 10.1 `Compute styles <#mod_1>`_
| 10.2 `Fix styles <#mod_2>`_
| 10.3 `Region styles <#mod_3>`_
| 10.4 `Collision styles <#mod_4>`_
| 10.5 `Surface collision styles <#mod_5>`_
| 10.6 `Chemistry styles <#mod_6>`_
| 10.7 `Dump styles <#mod_7>`_
| 10.8 `Input script commands <#mod_8>`_ 
| 

SPARTA is designed in a modular fashion so as to be easy to modify and
extend with new functionality.

In this section, changes and additions users can make are listed along
with minimal instructions.  If you add a new feature to SPARTA and
think it will be of general interest to users, please submit it to the
`developers <https://sparta.github.io/authors.html>`_ for inclusion in
the released version of SPARTA.

The best way to add a new feature is to find a similar feature in
SPARTA and look at the corresponding source and header files to figure
out what it does. You will need some knowledge of C++ to be able to
understand the hi-level structure of SPARTA and its class
organization, but functions (class methods) that do actual
computations are written in vanilla C-style code and operate on simple
C-style data structures (vectors, arrays, structs).

The new features described in this section require you to write a new
C++ derived class. Creating a new class requires 2 files, a source
code file (\*.cpp) and a header file (\*.h).  The derived class must
provide certain methods to work as a new option.  Depending on how
different your new feature is compared to existing features, you can
either derive from the base class itself, or from a derived class that
already exists.  Enabling SPARTA to invoke the new class is as simple
as putting the two source files in the src dir and re-building SPARTA.

The advantage of C++ and its object-orientation is that all the code
and variables needed to define the new feature are in the 2 files you
write, and thus shouldn't make the rest of SPARTA more complex or
cause side-effect bugs.

Here is a concrete example. Suppose you write 2 files collide\_foo.cpp
and collide\_foo.h that define a new class CollideFoo that computes
inter-particle collisions described in the classic 1997 paper by Foo,
et al. If you wish to invoke those potentials in a SPARTA input script
with a command like

collide foo mix-ID params.foo 3.0

then your collide\_foo.h file should be structured as follows:

#ifdef COLLIDE\_CLASS
CollideStyle(foo,CollideFoo)
#else
...
(class definition for CollideFoo)
...
#endif

where "foo" is the style keyword in the collid command, and CollideFoo
is the class name defined in your collide\_foo.cpp and collide\_foo.h
files.

When you re-build SPARTA, your new collision model becomes part of the
executable and can be invoked with a :doc:`collide <collide>` command
like the example above.  Arguments like a mixture ID, params.foo (a
file with collision parameters), and 3.0 can be defined and processed
by your new class.

As illustrated by this example, many kinds of options are referred to
in the SPARTA documentation as the "style" of a particular command.

The instructions below give the header file for the base class that
these styles are derived from.  Public variables in that file are ones
used and set by the derived classes which are also used by the base
class.  Sometimes they are also used by the rest of SPARTA.  Virtual
functions in the base class header file which are set = 0 are ones
that must be defined in the new derived class to give it the
functionality SPARTA expects.  Virtual functions that are not set to 0
are functions that can be optionally defined.

Here are additional guidelines for modifying SPARTA and adding new
functionality:

* Think about whether what you want to do would be better as a pre- or
  post-processing step. Many computations are more easily and more
  quickly done that way.
* Don't do anything within the timestepping of a run that isn't
  parallel.  E.g. don't accumulate a large volume of data on a single
  processor and analyze it.  This runs the risk of seriously degrading
  the parallel efficiency.

  If you have a question about how to compute something or about
  internal SPARTA data structures or algorithms, feel free to send an
  email to the `developers <https://sparta.github.io/authors.html>`_.

* If you add something you think is generally useful, also send an email
  to the `developers <https://sparta.github.io/authors.html>`_ so we can
  consider adding it to the SPARTA distribution.





.. _mod_1:

.. raw:: html

   <span id="mod_1"></span>

10.1 Compute styles
-------------------

:doc:`Compute style commands <compute>` calculate instantaneous
properties of the simulated system.  They can be global properties, or
per particle or per grid cell or per surface element properties.  The
result can be single value or multiple values (global or per particle
or per grid or per surf).

Here is a brief description of methods to define in a new derived
class.  See compute.h for details.  All of these methods are optional.

+------------------------+-----------------------------------------------------+
| init                   | initialization before a run                         |
+------------------------+-----------------------------------------------------+
| compute\_scalar        | compute a global scalar quantity                    |
+------------------------+-----------------------------------------------------+
| compute\_vector        | compute a global vector of quantities               |
+------------------------+-----------------------------------------------------+
| compute\_per\_particle | compute one or more quantities per particle         |
+------------------------+-----------------------------------------------------+
| compute\_per\_grid     | compute one or more quantities per grid cell        |
+------------------------+-----------------------------------------------------+
| compute\_per\_surf     | compute one or more quantities per surface element  |
+------------------------+-----------------------------------------------------+
| surf\_tally            | call when a particle hits a surface element         |
+------------------------+-----------------------------------------------------+
| boundary\_tally        | call when a particle hits a simulation box boundary |
+------------------------+-----------------------------------------------------+
| memory\_usage          | tally memory usage                                  |
+------------------------+-----------------------------------------------------+

Note that computes with "/particle" in their style name calculate per
particle quantities, with "/grid" in their name calculate per grid
cell quantities, and with "/surf" in their name calculate per surface
element properties.  All others calcuulate global quantities.

Flags may also need to be set by a compute to enable specific
properties.  See the compute.h header file for one-line descriptions.


----------


.. _mod_2:

.. raw:: html

   <span id="mod_2"></span>

10.2 Fix styles
---------------

:doc:`Fix style commands <fix>` perform operations during the
timestepping loop of a simulation.  They can define methods which are
invoked at different points within the timestep.  They can be used to
insert particles, perform load-balancing, or perform time-averaging of
various quantities.  They can also define and maintain new
per-particle vectors and arrays that define quantities that move with
particles when they migrate from processor to processor or when the
grid is rebalanced or adapated.  They can also produce output of
various kinds, similar to :doc:`compute <compute>` commands.

Here is a brief description of methods to define in a new derived
class.  See fix.h for details.  All of these methods are optional,
except setmask().

+-----------------+-------------------------------------------------------------------+
| setmask         | set flags that determine when the fix is called within a timestep |
+-----------------+-------------------------------------------------------------------+
| init            | initialization before a run                                       |
+-----------------+-------------------------------------------------------------------+
| start\_of\_step | called at beginning of timestep                                   |
+-----------------+-------------------------------------------------------------------+
| end\_of\_step   | called at end of timestep                                         |
+-----------------+-------------------------------------------------------------------+
| add\_particle   | called when a particle is created                                 |
+-----------------+-------------------------------------------------------------------+
| surf\_react     | called when a surface reaction occurs                             |
+-----------------+-------------------------------------------------------------------+
| memory\_usage   | tally memory usage                                                |
+-----------------+-------------------------------------------------------------------+

Flags may also need to be set by a fix to enable specific properties.
See the fix.h header file for one-line descriptions.

Fixes can interact with the Particle class to create new
per-particle vectors and arrays and access and update their
values.  These are the relevant Particle class methods:

+----------------+--------------------------------------------------+
| add\_custom    | add a new custom vector or array                 |
+----------------+--------------------------------------------------+
| find\_custom   | find a previously defined custom vector or array |
+----------------+--------------------------------------------------+
| remove\_custom | remove a custom vector or array                  |
+----------------+--------------------------------------------------+

See the :doc:`fix ambipolar <fix_ambipolar>` for an example of how these
are used.  It define an integer vector called "ionambi" to flag
particles as ambipolar ions, and a floatin-point array called
"velambi" to store the velocity vector for the associated electron.


----------


.. _mod_3:

.. raw:: html

   <span id="mod_3"></span>

10.3 Region styles
------------------

:doc:`Region style commands <region>` define geometric regions
within the simulation box.  Other commands use regions
to limit their computational scope.

Here is a brief description of methods to define in a new derived
class.  See region.h for details.  The inside() method is required.

inside: determine whether a point is inside/outside the region


----------


.. _mod_4:

.. raw:: html

   <span id="mod_4"></span>

10.4 Collision styles
---------------------

:doc:`Collision style commands <collide>` define collision models that
calculate interactions between particles in the same grid cell.

Here is a brief description of methods to define in a new derived
class.  See collide.h for details.  All of these methods are required
except init() and modify\_params().

+--------------------+---------------------------------------------------------------------------------------+
| init               | initialization before a run                                                           |
+--------------------+---------------------------------------------------------------------------------------+
| modify\_params     | process style-specific options of the :doc:`collide\_modify <collide_modify>` command |
+--------------------+---------------------------------------------------------------------------------------+
| vremax\_init       | estimate VREmax settings                                                              |
+--------------------+---------------------------------------------------------------------------------------+
| attempt\_collision | compute # of collisions to attempt for entire cell                                    |
+--------------------+---------------------------------------------------------------------------------------+
| attempt\_collision | compute # of collisions to attempt between 2 species groups                           |
+--------------------+---------------------------------------------------------------------------------------+
| test\_collision    | determine if a collision bewteen 2 particles occurs                                   |
+--------------------+---------------------------------------------------------------------------------------+
| setup\_collision   | pre-computation before a 2-particle collision                                         |
+--------------------+---------------------------------------------------------------------------------------+
| perform\_collision | calculate the outcome of a 2-particle collision                                       |
+--------------------+---------------------------------------------------------------------------------------+


----------


.. _mod_5:

.. raw:: html

   <span id="mod_5"></span>

10.5 Surface collision styles
-----------------------------

:doc:`Surface collision style commands <collide>` define collision
models that calculate interactions between a particle and surface
element.

Here is a brief description of methods to define in a new derived
class.  See surf\_collide.h for details.  All of these methods are
required except dynamic().

+---------+------------------------------------------------------+
| init    | initialization before a run                          |
+---------+------------------------------------------------------+
| collide | perform a particle/surface-element collision         |
+---------+------------------------------------------------------+
| dynamic | allow surface property to change during a simulation |
+---------+------------------------------------------------------+


----------


.. _mod_6:

.. raw:: html

   <span id="mod_6"></span>

10.6 Chemistry styles
---------------------

Particle/particle chemistry models in SPARTA are specified by
:doc:`reaction style commands <react>` which define lists of possible
reactions and their parameters.

Here is a brief description of methods to define in a new derived
class.  See react.h for details.  The init() method is optional;
the attempt() method is required.

+---------+---------------------------------------------------+
| init    | initialization before a run                       |
+---------+---------------------------------------------------+
| attempt | attempt a chemical reaction between two particles |
+---------+---------------------------------------------------+


----------


.. _mod_7:

.. raw:: html

   <span id="mod_7"></span>

10.7 Dump styles
----------------

:doc:`Dump commands <dump>` output snapshots of simulation data to a
file periodically during a simulation, in a particular file format.
Per particle, per grid cell, or per surface element data can be
output.

Here is a brief description of methods to define in a new derived
class.  See dump.h for details.  The init\_style(), modify\_param(), and
memory\_usage() methods are optional; all the others are required.

+---------------+---------------------------------------------------------------------------------+
| init\_style   | style-specific initialization before a run                                      |
+---------------+---------------------------------------------------------------------------------+
| modify\_param | process style-specific options of the :doc:`dump\_modify <dump_modify>` command |
+---------------+---------------------------------------------------------------------------------+
| write\_header | write the header of a snapshot to a file                                        |
+---------------+---------------------------------------------------------------------------------+
| count         | # of entities this processor will output                                        |
+---------------+---------------------------------------------------------------------------------+
| pack          | pack a processor's data into a buffer                                           |
+---------------+---------------------------------------------------------------------------------+
| write\_data   | write a buffer of data to a file                                                |
+---------------+---------------------------------------------------------------------------------+
| memory\_usage | tally memory usage                                                              |
+---------------+---------------------------------------------------------------------------------+


----------


.. _mod_8:

.. raw:: html

   <span id="mod_8"></span>

10.8 Input script commands
--------------------------

New commands can be added to SPARTA that will be recognized in input
scripts.  For example, the :doc:`create\_particles <create_particles>`,
:doc:`read\_surf <read_surf>`, and :doc:`run <run>` commands are all
implemented in this fashion.  When such a command is encountered in an
input script, SPARTA simply creates a class with the corresponding
name, invokes the "command" method of the class, and passes it the
arguments from the input script.  The command() method can perform
whatever operations it wishes on SPARTA data structures.

The single method the new class must define is as follows:

+---------+--------------------------------------------------+
| command | operations performed by the input script command |
+---------+--------------------------------------------------+

Of course, the new class can define other methods and variables as
needed.


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
