.. index:: remove\_surf

remove\_surf command
====================

**Syntax:**


.. parsed-literal::

   remove_surf surfID

* surfID = group ID for which surface elements to remove

**Examples:**


.. parsed-literal::

   remove_surf topsurf

**Description:**

Remove a group of surface elements that have previously been read-in
via the :doc:`read\_surf <read_surf>` command.  The :doc:`group surf <group>` or :doc:`read\_surf <read_surf>` can be used to assign
each surface element to one or more groups.  This command removes all
surface elements in the specified *surfID* group.

Note that the remaining surface elements must still constitute a
"watertight" surface or an error will be generated.  The definition of
watertight is explained in the Restrictions section of the
:doc:`read\_surf <read_surf>` doc page.

After surface elements have been deleted the IDs of the remaining
surface points and elements are renumbered so that the remaining N
elements have IDs from 1 to N.  The new list of surface elements can
be output via the :doc:`write\_surf <write_surf>` or :doc:`dump surf <dump>` commands.

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

:doc:`read\_surf <read_surf>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
