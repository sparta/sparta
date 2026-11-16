.. index:: undump

undump command
==============

**Syntax:**


.. parsed-literal::

   undump dump-ID

* dump-ID = ID of previously defined dump

**Examples:**


.. parsed-literal::

   undump mine
   undump 2

**Description:**

Delete a dump that was previously defined with a :doc:`dump <fix>`
command.  This also closes the file associated with the dump.

**Restrictions:** none

**Related commands:**

:doc:`dump <dump>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
