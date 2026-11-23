.. index:: suffix

.. _suffix-command:

suffix
======

SPARTA has no "suffix" input-script command.  Earlier versions of this
manual documented one, but SPARTA never accepted it: a "suffix" line in
an input script stops with an "Unknown command" error.

To use the accelerated variants of styles, such as those of the KOKKOS
package, do one of these:

* Launch SPARTA with the "-suffix" (or "-sf") `command-line switch <Section_start.html#start_7>`_, e.g. "-sf kk".  The suffix is then
  appended to every style that has a variant, without editing the input
  script.
* Name the variant explicitly in your input script, e.g.

.. code-block:: SPARTA

   collide vss/kk air ar.vss

See the :doc:`Accelerating SPARTA <Section_accelerate>` section of the
manual for more about the KOKKOS package.
