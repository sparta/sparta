13. Future and history
======================

This section lists features we are planning to add to SPARTA, features
of previous versions of SPARTA, and features of other parallel
molecular dynamics codes I've distributed.

| 13.1 `Coming attractions <#hist_1>`_
| 13.2 `Past versions <#hist_2>`_ 
| 


----------


.. _hist_1:

.. raw:: html

   <span id="hist_1"></span>

13.1 Coming attractions
---------------------------------------------------------------------------------

Features that have been requested but not yet implemented are tracked
as `issues <https://github.com/sparta/sparta/issues>`_ in the SPARTA
GitHub repository.  Please contact the
`developers <https://sparta.github.io/authors.html>`_ if you are
interested in contributing to any of those developments, or would be a
future user of that feature.

You can request a new feature by opening an issue on Github.


----------


.. _hist_2:

.. raw:: html

   <span id="hist_2"></span>

13.2 Past versions
----------------------------------------------------------------------------

Sandia's predecessor to SPARTA is a DSMC code called ICARUS.  It was
developed in the early 1990s by Tim Bartel and `Steve Plimpton <https://sjplimp.github.io>`_.  It was later modified and
extended by Michael Gallis.

ICARUS is a 2d code, written in Fortran, which models the flow
geometry around bodies with a collection of adjoining body-fitted grid
blocks.  The geometry of the grid cells within in a single block is
represented with analytic equations, which allows for fast particle
tracking.

Some details about ICARUS, including simulation snapshots and papers,
are discussed on `this page <https://sjplimp.github.io/dsmc.html>`_

Performance-wise ICARUS scaled quite well on several generations of
parallel machines, and is still used by Sandia researchers today.
ICARUS was export-controlled software, and so was not distributed
widely outside of Sandia.

SPARTA development began in late 2011.  In contrast to ICARUS, it is a
3d code, written in C++, and uses a hierarchical Cartesian grid to
track particles.  Surfaces are embedded in the grid, which cuts and
splits their flow volumes.

The `Authors link <https://sparta.github.io/authors.html>`_ on the SPARTA
web page gives a timeline of features added to the code since it's
initial open-source release.


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
