1. Introduction
===============

These sections provide an overview of what SPARTA can do, describe
what it means for SPARTA to be an open-source code, and acknowledge
the funding and people who have contributed to SPARTA.

| 1.1 `What is SPARTA <#intro_1>`_
| 1.2 `SPARTA features <#intro_2>`_
| 1.3 `Grids and surfaces in SPARTA <#intro_3>`_
| 1.4 `Open source distribution <#intro_4>`_
| 1.5 `Acknowledgments and citations <#intro_5>`_ 
| 


----------


.. _intro_1:

.. raw:: html

   <span id="intro_1"></span>

1.1 What is SPARTA
------------------

SPARTA is a Direct Simulation Montel Carlo code that models rarefied
gases, using collision, chemistry, and boundary condition models.  It
uses a hierarchical Cartesian grid to track and group particles for 3d
or 2d or axisymmetric models.  Objects emedded in the gas are
represented as triangulated surfaces and cut through grid cells.

For examples of SPARTA simulations, see the `SPARTA WWW Site <sws_>`_.

SPARTA runs efficiently on single-processor desktop or laptop
machines, but is designed for parallel computers.  It will run on any
parallel machine that compiles C++ and supports the `MPI <mpi_>`_
message-passing library.  This includes distributed- or shared-memory
parallel machines as well as commodity clusters.

.. _mpi: https://www.mcs.anl.gov/research/projects/mpi/



SPARTA can model systems with only a few particles up to millions or
billions.  See :doc:`Section 8 <Section_perf>` for information on SPARTA
performance and scalability, or the Benchmarks section of the `SPARTA WWW Site <sws_>`_.

SPARTA is a freely-available open-source code, distributed under the
terms of the `GNU Public License <gnu_>`_, or sometimes by request under
the terms of the `GNU Lesser General Public License (LGPL) <https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html>`_,
which means you can use or modify the code however you wish.  The only
restrictions imposed by the GPL or LGPL are on how you distribute the
code further.  See `Section 1.4 <#intro_4>`_ below for a brief discussion
of the open-source philosophy.

.. _gnu: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html



SPARTA is designed to be easy to modify or extend with new
capabilities, such as new collision or chemistry models, boundary
conditions, or diagnostics.  See :doc:`Section 10 <Section_modify>` for
more details.

SPARTA is written in C++ which is used at a hi-level to structure the
code and its options in an object-oriented fashion.  The kernel
computations use simple data structures and C-like code for effciency.
So SPARTA is really written in an object-oriented C style.

SPARTA was developed with internal funding at `Sandia National Laboratories <snl_>`_, a US Department of Energy lab.  See `Section 1.5 <#intro_5>`_ below for more information on SPARTA funding and
individuals who have contributed to SPARTA.

.. _snl: https://www.sandia.gov




----------


.. _intro_2:

.. raw:: html

   <span id="intro_2"></span>

1.2 SPARTA features
-------------------

This section highlights SPARTA features, with links to specific
commands which give more details.  The `next section <#intro_3>`_
illustrates the kinds of grid geometries and surface definitions which
SPARTA supports.

If SPARTA doesn't have your favorite collision model, boundary
condition, or diagnostic, see :doc:`Section 10 <Section_modify>` of the
manual, which describes how it can be added to SPARTA.

General features
----------------

* runs on a single processor or in parallel
* distributed-memory message-passing parallelism (MPI)
* spatial-decomposition of simulation domain for parallelism
* open-source distribution
* highly portable C++
* optional libraries used: MPI
* :doc:`easy to extend <Section_modify>` with new features and functionality
* runs from an :doc:`input script <Section_commands>`
* syntax for defining and using :doc:`variables and formulas <variable>`
* syntax for :doc:`looping over runs <jump>` and breaking out of loops
* run one or `multiple simulations simultaneously <Section_howto.html#howto_3>`_ (in parallel) from one script
* `build as library <Section_start.html#start_4>`_, invoke SPARTA thru `library interface <Section_howto.html#howto_6>`_ or provided :doc:`Python wrapper <Section_python>`
* `couple with other codes <Section_howto.html#howto_7>`_: SPARTA calls other code, other code calls SPARTA, umbrella code calls both

Models
------

* :doc:`3d or 2d <dimension>` or `2d-axisymmetric <Section_howto.html#howto_2>`_ domains
* variety of :doc:`global boundary conditions <boundary>`
* :doc:`create particles <create_particles>` within flow volume
* emit particles from simulation box faces due to :doc:`flow properties <fix_emit_face>`
* emit particles from simulation box faces due to :doc:`profile defined in file <fix_emit_face_file>`
* emit particles from surface elements due to :doc:`normal and flow properties <fix_emit_surf>`
* `ambipolar <Section_howto.html#howto_11>`_ approximation for ionized plasmas

Geometry
--------

* `Cartesian, heirarchical grids <#intro_3>`_ with multiple levels of local refinement
* :doc:`create grid from input script <create_grid>` or :doc:`read from file <read_grid>`
* embed `triangulated (3d) or line-segmented (2d) surfaces <#intro_3>`_ in grid, :doc:`read in from file <read_surf>`

Gas-phase collisions and chemistry
----------------------------------

* collisions between all particles or pairs of species groups within grid cells
* :doc:`collision models: <collide>` VSS (variable soft sphere), VHS (variable hard sphere), HS (hard sphere)
* :doc:`chemistry models: <react>` TCE, QK

Surface collisions and chemistry
--------------------------------

* for surface elements or global simulation box :doc:`boundaries <bound_modify>`
* :doc:`collisions: <surf_collide>` specular or diffuse
* :doc:`reactions <surf_react>`

Performance
-----------

* :doc:`grid cell weighting <global>` of particles
* :doc:`adaptation <adapt_grid>` of the grid cells between runs
* :doc:`on-the-fly adaptation <fix_adapt>` of the grid cells
* :doc:`static <balance_grid>` load-balancing of grid cells or particles
* :doc:`dynamic <fix_balance>` load-balancing of grid cells or particles

Diagnostics
-----------

* :doc:`global boundary statistics <compute_boundary>`
* :doc:`per grid cell statistics <compute_grid>`
* :doc:`per surface element statistics <compute_surf>`
* time-averaging of :doc:`global <fix_ave_time>`, :doc:`grid <fix_ave_grid>`, :doc:`surface <fix_ave_surf>` statistics

Output
------

* :doc:`log file of statistical info <stats_style>`
* :doc:`dump files <dump>` (text or binary) of per particle, per grid cell, per surface element values
* binary :doc:`restart files <restart>`
* on-the-fly :doc:`rendered images and movies <dump_image>` of particles, grid cells, surface elements

Pre- and post-processing
------------------------

* Various pre- and post-processing serial tools are packaged with
  SPARTA; see :doc:`Section 9 <Section_tools>` of the manual.
* Our group has also written and released a separate toolkit called
  `Pizza.py <pizza_>`_ which provides tools for doing setup, analysis,
  plotting, and visualization for SPARTA simulations.  Pizza.py is
  written in `Python <python_>`_ and is available for download from `the Pizza.py WWW site <pizza_>`_.

.. _pizza: https://lammps.github.io/pizza



.. _python: https://www.python.org




----------


.. _intro_3:

.. raw:: html

   <span id="intro_3"></span>

1.3 Grids and surfaces in SPARTA
--------------------------------

SPARTA overlays a grid over the simulation domain which is used to
track particles and to co-locate particles in the same grid cell for
performing collision and chemistry operations.  SPARTA uses a
Cartesian hierarchical grid.  Cartesian means that the faces of a grid
cell are aligned with the Cartesian xyz axes.  Hierarchical means that
individual grid cells can be sub-divided into smaller cells,
recursively.  This allows for flexible grid cell refinement in any
region of the simulation domain.  E.g. around a surface, or in a
high-density region of the gas flow.

An example 2d hierarchical grid is shown in the diagram, for a
circular surface object (in red) with the grid refined on the upwind
side of the object (flow from left to right).

.. image:: JPG/refine_grid.jpg
   :align: center

Objects represented with a surface triangulation (line segments in 2d)
can also be read in to define objects which particles flow around.
Individual surface elements are assigned to grid cells they intersect
with, so that particle/surface collisions can be efficiently computed.

As an example, here is coarsely triangulated representation of the
space shuttle (only 616 triangles!), which could be embedded in a
simulation box.  Click on the image for a larger picture.

.. image:: JPG/shuttle_small.jpg
   :target: JPG/shuttle.jpg
   :align: center

See `Sections 6.9 <Section_howto.html#howto_9>`_ and
`4.10 <Section_howto.html#>`_ for more details of both the grids and
surface objects that SPARTA supports and how to define them.


----------


.. _intro_4:

.. raw:: html

   <span id="intro_4"></span>

1.4 Open source distribution
----------------------------

SPARTA comes with no warranty of any kind.  As each source file states
in its header, it is a copyrighted code that is distributed free-of-
charge, under the terms of the `GNU Public License <gnu_>`_ (GPL).  This
is often referred to as open-source distribution - see
`www.gnu.org <gnuorg_>`_ or `www.opensource.org <opensource_>`_ for more
details.  The legal text of the GPL is in the LICENSE file that is
included in the SPARTA distribution.

.. _gnuorg: http://www.gnu.org



.. _opensource: https://opensource.org



Here is a summary of what the GPL means for SPARTA users:

\(1) Anyone is free to use, modify, or extend SPARTA in any way they
choose, including for commercial purposes.

\(2) If you distribute a modified version of SPARTA, it must remain
open-source, meaning you distribute it under the terms of the GPL.
You should clearly annotate such a code as a derivative version of
SPARTA.

\(3) If you release any code that includes SPARTA source code, then it
must also be open-sourced, meaning you distribute it under the terms
of the GPL.

\(4) If you give SPARTA files to someone else, the GPL LICENSE file and
source file headers (including the copyright and GPL notices) should
remain part of the code.

In the spirit of an open-source code, these are various ways you can
contribute to making SPARTA better.  You can send email to the
`developers <https://sparta.github.io/authors.html>`_ on any of these
topics.

* Point prospective users to the `SPARTA WWW Site <sws_>`_.  Mention it in
  talks or link to it from your WWW site.
* If you find an error or omission in this manual or on the `SPARTA WWW Site <sws_>`_, or have a suggestion for something to clarify or include,
  send an email to the
  `developers <https://sparta.github.io/authors.html>`_.
* If you find a bug, `Section 12.1 <Section_errors.html#err_2>`_ describes
  how to report it.
* If you publish a paper using SPARTA results, send the citation (and
  any cool pictures or movies) to add to the Publications, Pictures, and
  Movies pages of the `SPARTA WWW Site <sws_>`_, with links and attributions
  back to you.
* The tools sub-directory of the SPARTA distribution has various
  stand-alone codes for pre- and post-processing of SPARTA data.  More
  details are given in :doc:`Section 9 <Section_tools>`.  If you write a
  new tool that others will find useful, it can be added to the SPARTA
  distribution.
* SPARTA is designed to be easy to extend with new code for features
  like boundary conditions, collision or chemistry models, diagnostic
  computations, etc.  :doc:`Section 10 <Section_modify>` of the manual
  gives details.  If you add a feature of general interest, it can be
  added to the SPARTA distribution.
* The Benchmark page of the `SPARTA WWW Site <sws_>`_ lists SPARTA
  performance on various platforms.  The files needed to run the
  benchmarks are part of the SPARTA distribution.  If your machine is
  sufficiently different from those listed, your timing data can be
  added to the page.
* Cash.  Small denominations, unmarked bills preferred.  Paper sack OK.
  Leave on desk.  VISA also accepted.  Chocolate chip cookies
  encouraged.


----------


.. _intro_5:

.. raw:: html

   <span id="intro_5"></span>

1.5 Acknowledgments and citations
---------------------------------------------------------------------------------------------

SPARTA development has been funded by the `US Department of Energy <doe_>`_ (DOE).

.. _doe: https://www.energy.gov



If you use SPARTA results in your published work, please cite the
paper(s) listed under the `Citing SPARTA link <https://sparta.github.io/papers.html>`_ of the SPARTA WWW page, and
include a pointer to the `SPARTA WWW Site <sws_>`_
(https://sparta.github.io):

The `Publications link <https://sparta.github.io/papers.html>`_ on the
SPARTA WWW page lists papers that have cited SPARTA.  If your paper is
not listed there, feel free to send us the info.  If the simulations
in your paper produced cool pictures or animations, we'll be pleased
to add them to the `Pictures <https://sparta.github.io/pictures.html>`_
or `Movies <https://sparta.github.io/pictures.html>`_ pages of the SPARTA
WWW site.

The core group of SPARTA developers is at Sandia National Labs:

* Steve Plimpton, sjplimp at gmail.com
* Michael Gallis, magalli at sandia.gov


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
