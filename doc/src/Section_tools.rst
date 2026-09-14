9. Additional tools
===================

SPARTA is designed to be a computational kernel for performing DSMC
computations.  Additional pre- and post-processing steps are often
necessary to setup and analyze a simulation.  A few additional tools
are provided with the SPARTA distribution in the tools directory and
are described briefly below.

Our group has also written and released a separate toolkit called
`Pizza.py <pizza_>`_ which provides tools for doing setup, analysis,
plotting, and visualization for SPARTA simulations.  Pizza.py is
written in `Python <python_>`_ and is available for download from `the Pizza.py web site <pizza_>`_.

.. _pizza: https://lammps.github.io/pizza



.. _python: https://www.python.org



Some of the Pizza.py tools relevant to SPARTA are as follows:

* dump - read, write, manipulate particle dump files
* gl - 3d interactive visualization via OpenGL of dump or surface files
* sdata - read, write, manipulate surface files
* olog - read log files and extract columns of data
* vcr - VCR-style GUI for 3d interactive OpenGL visualization of dump or surface files

The dump, sdata, and olog tools are included in the SPARTA
distribution in the tools/pizza directory, and are used by some of the
scripts discussed below.

This is the list of tools included in the tools directory of the
SPARTA distribution.  Each is described in more detail below.

* `dump2cfg <#dump2cfg>`_ - convert a particle dump file to CFG format
* `dump2xyz <#dump2xyz>`_ - convert a particle dump file to XYZ format
* `grid\_refine <#gridrefine>`_ - refine a grid around a surface
* `implicit\_grid <#implicit>`_ - create a random porous region with implicit surfaces
* `jagged <#jagged>`_ - create jagged 2d/3d surfaces with explicit surfaces
* `log2txt <#log2txt>`_ - extract columns of info from a log file
* `logplot <#logplot>`_ - plot columns of info from a log file via GnuPlot
* `paraview <#paraview>`_ - converters of SPARTA data to `ParaView <ParaView website_>`_ format
* `stl2surf <#stl2surf>`_ - convert an STL text file into a SPARTA surface file
* `surf\_create <#surfcreate>`_ - create a surface file with simple objects
* `surf\_transform <#surftransform>`_ - transform surface via tranlate/scale/rotate operations

.. _ParaView website: https://www.paraview.org







.. raw:: html

   <span id="dump2cfg"></span>

dump2cfg tool
---------------------------------------------------------------------------

This is a Python script that converts a SPARTA particle dump file into
extended CFG format so that it can be visualized by the
`AtomEye <http://li.mit.edu/Archive/Graphics/A/>`_ visualization
program.  AtomEye is a very fast particle visualizer, capable of
interactive visualizations of millions of particles on a desktop
machine.  It is commonly used in the materials modeling community.

See the header of the script for the syntax used to run it.

This script uses one or more of the "Pizza.py" tools provided in the
tools/pizza directory.  See the tools/README file for info on how to
set an environment variable so that the Pizza.py tool files can be
found by Python, as well as instructions on various ways to run a
Python script.


----------


.. raw:: html

   <span id="dump2xyz"></span>

dump2xyz tool
---------------------------------------------------------------------------

This is a Python script that converts a SPARTA particle dump file into
XYZ format so that it can be visualized by various visualization
packages that read XYZ formatted files.  An example is
`VMD <https://www.ks.uiuc.edu/Research/vmd>`_ package, commonly used in
the molecular dynamics modeling community.

See the header of the script for the syntax used to run it.

This script uses one or more of the "Pizza.py" tools provided in the
tools/pizza directory.  See the tools/README file for info on how to
set an environment variable so that the Pizza.py tool files can be
found by Python, as well as instructions on various ways to run a
Python script.


----------


.. _gridrefine:

.. raw:: html

   <span id="gridrefine"></span>

grid\_refine tool
-----------------------------------------------------------------------------------

This is a Python script that creates a SPARTA grid file adapted
around the lines or triangles in a SPARTA surface file.  The resulting
grid file can be read by the :doc:`read\_grid <read_grid>` command.
The surface file can be read by the :doc:`read\_surf <read_surf>` command.

See the header of the script for the various adaptivity options that
are supported, and the syntax used to run it.


----------


.. _implicit:

.. raw:: html

   <span id="implicit"></span>

implicit\_grid tool
---------------------------------------------------------------------------------

This is a Python script which can be used to generate binary files
representing porous media samples, as read by the
:doc:`read\_isurf <read_isurf>` command.  The output files contain
randomized grid corner point values which induce implicit surfaces
which can contain huge numbers of surface elements.  They are useful
for stress testing the implicit surface options in SPARTA, as selected
by the :doc:`global surfs <global>` command.

See the header of the script for the syntax used to run it.

The examples/implicit directory uses these files as input.


----------


.. raw:: html

   <span id="jagged"></span>

jagged tools
----------------------------------------------------------------------

These are 2 Python scripts (jagged2d.py and jagged3d.py) which can be
used to generate SPARTA surface files in a pattern that can be very
jagged.  The surfaces can contain huge numbers of surface elements and
be read by the :doc:`read\_surf <read_surf>` command.  They are useful
for stress testing the explict surface options in SPARTA, including
distributed or non-distributed storage, as selected by the :doc:`global surfs <global>` command.

See the header of the scripts for the syntax used to run them.

The examples/jagged directory uses these files as input.


----------


.. raw:: html

   <span id="log2txt"></span>

log2txt tool
------------------------------------------------------------------------

This is a Python script that reads a SPARTA log file, extracts
selected columns of statistical output, and writes them to a text
file.  It knows how to concatenate log file info across multiple
successive runs.  The columnar output can then be read by various
plotting packages.

See the header of the script for the syntax used to run it.

This script uses one or more of the "Pizza.py" tools provided in the
tools/pizza directory.  See the tools/README file for info on how to
set an environment variable so that the Pizza.py tool files can be
found by Python, as well as instructions on various ways to run a
Python script.


----------


.. raw:: html

   <span id="logplot"></span>

logplot tool
------------------------------------------------------------------------

This is a Python script that reads a SPARTA log file, extracts the
selected columns of statistical output, and plots them via the GnuPlot
program.  It knows how to concatenate log file info across multiple
successive runs.

See the header of the script for the syntax used to run it.  You must
have GnuPlot installed on your system to use this script.  If you can
type "gnuplot" from the command line to start GnuPlot, it should work.
If not (e.g. because you need a path name), then edit these 2 lines as
needed in pizza/gnu.py:


.. parsed-literal::

   except: PIZZA_GNUPLOT = "gnuplot"
   except: PIZZA_GNUTERM = "x11"

For example, the first could become "/home/smith/bin/gnuplot".  The
second should only need changing if GnuPlot requires a different
setting to plot to your screen.

This script uses one or more of the "Pizza.py" tools provided in the
tools/pizza directory.  See the tools/README file for info on how to
set an environment variable so that the Pizza.py tool files can be
found by Python, as well as instructions on various ways to run a
Python script.


----------


.. raw:: html

   <span id="paraview"></span>

paraview tools
----------------------------------------------------------------------------

The tools/paraview directory has scripts which convert
SPARTA grid and surface data (input and output) to ParaView format.

`ParaView <ParaView website_>`_ is a popular, powerful, freely-available
visualization package.  You must have ParaView installed to use the
Python scripts.  See `Section 6.16 <Section_howto.html#howto_16>`_ for more details.

The scripts were developed by Tom Otahal (Sandia).


----------


.. raw:: html

   <span id="stl2surf"></span>

stl2surf tool
---------------------------------------------------------------------------

This is a Python script that reads a stereolithography (STL) text file
and converts it to a SPARTA surface file.  STL files contain a
collection of triangles and can be created by various mesh-generation
programs.  The format for SPARTA surface files is described on the
:doc:`read\_surf <read_surf>` command doc page.

See the header of the script for the syntax used to run it, e.g.


.. parsed-literal::

   % python stl2surf.py stlfile surffile

The script also checks the triangulated object to see if it is
"watertight" and issues a warning if it is not, since SPARTA will
perform the same check.  The :doc:`read\_surf <read_surf>` command doc
page explains what watertight means for 3d objects.


----------


.. _surfcreate:

.. raw:: html

   <span id="surfcreate"></span>

surf\_create tool
-----------------------------------------------------------------------------------

This is a Python script that creates a SPARTA surface file containing
one or more simple objects whose surface is represented as triangules
(3d) or line segments (2d).  Such files can be read by the
:doc:`read\_surf <read_surf>` command.  The 3d objects it supports are a
sphere, box, and spikysphere (randomized radius at each point).  The
2d objects it supports are a circle, rectangle, triangle, and
spikycircle (randomized radius at each point).

See the header of the script for the syntax used to run it.


----------


.. _surftransform:

.. raw:: html

   <span id="surftransform"></span>

surf\_transform tool
--------------------------------------------------------------------------------------------

This is a Python script that transforms a SPARTA surface file into a
new surface file using various operations supported by the
:doc:`read\_surf <read_surf>` command.  These operations include
translation, scaling, rotation, and inversion (changing which side of
the surface is inside vs outside).

See the header of the script for the syntax used to run it.


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
