.. index:: dump particle/vtk

dump particle/vtk command
=========================

dump grid/vtk command
=====================

dump surf/vtk command
=====================

**Syntax:**


.. parsed-literal::

   dump ID particle/vtk mix-ID N file args

   dump ID grid/vtk group-ID N file args

   dump ID surf/vtk group-ID N file args

* ID = user-assigned name for the dump
* particle/vtk or grid/vtk or surf/vtk = style of dump
* mix-ID = mixture ID for particle/vtk (which particles to dump)
* group-ID = grid or surface group ID for grid/vtk or surf/vtk
* N = dump every this many timesteps
* file = name of file to write dump info to
* args = list of attributes, same as for the corresponding :doc:`dump <dump>` style

**Examples:**


.. parsed-literal::

   dump 1 particle/vtk 100 dump.\*.vtu all x y z vx vy vz

   dump 2 grid/vtk 100 tmp.grid.\*.vtu all id proc vol

   dump 3 surf/vtk 100 tmp.surf.%.\*.vtp all id type

**Description:**

These dump styles are the VTK-format analogs of the :doc:`dump particle <dump>`, :doc:`dump grid <dump>`, and :doc:`dump surf <dump>`
styles.  They periodically write a snapshot of geometry and per-element
data in a native VTK file format that can be visualized directly in
ParaView or VisIt, without the offline conversion done by the scripts
in the tools/paraview directory.

They are part of the VTK package and require SPARTA to be built with an
external VTK library.  See the `Section 2.2 <Section_packages.html#VTK>`_ doc page for how to build SPARTA with
the VTK package.

The geometry written for each style is:

* particle/vtk = each selected particle is written as a VTK vertex (point)
* grid/vtk = each owned in-group grid cell is written as a VTK voxel (3d)   or pixel (2d) in an unstructured grid
* surf/vtk = each owned in-group surface element is written as a VTK   triangle (3d) or line (2d)

The list of attributes *args* is identical to the corresponding :doc:`dump particle <dump>`, :doc:`dump grid <dump>`, or :doc:`dump surf <dump>`
command, including *c\_ID*, *f\_ID*, *v\_name*, and custom-attribute
keywords.  Each requested attribute is written as a named VTK data
array: point data for particle/vtk, cell data for grid/vtk and
surf/vtk.

For dump particle/vtk the *x*\ , *y*\ , and *z* attributes are required;
they define the point coordinates.  Consecutive *vx vy vz* and *xs ys
zs* attributes are grouped into 3-component VTK vector arrays.

For dump grid/vtk and dump surf/vtk the element geometry (cell corner
coordinates or element vertices) is generated automatically from the
grid/surf data structures, so it does not need to be listed as an
attribute.

**File formats and names:**

The output file format is selected by the filename extension:

* .vtk = legacy VTK format (single file, ASCII or binary)
* .vtp = XML PolyData (particle/vtk and surf/vtk only)
* .vtu = XML UnstructuredGrid

Because VTK writes one file per snapshot, the filename must contain a
"\*" wildcard, which is replaced with the timestep (see :doc:`dump <dump>`
for "\*" and padding rules).

If the filename also contains a "%" character, each processor writes
its own piece file (with "%" replaced by the processor or cluster ID)
and processor 0 writes a parallel summary file (.pvtp for .vtp pieces,
.pvtu for .vtu pieces).  This requires the XML formats; the legacy
.vtk format does not support "%".  Grid cells are volumetric, so dump
grid/vtk supports only the .vtu / .pvtu (and legacy .vtk) formats, not
.vtp.

**Dump modify options:**

The :doc:`dump\_modify <dump_modify>` *binary* keyword is specific to
these styles:


.. parsed-literal::

   dump_modify ID binary yes

writes the VTK data in binary (yes) or ASCII (no, the default) form.

For dump particle/vtk the :doc:`dump\_modify <dump_modify>` *region* and
*thresh* keywords are also supported, exactly as for :doc:`dump particle <dump>`, to restrict which particles are written.

**Restrictions:**

These styles are part of the VTK package.  They are only enabled if
SPARTA was built with that package and an external VTK library.  See
the `Section 2.2 <Section_packages.html#VTK>`_ doc page for details.  The
VTK package must be built with CMake.

For dump grid/vtk, split grid cells (cells cut by a surface into
multiple sub-cells) are written as multiple VTK cells that share the
same bounding box, since each sub-cell occupies the same box as its
parent.  These coincident voxels/pixels overlap in a viewer; use the
per-cell data (e.g. vol or a per-sub-cell compute) to distinguish
them, or restrict the dump to a grid group that excludes split cells.

Dump particle/vtk and dump surf/vtk do not support string attributes.
Dump grid/vtk writes the *idstr* attribute as a VTK string array.

**Related commands:**

:doc:`dump <dump>`, :doc:`dump\_modify <dump_modify>`, :doc:`dump image <dump_image>`

**Default:**

binary no (ASCII).


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
