12. Errors
==========

This section describes the various kinds of errors you can encounter
when using SPARTA.

| 12.1 `Common problems <#err_1>`_
| 12.2 `Reporting bugs <#err_2>`_
| 12.3 `Error & warning messages <#err_3>`_ 
| 


----------


.. _err_1:

.. raw:: html

   <span id="err_1"></span>

12.1 Common problems
--------------------

If two SPARTA runs do not produce the same answer on different
machines or different numbers of processors, this is typically not a
bug.  On different machines, there can be numerical round-off in the
computations which causes slight differences in particle trajectories
or the number of particles, which will lead to numerical divergence of
the particle trajectores and averaged statistical quantities within a
few 100s or few 1000s of timesteps.  When running on different numbers
of processors, random numbers are used in different ways, so two
simulations can be immediately different.  However, the statistical
properties (e.g. overall particle temperature or per grid cell
temperature or surface energy flux) for the two runs on different
machines or on different numbers of processors should still be
similar.

A SPARTA simulation typically has two stages, setup and run.  Most
SPARTA errors are detected at setup time; others like running out of
memory may not occur until the middle of a run.

SPARTA tries to flag errors and print informative error messages so
you can fix the problem.  Of course, SPARTA cannot figure out physics
or numerical mistakes, like choosing too big a timestep or specifying
erroneous collision parameters.  If you run into errors that SPARTA
doesn't catch that you think it should flag, please send an email to
the `developers <https://sparta.github.io/authors.html>`_.

If you get an error message about an invalid command in your input
script, you can determine what command is causing the problem by
looking in the log.sparta file, or using the :doc:`echo command <echo>`
in your script or "-echo screen" as a `command-line argument <Section_start.html#start_7>`_ to see it on the screen.  For a
given command, SPARTA expects certain arguments in a specified order.
If you mess this up, SPARTA will often flag the error, but it may read
a bogus argument and assign a value that is valid, but not what you
wanted.

Generally, SPARTA will print a message to the screen and logfile and
exit gracefully when it encounters a fatal error.  Sometimes it will
print a WARNING to the screen and logfile and continue on; you can
decide if the WARNING is important or not.  A WARNING message that is
generated in the middle of a run is only printed to the screen, not to
the logfile, to avoid cluttering up statistical output.  If SPARTA
crashes or hangs without spitting out an error message first then it
could be a bug (see the `next section <#err_2>`_) or one of the following
cases:

SPARTA runs in the available memory a processor allows to be
allocated.  Most reasonable runs are compute limited, not memory
limited, so this shouldn't be a bottleneck on most platforms.  Almost
all large memory allocations in the code are done via C-style malloc's
which will generate an error message if you run out of memory.
Smaller chunks of memory are allocated via C++ "new" statements.  If
you are unlucky, you could run out of memory just when one of these
small requests is made, in which case the code will crash or hang (in
parallel), since SPARTA doesn't trap on those errors.

Illegal arithmetic can cause SPARTA to run slow or crash.  This is
typically due to invalid physics and numerics that your simulation is
computing.  If you see wild statistical values or NaN values in your
SPARTA output, something is wrong with your simulation.  If you
suspect this is happening, it is a good idea to print out statistical
info frequently (e.g. every timestep) via the :doc:`stats <stats>`
command so you can monitor what is happening.  Visualizing the
particle motion is also a good idea to insure your model is behaving
as you expect.

In parallel, one way SPARTA can hang is due to how different MPI
implementations handle buffering of messages.  If the code hangs
without an error message, it may be that you need to specify an MPI
setting or two (usually via an environment variable) to enable
buffering or boost the sizes of messages that can be buffered.


----------


.. _err_2:

.. raw:: html

   <span id="err_2"></span>

12.2 Reporting bugs
-------------------

If you are confident that you have found a bug in SPARTA, please
follow these steps.

Check the `New features and bug fixes <https://sparta.github.io/bug.html>`_ section of the `SPARTA web site <sws_>`_ to see if the bug has already been fixed.

If not, please email a description of the problem to the
`developers <https://sparta.github.io/authors.html>`_.

The most useful thing you can do to help us fix the bug is to isolate
the problem.  Run it on the smallest number of particles and grid
cells and fewest number of processors and with the simplest and
quick-to-run input script that reproduces the bug.  And try to
identify what command or combination of commands is causing the
problem.


----------


.. _err_3:

.. raw:: html

   <span id="err_3"></span>

12.3 Error & warning messages
-------------------------------------------------------------------------------------

These are two alphabetic lists of the `ERROR <#error>`_ and
`WARNING <#warn>`_ messages SPARTA prints out and the reason why.  If the
explanation here is not sufficient, the documentation for the
offending command may help.  Error and warning messages also list the
source file and line number where the error was generated.  For
example, this message

ERROR: Illegal create\_particles command (create\_particles.cpp:68)

means that line #68 in the file src/create\_particles.cpp generated the
error.  Looking in the source code may help you figure out what went
wrong.

.. raw:: html

   <span id="error"></span>

Errors:
---------------------------------------------------------------



*%d read\_surf point pairs are too close*
   A pair of points is very close together, relative to grid size,
   inidicating the grid is too large, or an ill-formed surface.

*%d read\_surf points are not inside simulation box*
   If clipping was not performed, all points in surf file
   must be inside (or on surface of) simulation box.

*%d surface elements not assigned to a collision model*
   All surface elements must be assigned to a surface collision model via
   the surf\_modify command before a simulation is perforemd.

*All universe/uloop variables must have same # of values*
   Self-explanatory.

*All variables in next command must be same style*
   Self-explanatory.

*Arccos of invalid value in variable formula*
   Argument of arccos() must be between -1 and 1.

*Arcsin of invalid value in variable formula*
   Argument of arcsin() must be between -1 and 1.

*Axi-symmetry is not yet supported in SPARTA*
   This error condition will be removed after axi-symmetry is
   fully implemented.

*Axi-symmetry only allowed for 2d simulation*
   Self-explanatory.

*BPG edge on more than 2 faces*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Bad grid of processors for balance\_grid block*
   Product of Px,Py,Pz must equal total number of processors.

*Bad grid of processors for create\_grid*
   For block style, product of Px,Py,Pz must equal total number of
   processors.

*Bigint setting in spatype.h is invalid*
   Size of bigint is less than size of smallint.

*Bigint setting in spatype.h is not compatible*
   Bigint size stored in restart file is not consistent with SPARTA
   version you are running.

*Both restart files must use % or neither*
   Self-explanatory.

*Both sides of boundary must be periodic*
   Cannot specify a boundary as periodic only on the lo or hi side.  Must
   be periodic on both sides.

*Bound\_modify surf requires wall be a surface*
   The box boundary must be of style "s" to be assigned a surface
   collision model.

*Bound\_modify surf\_collide ID is unknown*
   Self-explanatory.

*Boundary command after simulation box is defined*
   The boundary command cannot be used after a read\_data, read\_restart,
   or create\_box command.

*Box boundary not assigned a surf\_collide ID*
   Any box boundary of style "s" must be assigned to a surface collision
   model via the bound\_modify command, before a simulation is performed.

*Box bounds are invalid*
   The box boundaries specified in the read\_data file are invalid.  The
   lo value must be less than the hi value for all 3 dimensions.

*Box ylo must be 0.0 for axi-symmetric model*
   Self-explanatory.

*Can only use -plog with multiple partitions*
   Self-explanatory.  See doc page discussion of command-line switches.

*Can only use -pscreen with multiple partitions*
   Self-explanatory.  See doc page discussion of command-line switches.

*Cannot add new species to mixture all or species*
   This is done automatically for these 2 mixtures when
   each species is defined by the species command.

*Cannot balance grid before grid is defined*
   Self-explanatory.

*Cannot create grid before simulation box is defined*
   Self-explanatory.

*Cannot create grid when grid is already defined*
   Self-explanatory.

*Cannot create particles before grid is defined*
   Self-explanatory.

*Cannot create particles before simulation box is defined*
   Self-explanatory.

*Cannot create/grow a vector/array of pointers for %s*
   SPARTA code is making an illegal call to the templated memory
   allocaters, to create a vector or array of pointers.

*Cannot create\_box after simulation box is defined*
   A simulation box can only be defined once.

*Cannot open VSS parameter file %s*
   Self-explantory.

*Cannot open dir to search for restart file*
   Using a "\*" in the name of the restart file will open the current
   directory to search for matching file names.

*Cannot open dump file*
   The output file for the dump command cannot be opened.  Check that the
   path and name are correct.

*Cannot open file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct. If the file is a compressed file, also check that the gzip
   executable can be found and run.

*Cannot open file variable file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix ave/time file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix print file %s*
   The output file generated by the fix print command cannot be opened

*Cannot open gzipped file*
   SPARTA was compiled without support for reading and writing gzipped
   files through a pipeline to the gzip program with -DSPARTA\_GZIP.

*Cannot open input script %s*
   Self-explanatory.

*Cannot open log.sparta*
   The default SPARTA log file cannot be opened.  Check that the
   directory you are running in allows for files to be created.

*Cannot open logfile*
   The SPARTA log file named in a command-line argument cannot be opened.
   Check that the path and name are correct.

*Cannot open logfile %s*
   The SPARTA log file specified in the input script cannot be opened.
   Check that the path and name are correct.

*Cannot open print file %s*
   Self-explanatory.

*Cannot open reaction file %s*
   Self-explanatory.

*Cannot open restart file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.  If the file is a compressed file, also check that the gzip
   executable can be found and run.

*Cannot open screen file*
   The screen file specified as a command-line argument cannot be
   opened.  Check that the directory you are running in allows for files
   to be created.

*Cannot open species file %s*
   Self-explanatory.

*Cannot open universe log file*
   For a multi-partition run, the master log file cannot be opened.
   Check that the directory you are running in allows for files to be
   created.

*Cannot open universe screen file*
   For a multi-partition run, the master screen file cannot be opened.
   Check that the directory you are running in allows for files to be
   created.

*Cannot read grid before simulation box is defined*
   Self-explanatory.

*Cannot read grid when grid is already defined*
   Self-explanatory.

*Cannot read\_restart after simulation box is defined*
   The read\_restart command cannot be used after a read\_data,
   read\_restart, or create\_box command.

*Cannot read\_surf after particles are defined*
   This is because the newly read surface objects may enclose particles.

*Cannot read\_surf before grid ghost cells are defined*
   This needs to be documented if keep this restriction.

*Cannot read\_surf before grid is defined*
   Self-explantory.

*Cannot redefine variable as a different style*
   An equal-style variable can be re-defined but only if it was
   originally an equal-style variable.

*Cannot reset timestep with a time-dependent fix defined*
   The timestep cannot be reset when a fix that keeps track of elapsed
   time is in place.

*Cannot run 2d simulation with nonperiodic Z dimension*
   Use the boundary command to make the z dimension periodic in order to
   run a 2d simulation.

*Cannot set global surfmax when surfaces already exist*
   This setting must be made before any surfac elements are
   read via the read\_surf command.

*Cannot use collide\_modify with no collisions defined*
   A collision style must be specified first.

*Cannot use cwiggle in variable formula between runs*
   This is a function of elapsed time.

*Cannot use dump\_modify fileper without % in dump file name*
   Self-explanatory.

*Cannot use dump\_modify nfile without % in dump file name*
   Self-explanatory.

*Cannot use fix inflow in y dimension for axisymmetric*
   This is because the y dimension boundaries cannot be
   inflow boundaries for an axisymmetric model.

*Cannot use fix inflow in z dimension for 2d simulation*
   Self-explanatory.

*Cannot use fix inflow n > 0 with perspecies yes*
   This is because the perspecies option calculates the
   number of particles to insert itself.

*Cannot use fix inflow on periodic boundary*
   Self-explanatory.

*Cannot use group keyword with mixture all or species*
   This is because the groups for these 2 mixtures are
   pre-defined.

*Cannot use include command within an if command*
   Self-explanatory.

*Cannot use non-rcb fix balance with a grid cutoff*
   This is because the load-balancing will generate a partitioning
   of cells to processors that is dispersed and which will not work
   with a grid cutoff >= 0.0.

*Cannot use ramp in variable formula between runs*
   This is because the ramp() function is time dependent.

*Cannot use specified create\_grid options with more than one level*
   When defining a grid with more than one level, the other create\_grid
   keywords (stride, clump, block, etc) cannot be used.  The child grid
   cells will be assigned to processors in round-robin order as explained
   on the create\_grid doc page.

*Cannot use swiggle in variable formula between runs*
   This is a function of elapsed time.

*Cannot use vdisplace in variable formula between runs*
   This is a function of elapsed time.

*Cannot use weight cell radius unless axisymmetric*
   An axisymmetric model is required for this style of cell weighting.

*Cannot use write\_restart fileper without % in restart file name*
   Self-explanatory.

*Cannot use write\_restart nfile without % in restart file name*
   Self-explanatory.

*Cannot weight cells before grid is defined*
   Self-explanatory.

*Cannot write grid when grid is not defined*
   Self-explanatory.

*Cannot write restart file before grid is defined*
   Self-explanatory.

*Cell ID has too many bits*
   Cell IDs must fit in 32 bits (SPARTA small integer) or 64 bits (SPARTA
   big integer), as specified by the -DSPARTA\_SMALL, -DSPARTA\_BIG, or
   -DSPARTA\_BIGBIG options in the low-level Makefile used to build
   SPARTA.  See Section 2.2 of the manual for details.  And see Section
   4.8 for details on how cell IDs are formatted.

*Cell type mis-match when marking on neigh proc*
   Grid cell marking as inside, outside, or overlapping with surface
   elements failed.  Please report the issue to the SPARTA developers.

*Cell type mis-match when marking on self*
   Grid cell marking as inside, outside, or overlapping with surface
   elements failed.  Please report the issue to the SPARTA developers.

*Cellint setting in spatype.h is not compatible*
   Cellint size stored in restart file is not consistent with SPARTA
   version you are running.

*Collision mixture does not contain all species*
   The specified mixture must contain all species in the simulation so
   that they can be assigned to collision groups.

*Collision mixture does not exist*
   Self-explantory.

*Compute ID for compute reduce does not exist*
   Self-explanatory.

*Compute ID for fix ave/grid does not exist*
   Self-explanatory.

*Compute ID for fix ave/surf does not exist*
   Self-explanatory.

*Compute ID for fix ave/time does not exist*
   Self-explanatory.

*Compute ID must be alphanumeric or underscore characters*
   Self-explanatory.

*Compute boundary mixture ID does not exist*
   Self-explanatory.

*Compute grid mixture ID does not exist*
   Self-explanatory.

*Compute reduce compute array is accessed out-of-range*
   An index for the array is out of bounds.

*Compute reduce compute calculates global or surf values*
   The compute reduce command does not operate on this kind of values.
   The variable command has special functions that can reduce global
   values.

*Compute reduce compute does not calculate a per-grid array*
   This is necessary if a column index is used to specify the compute.

*Compute reduce compute does not calculate a per-grid vector*
   This is necessary if no column index is used to specify the compute.

*Compute reduce compute does not calculate a per-particle array*
   This is necessary if a column index is used to specify the compute.

*Compute reduce compute does not calculate a per-particle vector*
   This is necessary if no column index is used to specify the compute.

*Compute reduce fix array is accessed out-of-range*
   An index for the array is out of bounds.

*Compute reduce fix calculates global values*
   A fix that calculates peratom or local values is required.

*Compute reduce fix does not calculate a per-grid array*
   This is necessary if a column index is used to specify the fix.

*Compute reduce fix does not calculate a per-grid vector*
   This is necessary if no column index is used to specify the fix.

*Compute reduce fix does not calculate a per-particle array*
   This is necessary if a column index is used to specify the fix.

*Compute reduce fix does not calculate a per-particle vector*
   This is necessary if no column index is used to specify the fix.

*Compute reduce fix does not calculate a per-surf array*
   This is necessary if a column index is used to specify the fix.

*Compute reduce fix does not calculate a per-surf vector*
   This is necessary if no column index is used to specify the fix.

*Compute reduce replace requires min or max mode*
   Self-explanatory.

*Compute reduce variable is not particle-style variable*
   This is the only style of variable that can be reduced.

*Compute sonine/grid mixture ID does not exist*
   Self-explanatory.

*Compute surf mixture ID does not exist*
   Self-explanatory.

*Compute used in variable between runs is not current*
   Computes cannot be invoked by a variable in between runs.  Thus they
   must have been evaluated on the last timestep of the previous run in
   order for their value(s) to be accessed.  See the doc page for the
   variable command for more info.

*Could not create a single particle*
   The specified position was either not inside the simulation domain or
   not inside a grid cell with no intersections with any defined surface
   elements.

*Could not find compute ID to delete*
   Self-explanatory.

*Could not find dump grid compute ID*
   Self-explanatory.

*Could not find dump grid fix ID*
   Self-explanatory.

*Could not find dump grid variable name*
   Self-explanatory.

*Could not find dump image compute ID*
   Self-explanatory.

*Could not find dump image fix ID*
   Self-explanatory.

*Could not find dump modify compute ID*
   Self-explanatory.

*Could not find dump modify fix ID*
   Self-explanatory.

*Could not find dump modify variable name*
   Self-explanatory.

*Could not find dump particle compute ID*
   Self-explanatory.

*Could not find dump particle fix ID*
   Self-explanatory.

*Could not find dump particle variable name*
   Self-explanatory.

*Could not find dump surf compute ID*
   Self-explanatory.

*Could not find dump surf fix ID*
   Self-explanatory.

*Could not find dump surf variable name*
   Self-explanatory.

*Could not find fix ID to delete*
   Self-explanatory.

*Could not find split point in split cell*
   This is an error when calculating how a grid cell is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Could not find stats compute ID*
   Compute ID specified in stats\_style command does not exist.

*Could not find stats fix ID*
   Fix ID specified in stats\_style command does not exist.

*Could not find stats variable name*
   Self-explanatory.

*Could not find surf\_modify sc-ID*
   Self-explanatory.

*Could not find surf\_modify surf-ID*
   Self-explanatory.

*Could not find undump ID*
   A dump ID used in the undump command does not exist.

*Cound not find dump\_modify ID*
   Self-explanatory.

*Create\_box z box bounds must straddle 0.0 for 2d simulations*
   Self-explanatory.

*Create\_grid nz value must be 1 for a 2d simulation*
   Self-explanatory.

*Create\_particles global option not yet implemented*
   Self-explantory.

*Create\_particles mixture ID does not exist*
   Self-explanatory.

*Create\_particles single requires z = 0 for 2d simulation*
   Self-explanatory.

*Create\_particles species ID does not exist*
   Self-explanatory.

*Created incorrect # of particles: %ld versus %ld*
   The create\_particles command did not function
   properly.

*Delete region ID does not exist*
   Self-explanatory.

*Did not assign all restart particles correctly*
   One or more particles in the restart file were not assigned to a
   processor.  Please report the issue to the SPARTA developers.

*Did not assign all restart split grid cells correctly*
   One or more split grid cells in the restart file were not assigned
   to a processor.  Please report the issue to the SPARTA developers.

*Did not assign all restart sub grid cells correctly*
   One or more sub grid cells in the restart file were not assigned to a
   processor.  Please report the issue to the SPARTA developers.

*Did not assign all restart unsplit grid cells correctly*
   One or more unsplit grid cells in the restart file were not assigned
   to a processor.  Please report the issue to the SPARTA developers.

*Dimension command after simulation box is defined*
   The dimension command cannot be used after a read\_data,
   read\_restart, or create\_box command.

*Divide by 0 in variable formula*
   Self-explanatory.

*Dump every variable returned a bad timestep*
   The variable must return a timestep greater than the current timestep.

*Dump grid and fix not computed at compatible times*
   Fixes generate values on specific timesteps.  The dump grid output
   does not match these timesteps.

*Dump grid compute does not calculate per-grid array*
   Self-explanatory.

*Dump grid compute does not compute per-grid info*
   Self-explanatory.

*Dump grid compute vector is accessed out-of-range*
   Self-explanatory.

*Dump grid fix does not compute per-grid array*
   Self-explanatory.

*Dump grid fix does not compute per-grid info*
   Self-explanatory.

*Dump grid fix vector is accessed out-of-range*
   Self-explanatory.

*Dump grid variable is not grid-style variable*
   Self-explanatory.

*Dump image and fix not computed at compatible times*
   Fixes generate values on specific timesteps.  The dump image output
   does not match these timesteps.

*Dump image cannot use grid and gridx/gridy/gridz*
   Can only use grid option or one or more of grid x,y,z options
   by themselves, not together.

*Dump image compute does not have requested column*
   Self-explanatory.

*Dump image compute does not produce a vector*
   Self-explanatory.

*Dump image compute is not a per-grid compute*
   Self-explanatory.

*Dump image compute is not a per-surf compute*
   Self-explanatory.

*Dump image fix does not have requested column*
   Self-explanatory.

*Dump image fix does not produce a vector*
   Self-explanatory.

*Dump image fix does not produce per-grid values*
   Self-explanatory.

*Dump image fix does not produce per-surf values*
   Self-explanatory.

*Dump image requires one snapshot per file*
   Use a "\*" in the filename.

*Dump modify compute ID does not compute per-particle array*
   Self-explanatory.

*Dump modify compute ID does not compute per-particle info*
   Self-explanatory.

*Dump modify compute ID does not compute per-particle vector*
   Self-explanatory.

*Dump modify compute ID vector is not large enough*
   Self-explanatory.

*Dump modify fix ID does not compute per-particle array*
   Self-explanatory.

*Dump modify fix ID does not compute per-particle info*
   Self-explanatory.

*Dump modify fix ID does not compute per-particle vector*
   Self-explanatory.

*Dump modify fix ID vector is not large enough*
   Self-explanatory.

*Dump modify variable is not particle-style variable*
   Self-explanatory.

*Dump particle and fix not computed at compatible times*
   Fixes generate values on specific timesteps.  The dump particle output
   does not match these timesteps.

*Dump particle compute does not calculate per-particle array*
   Self-explanatory.

*Dump particle compute does not calculate per-particle vector*
   Self-explanatory.

*Dump particle compute does not compute per-particle info*
   Self-explanatory.

*Dump particle compute vector is accessed out-of-range*
   Self-explanatory.

*Dump particle fix does not compute per-particle array*
   Self-explanatory.

*Dump particle fix does not compute per-particle info*
   Self-explanatory.

*Dump particle fix does not compute per-particle vector*
   Self-explanatory.

*Dump particle fix vector is accessed out-of-range*
   Self-explanatory.

*Dump particle variable is not particle-style variable*
   Self-explanatory.

*Dump surf and fix not computed at compatible times*
   Fixes generate values on specific timesteps.  The dump surf output
   does not match these timesteps.

*Dump surf compute does not calculate per-surf array*
   Self-explanatory.

*Dump surf compute does not compute per-surf info*
   Self-explanatory.

*Dump surf compute vector is accessed out-of-range*
   Self-explanatory.

*Dump surf fix does not compute per-surf array*
   Self-explanatory.

*Dump surf fix does not compute per-surf info*
   Self-explanatory.

*Dump surf fix vector is accessed out-of-range*
   Self-explanatory.

*Dump surf variable is not surf-style variable*
   Self-explanatory.

*Dump\_modify buffer yes not allowed for this style*
   Not all dump styles allow dump\_modify buffer yes.  See the dump\_modify
   doc page.

*Dump\_modify region ID does not exist*
   Self-explanatory.

*Duplicate cell ID in grid file*
   Parent cell IDs must be unique.

*Edge not part of 2 vertices*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Edge part of invalid vertex*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Edge part of same vertex twice*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Empty brackets in variable*
   There is no variable syntax that uses empty brackets.  Check
   the variable doc page.

*Failed to allocate %ld bytes for array %s*
   The SPARTA simulation has run out of memory.  You need to run a
   smaller simulation or on more processors.

*Failed to open FFmpeg pipeline to file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct and writable and that the FFmpeg executable can be found and run.

*Failed to reallocate %ld bytes for array %s*
   The SPARTA simulation has run out of memory.  You need to run a
   smaller simulation or on more processors.

*File variable could not read value*
   Check the file assigned to the variable.

*Fix ID for compute reduce does not exist*
   Self-explanatory.

*Fix ID for fix ave/grid does not exist*
   Self-explanatory.

*Fix ID for fix ave/surf does not exist*
   Self-explanatory.

*Fix ID for fix ave/time does not exist*
   Self-explanatory.

*Fix ID must be alphanumeric or underscore characters*
   Self-explanatory.

*Fix ave/grid compute array is accessed out-of-range*
   Self-explanatory.

*Fix ave/grid compute does not calculate a per-grid array*
   Self-explanatory.

*Fix ave/grid compute does not calculate a per-grid vector*
   Self-explanatory.

*Fix ave/grid compute does not calculate per-grid values*
   Self-explanatory.

*Fix ave/grid fix array is accessed out-of-range*
   Self-explanatory.

*Fix ave/grid fix does not calculate a per-grid array*
   Self-explanatory.

*Fix ave/grid fix does not calculate a per-grid vector*
   Self-explanatory.

*Fix ave/grid fix does not calculate per-grid values*
   Self-explanatory.

*Fix ave/grid variable is not grid-style variable*
   Self-explanatory.

*Fix ave/surf compute array is accessed out-of-range*
   Self-explanatory.

*Fix ave/surf compute does not calculate a per-surf array*
   Self-explanatory.

*Fix ave/surf compute does not calculate a per-surf vector*
   Self-explanatory.

*Fix ave/surf compute does not calculate per-surf values*
   Self-explanatory.

*Fix ave/surf fix array is accessed out-of-range*
   Self-explanatory.

*Fix ave/surf fix does not calculate a per-surf array*
   Self-explanatory.

*Fix ave/surf fix does not calculate a per-surf vector*
   Self-explanatory.

*Fix ave/surf fix does not calculate per-surf values*
   Self-explanatory.

*Fix ave/surf variable is not surf-style variable*
   Self-explanatory.

*Fix ave/time cannot use variable with vector mode*
   Variables produce scalar values.

*Fix ave/time columns are inconsistent lengths*
   Self-explanatory.

*Fix ave/time compute array is accessed out-of-range*
   An index for the array is out of bounds.

*Fix ave/time compute does not calculate a scalar*
   Self-explantory.

*Fix ave/time compute does not calculate a vector*
   Self-explantory.

*Fix ave/time compute does not calculate an array*
   Self-explanatory.

*Fix ave/time compute vector is accessed out-of-range*
   The index for the vector is out of bounds.

*Fix ave/time fix array is accessed out-of-range*
   An index for the array is out of bounds.

*Fix ave/time fix does not calculate a scalar*
   Self-explanatory.

*Fix ave/time fix does not calculate a vector*
   Self-explanatory.

*Fix ave/time fix does not calculate an array*
   Self-explanatory.

*Fix ave/time fix vector is accessed out-of-range*
   The index for the vector is out of bounds.

*Fix ave/time variable is not equal-style variable*
   Self-explanatory.

*Fix command before simulation box is defined*
   The fix command cannot be used before a read\_data, read\_restart, or
   create\_box command.

*Fix for fix ave/grid not computed at compatible time*
   Fixes generate values on specific timesteps.  Fix ave/grid is
   requesting a value on a non-allowed timestep.

*Fix for fix ave/surf not computed at compatible time*
   Fixes generate their values on specific timesteps.  Fix ave/surf is
   requesting a value on a non-allowed timestep.

*Fix for fix ave/time not computed at compatible time*
   Fixes generate their values on specific timesteps.  Fix ave/time
   is requesting a value on a non-allowed timestep.

*Fix in variable not computed at compatible time*
   Fixes generate their values on specific timesteps.  The variable is
   requesting the values on a non-allowed timestep.

*Fix inflow mixture ID does not exist*
   Self-explanatory.

*Fix inflow used on outflow boundary*
   Self-explanatory.

*Fix used in compute reduce not computed at compatible time*
   Fixes generate their values on specific timesteps.  Compute reduce is
   requesting a value on a non-allowed timestep.

*Found edge in same direction*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Found no restart file matching pattern*
   When using a "\*" in the restart file name, no matching file was found.

*Gravity in y not allowed for axi-symmetric model*
   Self-explanatory.

*Gravity in z not allowed for 2d*
   Self-explanatory.

*Grid cell corner points on boundary marked as unknown = %d*
   Corner points of grid cells on the boundary of the simulation domain
   were not all marked successfully as inside, outside, or overlapping
   with surface elements.  Please report the issue to the SPARTA
   developers.

*Grid cells marked as unknown = %d*
   Grid cell marking as inside, outside, or overlapping with surface
   elements did not successfully mark all cells.  Please report the issue
   to the SPARTA developers.

*Grid cutoff is longer than box length in a periodic dimension*
   This is not allowed.  Reduce the size of the cutoff specified by the
   global gridcut command.

*Grid in/out other-mark error %d\\n*
   Grid cell marking as inside, outside, or overlapping with surface
   elements failed.  Please report the issue to the SPARTA developers.

*Grid in/out self-mark error %d for icell %d, icorner %d, connect %d %d, other cell %d, other corner %d, values %d %d\\n*
   A grid cell was incorrectly marked as inside, outside, or overlapping
   with surface elements.  Please report the issue to the SPARTA
   developers.

*Grid-style variables are not yet implemented*
   Self-explanatory.

*Illegal ... command*
   Self-explanatory.  Check the input script syntax and compare to the
   documentation for the command.  You can use -echo screen as a
   command-line option when running SPARTA to see the offending line.

*Inconsistent surface to grid mapping in read\_restart*
   When surface elements were mapped to grid cells after reading a
   restart file, an inconsitent count of elements in a grid cell was
   found, as compared to the original simulation, which should not
   happen.  Please report the issue to the SPARTA developers.

*Incorrect format of parent cell in grid file*
   Number of words in a parent cell line was not the expected number.

*Incorrect line format in VSS parameter file*
   Number of parameters in a line read from file is not valid.

*Incorrect line format in species file*
   Line read did not have expected number of fields.

*Incorrect line format in surf file*
   Self-explanatory.

*Incorrect point format in surf file*
   Self-explanatory.

*Incorrect triangle format in surf file*
   Self-explanatory.

*Index between variable brackets must be positive*
   Self-explanatory.

*Input line quote not followed by whitespace*
   An end quote must be followed by whitespace.

*Invalid Boolean syntax in if command*
   Self-explanatory.

*Invalid Nx,Ny,Nz values in grid file*
   A Nx or Ny or Nz value for a parent cell is <= 0.

*Invalid SPARTA restart file*
   The file does not appear to be a SPARTA restart file since it does not
   have the expected magic string at the beginning.

*Invalid attribute in dump grid command*
   Self-explanatory.

*Invalid attribute in dump modify command*
   Self-explantory.

*Invalid attribute in dump particle command*
   Self-explanatory.

*Invalid attribute in dump surf command*
   Self-explanatory.

*Invalid balance\_grid style for non-uniform grid*
   Some balance styles can only be used when the grid is uniform.  See
   the command doc page for details.

*Invalid call to ComputeGrid::post\_process\_grid()*
   This indicates a coding error.  Please report the issue to the SPARTA
   developers.

*Invalid call to ComputeSonineGrid::post\_process\_grid()*
   This indicates a coding error.  Please report the issue to the SPARTA
   developers.

*Invalid cell ID in grid file*
   A cell ID could not be converted into numeric format.

*Invalid character in species ID*
   The only allowed characters are alphanumeric, an underscore, a plus
   sign, or a minus sign.

*Invalid collide style*
   The choice of collision style is unknown.

*Invalid color in dump image command*
   The specified color name was not in the list of recognized colors.
   See the dump image doc page.

*Invalid color in dump\_modify command*
   The specified color name was not in the list of recognized colors.
   See the dump\_modify doc page.

*Invalid color map min/max values*
   The min/max values are not consistent with either each other or
   with values in the color map.

*Invalid command-line argument*
   One or more command-line arguments is invalid.  Check the syntax of
   the command you are using to launch SPARTA.

*Invalid compute ID in variable formula*
   The compute is not recognized.

*Invalid compute property/grid field for 2d simulation*
   Fields that reference z-dimension properties cannot be used
   in a 2d simulation.

*Invalid compute style*
   Self-explanatory.

*Invalid dump frequency*
   Dump frequency must be 1 or greater.

*Invalid dump grid field for 2d simulation*
   Self-explanatory.

*Invalid dump image filename*
   The file produced by dump image cannot be binary and must
   be for a single processor.

*Invalid dump image theta value*
   Theta must be between 0.0 and 180.0 inclusive.

*Invalid dump image zoom value*
   Zoom value must be > 0.0.

*Invalid dump movie filename*
   The file produced by dump movie cannot be binary or compressed
   and must be a single file for a single processor.

*Invalid dump style*
   The choice of dump style is unknown.

*Invalid dump surf field for 2d simulation*
   Self-explanatory.

*Invalid dump\_modify threshhold operator*
   Operator keyword used for threshold specification in not recognized.

*Invalid fix ID in variable formula*
   The fix is not recognized.

*Invalid fix ave/time off column*
   Self-explantory.

*Invalid fix style*
   The choice of fix style is unknown.

*Invalid flag in grid section of restart file*
   Unrecognized entry in restart file.

*Invalid flag in header section of restart file*
   Unrecognized entry in restart file.

*Invalid flag in layout section of restart file*
   Unrecognized entry in restart file.

*Invalid flag in particle section of restart file*
   Unrecognized entry in restart file.

*Invalid flag in peratom section of restart file*
   The format of this section of the file is not correct.

*Invalid flag in surf section of restart file*
   Unrecognized entry in restart file.

*Invalid image up vector*
   Up vector cannot be (0,0,0).

*Invalid immediate variable*
   Syntax of immediate value is incorrect.

*Invalid keyword in compute property/grid command*
   Self-explantory.

*Invalid keyword in stats\_style command*
   One or more specified keywords are not recognized.

*Invalid math function in variable formula*
   Self-explanatory.

*Invalid math/special function in variable formula*
   Self-explanatory.

*Invalid point index in line*
   Self-explanatory.

*Invalid point index in triangle*
   Self-explanatory.

*Invalid power expression in variable formula*
   A zero base cannot be raised to a negative power in a variable
   formula, since the result is infinite.

*Invalid react style*
   The choice of reaction style is unknown.

*Invalid reaction coefficients in file*
   Self-explanatory.

*Invalid reaction formula in file*
   Self-explanatory.

*Invalid reaction style in file*
   Self-explanatory.

*Invalid reaction type in file*
   Self-explanatory.

*Invalid read\_surf command*
   Self-explanatory.

*Invalid read\_surf geometry transformation for 2d simulation*
   Cannot perform a transformation that changes z cooridinates of points
   for a 2d simulation.

*Invalid region style*
   The choice of region style is unknown.

*Invalid replace values in compute reduce*
   Self-explanatory.

*Invalid reuse of surface ID in read\_surf command*
   Surface IDs must be unique.

*Invalid run command N value*
   The number of timesteps must fit in a 32-bit integer.  If you want to
   run for more steps than this, perform multiple shorter runs.

*Invalid run command start/stop value*
   Self-explanatory.

*Invalid run command upto value*
   Self-explanatory.

*Invalid special function in variable formula*
   Self-explanatory.

*Invalid species ID in species file*
   Species IDs are limited to 15 characters.

*Invalid stats keyword in variable formula*
   The keyword is not recognized.

*Invalid surf\_collide style*
   Self-explanatory.

*Invalid syntax in variable formula*
   Self-explanatory.

*Invalid use of library file() function*
   This function is called thru the library interface.  This
   error should not occur.  Contact the developers if it does.

*Invalid variable evaluation in variable formula*
   A variable used in a formula could not be evaluated.

*Invalid variable in next command*
   Self-explanatory.

*Invalid variable name*
   Variable name used in an input script line is invalid.

*Invalid variable name in variable formula*
   Variable name is not recognized.

*Invalid variable style in special function next*
   Only file-style or atomfile-style variables can be used with next().

*Invalid variable style with next command*
   Variable styles *equal* and *world* cannot be used in a next
   command.

*Ionization and recombination reactions are not yet implemented*
   This error conditions will be removed after those reaction styles are
   fully implemented.

*Irregular comm recv buffer exceeds 2 GB*
   MPI does not support a communication buffer that exceeds a 4-byte
   integer in size.

*Label wasn't found in input script*
   Self-explanatory.

*Log of zero/negative value in variable formula*
   Self-explanatory.

*MPI\_SPARTA\_BIGINT and bigint in spatype.h are not compatible*
   The size of the MPI datatype does not match the size of a bigint.

*Migrate cells send buffer exceeds 2 GB*
   MPI does not support a communication buffer that exceeds a 4-byte
   integer in size.

*Mismatched brackets in variable*
   Self-explanatory.

*Mismatched compute in variable formula*
   A compute is referenced incorrectly or a compute that produces per-atom
   values is used in an equal-style variable formula.

*Mismatched fix in variable formula*
   A fix is referenced incorrectly or a fix that produces per-atom
   values is used in an equal-style variable formula.

*Mismatched variable in variable formula*
   A variable is referenced incorrectly or an atom-style variable that
   produces per-atom values is used in an equal-style variable
   formula.

*Mixture %s fractions exceed 1.0*
   The sum of fractions must not be > 1.0.

*Mixture ID must be alphanumeric or underscore characters*
   Self-explanatory.

*Mixture group ID must be alphanumeric or underscore characters*
   Self-explanatory.

*Mixture species is not defined*
   One or more of the species ID is unknown.

*Modulo 0 in variable formula*
   Self-explanatory.

*More than one positive area with a negative area*
   SPARTA cannot determine which positive area the negative area is
   inside of, if a cell is so large that it includes both positive and
   negative areas.

*More than one positive volume with a negative volume*
   SPARTA cannot determine which positive volume the negative volume is
   inside of, if a cell is so large that it includes both positive and
   negative volumes.

*Must use -in switch with multiple partitions*
   A multi-partition simulation cannot read the input script from stdin.
   The -in command-line option must be used to specify a file.

*Next command must list all universe and uloop variables*
   This is to insure they stay in sync.

*No dump grid attributes specified*
   Self-explanatory.

*No dump particle attributes specified*
   Self-explanatory.

*No dump surf attributes specified*
   Self-explanatory.

*No positive areas in cell*
   This is an error when calculating how a 2d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*No positive volumes in cell*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Non digit character between brackets in variable*
   Self-explantory.

*Number of groups in compute boundary mixture has changed*
   This mixture property cannot be changed after this compute command is
   issued.

*Number of groups in compute grid mixture has changed*
   This mixture property cannot be changed after this compute command is
   issued.

*Number of groups in compute sonine/grid mixture has changed*
   This mixture property cannot be changed after this compute command is
   issued.

*Number of groups in compute surf mixture has changed*
   This mixture property cannot be changed after this compute command is
   issued.

*Number of groups in compute tvib/grid mixture has changed*
   This mixture property cannot be changed after this compute command is
   issued.

*Number of species in compute tvib/grid mixture has changed*
   This mixture property cannot be changed after this compute command is
   issued.

*Numeric index is out of bounds*
   A command with an argument that specifies an integer or range of
   integers is using a value that is less than 1 or greater than the
   maximum allowed limit.

*Nz value in read\_grid file must be 1 for a 2d simulation*
   Self-explanatory.

*Only ylo boundary can be axi-symmetric*
   Self-explanatory.  See the boundary doc page for more details.

*Owned cells with unknown neighbors = %d*
   One or more grid cells have unknown neighbors which will prevent
   particles from moving correctly.  Please report the issue to the
   SPARTA developers.

*Parent cell child missing*
   Hierarchical grid traversal failed.  Please report the issue to the
   SPARTA developers.

*Particle %d on proc %d hit inside of surf %d on step %ld*
   This error should not happen if particles start outside of physical
   objects.  Please report the issue to the SPARTA developers.

*Particle %d,%d on proc %d is in invalid cell  on timestep %ld*
   The particle is in a cell indexed by a value that is out-of-bounds for
   the cells owned by this processor.

*Particle %d,%d on proc %d is in split cell  on timestep %ld*
   This should not happend.  The particle should be in one
   of the sub-cells of the split cell.

*Particle %d,%d on proc %d is outside cell  on timestep %ld*
   The particle's coordinates are not within the grid cell
   it is supposed to be in.

*Particle vector in equal-style variable formula*
   Equal-style variables cannot use per-particle quantities.

*Particle-style variable in equal-style variable formula*
   Equal-style variables cannot use per-particle quantities.

*Partition numeric index is out of bounds*
   It must be an integer from 1 to the number of partitions.

*Per-particle compute in equal-style variable formula*
   Equal-style variables cannot use per-particle quantities.

*Per-particle fix in equal-style variable formula*
   Equal-style variables cannot use per-particle quantities.

*Per-processor particle count is too big*
   No processor can have more particle than fit in a 32-bit integer,
   approximately 2 billion.

*Point appears first in more than one CLINE*
   This is an error when calculating how a 2d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Point appears last in more than one CLINE*
   This is an error when calculating how a 2d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Processor partitions are inconsistent*
   The total number of processors in all partitions must match the number
   of processors SPARTA is running on.

*React tce can only be used with collide vss*
   Self-explanatory.

*Read\_grid did not find parents section of grid file*
   Expected Parents section but did not find keyword.

*Read\_surf did not find lines section of surf file*
   Expected Lines section but did not find keyword.

*Read\_surf did not find points section of surf file*
   Expected Parents section but did not find keyword.

*Read\_surf did not find triangles section of surf file*
   Expected Triangles section but did not find keyword.

*Region ID for dump custom does not exist*
   Self-explanatory.

*Region intersect region ID does not exist*
   One or more of the region IDs specified by the region intersect
   command does not exist.

*Region union region ID does not exist*
   One or more of the region IDs specified by the region union command
   does not exist.

*Replacing a fix, but new style != old style*
   A fix ID can be used a 2nd time, but only if the style matches the
   previous fix.  In this case it is assumed you with to reset a fix's
   parameters.  This error may mean you are mistakenly re-using a fix ID
   when you do not intend to.

*Request for unknown parameter from collide*
   VSS model does not have the parameter being requested.

*Restart file byte ordering is not recognized*
   The file does not appear to be a SPARTA restart file since it doesn't
   contain a recognized byte-ordering flag at the beginning.

*Restart file byte ordering is swapped*
   The file was written on a machine with different byte-ordering than
   the machine you are reading it on.

*Restart file incompatible with current version*
   This is probably because you are trying to read a file created with a
   version of SPARTA that is too old compared to the current version.

*Restart file is a multi-proc file*
   The file is inconsistent with the filename specified for it.

*Restart file is not a multi-proc file*
   The file is inconsistent with the filename specified for it.

*Restart variable returned a bad timestep*
   The variable must return a timestep greater than the current timestep.

*Reuse of compute ID*
   A compute ID cannot be used twice.

*Reuse of dump ID*
   A dump ID cannot be used twice.

*Reuse of region ID*
   A region ID cannot be used twice.

*Reuse of surf\_collide ID*
   A surface collision model ID cannot be used more than once.

*Run command before grid ghost cells are defined*
   Normally, ghost cells will be defined when the grid is created via the
   create\_grid or read\_grid commands.  However, if the global gridcut
   cutoff is set to a value >= 0.0, then ghost cells can only be defined
   if the partiioning of cells to processors is clumped, not dispersed.
   See the fix balance command for an explanation.  Invoking the fix
   balance command with a clumped option will trigger ghost cells to be
   defined.

*Run command before grid is defined*
   Self-explanatory.

*Run command start value is after start of run*
   Self-explanatory.

*Run command stop value is before end of run*
   Self-explanatory.

*Seed command has not been used*
   This command should appear near the beginning of your input script,
   before any random numbers are needed by other commands.

*Sending particle to self*
   This error should not occur.  Please report the issue to the SPARTA
   developers.

*Single area is negative, inverse donut*
   An inverse donut is a surface with a flow region interior to the donut
   hole and also exterior to the entire donut.  This means the flow
   regions are disconnected.  SPARTA cannot correctly compute the flow
   area of this kind of object.

*Single volume is negative, inverse donut*
   An inverse donut is a surface with a flow region interior to the donut
   hole and also exterior to the entire donut.  This means the flow
   regions are disconnected.  SPARTA cannot correctly compute the flow
   volume of this kind of object.

*Singlet BPG edge not on cell face*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Singlet CLINES point not on cell border*
   This is an error when calculating how a 2d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Small,big integers are not sized correctly*
   This error occurs whenthe sizes of smallint and bigint as defined in
   src/spatype.h are not what is expected.  Please report the issue
   to the SPARTA developers.

*Smallint setting in spatype.h is invalid*
   It has to be the size of an integer.

*Smallint setting in spatype.h is not compatible*
   Smallint size stored in restart file is not consistent with SPARTA
   version you are running.

*Species %s did not appear in VSS parameter file*
   Self-explanatory.

*Species ID does not appear in species file*
   Could not find the requested species in the specified file.

*Species ID is already defined*
   Species IDs must be unique.

*Sqrt of negative value in variable formula*
   Self-explanatory.

*Stats and fix not computed at compatible times*
   Fixes generate values on specific timesteps.  The stats output
   does not match these timesteps.

*Stats compute array is accessed out-of-range*
   Self-explanatory.

*Stats compute does not compute array*
   Self-explanatory.

*Stats compute does not compute scalar*
   Self-explanatory.

*Stats compute does not compute vector*
   Self-explanatory.

*Stats compute vector is accessed out-of-range*
   Self-explanatory.

*Stats every variable returned a bad timestep*
   The variable must return a timestep greater than the current timestep.

*Stats fix array is accessed out-of-range*
   Self-explanatory.

*Stats fix does not compute array*
   Self-explanatory.

*Stats fix does not compute scalar*
   Self-explanatory.

*Stats fix does not compute vector*
   Self-explanatory.

*Stats fix vector is accessed out-of-range*
   Self-explanatory.

*Stats variable cannot be indexed*
   A variable used as a stats keyword cannot be indexed.
   E.g. v\_foo must be used, not v\_foo\ **100**\ .

*Stats variable is not equal-style variable*
   Only equal-style variables can be output with stats output, not
   particle-style or grid-style or surf-style variables.

*Stats\_modify every variable returned a bad timestep*
   The variable must return a timestep greater than the current timestep.

*Stats\_modify int format does not contain d character*
   Self-explanatory.

*Substitution for illegal variable*
   Input script line contained a variable that could not be substituted
   for.

*Support for writing images in JPEG format not included*
   SPARTA was not built with the -DSPARTA\_JPEG switch in the Makefile.

*Support for writing images in PNG format not included*
   SPARTA was not built with the -DSPARTA\_PNG switch in the Makefile.

*Support for writing movies not included*
   SPARTA was not built with the -DSPARTA\_FFMPEG switch in the Makefile

*Surf file cannot contain lines for 3d simulation*
   Self-explanatory.

*Surf file cannot contain triangles for 2d simulation*
   Self-explanatory.

*Surf file does not contain lines*
   Required for a 2d simulation.

*Surf file does not contain points*
   Self-explanatory.

*Surf file does not contain triangles*
   Required for a 3d simulation.

*Surf-style variables are not yet implemented*
   Self-explanatory.

*Surf\_collide ID must be alphanumeric or underscore characters*
   Self-explanatory.

*Surf\_collide diffuse rotation invalid for 2d*
   Specified rotation vector must be in z-direction.

*Surf\_collide diffuse variable is invalid style*
   It must be an equal-style variable.

*Surf\_collide diffuse variable name does not exist*
   Self-explanatory.

*Surface check failed with %d duplicate edges*
   One or more edges appeared in more than 2 triangles.

*Surface check failed with %d duplicate points*
   One or more points appeared in more than 2 lines.

*Surface check failed with %d infinitely thin line pairs*
   Two adjacent lines have normals in opposite directions
   indicating the lines overlay each other.

*Surface check failed with %d infinitely thin triangle pairs*
   Two adjacent triangles have normals in opposite directions indicating
   the triangles overlay each other.

*Surface check failed with %d points on lines*
   One or more points are on a line they are not an end point of, which
   indicates an ill-formed surface.

*Surface check failed with %d points on triangles*
   One or more points are on a triangle they are not an end point of,
   which indicates an ill-formed surface.

*Surface check failed with %d unmatched edges*
   One or more edges did not appear in a triangle, or appeared only once
   and edge is not on surface of simulation box.

*Surface check failed with %d unmatched points*
   One or more points did not appear in a line, or appeared only once and
   point is not on surface of simulation box.

*Timestep must be >= 0*
   Reset\_timestep cannot be used to set a negative timestep.

*Too big a timestep*
   Reset\_timestep timestep value must fit in a SPARTA big integer, as
   specified by the -DSPARTA\_SMALL, -DSPARTA\_BIG, or -DSPARTA\_BIGBIG
   options in the low-level Makefile used to build SPARTA.  See
   Section 2.2 of the manual for details.

*Too many surfs in one cell*
   Use the global surfmax command to increase this max allowed number of
   surfs per grid cell.

*Too many timesteps*
   The cummulative timesteps must fit in a SPARTA big integer, as as
   specified by the -DSPARTA\_SMALL, -DSPARTA\_BIG, or -DSPARTA\_BIGBIG
   options in the low-level Makefile used to build SPARTA.  See Section
   2.2 of the manual for details.

*Too much buffered per-proc info for dump*
   Number of dumped values per processor cannot exceed a small integer
   (~2 billion values).

*Too much per-proc info for dump*
   Number of local atoms times number of columns must fit in a 32-bit
   integer for dump.

*Unbalanced quotes in input line*
   No matching end double quote was found following a leading double
   quote.

*Unexpected end of data file*
   SPARTA hit the end of the data file while attempting to read a
   section.  Something is wrong with the format of the data file.

*Unexpected end of grid file*
   Self-explantory.

*Unexpected end of surf file*
   Self-explanatory.

*Units command after simulation box is defined*
   The units command cannot be used after a read\_data, read\_restart, or
   create\_box command.

*Universe/uloop variable count < # of partitions*
   A universe or uloop style variable must specify a number of values >= to the
   number of processor partitions.

*Unknown command: %s*
   The command is not known to SPARTA.  Check the input script.

*Unknown outcome in reaction*
   The specified type of the reaction is not encoded in the reaction
   style.

*VSS parameters do not match current species*
   Species cannot be added after VSS colision file is read.

*Variable ID in variable formula does not exist*
   Self-explanatory.

*Variable evaluation before simulation box is defined*
   Cannot evaluate a compute or fix or atom-based value in a variable
   before the simulation has been setup.

*Variable for dump every is invalid style*
   Only equal-style variables can be used.

*Variable for dump image center is invalid style*
   Must be an equal-style variable.

*Variable for dump image phi is invalid style*
   Must be an equal-style variable.

*Variable for dump image theta is invalid style*
   Must be an equal-style variable.

*Variable for dump image zoom is invalid style*
   Must be an equal-style variable.

*Variable for restart is invalid style*
   It must be an equal-style variable.

*Variable for stats every is invalid style*
   It must be an equal-style variable.

*Variable formula compute array is accessed out-of-range*
   Self-explanatory.

*Variable formula compute vector is accessed out-of-range*
   Self-explanatory.

*Variable formula fix array is accessed out-of-range*
   Self-explanatory.

*Variable formula fix vector is accessed out-of-range*
   Self-explanatory.

*Variable has circular dependency*
   A circular dependency is when variable "a" in used by variable "b" and
   variable "b" is also used by varaible "a".  Circular dependencies with
   longer chains of dependence are also not allowed.

*Variable name between brackets must be alphanumeric or underscore characters*
   Self-explanatory.

*Variable name for compute reduce does not exist*
   Self-explanatory.

*Variable name for dump every does not exist*
   Self-explanatory.

*Variable name for dump image center does not exist*
   Self-explanatory.

*Variable name for dump image phi does not exist*
   Self-explanatory.

*Variable name for dump image theta does not exist*
   Self-explanatory.

*Variable name for dump image zoom does not exist*
   Self-explanatory.

*Variable name for fix ave/grid does not exist*
   Self-explanatory.

*Variable name for fix ave/surf does not exist*
   Self-explanatory.

*Variable name for fix ave/time does not exist*
   Self-explanatory.

*Variable name for restart does not exist*
   Self-explanatory.

*Variable name for stats every does not exist*
   Self-explanatory.

*Variable name must be alphanumeric or underscore characters*
   Self-explanatory.

*Variable stats keyword cannot be used between runs*
   Stats keywords that refer to time (such as cpu, elapsed) do not make
   sense in between runs.

*Vertex contains duplicate edge*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Vertex contains edge that doesn't point to it*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Vertex contains invalid edge*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Vertex has less than 3 edges*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*Vertex pointers to last edge are invalid*
   This is an error when calculating how a 3d grid is cut or split by
   surface elements.  It should not normally occur.  Please report the
   issue to the SPARTA developers.

*World variable count doesn't match # of partitions*
   A world-style variable must specify a number of values equal to the
   number of processor partitions.

*Y cannot be periodic for axi-symmetric*
   Self-explanatory.  See the boundary doc page for more details.

*Z dimension must be periodic for 2d simulation*
   Self-explanatory.



.. _warn:

.. raw:: html

   <span id="warn"></span>

Warnings:
---------------------------------------------------------------



*%d particles were in wrong cells on timestep %ld*
   This is the total number of particles that are incorrectly
   matched to their grid cell.

*Grid cell interior corner points marked as unknown = %d*
   Corner points of grid cells interior to the simulation domain were not
   all marked successfully as inside, outside, or overlapping with
   surface elements.  This should normally not happen, but does
   not affect simulations.

*More than one compute ke/particle*
   This may be inefficient since each such compute stores a vector
   of length equal to the number of particles.

*Restart file used different # of processors*
   The restart file was written out by a SPARTA simulation running on a
   different number of processors.  This means you will likely want to
   re-balance the grid cells and particles across processors.  This can
   be done using the balance or fix balance commands.

*Surface check found %d nearly infinitely thin line pairs*
   Two adjacent lines have normals in nearly opposite directions
   indicating the lines nearly overlay each other.

*Surface check found %d nearly infinitely thin triangle pairs*
   Two adjacent triangles have normals in nearly opposite directions
   indicating the triangles nearly overlay each other.

*Surface check found %d points nearly on lines*
   One or more points are nearly on a line they are not an end point of,
   which indicates an ill-formed surface.

*Surface check found %d points nearly on triangles*
   One or more points are nearly on a triangle they are not an end point
   of, which indicates an ill-formed surface.




.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
