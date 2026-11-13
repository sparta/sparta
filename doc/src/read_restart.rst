.. index:: read\_restart

read\_restart command
=====================

**Syntax:**


.. parsed-literal::

   read_restart file keyword args ...

* file = name of binary restart file to read in
* zero or one keyword/args pair may be listed
* keywords = *gridcut* or *balance*
  
  .. parsed-literal::
  
       gridcut arg = cutoff
         cutoff = acquire ghost cells up to this far away (distance units)
       balance args = same as for :doc:`balance_grid <balance_grid>` command



**Examples:**


.. parsed-literal::

   read_restart save.10000
   read_restart restart.\*
   read_restart flow.\*.%
   read_restart save.10000 gridcut -1.0
   read_restart save.10000 balance rcb cell

**Description:**

Read in a previously saved simulation from a restart file.  This
allows continuation of a previous run on the same or different number
of processors.  Information about what is stored in a restart file is
given below.  Basically this operation will re-create the simulation
box with all its particles, the hierarchical grid used to track
particles, and surface elements embedded in the grid, all with their
attributes at the point in time the information was written to the
restart file by a previous simluation.

Although restart files are saved in binary format to allow exact
regeneration of information, the random numbers used in the continued
run will not be identical to those used if the run had been continued.
Hence the new run will not be identical to the continued original run,
but should be statistically similar.

IMPORTANT NOTE: Because restart files are binary, they may not be
portable to other machines.  SPARTA will print an error message if
it cannot read a restart file for this reason.

If a restarted run is performed on the same number of processors as
the original run, then the assignment of grid cells (and their
particles) to processors will be the same as in the original
simulation.  If the processor count changes, then the assignment will
necessarily be different.  In particular, even if the original
assignment was "clumped", meaning each processor's cells were
geometrically compact, the new assignment will not be clumped; it will
be "dispersed".  See `Section 6.8 <Section_howto.html#howto_8>`_ of the
manual for an explanation of clumped and dispersed grid cell
assignments and their relative performance trade-offs.

Note that the restart file contains the setting for the :doc:`global gridcut <global>` command.  If it is >= 0.0 and the assignment of
grid cells to processors is "dispersed" (as described in the
preceeding paragraph), and there are surface elements defined in the
restart file, an error will be triggered.  This is because the
read\_restart command needs to mark all the grid cells as inside vs
outside the defined surface and cannot do this without ghost cell
information.  As explained on the doc page for the :doc:`global gridcut <global>` command, ghost cells cannot be setup with gridcut
>= 0.0 and "dispersed" grid cells.

The solution is to use one of the two keywords listed above, either
*gridcut* or *balance*\ .  The former allows you to reset the grid
cutoff to -1.0 so that ghost cells can be setup.  Note however that
this means each processor will own a copy of all grid cells (at least
until you change it later), which may be undesirable or even
impossible for large problems if it requires too much memory.  The
other solution is to use the *balance* keyword to trigger a re-balance
of the grid cells to processors as soon as the read\_restart command
reads them in.  The arguments for the *balance* keyword are identical
to those for the :doc:`balance\_grid <balance_grid>` command.  If you
choose a balancing style that results in a "clumped" assignment, then
ghost cells will be setup successfully.

NOTE: Only the *gridcut* or the *balance* keyword can be used, not
both of them.


----------


Similar to how restart files are written (see the
:doc:`write\_restart <write_restart>` and :doc:`restart <restart>`
commands), the restart filename can contain two wild-card characters.
If a "\*" appears in the filename, the directory is searched for all
filenames that match the pattern where "\*" is replaced with a timestep
value.  The file with the largest timestep value is read in.  Thus,
this effectively means, read the latest restart file.  It's useful if
you want your script to continue a run from where it left off.  See
the :doc:`run <run>` command and its "upto" option for how to specify
the run command so it doesn't need to be changed either.

If a "%" character appears in the restart filename, SPARTA expects a
set of multiple files to exist.  The :doc:`restart <restart>` and
:doc:`write\_restart <write_restart>` commands explain how such sets are
created.  Read\_restart will first read a filename where "%" is
replaced by "base".  This file tells SPARTA how many processors
created the set and how many files are in it.  Read\_restart then reads
the additional files.  For example, if the restart file was specified
as save.% when it was written, then read\_restart reads the files
save.base, save.0, save.1, ... save.P-1, where P is the number of
processors that created the restart file.

Note that P could be the total number of processors in the previous
simulation, or some subset of those processors, if the *fileper* or
*nfile* options were used when the restart file was written; see the
:doc:`restart <restart>` and :doc:`write\_restart <write_restart>` commands
for details.  The processors in the current SPARTA simulation share
the work of reading these files; each reads a roughly equal subset of
the files.  The number of processors which created the set can be
different than the number of processors in the current SPARTA
simulation.  This can be a fast mode of input on parallel machines
that support parallel I/O.


----------


A restart file stores only the following information about a
simulation, as specified by the associated commands:

* :doc:`units <units>`
* :doc:`dimension <dimension>`
* :doc:`simulation box size <create_box>` and :doc:`boundary conditions <boundary>`
* :doc:`global settings <global>`
* :doc:`particle species info <species>`
* :doc:`mixtures <mixture>`
* geometry of the hierarchical grid that overlays the simulation domain as :doc:`created <create_grid>` or :doc:`read from a file <read_grid>`
* geometry of all defined :doc:`surface elements <read_surf>`
* :doc:`group definitions <group>` for grid cells and surface elements
* :doc:`custom attributes <custom>` for particles, grid cells, or surface elements
* current simulation time
* current :doc:`timestep size <timestep>`
* current timestep number

No other information is stored in the restart file.  Specifically,
information about these simulation entities and their associated
commands is NOT stored:

* :doc:`random number seed <seed>`
* :doc:`computes <compute>`
* :doc:`fixes <fix>`
* :doc:`collision model <collide>`
* :doc:`chemistry (reaction) model <react>`
* :doc:`surface collision models <surf_collide>`
* :doc:`surface reaction models <surf_react>`
* assignment of surfaces/boundaries to surface models
* :doc:`variables <variable>`
* :doc:`regions <region>`
* output options for :doc:`stats <stats_style>`, :doc:`dump <dump>`, :doc:`restart <restart>` files

This means any information specified in the original input script by
these commands needs to be re-specified in the restart input script,
assuming the continued simulation needs the information.

Also note that many commands can be used after a restart file is read,
to override a setting that was stored in the restart file.  For
example, the :doc:`global <global>` command can be used to reset the
values of its specified keywords. If a global command is used in the
input file before the restart file is read, then it will be overriden
by values in the restart file. The only exception is the \*mem/limit\*
command, since it affects how the restart file is processed.

In particular, take note of the following issues:

The status of time-averaging fixes, such as :doc:`fix ave/time <fix_ave_time>`, :doc:`fix ave/grid <fix_ave_grid>`, :doc:`fix ave/surf <fix_ave_surf>`, does not carry over into the restarted
run.  E.g. if the *ave running* option is used with those commands in
the original script and again specified in the restart script, the
running averaged quantities do not persist into the new run.

The :doc:`surf\_modify <surf_modify>` command must be used in the restart
script to assign surface collision models, specified by the
:doc:`surf\_collide <surf_collide>` command, to all :doc:`global boundaries <boundary>` of type "s", and to any surfaces contained
in the restart file, as read in by the :doc:`read\_surf <read_surf>`
command.

If a collision model is specified in the restart script, and the
:doc:`collide\_modify vremax or remain <collide_modify>` command is used
to enable Vremax and fractional collision count to persist for many
timesteps, no information about these quantities persists from the
original simulation to the restarted simulation.  The initial run in
the restart script will re-initialize these data structures.

As noted above, custom attributes of particles, grid cells, or surface
elements defined in the previous input script and stored in the
restart file, will be re-assigned when the restart file is read.

If an input script command which normally defines a custom attribute
is specified, e.g. :doc:`fix ambipolar <fix_ambipolar>`, then if the
custom data for that attribute already exists, it will be re-used.  If
a corresponding input script command is not used, then the custom data
will be stored in the simulation (with particle in this case), but not
be used, which can be inefficient.  The :doc:`custom remove <custom>`
command can be used after the restart file is read, to delete unneded
custom attributes and their data.


----------


**Restrictions:** none

**Related commands:**

:doc:`read\_grid <read_grid>`, :doc:`read\_surf <read_surf>`,
:doc:`write\_restart <write_restart>`, :doc:`restart <restart>`

**Default:** none


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
