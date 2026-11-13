3. Commands
===========

This section describes how a SPARTA input script is formatted and what
commands are used to define a SPARTA simulation.

| 3.1 `SPARTA input script <#cmd_1>`_
| 3.2 `Parsing rules <#cmd_2>`_
| 3.3 `Input script structure <#cmd_3>`_
| 3.4 `Commands listed by category <#cmd_4>`_
| 3.5 `Commands listed alphabetically <#cmd_5>`_ 
| 


----------


.. _cmd_1:

.. raw:: html

   <span id="cmd_1"></span>

3.1 SPARTA input script
-----------------------

SPARTA executes by reading commands from a input script (text file),
one line at a time.  When the input script ends, SPARTA exits.  Each
command causes SPARTA to take some action.  It may set an internal
variable, read in a file, or run a simulation.  Most commands have
default settings, which means you only need to use the command if you
wish to change the default.

In many cases, the ordering of commands in an input script is not
important.  However the following rules apply:

\(1) SPARTA does not read your entire input script and then perform a
simulation with all the settings.  Rather, the input script is read
one line at a time and each command takes effect when it is read.
Thus this sequence of commands:


.. parsed-literal::

   timestep 0.5 
   run      100 
   run      100

does something different than this sequence:


.. parsed-literal::

   run      100 
   timestep 0.5 
   run      100

In the first case, the specified timestep (0.5 secs) is used for two
simulations of 100 timesteps each.  In the 2nd case, the default
timestep (1.0 sec is used for the 1st 100 step simulation and a 0.5
fmsec timestep is used for the 2nd one.

\(2) Some commands are only valid when they follow other commands.  For
example you cannot define the grid overlaying the simulation box until
the box itself has been defined.  Likewise you cannot read in
triangulated surfaces until a grid has been defined to store them.

Many input script errors are detected by SPARTA and an ERROR or
WARNING message is printed.  :doc:`Section 12 <Section_errors>` gives
more information on what errors mean.  The documentation for each
command lists restrictions on how the command can be used.


----------


.. _cmd_2:

.. raw:: html

   <span id="cmd_2"></span>

3.2 Parsing rules
-----------------

Each non-blank line in the input script is treated as a command.
SPARTA commands are case sensitive.  Command names are lower-case, as
are specified command arguments.  Upper case letters may be used in
file names or user-chosen ID strings.

Here is how each line in the input script is parsed by SPARTA:

\(1) If the last printable character on the line is a "&" character
(with no surrounding quotes), the command is assumed to continue on
the next line.  The next line is concatenated to the previous line by
removing the "&" character and newline.  This allows long commands to
be continued across two or more lines.

\(2) All characters from the first "#" character onward are treated as
comment and discarded.  See an exception in (6).  Note that a
comment after a trailing "&" character will prevent the command from
continuing on the next line.  Also note that for multi-line commands a
single leading "#" will comment out the entire command.

\(3) The line is searched repeatedly for $ characters, which indicate
variables that are replaced with a text string.  See an exception in
(6).

If the $ is followed by curly brackets, then the variable name is the
text inside the curly brackets.  If no curly brackets follow the $,
then the variable name is the single character immediately following
the $.  Thus ${myTemp} and $x refer to variable names "myTemp" and
"x".

How the variable is converted to a text string depends on what style
of variable it is; see the :doc:`variable <variable>` doc page for details.
It can be a variable that stores multiple text strings, and return one
of them.  The returned text string can be multiple "words" (space
separated) which will then be interpreted as multiple arguments in the
input command.  The variable can also store a numeric formula which
will be evaluated and its numeric result returned as a string.

As a special case, if the $ is followed by parenthesis, then the text
inside the parenthesis is treated as an "immediate" variable and
evaluated as an :doc:`equal-style variable <variable>`.  This is a way
to use numeric formulas in an input script without having to assign
them to variable names.  For example, these 3 input script lines:


.. parsed-literal::

   variable X equal (xlo+xhi)/2+sqrt(v_area)
   region 1 block $X 2 INF INF EDGE EDGE
   variable X delete

can be replaced by


.. parsed-literal::

   region 1 block $((xlo+xhi)/2+sqrt(v_area)) 2 INF INF EDGE EDGE

so that you do not have to define (or discard) a temporary variable X.

Note that neither the curly-bracket or immediate form of variables can
contain nested $ characters for other variables to substitute for.
Thus you cannot do this:


.. parsed-literal::

   variable        a equal 2
   variable        b2 equal 4
   print           "B2 = ${b$a}"

Nor can you specify this $($x-1.0) for an immediate variable, but
you could use $(v\_x-1.0), since the latter is valid syntax for an
:doc:`equal-style variable <variable>`.

See the :doc:`variable <variable>` command for more details of how
strings are assigned to variables and evaluated, and how they can be
used in input script commands.

\(4) The line is broken into "words" separated by whitespace (tabs,
spaces).  Note that words can thus contain letters, digits,
underscores, or punctuation characters.

\(5) The first word is the command name.  All successive words in the
line are arguments.

\(6) If you want text with spaces to be treated as a single argument,
it can be enclosed in either single (') or double (") or triple quotes
(""").  A long single argument enclosed in single or double quotes can
span multiple lines if the "&" character is used, as described above.
When the lines are concatenated together by SPARTA (and the "&"
characters and line breaks removed), the combined text will become a
single line.  If you want multiple lines of an argument to retain
their line breaks, the text can be enclosed in triple quotes, in which
case "&" characters are not needed and do not function as line
continuation character.

For example:

print "Volume = $v"
print 'Volume = $v'
print """
System volume = $v
System temperature = $t
"""
variable a string "red green blue &
purple orange cyan"
if "$\ *steps* > 1000" then quit

In each of these cases, the single, double, or triple quotes are
removed and the enclosed text stored internally as a single argument.

See the :doc:`dump modify format <dump_modify>`, :doc:`print <print>`,
:doc:`if <if>`, or :doc:`python <python>` commands for examples.

A "#" or "$" character that is between quotes will not be treated as a
comment indicator in (2) or substituted for as a variable in (3).

IMPORTANT NOTE: If the argument is itself a command that requires a
quoted argument (e.g. using a :doc:`print <print>` command as part of an
:doc:`if <if>` or :doc:`run every <run>` command), then single, double, or
triple quotes can be nested in the usual manner.  See the doc pages
for those commands for examples.  Only one of level of nesting is
allowed, but that should be sufficient for most use cases.


----------


.. _cmd_3:

.. raw:: html

   <span id="cmd_3"></span>

3.3 Input script structure
----------------------------------------------------------------------------------

This section describes the structure of a typical SPARTA input script.
The "examples" directory in the SPARTA distribution contains sample
input scripts; the corresponding problems are discussed in :doc:`Section 5 <Section_example>`, and animated on the `SPARTA WWW Site <sws_>`_.

A SPARTA input script typically has 4 parts:

1. Initialization
2. Problem definition
3. Settings
4. Run a simulation

The last 2 parts can be repeated as many times as desired.  I.e. run a
simulation, change some settings, run some more, etc.  Each of the 4
parts is now described in more detail.  Remember that almost all the
commands need only be used if a non-default value is desired.

\(1) Initialization

Set parameters that need to be defined before the simulation domain,
particles, grid cells, and surfaces are defined.

Relevant commands include :doc:`dimension <dimension>`,
:doc:`units <units>`, and :doc:`seed <seed>`.

\(2) Problem definition

These items must be defined before running a SPARTA calculation, and
typically in this order:

* :doc:`create\_box <create_box>` for the simulation box
* :doc:`create\_grid <create_grid>` or :doc:`read\_grid <read_grid>` for grid cells
* :doc:`read\_surf <read_surf>` or :doc:`read\_isurf <read_isurf>` for surfaces
* :doc:`species <species>` for particle species properties
* :doc:`create\_particles <create_particles>` for particles

The first two are required.  Surfaces are optional.  Particles are also
optional in the setup stage, since they can be added as the simulation
runs.

The system can also be load-balanced after the grid and/or particles
are defined in the setup stage using the
:doc:`balance\_grid <balance_grid>` command.  The grid can also be
adapted before or betwee simulations using the
:doc:`adapt\_grid <adapt_grid>` command.

\(3) Settings

Once the problem geometry, grid cells, surfaces, and particles are
defined, a variety of settings can be specified, which include
simulation parameters, output options, etc.

Commands that do this include

:doc:`global <global>`
:doc:`timestep <timestep>`
:doc:`collide <collide>` for a collision model
:doc:`react <react>` for a chemisty model
:doc:`fix <fix>` for boundary conditions, time-averaging, load-balancing, etc
:doc:`compute <compute>` for diagnostic computations
:doc:`stats\_style <stats_style>` for screen output
:doc:`dump <dump>` for snapshots of particle, grid, and surface info
:doc:`dump image <dump>` for on-the-fly images of the simulation
:doc:`dump vtk <dump_vtk>` for native VTK-format snapshots (VTK package)

\(4) Run a simulation

A simulation is run using the :doc:`run <run>` command.


----------


.. _cmd_4:

.. raw:: html

   <span id="cmd_4"></span>

3.4 Commands listed by category
-------------------------------

This section lists many SPARTA commands, grouped by category.  The
`next section <#cmd_5>`_ lists all commands alphabetically.

Initialization:

:doc:`dimension <dimension>`, :doc:`package <package>`, :doc:`seed <seed>`,
:doc:`suffix <suffix>`, :doc:`units <units>`

Problem definition:

:doc:`boundary <boundary>`, :doc:`bound\_modify <bound_modify>`,
:doc:`create\_box <create_box>`, :doc:`create\_grid <create_grid>`,
:doc:`create\_particles <create_particles>`, :doc:`mixture <mixture>`,
:doc:`read\_grid <read_grid>`, :doc:`read\_isurf <read_isurf>`,
:doc:`read\_particles <read_particles>`, :doc:`read\_surf <read_surf>`,
:doc:`read\_restart <read_restart>`, :doc:`species <species>`,

Settings:

:doc:`collide <collide>`, :doc:`collide\_modify <collide_modify>`,
:doc:`compute <compute>`, :doc:`fix <fix>`, :doc:`global <global>`,
:doc:`react <react>`, :doc:`react\_modify <react_modify>`,
:doc:`region <region>`, :doc:`surf\_collide <surf_collide>`,
:doc:`surf\_modify <surf_modify>`, :doc:`surf\_react <surf_react>`,
:doc:`timestep <timestep>`, :doc:`uncompute <uncompute>`,
:doc:`unfix <unfix>`

Output:

:doc:`dump <dump>`, :doc:`dump\_image <dump_image>`,
:doc:`dump\_modify <dump_modify>`, :doc:`restart <restart>`,
:doc:`stats <stats>`, :doc:`stats\_modify <stats_modify>`,
:doc:`stats\_style <stats_style>`, :doc:`undump <undump>`,
:doc:`write\_grid <write_grid>`, :doc:`write\_isurf <write_isurf>`,
:doc:`write\_surf <write_surf>`, :doc:`write\_restart <write_restart>`

Actions:

:doc:`adapt\_grid <adapt_grid>`, :doc:`balance\_grid <balance_grid>`,
:doc:`run <run>`, :doc:`scale\_particles <scale_particles>`

Miscellaneous:

:doc:`clear <clear>`, :doc:`echo <echo>`, :doc:`if <if>`,
:doc:`include <include>`, :doc:`jump <jump>`, :doc:`label <label>`,
:doc:`log <log>`, :doc:`next <next>`, :doc:`partition <partition>`,
:doc:`print <print>`, :doc:`python <python>`, :doc:`quit <quit>`,
:doc:`shell <shell>`, :doc:`variable <variable>`


----------


.. _cmd_5:

.. raw:: html

   <span id="cmd_5"></span>

.. _comm:

.. raw:: html

   <span id="comm"></span>

3.5 Individual commands
-------------------------------------------------------------------------------------------------------------------------------------

This section lists all SPARTA commands alphabetically, with a separate
listing below of styles within certain commands.  The `previous section <#cmd_4>`_ lists many of the same commands, grouped by category.

+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`adapt\_grid <adapt_grid>`           | :doc:`balance\_grid <balance_grid>` | :doc:`boundary <boundary>`            | :doc:`bound\_modify <bound_modify>`     | :doc:`clear <clear>`                    | :doc:`collide <collide>`                        |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`collide\_modify <collide_modify>`   | :doc:`compute <compute>`            | :doc:`create\_box <create_box>`       | :doc:`create\_grid <create_grid>`       | :doc:`create\_isurf <create_isurf>`     | :doc:`create\_particles (k) <create_particles>` |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`custom <custom>`                    | :doc:`dimension <dimension>`        | :doc:`dump <dump>`                    | :doc:`dump image <dump_image>`          | :doc:`dump\_modify <dump_modify>`       | :doc:`dump movie <dump_image>`                  |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`dump vtk <dump_vtk>`                | :doc:`echo <echo>`                  | :doc:`fix <fix>`                      | :doc:`global <global>`                  | :doc:`group <group>`                    | :doc:`if <if>`                                  |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`include <include>`                  | :doc:`jump <jump>`                  | :doc:`label <label>`                  | :doc:`log <log>`                        | :doc:`mixture <mixture>`                | :doc:`move\_surf <move_surf>`                   |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`next <next>`                        | :doc:`package <package>`            | :doc:`partition <partition>`          | :doc:`print <print>`                    | :doc:`python <python>`                  | :doc:`quit <quit>`                              |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`react <react>`                      | :doc:`react\_modify <react_modify>` | :doc:`read\_grid <read_grid>`         | :doc:`read\_isurf <read_isurf>`         | :doc:`read\_particles <read_particles>` | :doc:`read\_restart <read_restart>`             |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`read\_surf (k) <read_surf>`         | :doc:`region <region>`              | :doc:`remove\_surf (k) <remove_surf>` | :doc:`reset\_timestep <reset_timestep>` | :doc:`restart <restart>`                | :doc:`run <run>`                                |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`scale\_particles <scale_particles>` | :doc:`seed <seed>`                  | :doc:`shell <shell>`                  | :doc:`species <species>`                | :doc:`species\_modify <species_modify>` | :doc:`stats <stats>`                            |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`stats\_modify <stats_modify>`       | :doc:`stats\_style <stats_style>`   | :doc:`suffix <suffix>`                | :doc:`surf\_collide <surf_collide>`     | :doc:`surf\_react <surf_react>`         | :doc:`surf\_modify <surf_modify>`               |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`timestep <timestep>`                | :doc:`uncompute <uncompute>`        | :doc:`undump <undump>`                | :doc:`unfix <unfix>`                    | :doc:`units <units>`                    | :doc:`variable <variable>`                      |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+
| :doc:`write\_grid <write_grid>`           | :doc:`write\_isurf <write_isurf>`   | :doc:`write\_restart <write_restart>` | :doc:`write\_surf <write_surf>`         |                                         |                                                 |
+-------------------------------------------+-------------------------------------+---------------------------------------+-----------------------------------------+-----------------------------------------+-------------------------------------------------+


----------


Fix styles
----------

See the :doc:`fix <fix>` command for one-line descriptions of each style
or click on the style itself for a full description.  Some of the
styles have accelerated versions, which can be used if SPARTA is built
with the :doc:`appropriate accelerated package <Section_accelerate>`.
This is indicated by additional letters in parenthesis: k = KOKKOS.

+--------------------------------------+------------------------------------------------+--------------------------------------+----------------------------------------+----------------------------------------------------------+---------------------------------------------+
| :doc:`ablate <fix_ablate>`           | :doc:`adapt (k) <fix_adapt>`                   | :doc:`ambipolar (k) <fix_ambipolar>` | :doc:`ave/grid (k) <fix_ave_grid>`     | :doc:`ave/histo (k) <fix_ave_histo>`                     | :doc:`ave/histo/weight (k) <fix_ave_histo>` |
+--------------------------------------+------------------------------------------------+--------------------------------------+----------------------------------------+----------------------------------------------------------+---------------------------------------------+
| :doc:`ave/surf <fix_ave_surf>`       | :doc:`ave/time <fix_ave_time>`                 | :doc:`balance (k) <fix_balance>`     | :doc:`controller <fix_controller>`     | :doc:`custom <fix_custom>`                               | :doc:`dt/reset (k) <fix_dt_reset>`          |
+--------------------------------------+------------------------------------------------+--------------------------------------+----------------------------------------+----------------------------------------------------------+---------------------------------------------+
| :doc:`emit/face (k) <fix_emit_face>` | :doc:`emit/face/file (k) <fix_emit_face_file>` | :doc:`emit/surf (k) <fix_emit_surf>` | :doc:`field/grid (k) <fix_field_grid>` | :doc:`field/particle (k) <fix_field_particle>`           | :doc:`grid/check (k) <fix_grid_check>`      |
+--------------------------------------+------------------------------------------------+--------------------------------------+----------------------------------------+----------------------------------------------------------+---------------------------------------------+
| :doc:`halt <fix_halt>`               | :doc:`move/surf (k) <fix_move_surf>`           | :doc:`print <fix_print>`             | :doc:`surf/temp (k) <fix_surf_temp>`   | :doc:`temp/global/rescale (k) <fix_temp_global_rescale>` | :doc:`temp/rescale (k) <fix_temp_rescale>`  |
+--------------------------------------+------------------------------------------------+--------------------------------------+----------------------------------------+----------------------------------------------------------+---------------------------------------------+
| :doc:`vibmode (k) <fix_vibmode>`     |                                                |                                      |                                        |                                                          |                                             |
+--------------------------------------+------------------------------------------------+--------------------------------------+----------------------------------------+----------------------------------------------------------+---------------------------------------------+


----------


Compute styles
--------------

See the :doc:`compute <compute>` command for one-line descriptions of
each style or click on the style itself for a full description.  Some
of the styles have accelerated versions, which can be used if SPARTA
is built with the :doc:`appropriate accelerated package <Section_accelerate>`.  This is indicated by additional
letters in parenthesis: k = KOKKOS.

+--------------------------------------------------------------+--------------------------------------------------------------+----------------------------------------------------------+------------------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------+
| :doc:`boundary (k) <compute_boundary>`                       | :doc:`count (k) <compute_count>`                             | :doc:`distsurf/grid (k) <compute_distsurf_grid>`         | :doc:`dt/grid (k) <compute_dt_grid>`                       | :doc:`eflux/grid (k) <compute_eflux_grid>`       | :doc:`fft/grid (k) <compute_fft_grid>`                         |
+--------------------------------------------------------------+--------------------------------------------------------------+----------------------------------------------------------+------------------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------+
| :doc:`gas/collision/grid (k) <compute_gas_collision_grid>`   | :doc:`gas/collision/tally (k) <compute_gas_collision_tally>` | :doc:`gas/reaction/grid (k) <compute_gas_reaction_grid>` | :doc:`gas/reaction/tally (k) <compute_gas_reaction_tally>` | :doc:`grid (k) <compute_grid>`                   | :doc:`isurf/grid (k) <compute_isurf_grid>`                     |
+--------------------------------------------------------------+--------------------------------------------------------------+----------------------------------------------------------+------------------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------+
| :doc:`ke/particle (k) <compute_ke_particle>`                 | :doc:`lambda/grid (k) <compute_lambda_grid>`                 | :doc:`pflux/grid (k) <compute_pflux_grid>`               | :doc:`property/grid (k) <compute_property_grid>`           | :doc:`property/surf (k) <compute_property_surf>` | :doc:`react/boundary (k) <compute_react_boundary>`             |
+--------------------------------------------------------------+--------------------------------------------------------------+----------------------------------------------------------+------------------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------+
| :doc:`react/surf (k) <compute_react_surf>`                   | :doc:`react/isurf/grid (k) <compute_react_isurf_grid>`       | :doc:`reduce (k) <compute_reduce>`                       | :doc:`sonine/grid (k) <compute_sonine_grid>`               | :doc:`surf (k) <compute_surf>`                   | :doc:`surf/collision/tally (k) <compute_surf_collision_tally>` |
+--------------------------------------------------------------+--------------------------------------------------------------+----------------------------------------------------------+------------------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------+
| :doc:`surf/reaction/tally (k) <compute_surf_reaction_tally>` | :doc:`temp (k) <compute_temp>`                               | :doc:`thermal/grid (k) <compute_thermal_grid>`           | :doc:`tvib/grid (k) <compute_tvib_grid>`                   |                                                  |                                                                |
+--------------------------------------------------------------+--------------------------------------------------------------+----------------------------------------------------------+------------------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------+


----------


Collide styles
--------------

See the :doc:`collide <collide>` command for details of each style.
Some of the styles have accelerated versions, which can be used if
SPARTA is built with the :doc:`appropriate accelerated package <Section_accelerate>`.  This is indicated by additional
letters in parenthesis: k = KOKKOS.

+--------------------------+
| :doc:`vss (k) <collide>` |
+--------------------------+


----------


Surface collide styles
----------------------

See the :doc:`surf\_collide <surf_collide>` command for details of each
style.  Some of the styles have accelerated versions, which can be
used if SPARTA is built with the :doc:`appropriate accelerated package <Section_accelerate>`.  This is indicated by additional
letters in parenthesis: k = KOKKOS.

+-------------------------------------+---------------------------------------+------------------------------------+
| :doc:`adiabatic (k) <surf_collide>` | :doc:`cll (k) <surf_collide>`         | :doc:`diffuse (k) <surf_collide>`  |
+-------------------------------------+---------------------------------------+------------------------------------+
| :doc:`impulsive (k) <surf_collide>` | :doc:`piston (k) <surf_collide>`      | :doc:`specular (k) <surf_collide>` |
+-------------------------------------+---------------------------------------+------------------------------------+
| :doc:`td (k) <surf_collide>`        | :doc:`transparent (k) <surf_collide>` | :doc:`vanish (k) <surf_collide>`   |
+-------------------------------------+---------------------------------------+------------------------------------+


----------


Surface reaction styles
-----------------------

See the :doc:`surf\_react <surf_react>` command for details of each
style. Some of the styles have accelerated versions, which can be
used if SPARTA is built with the :doc:`appropriate accelerated package <Section_accelerate>`.  This is indicated by additional
letters in parenthesis: k = KOKKOS.

+---------------------------------------+--------------------------------+
| :doc:`adsorb (k) <surf_react_adsorb>` | :doc:`global (k) <surf_react>` |
+---------------------------------------+--------------------------------+
| :doc:`prob (k) <surf_react>`          |                                |
+---------------------------------------+--------------------------------+


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html
