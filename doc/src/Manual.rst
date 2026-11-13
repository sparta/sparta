.. raw:: html

   <H1></H1>

SPARTA Documentation
====================

27 Aug 2026 version
-------------------

Version info:
-------------

The SPARTA "version" is the date when it was released, such as 3 Mar
2014. SPARTA is updated continuously.  Whenever we fix a bug or add a
feature, we release it immediately, and post a notice on `this page of the WWW site <bug_>`_.  Each dated copy of SPARTA contains all the
features and bug-fixes up to and including that version date. The
version date is printed to the screen and logfile every time you run
SPARTA. It is also in the file src/version.h and in the SPARTA
directory name created when you unpack a tarball, and at the top of
the first page of the manual (this page).

* If you browse the HTML doc pages on the SPARTA WWW site, they always
  describe the most current version of SPARTA.
* The doc pages are not shipped built.  They are generated from the
  reStructuredText sources in doc/src, so they always describe the
  version you have.  "make -C doc html" builds them into doc/html, and
  "make -C doc pdf" builds the `PDF file <Manual.pdf>`_ of the whole
  manual beside them.  "make -C doc" on its own lists the targets.
* At some point, there also will be a Developer.pdf file
  in the doc directory, which describes the internal structure and
  algorithms of SPARTA.

SPARTA stands for Stochastic PArallel Rarefied-gas Time-accurate
Analyzer.

SPARTA is a Direct Simulation Montel Carlo (DSMC) simulator designed
to run efficiently on parallel computers.  It was developed at Sandia
National Laboratories, a US Department of Energy facility, with
funding from the DOE.  It is an open-source code, distributed freely
under the terms of the GNU Public License (GPL), or sometimes by
request under the terms of the GNU Lesser General Public License
(LGPL).

The primary developers of SPARTA are `Steve Plimpton <sjp_>`_, and Michael
Gallis who can be contacted at sjplimp at gmail.com and magalli at
sandia.gov.  The `SPARTA WWW Site <sws_>`_ at https://sparta.github.io has
more information about the code and its uses.

.. _bug: https://sparta.github.io/bug.html



.. _sjp: https://sjplimp.github.io




----------


The SPARTA documentation is organized into the following sections.  If
you find errors or omissions in this manual or have suggestions for
useful information to add, please send an email to the developers so
we can improve the SPARTA documentation.

Once you are familiar with SPARTA, you may want to bookmark `this page <Section_commands.html#comm>`_ at Section\_commands.html#comm since
it gives quick access to documentation for all SPARTA commands.

`PDF file <Manual.pdf>`_ of the entire manual, built by
"make -C doc pdf"

* :doc:`Introduction <Section_intro>`
  1.1 `What is SPARTA <intro_1_>`_
  1.2 `SPARTA features <intro_2_>`_
  1.3 `Grids and surfaces in SPARTA <intro_3_>`_
  1.4 `Open source distribution <intro_4_>`_
  1.5 `Acknowledgments and citations <intro_5_>`_

* :doc:`Getting started <Section_start>`
  2.1 `What's in the SPARTA distribution <start_1_>`_
  2.2 `Making SPARTA <start_2_>`_
  2.3 `Building SPARTA with optional packages <start_3_>`_
  2.4 `Building SPARTA as a library <start_4_>`_
  2.5 `Testing SPARTA <start_5_>`_
  2.6 `Running SPARTA <start_6_>`_
  2.7 `Command-line options <start_7_>`_
  2.8 `Screen output <start_8_>`_

* :doc:`Commands <Section_commands>`
  3.1 `SPARTA input script <cmd_1_>`_
  3.2 `Parsing rules <cmd_2_>`_
  3.3 `Input script structure <cmd_3_>`_
  3.4 `Commands listed by category <cmd_4_>`_
  3.5 `Commands listed alphabetically <cmd_5_>`_

* :doc:`Packages <Section_packages>`
* :doc:`Accelerating SPARTA performance <Section_accelerate>`
  5.1 `Measuring performance <acc_1_>`_
  5.2 `Packages with optimized styles <acc_2_>`_
  5.3 `KOKKOS package <acc_3_>`_

* :doc:`How-to discussions <Section_howto>`
  6.1 `2d simulations <howto_1_>`_
  6.2 `Axisymmetric simulations <howto_2_>`_
  6.3 `Running multiple simulations from one input script <howto_3_>`_
  6.4 `Output from SPARTA <howto_4_>`_
  6.5 `Visualizing SPARTA snapshots <howto_5_>`_
  6.6 `Library interface to SPARTA <howto_6_>`_
  6.7 `Coupling SPARTA to other codes <howto_7_>`_
  6.8 `Details of grid geometry in SPARTA <howto_8_>`_
  6.9 `Details of surfaces in SPARTA <howto_9_>`_
  6.10 `Restarting a simulation <howto_10_>`_
  6.11 `Using the ambipolar approximation <howto_11_>`_
  6.12 `Using multiple vibrational energy levels <howto_12_>`_
  6.13 `Surface elements: explicit, implicit, distributed <howto_13_>`_
  6.14 `Implicit surface ablation <howto_14_>`_
  6.15 `Transparent surface elements <howto_15_>`_
  6.16 `Visualizing SPARTA output with ParaView <howto_16_>`_
  6.17 `Custom per-particle, per-grid, per-surf attributes <howto_17_>`_
  6.18 `Variable timestep simulations <howto_18_>`_
  6.19 `Details of particles in SPARTA <howto_19_>`_

* :doc:`Example problems <Section_example>`
* :doc:`Performance & scalability <Section_perf>`
* :doc:`Additional tools <Section_tools>`
* :doc:`Modifying & extending SPARTA <Section_modify>`
  10.1 `Compute styles <mod_1_>`_
  10.2 `Fix styles <mod_2_>`_
  10.3 `Region styles <mod_3_>`_
  10.4 `Collision styles <mod_4_>`_
  10.5 `Surface collision styles <mod_5_>`_
  10.6 `Chemistry styles <mod_6_>`_
  10.7 `Dump styles <mod_7_>`_
  10.8 `Input script commands <mod_8_>`_

* :doc:`Python interface <Section_python>`
  11.1 `Creating a shared MPI library <py_1_>`_
  11.2 `Extending Python with a parallel version of SPARTA <py_2_>`_
  11.3 `Extending Python with MPI <py_3_>`_
  11.4 `Testing the Python-SPARTA interface <py_4_>`_
  11.5 `Using SPARTA from Python <py_5_>`_
  11.6 `Example Python scripts that use SPARTA <py_6_>`_

* :doc:`Errors <Section_errors>`
  12.1 `Common problems <err_1_>`_
  12.2 `Reporting bugs <err_2_>`_
  12.3 `Error & warning messages <err_3_>`_

* :doc:`Future and history <Section_history>`
  13.1 `Coming attractions <hist_1_>`_
  13.2 `Past versions <hist_2_>`_



.. _intro_1: Section_intro.html#intro_1



.. _intro_2: Section_intro.html#intro_2



.. _intro_3: Section_intro.html#intro_3



.. _intro_4: Section_intro.html#intro_4



.. _intro_5: Section_intro.html#intro_5



.. _start_1: Section_start.html#start_1



.. _start_2: Section_start.html#start_2



.. _start_3: Section_start.html#start_3



.. _start_4: Section_start.html#start_4



.. _start_5: Section_start.html#start_5



.. _start_6: Section_start.html#start_6



.. _start_7: Section_start.html#start_7



.. _start_8: Section_start.html#start_8



.. _cmd_1: Section_commands.html#cmd_1



.. _cmd_2: Section_commands.html#cmd_2



.. _cmd_3: Section_commands.html#cmd_3



.. _cmd_4: Section_commands.html#cmd_4



.. _cmd_5: Section_commands.html#cmd_5



.. _acc_1: Section_accelerate.html#acc_1



.. _acc_2: Section_accelerate.html#acc_2



.. _acc_3: Section_accelerate.html#acc_3



.. _howto_1: Section_howto.html#howto_1



.. _howto_2: Section_howto.html#howto_2



.. _howto_3: Section_howto.html#howto_3



.. _howto_4: Section_howto.html#howto_4



.. _howto_5: Section_howto.html#howto_5



.. _howto_6: Section_howto.html#howto_6



.. _howto_7: Section_howto.html#howto_7



.. _howto_8: Section_howto.html#howto_8



.. _howto_9: Section_howto.html#howto_9



.. _howto_10: Section_howto.html#howto_10



.. _howto_11: Section_howto.html#howto_11



.. _howto_12: Section_howto.html#howto_12



.. _howto_13: Section_howto.html#howto_13



.. _howto_14: Section_howto.html#howto_14



.. _howto_15: Section_howto.html#howto_15



.. _howto_16: Section_howto.html#howto_16



.. _howto_17: Section_howto.html#howto_17



.. _howto_18: Section_howto.html#howto_18



.. _howto_19: Section_howto.html#howto_19



.. _mod_1: Section_modify.html#mod_1



.. _mod_2: Section_modify.html#mod_2



.. _mod_3: Section_modify.html#mod_3



.. _mod_4: Section_modify.html#mod_4



.. _mod_5: Section_modify.html#mod_5



.. _mod_6: Section_modify.html#mod_6



.. _mod_7: Section_modify.html#mod_7



.. _mod_8: Section_modify.html#mod_8



.. _py_1: Section_python.html#py_1



.. _py_2: Section_python.html#py_2



.. _py_3: Section_python.html#py_3



.. _py_4: Section_python.html#py_4



.. _py_5: Section_python.html#py_5



.. _py_6: Section_python.html#py_6



.. _err_1: Section_errors.html#err_1



.. _err_2: Section_errors.html#err_2



.. _err_3: Section_errors.html#err_3



.. _hist_1: Section_history.html#hist_1



.. _hist_2: Section_history.html#hist_2



.. raw:: html

   </BODY>


.. _sws: https://sparta.github.io
.. _sd: Manual.html
.. _sc: Section_commands.html

.. toctree::
   :hidden:
   :maxdepth: 2

   Section_intro
   Section_start
   Section_commands
   Section_packages
   Section_accelerate
   Section_howto
   Section_example
   Section_perf
   Section_tools
   Section_modify
   Section_python
   Section_errors
   Section_history
   adapt_grid
   balance_grid
   bound_modify
   boundary
   clear
   collide
   collide_modify
   compute
   compute_boundary
   compute_count
   compute_distsurf_grid
   compute_dt_grid
   compute_eflux_grid
   compute_fft_grid
   compute_gas_collision_grid
   compute_gas_collision_tally
   compute_gas_reaction_grid
   compute_gas_reaction_tally
   compute_grid
   compute_isurf_grid
   compute_ke_particle
   compute_lambda_grid
   compute_pflux_grid
   compute_property_grid
   compute_property_surf
   compute_react_boundary
   compute_react_isurf_grid
   compute_react_surf
   compute_reduce
   compute_sonine_grid
   compute_surf
   compute_surf_collision_tally
   compute_surf_reaction_tally
   compute_temp
   compute_thermal_grid
   compute_tvib_grid
   create_box
   create_grid
   create_isurf
   create_particles
   custom
   dimension
   dump
   dump_image
   dump_modify
   dump_vtk
   echo
   fix
   fix_ablate
   fix_adapt
   fix_ambipolar
   fix_ave_grid
   fix_ave_histo
   fix_ave_surf
   fix_ave_time
   fix_balance
   fix_controller
   fix_custom
   fix_dt_reset
   fix_emit_face
   fix_emit_face_file
   fix_emit_surf
   fix_field_grid
   fix_field_particle
   fix_grid_check
   fix_halt
   fix_move_surf
   fix_print
   fix_surf_temp
   fix_temp_global_rescale
   fix_temp_rescale
   fix_vibmode
   global
   group
   if
   include
   jump
   label
   log
   mixture
   move_surf
   next
   package
   partition
   print
   python
   quit
   react
   react_modify
   read_grid
   read_isurf
   read_particles
   read_restart
   read_surf
   region
   remove_surf
   reset_timestep
   restart
   run
   scale_particles
   seed
   shell
   species
   species_modify
   stats
   stats_modify
   stats_style
   suffix
   surf_collide
   surf_modify
   surf_react
   surf_react_adsorb
   timestep
   uncompute
   undump
   unfix
   units
   variable
   write_grid
   write_isurf
   write_restart
   write_surf
