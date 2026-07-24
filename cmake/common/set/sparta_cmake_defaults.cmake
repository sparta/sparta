# ##############################################################################
# This file sets common default options that all sparta builds use. These
# options can be overridden at configure time via `cmake -DVAR=VAL` or `cmake -C
# /path/to/preset/presets.cmake`
# ##############################################################################
set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS
    -DSPARTA_GZIP
    CACHE
      STRING
      "Compiler flags used when building object files for the \"spa_\" executable"
)

set(SPARTA_MACHINE
    ""
    CACHE
      STRING
      "Suffix to append to spa binary (WON'T enable any features automatically)"
)

if(SPARTA_ENABLE_TESTING)
  set(SPARTA_ENABLED_TEST_SUITES
      "ablation"
      "adapt"
      "vibrate"
      "surf_collide"
      "surf"
      "surf_react_adsorb"
      "step"
      "spiky"
      "sphere"
      "jagged"
      # FAILING."implicit"
      "free"
      "flowfile"
      "emit"
      "collide"
      "circle"
      "chem"
      "cylinder"
      "axi"
      "ambi"
      "relax_const"
      "relax_variable"
      "thermostat"
      "bfield"
      "adjust_temp"
      "shock_tube"
      "variable_timestep"
      "surf_react_heatflux"
      "chem_rates"
      "custom"
      "explicit2implicit"
      "mfp_mct"
      "torque")

  set(SPARTA_DISABLED_TESTS
      "in.ablation.3d.reactions" # Failing
      "in.axi" # Failing
      "in.collide" # Failing
      "in.ambi" # Failing
      "in.cylinder" # Long runtime
      "in.jagged.3d" # Long runtime
      "in.jagged.3d.distributed" # Long runtime
      "in.custom.cube.read.restart" # Failing
      "in.custom.cube.set.restart" # Failing
      "in.custom.step.read.restart" # Failing
      "in.custom.step.set.restart" # Failing
      # Stage 2 of a write_isurf/read_isurf round trip: requires the corner
      # point file written by in.exp2imp.axi.spherecone.readback, so it
      # cannot be run standalone.
      "in.exp2imp.axi.spherecone.readback2"
  )

  # When running the KOKKOS regression tests (SPARTA_KOKKOS_EXACT, run with
  # "-k on -sf kk"), skip the inputs that use features which are not yet
  # KOKKOS-enabled and would error out at run time. These tests still run in
  # the non-KOKKOS configurations.
  if(SPARTA_KOKKOS_EXACT)
    list(APPEND SPARTA_DISABLED_TESTS
        # fix ave/grid for grid/surf inputs not yet supported in KOKKOS
        "in.ablation.2d"
        "in.ablation.3d"
        "in.ablate.axi.spherecone"
        # surf_collide adiabatic/cll/td/impulsive styles not KOKKOS-enabled
        "in.beam.adiabatic"
        "in.beam.cll"
        "in.beam.impulsive"
        "in.beam.td"
        "in.circle.adiabatic"
        "in.circle.cll"
        "in.circle.impulsive"
        "in.circle.td"
        # surf_react gs/ps styles use a non-KOKKOS-enabled surf_collide method
        "in.beam.face.gs"
        "in.beam.face.gs_ps"
        "in.beam.face.ps"
        "in.beam.surf.gs"
        "in.beam.surf.gs_ps"
        "in.beam.surf.ps"
        "in.circle.gs"
        "in.circle.gs_ps"
        "in.circle.ps"
        # external field fix not KOKKOS-enabled
        "in.bfield"
        "in.bfield.grid"
        # VTK dump styles have no KOKKOS variant
        "in.vtk"
        "in.vtk.3d"
    )
  endif()

  list(APPEND __DEFAULT_MPI_RANKS "1")
  list(APPEND __DEFAULT_MPI_RANKS "4")
endif()
