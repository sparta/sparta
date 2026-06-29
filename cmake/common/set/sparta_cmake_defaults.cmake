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
  )

  # When running the KOKKOS regression tests (SPARTA_KOKKOS_EXACT, run with
  # "-k on -sf kk"), skip the inputs that use features which are not yet
  # KOKKOS-enabled and would error out at run time. These tests still run in
  # the non-KOKKOS configurations.
  if(SPARTA_KOKKOS_EXACT)
    list(APPEND SPARTA_DISABLED_TESTS
        # implicit-surface ablation in 3D errors under KOKKOS (zero collision
        # cell volume); the 2D case runs bit-for-bit and is enabled
        "in.ablation.3d"
        # external field fix not KOKKOS-enabled
        "in.bfield"
        "in.bfield.grid"
    )
  endif()

  list(APPEND __DEFAULT_MPI_RANKS "1")
  list(APPEND __DEFAULT_MPI_RANKS "4")
endif()
