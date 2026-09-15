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
      "tally_computes"
      "ambi_3body"
      "regions"
      "qk"
      "explicit2implicit"
      "mfp_mct"
      "optmove"
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
      # The two stages of the write_isurf/read_isurf round trip. Stage 1
      # writes a corner point file into the shared suite run directory under
      # a fixed name, which would race with itself when ctest runs the mpi_1
      # and mpi_4 copies of the test concurrently, and stage 2 consumes that
      # file so it cannot run standalone. They are kept as a documented
      # example pair; create_isurf itself is covered by the two
      # in.exp2imp.axi.* tests above.
      "in.exp2imp.axi.spherecone.readback"
      "in.exp2imp.axi.spherecone.readback2"
      # Ablation driven by particle flux over 500 steps. The surface geometry
      # is deterministic, but the incident flux that drives it is stochastic
      # and amplifies last-bit floating point differences between platforms,
      # so a gold standard log for it is not portable (compare in.axi above).
      # Kept as a runnable example.
      "in.ablate.axi.spherecone"
  )

  # When running the KOKKOS regression tests (SPARTA_KOKKOS_EXACT, run with
  # "-k on -sf kk"), skip the inputs that use features which are not yet
  # KOKKOS-enabled and would error out at run time. These tests still run in
  # the non-KOKKOS configurations.
  if(SPARTA_KOKKOS_EXACT)
    list(APPEND SPARTA_DISABLED_TESTS
        # VTK dump styles have no KOKKOS variant
        "in.vtk"
        "in.vtk.3d"
    )
  endif()

  list(APPEND __DEFAULT_MPI_RANKS "1")
  list(APPEND __DEFAULT_MPI_RANKS "4")
endif()
