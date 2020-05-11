################################################################################
# This file sets sparta boolean options listed via `cmake -L`
################################################################################

########## BEGIN SPARTA OPTIONAL DEPENDENCIES ##########
option(PKG_MPI_STUBS
        "Enable or disable sparta mpi stubs package. Default: OFF."
        OFF
      )
option(PKG_FFT
        "Enable or disable sparta fft package. Default: OFF."
        OFF
      )
option(PKG_KOKKOS
        "Enable or disable sparta kokkos package. Default: OFF."
        OFF
      )
########## END   SPARTA OPTIONAL DEPENDENCIES ##########

########## BEGIN SPARTA TPL DEPENDENCIES ##########
option(BUILD_MPI_TPL
        "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
        ON
      )

option(BUILD_KOKKOS_TPL
        "Enable or disable KOKKOS TPL. Default: OFF."
        OFF
      )
########## END   SPARTA TPL DEPENDENCIES ##########
