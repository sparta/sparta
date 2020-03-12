# SPARTA_ENABLE_MPICH_TPL:BOOL
# SPARTA_ENABLE_OMPI_TPL:BOOL
# SPARTA_ENABLE_MPICH_TPL:BOOL
# SPARTA_ENABLE_OPENMP_TPL:BOOL

# SPARTA_ENABLE_KOKKOS:BOOL
# SPARTA_ENABLE_FFT:BOOL

########## BEGIN SPARTA OPTIONAL DEPENDENCIES ##########
option(SPARTA_ENABLE_MPI_STUBS
        "Enable or disable sparta mpi stubs package. Default: ON."
        ON
      )
option(SPARTA_ENABLE_FFT
        "Enable or disable sparta fft package. Default: OFF."
        OFF
      )
option(SPARTA_ENABLE_KOKKOS
        "Enable or disable sparta kokkos package. Default: OFF."
        OFF
      )
########## END   SPARTA OPTIONAL DEPENDENCIES ##########

########## BEGIN SPARTA TPL DEPENDENCIES ##########
option(SPARTA_ENABLE_MPI_TPL
        "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: OFF."
        OFF
      )
########## END   SPARTA TPL DEPENDENCIES ##########
