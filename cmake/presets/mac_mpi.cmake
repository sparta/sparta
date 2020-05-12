include(${CMAKE_SOURCE_DIR}/cmake/presets/mac.cmake)

#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_CMAKE_CONFIG_STRING
    mac_mpi
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
set(BUILD_MPI_TPL
    ON
    CACHE
    BOOL
    "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
    FORCE
    )
set(PKG_MPI_STUBS
    OFF
    CACHE
    BOOL
    "Enable or disable sparta mpi stubs package. Default: OFF."
    FORCE
    )
#################### END SPARTA OPTIONS ####################
