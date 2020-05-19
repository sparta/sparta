# mac_mpi = Apple PowerBook laptop, mpic++

include(${CMAKE_SOURCE_DIR}/cmake/presets/mac.cmake)

#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_CMAKE_CONFIG_STRING
    mac_mpi
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
#################### END SPARTA OPTIONS ####################

#################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_CXX_COMPILER "mpic++" CACHE STRING "" FORCE)
#################### END CMAKE OPTIONS ####################

#################### BEGIN TPL OPTIONS ####################
set(BUILD_MPI
    ON
    CACHE
    BOOL
    "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
    FORCE
    )
#################### END   TPL OPTIONS ####################

#################### BEGIN PACKAGE OPTIONS ####################
set(PKG_MPI_STUBS
    OFF
    CACHE
    BOOL
    "Enable or disable sparta mpi stubs package. Default: OFF."
    FORCE
    )
#################### END   PACKAGE OPTIONS ####################
