# mpi = Linux box/cluster, g++, installed MPI

include(${CMAKE_SOURCE_DIR}/cmake/presets/mpi.cmake)
#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_CMAKE_CONFIG_STRING
    mpi_debug
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
#################### END   SPARTA OPTIONS ####################

#################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_CXX_FLAGS "-O -g" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-O -g" CACHE STRING "" FORCE)
#################### END CMAKE OPTIONS ####################
