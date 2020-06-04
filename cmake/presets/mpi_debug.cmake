# mpi = Linux box/cluster, g++, installed MPI

include(${CMAKE_CURRENT_LIST_DIR}/mpi.cmake)
#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
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
