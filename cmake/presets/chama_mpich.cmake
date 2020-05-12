# chama_mpich = mpicxx, MPICH (with MVAPICH module loaded)

include(${CMAKE_SOURCE_DIR}/cmake/presets/chama.cmake)
#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_CMAKE_CONFIG_STRING
    chama_mpich
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
#################### END   SPARTA OPTIONS ####################

#################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
#################### END CMAKE OPTIONS ####################
