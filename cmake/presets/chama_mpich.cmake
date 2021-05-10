# chama_mpich = mpicxx, MPICH (with MVAPICH module loaded)

include(${CMAKE_CURRENT_LIST_DIR}/chama.cmake)
# ################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    chama_mpich
    CACHE STRING
          "Descriptive string to describe \"spa_\" executable configuration"
          FORCE)
# ################### END   SPARTA OPTIONS ####################

# ################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_CXX_COMPILER
    "mpicxx"
    CACHE STRING "" FORCE)
# ################### END CMAKE OPTIONS ####################
