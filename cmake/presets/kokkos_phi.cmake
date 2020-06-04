# kokkos_phi = KOKKOS package with OpenMP backend, KNL Xeon Phi support, Cray
# System, default MPI

include(${CMAKE_CURRENT_LIST_DIR}/kokkos_common.cmake)
# ################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    kokkos_phi
    CACHE STRING
          "Descriptive string to describe \"spa_\" executable configuration"
          FORCE)
# ################### END   SPARTA OPTIONS ####################

# ################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_C_COMPILER
    "cc"
    CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER
    "CC"
    CACHE STRING "" FORCE)
# ################### END CMAKE OPTIONS ####################

# ################### BEGIN KOKKOS OPTIONS ####################
set(Kokkos_ARCH_KNL
    ON
    CACHE STRING "")
# ################### END   KOKKOS OPTIONS ####################
