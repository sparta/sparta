# kokkos_omp = KOKKOS package with OpenMP backend, MPI compiler, default MPI

include(${CMAKE_CURRENT_LIST_DIR}/kokkos_common.cmake)
# ################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    kokkos_omp
    CACHE STRING
          "Descriptive string to describe \"spa_\" executable configuration"
          FORCE)
# ################### END   SPARTA OPTIONS ####################

# ################### BEGIN KOKKOS OPTIONS ####################
# ################### END   KOKKOS OPTIONS ####################
