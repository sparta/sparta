# kokkos_omp = KOKKOS package with OpenMP backend, MPI compiler, default MPI

include(${CMAKE_SOURCE_DIR}/cmake/presets/astra_kokkos.cmake)
#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_CMAKE_CONFIG_STRING
    kokkos_omp
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
#################### END   SPARTA OPTIONS ####################

#################### BEGIN KOKKOS OPTIONS ####################
set(KOKKOS_ARCH "" CACHE STRING "" FORCE)
#################### END   KOKKOS OPTIONS ####################
