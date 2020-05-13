# astra_kokkos = Arm ThunderX2 architecture, KOKKOS package with OpenMP backend, MPI compiler, default MPI

include(${CMAKE_SOURCE_DIR}/cmake/presets/kokkos_common.cmake)
#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_CMAKE_CONFIG_STRING
    astra_kokkos
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
#################### END   SPARTA OPTIONS ####################

#################### BEGIN KOKKOS OPTIONS ####################
set(KOKKOS_ARCH "ARMv8-TX2" CACHE STRING "")
#################### END   KOKKOS OPTIONS ####################
