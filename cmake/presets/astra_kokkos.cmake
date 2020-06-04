# astra_kokkos = Arm ThunderX2 architecture, KOKKOS package with OpenMP backend, MPI compiler, default MPI

include(${CMAKE_CURRENT_LIST_DIR}/kokkos_common.cmake)
#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    astra_kokkos
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
#################### END   SPARTA OPTIONS ####################

#################### BEGIN KOKKOS OPTIONS ####################
set(Kokkos_ARCH_ARMV8_THUNDERX2 ON CACHE STRING "")
#################### END   KOKKOS OPTIONS ####################
