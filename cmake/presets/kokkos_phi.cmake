# kokkos_phi = KOKKOS package with OpenMP backend, KNL Xeon Phi support, Cray System, default MPI

include(${CMAKE_SOURCE_DIR}/cmake/presets/kokkos_common.cmake)
#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_CMAKE_CONFIG_STRING
    kokkos_phi
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
#################### END   SPARTA OPTIONS ####################

#################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_C_COMPILER "cc" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "CC" CACHE STRING "" FORCE)
#################### END CMAKE OPTIONS ####################

#################### BEGIN KOKKOS OPTIONS ####################
set(KOKKOS_ARCH "KNL" CACHE STRING "")
#################### END   KOKKOS OPTIONS ####################
