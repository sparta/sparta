# kokkos_cuda = KOKKOS package with CUDA backend, default MPI, nvcc/mpicxx
# compiler with OpenMPI or MPICH

include(${CMAKE_CURRENT_LIST_DIR}/kokkos_common.cmake)
# ################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    kokkos_cuda
    CACHE STRING
          "Descriptive string to describe \"spa_\" executable configuration"
          FORCE)
# ################### END   SPARTA OPTIONS ####################

# ################### BEGIN CMAKE OPTIONS ####################
# TODO: Should CMAKE_CXX_COMPILER be set to nvcc_wrapper
# src/KOKKOS/CMakeLists.txt? set(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "")
set(CMAKE_CXX_COMPILER
    ${CMAKE_SOURCE_DIR}lib/kokkos/bin/nvcc_wrapper
    CACHE STRING "" FORCE)
# ################### END CMAKE OPTIONS ####################

# ################### BEGIN KOKKOS OPTIONS ####################
set(Kokkos_ENABLE_Cuda
    ON
    CACHE STRING "")
set(Kokkos_ARCH_KEPLER35
    ON
    CACHE STRING "")
# ################### END   KOKKOS OPTIONS ####################
