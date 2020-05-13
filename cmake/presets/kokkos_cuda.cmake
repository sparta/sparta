# kokkos_cuda = KOKKOS package with CUDA backend, default MPI, nvcc/mpicxx compiler with OpenMPI or MPICH

include(${CMAKE_SOURCE_DIR}/cmake/presets/kokkos_common.cmake)
#################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_CMAKE_CONFIG_STRING
    kokkos_cuda
    CACHE
    STRING "Descriptive string to describe \"spa_\" executable configuration"
    FORCE
    )
#################### END   SPARTA OPTIONS ####################

#################### BEGIN CMAKE OPTIONS ####################
# TODO: Should CMAKE_CXX_COMPILER be set to nvcc_wrapper src/KOKKOS/CMakeLists.txt?
#set(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "")
set(CMAKE_CXX_COMPILER ${CMAKE_SOURCE_DIR}lib/kokkos/bin/nvcc_wrapper CACHE STRING "" FORCE)
#################### END CMAKE OPTIONS ####################

#################### BEGIN KOKKOS OPTIONS ####################
set(Kokkos_ENABLE_Cuda ON CACHE STRING "")
set(KOKKOS_ARCH "Kepler35" CACHE STRING "")
#################### END   KOKKOS OPTIONS ####################
