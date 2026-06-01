# vortex_kokkos = KOKKOS package with CUDA backend, NVIDIA V100 GPU + IBM Power9
# CPU, OpenMPI, nvcc_wrapper compiler, cuFFT
#
# This preset targets a Power9 + NVIDIA V100 system using OpenMPI.
# MY_MPI_PATH should point to the MPI installation's bin directory.

include(${CMAKE_CURRENT_LIST_DIR}/kokkos_common.cmake)
# ################### BEGIN SPARTA OPTIONS ####################
set(SPARTA_MACHINE
    vortex_kokkos
    CACHE STRING
          "Descriptive string to describe \"spa_\" executable configuration"
          FORCE)
# ################### END   SPARTA OPTIONS ####################

# ################### BEGIN CMAKE OPTIONS ####################
set(CMAKE_CXX_COMPILER
    ${CMAKE_CURRENT_LIST_DIR}/../../lib/kokkos/bin/nvcc_wrapper
    CACHE STRING "" FORCE)
# ################### END CMAKE OPTIONS ####################

# ################### BEGIN KOKKOS OPTIONS ####################
set(Kokkos_ENABLE_CUDA
    ON
    CACHE STRING "")
set(Kokkos_ARCH_POWER9
    ON
    CACHE STRING "")
set(Kokkos_ARCH_VOLTA70
    ON
    CACHE STRING "")
# ################### END   KOKKOS OPTIONS ####################

# If FFT package is also enabled, use CUFFT for FFTs
set(FFT_KOKKOS "CUFFT" CACHE STRING "" FORCE)

# Vortex uses OpenMPI. Uncomment and set MY_MPI_PATH if cmake cannot find MPI:
# set(MPI_CXX_INCLUDE_PATH "$ENV{MY_MPI_PATH}../include" CACHE STRING "" FORCE)
# set(CMAKE_EXE_LINKER_FLAGS
#     "-L$ENV{MY_MPI_PATH}../lib -lmpi_ibm -lopen-pal -lopen-rte -lhwloc_ompi -levent_core -levent_pthreads -levent"
#     CACHE STRING "" FORCE)
