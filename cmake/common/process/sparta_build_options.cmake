################################################################################
# This file process boolean options set by 
# cmake/common/set/sparta_build_options.cmake and sets up the correct SPARTA
# targets based on those options.
# This file is also responsible for handeling depencies among boolean options.
################################################################################

#################### BEGIN COMBINE CXX FLAGS ####################
set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS ${SPARTA_CXX_COMPILE_FLAGS} ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
#################### BEGIN COMBINE CXX FLAGS ####################

#################### BEGIN PROCESS MPI TPL/PKG ####################
if(BUILD_MPI AND NOT PKG_MPI_STUBS)
    find_package(MPI REQUIRED)
    set(TARGET_SPARTA_BUILD_MPI MPI::MPI_CXX)
    target_compile_options(${TARGET_SPARTA_BUILD_MPI} INTERFACE ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
else()
    set(PKG_MPI_STUBS ON)
    set(BUILD_MPI OFF)
    set(TARGET_SPARTA_BUILD_MPI pkg_mpi_stubs)
endif()
list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_MPI})
#################### END PROCESS MPI TPL/PKG ####################

#################### BEGIN PROCESS PKGS ####################
if(PKG_FFT)
    set(TARGET_SPARTA_PKG_FFT pkg_fft)
    list(APPEND TARGET_SPARTA_PKGS ${TARGET_SPARTA_PKG_FFT})
endif()

if(PKG_KOKKOS)
    set(TARGET_SPARTA_PKG_KOKKOS pkg_kokkos)
    list(APPEND TARGET_SPARTA_PKGS ${TARGET_SPARTA_PKG_KOKKOS})
    # PKG_KOKKOS depends on BUILD_KOKKOS
    set(BUILD_KOKKOS ON)
endif()
#################### END PROCESS PKGS ####################

#################### BEGIN PROCESS TPLS ####################
if(BUILD_KOKKOS)
    # The 2 lines below find an install version of Kokkos. We are using and
    # in-tree copy of kokkos in sparta for now.
    #find_package(KOKKOS REQUIRED)
    #set(TARGET_SPARTA_BUILD_KOKKOS Kokkos::kokkos)
    set(TARGET_SPARTA_BUILD_KOKKOS kokkos)
    list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_KOKKOS})
    # BUILD_KOKKOS does not depend on PKG_KOKKOS, do not attempt to resolve dependency
endif()
#################### END PROCESS TPLS ####################
