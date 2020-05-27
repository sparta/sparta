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
if(BUILD_MPI AND PKG_MPI_STUBS)
    message(WARNING "Both BUILD_MPI: ${BUILD_MPI} and PKG_MPI_STUBS: "
                    "${PKG_MPI_STUBS} are selected. Defaulting to PKG_MPI_STUBS.")
endif()

if(BUILD_MPI AND NOT PKG_MPI_STUBS)
    find_package(MPI REQUIRED)
    # TODO: if NOT MPI_FOUND, handle finding mpi installs
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
if(FFT AND PKG_FFT)
    message(WARNING "Both FFT: ${FFT} and PKG_FFT: ${PKG_FFT} are selected. "
                    "Defaulting to PKG_FFT.")
endif()

if(PKG_FFT)
    set(FFT OFF CACHE STRING "" FORCE)
    set(TARGET_SPARTA_PKG_FFT pkg_fft)
    list(APPEND TARGET_SPARTA_PKGS ${TARGET_SPARTA_PKG_FFT})
    set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DFFT_NONE ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
endif()

if(FFT AND NOT PKG_FFT)
    find_package(${FFT} REQUIRED)
    set(TARGET_SPARTA_BUILD_FFT ${FFT}::${FFT})
    set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DFFT_${FFT} ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
    list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_FFT})
    
    set(PKG_FFT ON CACHE STRING "" FORCE)
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

#################### BEGIN PROCESS EXTRA OPTIONS ####################

#################### END   PROCESS EXTRA OPTIONS ####################
