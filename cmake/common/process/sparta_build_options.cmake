# ##############################################################################
# This file process boolean options set by
# cmake/common/set/sparta_build_options.cmake and sets up the correct SPARTA
# targets based on those options. This file is also responsible for handeling
# dependencies among boolean options.
# ##############################################################################

# ################### BEGIN PROCESS EXTRA OPTIONS ####################
if(SPARTA_ENABLE_ALL_PKGS AND SPARTA_DISABLE_ALL_PKGS)
  message(FATAL_ERROR "SPARTA_ENABLE_ALL_PKGS=SPARTA_DISABLE_ALL_PKGS=ON")
endif()

if(SPARTA_ENABLE_ALL_PKGS)
  foreach(pkg IN LISTS SPARTA_PKG_LIST)
    get_property(
      pkg_help
      CACHE ${pkg}
      PROPERTY HELPSTRING)
    set(${pkg}
        ON
        CACHE STRING "${pkg_help}" FORCE)
  endforeach()
endif()

if(SPARTA_DISABLE_ALL_PKGS)
  foreach(pkg IN LISTS SPARTA_PKG_LIST)
    get_property(
      pkg_help
      CACHE ${pkg}
      PROPERTY HELPSTRING)
    set(${pkg}
        OFF
        CACHE STRING "${pkg_help}" FORCE)
  endforeach()
endif()
# ################### END   PROCESS EXTRA OPTIONS ####################

# ################### BEGIN PROCESS MPI TPL/PKG ####################
if(BUILD_MPI AND PKG_MPI_STUBS)
  message(WARNING "Both BUILD_MPI: ${BUILD_MPI} and PKG_MPI_STUBS: "
                  "${PKG_MPI_STUBS} are selected. Defaulting to PKG_MPI_STUBS.")
endif()

if(BUILD_MPI AND NOT PKG_MPI_STUBS)
  find_package(MPI REQUIRED)
  # TODO: if NOT MPI_FOUND, handle finding mpi installs
  set(TARGET_SPARTA_BUILD_MPI MPI::MPI_CXX)
else()
  set(PKG_MPI_STUBS ON)
  set(BUILD_MPI OFF)
  set(TARGET_SPARTA_BUILD_MPI pkg_mpi_stubs)
endif()
list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_MPI})
# ################### END PROCESS MPI TPL/PKG ####################

# ################### BEGIN PROCESS PKGS ####################
if(FFT AND PKG_FFT)
  message(WARNING "Both FFT: ${FFT} and PKG_FFT: ${PKG_FFT} are selected. "
                  "Defaulting to PKG_FFT.")
endif()

if(PKG_FFT)
  set(FFT
      OFF
      CACHE STRING "" FORCE)
  set(TARGET_SPARTA_PKG_FFT pkg_fft)
  list(APPEND TARGET_SPARTA_PKGS ${TARGET_SPARTA_PKG_FFT})
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DFFT_NONE
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
endif()

if(FFT AND NOT PKG_FFT)
  find_package(${FFT} REQUIRED)
  set(TARGET_SPARTA_BUILD_FFT ${FFT}::${FFT})
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DFFT_${FFT}
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
  list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_FFT})

  set(PKG_FFT
      ON
      CACHE STRING "" FORCE)
  set(TARGET_SPARTA_PKG_FFT pkg_fft)
  list(APPEND TARGET_SPARTA_PKGS ${TARGET_SPARTA_PKG_FFT})
  if(SPARTA_ENABLE_TESTING)
    set(SPARTA_ENABLED_TEST_SUITES ${SPARTA_ENABLED_TEST_SUITES} "fft")
  endif()
endif()

if(PKG_KOKKOS)
  set(TARGET_SPARTA_PKG_KOKKOS pkg_kokkos)
  list(APPEND TARGET_SPARTA_PKGS ${TARGET_SPARTA_PKG_KOKKOS})
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DSPARTA_KOKKOS
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
  # PKG_KOKKOS depends on BUILD_KOKKOS
  set(BUILD_KOKKOS ON)
endif()
# ################### END PROCESS PKGS ####################

# ################### BEGIN PROCESS TPLS ####################
if(BUILD_KOKKOS)
  # The 2 lines below find an install version of Kokkos. We are using and
  # in-tree copy of kokkos in sparta for now. find_package(KOKKOS REQUIRED)
  # set(TARGET_SPARTA_BUILD_KOKKOS Kokkos::kokkos)
  set(TARGET_SPARTA_BUILD_KOKKOS kokkos)
  list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_KOKKOS})
  # BUILD_KOKKOS does not depend on PKG_KOKKOS, do not attempt to resolve
  # dependency
endif()

if(BUILD_JPEG)
  find_package(JPEG REQUIRED)
  set(TARGET_SPARTA_BUILD_JPEG JPEG::JPEG)
  list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_JPEG})
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DSPARTA_JPEG
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
endif()

if(BUILD_PNG)
  find_package(PNG REQUIRED)
  set(TARGET_SPARTA_BUILD_PNG PNG::PNG)
  list(APPEND TARGET_SPARTA_BUILD_TPLS ${TARGET_SPARTA_BUILD_PNG})
  set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS -DSPARTA_PNG
                                       ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
endif()
# ################### END PROCESS TPLS ####################

# ################### BEGIN COMBINE CXX FLAGS ####################
set(SPARTA_DEFAULT_CXX_COMPILE_FLAGS ${SPARTA_CXX_COMPILE_FLAGS}
                                     ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
# ################### BEGIN COMBINE CXX FLAGS ####################

if(BUILD_MPI)
  set(CRAYPE_VERSION $ENV{CRAYPE_VERSION})
  if(NOT CRAYPE_VERSION)
    set_property(
      TARGET ${TARGET_SPARTA_BUILD_MPI}
      PROPERTY INTERFACE_COMPILE_OPTIONS ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})
  endif()
endif()

if(SPARTA_CTEST_CONFIGS)
  string(REPLACE " " ";" SPARTA_CTEST_CONFIGS "${SPARTA_CTEST_CONFIGS}")
endif()
