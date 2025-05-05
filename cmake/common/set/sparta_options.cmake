# ##############################################################################
# This file sets sparta boolean options listed via `cmake -L`
# ##############################################################################

# -------------------------------------------------------------------------------
# macro(sparta_option) Adds a sparta package option by invoking cmake's option
# routine and appends the option to OPTIONS_LIST
#
# Inputs: OPTION_NAME OPTION_STRING OPTION_VALUE OPTION_LIST
#
# Outputs: OPTION_LIST OPTION_NAME
macro(sparta_option OPTION_NAME OPTION_STRING OPTION_VALUE OPTION_LIST)
  option(${OPTION_NAME} ${OPTION_STRING} ${OPTION_VALUE})
  list(APPEND ${OPTION_LIST} ${OPTION_NAME})
endmacro()

# ######### BEGIN SPARTA OPTIONAL DEPENDENCIES ##########
sparta_option(
  PKG_MPI_STUBS "Enable or disable sparta mpi stubs package. Default: OFF." OFF
  SPARTA_PKG_LIST)

sparta_option(PKG_FFT "Enable or disable sparta fft package. Default: OFF." OFF
              SPARTA_PKG_LIST)
sparta_option(
  PKG_KOKKOS "Enable or disable sparta kokkos package. Default: OFF." OFF
  SPARTA_PKG_LIST)

sparta_option(
  PKG_PYTHON "Enable or disable sparta python package. Default: OFF." OFF
  SPARTA_PKG_LIST)
# ######### END   SPARTA OPTIONAL DEPENDENCIES ##########

# ######### BEGIN SPARTA TPL DEPENDENCIES ##########
sparta_option(
  BUILD_MPI
  "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
  ON
  SPARTA_BUILD_TPL_LIST)

sparta_option(BUILD_KOKKOS "Enable or disable KOKKOS TPL. Default: OFF." OFF
              SPARTA_BUILD_TPL_LIST)

sparta_option(BUILD_JPEG "Enable or disable JPEG TPL. Default: OFF." OFF
              SPARTA_BUILD_TPL_LIST)

sparta_option(BUILD_PNG "Enable or disable PNG TPL. Default: OFF." OFF
              SPARTA_BUILD_TPL_LIST)

set(FFT "OFF" CACHE STRING "Select a FFT TPL from FFTW3, MKL, or KISS. Default: KISS.")
set(FFT_KOKKOS "OFF" CACHE STRING "Select a FFT TPL for Kokkos from CUFFT, HIPFFT, FFTW3, MKL, or KISS. Default: KISS.")
# ######### END   SPARTA TPL DEPENDENCIES ##########

# ######### BEGIN SPARTA EXTRA OPTIONS ##########
sparta_option(
  BUILD_SHARED_LIBS
  "Enable or disable building of sparta as a shared library. Default: OFF." OFF
  SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(
  SPARTA_MACHINE
  "Suffix to append to spa binary (WON'T enable any features automatically)" ""
  SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(SPARTA_ENABLE_ALL_PKGS "Enable all sparta packages. Default: OFF"
              "" SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(
  SPARTA_DISABLE_ALL_PKGS "Disable all sparta packages. Default: OFF" ""
  SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(
  SPARTA_LIST_PKGS "List available packages and quit. Default: OFF." ""
  SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(
  SPARTA_LIST_TPLS "List available packages and quit. Default: OFF." ""
  SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(SPARTA_ENABLE_TESTING "Enable sparta testing. Default: OFF" OFF
              SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(
  SPARTA_DSMC_TESTING_PATH "Enable sparta dsmc_testing. Default: OFF" OFF
  SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(
  SPARTA_CTEST_CONFIGS
  "Additional ctest configurations, separtaed by \"\;\", that allow SPARTA_SPA_ARGS_<CONFIG_NAME> or SPARTA_DSMC_TESTING_DRIVER_ARGS_<CONFIG_NAME> to be specified. Default: \"\""
  ""
  SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(
  SPARTA_SPA_ARGS "Additional arguments for the sparta binary. Default: OFF"
  OFF SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(
  SPARTA_DSMC_TESTING_DRIVER_ARGS
  "Additional arguments for ${SPARTA_DSMC_TESTING_PATH}/regression.py. Default: OFF"
  OFF
  SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(
  SPARTA_MULTIBUILD_CONFIGS
  "Additional sparta build configurations, separtaed by \"\;\", that allow multiple build configurations to be built from a single build directory. Default: \"\""
  ""
  SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(
  SPARTA_MULTIBUILD_PRESET_DIR
  "The path to your preset files at ${SPARTA_MULTIBUILD_PRESET_DIR}/${SPARTA_MULTIBUILD_CONFIG}.cmake. Default: \"\""
  ""
  SPARTA_EXTRA_OPTIONS_LIST)

sparta_option(SPARTA_ENABLE_PARAVIEW_TESTING "Enable ParaView testing. Default: OFF" OFF
              SPARTA_EXTRA_OPTIONS_LIST)

if(SPARTA_CTEST_CONFIGS)
  foreach(config ${SPARTA_CTEST_CONFIGS})
    list(APPEND SPARTA_EXTRA_OPTIONS_LIST SPARTA_SPA_ARGS_${config})
    list(APPEND SPARTA_EXTRA_OPTIONS_LIST
         SPARTA_DSMC_TESTING_DRIVER_ARGS_${config})
  endforeach()
endif()

# ######### END   SPARTA EXTRA OPTIONS ##########
