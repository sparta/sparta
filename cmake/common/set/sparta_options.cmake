################################################################################
# This file sets sparta boolean options listed via `cmake -L`
################################################################################

#-------------------------------------------------------------------------------
# macro(sparta_option)
# Adds a sparta package option by invoking cmake's option routine and appends the option
# to OPTIONS_LIST
#
# Inputs:
#   OPTION_NAME
#   OPTION_STRING
#   OPTION_VALUE
#   OPTION_LIST
#
# Outputs:
#  OPTION_LIST
#  OPTION_NAME
macro(sparta_option OPTION_NAME OPTION_STRING OPTION_VALUE OPTION_LIST)
  option(${OPTION_NAME}
         ${OPTION_STRING}
         ${OPTION_VALUE}
        )
  list(APPEND ${OPTION_LIST} ${OPTION_NAME})
endmacro()


########## BEGIN SPARTA OPTIONAL DEPENDENCIES ##########
sparta_option(PKG_MPI_STUBS
              "Enable or disable sparta mpi stubs package. Default: OFF."
              OFF
              SPARTA_PKG_LIST
              )

sparta_option(PKG_FFT
              "Enable or disable sparta fft package. Default: OFF."
              OFF
              SPARTA_PKG_LIST
             )
sparta_option(PKG_KOKKOS
              "Enable or disable sparta kokkos package. Default: OFF."
              OFF
              SPARTA_PKG_LIST
             )
########## END   SPARTA OPTIONAL DEPENDENCIES ##########

########## BEGIN SPARTA TPL DEPENDENCIES ##########
sparta_option(BUILD_MPI
              "Enable or disable MPI TPL. Assumes environment has MPI_ROOT set. Default: ON."
              ON
              SPARTA_BUILD_TPL_LIST
             )

sparta_option(BUILD_KOKKOS
              "Enable or disable KOKKOS TPL. Default: OFF."
              OFF
              SPARTA_BUILD_TPL_LIST
             )

sparta_option(BUILD_JPEG
              "Enable or disable JPEG TPL. Default: OFF."
              OFF
              SPARTA_BUILD_TPL_LIST
             )

sparta_option(BUILD_PNG
              "Enable or disable PNG TPL. Default: OFF."
              OFF
              SPARTA_BUILD_TPL_LIST
             )

option(FFT
       "Select a FFT TPL from FFTW2, FFTW3, and MKL. Default: OFF."
       OFF
      )
########## END   SPARTA TPL DEPENDENCIES ##########

########## BEGIN SPARTA EXTRA OPTIONS ##########
sparta_option(BUILD_SHARED_LIBS
              "Enable or disable building of sparta as a shared library. Default: OFF."
              OFF
              SPARTA_EXTRA_OPTIONS_LIST
             )

sparta_option(SPARTA_MACHINE
              "Suffix to append to spa binary (WON'T enable any features automatically)"
              ""
              SPARTA_EXTRA_OPTIONS_LIST
              )
########## END   SPARTA EXTRA OPTIONS ##########
