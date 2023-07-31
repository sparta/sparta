# Copied from
# https://github.com/lammps/lammps/tree/3b77546eb9b0b5f7a402531b130e3b6dfad98a0c
# - Find mkl Find the native MKL headers and libraries.
#
# MKL_INCLUDE_DIRS - where to find mkl.h, etc. MKL_LIBRARIES    - List of
# libraries when using mkl. MKL_FOUND        - True if mkl found.
#

find_path(MKL_INCLUDE_DIR mkl_dfti.h HINTS $ENV{MKLROOT}/include)

find_library(
  MKL_LIBRARY
  NAMES mkl_rt
  HINTS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if all
# listed variables are TRUE

find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARY MKL_INCLUDE_DIR)

if(MKL_FOUND)
  set(MKL_LIBRARIES ${MKL_LIBRARY})
  set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})

  if(NOT TARGET MKL::MKL)
    add_library(MKL::MKL UNKNOWN IMPORTED)
    set_target_properties(
      MKL::MKL PROPERTIES IMPORTED_LOCATION "${MKL_LIBRARY}"
                          INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}")
  endif()
endif()

mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARY)
