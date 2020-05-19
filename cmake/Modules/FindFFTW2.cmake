# Copied from https://github.com/lammps/lammps/tree/3b77546eb9b0b5f7a402531b130e3b6dfad98a0c
# - Find fftw2
# Find the native double precision FFTW2 headers and libraries.
#
#  FFTW2_INCLUDE_DIRS  - where to find fftw2.h, etc.
#  FFTW2_LIBRARIES     - List of libraries when using fftw2.
#  FFTW2_OMP_LIBRARIES - List of libraries when using fftw2.
#  FFTW2_FOUND         - True if fftw2 found.
#

find_package(PkgConfig)

pkg_check_modules(PC_FFTW2 fftw2)
find_path(FFTW2_INCLUDE_DIR fftw2.h HINTS ${PC_FFTW2_INCLUDE_DIRS})
find_library(FFTW2_LIBRARY NAMES fftw2 HINTS ${PC_FFTW2_LIBRARY_DIRS})
find_library(FFTW2_OMP_LIBRARY NAMES fftw2_omp HINTS ${PC_FFTW2_LIBRARY_DIRS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW2_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(FFTW2 DEFAULT_MSG FFTW2_LIBRARY FFTW2_INCLUDE_DIR)

# Copy the results to the output variables and target.
if(FFTW2_FOUND)
  set(FFTW2_LIBRARIES ${FFTW2_LIBRARY} )
  set(FFTW2_INCLUDE_DIRS ${FFTW2_INCLUDE_DIR} )

  if(NOT TARGET FFTW2::FFTW2)
    add_library(FFTW2::FFTW2 UNKNOWN IMPORTED)
    set_target_properties(FFTW2::FFTW2 PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      IMPORTED_LOCATION "${FFTW2_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW2_INCLUDE_DIRS}")
  endif()
  if(FFTW2_OMP_LIBRARY)
    set(FFTW2_OMP_LIBRARIES ${FFTW2_OMP_LIBRARY})
    if(NOT TARGET FFTW2::FFTW2_OMP)
      add_library(FFTW2::FFTW2_OMP UNKNOWN IMPORTED)
      set_target_properties(FFTW2::FFTW2_OMP PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
	IMPORTED_LOCATION "${FFTW2_OMP_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW2_INCLUDE_DIRS}")
    endif()
  endif()
endif()

mark_as_advanced(FFTW2_INCLUDE_DIR FFTW2_LIBRARY FFTW2_OMP_LIBRARY)
