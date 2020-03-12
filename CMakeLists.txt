# TODO: Is 3.13 widley available??
cmake_minimum_required(VERSION 3.13.0)
project(SPARTA)

########## BEGIN SET CMAKE DEFAULTS ##########

########## END   SET CMAKE DEFAULTS ##########

set(SPARTA_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(SPARTA_CMAKE_COMMON_DIR ${SPARTA_CMAKE_DIR}/common)

# Load in sparta's ENABLE/DISABLE options
include(${SPARTA_CMAKE_COMMON_DIR}/sparta_options.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/common/sparta_compiler_flags.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/machine_files/default.cmake)
# TODO: process SpArTA variables?

## ENABLE variables
## ENABLE variables will trigger a find package or local CMakeLists inclusion.
## If the find package fails, it
## should suggest setting env vars to the user's local installation of using a
## module to set env vars.

# TODO: process enable variables and set SPARTA_CMAKE_CONFIG_STRING string?
# TODO: Use TriBits? Would this help with TPL enables?

## Machine files
## Machine files can be specified via cmake -C machine.file, they will cache
## variables (as defaults) that can be optionally overriden. If a set comes
## after set-cache, the cache is ignored.
## EXAMPLES:
# chama.cmake
# cmake -C chama.cmake
# cmake -C chama.cmake -DSPARTA_ENABLE_MPICH
# cmake -C chama.cmake -DSPARTA_BIGBIG

# Setup sparta targets
set(TARGET_SPARTA spa_${SPARTA_CMAKE_CONFIG_STRING})
set(TARGET_SPARTA_MPI_STUBS mpi_stubs)

add_subdirectory(src)
#ADD_SUBDIRECTORY(lib)?
