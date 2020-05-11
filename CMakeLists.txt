# TODO: Is 3.13 widley available??
cmake_minimum_required(VERSION 3.13.0)
project(SPARTA)

########## BEGIN SET COMMON DIRS ##########
set(SPARTA_TPL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(SPARTA_PKG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SPARTA_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(SPARTA_CMAKE_COMMON_DIR ${SPARTA_CMAKE_DIR}/common)
########## END   SET COMMON DIRS ##########

########## BEGIN SET CMAKE DEFAULTS ##########
include(${SPARTA_CMAKE_COMMON_DIR}/set/sparta_cmake_defaults.cmake)
########## END   SET CMAKE DEFAULTS ##########

########## BEGIN SET SPARTA OPTIONS ##########
include(${SPARTA_CMAKE_COMMON_DIR}/set/sparta_options.cmake)
########## END   SET SPARTA OPTIONS ##########

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
#set(CMAKE_CXX_FLAGS CACHE "${CMAKE_CXX_FLAGS} ${SPARTA_DEFAULT_CXX_COMPILE_FLAGS}" FORCE)
#add_compile_options(${SPARTA_DEFAULT_CXX_COMPILE_FLAGS})

include(${SPARTA_CMAKE_COMMON_DIR}/process/sparta_build_options.cmake)

add_subdirectory(src)
