cmake_minimum_required(VERSION 3.10.0)
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

set(TARGET_SPARTA spa_${SPARTA_CMAKE_CONFIG_STRING})

include(${SPARTA_CMAKE_COMMON_DIR}/process/sparta_build_options.cmake)
include(${SPARTA_CMAKE_COMMON_DIR}/print/sparta_options.cmake)

add_subdirectory(src)
