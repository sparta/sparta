cmake_minimum_required(VERSION 3.10.0)

project(SPARTA)

set(SPARTA_BUILD_VER "1.0.0")
set(SPARTA_API_VER "1")

# ######### BEGIN CMAKE POLICY SETTINGS ##########
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  message(STATUS "Setting policy CMP0074 to use <TPL>_ROOT variables")
  cmake_policy(SET CMP0074 NEW)
endif()
# ######### END CMAKE POLICY SETTINGS ##########

# ######### BEGIN SET COMMON DIRS ##########
set(SPARTA_TPL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../lib)
set(SPARTA_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)
set(SPARTA_EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../examples)
set(SPARTA_DSMC_TESTING_EXAMPLES_DIR ${SPARTA_DSMC_TESTING_PATH}/examples)
set(SPARTA_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Modules)
set(SPARTA_CMAKE_COMMON_DIR ${SPARTA_CMAKE_DIR}/common)
# ######### END   SET COMMON DIRS ##########

# ######### BEGIN SET CMAKE DEFAULTS ##########
include(${SPARTA_CMAKE_COMMON_DIR}/set/sparta_cmake_defaults.cmake)
# ######### END   SET CMAKE DEFAULTS ##########

# ######### BEGIN SET SPARTA OPTIONS ##########
include(${SPARTA_CMAKE_COMMON_DIR}/set/sparta_options.cmake)
# ######### END   SET SPARTA OPTIONS ##########

# Include CTest in top-level CMakeLists. Tests are configured and added in
# examples/CMakeLists.txt
if(SPARTA_ENABLE_TESTING OR SPARTA_DSMC_TESTING_PATH)
  if(SPARTA_DSMC_TESTING_PATH)
    file(COPY ${SPARTA_DSMC_TESTING_PATH}/cmake/CTestConfig.cmake
         DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  include(CTest)
endif()

if(SPARTA_MULTIBUILD_CONFIGS)
  include(ExternalProject)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
       "# CMake generated Testfile for\n")
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
       "# Source directory: ${SPARTA_CMAKE_DIR}\n")
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
       "# Build directory: ${CMAKE_CURRENT_BINARY_DIR}\n")
  get_filename_component(SPARTA_MULTIBUILD_PRESET_DIR_ABSOLUTE
                         ${SPARTA_MULTIBUILD_PRESET_DIR} ABSOLUTE)

  foreach(BUILD_TYPE ${SPARTA_MULTIBUILD_CONFIGS})
    message("Adding ${BUILD_TYPE}...")
    ExternalProject_Add(
      ${BUILD_TYPE}
      SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
      CONFIGURE_COMMAND
        ${CMAKE_COMMAND} -DSPARTA_CMAKE_DIR=${SPARTA_CMAKE_DIR} -C
        ${SPARTA_MULTIBUILD_PRESET_DIR_ABSOLUTE}/${BUILD_TYPE}.cmake
        ${SPARTA_CMAKE_DIR}
      INSTALL_COMMAND "")
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
         "subdirs(\"${BUILD_TYPE}-prefix/src/${BUILD_TYPE}-build\")\n")
  endforeach()

  return()
endif()

set(TARGET_SPARTA spa_${SPARTA_MACHINE})

set(TARGET_SPARTA_LIB sparta_${SPARTA_MACHINE})
if(BUILD_SHARED_LIBS)
  set(TARGET_SPARTA_LIB_LINK "libsparta${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(TARGET_SPARTA_LIB_NAME
      "libsparta_${SPARTA_MACHINE}${CMAKE_SHARED_LIBRARY_SUFFIX}")
else()
  set(TARGET_SPARTA_LIB_LINK "libsparta${CMAKE_STATIC_LIBRARY_SUFFIX}")
  set(TARGET_SPARTA_LIB_NAME
      "libsparta_${SPARTA_MACHINE}${CMAKE_STATIC_LIBRARY_SUFFIX}")
endif()

# ######### BEGIN LIST AND EXIT ##########
if(SPARTA_LIST_PKGS OR SPARTA_LIST_TPLS)
  if(SPARTA_LIST_PKGS)
    message(STATUS "SPARTA Packages:")
    foreach(opt IN LISTS SPARTA_PKG_LIST)
      get_property(
        opt_help
        CACHE ${opt}
        PROPERTY HELPSTRING)
      message(STATUS "  // ${opt_help}")
      message(STATUS "  ${opt}=${${opt}}")
    endforeach()
    message(STATUS "")
  endif()

  if(SPARTA_LIST_TPLS)
    message(STATUS "SPARTA TPLs:")
    foreach(opt IN LISTS SPARTA_BUILD_TPL_LIST)
      get_property(
        opt_help
        CACHE ${opt}
        PROPERTY HELPSTRING)
      message(STATUS "  // ${opt_help}")
      message(STATUS "  ${opt}=${${opt}}")
    endforeach()
  endif()
  # EXIT
  return()
endif()
# ######### END LIST AND EXIT ##########

# ######### BEGIN PROCESS SPARTA OPTIONS ##########
include(${SPARTA_CMAKE_COMMON_DIR}/process/sparta_build_options.cmake)
# ######### END   PROCESS SPARTA OPTIONS ##########

# ######### BEGIN PRINT SPARTA OPTIONS ##########
include(${SPARTA_CMAKE_COMMON_DIR}/print/sparta_options.cmake)
# ######### END   PRINT SPARTA OPTIONS ##########

add_subdirectory(${SPARTA_SRC_DIR} ${CMAKE_CURRENT_BINARY_DIR}/src)
add_subdirectory(${SPARTA_EXAMPLES_DIR} ${CMAKE_CURRENT_BINARY_DIR}/examples)
if(SPARTA_DSMC_TESTING_PATH)
  add_subdirectory(${SPARTA_DSMC_TESTING_EXAMPLES_DIR}
                   ${CMAKE_CURRENT_BINARY_DIR}/dsmc_testing/examples)
endif()
